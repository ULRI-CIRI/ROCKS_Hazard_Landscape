#' Install packages if not already installed
#'
#' Checks whether each package in \var{packages} is installed and installs
#' any that are missing using \code{\link[utils]{install.packages}}.
#'
#' @param packages A character vector of package names to check and install.
#' @param ... Additional arguments passed to \code{\link[utils]{install.packages}}.
#' @return Invisibly returns \code{NULL}.
#' @examples
#' install_if_missing(c("dplyr", "ggplot2", "tibble"))
install_if_missing <- function(packages, ...) {
  missing <- packages[!packages %in% rownames(utils::installed.packages())]
  if (length(missing) > 0) {
    message("Installing missing packages: ", paste(missing, collapse = ", "))
    utils::install.packages(missing, ...)
  }
  invisible(NULL)
}

#' Get a PostgreSQL database connection from the environment variable specified
#'
#' Returns a database connection from the connection parameters in the environment varialbe passed.
#' The environment variable is expected to have a valid PostgreSQL connection string. As an example,
#' it could look something like this: PIPELINES_CONN_STR=postgres://username:password@host:port/db
#' All five fields (username, password, host, port, database) must be specified.  See details on the
#' connection string here: https://www.postgresql.org/docs/current/libpq-connect.html
#'
#' @param enviro_variable The environments variable which contains the database connection
#' @return An active database connection that will be closed when the application exists
pg_connection_from_env <- function(enviro_variable) {
  conn_str <- Sys.getenv(enviro_variable, unset = "empty")
  if (conn_str == "empty") {
    stop(
      "Environment variable ",
      enviro_variable,
      " not found"
    )
  }
  parts <- regmatches(
    conn_str,
    regexec(
      "postgres://([^:]+):([^@]+)@([^:/]+)(?::(\\d+))?/(.+)",
      conn_str,
      perl = TRUE
    )
  )[[1]]
  conn <- DBI::dbConnect(
    RPostgres::Postgres(),
    user = parts[2],
    password = parts[3],
    host = parts[4],
    port = if (nchar(parts[5]) > 0) as.integer(parts[5]) else 5432L,
    dbname = parts[6]
  )
  # ensure the database connection is closed at the end of the R session
  reg.finalizer(globalenv(), function(e) DBI::dbDisconnect(conn), onexit = TRUE)
  return(conn)
}

#' Escape a vector of values for use in a SQL IN clause or array literal
#'
#' Quotes and escapes each element using ANSI SQL literal escaping, then
#' collapses them into a single comma-separated string safe for interpolation
#' into a SQL query. The result can be dropped directly into an IN clause:
#' \code{WHERE id IN (pg_escape_vector(ids))}.
#'
#' When \code{array_type} is specified, the result is wrapped in a PostgreSQL
#' array literal suitable for passing to functions and stored procedures that
#' expect an array argument, e.g. \code{TEXT[]}.
#'
#' @param vector A vector of values to escape and collapse.
#' @param array_type Optional. When provided, wraps the result in a PostgreSQL
#'   array literal cast to this type, e.g. \code{"TEXT"} produces
#'   \code{ARRAY['a', 'b']::TEXT[]}.
#' @return A \code{DBI::SQL} object containing the escaped, comma-separated
#'   values, or an array literal if \code{array_type} is specified.
#' @examples
#' pg_escape_vector(c("red", "blue", "green"))
#' # SQL: 'red', 'blue', 'green'
#'
#' pg_escape_vector(c("CICHEM000200014", "CICHEM000200048"), array_type = "TEXT")
#' # SQL: ARRAY['CICHEM000200014', 'CICHEM000200048']::TEXT[]
pg_escape_vector <- function(vector, array_type = NULL) {
  ansi <- DBI::ANSI()
  collapsed <- paste(
    sapply(
      vector,
      function(x) DBI::dbQuoteLiteral(ansi, x)
    ),
    collapse = ", "
  )
  if (!is.null(array_type)) {
    return(DBI::SQL(paste0("ARRAY[", collapsed, "]::", array_type, "[]")))
  }
  return(DBI::SQL(collapsed))
}

#' Query PostgreSQL database and return a tibble
#'
#' Interpolates named parameters into a SQL query string using
#' \code{DBI::sqlInterpolate}, executes it against the provided connection,
#' and returns the result as a \code{\link[tibble]{tibble}}.
#'
#' Use \code{?name} placeholders in \var{query_str} for each named argument
#' passed via \code{...}. Values are safely escaped via \code{DBI::ANSI()},
#' so the function is safe against SQL injection as long as all user-supplied
#' values are passed as named arguments rather than interpolated into the
#' query string directly.
#'
#' By default values are interpolated using dbQuoteLiteral - this should be used
#' when quoting content that would be in rows in a database. If you'd instead be
#' passing database identifiers (tables columns, schema) as parameters, use
#' DBI::dbQuoteIdentifier to the argument you're passing into this function.
#'
#' @param conn A \code{DBIConnection} object.
#' @param query_str A SQL query string containing \code{?name} placeholders.
#' @param schema An optional character string specifying the PostgreSQL schema
#'   to search. If provided, the session \code{search_path} is temporarily set
#'   to this schema for the duration of the query and then restored afterward.
#' @param ... Named arguments whose values will be substituted for the
#'   corresponding \code{?name} placeholders in \var{query_str}. Use
#'   \code{\link{pg_escape_vector}} to pass multi-value IN clause parameters.
#' @return A \code{\link[tibble]{tibble}} containing the query results.
#' @examples
#' pg_query(conn, "SELECT * FROM t WHERE id = ?id AND ver = ?ver", id = 42, ver = "2.9")
#' pg_query(conn, "SELECT * FROM t", schema = "my_schema")
pg_query <- function(conn, query_str, schema = NULL, ...) {
  ansi <- DBI::ANSI()
  # change schema, just for this query, if specified
  if (!is.null(schema)) {
    original_schema <- DBI::dbGetQuery(conn, "SHOW search_path")[[1]]
    if (original_schema == "") {
      original_schema <- "public"
    }
    schema_query <- "SET search_path TO ?"
    old_schema <- DBI::sqlInterpolate(
      ansi,
      schema_query,
      DBI::dbQuoteIdentifier(ansi, original_schema)
    )
    new_schema <- DBI::sqlInterpolate(
      ansi,
      schema_query,
      DBI::dbQuoteIdentifier(ansi, schema)
    )
    DBI::dbExecute(conn, new_schema)
    on.exit(DBI::dbExecute(conn, old_schema), add = TRUE)
  }
  # prepare query using proper interpolation
  query <- DBI::sqlInterpolate(ansi, query_str, ...)
  # return a tibble
  return(tibble::as_tibble(DBI::dbGetQuery(conn, query)))
}