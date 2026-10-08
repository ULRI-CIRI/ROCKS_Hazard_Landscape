
library(DBI)
library(RPostgres)

connect_tamdb <- function() {
  DBI::dbConnect(
    RPostgres::Postgres(),
    dbname   = Sys.getenv("DB_NAME"),
    host     = Sys.getenv("DB_HOST"),
    port     = as.integer(Sys.getenv("DB_PORT")),
    user     = Sys.getenv("DB_USER"),
    password = Sys.getenv("DB_PASSWORD")
  )
}
