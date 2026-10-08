#---------------------------------------------------------#
# Calculation of AEDs for Profiles Summary PODs
#---------------------------------------------------------#
# Katie Paul Friedman, Katie.PaulFriedman@UL.org
# August 29, 2026
# Updated September 8, 2026

#---------------------------------------------------------#
## Prepare environment
#---------------------------------------------------------#

library(data.table)
library(dplyr)
library(openxlsx)
library(tictoc)

packageVersion('httk') # was 2.7.4
# need to reverse the httk version to 2.6.1 to match AEDs from invitrodb v4.3
# update renv
# renv::remove("httk")
# renv::settings$use.cache(FALSE)
# renv::install("httk@2.6.1") # installed earlier version of httk to match httk version used with invitrodb v4.3
# renv::snapshot()
packageVersion('httk') # now 2.6.1

library(httk)
library(httr2)
library(purrr)
library(tibble)

#---------------------------------------------------------#
## Load data
#---------------------------------------------------------#

#rocks_profile_pod_table <- fread("./data/raw/rocks_profile_pod_table.csv")
rocks_profile_pod_table <- fread("./data/raw/rocks_profile_pod_table_08sept2026.csv")
# 1997 rows connected to cichem_id, which will not be understood by httk v2.6.1

# convert to dtxsid
# CIIW UI used to retrieve DTXSID from CICHEM ID
# cichem_map <- fread("./data/raw/identifier-tools_2026-08-29_07-03-21.csv")
cichem_map <- fread("./data/raw/identifier-tools_2026-09-08_20-08-06.csv")
# rename "output" to dtxsid
setnames(cichem_map, c('input','output'), c('cichem_id','dtxsid'))
# merge rocks_profile_pod_table with cichem_map to get the dsstox_substance_id for each chemical
rocks_profile_pod_table <- merge.data.table(rocks_profile_pod_table, 
                                            cichem_map, 
                                            by = "cichem_id",
                                            all.x = TRUE)
head(rocks_profile_pod_table)
# check for any NA dtxsid
rocks_profile_pod_table[is.na(dtxsid), .N] # 0, good

#---------------------------------------------------------#
## Run AEDs
#---------------------------------------------------------#

# load all the QSPRs/QSARs so that we can have a more complete chemical picture
httk::load_dawson2021()
httk::load_pradeep2020()
httk::load_sipes2017()
httk::load_honda2023()

# Operationalize calc_mc_oral_equiv()
set.seed(12345)

## Css50 first
tictoc::tic()

bioactive_cols <- c('min_pod', 'min_pod_httr')

for (col in bioactive_cols) {
  col_safe <- gsub("\\.", "", col) # For cleaner new column names
  new_col <- paste0(col_safe, "_aed_sumcl50")
  rocks_profile_pod_table[[new_col]] <- NA_real_ # Initialize
  
  for (i in seq_len(nrow(rocks_profile_pod_table))) {
    conc <- rocks_profile_pod_table[[col]][i]
    dtxsid <- rocks_profile_pod_table$dtxsid[i]
    
    if (!is.na(conc) && conc > 0 && !is.na(dtxsid) && nchar(dtxsid) > 0) {
      rocks_profile_pod_table[[new_col]][i] <- tryCatch({
        calc_mc_oral_equiv(
          conc = conc,
          dtxsid = dtxsid,
          which.quantile = 0.50,
          species = 'Human',
          restrictive.clearance = TRUE,
          output.units = 'mgpkgpday',
          model = 'sumclearances',
          Caco2.options = list(
            Caco2.Pab.default = 1.6,
            Caco2.Fabs = TRUE,
            Caco2.Fgut = TRUE,
            overwrite.invivo = FALSE,
            keepit100 = FALSE
          )
        )
      }, error = function(e) NA_real_)
    }
    # If NA or conc <= 0, skip (remains as NA)
  }
}

# Calc with Css95 

for (col in bioactive_cols) {
  col_safe <- gsub("\\.", "", col)
  new_col <- paste0(col_safe, "_aed_sumcl95")
  rocks_profile_pod_table[[new_col]] <- NA_real_ # Initialize
  
  for (i in seq_len(nrow(rocks_profile_pod_table))) {
    conc <- rocks_profile_pod_table[[col]][i]
    dtxsid <- rocks_profile_pod_table$dtxsid[i]
    
    if (!is.na(conc) && conc > 0 && !is.na(dtxsid) && nchar(dtxsid) > 0) {
      rocks_profile_pod_table[[new_col]][i] <- tryCatch({
        calc_mc_oral_equiv(
          conc = conc,
          dtxsid = dtxsid,
          which.quantile = 0.95,
          species = 'Human',
          restrictive.clearance = TRUE,
          output.units = 'mgpkgpday',
          model = 'sumclearances',
          Caco2.options = list(
            Caco2.Pab.default = 1.6,
            Caco2.Fabs = TRUE,
            Caco2.Fgut = TRUE,
            overwrite.invivo = FALSE,
            keepit100 = FALSE
          )
        )
      }, error = function(e) NA_real_)
    }
    # If NA or conc <= 0, skip (remains as NA)
  }
}

tictoc::toc() 

save(rocks_profile_pod_table, file = "./data/raw/rocks_profile_pod_table_aeds_08sept26.RData")
write.csv(rocks_profile_pod_table, file = "./data/raw/rocks_profile_pod_table_aeds_08sept26.csv", row.names = FALSE)


#---------------------------------------------------------#
## Verify AEDs
#---------------------------------------------------------#
# previous ROCKS profile table
rocks_profile_pod_table_old <- fread("./data/raw/rocks_profile_pod_table_aeds.csv")

# check rocks_profile_pod_table to see if columns are numeric
rocks_profile_pod_table %>% select(min_pod_aed_sumcl50, min_pod_aed_sumcl95, min_pod_httr_aed_sumcl50, min_pod_httr_aed_sumcl95) %>% str()
rocks_profile_pod_table_old %>% select(min_pod_aed_sumcl50, min_pod_aed_sumcl95) %>% str()

# merge to allow plot
m <- merge(
  rocks_profile_pod_table,
  rocks_profile_pod_table_old,
  by = "cichem_id",          
  suffixes = c("_new", "_old")
)

# verify small changes only with plot
plot(log10(m$min_pod_aed_sumcl50_new),
     log10(m$min_pod_aed_sumcl50_old),
     xlab = "New AEDs (sumcl50) [log10]",
     ylab = "Old AEDs (sumcl50) [log10]",
     main = "Comparison of New vs Old AEDs (sumcl50) [log10]")
abline(0, 1, col = "red", lty = 2)   # reference line for easy comparison


