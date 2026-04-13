source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

### NEO_DIR <- INPUT DIRECTORY SPECIFIED IN map.r
### SOM_DIR <- INPUT DIRECTORY SPECIFIED IN map.r
### TMP <- OUTPUT DIRECTORY SPECIFIED IN map.r

get_fp <- function (i, type = "purity") {
  if (type == "neo_pep") { ac(paste0(NEO_DIR, i, "/", i, ".neo.peptide_scores.tsv.gz"))} 
}

reader <- function (i_file = "ab", sample = "cd") { fread(i_file) %>% mutate(sampleId = sample)}

patients <- list.files(NEO_DIR)

lr_threshold <- .01

out_file <- file.path(TMP_DIR, "neo_pep.csv")
first <- TRUE

system.time(
for (i in patients[1:10]) {
  
  print(i); flush.console();
  i_file <- get_fp(i, type = "neo_pep")

  if (!file.exists(i_file)) next

  oo <- tryCatch({
    oo <- reader(i_file, sample = i)
    oo <- oo[oo$LikelihoodRank < lr_threshold, ]
    oo
  }, error = function(e) NULL)

  if (!is.null(oo) && nrow(oo) > 0) {
    fwrite(oo, out_file, append = !first)
    first <- FALSE
  }
})
