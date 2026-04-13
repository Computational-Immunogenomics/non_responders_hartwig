source(paste0(dirname(getwd()),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))

neo <- fread( paste0(TMP_DIR, "neo.csv"))

neo_ready <- 
neo %>% 
 fi( VariantCopyNumber > .5, SubclonalLikelihood < .2) %>% 
 gb( sampleId ) %>% 
 su( neo_ct = n(), neo_rna_ct = sum(RnaFrags > 0) ) %>%
 mu( neo_rna_ct = ifelse(neo_rna_ct == 0, NA, neo_rna_ct))

neo_pep <- fread( paste0(TMP_DIR, "neo_pep.csv"))

neo_pep_ready <- 
neo_pep %>% 
 gb( sampleId ) %>% 
 su( neo_pep_ct = n(), neo_pep_rna_ct = sum(ExpLikelihoodRank < .02) ) 

neo_go <- neo_ready %>% lj(neo_pep_ready, by = "sampleId")

fwrite(neo_go, paste0(READY_DIR, "neo_ready.csv"))
