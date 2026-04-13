source(paste0(dirname(dirname(getwd())),'/map.r'))
source(paste0(HELP_DIR, "shortcuts.r"))
source(paste0(HELP_DIR, "helpers.r"))

library(patchwork)

library(gridExtra)

library(cowplot)

main <- readRDS(paste0(FIG_DIR, "volcano_main.rds"))
highlights <- readRDS(paste0(FIG_DIR, "highlights.rds"))
summaries <- readRDS(paste0(FIG_DIR, "cohort_summaries.rds"))
cohorts <- readRDS(paste0(FIG_DIR, "0_cohorts_full.rds"))

highlights <- plot_grid(highlights, labels = c("c"), label_size = 18)

main2 <- plot_grid(main, highlights, ncol = 1) 

main2_annotated <- main2 + plot_annotation(title = "Figure 2",theme = theme(plot.title = element_text(hjust = 0, size = 20, face = "bold")))

options(repr.plot.width = 15, repr.plot.height = 14)

main2_annotated

ggsave( paste0(FIG_DIR, "Main_Figure_2.png"), plot = main2, width = 15, height = 14)
ggsave( paste0(FIG_DIR, "Main_Figure_2.pdf"), plot = main2_annotated, width = 15, height = 14)

options(repr.plot.width = 26, repr.plot.height = 12)

supplementary1 <- 
plot_grid(cohorts, summaries, ncol = 2, 
          labels = c("a", "b"),
          label_size = 26,     # Font size
          label_x = 0,         # Left-align label
          label_y = 1,         # Top-align label
          hjust = -0.1,        # Slight shift left
          vjust = 1.2) +
plot_annotation(title = "Supplementary Figure 1",theme = theme(plot.title = element_text(hjust = 0, size = 26, face = "bold")))

options(repr.plot.width = 26, repr.plot.height = 12)
supplementary1
ggsave( paste0(FIG_DIR, "Supplementary_Figure_1.png"), plot = supplementary1, width = 26, height = 12)
ggsave( paste0(FIG_DIR, "Supplementary_Figure_1.pdf"), plot = supplementary1, width = 26, height = 12)
