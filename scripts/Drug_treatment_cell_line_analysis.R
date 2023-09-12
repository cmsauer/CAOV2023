# script to analyse drug treatment data for ovarian cancer cell lines
# Carolin Sauer

# load libraries
require(tidyverse)
require(GRmetrics)

source("cms_pal.R")

# Load data

GRmet <- read.csv("../data/cell_line_drug_data/combined_GRmetrics_output_CMS20211214.csv")

CLs_select <- read.csv("../data/cell_line_drug_data/cell_line_datasheet.csv")  %>%
  filter(str_detect(use_experiment, "exp")) %>% 
  mutate(identifier = paste0(Cell.Line, "_", use_experiment))

Drug_info <- read.csv("../data/cell_line_drug_data/Drug_info.csv")

drug_cols <- c(CMSroyalblue, CMSdarkgrey, CMSlightred, CMSteal, CMSorange, CMSpurple)

# process/prep data
analysis_data <- GRmet %>% 
  filter(identifier %in% CLs_select$identifier) %>% 
  select(cell_line, identifier, fit_GR, agent, ctrl_cell_doublings, GR50, GRmax) %>% 
  left_join(CLs_select, by = "identifier") %>% 
  left_join(Drug_info, by = "agent") %>% 
  arrange(cell_line, group)

analysis_data$CA_group <- factor(analysis_data$CA_group, levels = c("low", "high"))
analysis_data$agent <- str_replace(analysis_data$agent, "Paclitaxel", "paclitaxel")
analysis_data$agent <- factor(analysis_data$agent, levels = c("oxaliplatin", "paclitaxel", "tozasertib", "alisertib", "barasertib", "cw069", "bi2536", "volasertib", "cfi400945", "az3146", "gsk923295"))


# PLOTTING

## Comparison of all drugs across all cell lines

format <- scales::format_format(big.mark = ",", decimal.mark = ".", scientific = FALSE)

gr50_all <- analysis_data %>%
  filter(cell_line != "A2780-ADR",
         cell_line != "OV90(mesmasson)",
         GR50 > 0) %>%
  ggplot(aes(x = agent , y = GR50, colour = group))  +
  geom_jitter(size = 0.5, alpha = 0.6, position = position_jitterdodge(jitter.width = 0.5)) +
  scale_y_continuous(trans = "log10", labels = format) +
  # scale_y_continuous(trans = "log10", labels = scales::comma) +
  geom_boxplot(notch = T, outlier.colour = NA, fill = NA, size = 0.25) +
  scale_colour_manual(values = drug_cols, name = "Drug group/target") +
  ggpubr::stat_compare_means(size = 1.5, label.x = 8.25, label.y = log10(10000)) +
  labs(x = "", y = "Potency (GR50 (µM))") +
  theme_bw(7) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none", 
        legend.key.size = unit(0.1, "in"))

grmax_all <- analysis_data %>% 
  filter(cell_line != "A2780-ADR",
         cell_line != "OV90(mesmasson)") %>%
  ggplot(aes(x = agent , y = GRmax, colour = group))  +
  geom_hline(yintercept = 0, colour = CMSgrey, size = 0.2, linetype = "dashed") +
  geom_jitter(size = 0.5, alpha = 0.6, position = position_jitterdodge(jitter.width = 0.5)) +
  geom_boxplot(notch = T, outlier.colour = NA, fill = NA, size = 0.25) +
  scale_colour_manual(values = drug_cols, name = "Drug group/target") +
  ggpubr::stat_compare_means(size = 1.5, label.x = 8.25, label.y = 0.75) +
  labs(x = "", y = "Efficacy (GRmax)") +
  ylim(-1, 1) +
  theme_bw(7) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        legend.position = "none", 
        legend.key.size = unit(0.1, "in"))

cowplot::plot_grid(gr50_all, grmax_all, nrow = 2, align = "v")
ggsave("../figures/CAOV_paper_figure9b.pdf",width = 2.75, height = 3)


# Correlation with CA

gr50 <- analysis_data %>% 
  filter(cell_line != "A2780-ADR") %>%
  filter(GR50 != Inf) %>% 
  ggplot(aes(x = X..CA.cells, y = (GR50))) +
  scale_y_continuous(trans = "log10", labels = format) +
  # scale_y_continuous(trans = "log10", labels = scales::comma) +
  geom_point(aes(colour = group), size = 0.5, alpha = 0.6) +
  geom_smooth(method = "lm", colour = CMSdarkgrey, size = 0.5) +
  scale_colour_manual(values = drug_cols, name = "Drug group/target") +
  facet_wrap(.~agent, nrow = 1) +
  ggpubr::stat_cor(method = "spearman", size = 1.5) +
  labs(x = "\n% non-mitotic cells with 2 or more centrosomes", y = "Potency (GR50 (µM))") +
  theme_bw(7) +
  theme(strip.background = element_blank(),
        legend.position = "none", 
        legend.key.size = unit(0.1, "in"))

grmax <- analysis_data %>% 
  filter(cell_line != "A2780-ADR",
         cell_line != "OV90(mesmasson)") %>%
  ggplot(aes(x = X..CA.cells, y = GRmax)) +
  geom_hline(yintercept = 0, colour = CMSgrey, size = 0.2, linetype = "dashed")  +
  geom_point(aes(colour = group), size = 0.5, alpha = 0.6) +
  geom_smooth(method = "lm", colour = CMSdarkgrey, size = 0.5) +
  scale_colour_manual(values = drug_cols, name = "Drug group/target") +
  facet_wrap(.~agent, nrow = 1) +
  ggpubr::stat_cor(method = "spearman", size = 1.5) +
  labs(x = "\n% non-mitotic cells with 2 or more centrosomes", y = "Efficacy (GRmax)") +
  ylim(-1, 1) +
  theme_bw(7) +
  theme(strip.background = element_blank(),
        legend.position = "none", 
        legend.key.size = unit(0.1, "in"))

cowplot::plot_grid(gr50, grmax, nrow = 2, align = "v")
# ggsave("../figures/CAOV_paper_figure9c.pdf",width = 7.25, height = 3)
# write.csv(analysis_data, "../source_data_figures/source_data_figure9.csv", row.names = F)


