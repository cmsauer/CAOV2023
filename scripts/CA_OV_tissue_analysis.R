# script to analyse CA imaging data for OV04 and BriTROC-1 tissue samples
# Carolin Sauer

# load libraries
require(tidyverse)
require(broom)
require(survival)
require(survminer)

source("cms_pal.R")

# load data
### CA data from BriTROC and OV04 (CAOV) samples have been normalised as described and outliers have been removed prior to this analysis.
CAOV_BriTROC_norm_joined <- read.csv("../data/tissue_imaging_data/Joined_normDATA_BriTROC_CAOV_20220107.csv")

### Raw data per object (centrosome) detected for each imaging field and tissue. Data needed for spatial analyses
OV04_eval <- read.csv("../data/tissue_imaging_data/OV04_all_Eval_RmOutliers_FTnorm_20220106.csv") %>% 
  unite("identifier", c(Tissue.Section, Field), sep = "_", remove = T)
OV04_spatials <- read.table("../data/tissue_imaging_data/OV04_spatials.txt", sep = "\t", header = T )
OV04_spatials <- OV04_eval %>% 
  left_join(OV04_spatials, by = "identifier") 

BriTROC_eval <- read_csv("../data/tissue_imaging_data/BriTROC_data_all_joined_RmOutliers_CTRLnorm_20220107.csv") %>% 
  filter(SampleType == "QCpassed_Sample") %>% 
  rename(TMA = "Tissue Section") %>% 
  unite("identifier", c(BriTROC_core, TMA, Field), sep = "_", remove = F)
BriTROC_spatials <- read.table("../data/tissue_imaging_data/BriTROC_spatials.txt", sep = "\t", header = T) 
BriTROC_spatials <- BriTROC_eval %>% 
  left_join(BriTROC_spatials, by = "identifier")

### Clinical data and metadata information
## Sample data
caov_samples <- read.csv("../data/clinical_and_meta_data/caov_samples_tissueSites_20220110.csv") %>% mutate(cohort = "OV04")
britroc_samples <- read.csv("../data/clinical_and_meta_data/britroc_samples_tissueSites_20221122.csv") %>% mutate(cohort="BriTROC")
samples_all <- bind_rows(caov_samples, britroc_samples)

## Per patient Data
britroc_patients <- read.csv("../data/clinical_and_meta_data/britroc_clinDat_perPatient_20221122.csv")
caov_patients <- read.csv("../data/clinical_and_meta_data/caov_clinDat_perPatient_20220110.csv")
patients_all <- bind_rows(caov_patients, britroc_patients)

### CA threshold based on the 95% confidence interval of the FT sample with the highest CA score
upperCI <- 1.833243 


# FIGURE 1
## Fig 1d
CA_tissue <- CAOV_BriTROC_norm_joined %>% 
  select(SampleID, SampleType, Cohort, CAScoreMedian, CSizeMedian, heterosce) %>% 
  unique() %>% 
  mutate(SampleType = case_when(SampleType == "Control" ~ "Normal (other)",
                                SampleType == "FT" ~ "Normal (FT)",
                                SampleType == "Liver_Control" ~ "Liver",
                                SampleType == "QCpassed_Sample" ~ "HGSOC (BriTROC)",
                                SampleType == "T" ~ "HGSOC (OV04)")) %>% 
  mutate(SampleType2 = case_when(SampleType %in% c("Normal (other)","Normal (FT)") ~ "Normal",
                                 SampleType == "Liver" ~ "Liver",
                                 SampleType == "HGSOC (BriTROC)" ~ "Tumour \n(BriTROC)",
                                 SampleType == "HGSOC (OV04)" ~ "Tumour \n(OV04)"))

CA_tissue$SampleType2 = factor(CA_tissue$SampleType2, levels = c("Normal", "Tumour \n(OV04)", "Tumour \n(BriTROC)", "Liver"))

boxplot <- CA_tissue %>%
  ggplot(aes(x = SampleType2 , y = CAScoreMedian)) +
  geom_jitter(aes(colour = SampleType), width = 0.2, size = 0.5) +
  scale_colour_manual(values = c("Normal (other)" = CMSpink, "Normal (FT)" = CMSorange, 
                                 "Liver" = CMSpurple, "HGSOC (BriTROC)" = CMSblue, "HGSOC (OV04)" = CMSlightblue), name = "Sample type") +
  geom_boxplot(notch = T, outlier.colour = NA, fill = NA, size = 0.5) +
  geom_hline(yintercept = upperCI, colour = CMSdarkgrey, linetype = "dashed") +
  annotate("text", x = 1, y = upperCI + 0.4, label = "CA cutoff", size = 2) +
  theme_classic(7) +
  ylim(0,10) +
  labs(x = "Sample type", y = "CA score (median)") +
  ggpubr::stat_compare_means(paired = F, var.equal = FALSE, size = 1.75, label.y = 10, label.x = 1.25) +
  ggpubr::stat_compare_means(comparisons = list(c("Normal", "Tumour \n(OV04)"), 
                                                c("Normal", "Tumour \n(BriTROC)"),
                                                c("Normal", "Liver")), size = 1.75, label.y = c(7,8,9)) +
  theme(legend.position = "none",
        legend.key.size = unit(0.1, "in"),
        axis.title.x = element_blank())

dens <- CA_tissue %>%
  ggplot(aes(y = CAScoreMedian)) +
  # scale_y_continuous(breaks=seq(0, 55, 10), limits=c(0, 55)) +
  geom_density(aes(fill = SampleType, colour = SampleType), alpha = 0.5) +
  scale_fill_manual(values = c("Normal (other)" = CMSpink, "Normal (FT)" = CMSorange, 
                               "Liver" = CMSpurple, "HGSOC (BriTROC)" = CMSblue, "HGSOC (OV04)" = CMSlightblue), name = "Sample type") +
  scale_colour_manual(values = c("Normal (other)" = CMSpink, "Normal (FT)" = CMSorange, 
                                 "Liver" = CMSpurple, "HGSOC (BriTROC)" = CMSblue, "HGSOC (OV04)" = CMSlightblue), name = "Sample type") +
  theme_void(7) +
  ylim(0,10) +
  theme(legend.position = "none")

cowplot::plot_grid(boxplot, dens, rows = 1, align = "h", rel_widths = c(1,0.25))
# ggsave("../figures/CAOV_paper_figure1d.pdf", height = 2, width = 3)

## Fig 1e
CA_tissue %>%
  ggplot(aes(x = SampleType2 , y = CSizeMedian)) +
  geom_jitter(aes(colour = SampleType), width = 0.2, size = 0.5) +
  scale_colour_manual(values = c("Normal (other)" = CMSpink, "Normal (FT)" = CMSorange, 
                                 "Liver" = CMSpurple, "HGSOC (BriTROC)" = CMSblue, "HGSOC (OV04)" = CMSlightblue), name = "Sample type") +
  geom_boxplot(notch = T, outlier.colour = NA, fill = NA, size = 0.5) +
  theme_classic(7) +
  labs(x = "Sample type", y = "Centrosome size (median)") +
  # ylim(0.5,1.5)+
  ggpubr::stat_compare_means(paired = F, var.equal = FALSE, size = 1.75, label.x = 1.25, label.y = 1.4) +
  ggpubr::stat_compare_means(comparisons = list(c("Normal", "Tumour \n(OV04)"), 
                                                c("Normal", "Tumour \n(BriTROC)"),
                                                c("Normal", "Liver")), size = 1.75, label.y = c(1.25,1.3,1.35)) +
  theme(legend.position = "right",
        legend.key.size = unit(0.1, "in"),
        axis.title.x = element_blank())

# ggsave("../figures/CAOV_paper_figure1e.pdf", height = 2, width = 3)
# write.csv(CA_tissue, "../source_data_figures/source_data_figure1.csv",row.names = F)


# FIGURE 2
CAOV_BriTROC_norm_joined %>% 
  group_by(SampleID) %>%
  mutate(SampleType = case_when(SampleType == "Control" ~ "Normal (other)",
                                SampleType == "FT" ~ "Normal (FT)",
                                SampleType == "Liver_Control" ~ "Liver",
                                SampleType == "QCpassed_Sample" ~ "HGSOC (BriTROC)",
                                SampleType == "T" ~ "HGSOC (OV04)")) %>% 
  ggplot(aes(y = reorder(SampleID, -CAScoreMedian) , x = CAscore_norm)) +
  geom_boxplot(notch = T, outlier.shape = NA, aes(fill = SampleType), lwd = 0.2) +
  scale_fill_manual(values = c("Normal (other)" = CMSpink, "Normal (FT)" = CMSorange, 
                               "Liver" = CMSpurple, "HGSOC (BriTROC)" = CMSblue, "HGSOC (OV04)" = CMSlightblue), name = "Sample type") +
  geom_vline(xintercept = upperCI, colour = CMSdarkgrey, linetype = "dashed") +
  annotate("text", y = 310, x = 2.25, label = "CA cutoff\n \n63.8% with CA \n(183/287 HGSOC)", size = 2, hjust = 0) +
  scale_x_continuous(limits=c(0, 16), breaks=c(0,1,2,3,4,5,6,8,10,12,14,16)) +
  theme_classic(8) +
  labs(y = "Samples", x = "CA Score") +
  theme(axis.text.y = element_blank(),
        legend.position = "bottom",
        legend.key.size = unit(0.1, "in"))
# ggsave("../figures/CAOV_paper_figure2.pdf", height = 10, width = 5.5)
# write.csv(CAOV_BriTROC_norm_joined, "../source_data/source_data_figure2.csv",row.names = F)


# FIGURE 3
## Fig3a-b and Supplementary Figure 1
## Spatial Heatmaps - OV04 (CAOV)

for (i in unique(OV04_spatials$Tissue.Section)) {
  temp = OV04_spatials %>% 
    filter(Tissue.Section == i) %>% 
    unique() %>% 
    mutate(lim_CAscore_norm = if_else(CAscore_norm < 0.5, 0.5, CAscore_norm)) %>% 
    mutate(lim_CAscore_norm = if_else(CAscore_norm > 8, 8, lim_CAscore_norm)) 
  
  #If want to normalise/re-scale positions    
  # group_by(Tissue.Section) %>%
  #     mutate(minX = min(Position.X..µm.),
  #          minY = min(Position.Y..µm.),
  #          X = Position.X..µm.-minX,
  #          Y = Position.Y..µm.-minY)
  
  temp %>% 
    ggplot(aes(x = Position.X..µm., 
               y = Position.Y..µm., 
               colour = lim_CAscore_norm)) +
    # scale_x_continuous(limits = c(-6550, 8250)) +
    # scale_y_continuous(limits = c(-5000, 13750)) +
    geom_point(size = 0.001) +
    scale_colour_gradient(low=CMSpurple,
                          high=CMSyellow, 
                          "Imaging field \nCA score",
                          breaks = c(2, 4, 6, 8),
                          limits = c(0,8),
                          labels = c(" 2"," 4"," 6",">8")) +
    labs(title = paste0(temp[1,]$Tissue.Section, "-", temp[1,]$Tissue.Type, "  ", "CA score = ", round(temp[1,]$CAScoreMedian, digits = 2), "\nHeterog. score = ", round(temp[1,]$heterosce, digits = 2)),
         x = "Position X (µm)",
         y = "Position Y (µm)") +
    theme_classic(6) +
    theme(legend.key.size = unit(0.1, "in"))
  
  # ggsave(paste0("../figures/spatialHeatmaps/spatialHeatmaps_OV04/spatialHeatmap_OV04_", i, ".png"), height = 2.5, width = 3.5)
  cat("Plot generated for ", i, "\n", sep = "")  
}

## Spatial Heatmaps - BriTROC
for (i in unique(BriTROC_spatials$BriTROC_core)) {
  temp = BriTROC_spatials %>% 
    filter(BriTROC_core == i) %>% 
    unique() %>% 
    mutate(lim_CAScoreMedian = if_else(CAscore_norm < 0.5, 0.5, CAscore_norm)) %>% 
    mutate(lim_CAScoreMedian = if_else(CAscore_norm > 8, 8, lim_CAScoreMedian)) %>% 
    group_by(Core) %>% 
    mutate(minX = min(Position.X..µm.),
           minY = min(Position.Y..µm.)) %>% 
    mutate(X = case_when(Core == "Core 1" ~ Position.X..µm.-minX - 3000,
                         Core == "Core 2" ~ Position.X..µm.-minX ,
                         Core == "Core 3" ~ Position.X..µm.-minX + 3000),
           Y = case_when(Core == "Core 1" ~ Position.Y..µm.-minY ,
                         Core == "Core 2" ~ Position.Y..µm.-minY ,
                         Core == "Core 3" ~ Position.Y..µm.-minY )) %>% 
    mutate(BriTROC_core = str_remove(BriTROC_core, "_"))
  
  temp  %>% 
    ggplot(aes(x = X, 
               y = Y, 
               colour = lim_CAScoreMedian)) +
    scale_x_continuous(limits = c(-3500, 4000)) +
    scale_y_continuous(limits = c(-2000, 2000)) +
    geom_point(size = 0.001) +
    annotate("text", x = c(-2600, 400, 3400), y = 1250, label = c("Core 1", "Core 2", "Core 3"), size = 2) +
    scale_colour_gradient(low=CMSpurple,
                          high=CMSyellow, 
                          "Imaging field \nCA score",
                          breaks = c(2, 4, 6, 8),
                          limits = c(0,8),
                          labels = c(" 2"," 4"," 6",">8")) +
    labs(title = paste0(temp[1,]$BriTROC_core, "-T", temp[1,]$Tissue.Type, "  ", "CA score = ", round(temp[1,]$CAScoreMedian, digits = 2), 
                        "\nHeterog. score = ", round(temp[1,]$heterosce, digits = 2)),
         x = "Position X (µm)",
         y = "Position Y (µm)") +
    theme_classic(6) +
    theme(legend.key.size = unit(0.1, "in"))
  
  # ggsave(paste0("../figures/spatialHeatmaps/spatialHeatmaps_BriTROC/spatialHeatmap_scaled_BriTROC_", i, ".png"), height = 2.5, width = 3.5)s
  cat("Plot generated for ", i, "\n", sep = "")   
}


## Fig3c - Compare CA scores in patients with multiple tissues collected/imaged
mlt_samp_patients <- samples_all %>%
  filter(n > 1) 

mlt_samp_patients  %>% 
  mutate(site_categ = factor(site_categ, levels = c("ovary/fallopian tube","omentum","multiple tissues \n(primary site)",
                                                    "metastatic","unknown"))) %>% 
  ggplot(aes(x = factor(patient_id), y = CAScoreMedian)) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = upperCI, ymax = Inf,
           alpha = 0.3, fill = CMSgrey) +
  geom_hline(yintercept = upperCI, linetype = "dashed", colour = "black") +
  geom_line(aes(group = patient_id), colour = CMSdarkgrey, size =0.25) +
  geom_point(aes(shape = site_categ,colour = cohort),  size = 1) +
  geom_point(data = unique(mlt_samp_patients[,c("patient_id", "av_CAScore")]), 
             aes(x = factor(patient_id), y = av_CAScore), colour = CMSlightred, shape = 8, size = 0.7) +
  scale_shape_discrete(name = "Tissue site") +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  # ylim(0,8) +
  labs(x = "patients (OV04+BriTROC)", y = "CA score") +
  theme_classic(7) +
  theme(legend.key.size = unit(0.05, "in"),
        legend.position = "bottom",
        legend.box = "vertical",
        axis.text.x = element_blank()) +
  guides(shape=guide_legend(nrow=2,byrow=TRUE),
         colour=guide_legend(nrow=1,byrow=TRUE))
# ggsave("../figures/CAOV_paper_figure3c.pdf", height = 2.5, width = 3)

#All samples tissue site comparison
## Fig3d - CA score
samples_all %>%
  mutate(site_categ = factor(site_categ, levels = c("ovary/fallopian tube","omentum","multiple tissues \n(primary site)",
                                                    "metastatic","unknown"))) %>% 
  ggplot(aes(x = site_categ, y = CAScoreMedian)) +
  geom_jitter(aes(shape = site_categ, colour = cohort),width = 0.2,  size = 1) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  # ylim(0,8) +
  geom_hline(yintercept = upperCI, linetype = "dashed", colour = "black") +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  theme_classic(7) +
  guides(x = guide_axis(angle = 90)) +
  ylab("CA score") +
  theme(axis.title.x = element_blank(), 
        legend.key.size = unit(0.05, "in"),
        legend.position = c(0.9,0.9)) +
  ggpubr::stat_compare_means(size = 1.5) +
  guides(shape="none")
# ggsave("../figures/CAOV_paper_figure3d.pdf", height = 2, width = 2)

## Fig3e - Heterosc score
samples_all %>%
  mutate(site_categ = factor(site_categ, levels = c("ovary/fallopian tube","omentum","multiple tissues \n(primary site)",
                                                    "metastatic","unknown"))) %>% 
  ggplot(aes(x = site_categ, y = heterosce)) +
  geom_jitter(aes(shape = site_categ, colour = cohort),width = 0.2,  size = 1) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  theme_classic(7) +
  guides(x = guide_axis(angle = 90)) +
  ylab("Tissue heterogeneity \n(CA score)") +
  theme(axis.title.x = element_blank(), 
        legend.key.size = unit(0.05, "in"),
        legend.position = c(0.9,0.9)) +
  ggpubr::stat_compare_means(size = 1.5) +
  guides(shape="none")
# ggsave("../figures/CAOV_paper_figure3e.pdf", height = 2, width = 2)
# write.csv(samples_all, "../source_data_figures/source_data_figure3cde.csv", row.names = F)


# FIGURE 4
## Fig4a - Histology
patients_all %>% 
  filter(histological_type != "") %>% 
  mutate(histological_type = case_when(str_detect(histological_type, "serous") ~ "HGSOC",
                                       str_detect(histological_type, "endom") ~ "endometrioid")) %>% 
  ggplot(aes(x = histological_type, y = CAScore)) +
  geom_jitter(width = 0.2, aes(colour = cohort), size = 0.5) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  theme_classic(7) +
  ylab("CA score") +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  ggpubr::stat_compare_means(size = 2) +
  guides(x = guide_axis(angle = 90)) 
# ggsave("../figures/CAOV_paper_figure4a.pdf", height = 1.925, width = 1.25)

## Fig4b - Stage
patients_all %>% 
  filter(!is.na(tumour_stage_at_diagnosis)) %>% 
  mutate(tumour_stage_at_diagnosis = factor(tumour_stage_at_diagnosis, levels = c("1","2","3","4"))) %>% 
  ggplot(aes(x = tumour_stage_at_diagnosis, y = CAScore)) +
  geom_jitter(width = 0.2, aes(colour = cohort), size = 0.5) +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  theme_classic(7) +
  ylab("CA score") +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  ggpubr::stat_compare_means(size = 2)
# ggsave("../figures/CAOV_paper_figure4b.pdf", height = 1.5, width = 2)

## Fig4c - BRCA germline
patients_all %>% 
  mutate(germline = if_else(germline == "", "unknown", as.character(germline))) %>% 
  ggplot(aes(x = germline, y = CAScore)) +
  geom_jitter(width = 0.2, aes(colour = cohort), size = 0.5) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  theme_classic(7) +
  ylab("CA score") +
  theme(axis.title.x = element_blank(), legend.position = "none") +
  ggpubr::stat_compare_means(size = 2) +
  guides(x = guide_axis(angle = 90)) 
# ggsave("../figures/CAOV_paper_figure4c.pdf", height = 1.925, width = 1.5)

## Fig4d-f - Surgery type (IPS vs DPS) 
### CA score
patients_all %>% 
  filter(debulking_surgery_type %in% c("DPS", "IPS")) %>% 
  ggplot(aes(x = factor(debulking_surgery_type, levels = c("IPS", "DPS")), y = CAScore)) +
  geom_jitter(aes(colour = cohort), width = 0.2, size = 0.5) +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  theme_classic(7) +
  ylab("CA score") +
  theme(axis.title.x = element_blank(), 
        legend.key.size = unit(0.05, "in"),
        legend.position = "none") +
  ggpubr::stat_compare_means(size = 2)
# ggsave("../figures/CAOV_paper_figure4d.pdf", height = 1.5, width = 1.25)

### Heterosc
patients_all %>% 
  filter(debulking_surgery_type %in% c("DPS", "IPS")) %>% 
  ggplot(aes(x = factor(debulking_surgery_type, levels = c("IPS", "DPS")), y = Heterosc)) +
  geom_jitter(aes(colour = cohort), width = 0.2, size = 0.5) +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  theme_classic(7) +
  ylab("Tissue heterogeneity") +
  theme(axis.title.x = element_blank(), 
        legend.key.size = unit(0.05, "in"),
        legend.position = "none") +
  ggpubr::stat_compare_means(size = 2)
# ggsave("../figures/CAOV_paper_figure4e.pdf", height = 1.5, width = 1.25)

### Centrosome Size
patients_all %>% 
  filter(debulking_surgery_type %in% c("DPS", "IPS")) %>% 
  ggplot(aes(x = factor(debulking_surgery_type, levels = c("IPS", "DPS")), y = CSize)) +
  geom_jitter(aes(colour = cohort), width = 0.2, size = 0.5) +
  scale_colour_manual(values = c("BriTROC"=CMSblue, "OV04"=CMSlightblue)) +
  geom_boxplot(fill = NA, outlier.colour = NA, notch = T) +
  theme_classic(7) +
  ylab("Centrosome size") +
  theme(axis.title.x = element_blank(), 
        legend.key.size = unit(0.05, "in"),
        legend.position = "none") +
  ggpubr::stat_compare_means(size = 2)
# ggsave("../figures/CAOV_paper_figure4f.pdf", height = 1.5, width = 1.25)

## Fig4g-h - Cox PH - Overall free survival
### OV04
surv_data_pretty <- patients_all %>% 
  filter(cohort == "OV04") %>%
  mutate(CA_status = str_replace(CA_status, "^CA$", "significant CA")) %>% 
  select(patient_id, Stage = tumour_stage_at_diagnosis, "Age (standardised)" = StandardisedAge, Surgery = debulking_surgery_type, 
         "CA status" = CA_status, days_OS, surv_status, csize_stat, Hetero_stat, CA_status2, status, status2)
surv_data_pretty$`CA status` <- factor(surv_data_pretty$`CA status`, levels = c("no CA", "significant CA"))
surv_data_pretty$Surgery <- factor(surv_data_pretty$Surgery, levels = c("IPS", "DPS"))

res.cox.CAScore.os <- coxph(Surv(days_OS, surv_status) ~ `CA status` + Surgery + Stage + `Age (standardised)`, data = surv_data_pretty)
summary(res.cox.CAScore.os)

ggforest(res.cox.CAScore.os, data = surv_data_pretty)
# ggsave("../figures/CAOV_paper_figure4g.pdf", width = 7, height = 4)

### BriTROC
surv_data_pretty <- patients_all %>% 
  filter(cohort == "BriTROC") %>%
  mutate(CA_status = str_replace(CA_status, "^CA$", "significant CA")) %>% 
  select(patient_id, Stage = tumour_stage_at_diagnosis, "Age (standardised)" = StandardisedAge, Surgery = debulking_surgery_type, 
         "CA status" = CA_status, days_OS, surv_status, csize_stat, Hetero_stat, CA_status2, status, status2)
surv_data_pretty$`CA status` <- factor(surv_data_pretty$`CA status`, levels = c("no CA", "significant CA"))
surv_data_pretty$Surgery <- factor(surv_data_pretty$Surgery, levels = c("IPS", "DPS"))

res.cox.CAScore.os <- coxph(Surv(days_OS, surv_status) ~ `CA status` + Surgery + Stage + `Age (standardised)`, data = surv_data_pretty)
summary(res.cox.CAScore.os)

ggforest(res.cox.CAScore.os, data = surv_data_pretty)
# ggsave("../figures/CAOV_paper_figure4h.pdf", width = 7, height = 4)
# write.csv(patients_all, "../source_data_figures/source_data_figure4.csv", row.names = F)

