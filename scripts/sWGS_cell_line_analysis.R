# script to analyse sWGS and copy number signature data for ovarian cancer cell lines
# Carolin Sauer

# load libraries
require(tidyverse)
require(data.table)
require(shades)

source("cms_pal.R")

# Load data
## meta data and cell line characterisation (imaging) data
meta <- read_csv("../data/cell_line_genomics_data/cell-info_metadata_CMS20211128.csv") %>% 
  rename(sample_id = DNA_JBLAB_number)

eval_data <- read.csv("../data/cell_line_imaging_data/CellLine_summary_data_stain1-4_V2_CMS20220523.csv") 

cns <- read_csv("../data/cell_line_genomics_data/signature_exposures_CMS20211130.csv")

# Correlations with ploidy and tMAD score (Fig8a-d)
## Read in genomic features
gen_feat <- read.csv("../data/cell_line_genomics_data/cell_line_genomic_features.csv")
gen_feat <- gen_feat %>% 
  left_join(meta, by = "sample_id") %>% 
  left_join(eval_data, by = "Cell.Type")

gen_feat %>% 
  ggplot(aes(x = X..CA.cells, y = ploidy)) +
  geom_point(colour = CMSgrey, size = 0.5, alpha = 0.85) +
  geom_smooth(method = "lm", size = 0.5, colour = CMSdarkgrey) +
  ggpubr::stat_cor(method = "spearman", size = 2) +
  theme_bw(7) +
  labs(x = "CA frequency", y = "Ploidy") 
ggsave("../figures/CAOV_paper_figure8a.pdf", height = 2, width = 1.75)
write.csv(gen_feat, "../source_data_figures/source_data_figure8a.csv", row.names = F)

# tMAD scores
tmad <- read.csv("../data/cell_line_genomics_data/tMAD_3Million_500_control_CONTROL_noMedianNorm.csv") 

tmad_comb <- tmad %>% 
  filter(sample_id %in% gen_feat$sample_id) %>% # those samples that passed sequencing (and CNS analysis) QC
  left_join(meta[c(1:3,8)], by = "sample_id") %>%
  left_join(eval_data, by = "Cell.Type") 

tmad_comb %>% 
  filter(subtype != "normal") %>%
  ggplot(aes(x=X..CA.cells, y=tMAD)) +
  geom_point(size = 0.5, colour = CMSgrey, alpha = 0.85) +
  geom_smooth(method = "lm", colour = CMSgreen, size = 0.5) +
  ggpubr::stat_cor(size = 2 ,method = "spearman") +
  labs(x = "CA frequency", y = "tMAD score") +
  theme_bw(6)
# ggsave("../figures/CAOV_paper_figure8b.pdf", height = 2, width = 1.75)

tmad_comb %>% 
  filter(subtype != "normal") %>%
  ggplot(aes(x=X..MN.cells, y=tMAD)) +
  scale_x_continuous(trans = "log1p") +
  geom_point(size = 0.5, colour = CMSgrey, alpha = 0.85) +
  geom_smooth(method = "lm", colour = CMSgreen, size = 0.5) +
  ggpubr::stat_cor(size = 2 ,method = "spearman") +
  labs(x = "MN frequency", y = "tMAD score") +
  theme_bw(6)
# ggsave("../figures/CAOV_paper_figure8c.pdf", height = 2, width = 1.75)

tmad_comb %>% 
  filter(subtype != "normal") %>%
  ggplot(aes(x=X..CA.cells*X..MN.cells, y=tMAD)) +
  scale_x_continuous(trans = "log1p",
                     labels = scales::unit_format(unit = "k",scale = 1/1000, accuracy = 0.5,sep = "")) +
  geom_point(size = 0.5, colour = CMSgrey, alpha = 0.85) +
  geom_smooth(method = "lm", colour = CMSgreen, size = 0.5) +
  ggpubr::stat_cor(size = 2 ,method = "spearman") +
  labs(x = "CA × MN frequency", y = "tMAD score") +
  theme_bw(6) +
  theme(axis.text.x = element_text(colour = "white"))
# ggsave("../figures/CAOV_paper_figure8d.pdf", height = 2, width = 1.75)
# write.csv(tmad_comb, "../source_data_figures/source_data_figure8b-d.csv", row.names = F)

# Copy Number Signatures
## combine all data for CNS analysis
comb <- meta %>% 
  left_join(eval_data, by = "Cell.Type") %>% 
  left_join(cns, by = "sample_id")

## Prep data for ALR model run on cluster (see separate script)

cns_IF_data <- comb %>% 
  filter(!is.na(s1)) %>% 
  mutate(CA_group = case_when(X..CA.cells > median(eval_data$X..CA.cells) ~ "high",
                              X..CA.cells < median(eval_data$X..CA.cells) ~ "low"),
         MN_group = case_when(X..MN.cells > median(eval_data$X..MN.cells) ~ "high",
                              X..MN.cells < median(eval_data$X..MN.cells) ~ "low")) %>% 
  select(sample_id, s1,s2,s3,s4,s5,s6,s7, X..CA.cells, CA_group, X..MN.cells, MN_group)
# write.csv(cns_IF_data, "../data/cell_line_genomics_data/cns_IF_data_forALRmodel_CMS20220912.csv", row.names = F, na = "")



# Correlations with CNS (Fig8e-f)
### mean CA Score
comb %>% 
  filter(!is.na(Cell.Type),
         !is.na(s1)) %>% 
  gather(17:23, key = "Signature", value = "Exposure")  %>% 
  unique() %>% 
  group_by(Signature) %>% 
  ggplot(aes(x = X..CA.cells, y = Exposure)) +
  geom_point(aes(colour = Signature), size = 0.75) +
  scale_colour_manual(values = saturation(CMS_pal, 0.5)) +
  geom_smooth(method = "lm", colour = CMSdarkgrey, size = 0.5) +
  facet_wrap(~Signature, nrow = 1) +
  ggpubr::stat_cor(size = 2, method = "spearman", label.sep = "\n") +
  labs(x = "Mean % non-mitotic cells with 2 or more centrosomes", y = "Copy number signature activity") +
  theme_bw(7) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))
# ggsave("../figures/CAOV_paper_figure8e.pdf", height = 2, width = 4)

### mean MN fraction
comb %>% 
  filter(!is.na(Cell.Type),
         !is.na(s1)) %>% 
  gather(17:23, key = "Signature", value = "Exposure")  %>% 
  unique() %>% 
  group_by(Signature) %>% 
  ggplot(aes(x = X..MN.cells, y = Exposure)) +
  geom_point(aes(colour = Signature), size = 0.75) +
  scale_colour_manual(values = saturation(CMS_pal, 0.5)) +
  scale_y_continuous(breaks = seq(0,1,0.25)) +
  geom_smooth(method = "lm", colour = CMSdarkgrey, size = 0.5) +
  facet_wrap(~Signature, nrow = 1) +
  ggpubr::stat_cor(size = 2, method = "spearman", label.sep = "\n") +
  labs(x = "Mean % non-mitotic cells with micronuclei", y = "Copy number signature activity") +
  theme_bw(7) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"))
# ggsave("../figures/CAOV_paper_figure8f.pdf", height = 2, width = 4)
# write.csv(comb, "../source_data_figures/source_data_figure8e-f.csv", row.names = F)

## code for ALR model (Fig8g-h) provided in separate script. (see scripts/ALR_analysis/CNS_analyses_FE_no_zeros_CAandMN_CMS20211227.R)

# Genome subclonality
## Calculate genome subclonality (also see rascal paper)
rcn_data <- read.table("../data/cell_line_genomics_data/copyNumberSegmented_30kb_HGSOC_CellLines.txt", header = T, sep = "\t") 

rcn <- rcn_data %>% 
  filter(! chromosome %in% c("X", "Y")) %>%
  mutate(chromosome = factor(chromosome, levels = unique(chromosome))) %>%
  mutate_at(vars(start, end), as.integer) %>%
  mutate(position = round((start + end) / 2))

library(rascal)
# copy number values should be relative to the average copy number within the sample so
# we scale this by the median value
rcn <- rcn %>%
  group_by(sample) %>% 
  mutate_at(vars(copy_number, segmented), ~ . / median(segmented, na.rm = TRUE))

segments <- rcn %>% 
  group_by(sample) %>% 
  filter(!is.na(segmented)) %>%
  mutate(length = end - start + 1) %>% 
  arrange(sample, chromosome, start) %>% 
  mutate(new_segment = row_number() == 1 | !(sample == lag(sample) & chromosome == lag(chromosome) & segmented == lag(segmented))) %>% 
  mutate(segment = cumsum(new_segment)) %>% 
  group_by(sample, segment) %>% #otherwise don't get all samples returned (as not coded into function)
  summarize(sample = first(sample), 
            chromosome = first(chromosome), start = first(start), 
            end = last(end), copy_number = first(segmented), bin_count = n(), 
            sum_of_bin_lengths = sum(length), weight = sum(length)/median(length))

chromosomes <- chromosome_offsets(rcn)

rcn <- convert_to_genomic_coordinates(rcn, "position", chromosomes)

segments <- convert_to_genomic_coordinates(segments, c("start", "end"), chromosomes)

## convert to absolute using previously estimated fits
fits <-  gen_feat %>% select(sample, sample_id, ploidy, cellularity, distance) %>%
  mutate(sample = str_replace(sample, "IGROV_1", "IGROV-1")) # fix sample name
acn_segments <-  fits %>% 
  left_join(segments, by = "sample") %>% 
  mutate(acn_segment = relative_to_absolute_copy_number(relative_copy_numbers = copy_number, ploidy = ploidy, cellularity = cellularity))

#determine subclonality
subclonal_segs_cl <- acn_segments %>%
  filter(!is.na(acn_segment)) %>% 
  mutate(deviation = abs(acn_segment - round(acn_segment))) %>% 
  mutate(clonality = if_else(deviation > 0.25, "subclonal", "clonal")) 

subclonal_summary_cl <- subclonal_segs_cl %>% 
  group_by(sample) %>% 
  mutate(genome_length = sum(sum_of_bin_lengths),
         seg_count =  max(segment)) %>%
  filter(clonality == "subclonal") %>% 
  add_count(sample) %>% 
  mutate(clonality_lenght = sum(sum_of_bin_lengths),
         fraction_genome_subclonal = (clonality_lenght/genome_length)*100,
         fraction_segments_subclonal = (n/seg_count)*100) %>% 
  select(sample, ploidy, cellularity, distance, clonality, fraction_genome_subclonal, fraction_segments_subclonal) %>% 
  unique() %>% 
  ungroup()

# plotting data 
subclon_plotting <- tmad_comb %>% left_join(subclonal_summary_cl[,c("sample", "fraction_genome_subclonal")], by = c("sample_id" = "sample"))

subclon_plotting <- gen_feat %>% 
  select(sample)
  left_join(subclonal_summary_cl[,c("sample", "fraction_genome_subclonal")], by = "sample")

subclon_plotting <- gen_feat %>% 
  mutate(sample = str_replace(sample, "IGROV_1", "IGROV-1")) %>%
  select(sample, sample_id, cell_line, Cell.Type, X..CA.cells, X..MN.cells) %>%
  left_join(subclonal_summary_cl[,c("sample", "fraction_genome_subclonal")], by = "sample")

cor.test(subclon_plotting$X..CA.cells, subclon_plotting$fraction_genome_subclonal, method = "spearman")
cor.test(subclon_plotting$X..MN.cells, subclon_plotting$fraction_genome_subclonal, method = "spearman")

# plotting

#CA
subclon_plotting %>% 
  ggplot(aes(x = X..CA.cells, y = fraction_genome_subclonal)) +
  geom_point(colour = CMSgrey, size = 0.5, alpha = 0.85) +
  geom_smooth(method = "lm", size = 0.5, colour = CMSorange) +
  ggpubr::stat_cor(method = "spearman", size = 2) +
  theme_bw(7) +
  labs(x = "% Non-mitotic cells \nwith 2 or more centrosomes", y = "% Genome subclonality") 
# ggsave("../figures/CAOV_paper_figure8i_p1.pdf", height = 2, width = 1.5)

subclon_plotting %>% 
  filter(!is.na(X..CA.cells)) %>% 
  mutate(CA_group = case_when(X..CA.cells > 32 ~ "high", 
                              X..CA.cells < 21 ~"low",
                              X..CA.cells <= 32 & X..CA.cells >= 21 ~ "mid")) %>% 
  ggplot(aes(x = factor(CA_group, levels = c("low", "mid", "high")), y = fraction_genome_subclonal)) +
  geom_jitter(colour = CMSgrey, size = 0.5, width = 0.25, alpha = 0.85) + 
  geom_boxplot(notch = T, fill = NA, outlier.colour = NA, colour = CMSdarkgrey)+
  ggpubr::stat_compare_means(size = 1.75, label.y = 45, label.x = 0.9) +
  ggpubr::stat_compare_means(size = 1.75, comparisons = list(c("low", "mid"), c("high", "mid"), c("high", "low")), label = "p.signif") +
  labs(x = "Centrosome amplification", y = "% Genome subclonality") +
  theme_bw(7) +
  theme(legend.key.size = unit(0.1, "in"),
        legend.position = "right")
# ggsave("../figures/CAOV_paper_figure8i_p2.pdf", height = 2, width = 1.5)


#MN
subclon_plotting %>% 
  ggplot(aes(x = X..MN.cells, y = fraction_genome_subclonal)) +
  geom_point(colour = CMSgrey, size = 0.5, alpha = 0.85) +
  geom_smooth(method = "lm", size = 0.5, colour = CMSorange) +
  ggpubr::stat_cor(method = "spearman", size = 2) +
  scale_x_continuous(trans = "log1p") +
  theme_bw(7) +
  labs(x = "% Non-mitotic cells \nwith micronuclei", y = "% Genome subclonality") 
# ggsave("../figures/CAOV_paper_figure8j_p1.pdf", height = 2, width = 1.5)

subclon_plotting %>% 
  filter(!is.na(X..MN.cells)) %>% 
  mutate(MN_group = case_when(X..MN.cells > 19 ~ "high", 
                              X..MN.cells < 10 ~"low",
                              X..MN.cells <= 19 & X..MN.cells >= 10 ~ "mid")) %>% 
  ggplot(aes(x = factor(MN_group, levels = c("low", "mid", "high")), y = fraction_genome_subclonal)) +
  geom_jitter(colour = CMSgrey, size = 0.5, width = 0.25, alpha = 0.85) + 
  geom_boxplot(notch = T, fill = NA, outlier.colour = NA, colour = CMSdarkgrey)+
  ggpubr::stat_compare_means(size = 1.75, label.y = 45, label.x = 0.9) +
  ggpubr::stat_compare_means(size = 1.75, comparisons = list(c("low", "mid"), c("high", "mid"), c("high", "low")), label = "p.signif") +
  labs(x = "Micronuclei", y = "% Genome subclonality") +
  theme_bw(7) +
  theme(legend.key.size = unit(0.1, "in"),
        legend.position = "right")
# ggsave("../figures/CAOV_paper_figure8j_p2.pdf", height = 2, width = 1.5)
# write.csv(subclon_plotting, "../source_data_figures/source_data_figure8i-j.csv", row.names = F)







