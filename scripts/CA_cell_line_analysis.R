# script to analyse CA imaging data for ovarian cancer cell lines
# Carolin Sauer

# load libraries
require(tidyverse)
require(GGally)

source("cms_pal.R")

# load and prepare data

## Stain 1 results - main evaluation data
EvalData_FieldData <- read.csv("../data/cell_line_imaging_data/Stain1_EvalData_FieldData.csv") 
Field_info <- read.csv("../data/cell_line_imaging_data/Stain1_EvalData_FieldData_FieldSelect_NumberNuclei.csv") 
### Mitotic Index
Mitotic <- read.csv("../data/cell_line_imaging_data/Stain1_MitosisData_FieldData.csv") 
mitoticIndex_cell <- Mitotic %>% 
  group_by(Cell.Type, MitosisStatus) %>% 
  mutate(Index_mean = mean(Index)) %>% 
  select(Cell.Type, MitosisStatus, Index_mean) %>%
  unique() %>% 
  spread(key = MitosisStatus, value = Index_mean) %>% 
  mutate(mitotic = if_else(is.na(mitotic), 0, mitotic))

### CA Score
CellCAscore <- EvalData_FieldData %>% 
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells),
         meanCAScore = mean(CAScore)) %>%
  select(Cell.Type, med_Fraction_CAcells, meanCAScore) %>% 
  unique()

s1_n <- sum(EvalData_FieldData$Nuclei.Number)
nucs_s1 <- Field_info[, c("Cell.Type", "Nuclei.Number")] %>% rename(numberNuclei = Nuclei.Number) %>% mutate(stain = "s1")
s1_c <- sum(EvalData_FieldData$Centrosome.Number)

## Stain 2 results
stain2_eval <- read.csv("../data/cell_line_imaging_data/all_PlateResults_stain2_CMS20211123.csv")  

stain2_eval_edit <- stain2_eval %>% 
  filter(Nuclei.Selected...Number.of.Objects >= 150) %>% #well results (total of 30 cells per field minimum)
  unite("Well_ID", c(Plate,Row,Column), sep = ".", remove = F) %>%
  group_by(Cell.Type) %>% 
  select(Well_ID, 
         Cell.Type,
         numberNuclei = Nuclei.Selected...Number.of.Objects,
         numberCentrosomes = Centrosomes...Number.of.Objects,
         numberCETN3posCentrosomes = Centrin3.pos.Centrosomes...Number.of.Objects,
         fracCETN3pos = Centrosomes...Centrin3.pos.Centrosomes...Mean.per.Well) %>% 
  mutate(FracPCMfragmentation = 1-fracCETN3pos) %>% 
  unique()

stain2_eval_edit_cell <- stain2_eval_edit %>% 
  group_by(Cell.Type) %>% 
  mutate(mean_fracCETN3pos = mean(fracCETN3pos),
         mean_FracPCMfragmentation = mean(FracPCMfragmentation)) %>% 
  select(Cell.Type, "% Centrin3 positive centrosomes"=mean_fracCETN3pos, "% Acentriolar centrosomes"=mean_FracPCMfragmentation) %>% 
  unique() 

CETN3 <- EvalData_FieldData %>%  
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells)) %>% 
  select(Cell.Type, med_Fraction_CAcells) %>% 
  unique() %>% 
  left_join(stain2_eval_edit_cell, by = "Cell.Type")

s2_n <- sum(stain2_eval_edit$numberNuclei)
nucs_s2 <- stain2_eval_edit[, c("Cell.Type", "numberNuclei")] %>% mutate(stain = "s2")
s2_c <- sum(stain2_eval_edit$numberCentrosomes)

## Stain 3 results
stain3_eval <- read.csv("../data/cell_line_imaging_data/all_PlateResults_stain3_CMS20211126.csv") 

stain3_eval_edit <- stain3_eval %>% 
  filter(Nuclei.Selected...Number.of.Objects >= 75) %>% 
  unite("Well_ID", c(Plate,Row,Column), sep = ".", remove = F) %>%
  group_by(Cell.Type) %>% 
  select(Well_ID, 
         Cell.Type,
         numberNuclei = Nuclei.Selected...Number.of.Objects,
         numberCentrosomes = Centrosomes...Number.of.Objects,
         numberCep164posCentrosomes = Cep164.pos.Centrosomes...Number.of.Objects,
         fracCep164pos = Centrosomes...Cep164.pos.Centrosomes...Mean.per.Well) %>% 
  mutate(fracCep164neg = 1-fracCep164pos) %>% 
  unique()

stain3_eval_edit_cell <- stain3_eval_edit %>% 
  group_by(Cell.Type) %>% 
  mutate(mean_fracCep164pos = mean(fracCep164pos),
         mean_fracCep164neg = mean(fracCep164neg)) %>% 
  select(Cell.Type, "% Cep164 positive centrosomes"=mean_fracCep164pos, "% Acentriolar centrosomes"=mean_fracCep164neg) %>% 
  unique() 

Cep164 <- EvalData_FieldData %>%  
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells)) %>% 
  select(Cell.Type, med_Fraction_CAcells) %>% 
  unique() %>% 
  left_join(stain3_eval_edit_cell, by = "Cell.Type")

s3_n <- sum(stain3_eval_edit$numberNuclei)
nucs_s3 <- stain3_eval_edit[, c("Cell.Type", "numberNuclei")] %>% mutate(stain = "s3")
s3_c <- sum(stain3_eval_edit$numberCentrosomes)

## Stain 4 results
stain4_eval <- read.csv("../data/cell_line_imaging_data/all_PlateResults_stain4_CMS20211127.csv") 
stain4_eval_edit <- stain4_eval %>% 
  filter(Nuclei.Selected...Number.of.Objects >= 75) %>% 
  unite("Well_ID", c(Plate,Row,Column), sep = ".", remove = F) %>%
  group_by(Cell.Type) %>% 
  select(Well_ID, 
         Cell.Type,
         numberNuclei = Nuclei.Selected...Number.of.Objects,
         numberCentrosomes = Centrosomes...Number.of.Objects,
         gH2AX_meanIntensity = Nuclei.Selected...Intensity.gH2AX.Mean...Mean.per.Well,
         numbergH2AXposNuclei = gH2AX.positive...Number.of.Objects,
         fracgH2AXpos = Nuclei.Selected...gH2AX.positive...Mean.per.Well) %>%  
  unique()

stain4_eval_edit_cell <- stain4_eval_edit %>% 
  group_by(Cell.Type) %>% 
  mutate(mean_fracgH2AXpos = mean(fracgH2AXpos),
         mean_gH2AXIntensity = mean(gH2AX_meanIntensity)) %>% 
  select(Cell.Type, mean_fracgH2AXpos, mean_gH2AXIntensity) %>% 
  unique() 

gH2AX <- EvalData_FieldData %>%  
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells)) %>% 
  select(Cell.Type, med_Fraction_CAcells) %>% 
  unique() %>% 
  left_join(stain4_eval_edit_cell, by = "Cell.Type")

s4_n <- sum(stain4_eval_edit$numberNuclei)
nucs_s4 <- stain4_eval_edit[, c("Cell.Type", "numberNuclei")] %>% mutate(stain = "s4")
s4_c <- sum(stain4_eval_edit$numberCentrosomes)


# Total number of nuclei and centrosomes
tot_n <- sum(s1_n, s2_n, s3_n, s4_n)
tot_c <- sum(s1_c, s2_c, s3_c, s4_c)

all_nucs <- bind_rows(nucs_s1, nucs_s2, nucs_s3, nucs_s4) %>%
  arrange(Cell.Type) %>%
  group_by(Cell.Type, stain) %>%
  mutate(Nuclei.sum_perStain = sum(numberNuclei)) %>%
  select(-numberNuclei) %>%
  unique() %>%
  group_by(Cell.Type) %>%
  mutate(mean_Nuclei.sum = round(mean(Nuclei.sum_perStain), digits = 0)) %>%
  select(-c(stain, Nuclei.sum_perStain)) %>% unique()

# Combined cell line data
Combined_perCL <- EvalData_FieldData %>% 
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells),
         meanFraction_CAcells = mean(Fraction_CAcells),
         meanCAScore = mean(CAScore),
         meanFraction_MNcells = mean(Fraction_MNcells)) %>% 
  select(Cell.Type, med_Fraction_CAcells, meanFraction_CAcells, meanCAScore, meanFraction_MNcells) %>% 
  unique()  %>% 
  left_join(mitoticIndex_cell[,1:2], by = "Cell.Type") %>% 
  left_join(Cep164[,c(1,3)], by = "Cell.Type") %>% 
  left_join(gH2AX[,c(1,3,4)], by = "Cell.Type")

# cell line meta data
meta_info <- read.csv("../data/cell_line_imaging_data/cell-info_metadata_CMS20211128.csv") %>% 
  mutate(Cell.Type = as.character(Cell.Type)) 


# FIGURE 5
CAfrac <- EvalData_FieldData %>% 
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells),
         meanCAScore = mean(CAScore)) %>%
  ggplot(aes( y = reorder(Cell.Type, -med_Fraction_CAcells), x = Fraction_CAcells)) +
  geom_bar(CellCAscore, mapping = aes(x = meanCAScore*40, y = reorder(Cell.Type, -med_Fraction_CAcells)), 
           fill = CMSlightblue, alpha = 0.5, width = 0.75, stat = "identity") +
  scale_x_continuous(sec.axis = sec_axis(~ ./40, name = "Mean CA score")) +
  geom_jitter(colour = CMSgrey, size = 0.5) +
  geom_boxplot(fill = NA, outlier.colour = NA) +
  theme_classic(6) +
  labs(x = "% Non-mitotic cells \nwith 2 or more centrosomes") +
  theme(legend.position = "none",
        axis.title.x.top = element_text(color = CMSlightblue),
        axis.text.x.top = element_text(color = CMSlightblue),
        axis.ticks.x.top = element_line(color = CMSlightblue),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(colour = "#d9d9d9", size = 0.25)) 

#number of Nuclei
NoNuclei <- EvalData_FieldData %>%  
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells)) %>% 
  select(Cell.Type, med_Fraction_CAcells) %>% 
  unique() %>% 
  # left_join(Field_info, by = "Cell.Type") %>%  #old...
  left_join(all_nucs, by = "Cell.Type") %>%
  select(Cell.Type, med_Fraction_CAcells, mean_Nuclei.sum) %>% 
  unique() %>% 
  left_join(meta_info, by = "Cell.Type") %>%
  ggplot(aes( y = reorder(Cell.Type, -med_Fraction_CAcells), x = mean_Nuclei.sum)) +
  geom_bar(stat = "identity", aes(fill = subtype2), width = 0.75) +
  scale_fill_manual(values = c("HGSOC"=CMSteal, "unknown"=CMSroyalblue, "normal"=CMSdarkgrey, "LGSOC" = CMSorange, "endometrioid" = CMSpink, "mucinous" = CMSlightred), name = "subtype") +
  scale_x_continuous(breaks=seq(0,6000,2000),
                     labels = scales::unit_format(unit = "k",scale = 1/1000, accuracy = 1, sep = "")) +
  theme_classic(6) +
  labs(x = "Avg. # nuclei \nanalysed") +
  theme(legend.position = "none",
        legend.key.size = unit(0.1, "in"),
        axis.title.y = element_blank(),
        panel.grid.major = element_line(colour = "#d9d9d9", size = 0.25),
        panel.grid.minor.x = element_line(colour = "#d9d9d9", size = 0.25)) 

#Fraction MN
MNfrac <- EvalData_FieldData %>% 
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells)) %>% 
  mutate(meanFraction_MNcells = mean(Fraction_MNcells)) %>%
  select(Cell.Type, med_Fraction_CAcells, meanFraction_MNcells) %>% 
  unique() %>% 
  ggplot(aes( y = reorder(Cell.Type, -med_Fraction_CAcells), x = meanFraction_MNcells)) +
  geom_bar(stat = "identity", fill = CMSblue, width = 0.75) +
  theme_classic(6) +
  labs(x = "% Non-mitotic cells \nwith micronuclei") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(colour = "#d9d9d9", size = 0.25),
        panel.grid.minor.x = element_line(colour = "#d9d9d9", size = 0.25)) 

#Mitotic Index
MI <- EvalData_FieldData %>%  
  group_by(Cell.Type) %>% 
  mutate(med_Fraction_CAcells = median(Fraction_CAcells)) %>% 
  select(Cell.Type, med_Fraction_CAcells) %>% 
  unique() %>% 
  left_join(mitoticIndex_cell, by = "Cell.Type") %>% 
  ggplot(aes( y = reorder(Cell.Type, -med_Fraction_CAcells), x = mitotic)) +
  geom_bar(stat = "identity", fill = CMSred, width = 0.75) +
  theme_classic(6) +
  labs(x = "Mitotic \nindex") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(colour = "#d9d9d9", size = 0.25),
        panel.grid.minor.x = element_line(colour = "#d9d9d9", size = 0.25))

# PCM fragmentation/Centriolar Centrosome
Cep164Frac <- Cep164 %>%
  gather(3:4, key = "fraction.type", value = "fraction") %>%
  mutate(fraction = fraction*100) %>%
  filter(fraction.type == "% Cep164 positive centrosomes") %>% 
  ggplot(aes(y = reorder(Cell.Type, -med_Fraction_CAcells), x = fraction)) +
  geom_bar(position = "stack", stat="identity", width = 0.75, fill = CMSgrey) +
  scale_x_continuous(breaks=seq(0,100,25)) +
  theme_classic(6) +
  labs(x = "% Cep164 \npositive") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(colour = "#d9d9d9", size = 0.25),
        panel.grid.minor.x = element_line(colour = "#d9d9d9", size = 0.25))

#Centrosome Size
CSize <- EvalData_FieldData %>% 
  group_by(Cell.Type) %>% 
  mutate(med_Csize = median(avCentrosome.Size, na.rm = T),
         mean_CSize = mean(avCentrosome.Size, na.rm = T),
         med_Fraction_CAcells = median(Fraction_CAcells)) %>%
  select(Cell.Type, mean_CSize, med_Csize, med_Fraction_CAcells) %>% 
  unique() %>% 
  ggplot(aes( y = reorder(Cell.Type, -med_Fraction_CAcells), x = mean_CSize)) +
  geom_bar(fill = CMSdarkgrey, width = 0.75, stat = "identity") +
  theme_classic(6) +
  labs(x = "Centrosome \nsize (µm²)") +
  theme(legend.position = "none",
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.grid.major = element_line(colour = "#d9d9d9", size = 0.25)) 

cowplot::plot_grid(NoNuclei, MI, CAfrac, CSize, Cep164Frac, MNfrac, nrow = 1, align = "h", rel_widths = c(1.1,0.5,2.5,0.5,0.5,0.75))
# ggsave("../figures/CAOV_paper_figure5.pdf", height = 8, width = 6)
# write.csv(EvalData_FieldData, "../source_data_figures/source_data_figure5.csv", row.names = F)


# FIGURE 6
# Fig6c - multi-cor plot
testing_data <- EvalData_FieldData %>% 
  select(Field_ID, Experiment, Cell.Type, AV_Well_NucleusSize, AV_Well_NucleusIntensity, avCentrosome.Size, Fraction_CAcells, Fraction_MNcells)

testing_data <- EvalData_FieldData %>% 
  group_by(Cell.Type) %>% 
  mutate(meanNucleiSize = mean(AV_Well_NucleusSize),
         meanNucleiIntensity = mean(AV_Well_NucleusIntensity),
         meanCentrosomeSize = mean(avCentrosome.Size)) %>% 
  select(Cell.Type, meanNucleiSize, meanNucleiIntensity, meanCentrosomeSize) %>% 
  unique() %>% 
  left_join(Combined_perCL[c("Cell.Type", "meanFraction_CAcells", "meanFraction_MNcells", "mitotic", "mean_fracgH2AXpos")], by = "Cell.Type") %>% 
  select(Cell.Type, 
         '% CA cells' = meanFraction_CAcells, 
         'Nuclei intensity' = meanNucleiIntensity,
         '% MN cells' = meanFraction_MNcells,
         'Mitotic index' = mitotic,
         '% ɣH2AX cells' = mean_fracgH2AXpos)

# function for scatter plots with smoothed trend line
lower_plots <- function(data, mapping, ...) {
  ggplot(data = data, mapping = mapping) +
    geom_point(color = CMSdarkgrey, size = 0.5, alpha = 0.8) +
    geom_smooth(method = "lm", size = 0.75,...) 
} 

# working with a trick of global assignment
diag_plots <- function(data, mapping, ...) {
  # increase counter each run globally so outside the function as well and this does the trick!
  x <<- x + 1
  ggplot(data = data, mapping = mapping) +
    # choose color by counter and send bin width argument in
    geom_histogram(fill = clrs[x], ...)
} 

# set the color and counter
clrs <- c(CMSlightblue, CMSteal, CMSblue, CMSred, CMSgreen, CMSlightblue, CMSteal, CMSblue, CMSred, CMSgreen, CMSlightblue, CMSteal, CMSblue, CMSred, CMSgreen)
x <- 0

ggpairs(testing_data,
        columns = 2:ncol(testing_data), 
        columnLabels = colnames(testing_data[2:ncol(testing_data)]),
        upper = list(continuous = wrap('cor',method = "spearman", size = 2.5, colour = CMSdarkgrey)),
        diag = list(continuous = wrap(diag_plots, colour = CMSdarkgrey, alpha = 0.9, size = 0.1)),
        lower = list(continuous = wrap(lower_plots, color=CMSdarkgrey, se=T))) +
  theme_bw(7) +
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
# ggsave("../figures/CAOV_paper_figure6c.pdf", height = 5, width = 5)
# write.csv(testing_data, "../source_data_figures/source_data_figure6c.csv", row.names = F)


# Fig6d-g
CellLine_meta_compare <- EvalData_FieldData %>% 
  group_by(Cell.Type) %>% 
  mutate(meanCAfra = mean(Fraction_CAcells),
         meanCAScore = mean(CAScore),
         meanMNfra = mean(Fraction_MNcells),
         meanCSize = mean(avCentrosome.Size)) %>% 
  select(Cell.Type, meanCAfra, meanCSize, meanMNfra) %>% 
  unique() %>% 
  left_join(meta_info, by = "Cell.Type")
CellLine_meta_compare$OxygenConditions <- factor(CellLine_meta_compare$OxygenConditions, levels=c("21%", "5%"))

## Suptype comparison
### Fraction CA
CellLine_meta_compare %>% 
  mutate(overallGroup = if_else(subtype == "HGSOC", "HGSOC", "other")) %>% 
  ggplot(aes(x = overallGroup, y = meanMNfra)) +
  geom_jitter(aes(colour = subtype2), width = 0.25, size = 0.75) +
  geom_boxplot(colour = CMSdarkgrey, notch = T, outlier.colour = NA, fill = NA) +
  scale_colour_manual(values = c("HGSOC"=CMSteal, "unknown"=CMSroyalblue, "normal"=CMSdarkgrey, 
                                 "LGSOC" = CMSorange, "endometrioid" = CMSpink, "mucinous" = CMSlightred), name = "subtype") +
  ggpubr::stat_compare_means(size = 2) +
  theme_bw(7) +
  labs(x = "", y = "% Non-mitotic cells \nwith micronuclei") +
  theme(legend.key.size = unit(0.1, "in"),
        legend.position = "none")+
  guides(colour=guide_legend(nrow=3))
# ggsave("../figures/CAOV_paper_figure6d.pdf", height = 2, width = 1.5)

### Fraction MN
CellLine_meta_compare %>% 
  mutate(overallGroup = if_else(subtype == "HGSOC", "HGSOC", "other")) %>% 
  ggplot(aes(x = overallGroup, y = meanMNfra)) +
  geom_jitter(aes(colour = subtype2), width = 0.25, size = 0.75) +
  geom_boxplot(colour = CMSdarkgrey, notch = T, outlier.colour = NA, fill = NA) +
  scale_colour_manual(values = c("HGSOC"=CMSteal, "unknown"=CMSroyalblue, "normal"=CMSdarkgrey, 
                                 "LGSOC" = CMSorange, "endometrioid" = CMSpink, "mucinous" = CMSlightred), name = "subtype") +
  ggpubr::stat_compare_means(size = 2) +
  theme_bw(7) +
  labs(x = "", y = "% Non-mitotic cells \nwith micronuclei") +
  theme(legend.key.size = unit(0.1, "in"),
        legend.position = "none")+
  guides(colour=guide_legend(nrow=3))
# ggsave("../figures/CAOV_paper_figure6e.pdf", height = 2, width = 2.25)

## Oxygen growth conditions
CellLine_meta_compare %>% 
  ggplot(aes(x = OxygenConditions, y = meanCAfra)) +
  geom_jitter(colour = CMSdarkgrey, width = 0.25, size = 0.75) +
  geom_boxplot(aes(colour = OxygenConditions), notch = T, outlier.colour = NA, fill = NA) +
  scale_colour_manual(values = c(CMSgrey, CMSpurple)) +
  ggpubr::stat_compare_means(size = 2) +
  theme_bw(7) +
  labs(x = "Oxygen growth conditions", y = "% Non-mitotic cells \nwith 2 or more centrosomes") +
  theme(legend.position = "none")
# ggsave("../figures/CAOV_paper_figure6f.pdf", height = 2, width = 1.5)

CellLine_meta_compare %>% 
  ggplot(aes(x = OxygenConditions, y = meanMNfra)) +
  geom_jitter(colour = CMSdarkgrey, width = 0.25, size = 0.75) +
  geom_boxplot(aes(colour = OxygenConditions), notch = T, outlier.colour = NA, fill = NA) +
  scale_colour_manual(values = c(CMSgrey, CMSpurple)) +
  ggpubr::stat_compare_means(size = 2) +
  theme_bw(7) +
  labs(x = "Oxygen growth conditions", y = "% Non-mitotic cells \nwith micronuclei") +
  theme(legend.position = "none")
# ggsave("../figures/CAOV_paper_figure6g.pdf", height = 2, width = 1.5)
# write.csv(CellLine_meta_compare, "../source_data_figures/source_data_figure6d-g.csv")



# Fig6h-j
Normoxia <- read.table('../data/cell_line_imaging_data/oxygen_growth_conds_experiment/Normal_oxygen_raw_results.txt',  skip = 9, header = TRUE, sep = "\t")
Hypoxia <- read.table('../data/cell_line_imaging_data/oxygen_growth_conds_experiment/Low_oxygen_raw_results.txt',  skip = 9, header = TRUE, sep = "\t")  

Joined_all <- bind_rows(Normoxia, Hypoxia)
colnames(Joined_all)

Norm_vs_Hypoxia <- Joined_all %>% 
  select(Row, Column, Field, Object.No, Compound, Concentration, Cell.Type, 
         Valid.Cells...Object.No.in.Nuclei, Valid.Cells...Nucleus.Area..µm.., 
         Valid.Cells...Intensity.Nuclei.Mean, Valid.Cells...Nucleus.Roundness, 
         Valid.Cells...Intensity.Pax8.Mean, Valid.Cells...Number.of.Micronuclei..per.Cell, 
         Valid.Cells...Valid.with.MN, Valid.Cells...Total.Spot.Area..Mean.per.Cell, 
         Valid.Cells...Number.of.Centrosomes..per.Cell...Mean.per.Cell, 
         Valid.Cells...Valid.cell.with.2.or.more.Centrosomes, Valid.Cells...Valid.cells.with.MN.and.CA) %>% 
  unite("Well_ID", 1:2, sep = ".")

PlateLayout_both <- Norm_vs_Hypoxia %>% 
  select(Well_ID, Cell.Type) %>% 
  unique()

NormVsHypoxia_WellData <- Norm_vs_Hypoxia %>% 
  group_by(Well_ID, Compound) %>% 
  add_count(Well_ID, Compound) %>% 
  rename(Nuclei.Number = n) %>% 
  group_by(Well_ID, Compound) %>% 
  mutate(AV_Well_NucleusSize = mean(Valid.Cells...Nucleus.Area..µm..),
         AV_Well_NucleusIntensity = mean(Valid.Cells...Intensity.Nuclei.Mean),
         AV_Well_NucleusRoundness = mean(Valid.Cells...Nucleus.Roundness),
         AV_Well_Pax8Intensity = mean(Valid.Cells...Intensity.Pax8.Mean),
         Centrosome.Number = sum(Valid.Cells...Number.of.Centrosomes..per.Cell...Mean.per.Cell, na.rm=TRUE),
         AV_Well_Centrosome.Size = mean((Valid.Cells...Total.Spot.Area..Mean.per.Cell/Valid.Cells...Number.of.Centrosomes..per.Cell...Mean.per.Cell), na.rm = TRUE),
         MN.Number = sum(Valid.Cells...Number.of.Micronuclei..per.Cell, na.rm=TRUE),
         CAScore = (Centrosome.Number/Nuclei.Number),
         MNScore = (MN.Number/Nuclei.Number),
         MN_posCells = sum(Valid.Cells...Valid.with.MN),
         CentrosomeAmp_Cells = sum(Valid.Cells...Valid.cell.with.2.or.more.Centrosomes, na.rm=TRUE),
         Fraction_MNcells = (MN_posCells/Nuclei.Number)*100,
         Fraction_CAcells = (CentrosomeAmp_Cells/Nuclei.Number)*100) %>% 
  select(Well_ID, Compound, Concentration, Cell.Type, Nuclei.Number, AV_Well_NucleusSize, AV_Well_NucleusIntensity, AV_Well_NucleusRoundness, AV_Well_Pax8Intensity, Centrosome.Number, AV_Well_Centrosome.Size, CAScore, CentrosomeAmp_Cells, Fraction_CAcells, MN.Number, MNScore, MN_posCells, Fraction_MNcells) %>% 
  unique()

NormVsHypoxia_WellData$Compound <- factor(NormVsHypoxia_WellData$Compound,
                                          levels = c('Normoxia','Hypoxia'),ordered = TRUE)

#Data per Cell Line
NormVsHypoxia_CellLineData <- NormVsHypoxia_WellData %>% 
  group_by(Cell.Type, Compound) %>% 
  mutate(CellLine_Fraction_CAcells = mean(Fraction_CAcells),
         CellLine_Fraction_MNcells = mean(Fraction_MNcells),
         CellLine_CentrosomeSize = mean(AV_Well_Centrosome.Size),
         CellLine_NucleiSize = mean(AV_Well_NucleusSize),
         CellLine_NucleiIntensity = mean(AV_Well_NucleusIntensity)) %>% 
  select(Cell.Type, Compound, CellLine_Fraction_CAcells, CellLine_Fraction_MNcells, CellLine_CentrosomeSize, CellLine_NucleiSize, CellLine_NucleiIntensity) %>%
  unique()

## CA fraction vs growth conditions
NormVsHypoxia_CellLineData %>% 
  group_by(Cell.Type, Compound) %>% 
  mutate(Compound = case_when(Compound == "Hypoxia" ~ "5%",
                              Compound == "Normoxia" ~ "21%")) %>% 
  ggplot(aes( x = Compound, y = CellLine_Fraction_CAcells)) +
  geom_boxplot(colour = CMSdarkgrey, fill = NA, notch = T, outlier.colour = NA) +
  geom_point(aes(colour = Cell.Type), size = 0.75) +
  geom_line(aes(group = Cell.Type, colour = Cell.Type), size = 0.5) +
  scale_colour_manual(values =  CMS_pal, name = "cell line") +
  theme_bw(7) +
  labs(x = "Oxygen growth conditions",
       y = "% Cells with more than 2 centrosomes") +
  theme(legend.position = "right",
        legend.key.size = unit(0.1, "in")) + 
  stat_compare_means(paired = TRUE, size = 2)
# ggsave("../figures/CAOV_paper_figure6h.pdf", height = 2, width = 2)

## MN fraction vs growth conditions
NormVsHypoxia_CellLineData %>% 
  group_by(Cell.Type, Compound) %>% 
  mutate(Compound = case_when(Compound == "Hypoxia" ~ "5%",
                              Compound == "Normoxia" ~ "21%")) %>% 
  ggplot(aes( x = Compound, y = CellLine_Fraction_MNcells)) +
  geom_boxplot(colour = CMSdarkgrey, fill = NA, notch = T, outlier.colour = NA) +
  geom_point(aes(colour = Cell.Type), size = 0.75) +
  geom_line(aes(group = Cell.Type, colour = Cell.Type), size = 0.5) +
  scale_colour_manual(values =  CMS_pal, name = "cell line") +
  theme_bw(7) +
  labs(x = "Oxygen growth conditions",
       y = "% Cells with micronuclei") +
  theme(legend.position = "right",
        legend.key.size = unit(0.1, "in")) + 
  stat_compare_means(paired = TRUE, size = 2)
# ggsave("../figures/CAOV_paper_figure6i.pdf", height = 2, width = 2)
# write.csv(NormVsHypoxia_CellLineData, "../source_data_figures/source_data_figure6h-i.csv")

