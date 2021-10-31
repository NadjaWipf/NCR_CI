####################################################################################################################################################
## This R code is to replicate the analyses, figures and tables from the research article entitled:
## "Multi-insecticide resistant malaria vectors in the field remain susceptible to malathion, despite the presence of Ace1 point mutations", 2021
## Analyses: DNA genotyping of mosquito legs, calculate mutant allele frequencies
## Figures: Fig 6, S4 Fig
## Code written by NW
####################################################################################################################################################

# --------------------------------------------------------------------------------------------------------------
# load required packages
library(cowplot)
library(ggpubr)
library(tidyverse)

# set working directory to a "data" folder where you have saved all relevant STables "input_data" as .txt files.
# setwd("~/data")

# --------------------------------------------------------------------------------------------------------------

# load data frame
legs <- read.table("S6Table.txt", header = T, sep= "\t")

# reorder and rename variables
colnames(legs)
legs <- legs%>%
  select(ID, Site, Insecticide, RNAseq, qPCR_Vgsc_L995ForS, qPCR_Ace1_G280S)%>%
  # create combined variable for the Vgsc-L995F/S and Ace1-G280S locus
  tidyr::unite("Vgsc_Ace1", qPCR_Vgsc_L995ForS:qPCR_Ace1_G280S, remove=F)%>%
  mutate_if(is.character, as.factor)

levels(legs$Site)
legs$Site <- factor(legs$Site, levels = c("Agboville", "Dabou", "Tiassale", "Mali-NIH", "Ngousso"),
                    labels = c("Agboville", "Dabou", "Tiassalé", "Mali-NIH", "Ngousso"))

levels(legs$Insecticide)
legs$Insecticide <- factor(legs$Insecticide, levels = c("Control", "Deltamethrin", "Malathion"),
                           labels = c("unexposed controls (C)", "6.4% deltamethrin survivors (D)", "2.5% malathion survivors (M)"))

levels(legs$qPCR_Vgsc_L995ForS)
legs$qPCR_Vgsc_L995ForS <- fct_relevel(legs$qPCR_Vgsc_L995ForS, c("FF", "LF", "LL"))

levels(legs$qPCR_Ace1_G280S)
legs$qPCR_Ace1_G280S <- fct_relevel(legs$qPCR_Ace1_G280S, c("SS", "GS", "GG"))

levels(legs$Vgsc_Ace1)
legs$Vgsc_Ace1 <- fct_relevel(legs$Vgsc_Ace1, c("FF_SS", "FF_GS", "LF_SS", "FF_GG", "LF_GS", "LF_GG", "LL_GS", "LL_GG"))

# --------------------------------------------------------------------------------------------------------------

## Fig 6. (A) Vgsc-L995F and (B) Ace1-G280S mutation in 40 insecticide insecticide-unexposed individuals per field population

# subset the 120 specimens
legs_120 <- legs%>%
  filter(Site %in% c("Agboville", "Dabou", "Tiassalé") & Insecticide == "unexposed controls (C)")%>%
  droplevels()


#################################################################
###   calculate allele frequency for 120 control specimens   ####
#################################################################

# calculate kdr = Vgsc-995F allele frequency per Site + experimental condition

kdr1 <- legs_120 %>%
  group_by(Site, Insecticide, qPCR_Vgsc_L995ForS) %>%
  summarise(n=n()) %>%
  droplevels()

kdr2 <- spread(kdr1, qPCR_Vgsc_L995ForS, n, fill = 0) # spreads LL/LF/FF into 3 columns, na = 0
kdr2$kdr_total <- (2*(kdr2$LL + kdr2$LF + kdr2$FF)) # calc all alleles per obs
kdr2$F_count <- ((kdr2$LF) + (2*kdr2$FF)) # calc F alleles per obs
kdr2$"kdrF_freq" <- (kdr2$F_count/kdr2$kdr_total*100) # calc F allele freq per obs

## calculate Ace1 = Ace1-280S allele frequency per Site + experimental condition
ace1 <- legs_120 %>%
  group_by(Site, Insecticide, qPCR_Ace1_G280S) %>%
  summarise(n=n())  %>%
  droplevels()

ace2 <- spread(ace1, qPCR_Ace1_G280S, n, fill = 0)
ace2$ace_total <- (2*(ace2$GG + ace2$GS + ace2$SS))
ace2$S_count <- ((ace2$GS) + (2*ace2$SS))
ace2$"aceS_freq" <- (ace2$S_count/ace2$ace_total*100)

# join the 2 tables
freq_table <- full_join(kdr2, ace2, by = c("Site", "Insecticide"))

## add super-kdr = Vgsc-1570Y frequency (0%) - 100% NN
freq_table <- add_column(freq_table, "kdrS_1570Y_freq" = 0)

# round to full numbers, remove unnecessary columns
freq <- freq_table %>%
  select(-c(F_count, S_count)) %>%
  mutate_if (is.numeric, round, 0)

# optional: to save table un-comment next line
# write.table(freq, file="Allele_Freq_120controls.txt", row.names=F, dec=".", sep="\t")

# --------------------------------------------------------------------------------------------------------------
#####################
###   Figure 6   ####
#####################

# Fig 6A - BAR PLOT for kdr Vgsc_L995F of 120 unexposed controls
kdr_bar <- ggplot(data=legs_120, aes(Site))+
  theme_bw()+
  facet_grid(~ Insecticide)+
  geom_bar(aes(fill=qPCR_Vgsc_L995ForS), colour="white", size=0.1)+
  geom_text(stat='count', aes(label=qPCR_Vgsc_L995ForS, fill=qPCR_Vgsc_L995ForS), position = position_stack(vjust = 0.5))+
  scale_y_continuous(breaks = seq(0, 40, 5), limits = c(-3, 41), expand = c(0, 0))+
  scale_fill_manual(values = c("#238b45", "#66c2a4", "#ccece6"),
                    name = "Vgsc-L995F alleles\nbelow bar:\nmutant allele frequency (%)",
                    labels= c("FF = homozygous mutant", "LF = heterozygous", "LL = homozygous wildtype"))+
  labs(y = "Number of mosquitoes\n")+
  theme(axis.title.y = element_text(size=12, colour="black", face= "bold"),
        axis.title.x = element_blank(),
        axis.text  = element_text(size=12, colour="black"),
        strip.text.x = element_text(colour="black", size=12, face="plain"))

Fig6A <- kdr_bar +
  geom_text(data = freq, aes(label = paste(kdrF_freq, "%", sep=""), fontface="bold"), y=-1.3)

# --------------------------------------------------------------------------------------------------------------

# Fig 6B - BAR PLOT for Ace1_G280S of 120 unexposed controls
ace_bar <- ggplot(data=legs_120, aes(Site))+
  theme_bw()+
  facet_grid(~ Insecticide)+
  geom_bar(aes(fill=qPCR_Ace1_G280S), colour="white", size=0.1)+
  geom_text(stat='count', aes(label=qPCR_Ace1_G280S, fill=qPCR_Ace1_G280S), position = position_stack(vjust = 0.5))+
  scale_y_continuous(breaks = seq(0, 40, 5), limits = c(-3, 41), expand = c(0, 0))+
  scale_fill_manual(values = c("#0570b0", "#74a9cf", "#d0d1e6"),
                    name = "Ace1-G280S alleles\nbelow bar:\nmutant allele frequency (%)",
                    labels= c("SS = homozygous mutant", "GS = heterozygous", "GG = homozygous wildtype"))+
  labs(y = "Number of mosquitoes\n")+
  theme(axis.title.y = element_text(size=12, colour="black", face ="bold"),
        axis.title.x = element_blank(),
        axis.text  = element_text(size=12, colour="black"),
        strip.text.x = element_text(colour="black", size=12, face="plain"))

Fig6B <- ace_bar +
  geom_text(data = freq, aes(label = paste(aceS_freq, "%", sep=""), fontface="bold"), y=-1.3)

# arrange Fig6A and Fig6B next to eachother 
ggpubr::ggarrange(Fig6A, Fig6B,
          labels = c("A", "B"),
          ncol = 2)
# --------------------------------------------------------------------------------------------------------------

# subset RNA-sequenced specimens
Rlegs <- legs%>%filter(RNAseq == 1)%>%droplevels()

# Fige 6C - BAR PLOT showing in which combination the two resistance-associated alleles occurred in the 55 RNA -sequenced individuals
Fig6C <- ggplot(data=Rlegs, aes(Site, fill=Vgsc_Ace1))+
  theme_bw()+
  geom_bar(colour= "white")+
  facet_grid(~ Insecticide)+
  geom_text(stat='count', aes(label= Vgsc_Ace1), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a8ddb5", "#ccebc5", "#e0f3db"),
                    name = "Alleles of Vgsc & Ace1")+
  labs(
    #title = "Combined Vgsc & Ace1 genotype of 55 RNA-sequenced individuals", 
    y = "Number of mosquitoes\n", x = "\nMosquito population")+   # label x and y axix, \n gives more space btw label and axis
  theme(#plot.title = element_text(colour="black", size=16, face="plain"), # change text size & colour
    axis.text = element_text(colour="black", size=10, face="plain"), 
    axis.title = element_text(colour="black", size=12, face="bold"),
    strip.text.x = element_text(colour="black", size=12, face="plain"))

# --------------------------------------------------------------------------------------------------------------

# Figure 6 => arrange the three plots in one final figure with three panels

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()

sep2 <-ggarrange(Fig6A, blankPlot, Fig6B,
                 ncol = 3, labels = c("A",  NA, "B"),
                 widths = c(4.4, 0.3, 4.4)
                 #, common.legend = TRUE, legend = "right"
                 )

Fig6.final <- ggarrange(sep2, blankPlot, Fig6C,
                   ncol = 1, nrow = 3,  labels = c( NA , NA, "C"),
                   heights = c(6, 0.3, 5.5))


# ggsave("Fig6.tiff", plot= Fig6.final, device = "tiff",
       # scale =1, width = 30, height = 22, units = "cm", dpi = 300)


# --------------------------------------------------------------------------------------------------------------

# S4 Fig. Bar plot showing combined Vgsc-L995F and Ace1-G280S genotypes of Ivorian field populations. 
S4Fig.legs_120 <- ggplot(data=legs_120, aes(Site, fill=Vgsc_Ace1))+
  theme_bw()+
  geom_bar(colour= "white")+
  facet_grid(~ Insecticide)+
  geom_text(stat='count', aes(label= Vgsc_Ace1), position = position_stack(vjust = 0.5))+
  scale_fill_manual(values = c("#014636", "#016c59", "#02818a", "#3690c0", "#67a9cf", "#a8ddb5", "#ccebc5", "#e0f3db"),
                    name = "Alleles of Vgsc & Ace1")+
  labs(title = "Combined Vgsc & Ace1 genotype of field populations",   # give a title
       y = "Number of mosquitoes\n", x = "\nMosquito population")+   # label x and y axix, \n gives more space btw label and axis
  theme(plot.title = element_text(colour="black", size=16, face="plain"), # change text size & colour
        axis.title = element_text(size=12, colour="black", face ="bold"),
        axis.text  = element_text(size=12, colour="black"),
        strip.text.x = element_text(colour="black", size=12, face="plain"))

# ggsave("S4Fig.tiff", plot = S4Fig.legs_120, device = "tiff", scale =1, width = 25, height = 21, units = "cm", dpi = 300)

# --------------------------------------------------------------------------------------------------------------

# optional.
#######################################################
###  calculate allele frequency of RNAseq samples  ####
#######################################################

colnames(Rlegs)

## calculate mutant allele frequency
# calculate kdr = Vgsc-995F allele frequency per Site + experimental condition
kdr1 <- Rlegs %>%
  group_by(Site, Insecticide, qPCR_Vgsc_L995ForS) %>%
  summarise(n=n()) %>%
  droplevels()

kdr2 <- spread(kdr1, qPCR_Vgsc_L995ForS, n, fill = 0) # spreads LL/LF/FF into 3 columns, na = 0
kdr2$kdr_total <- (2*(kdr2$LL + kdr2$LF + kdr2$FF)) # calc all alleles per obs
kdr2$F_count <- ((kdr2$LF) + (2*kdr2$FF)) # calc F alleles per obs
kdr2$"kdrF_freq" <- (kdr2$F_count/kdr2$kdr_total*100) # calc F allele freq per obs

## calculate Ace1 = Ace1-280S allele frequency per Site + experimental condition
ace1 <- Rlegs %>%
  group_by(Site, Insecticide, qPCR_Ace1_G280S) %>%
  summarise(n=n())  %>%
  droplevels()

ace2 <- spread(ace1, qPCR_Ace1_G280S, n, fill = 0)
ace2$ace_total <- (2*(ace2$GG + ace2$GS + ace2$SS))
ace2$S_count <- ((ace2$GS) + (2*ace2$SS))
ace2$"aceS_freq" <- (ace2$S_count/ace2$ace_total*100)

# join the 2 tables
freq_table <- full_join(kdr2, ace2, by = c("Site", "Insecticide"))

## add super-kdr = Vgsc-1570Y frequency (0%) - 100% NN
freq_table <- add_column(freq_table, "kdrS_1570Y_freq" = 0)

# round to full numbers, remove unnecessary columns
freq <- freq_table %>%
  select(-c(F_count, S_count)) %>%
  mutate_if (is.numeric, round, 0)

# save it as table if you like (un-comment next line)
# write.table(freq, file="Allele_Freq_55RNAseq.txt", row.names=F, dec=".", sep="\t")

# --------------------------------------------------------------------------------------------------------------
# End of R script.