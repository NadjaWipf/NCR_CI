####################################################################################################################################################
## This R code is to replicate the analyses, figures and tables from the research article entitled:
## "Multi-insecticide resistant malaria vectors in the field remain susceptible to malathion, despite the presence of Ace1 point mutations", 2021
## Analyses: DNA genotyping of mosquito legs, calculate mutant allele frequencies
## Figures: Fig 3, Fig 4, Fig 5, S3 Fig
## Code written by NW
####################################################################################################################################################

# --------------------------------------------------------------------------------------------------------------
# load required packages
library(cowplot)
library(ggpubr)
library(ggrepel)
library(lemon)
library(tidyverse)
library(scales)

# set working directory to a "data" folder where you have saved all relevant STables "input_data" as .txt files.
# setwd("~/data")

# --------------------------------------------------------------------------------------------------------------

# load data frame
fieldvslab <- read.table("S4Table.txt", header = T, sep= "\t", fill=TRUE, quote="")

colnames(fieldvslab)
levels(as.factor(fieldvslab$treatment_comparison))

# subset comparisons with insecticide-unexposed controls
fl.c <- fieldvslab%>%
  filter(treatment_comparison == "Field_C_vs.Lab2_C")%>%droplevels()%>%
  mutate_if(is.character, as.factor)

# new name for MS figures
levels(fl.c$comparison_full)
fl.c$comparison_full <- str_replace_all(fl.c$comparison_full,
                                        c("Agb_C_vs.Lab2_C"="Agb_C vs. Lab2_C", "Dab_C_vs.Lab2_C"="Dab_C vs. Lab2_C", "Tia_C_vs.Lab2_C"="Tia_C vs. Lab2_C"))

# rename some column names
ms <- fl.c%>%
  dplyr::rename(protein_family = description, comparison = comparison_full, stelle = significance_level, logFC = log2FC)%>%
  mutate_if(is.character, as.factor)

# --------------------------------------------------------------------------------------------------------------

##################################
###    Fig 3. Volcano plots    ###
##################################

### create a variable with the protein_families to highlight in the volcano plots
levels(ms$protein_family)

ms.v <- ms%>%
  mutate(IR_volcano_small = case_when(str_detect(protein_family,"ABC_transporter") ~ "ABC transporter",
                                      str_detect(protein_family,"COE")  ~ "COE",
                                      str_detect(gene_name,"ACE1")  ~ "COE",
                                      str_detect(protein_family,"^cuticular_protein") ~ "cuticular protein",
                                      str_detect(protein_family,"GST")  ~ "GST",
                                      str_detect(protein_family,"UGT")  ~ "UGT",
                                      str_detect(protein_family,"^P450")  ~ "P450",
                                      str_detect(protein_family,"ATP_synthase$")  ~ "H+-transporting ATP synthase"
                                      ))

# rename protein_family_detail names to have same variables 
anno.fl <- ms.v%>%
  dplyr::rename(protein_family_detail = protein_family,
                protein_family_with_NA = IR_volcano_small)%>%
  mutate_if(is.character, as.factor)

levels(anno.fl$protein_family_detail)
levels(anno.fl$protein_family_with_NA)

# How many gene_names per protein_family for the volcano plot?
n <- anno.fl %>%
  filter(comparison == "Agb_C vs. Lab2_C")%>%
  #group_by(protein_family_detail, protein_family_with_NA)%>%
  group_by(protein_family_with_NA)%>%
  summarize(n =n())

# define the gene_names to lable in the volcano plot
gene_labels <- c("GSTD3", "GSTD12", "GSTMS2", "COEAE5G", "COEAE6G", "ACE1",
                 "CYP4H17", "CYP6AA1", "CYP6AA2", "CYP6AG1", "CYP6AG2", "CYP6AK1", "CYP6M2", "CYP6N1", "CYP6P1", "CYP6P2", "CYP6P3", "CYP6P4", "CYP6P5", "CYP6Z1", "CYP6Z4", "CYP9K1", "CYP9L1",
                 "COEAE6O", "COEAE8O","CPLCP25","LYSC5", "TRYP6", "TRYP7",
                 "COEAE3H","COEAE6G","CPLCP25",
                 "CYP303A1", "CLIPA5", "HPX10", "TEP11",
                 "CYP12F3", "TSF1")

# duplicate protein_family variable and name NAs as "other" - otherwise NA points will not be drawn in volcano plot
anno.fl <- anno.fl%>%
  mutate(protein_family = protein_family_with_NA)%>%
  mutate(IR_enz_names = ifelse(gene_name %in% gene_labels, as.character(gene_name), NA))%>%  # label selected  genes
  mutate_if(is.character, as.factor)

anno.fl$protein_family <- anno.fl$protein_family %>% fct_explicit_na(na_level = "other")

## Fig 3. create function for volcano plot with labels
volcano_fct_FvsL.lab <- function(df.all, comp , pval, fc_up, fc_down) {
  
  # filter run_comparison_no
  mRNA <- df.all %>%
    filter(comparison %in% comp)%>%droplevels()
  
  col_vector <- c("orchid", "blue", "gold", "deepskyblue", "slategray", "red", "limegreen", "lightgray")
  myColors <- col_vector
  names(myColors) <- levels(mRNA$protein_family)
  
  ## volcano plot
  volcano_plot <- ggplot(mRNA, aes(x = logFC, y = -log10(FDR))) +
    geom_point(aes(color = factor(protein_family)), alpha= 0.8) +
    geom_point(aes(color = factor(protein_family_with_NA))) +
    scale_color_manual(name = "Protein family", values = myColors) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom") +
    #geom_vline(xintercept = c(-1, 1), linetype=3) +
    geom_vline(xintercept = c(0), linetype=3) +
    geom_hline(yintercept = -log10(0.01), linetype = 3) +
    labs(x = expression(log[2](FC)), y = expression(-log[10](FDR)), 14)+
    #facet_grid(. ~ comparison)+
    #facet_wrap( ~ comparison, nrow = 2, ncol=2)+
    facet_rep_wrap(~ comparison, ncol=2, repeat.tick.labels= T)+
    theme(axis.title = element_text(colour="#000000"),
          axis.text  = element_text(colour="#000000"),
          strip.text = element_text(size = 12, colour="#000000"))+
    
    # label gene_names
    geom_text_repel(
      data = subset(mRNA, logFC <= fc_down),
      aes(label = ifelse(FDR <= pval, str_c(IR_enz_names), "")),
      nudge_x = -8,
      size= 3.5,
      segment.alpha = 0.5,
      direction = "y",
      hjust = 1
    ) +
    geom_text_repel(
      data = subset(mRNA, logFC >= fc_up),
      aes(label = ifelse(FDR <= pval, str_c(IR_enz_names), "")),
      nudge_x = 8,
      size = 3.5,
      segment.alpha = 0.5,
      direction = "y",
      hjust = 0)
  
  return(volcano_plot)
}

anno.fl$comparison <- fct_relevel(anno.fl$comparison,c("Agb_C vs. Lab2_C", "Tia_C vs. Lab2_C", "Dab_C vs. Lab2_C"))

Fig3 <- volcano_fct_FvsL.lab(anno.fl, comp = c("Agb_C vs. Lab2_C", "Dab_C vs. Lab2_C", "Tia_C vs. Lab2_C"),
                                    pval = 0.01, fc_up = 0, fc_down = 0)

# un-comment to save / export figure
# ggsave("Fig3_volcano_lab.tiff", plot = Fig3, device = "tiff", scale = 1, width = 30, height  = 35, units = "cm", dpi = 300)
# ggsave("Fig3_volcano_lab.pdf", plot = Fig3, device = "pdf", scale = 1, width = 31, height  = 35, units = "cm", dpi = 300)

# --------------------------------------------------------------------------------------------------------------

#############################
###   Fig 4. HEAT MAPs    ###
#############################

## create heat plot function
    # no title / no axis titles / no font size => for figures in manuscript

HEAT.map.fct.m <- function(df) {
  
  # set font to bold only gene_name, but keep it regular for gene_id
  # bold.labels <- levels(df$gene_name)
  
  ## heat plot
  heat_map <- ggplot(df, aes(x=comparison, y= gene_name:gene_id, fill=round(logFC, digits=1))) +
    theme_minimal()+        # plain background
    geom_tile() +           # creates the heat map
    geom_text(aes(label = paste(format(round(logFC, digits=1), nsmall = 1), stelle))) +  # round FC to one digit after comma (keep 0) + add sign.*-***
    scale_fill_gradient2(low = muted("blue"), mid = "white", high = muted("red"), midpoint = 0, limits = c(-7.6,6.5), name = expression(log[2](FC)))+
    theme(axis.title = element_blank(),
          # axis.text.y = element_text(colour="black", face = bold.labels), 
          axis.text.y = element_text(colour="black"), 
          axis.text.x = element_text(colour="black", face="bold"),
          legend.position = "bottom")
  
  return(heat_map)
  
}

levels(ms$protein_family)

# Fig.4A : heat map for Ace-1 & COEs
esti.sig <- ms%>%
  #filter(logFC > 0)%>%
  filter(gene_name == "ACE1" | protein_family == "COE")%>%
  filter(stelle != "ns")%>%droplevels()

esti.a <- esti.sig%>%
  group_by(gene_id)%>%
  filter(stelle != "*")%>%
  droplevels()

gene.id.a <- levels(esti.a$gene_id)

esti.aa <- esti.sig%>%
  filter(gene_id %in% c(gene.id.a))%>%droplevels()

hm.COE.a <- HEAT.map.fct.m(esti.aa)+theme(legend.position = "none")
hm.COE.a

# Fig.4B: heat map for GSTs
gst.sig <- ms%>%
  filter(protein_family == "GST")%>%
  filter(stelle != "ns")%>%
  droplevels()

gst.b <- gst.sig%>%
  group_by(gene_id)%>%
  filter(stelle != "*")%>%
  droplevels()

gene.id.b <- levels(gst.b$gene_id)

gst.bb <- gst.sig%>%
  filter(gene_id %in% c(gene.id.b))%>%droplevels()

hm.GST.b <- HEAT.map.fct.m(gst.bb)
hm.GST.b <- HEAT.map.fct.m(gst.bb)+theme(legend.position = "none")
hm.GST.b

# Fig.4C: heat map for P450s
cyp.sig <- ms%>%
  filter(protein_family == "P450")%>%
  filter(stelle != "ns")%>%
  #filter(logFC > 0)%>%  # if only UP regulated genes shown
  droplevels()

cyp.c <- cyp.sig%>%
  group_by(gene_id)%>%
  filter(stelle != "*")%>%
  droplevels()

gene.id.c <- levels(cyp.c$gene_id)

cyp.cc <- cyp.sig%>%
  filter(gene_id %in% c(gene.id.c))%>%droplevels()

nlevels(cyp.sig$gene_name)
hm.p450.c <- HEAT.map.fct.m(cyp.cc)
hm.p450.c

# Fig.4. one legend for all - to have correct colour range
bind_rows_keep_factors <- function(...) {
  ## Identify all factors
  factors <- unique(unlist(
    map(list(...), ~ select_if(..., is.factor) %>% names())
  ))
  ## Bind dataframes, convert characters back to factors
  suppressWarnings(bind_rows(...)) %>% 
    mutate_at(vars(one_of(factors)), factor)
  #Source: https://stackoverflow.com/questions/42278020/bind-rows-of-data-frames-with-some-factor-columns
}

all.3 <- bind_rows_keep_factors(esti.aa, gst.bb, cyp.cc)

hm.all <- HEAT.map.fct.m(all.3)

# Fig4. arrange the three plots in one

blankPlot <- ggplot()+geom_blank(aes(1,1)) + 
  cowplot::theme_nothing()
legend.cyp <- get_legend(hm.p450.c)
legend.all <- get_legend(hm.all)

ab <- ggpubr::ggarrange(blankPlot, hm.COE.a, blankPlot, hm.GST.b, legend.all,
                        nrow = 5, labels = c(NA, "A",  NA, "B", NA), vjust = 0.3, heights = c(0.2, 3, 0.3, 2.2, 0.8)
)

c <- ggpubr::ggarrange(blankPlot, hm.p450.c,
                       nrow = 2,  labels = c( NA , "C"), vjust = 0.3, heights = c(0.2, 6.3),
                       legend= "none"
) 

Fig4 <- ggpubr::ggarrange(ab, blankPlot, c,
                         ncol = 3, nrow = 1, widths = c(3.8, 0.2, 4)
) 

# un-comment to save / export figure
# Fig4 + ggsave("Fig4.tiff", device = "tiff", scale = 1, width = 30, height = 22, units = "cm", dpi = 300)
# Fig4 + ggsave("Fig4.pdf", device = "pdf", scale = 1, width = 30, height = 22, units = "cm", dpi = 300)


# --------------------------------------------------------------------------------------------------------------

######################################################
###     Fig 5. lm detox RT-qPCR vs. RNA-seq       ####
######################################################

# pull gene_id  s for transcripts for which we have detox qPCR data

detox_seq <- ms%>%
  filter(gene_name %in% c("CYP6M2","CYP6P1","CYP6P3","CYP6P4","CYP6Z1","CYP9K1"))%>%droplevels()%>%
  rename(logFC_seq = logFC)%>%
  mutate(seq_FC = 2^(logFC_seq))%>%
  select(comparison, population_comparison, gene_id, gene_name, logFC_seq, seq_FC, FDR)
  
  
detox_pcr.1 <- read.table("S5Table.txt", sep="\t", header=T)
detox_pcr <- detox_pcr.1%>%
  filter(gene_name %in% c("CYP6M2","CYP6P1","CYP6P3","CYP6P4","CYP6Z1","CYP9K1"))%>%droplevels()%>%
  mutate(logFC_pcr = log2(qPCR_FC))%>%
  select(comparison, gene_name, logFC_pcr, qPCR_FC, P_H1)%>%
  mutate_if(is.character, as.factor)
  
detox_pcr$comparison <- str_replace_all(detox_pcr$comparison, c("_Lab2"="_C vs. Lab2_C"))

detox.2 <- full_join(detox_seq, detox_pcr, by = c("comparison", "gene_name"))

# separated by site => shows that slopes are similar enough to pool
ggplot(detox.2, aes(x = logFC_pcr, y = logFC_seq, colour = comparison))+
  geom_point()+
  geom_text(aes(label = gene_name), hjust=0, vjust=0, nudge_x = .05)+
  labs(x = "RT-qPCR log2 FC", y = "RNA-seq log2 FC")+
  stat_smooth(method = "lm")+
  theme_bw()

# correlation test, r = 0.9609266 (ADT)
cor.test(formula = ~ logFC_pcr + logFC_seq,
         data = detox.2)

# linear regression model
fit <- lm(logFC_seq ~ logFC_pcr , data = detox.2)
summary(fit)
summary(fit)$r.squared 
# r.squared ADT = 0.92338

# nice one with equation + R2
# https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph

lm_eqn <- function(detox.2){
  m <- lm(logFC_seq ~ logFC_pcr, detox.2);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2,
                   # italic(y) == a + b %.% italic(x)*","~~italic(R)^2~"="~r2*","~italic(p)~"="~pval, 
                   list(a = format(unname(coef(m)[1]), digits = 2),
                        b = format(unname(coef(m)[2]), digits = 2),
                        #pval = format(summary(m)$coefficients[2,4], digits = 2), # put r + p-value OR R2 + equ
                        r2 = format(summary(m)$r.squared, digits = 3)))
  as.character(as.expression(eq));
}

# Fig 5.Linear relationship in fold changes measured in P450s between RNA-seq and RT-qPCR. 
Fig5 <- ggplot(detox.2, aes(x = logFC_pcr, y = logFC_seq))+
  stat_smooth(method = "lm", colour = "black", size = 0.5)+
  geom_point(aes(colour = gene_name, shape = comparison), size= 2)+
  labs(x = expression(`RT-qPCR`~log[2](FC)), y = expression(`RNA-seq`~log[2](FC)), size = 14,
       colour = "gene name")+
  scale_shape_discrete(labels = c("Agb_C vs. Lab2_C", "Dab_C vs. Lab2_C", "Tia_C vs. Lab2_C"))+
  xlim(0.8,5.7)+ylim(0.8,5.7)+
  geom_text(x = 2.2, y = 5.2, label = lm_eqn(detox.2), parse = TRUE)+
  theme_bw()

# un-comment to save / export figure
#ggsave("Fig5.tiff", plot =Fig5, device = "tiff", scale = 1, width = 16, height  = 11, units = "cm", dpi = 300)

########################################################################################################################
# End of R script.