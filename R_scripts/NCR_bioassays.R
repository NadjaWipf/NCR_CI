####################################################################################################################################################
## This R code is to replicate the analyses, figures and tables from the research article entitled:
## "Multi-insecticide resistant malaria vectors in the field remain susceptible to malathion, despite the presence of Ace1 point mutations", 2021
## Analysis: insecticide bioassay
## Figures: Fig 2 (A,B,C), Table 1, STable 
## Code written by NW and PM
####################################################################################################################################################

# --------------------------------------------------------------------------------------------------------------
# load required packages
library(ggpubr)
library(lme4)
library(MASS)
library(tidyverse)

# useful function to keep factors as factors after merging different dfs using dyplr::bind_rows()
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

# set working directory to a "data" folder where you have saved all relevant STables "input_data" as .txt files.
# setwd("~/data")

# --------------------------------------------------------------------------------------------------------------
### Fig 2. Phenotypic insecticide resistance in Anopheles coluzzii from southern Côte d'Ivoire assessed with WHO discriminating concentrations and dose-response bioassays.

##############################
##  Data analysis  Fig. 2A  ##
##############################

# import bioassay data "WHO insecticide susceptibility assays using  discriminating concentration" (Fig. 2A)

# Note: Because we had 100% mortality for one bioassay (Malathion in Tiassalé, 2018),
# the GLM could not properly estimate the confidence interval with the original dataset,
# therefore we artificially added one survivor into the data set and stated this in the figure legend.

dc <- read.table("discriminating_conc_input_data.txt", sep="\t", header=T)

dc <- dc%>%mutate_if(is.character, as.factor)

## subset Bendiocarb data
dc_ben <- filter(dc, Insecticide == "Bendiocarb")

# GLM for Bendiocarb subset
ben.s <- glm(cbind(Dead,Alive) ~ 0 + Site, data = dc_ben, family=binomial(link = "logit"))
summary(ben.s)

# compute dispersion for binomial model
N=nrow(dc_ben)
p=length(coef(ben.s))+1
sum(resid(ben.s, type="pearson")^2)/(N-p)

# calculate average mortality and 95% confidence interval and create table
ci.ben.s <- data.frame(
  Outcome="Proportion dead",
  Insecticide="0.1% Bendiocarb (CB)",
  Site=c("Agboville", "Dabou", "Tiassalé"),
  Estimate=plogis(summary(ben.s)$coefficients[,"Estimate"]),
  CI.low=plogis(summary(ben.s)$coefficients[,"Estimate"]-1.96*summary(ben.s)$coefficients[,"Std. Error"]),
  CI.high=plogis(summary(ben.s)$coefficients[,"Estimate"]+1.96*summary(ben.s)$coefficients[,"Std. Error"])
)

## subset DDT data
dc_ddt <- filter(dc, Insecticide == "DDT")

# GLM for DDT subset
ddt.s <- glm(cbind(Dead,Alive) ~ 0 + Site, data = dc_ddt, family=binomial(link = "logit"))
summary(ddt.s)

# compute dispersion for binomial model
N=nrow(dc_ddt)
p=length(coef(ddt.s))+1
sum(resid(ddt.s, type="pearson")^2)/(N-p)

# calculate average mortality and 95% confidence interval and create table
ci.ddt.s <- data.frame(
  Outcome="Proportion dead",
  Insecticide="4% DDT (OC)",
  Site=c("Agboville", "Dabou", "Tiassalé"),
  Estimate=plogis(summary(ddt.s)$coefficients[,"Estimate"]),
  CI.low=plogis(summary(ddt.s)$coefficients[,"Estimate"]-1.96*summary(ddt.s)$coefficients[,"Std. Error"]),
  CI.high=plogis(summary(ddt.s)$coefficients[,"Estimate"]+1.96*summary(ddt.s)$coefficients[,"Std. Error"])
)

## subset deltamethrin data
dc_del <- filter(dc, Insecticide == "Deltamethrin")

# GLM for deltamethrin subset
del.s <- glm(cbind(Dead,Alive) ~ 0 + Site, data = dc_del, family=binomial(link = "logit"))
summary(del.s)

# compute dispersion for binomial model
N=nrow(dc_del)
p=length(coef(del.s))+1
sum(resid(del.s, type="pearson")^2)/(N-p)

# calculate average mortality and 95% confidence interval and create table
ci.del.s <- data.frame(
  Outcome="Proportion dead",
  Insecticide="0.05% Deltamethrin (PY)",
  Site=c("Agboville", "Dabou", "Tiassalé"),
  Estimate=plogis(summary(del.s)$coefficients[,"Estimate"]),
  CI.low=plogis(summary(del.s)$coefficients[,"Estimate"]-1.96*summary(del.s)$coefficients[,"Std. Error"]),
  CI.high=plogis(summary(del.s)$coefficients[,"Estimate"]+1.96*summary(del.s)$coefficients[,"Std. Error"])
)

## subset Malathion data
dc_mal <- filter(dc, Insecticide == "Malathion")

# GLM for Malathion subset
mal.s <- glm(cbind(Dead,Alive) ~ 0 + Site, data = dc_mal, family=binomial(link = "logit"))
summary(mal.s)

# compute dispersion for binomial model
N=nrow(dc_mal)
p=length(coef(mal.s))+1
sum(resid(mal.s, type="pearson")^2)/(N-p)

# calculate average mortality and 95% confidence interval and create table
ci.mal.s <- data.frame(
  Outcome="Proportion dead",
  Insecticide="5% Malathion (OP)",
  Site=c("Agboville", "Dabou", "Tiassalé"),
  Estimate=plogis(summary(mal.s)$coefficients[,"Estimate"]),
  CI.low=plogis(summary(mal.s)$coefficients[,"Estimate"]-1.96*summary(mal.s)$coefficients[,"Std. Error"]),
  CI.high=plogis(summary(mal.s)$coefficients[,"Estimate"]+1.96*summary(mal.s)$coefficients[,"Std. Error"])
)

## combine results of all 4 GLMs
ci.all.glm <- bind_rows_keep_factors(ci.ben.s, ci.ddt.s, ci.del.s, ci.mal.s)

# summarize total number of exposed  mosquitoes (n_exp) per Insecticide & Site
Exp <- dc %>%
  group_by(Insecticide, Site)%>%
  summarise(n_exp = sum(Exposed))

# add n_exp to results table, rename, summarise and export results
colnames(ci.all.glm)
ci.all.glm <- ci.all.glm %>%
  mutate(Exposed = Exp$n_exp)%>%
  rename(AVG.mortality = Estimate)%>%
  select(Insecticide, Site, Exposed, AVG.mortality, CI.low, CI.high)%>%
  mutate_if(is.character, as.factor)

# write.table(ci.all.glm, file = "disc_concentration_glm_output.txt", row.names = F, dec=".", sep="\t")


#########################
##  PLOT GLM Fig. 2A   ##
#########################

# reorder levels alphabetically
levels(ci.all.glm$Insecticide)
ci.all.glm$Insecticide <- fct_relevel(ci.all.glm$Insecticide,c("0.1% Bendiocarb (CB)", "4% DDT (OC)", "0.05% Deltamethrin (PY)", "5% Malathion (OP)"))

# plot discriminating concentration bioassay (Fig. 2A)
p <- ggplot(ci.all.glm, aes(y=100*AVG.mortality, x=Site, shape=Site, colour=Site))+
  geom_hline(yintercept = 98, linetype = "dotted", alpha = .5)+
  geom_hline(yintercept = 90, linetype = "dashed", alpha = .5)+
  geom_hline(yintercept = c(0,50,100), linetype=1, colour="gray60", alpha=.4)+
  geom_point(size=3)+
  geom_errorbar(aes(ymin=100*CI.low, ymax=100*CI.high), width=0.15, size=0.7)+
  scale_colour_manual(values=c("gray60", "gray40", "gray15"))+
  facet_grid(.~Insecticide)+
  scale_y_continuous(limits = c(-3, 100))+
  labs(y = "Mortality [%]",
       shape = "Mosquito field population",
       colour = "Mosquito field population")+
  theme_bw()+
  theme(axis.title.y = element_text(size=13, colour="#000000"),
        axis.title.x = element_blank(),
        axis.text  = element_text(size=11, colour="#000000"),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=10, face="bold", colour="#000000"),
        legend.position = "bottom",
        legend.text = element_text(size=11, colour="#000000"),
        legend.title = element_text(size=11, face="bold", colour="#000000"))
p

# add number of tested mosquitoes
p2 <- p + geom_text(data=ci.all.glm, aes(label = paste("n=", Exposed, sep="")), y= -3, colour="#000000")
p2

# t create multiple colours for strip background. Source: https://github.com/tidyverse/ggplot2/issues/2096
g <- ggplot_gtable(ggplot_build(p2))
# to colour in just the TOP rectangles
strip_t <- which(grepl('strip-t', g$layout$name))
# set the 4 insecticide colours
fills <- c("#E69F00", "#339933", "#D82E2E", "#0483BA")
# some fancy code I don't understand, but it works
k <- 1
for (i in strip_t) {
  j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
  g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}
ggarrange(g)
Fig2A.dc <- g


#####################################
##  Data analysis  Fig. 2B and 2C  ##
#####################################

### Import dose-response input data into R

do.re <- read.table("dose-response_input_data.txt", sep="\t", header=T)

# check numbers
with(do.re,Exposed-(Dead+Alive)) # all 0 = OK
with(do.re,aggregate(do.re$Exposed,
                     by=list(Insecticide=Insecticide,Conc=Conc,Site=Site),
                     sum))
summary(do.re)

# tell R to recognise Date as "dd.mm.yyyy"
do.re$Date <- as.Date(do.re$Date,"%d.%m.%Y") 

# subset data set for each insecticide
del <- filter(do.re, Insecticide == "Deltamethrin")
mal <- filter(do.re, Insecticide == "Malathion")

# fit GLM for deltamethrin
del.fit <- glm(cbind(Dead,Alive) ~ 0 + Site/log(Conc), data=del, family=binomial(link="logit"))
summary(del.fit)
# compute dispersion for binomial del model
sum(resid(del.fit, type="pearson")^2)/((nrow(del))-(length(coef(del.fit))+1))

# fit GLM for malathion
mal.fit <- glm(cbind(Dead,Alive) ~ 0 + Site/log(Conc), data=mal, family=binomial(link="logit"))
summary(mal.fit)
# compute dispersion for binomial malmodel
sum(resid(mal.fit, type="pearson")^2)/((nrow(del))-(length(coef(mal.fit))+1))

# estimate mortality at other insecticide concentrations
new.data <- expand.grid(Conc=seq(0.00005,100,0.0001), # calculate till conc. 100%
                        Site=c("Agboville", "Dabou", "Ngousso", "Tiassale")
                        )

del.pred <- predict(del.fit, new.data, se.fit=T)
est.del <- new.data%>%
  mutate(Insecticide = as.factor("Deltamethrin"))%>%
  mutate(Mortality = plogis(del.pred$fit))%>%
  mutate(CI.lo = plogis(del.pred$fit-1.96*del.pred$se.fit))%>%
  mutate(CI.hi = plogis(del.pred$fit+1.96*del.pred$se.fit))

mal.pred <- predict(mal.fit, new.data, se.fit=T)
est.mal <- new.data%>%
  mutate(Insecticide = as.factor("Malathion"))%>%
  mutate(Mortality = plogis(mal.pred$fit))%>%
  mutate(CI.lo = plogis(mal.pred$fit-1.96*mal.pred$se.fit))%>%
  mutate(CI.hi = plogis(mal.pred$fit+1.96*mal.pred$se.fit))

est.all <- bind_rows_keep_factors(est.del,est.mal)


##############################################      
######   > one plot 2 insecticides <    ######
##############################################

# summarise actually measured mortalities from original data set (for plot)
mort <- do.re %>%
  group_by(Insecticide, Conc, Site)%>%
  summarise(Exposed = sum(Exposed), Dead = sum(Dead))%>%
  mutate(Mortality = Dead/Exposed)

## plot dose-response curve (drc) plot for each insecticide separately

# subset for deltamethrin drc
mort.del <- mort%>%
  mutate_if(is.character, as.factor)%>%
  filter(Insecticide == "Deltamethrin")%>%
  droplevels()

# reorder levels
mort.del$Site <- fct_relevel(mort.del$Site, c("Ngousso","Agboville","Dabou","Tiassale"))
est.del$Site <- fct_relevel(est.del$Site, c("Ngousso","Agboville","Dabou","Tiassale"))

# plot deltamethrin drc (Fig. 2B)
Fig2B.drc.del <- ggplot(data=mort.del, aes(x=Conc, y=Mortality, colour=Site, shape=Site))+
  facet_grid(.~Insecticide)+
  geom_vline(xintercept=0.05, colour="#D82E2E", linetype=3, size=.9)+ #add dotted line to mark diagn.conc.
  geom_hline(yintercept = c(0,.5,1), colour="gray60", linetype=1, alpha=.4)+ #add thicker line for 50% + 100% mortality
  geom_smooth(data=est.del, stat="identity", aes(x=Conc, y=Mortality, ymin=CI.lo, ymax=CI.hi, colour=Site, linetype=Site), size=0.7)+
  geom_point(size=2, alpha =.9)+
  scale_colour_manual(
    name="Mosquito population",
    values=c("gray75", "gray60", "gray40", "gray15"),
    labels=c("Ngousso (lab)","Agboville (field)","Dabou (field)", "Tiassalé (field)"))+
  scale_shape_manual(
    name="Mosquito population",
    values=c(9,19,17,15),
    labels=c("Ngousso (lab)","Agboville (field)","Dabou (field)", "Tiassalé (field)"))+
  scale_linetype_manual(
    name="Mosquito population",
    values=c(6,2,5,1),
    labels=c("Ngousso (lab)","Agboville (field)","Dabou (field)", "Tiassalé (field)"))+
  ylab("Mortality [%]")+
  xlab("Insecticide concentration [%]")+
  scale_x_continuous(breaks=c(0.0001,0.0004, 0.0008, 0.0016, 0.0063, 0.0125, 0.025, 0.05, 0.1, 0.2, 0.4, 0.8, 1.6, 6.4, 12.8, 100),
                     labels=c("0.0001","0.0004", "0.0008", "0.0016", "0.0063", "0.0125", "0.025", "0.05", "0.1", "0.2", "0.4", "0.8", "1.6", "6.4", "12.8", "100"),
                     limits = c(0.00005,100.1), trans="log2")+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels=c("0","25", "50", "75", "100"))+
  theme_bw()+
  theme(axis.title = element_text(colour="#000000", size=13),
        axis.text.y  = element_text(colour="#000000", size= 11),
        axis.text.x = element_text(colour="#000000", angle=40, hjust=1, vjust=1, size = 11),
        legend.title = element_text(colour="#000000", face="bold", size = 11),
        legend.text = element_text(colour="#000000", size = 11),
        legend.position = "bottom",
        strip.text = element_text(colour="#000000", face="bold", size=13),
        strip.background = element_rect(fill="#D82E2E"))

#######################
###  Malathion drc  ### 
#######################

# subset for malathion drc
mort.mal <- mort%>%
  mutate_if(is.character, as.factor)%>%
  filter(Insecticide == "Malathion")%>%
  droplevels()

# reorder levels
mort.mal$Site <- fct_relevel(mort.mal$Site, c("Ngousso","Tiassale","Agboville","Dabou"))
est.mal$Site <- fct_relevel(est.mal$Site, c("Ngousso","Tiassale","Agboville","Dabou"))

# plot deltamethrin drc (Fig. 2C)
Fig2C.drc.mal <- ggplot(data=mort.mal, aes(x=Conc, y=Mortality, colour=Site, shape=Site))+
  facet_grid(.~Insecticide)+
  geom_vline(xintercept=5, colour="#0483BA", linetype=3, size=.9)+ #add dotted line to mark diagn.conc.
  geom_hline(yintercept = c(0,.5,1), colour="gray60", linetype=1, alpha=.4)+ #add thicker line for 50% + 100% mortality
  geom_smooth(data=est.mal, stat="identity", aes(x=Conc, y=Mortality, ymin=CI.lo, ymax=CI.hi, colour=Site, linetype=Site), size=0.7)+
  geom_point(size=2, alpha =.9)+
  scale_colour_manual(
    name="Mosquito population",
    values=c("gray75","gray15", "gray60", "gray40"),
    labels=c("Ngousso (lab)", "Tiassalé (field)", "Agboville (field)", "Dabou (field)"))+
  scale_shape_manual(
    name="Mosquito population",
    values=c(9,15,19,17),
    labels=c("Ngousso (lab)", "Tiassalé (field)", "Agboville (field)", "Dabou (field)"))+
  scale_linetype_manual(
    name="Mosquito population",
    values=c(6,1,2,5),
    labels=c("Ngousso (lab)", "Tiassalé (field)", "Agboville (field)", "Dabou (field)"))+
  ylab("Mortality [%]")+
  xlab("Insecticide concentration [%]")+
  scale_x_continuous(breaks=c(0.0391,0.1563, 0.3125, 0.625, 1.25, 2.5, 5, 20, 100),
                     labels=c("0.0391","0.1563", "0.3125", "0.625", "1.25", "2.5", "5", "20", "100"),
                     limits = c(0.00005,100.1), trans="log2")+
  scale_y_continuous(breaks=c(0, 0.25, 0.5, 0.75, 1),
                     labels=c("0","25", "50", "75", "100"))+
  theme_bw()+
  theme(axis.title = element_text(colour="#000000", size=13),
        axis.text.y  = element_text(colour="#000000", size= 11),
        axis.text.x = element_text(colour="#000000", angle=40, hjust=1, vjust=1, size = 11),
        legend.title = element_text(colour="#000000", face="bold", size = 11),
        legend.text = element_text(colour="#000000", size = 11),
        legend.position = "bottom",
        strip.text = element_text(colour="#000000", face="bold", size=13),
        strip.background = element_rect(fill="#0483BA"))


# arrange the three plots below each other (Fig. 2 with panels A, B and C)
Fig2.ABC <- ggarrange(Fig2A.dc, Fig2B.drc.del, Fig2C.drc.mal, nrow = 3, labels = c('A', 'B', 'C'))

#Fig2.ABC + ggsave("Fig2.ABC.pdf", device = "pdf", scale = 1, width = 16, height = 22, units = "cm", dpi = 300)
#Fig2.ABC + ggsave("Fig2.ABC.tiff", device = "tiff", scale = 1, width = 20, height = 30, units = "cm", dpi = 300)

###############################      
###  Get the LC50 & LC90    ###
###############################

# Create data set with lethal concentrations (LCs 50 and 90) and 95% confidence intervals (CIs)
lds=function(mod, cf) {
  tmp=dose.p(mod,p=c(.5,.9), cf=cf)
  LC50=as.vector(tmp[1])
  SE50=attr(tmp,"SE")[1]
  LC50.LCI=as.vector(tmp[1])-1.96*attr(tmp,"SE")[1]
  LC50.UCI=as.vector(tmp[1])+1.96*attr(tmp,"SE")[1]
  LC90=as.vector(tmp[2])
  SE90=attr(tmp,"SE")[2]
  LC90.LCI=as.vector(tmp[2])-1.96*attr(tmp,"SE")[2]
  LC90.UCI=as.vector(tmp[2])+1.96*attr(tmp,"SE")[2]
  cbind(LC50, SE50, LC50.LCI, LC50.UCI, LC90, SE90, LC90.LCI, LC90.UCI)
}
lcs=as.data.frame(
  rbind(
    lds(del.fit, cf=c(1,5)),
    lds(del.fit, cf=c(2,6)),
    lds(del.fit, cf=c(3,7)),
    lds(del.fit, cf=c(4,8)),
    lds(mal.fit, cf=c(1,5)),
    lds(mal.fit, cf=c(2,6)),
    lds(mal.fit, cf=c(3,7)),
    lds(mal.fit, cf=c(4,8))
  ))
lcs=exp(lcs)  
lcs=cbind(
  Insecticide=rep(c("Deltamethrin","Malathion"), each=4),
  Population=rep(c("Agboville","Dabou","Ngousso", "Tiassale"),2),
  N=c(
    sum(del$Exposed[del$Site=="Agboville"]),
    sum(del$Exposed[del$Site=="Dabou"]),
    sum(del$Exposed[del$Site=="Ngousso"]),
    sum(del$Exposed[del$Site=="Tiassale"]),
    sum(mal$Exposed[mal$Site=="Agboville"]),
    sum(mal$Exposed[mal$Site=="Dabou"]),
    sum(mal$Exposed[mal$Site=="Ngousso"]),
    sum(mal$Exposed[mal$Site=="Tiassale"])
  ),
  lcs)

# export dose-response curves lethal concentrations glm output (LC 50 for Table 1)
# write.table(lcs,file="drc_lethal_conc_glm_output.txt",row.names=F, dec=".", sep="\t")

# --------------------------------------------------------------------------------------------------------------
# End of R script.