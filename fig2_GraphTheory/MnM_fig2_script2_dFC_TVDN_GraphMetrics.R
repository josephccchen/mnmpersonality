## EDIT THIS SECTION -----------------------------------------------------------

inputdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/MnM_dFC/r_0.5/k_1.6'
ratingsdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/events'
outputdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/fig2_GraphTheory'
personalitydf = read.csv("/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/MnM_PersonalityInfo.csv")

setwd(inputdir)

# font sizes
fontsize = 8
smallfontsize = 3

## Set-Up ----------------------------------------------------------------------
# libraries
library(ggplot2)
library(reshape2)
library(cowplot)
library(gridExtra)
library(dplyr)
library(DescTools)
library(plyr)
library(lme4)
library(multcomp)
library(car)
library(emmeans)
library(MKinfer)
library(compute.es)
library(plotrix)
library(tidyr)
library(abind)



## read in graph metrics from matlab script ----
graphdfname = paste(inputdir, '/graphmatrics_summary.csv', sep='')
graphdf = read.csv(graphdfname)
  
## read in state switching data from python script ----
df1name = paste(inputdir, '/df1.csv',sep='')
df1 = read.csv(df1name)
df1 = df1[, !(names(df1) %in% "X")] #drop the first index column
df1 = data.frame(df1, 
                   Extraversion = NA,
                   Extravert = NA,
                   TrialRatings = NA,
                   TrialAvgRating = NA)

# remove sub-2003 sub-2020
df1 = subset(df1, ID != "sub-2003" | ID != "sub-2020")

  
## read in behavioural focus data ----
for (iRow in seq(1,dim(df1)[1])) {
    tempSubj = df1$ID[iRow]
    tempTask = df1$Task[iRow]
    tempAttnFile = paste(ratingsdir, '/', tempSubj, '_task-', tempTask, '_events.csv', sep='')
    tempattndf = read.csv(tempAttnFile)
    tempratings = tempattndf$rating[tempattndf$trial_type == 'response']
    tempallratingstext = paste(tempratings, collapse = ", ")
    tempavg = mean(as.numeric(tempratings))
    df1$TrialRatings[iRow] = tempallratingstext
    df1$TrialAvgRating[iRow] = tempavg
    
}
  
# add personality data to the dataframe ----
for (iRow in seq(1, dim(df1)[1])) {
  tempSubj = df1$ID[iRow]
  tempIndex = match(tempSubj, personalitydf$ID, nomatch = 0)
  if (tempIndex == 0) {next}
  df1$Extraversion[iRow] = personalitydf$Personality_ExtraversionRaw[tempIndex]
}
  
for (iRow in seq(1,dim(df1)[1])) {
  if (is.na(df1$Extraversion[iRow])) { next}
  if (df1$Extraversion[iRow] < 24) {df1$Extravert[iRow] = 0}
  if (df1$Extraversion[iRow] ==24) {df1$Extravert[iRow] = NA}
  if (df1$Extraversion[iRow] > 24) {df1$Extravert[iRow] = 1}
}
  

## combine all the data into one df ----
alldf = data.frame(df1,
                   GEFF = graphdf$GEFF,
                   Modularity = graphdf$Modularity)

## remove rows with NAN in Extravert from df1 ----
alldf2 = alldf[!is.na(df1$Extravert),]
alldf2$Extravert = as.factor(alldf2$Extravert)
alldf$Extravert = as.factor(alldf$Extravert)

## make summarydf ----

alldf_summ = ddply(alldf2, ~ Task * Extravert, summarise,
                   NumECPTS.mean = mean(NumECPTS, na.rm = T), NumECPTS.se = std.error(NumECPTS, na.rm = T),
                   GEFF.mean = mean(GEFF), GEFF.se = std.error(GEFF),
                   Modularity.mean = mean(Modularity), Modularity.se = std.error(Modularity))

alldf_summ$Extravert = as.factor(alldf_summ$Extravert)

alldfIntro = subset(alldf2, alldf2$Extravert == "0")
alldfExtra = subset(alldf2, alldf2$Extravert == "1")

## Legend plot ----

legenddf = data.frame(value = runif(12),
                      task = c("a", "a", "b", "b"),
                      group = c("1", "2"))
legend_plot = ggplot(legenddf, aes(x = task, y = value, group = group, colour = group, shape = group)) + 
  geom_point() + geom_line() + theme_classic() +
  scale_shape_manual(values = c(16,17), labels = c("Introvert", "Extravert")) +
  scale_colour_manual(values = c( "black", "grey"), 
                      label= c("Introvert", "Extravert")) +
  theme(legend.position = "top", legend.title = element_blank())
legend_g = get_legend(legend_plot)


## GEFF ~ Attn ---------

# lme
GEFFxAttn_lme = lmer(GEFF ~ TrialAvgRating * Extravert + (1|ID), alldf)
GEFFxAttn_ANOVA = Anova(GEFFxAttn_lme,type=3)
GEFFxAttn_summary = summary(GEFFxAttn_lme)
GEFFxAttn_InteractionChi = GEFFxAttn_ANOVA$Chisq[4]
GEFFxAttn_Interactionp = GEFFxAttn_ANOVA$`Pr(>Chisq)`[4]
GEFFxAttn_IntroB0 = GEFFxAttn_summary$coefficients[1,1]
GEFFxAttn_IntroB1 = GEFFxAttn_summary$coefficients[2,1]
GEFFxAttn_ExtraB0 = GEFFxAttn_summary$coefficients[3,1] + GEFFxAttn_IntroB0
GEFFxAttn_ExtraB1 = GEFFxAttn_summary$coefficients[4,1] + GEFFxAttn_IntroB1
GEFFxAttn_label = paste0("Interaction X2 = ", sprintf("%.1f", GEFFxAttn_InteractionChi) , ", p = ", sprintf("%.3f", GEFFxAttn_Interactionp))
GEFFxAttn_Intro_cor = cor(alldfIntro$TrialAvgRating, alldfIntro$GEFF, use = 'pairwise.complete.obs')
GEFFxAttn_Extra_cor = cor(alldfExtra$TrialAvgRating, alldfExtra$GEFF, use = 'pairwise.complete.obs')
GEFFxAttn_Intro_label = paste("Introvert r = ", round(GEFFxAttn_Intro_cor,2))
GEFFxAttn_Extra_label = paste("Extravert r = ", round(GEFFxAttn_Extra_cor,2))

# now plot
GEFFxlimits = c(0.40,0.60,0.04)
GEFFxAttnplot = ggplot(alldf2, aes(x = TrialAvgRating, y = GEFF , group = Extravert, colour = Extravert, shape = Extravert)) + 
  geom_point(alpha = 0.5) +
  scale_x_continuous(limits = c(0.4,1), breaks = seq(0.4,1,0.2), name = "Attention") +
  scale_y_continuous(name = "Global Efficiency",limits = c(GEFFxlimits[1], GEFFxlimits[2]), breaks = seq(GEFFxlimits[1],GEFFxlimits[2],GEFFxlimits[3])) +
  scale_colour_manual(values = c("black", "grey"), labels = c("Introvert", "Extraverts")) +
  scale_shape_manual(values = c(16,17), labels = c("Introverts", "Extraverts")) +
  theme_classic() +
  theme(legend.position = "none", legend.title = element_blank()) +
  annotate('text', x=0.7, y= 0.4, label=GEFFxAttn_label, size = smallfontsize, hjust = 0.5, vjust = 0, fontface = 'italic') + 
  geom_abline(slope = GEFFxAttn_IntroB1, intercept = GEFFxAttn_IntroB0, colour = 'black') +
  annotate('text', x=0.7, y= 0.42, label=GEFFxAttn_Intro_label, colour = "black", size = smallfontsize, hjust = 0.5, vjust = 0, fontface = 'italic') + 
  geom_abline(slope = GEFFxAttn_ExtraB1, intercept = GEFFxAttn_ExtraB0, colour = 'grey') +
  annotate('text', x=0.7, y= 0.44, label=GEFFxAttn_Extra_label, colour = "grey", size = smallfontsize, hjust = 0.5, vjust =0, fontface = 'italic')



## Modularity ~ Attn -------
ModularityxAttn_lme = lmer(Modularity ~ TrialAvgRating * Extravert + (1|ID), alldf)
ModularityxAttn_ANOVA = Anova(ModularityxAttn_lme,type=3)
ModularityxAttn_summary = summary(ModularityxAttn_lme)
ModularityxAttn_InteractionChi = ModularityxAttn_ANOVA$Chisq[4]
ModularityxAttn_Interactionp = ModularityxAttn_ANOVA$`Pr(>Chisq)`[4]
ModularityxAttn_IntroB0 = ModularityxAttn_summary$coefficients[1,1]
ModularityxAttn_IntroB1 = ModularityxAttn_summary$coefficients[2,1]
ModularityxAttn_ExtraB0 = ModularityxAttn_summary$coefficients[3,1] + ModularityxAttn_IntroB0
ModularityxAttn_ExtraB1 = ModularityxAttn_summary$coefficients[4,1] + ModularityxAttn_IntroB1
ModularityxAttn_label = paste0("Interaction X2 = ", sprintf("%.1f", ModularityxAttn_InteractionChi) , ", p = ", sprintf("%.3f", GEFFxAttn_Interactionp))
ModularityxAttn_Intro_cor = cor(alldfIntro$TrialAvgRating, alldfIntro$Modularity, use = 'pairwise.complete.obs')
ModularityxAttn_Extra_cor = cor(alldfExtra$TrialAvgRating, alldfExtra$Modularity, use = 'pairwise.complete.obs')
ModularityxAttn_Intro_label = paste("Introvert r = ", round(ModularityxAttn_Intro_cor,2))
ModularityxAttn_Extra_label = paste("Extravert r = ", round(ModularityxAttn_Extra_cor,2))

# now plot
Modularityxlimits = c(0.18,0.44,0.04)
ModularityxAttnplot = ggplot(alldf2, aes(y = Modularity, x = TrialAvgRating, group = Extravert, colour = Extravert, shape = Extravert)) + 
  geom_point(alpha = 0.5) + 
  scale_x_continuous(limits = c(0.4,1), breaks = seq(0.4,1,0.2), name = "Attention") +
  scale_y_continuous(name = "Modularity",limits = c(Modularityxlimits[1], Modularityxlimits[2]), breaks = seq(Modularityxlimits[1],Modularityxlimits[2],Modularityxlimits[3])) +
  scale_colour_manual(values = c("black", "grey"), labels = c("Introvert", "Extravert")) +
  scale_shape_manual(values = c(16,17), labels = c("Introverts", "Extraverts")) +
  theme_classic() +
  theme(legend.position = "none", legend.title = element_blank()) +
  annotate('text', x=0.7, y= 0.18, label=ModularityxAttn_label, size = smallfontsize, hjust = 0.5, vjust = 0, fontface = 'italic') + 
  geom_abline(slope = ModularityxAttn_IntroB1, intercept = ModularityxAttn_IntroB0, colour = 'black') +
  annotate('text', x=0.7, y= 0.2, label=ModularityxAttn_Intro_label, colour = "black", size = smallfontsize, hjust = 0.5, vjust = 0, fontface = 'italic') + 
  geom_abline(slope = ModularityxAttn_ExtraB1, intercept = ModularityxAttn_ExtraB0, colour = 'grey') +
  annotate('text', x=0.7, y= 0.22, label=ModularityxAttn_Extra_label, colour = "grey", size = smallfontsize, hjust = 0.5, vjust =0, fontface = 'italic')


## megaplot ----
bigplot = plot_grid(GEFFxAttnplot, ModularityxAttnplot,
                    nrow = 1, ncol = 2, labels = "AUTO")
bigbigplot = plot_grid(legend_g, bigplot, ncol = 1, rel_heights = c(0.05,1))

ggsave(paste0(outputdir, '/MnM_fig2_dFC_TVDN_connectivity2networks.pdf'), plot = bigbigplot, width = 15, height = 6, units = 'cm', device = 'pdf')

