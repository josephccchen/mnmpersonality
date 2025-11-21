## EDIT THIS SECTION -----------------------------------------------------------

ratingsdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/events'
scriptsdir = '/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/fig3_and_4_NumEdges_and_StrengthEdges'
df1 = read.csv(paste0(scriptsdir, '/MnM_fig34_dFC_EdgeAnalyses.csv'))

# font sizes
fontsize = 8
smallfontsize = 2

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
  
df2 = ddply(df1, ~ ID * Task, summarise,
            Extraversion = mean(Extraversion, na.rm=T),
            Extravert = mean(Extravert, na.rm=T),
            AllConn_CEN.mean = mean(AllConn_CEN, na.rm=T),
            AllConn_DMN.mean = mean(AllConn_DMN, na.rm=T),
            AllConn_SN.mean = mean(AllConn_SN, na.rm=T),
            Intra_CEN.mean = mean(Intra_CEN, na.rm=T),
            Intra_DMN.mean = mean(Intra_DMN, na.rm=T),
            Intra_SN.mean = mean(Intra_SN, na.rm=T),
            Extra_CEN.mean = mean(Extra_CEN, na.rm=T),
            Extra_DMN.mean = mean(Extra_DMN, na.rm=T),
            Extra_SN.mean = mean(Extra_SN, na.rm=T),
            CEN_SN.mean = mean(CEN_SN, na.rm=T),
            SN_DMN.mean = mean(SN_DMN, na.rm=T),
            DMN_CEN.mean = mean(DMN_CEN, na.rm=T),
            Intra_CEN_z.mean = mean(Intra_CEN_z, na.rm=T),
            Intra_DMN_z.mean = mean(Intra_DMN_z, na.rm=T),
            Intra_SN_z.mean = mean(Intra_SN_z, na.rm=T),
            Extra_CEN_z.mean = mean(Extra_CEN_z, na.rm=T),
            Extra_DMN_z.mean = mean(Extra_DMN_z, na.rm=T),
            Extra_SN_z.mean = mean(Extra_SN_z, na.rm=T),
            CEN_SN_z.mean = mean(CEN_SN_z, na.rm=T),
            SN_DMN_z.mean = mean(SN_DMN_z, na.rm=T),
            DMN_CEN_z.mean = mean(DMN_CEN_z, na.rm=T))

df3 = data.frame(df2, TrialAvgRating = NA)
AllSubj = unique(df3$ID)
NumSubj = length(AllSubj)
  
## read in behavioural focus data ----
for (iRow in seq(1,dim(df3)[1])) {
    tempSubj = df3$ID[iRow]
    tempTask = df3$Task[iRow]
    tempAttnFile = paste(ratingsdir, '/', tempSubj, '_task-', tempTask, '_events.csv', sep='')
    tempattndf = read.csv(tempAttnFile)
    tempratings = tempattndf$rating[tempattndf$trial_type == 'response']
    tempallratingstext = paste(tempratings, collapse = ", ")
    tempavg = mean(as.numeric(tempratings))
    df3$TrialAvgRating[iRow] = tempavg
}

## stats
interestcomps = c("BreathNVHA Extrovert0 - BreathNVHA Extrovert1",
                  "BreathNVLA Extrovert0 - BreathNVLA Extrovert1",
                  "Meditation Extrovert0 - Meditation Extrovert1")

## Fig3 NumEdges ----------
PublicationColumnsOfInterest = c("Intra_CEN.mean", "Intra_DMN.mean", "Intra_SN.mean",
                      "Extra_CEN.mean", "Extra_DMN.mean", "Extra_SN.mean",
                      "CEN_SN.mean", "SN_DMN.mean", "DMN_CEN.mean")
PrettyNames = c("Internal CEN Edges", "Internal DMN Edges", "Internal SN Edges",
                "External CEN Edges", "External DMN Edges", "External SN Edges",
                "CEN-DMN Edges", "DMN-SN Edges", "SN-CEN Edges")

PublicationNumColsOfInterest = length(PublicationColumnsOfInterest)
PublicationGraphs = list(length = PublicationNumColsOfInterest)
for (iCol in seq(1,PublicationNumColsOfInterest)) {
  tempColName = PublicationColumnsOfInterest[iCol]
  tempColIndex = match(tempColName, colnames(df3))
  
  # create tempdf
  tempdf = data.frame(ID = df3$ID,
                      Task = df3$Task,
                      Extraversion = df3$Extraversion,
                      Extravert = df3$Extravert,
                      measure = df3[,tempColIndex],
                      TrialAvgRating = df3$TrialAvgRating)
  tempdf$Extravert = as.factor(tempdf$Extravert)
  
  # LMM stats
  tempmodel = lmer(TrialAvgRating ~ measure*Extravert + (1|ID), tempdf)
  tempmodelsummary = summary(tempmodel)
  tempanova = Anova(tempmodel, type = 3)
  tempb0 = tempmodelsummary$coefficients[1,1]
  tempmeasure = tempmodelsummary$coefficients[2,1]
  tempExtravert1 = tempmodelsummary$coefficients[3,1]
  tempInteractionChi = tempanova$Chisq[4]
  tempInteraction = tempmodelsummary$coefficients[4,1]
  temppvalue = tempanova$`Pr(>Chisq)`[4]
  
  # lines of best fit
  tempIntrovertb0 = tempb0
  tempIntrovertb1 = tempmeasure
  tempExtravertb0 = tempb0 + tempExtravert1
  tempExtravertb1 = tempmeasure + tempInteraction
  
  # Pearson Correlation
  tempcor_Introvert = cor.test(tempdf$measure[tempdf$Extravert == "0"], tempdf$TrialAvgRating[tempdf$Extravert == "0"])
  tempcor_Introvert_r = tempcor_Introvert$estimate
  #tempcor_Introvert_p = tempcor_Introvert$p.value
  tempcor_Extravert = cor.test(tempdf$measure[tempdf$Extravert == "1"], tempdf$TrialAvgRating[tempdf$Extravert == "1"])
  tempcor_Extravert_r = tempcor_Extravert$estimate
  #tempcor_Extravert_p = tempcor_Extravert$p.value
  
  # plot details
  #threshold99 = quantile(tempdf$measure, 0.99, na.rm=T)
  #tempxmax = threshold99
  tempxmax = max(tempdf$measure, na.rm=T)
  tempxmin = min(tempdf$measure, na.rm=T)
  tempxmid = (tempxmax+tempxmin) /2
  
  # labels to annotate plot
  LMMlabel = paste0("X2 = ", sprintf("%.1f", tempInteractionChi), ", p = ", sprintf("%.3f", temppvalue))
  PearsonIntroLabel = paste0("Introvert r = ", sprintf("%.2f",tempcor_Introvert_r))
  PearsonExtraLabel = paste0("Extravert r = ", sprintf("%.2f",tempcor_Extravert_r))

  # actual plot
  tempplot = ggplot(tempdf, aes(x = measure, y = TrialAvgRating, group = Extravert, colour = Extravert, shape = Extravert)) +
    geom_point(alpha = 0.3) +
    theme_classic() +
    geom_abline(intercept = tempIntrovertb0, slope = tempIntrovertb1, colour = "black") +
    geom_abline(intercept = tempExtravertb0, slope = tempExtravertb1, colour = "grey") +
    scale_x_continuous(name = PrettyNames[iCol], limits = c(tempxmin, tempxmax)) +
    scale_y_continuous(name = "Attention", limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_colour_manual(values = c("black", "grey")) +
    scale_shape_manual(values = c(16,17)) +
    annotate('text', x = tempxmid, y = 0, label = LMMlabel, hjust = 0.5, vjust = 0, size = 3) +
    #annotate('text', x = tempxmid, y = 0.10, label = PearsonIntroLabel, hjust = 0.5, vjust=0, size = 3, colour = "blue") +
    #annotate('text', x = tempxmid, y = 0.20, label = PearsonExtraLabel, hjust = 0.5, vjust=0, size = 3, colour = "red") +
    theme(legend.title = element_blank(),
          legend.position = "none")
  PublicationGraphs[[iCol]] = tempplot
}

legenddf = data.frame(value = runif(12),
                      task = c("a", "a", "b", "b"),
                      group = c("1", "2"))
legend_plot = ggplot(legenddf, aes(x = task, y = value, group = group, colour = group, shape = group)) + 
  geom_point() + geom_line() + theme_classic() +
  scale_shape_manual(values = c(16,17), label = c("Introvert", "Extravert")) +
  scale_colour_manual(values = c( "black", "grey"), 
                      label= c("Introvert", "Extravert")) +
  theme(legend.position = "top", legend.title = element_blank())
legend_g = get_legend(legend_plot)

PublicationBigPlot = plot_grid(PublicationGraphs[[1]], PublicationGraphs[[2]], PublicationGraphs[[3]],
                               PublicationGraphs[[4]], PublicationGraphs[[5]], PublicationGraphs[[6]], 
                               ncol = 3, labels = "AUTO")
PublicationBigPlot2 = plot_grid(legend_g, PublicationBigPlot, nrow = 2, rel_heights = c(0.05,1))
PublicationGraphs_savename = paste0(scriptsdir, '/MnM_fig3_dFC_NumEdges.pdf')
ggsave(PublicationGraphs_savename, plot = PublicationBigPlot2, device = 'pdf', width = 15, height = 9, units = 'cm')


## Fig4 StrengthEdges ----------
SecondColumnsOfInterest = c("Intra_CEN_z.mean", "Intra_DMN_z.mean", "Intra_SN_z.mean",
                           "Extra_CEN_z.mean", "Extra_DMN_z.mean", "Extra_SN_z.mean",
                           "CEN_SN_z.mean", "SN_DMN_z.mean", "DMN_CEN_z.mean")

PrettyNames = c("Internal CEN FC (z)", "Internal DMN FC (z)", "Internal SN FC (z)",
                "External CEN FC (z)", "External DMN FC (z)", "External SN FC (z)",
                "CEN-DMN FC (z)", "DMN-SN FC (z)", "SN-CEN FC (z)")

SecondNumColsOfInterest = length(SecondColumnsOfInterest)
SecondGraphs = list(length = SecondNumColsOfInterest)
for (iCol in seq(1,SecondNumColsOfInterest)) {
  tempColName = SecondColumnsOfInterest[iCol]
  tempColIndex = match(tempColName, colnames(df3))
  
  # create tempdf
  tempdf = data.frame(ID = df3$ID,
                      Task = df3$Task,
                      Extraversion = df3$Extraversion,
                      Extravert = df3$Extravert,
                      measure = df3[,tempColIndex],
                      TrialAvgRating = df3$TrialAvgRating)
  tempdf$Extravert = as.factor(tempdf$Extravert)
  
  # LMM stats
  tempmodel = lmer(TrialAvgRating ~ measure*Extravert + (1|ID), tempdf)
  tempmodelsummary = summary(tempmodel)
  tempanova = Anova(tempmodel, type = 3)
  tempb0 = tempmodelsummary$coefficients[1,1]
  tempmeasure = tempmodelsummary$coefficients[2,1]
  tempExtravert1 = tempmodelsummary$coefficients[3,1]
  tempInteractionChi = tempanova$Chisq[4]
  tempInteraction = tempmodelsummary$coefficients[4,1]
  temppvalue = tempanova$`Pr(>Chisq)`[4]
  
  # lines of best fit
  tempIntrovertb0 = tempb0
  tempIntrovertb1 = tempmeasure
  tempExtravertb0 = tempb0 + tempExtravert1
  tempExtravertb1 = tempmeasure + tempInteraction
  
  # Pearson Correlation
  tempcor_Introvert = cor.test(tempdf$measure[tempdf$Extravert == "0"], tempdf$TrialAvgRating[tempdf$Extravert == "0"])
  tempcor_Introvert_r = tempcor_Introvert$estimate
  #tempcor_Introvert_p = tempcor_Introvert$p.value
  tempcor_Extravert = cor.test(tempdf$measure[tempdf$Extravert == "1"], tempdf$TrialAvgRating[tempdf$Extravert == "1"])
  tempcor_Extravert_r = tempcor_Extravert$estimate
  #tempcor_Extravert_p = tempcor_Extravert$p.value
  
  # plot details
  #threshold99 = quantile(tempdf$measure, 0.99, na.rm=T)
  #tempxmax = threshold99
  tempxmax = max(tempdf$measure, na.rm=T)
  tempxmax = 2.2
  tempxmin = min(tempdf$measure, na.rm=T)
  tempxmin = 0.8
  tempxmid = (tempxmax+tempxmin) /2
  
  # labels to annotate plot
  LMMlabel = paste0("X2 = ", sprintf("%.1f", tempInteractionChi), ", p = ", sprintf("%.3f", temppvalue))
  PearsonIntroLabel = paste0("Introvert r = ", signif(tempcor_Introvert_r,3))
  PearsonExtraLabel = paste0("Extravert r = ", signif(tempcor_Extravert_r,3))
  
  # actual plot
  tempplot = ggplot(tempdf, aes(x = measure, y = TrialAvgRating, group = Extravert, colour = Extravert, shape = Extravert)) +
    geom_point(alpha = 0.3) +
    theme_classic() +
    geom_abline(intercept = tempIntrovertb0, slope = tempIntrovertb1, colour = "black") +
    geom_abline(intercept = tempExtravertb0, slope = tempExtravertb1, colour = "grey") +
    scale_x_continuous(name = PrettyNames[iCol], limits = c(tempxmin, tempxmax), 
                       breaks = seq(tempxmin, tempxmax, 0.4)) +
    scale_y_continuous(name = "Attention", limits = c(0,1), breaks = seq(0,1,0.2)) +
    scale_colour_manual(values = c("black", "grey")) +
    scale_shape_manual(values = c(16,17)) +
    annotate('text', x = tempxmid, y = 0.00, label = LMMlabel, hjust = 0.5, vjust = 0, size = 3) +
    #annotate('text', x = tempxmid, y = 0.10, label = PearsonIntroLabel, hjust = 0.5, vjust=0, size = 3, colour = "blue") +
    #annotate('text', x = tempxmid, y = 0.20, label = PearsonExtraLabel, hjust = 0.5, vjust=0, size = 3, colour = "red") +
    theme(legend.title = element_blank(),
          legend.position = "none")
  tempplot
  SecondGraphs[[iCol]] = tempplot
}

SecondBigPlot = plot_grid(SecondGraphs[[1]], SecondGraphs[[2]], SecondGraphs[[3]],
                          SecondGraphs[[4]], SecondGraphs[[5]], SecondGraphs[[6]], 
                          ncol = 3, labels = "AUTO")
SecondBigPlot2 = plot_grid(legend_g, SecondBigPlot, nrow = 2, rel_heights = c(0.05,1))
SecondGraphs_savename = paste0(scriptsdir, '/MnM_fig4_dFC_StrengthEdges.pdf')
ggsave(SecondGraphs_savename, plot = SecondBigPlot2, device = 'pdf', width = 15, height = 9, units = 'cm')

