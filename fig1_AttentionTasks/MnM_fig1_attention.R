## Script settings (EDIT ME) ----
datadir = "/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data"
outputdir = "/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/fig1_AttentionTasks"
personalitydf = read.csv("/Users/josephchen/Library/CloudStorage/Box-Box/josephchen/MnM/MnMscripts_forpub/data/MnM_PersonalityInfo.csv")

datafiles = list.files(datadir)
fontsize = 8
smallfontsize = 2

InterestComparisons = c("BreathNVHA Extrovert0 - BreathNVHA Extrovert1",
                        "BreathNVLA Extrovert0 - BreathNVLA Extrovert1",
                        "BreathNVHA Extrovert0 - BreathNVLA Extrovert0",
                        "BreathNVHA Extrovert1 - BreathNVLA Extrovert1",
                        "BreathNVLA Extrovert0 - Meditation Extrovert0",
                        "BreathNVHA Extrovert0 - Meditation Extrovert0",
                        "BreathNVLA Extrovert1 - Meditation Extrovert1",
                        "BreathNVHA Extrovert1 - Meditation Extrovert1")

AllSubj = c("sub-2000", "sub-2001", "sub-2005", "sub-2007",
            "sub-2010", "sub-2011", "sub-2012", "sub-2013", "sub-2016",
            "sub-2017", "sub-2018", "sub-2019", "sub-2022",
            "sub-2024", "sub-2025", "sub-2026", "sub-2028", "sub-2029",
            "sub-2031", "sub-2032", "sub-2033", "sub-2034", "sub-2035",
            "sub-2036", "sub-2037", "sub-2038", "sub-2039", "sub-2040",
            "sub-2041", "sub-2042", "sub-2043", "sub-2044", "sub-2045")

NumSubj = length(AllSubj)

## libraries ----
setwd(outputdir)

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

# read the personality data ----
alldf = data.frame(ID = NA,
                   Age = NA,
                   Gender = NA,
                   task = NA,
                   trial = NA,
                   attention = NA,
                   extraversion = NA,
                   extravert = NA)

# Fill data frame with data ----
for (iSubj in seq(1, NumSubj)) {
  tempSubj = AllSubj[iSubj]
  
  # extract personality + demographics data
  tempindex = match(tempSubj, personalitydf$ID)
  tempAge = personalitydf$Age[tempindex]
  tempGender = personalitydf$Gender[tempindex]
  tempextraversion = personalitydf$Personality_ExtraversionRaw[tempindex]
  tempextravert = NA
  if (tempextraversion > 24) {tempextravert = 1}
  if (tempextraversion < 24) {tempextravert = 0}
  
  # retrieve attention data
  tempNoMusic = read.csv(paste0(datadir, "/events/", tempSubj, "_task-Meditation_events.csv"))
  tempLAMusic = read.csv(paste0(datadir, "/events/", tempSubj, "_task-BreathNVLA_events.csv"))
  tempHAMusic = read.csv(paste0(datadir, "/events/", tempSubj, "_task-BreathNVHA_events.csv"))
  
  # ratings
  tempNoMusicRating = suppressWarnings(as.vector(na.omit(as.numeric(tempNoMusic$rating))))
  tempLAMusicRating = suppressWarnings(as.vector(na.omit(as.numeric(tempLAMusic$rating))))
  tempHAMusicRating = suppressWarnings(as.vector(na.omit(as.numeric(tempHAMusic$rating))))
  
  # create temporary objects
  tempSubj = rep(tempSubj, 36)
  tempAge = rep(tempAge, 36)
  tempGender = rep(tempGender,36)
  tempextraversion = rep(tempextraversion, 36)
  tempextravert = rep(tempextravert, 36)
  temptask = rep(c("NoMusic", "LAMusic", "HAMusic"), each = 12)
  temptrial = seq(1,36)
  tempattention = c(tempNoMusicRating, tempLAMusicRating, tempHAMusicRating)
  
  # save into alldf
  alldf <- rbind(alldf, list(tempSubj, tempAge, tempGender, temptask, temptrial, tempattention, tempextraversion, tempextravert))
  
}
alldf <- alldf[-1,]
rownames(alldf) <- NULL
alldf = alldf[complete.cases(alldf),]
alldf$extravert = as.factor(alldf$extravert)
alldf$task = as.factor(alldf$task)
alldf$ID = as.factor(alldf$ID)
alldf$trial = as.factor(alldf$trial)

# fMRIfocus ~ Task ----
alldf_summ = ddply(alldf, ~ task * extravert, summarise,
                   attn.mean = mean(attention),
                   attn.se = std.error(attention))
BetweenComps = c("HAMusic extravert0 - HAMusic extravert1",
                 "LAMusic extravert0 - LAMusic extravert1",
                 "NoMusic extravert0 - NoMusic extravert1")
WithinComps_Intro = c("HAMusic extravert0 - LAMusic extravert0",
                      "LAMusic extravert0 - NoMusic extravert0")
WithinComps_Extro = c("HAMusic extravert1 - LAMusic extravert1",
                      "LAMusic extravert1 - NoMusic extravert1")
model1 = lmer(attention ~ task * extravert + (1|ID), alldf)
Anova(model1, type = 3)
comp = pairs(emmeans(model1, ~ task * extravert), adjust = "none")
compdf = data.frame(comp)
comp_between = which(compdf$contrast %in% BetweenComps)
comp_within_intro = which(compdf$contrast %in% WithinComps_Intro)
comp_within_extro = which(compdf$contrast %in% WithinComps_Extro)
tvals_between = compdf$t.ratio[comp_between]
pvals_between = compdf$p.value[comp_between]
tvals_within_intro = compdf$t.ratio[comp_within_intro]
pvals_within_intro = compdf$p.value[comp_within_intro]
tvals_within_extro = compdf$t.ratio[comp_within_extro]
pvals_within_extro = compdf$p.value[comp_within_extro]


stars_between = vector('character',length = length(pvals_between))
stars_within_intro = vector('character', length = length(pvals_within_intro))
stars_within_extro = vector('character', length = length(pvals_within_extro))
stars_between[pvals_between < 0.05] = 'x'
stars_within_intro[pvals_within_intro < 0.05] = '*'
stars_within_extro[pvals_within_extro < 0.05] = '*'

# plot
focus_g = ggplot(alldf_summ, aes(x = task, y = attn.mean, group = extravert, shape = extravert, color = extravert)) +
  geom_jitter(data = alldf, mapping = aes(x = task, y = attention, group = extravert, shape = extravert), 
              alpha = 0.1, position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.4)) +
  geom_errorbar(aes(x = task, ymin = attn.mean - attn.se, ymax = attn.mean + attn.se),
                width = 0.3, position = position_dodge(width = 0.4)) +
  geom_point(position = position_dodge(width = 0.4), size = 2) + 
  scale_shape_manual(values = c(16,17), labels = c("Introverts", "Extraverts")) +
  scale_colour_manual(values = c("black", "grey"), labels = c("Introverts", "Extraverts")) + 
  scale_x_discrete(labels = c("HA Music", "LA Music", "No Music")) +
  scale_y_continuous(name = "Focus on Breath",
                     limits = c(-0.02,1.06), breaks = seq(0,1,.2)) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_text(angle = 90),
        legend.position = "top",
        legend.title = element_blank()) +
  annotate("text", x = 1, y = 1.02, label = stars_between[1], hjust = 0.5, vjust = 1, fontface = 'bold', size = fontsize-4, colour = 'black') +
  annotate("text", x = 2, y = 1.02, label = stars_between[2], hjust = 0.5, vjust = 1, fontface = 'bold', size = fontsize, colour = 'black') +
  annotate("text", x = 3, y = 1.02, label = stars_between[3], hjust = 0.5, vjust = 1, fontface = 'bold', size = fontsize, colour = 'black') +
  #annotate("text", x = (1.06+2.06)/2, y = 0.85, label = stars_within_extro[1], hjust = 0.5, vjust = 1, fontface = 'bold', size = fontsize, colour = 'red') +
  #annotate("segment", x = 1.06, xend = 2.06, y = 0.84, yend = 0.84, colour = 'red') +
  annotate("text", x = (2.10+3.10)/2, y = 0.84, label = stars_within_extro[2], hjust = 0.5, vjust = 0, fontface = 'bold', size = fontsize, colour = 'grey') +
  annotate("segment", x = 2.10, xend = 3.10, y = 0.83, yend = 0.83, colour = 'grey') + 
  annotate("text", x = (0.90+1.9)/2, y = 0.59, label = stars_within_intro[1], hjust = 0.5, vjust = 1, fontface = 'bold', size = fontsize, colour = 'black') +
  annotate("segment", x = 0.9, xend = 1.9, y = 0.6, yend = 0.6, colour = 'black') + 
  annotate("text", x = (1.9+2.9)/2, y = 0.57, label = stars_within_intro[2], hjust = 0.5, vjust = 1, fontface = 'bold', size = fontsize, colour = 'black') +
  annotate("segment", x = 1.9, xend = 2.9, y = 0.58, yend = 0.58, colour = 'black') +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
focus_g

ggsave('MnM_FocusOnTask_IntrovertsVsExtroverts_Task.png', focus_g, units = "cm", width = 8, height = 8, dpi = 300, device = 'png')
ggsave('MnM_FocusOnTask_IntrovertsVsExtroverts_Task.pdf', focus_g, units = "cm", width = 8, height = 8, dpi = 300, device = 'pdf')
