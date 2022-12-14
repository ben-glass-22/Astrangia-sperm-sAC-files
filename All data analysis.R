#Astrangia sperm sAC project data analysis and presentation

library(ggplot2)
library(plyr)
library(ggh4x)
library(MuMIn)
library(lme4)
library(car)
library(emmeans)
library(multcomp)
library(multcompView)

summarySE <- function(data = NULL,
                      measurevar,
                      groupvars = NULL,
                      na.rm=FALSE,
                      conf.interval=.95,
                      .drop=TRUE) {
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)
  
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

#bulk cAMP assays

bcadata <-
  read.csv(file.choose(),
           header = T,
           stringsAsFactors = T)

bcadatasum <- summarySE(bcadata,
                        measurevar = "cAMP_nmol_ng",
                        groupvars = c("Time",
                                      "Treatment",
                                      "Stim"))

bcaplot <-
  ggplot(data = bcadatasum,
         aes(x = Time,
             y = cAMP_nmol_ng,
             color = Stim)) +
  geom_point() +
  geom_line(aes(group = Stim)) +
  geom_errorbar(aes(ymin = cAMP_nmol_ng - se,
                    ymax = cAMP_nmol_ng + se),
                width = .5) +
  facet_wrap(~ Treatment) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = c(0.88, 0.76),
        legend.key = element_rect(fill = NA),
        axis.text.x = element_text(size = 8,
                                   color = "black"), 
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 8)) +
  xlab("Time (min)") +
  ylab(expression(cAMP~(nmol~ng~protein^{-1}))) +
  scale_y_continuous(limits = c(0, 100)) +
  scale_color_manual(values = c("No" = "#54acb9",
                                "Yes" = "#3f7c8b"))

bcaplot

ggsave("Bulk cAMP figure.png",
       width = 2.83,
       height = 2.38,
       units = c("in"))

bcadata$Time <- as.factor(bcadata$Time)

bcamodel <- lm(cAMP_nmol_mg ~ Time + Treatment*Stim,
                 data = bcadata)
bcaanova <- Anova(bcamodel,
                  type = "III")
bcaanova

#PKA substrate assays

PKAdata <-
  read.csv(file.choose(),
           header = T,
           stringsAsFactors = T)

PKAdatasum <- summarySE(PKAdata,
                        measurevar = "Fold_increase",
                        groupvars = c("Time",
                                      "Treatment"))

PKAplot <-
  ggplot(data = PKAdatasum,
         aes(x = Time,
             y = Fold_increase,
             color = Treatment)) +
  geom_point() +
  geom_line() +
  geom_errorbar(aes(ymin = Fold_increase - se,
                    ymax = Fold_increase + se),
                width = .2) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = c(0.85, 0.8),
        legend.key = element_rect(fill = NA),
        axis.text.x = element_text(size = 8,
                                   color = "black"), 
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.x = element_text(size = 8), 
        axis.title.y = element_text(size = 8),
        strip.text.x = element_text(size = 8)) +
  xlab("Time (min)") +
  ylab("PKA substrate phosphorylation (fold increase)") +
  scale_color_manual(values = c("DMSO" = "#54acb9",
                                "KH7" = "#826370",
                                "H-89" = "#e2876b"))

PKAplot

ggsave("PKA sub figure.png",
       width = 2.63,
       height = 2.58,
       units = c("in"))

PKAdata$Time <- as.factor(PKAdata$Time)

PKAmodel <- lm(Fold_increase ~ Time*Treatment,
                 data = PKAdata)
PKAanova <- Anova(PKAmodel,
                  type = "III")
PKAanova

#motility videos

motilitydata <- read.csv(file.choose(),
                         header = T,
                         stringsAsFactors = T)

motilitydatasum <- summarySE(motilitydata,
                        measurevar = "Prct_motile",
                        groupvars = c("Medium",
                                      "Treatment",
                                      "Stim"))

motilityplot <-
  ggplot(data = motilitydatasum,
         aes(x = factor(Treatment,
                        levels = c("None",
                                   "DMSO",
                                   "KH7",
                                   "H-89")),
             y = Prct_motile,
             fill = Stim)) +
  geom_bar(stat = "identity",
           position = position_dodge()) +
  geom_errorbar(aes(ymin = Prct_motile - se,
                    ymax = Prct_motile + se),
                position = position_dodge2()) +
  facet_wrap(~ Medium,
             scales = "free_x") +
  force_panelsizes(cols = c(.8, 0.2)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = c(0.9, 0.85),
        axis.text.x = element_text(size = 8,
                                   color = "black"), 
        axis.text.y = element_text(size = 8,
                                   color = "black"),
        axis.title.x = element_text(size = 10), 
        axis.title.y = element_text(size = 10),
        strip.text.x = element_text(size = 8)) +
  xlab("Treatment") +
  ylab("Motility (%)") +
  scale_fill_manual(values = c("No" = "#54acb9",
                               "Yes" = "#3f7c8b")) +
  ylim(0, 88)

motilityplot

ggsave("Motility figure.png",
       width = 4.5,
       height = 4,
       units = c("in"))

motilitydatanafsw <- subset(motilitydata,
                            Medium == "NaFSW")
motilitydatasw <- subset(motilitydata,
                         Medium == "SW")

nafswmod <- lm(Prct_motile ~ Treatment*Stim,
                data = motilitydatanafsw)
nafswanova <- Anova(nafswmod,
                    type = "III")
nafswanova
nafswtukey <- emmeans(nafswmod,
                     list(pairwise ~ Treatment*Stim),
                     adjust = "Tukey")
nafswtukey

swmod <- lm(Prct_motile ~ Stim,
               data = motilitydatasw)
swanova <- Anova(swmod,
                    type = "II")
swanova