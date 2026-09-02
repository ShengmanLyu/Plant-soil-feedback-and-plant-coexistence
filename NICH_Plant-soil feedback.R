#*********************************************************************************
# Plant-soil feedbacks pilot
# 28.01.2025
#*********************************************************************************
rm(list=ls())

library(tidyverse)
library(mgcv)
library(readxl)
library(lme4)
library(lmerTest)
library(car)
library(stringr)
library(viridis)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Biomass ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
d <- read_excel("/Users/slyu1/LVSM/NICH/Data/plant-soil/plant_soil feedback.xlsx", na="NA")  

# sample size
d %>%
  filter(comment != "dying") %>%
  filter(!is.na(biomass.mg)) %>%
  group_by(focal.species, background.species, site) %>%
  summarize(sample.size=n()) %>%
  filter(sample.size<3)

# Plot biomass two sites (high site)
d.none <- d %>%
  filter(background.species != "none")

d.az <- d %>%
  filter(background.species == "none") %>%
  mutate(site = "Anzeindaz")

d.sl <- d %>%
  filter(background.species == "none") %>%
  mutate(site = "Solalex")

bind_rows(d.none, d.az, d.sl)  %>%
  #filter(!is.na(site)) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.community=ifelse(background.species == "none", "Sterilized", Soil.community)) %>%
  mutate(focal.species=factor(focal.species, levels=c("Anal", "Asal", "Plal", "Melu", "Sapr", "Plla"))) %>%
  mutate(Soil.inoculum = str_replace(background.species, "none", "Sterilized")) %>%
  mutate(Soil.inoculum = factor(Soil.inoculum, levels=c("Sterilized", "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  ggplot(aes(y=biomass.mg, x=Soil.inoculum, col=Soil.community)) +
  #geom_point() +
  stat_summary(size=0.4) +
  #geom_boxplot(alpha=0) +
  facet_grid(rows=vars(Elevation), cols = vars(focal.species), scales="free_y") +
  #facet_wrap(~focal.species, nrow=1, scales="free_y") +
  scale_x_discrete(name="Soil inoculum") +
  scale_y_continuous(name="Plant biomass (mg)") +
  scale_color_manual(values = c("Sterilized" = col.viridis[3], "Interspecific" = col.viridis[2], "Intraspecific" = col.viridis[1])) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=60, hjust=1, size=9),
        legend.position = "null")

# test with full dataset
# "sterilized" treatment were excluded because no sites
d %>%
  filter(background.species != "sterilized") %>%
  lmer(biomass.mg ~ focal.species+ background.species + site + focal.species:background.species + background.species:site + (1|replicate) , data=.) %>%
  #summary()
  Anova()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### PSF for each species ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# calculate PSF
# Match microbial interactions
for(i in 1:nrow(d)) {
  # i = 1
  di <- d[i,]
  site.i <- di$site
  focal.i <- di$focal.species
  bg.i <- di$background.species
  biomass.i <- di$biomass.mg
  if(is.na(biomass.i) | bg.i == "none") next
  
  # PSF using sterilized soil as a reference
  biomass.none <- mean(filter(d, focal.species == focal.i & background.species == "none")$biomass.mg, na.rm=TRUE)
  d[i, "PSF.none"] <- log(biomass.i/biomass.none)
}

# plot PSF for two sites
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific inoculation", "Interspecific inoculation")) %>%
  mutate(focal.species=factor(focal.species, levels=c("Anal", "Asal", "Plal", "Melu", "Sapr", "Plla"))) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  ggplot(aes(y=PSF.none, x=Soil.inoculum, col=Soil.community)) +
  geom_hline(yintercept = 0, linetype="dashed") +
  #geom_point() +
  stat_summary(size=0.4) +
  #geom_boxplot(alpha=0) +
  facet_grid(rows=vars(Elevation), cols = vars(focal.species), scales="free_y") +
  scale_x_discrete(name="Soil inoculum") +
  scale_y_continuous(name="Microbial effects (m)") +
  scale_color_manual(values = c("Interspecific inoculation" = col.viridis[2], "Intraspecific inoculation" = col.viridis[1])) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=60, hjust=1, size=9.5),
        legend.position = "null")

# test with full model
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  # the full model
  lmer(PSF.none ~ focal.species * Soil.inoculum * Elevation + (1|replicate) , data=.) %>% 
  #lmer(PSF.none ~ focal.species + Soil.inoculum + Elevation + Elevation:Soil.inoculum + (1|replicate) , data=.) %>% 
  #summary()
  Anova()

# test with simplified model and soil community
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  # the full model
  #lmer(PSF.none ~ focal.species * Soil.inoculum * Elevation + Soil.community + (1|replicate) , data=.) %>% 
  lmer(PSF.none ~ focal.species + Soil.inoculum + Soil.community + Elevation + Elevation:Soil.inoculum + (1|replicate) , data=.) %>% 
  #summary()
  Anova()

# test with simplified model
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  # the full model
  #lmer(PSF.none ~ focal.species * Soil.inoculum * Elevation + Soil.community + (1|replicate) , data=.) %>% 
  lmer(PSF.none ~ focal.species + Soil.inoculum + Elevation + Elevation:Soil.inoculum + (1|replicate) , data=.) %>% 
  #summary()
  Anova()


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot by focal species
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  ggplot(aes(y=exp(PSF.none), x=focal.species)) +
  stat_summary(size=0.3) +
  scale_x_discrete(name="Focal species") +
  scale_y_continuous(name="PSF \n ln(inoculated/sterilized)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=60, hjust=1, size=10),
        legend.position = "null")

# plot by conditioning species and site
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  ggplot(aes(y=exp(PSF.none), x=Soil.inoculum, col=Elevation)) +
  stat_summary(size=0.3) +
  scale_x_discrete(name="Soil inoculum") +
  scale_y_continuous(name="PSF \n ln(inoculated/sterilized)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=60, hjust=1, size=10))

# plot by sites
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  ggplot(aes(y=exp(PSF.none), x=Elevation)) +
  stat_summary(size=0.3) +
  scale_x_discrete(name="Elevation") +
  scale_y_continuous(name="PSF \n ln(inoculated/sterilized)") +
  facet_wrap(~Soil.inoculum) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=60, hjust=1, size=10),
        legend.position = "null")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot by intra vs inter
d %>%
  filter(!is.na(PSF.none)) %>%
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Soil.community=ifelse(focal.species == background.species, "Intraspecific", "Interspecific")) %>%
  mutate(Soil.inoculum = factor(background.species, levels=c( "Anal", "Asal", "Melu", "Plal", "Plla", "Sapr"))) %>%
  ggplot(aes(y=exp(PSF.none), x=Soil.community)) +
  stat_summary(size=0.3) +
  scale_x_discrete(name="Soil community") +
  scale_y_continuous(name="PSF \n ln(inoculated/sterilized)") +
  facet_wrap(~Soil.inoculum) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle=60, hjust=1, size=10),
        legend.position = "null")

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### PSF-driven coexistence ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
outcome <- read_csv("/Users/slyu1/LVSM/NICH/Analysis/NICH-CH3/Coexistence.csv", na="NA")
outcome <- outcome %>% 
  mutate(Elevation = ifelse(site == "Anzeindaz", "High site", "Low site")) %>%
  mutate(Elevation = factor(Elevation, levels=c("Low site", "High site")))

# Match microbial interactions
for(i in 1:nrow(outcome)) {
  # i = 126
  print(i)
  oi <- outcome[i,]
  site.i <- oi$site
  sps1 <- oi$sps1
  sps2 <- oi$sps2
  
  # skip Les Posses
  if(site.i == "Les Posses") next
  # skip species not in the plant-soil experiment
  if(!(sps1 %in% unique(d$focal.species)) | !(sps2 %in% unique(d$focal.species))) next
  if(!(sps1 %in% unique(d$background.species)) | !(sps2 %in% unique(d$background.species))) next
  # match data
  #sps1.none <- mean(filter(d, focal.species == sps1 & background.species == "none")$biomass.mg, na.rm=TRUE)
  #sps2.none <- mean(filter(d, focal.species == sps2 & background.species == "none")$biomass.mg, na.rm=TRUE)
  
  # Pairwise PSF
  print(i)
  sps1.none <- mean(filter(d, focal.species == sps1 & background.species == "none")$biomass.mg, na.rm=TRUE)
  sps2.none <- mean(filter(d, focal.species == sps2 & background.species == "none")$biomass.mg, na.rm=TRUE)
  sps1.intra <- mean(filter(d, site == site.i & focal.species == sps1 & background.species == sps1)$biomass.mg, na.rm=TRUE)
  sps2.intra <- mean(filter(d, site == site.i & focal.species == sps2 & background.species == sps2)$biomass.mg, na.rm=TRUE)
  sps1.inter <-  mean(filter(d, site == site.i & focal.species == sps1 & background.species == sps2)$biomass.mg, na.rm=TRUE)
  sps2.inter <-  mean(filter(d, site == site.i & focal.species == sps2 & background.species == sps1)$biomass.mg, na.rm=TRUE)
  
  # coexistence components using log response ratio
  m11 <- log(sps1.intra/sps1.none)
  m12 <- log(sps1.inter/sps1.none)
  m22 <- log(sps2.intra/sps2.none)
  m21 <- log(sps2.inter/sps2.none)
  
  # coexistence components using response ratio
  #m11 <- sps1.intra/sps1.none
  #m12 <- sps1.inter/sps1.none
  #m22 <- sps2.intra/sps2.none
  #m21 <- sps2.inter/sps2.none
  
  # Coexistence based on PSF
  PSF.RFD <- 0.5*(m11 + m12) - 0.5*(m22 + m21)
  PSF.ND <- -0.5*(m11 + m22 - m12 - m21)
  #PSF.total1 <-(sps1.intra * sps1.inter)/(sps1.none^2)
  #PSF.total2 <-(sps2.intra * sps2.inter)/(sps2.none^2)
  #PSF.total1 <-log((sps1.intra + sps1.inter)/(2*sps1.none))
  #PSF.total2 <-log((sps2.intra + sps2.inter)/(2*2sps2.none))
  #PSF.RFD <- 1/2*log(PSF.total1/PSF.total2)
  #PSF.ND <- -1/2*(log(sps1.intra/sps1.inter) + log(sps2.intra/sps2.inter))
  
  # outcomes based on PSF
  if(abs(PSF.RFD) > abs(PSF.ND)) outcome.PSF = "Competitive exclusion"
  else if(PSF.ND > 0) outcome.PSF = "Coexistence"
  else if(PSF.ND < 0) outcome.PSF = "Priority effects"
  
  print(c(sps1.intra, sps1.inter, sps1.none))
  print(c(sps2.intra, sps2.inter, sps2.none))
  
  
  # Trait differences
  # height
  #trait.sps1 <- filter(trait, site == site.i & species == sps1)$reproductive.height
  #trait.sps2 <- filter(trait, site == site.i & species == sps2)$reproductive.height
  
  # light
  trait.sps1 <- filter(trait, site == site.i & species == sps1)$light.intercept.10cm
  trait.sps2 <- filter(trait, site == site.i & species == sps2)$light.intercept.10cm
  
  # vegetative.height
  #trait.sps1 <- filter(trait, site == site.i & species == sps1)$vegetative.height
  #trait.sps2 <- filter(trait, site == site.i & species == sps2)$vegetative.height
  
  if(length(trait.sps1) ==0 | length(trait.sps2) ==0) trait.dif = NA 
  else trait.dif = -trait.sps1 - (-trait.sps2)
  
  # combine
  outcome[i, "PSF.RFD"] <- PSF.RFD
  outcome[i, "PSF.ND"] <- PSF.ND
  outcome[i, "Trait.dif"] <- trait.dif
  outcome[i, "outcome.PSF"] <- outcome.PSF
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Coexistence in the field
outcome %>%
  mutate(FD2 = sqrt(sensitivity.21/sensitivity.12)) %>%
  ggplot(aes(x=FD2, y=fd)) +
  geom_point()

# outcome
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  ggplot(aes(x=-log(1-nd), y=log(fd), label=paste(sps1, sps2, sep=""))) + 
  # coexistence area
  geom_polygon(data=data.frame(x=c(0, 1.6, 1.6, 0),
                               y= c(0, 1.6, -1.6, 0)), aes(x=x, y=y), fill=col.viridis[2], alpha=0.4) +
  # priority area
  geom_polygon(data=data.frame(x=c(0, -0.6, -0.6, 0),
                               y= c(0, 0.6, -0.6, 0)), aes(x=x, y=y), fill=col.viridis[3], alpha=0.4) +
  geom_hline(yintercept = 0, size=0.1) +
  geom_vline(xintercept = 0, size=0.1) +
  geom_abline(intercept = 0, slope=1,  linetype="dashed") +
  geom_abline(intercept = 0, slope=-1, linetype="dashed") +
  geom_point(aes(shape=outcome.ndfd), size=2) + 
  #geom_text() +
  scale_x_continuous(name="Total niche differences in the field") +
  scale_y_continuous(name="Total fitness differences in the field") +
  facet_wrap(~Elevation) +
  coord_cartesian(xlim=c(-0.5, 1.5), ylim=c(-1.5, 1.5)) + 
  scale_shape_manual(values = c(16, 21)) +
  theme_bw() +
  theme(legend.position = "NULL",
        strip.text.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank()
  )

# RFD ~ site
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  lmer(log(fd) ~ Elevation +(1|sps1) + (1|sps2), data=.) %>%
  Anova()

# ND ~ site
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  lmer(-log(1-nd) ~ site +(1|sps1) + (1|sps2), data=.) %>%
  Anova()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PSF-driven coexistence
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  ggplot(aes(x=PSF.ND, y=PSF.RFD, shape=outcome.PSF, label=paste(sps1, sps2, sep=""))) + 
  # coexistence area
  geom_polygon(data=data.frame(x=c(0, 1.6, 1.6, 0),
                               y= c(0, 1.6, -1.6, 0)), aes(x=x, y=y), fill=col.viridis[2], alpha=0.4) +
  # priority area
  geom_polygon(data=data.frame(x=c(0, -0.6, -0.6, 0),
                               y= c(0, 0.6, -0.6, 0)), aes(x=x, y=y), fill=col.viridis[3], alpha=0.4) +
  geom_hline(yintercept = 0, size=0.1) +
  geom_vline(xintercept = 0, size=0.1) +
  geom_abline(intercept = 0, slope=1, linetype="dashed") +
  geom_abline(intercept = 0, slope=-1, linetype="dashed") +
  geom_point(size=2) + 
  #geom_text(size=3) +
  scale_x_continuous(name="PSF-driven niche differences") +
  scale_y_continuous(name="PSF-driven fitness differences") +
  coord_cartesian(xlim=c(-0.5, 1.5), ylim=c(-1.5, 1.5)) + 
  scale_shape_manual(values = c(16, 21)) +
  facet_wrap(~Elevation) +
  theme_bw() +
  theme(legend.position = "NULL",
        strip.text = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank()
  )

# RFD ~ site
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  lmer(abs(PSF.RFD) ~ Elevation +(1|sps1) + (1|sps2), data=.) %>%
  Anova()

# ND ~ site
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  lmer(PSF.ND ~ site +(1|sps1) + (1|sps2), data=.) %>%
  Anova()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Comparison in outcomes
# count outcomes
outcome %>%
  filter(!is.na(outcome.PSF)) %>%
  group_by(site, outcome.PSF) %>%
  summarize(n=n()) 

outcome %>%
  filter(!is.na(outcome.PSF)) %>%
  group_by(site, outcome.ndfd) %>%
  summarize(n=n()) 

# Chisq test: differed significantly
# high site
chisq.test(matrix(c(5, 4, 9, 0), ncol = 2))
# low site
chisq.test(matrix(c(8, 1, 1, 8), ncol = 2))

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PSF ~ RFD
# PSF differences
outcome %>%
  filter(site != "Les Posses") %>%
  filter(!is.na(PSF.RFD)) %>%
  ggplot(aes(x=PSF.RFD, y=log(fd), col=Elevation)) +
  geom_hline(yintercept = 0, linetype="dashed", size=0.2) + 
  geom_vline(xintercept = 0, linetype="dashed", size=0.2) +
  geom_point(aes(shape=outcome.ndfd), size=2) +
  geom_smooth(data=filter(outcome, Elevation=="High site"), method = lm) + 
  geom_smooth(data=filter(outcome, Elevation=="Low site"), method = lm, linetype="dashed") + 
  #geom_text(aes(label = pair)) +
  scale_x_continuous(name="PSF-driven fitness differences") +
  scale_y_continuous(name="Total fitness differences in the field") +
  scale_color_manual(values = c("High site" = "blue", "Low site" = "orange")) +
  scale_shape_manual(values=c(16,21)) +
  facet_wrap(~Elevation) +
  theme_bw() +
  theme(legend.position = "NULL",
        strip.text.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank())

# test
outcome %>%
  filter(site != "Les Posses") %>%
  filter(!is.na(PSF.RFD)) %>%
  #filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(log(fd) ~ PSF.RFD*site + (1|sps1) + (1|sps2), data=.) %>%
  Anova()

# test in high site
outcome %>%
  filter(site == "Anzeindaz") %>%
  filter(!is.na(PSF.RFD)) %>%
  #filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(log(fd) ~ PSF.RFD + (1|sps1) + (1|sps2), data=.) %>%
  #summary()
  Anova()

# test in low site
outcome %>%
  filter(site == "Solalex") %>%
  filter(!is.na(PSF.RFD)) %>%
  #filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(log(fd) ~ PSF.RFD + (1|sps1) + (1|sps2), data=.) %>%
  #summary()
  Anova()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PSF ~ ND
# Plot
outcome %>%
  filter(site != "Les Posses") %>%
  filter(!is.na(PSF.ND)) %>%
  ggplot(aes(x=PSF.ND, y=-log(1-nd), col=Elevation)) +
  geom_hline(yintercept = 0, linetype="dashed", size=0.2) + 
  geom_vline(xintercept = 0, linetype="dashed", size=0.2) +
  geom_point(aes(shape=outcome.ndfd), size=2) +
  #geom_text(aes(label = pair)) +
  geom_smooth(data=filter(outcome, Elevation=="Low site"), method = lm) + 
  geom_smooth(data=filter(outcome, Elevation=="High site"), method = lm, linetype="dashed") + 
  scale_x_continuous(name="PSF-driven niche differences") +
  scale_y_continuous(name="Total niche differences in the field") +
  scale_color_manual(values = c("High site" = "blue", "Low site" = "orange")) +
  scale_shape_manual(values=c(16,21)) +
  facet_wrap(~Elevation) +
  theme_bw() +
  theme(legend.position = "NULL",
        strip.text.x = element_text(size = 10),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank())

# test
outcome %>%
  filter(site != "Les Posses") %>%
  filter(!is.na(PSF.ND)) %>%
  #filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(-log(1-nd) ~ PSF.ND*site + (1|sps1) + (1|sps2), data=.) %>%
  Anova()

# test
outcome %>%
  filter(site == "Solalex") %>%
  filter(!is.na(PSF.ND)) %>%
  #filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  lmer(-log(1-nd) ~ PSF.ND + (1|sps1) + (1|sps2), data=.) %>%
  Anova()

# test
outcome %>%
  filter(site == "Anzeindaz") %>%
  filter(!is.na(PSF.ND)) %>%
  #filter(sensitivity.12 != 0.1 & sensitivity.21 != 0.1) %>%
  #lm(-log(1-nd) ~ PSF.ND, data=.) %>%
  lmer(-log(1-nd) ~ PSF.ND + (1|sps1) + (1|sps2), data=.) %>%
  Anova()

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
##### Figure intro ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PSF interacts with other processes jointly to shape coexistence

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Competition-driven coexistence
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  ggplot(aes(x=PSF.ND, y=PSF.RFD, shape=outcome.PSF, label=paste(sps1, sps2, sep=""))) + 
  # coexistence area
  geom_polygon(data=data.frame(x=c(0, 1.6, 1.6, 0),
                               y= c(0, 1.6, -1.6, 0)), aes(x=x, y=y), fill=col.viridis[2], alpha=0.4) +
  # priority area
  geom_polygon(data=data.frame(x=c(0, -1.1, -1.1, 0),
                               y= c(0, 1.1, -1.1, 0)), aes(x=x, y=y), fill=col.viridis[3], alpha=0.4) +
  geom_hline(yintercept = 0, size=0.1) +
  geom_vline(xintercept = 0, size=0.1) +
  geom_abline(intercept = 0, slope=1, linetype="dashed") +
  geom_abline(intercept = 0, slope=-1, linetype="dashed") +
  #geom_point(size=2) + 
  #geom_text(size=3) +
  scale_x_continuous(name="Competition-driven niche differences") +
  scale_y_continuous(name="Competition-driven fitness differences") +
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1, 1)) + 
  scale_shape_manual(values = c(16, 21)) +
  theme_bw() +
  theme(legend.position = "NULL",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank()
  )

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# PSF-driven coexistence
outcome %>% 
  filter(site != "Les Posses") %>%
  filter(!is.na(Trait.dif)) %>%
  ggplot(aes(x=PSF.ND, y=PSF.RFD, shape=outcome.PSF, label=paste(sps1, sps2, sep=""))) + 
  # coexistence area
  geom_polygon(data=data.frame(x=c(0, 1.5, 1.5, 0),
                               y= c(0, 1.5, -1.5, 0)), aes(x=x, y=y), fill=col.viridis[2], alpha=0.4) +
  # priority area
  geom_polygon(data=data.frame(x=c(0, -1.1, -1.1, 0),
                               y= c(0, 1.1, -1.1, 0)), aes(x=x, y=y), fill=col.viridis[3], alpha=0.4) +
  geom_hline(yintercept = 0, size=0.1) +
  geom_vline(xintercept = 0, size=0.1) +
  geom_abline(intercept = 0, slope=1, linetype="dashed") +
  geom_abline(intercept = 0, slope=-1, linetype="dashed") +
  #geom_point(size=2) + 
  #geom_text(size=3) +
  scale_x_continuous(name="PSF-driven niche differences") +
  scale_y_continuous(name="PSF-driven fitness differences") +
  coord_cartesian(xlim=c(-1, 1), ylim=c(-1, 1)) + 
  scale_shape_manual(values = c(16, 21)) +
  theme_bw() +
  theme(legend.position = "NULL",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.background = element_blank()
  )
