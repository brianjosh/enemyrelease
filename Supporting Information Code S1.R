#Supporting Information Code S1

#Code to support "Targeted aboveground enemy release increases seedling survival in grasslands"
#Code prepared by:
#Joshua Brian (Corresponding: joshua.brian@kcl.ac.uk; jshbrian@gmail.com)
#Maria Angeles Perez Navarro
#Harry Shepherd

###########################################################################################

#Contents of script:

#Part 1: Libraries and notes   --- Line 24
#Part 2: Damage analysis       --- Line 61
#Part 3: Phylogenetic analysis --- Line 159
#Part 4: Traits                --- Line 319
#Part 5: Survival analysis     --- Line 414
#Part 6: Main figures          --- Line 855
#Part 7: Supplementary figures --- Line 1047
#A joke to finish              --- Line 1186

###########################################################################################

#Part 1: Libraries and notes

#Read in required libraries
library(lme4)
library(glmmTMB)
library(car)
library(MuMIn)
library(tidyverse)
library(DHARMa)
library(pdp)
library(visreg)
library(rdacca.hp)
library(ggpubr)
library(cowplot)
library(readr)
library(tidyr)
library(stringr)
library(radiant.data)
library(taxize)
library(Taxonstand)
library(phytools)
library(V.PhyloMaker)
library(dendextend)
library(phylogram)
library(circlize)
library(phylosignal)

#IMPORTANT NOTE: Species in this code and in data files are referred to by a 
#five letter abbreviation (e.g. Achillea millefolium = Achmi). However, 
#Cedar Creek has its own name for species based on earlier
#taxonomy. Two species (S. oolentangiense and C. gracile) are known as 
#Aster azureus and Bouteloua curtipendula, hence abbreviations in this code and 
#input files are Astaz and Bougr. These were updated to currently recognised taxonomy
#for figures and in the final paper.

###########################################################################################

#Part 2: Damage analysis

#Descriptive statistics
damageE <- read.csv("damage_E.csv", header=T, stringsAsFactors = T)
damageE$plot <- as.factor(damageE$plot)
hist(damageE$meaninsect)
hist(damageE$meanfungal)
damageM <- read.csv("damage_M.csv", header=T, stringsAsFactors = T)
damageL$plot <- as.factor(damageM$plot)
hist(damageM$meaninsect)
hist(damageM$meanfungal)

aggregate(damageE$meaninsect, list(damageE$treatment), mean)
EinsectSD <- aggregate(damageE$meaninsect, list(damageE$treatment), sd)
(EinsectSE <- EinsectSD$x/sqrt(length(damageE$meaninsect)))

aggregate(damageE$meanfungal, list(damageE$treatment), mean)
EfungalSD <- aggregate(damageE$meanfungal, list(damageE$treatment), sd)
(EfungalSE <- EfungalSD$x/sqrt(length(damageE$meanfungal)))

aggregate(damageM$meaninsect, list(damageM$treatment), mean)
MinsectSD <- aggregate(damageM$meaninsect, list(damageM$treatment), sd)
(MinsectSE <- MinsectSD$x/sqrt(length(damageM$meaninsect)))

aggregate(damageM$meanfungal, list(damageM$treatment), mean)
MfungalSD <- aggregate(damageM$meanfungal, list(damageM$treatment), sd)
(MfungalSE <- MfungalSD$x/sqrt(length(damageM$meanfungal)))

#Pairwise Wilcoxon tests comparing different treatments and sites

overall <- read.csv("damage_overall.csv", header=T, stringsAsFactors = T)
aggregate(overall$damage, list(overall$combo), mean)
overall <- split(overall, overall$type)
overallinsect <- overall$ins
overallinsect$combo <- factor(overallinsect$combo)
overallinsect$plot <- as.factor(overallinsect$plot)
overallfungi <- overall$fun
overallfungi$combo <- factor(overallfungi$combo)
overallfungi$plot <- as.factor(overallfungi$plot)

kruskal.test(damage ~ treatment, data=overallinsect)
pairwise.wilcox.test(overallinsect$damage, overallinsect$combo, p.adjust.method = "BH")

kruskal.test(damage ~ treatment, data=overallfungi)
pairwise.wilcox.test(overallfungi$damage, overallfungi$combo, p.adjust.method = "BH")

#Confirm that treatments differ, and test effect of richness and site,
#using zero-inflated generalised linear models

insectmodel <- glmmTMB(damage ~ treatment + stage + richness +
                         (1|species/plot),
                       family=poisson, ziformula=~1,
                       data=overallinsect)
summary(insectmodel)
Anova(insectmodel, type="III")
simulationOutput <- simulateResiduals(fittedModel = insectmodel, plot=F)
plot(simulationOutput)

fungimodel <- glmmTMB(damage ~ treatment + stage + richness +
                        (1|species/plot),
                      family=poisson, ziformula=~1,
                      data=overallfungi)
summary(fungimodel)
Anova(fungimodel, type="III")

#Compare effect sizes per species to damage from control only

damageeffect <- read.csv("effect_size_species.csv", header=T, stringsAsFactors = T)

treatmentdamage <- split(damageeffect, damageeffect$treatment)
insecticide <- treatmentdamage$I
fungicide <- treatmentdamage$IF

#Statistically test these

insectmod <- lm(ef ~ meaninsectdamageC, data=insecticide)
summary(insectmod)
#without outlier
insecticideadj <- insecticide[-12, ]
insectmodadj <- lm(ef ~ meaninsectdamageC, data=insecticideadj)
summary(insectmodadj)

fungalmod <- lm(ef ~ meaninsectfungaldamageC, data=fungicide)
summary(fungalmod)
#without outlier
fungicideadj <- fungicide[-10, ]
fungalmodadj <- lm(ef ~ meaninsectfungaldamageC, data=fungicideadj)
summary(fungalmodadj)

fungalmod <- lm(ef ~ meanfungaldamageC, data=fungicide)
summary(fungalmod)
#without outlier
fungicideadj <- fungicide[-10, ]
fungalmodadj <- lm(ef ~ meanfungaldamageC, data=fungicideadj)
summary(fungalmodadj)

###########################################################################################

#Part 3: Phylogenetic analysis

#Read in data file containing percent abundances for every species in every plot
abun_traits <- read.csv("abundances.csv", header=T)
str(abun_traits)

#get species names
species <- unique(abun_traits$species)

sp_file <- abun_traits%>%
  select(species, genus, family, species.relative, genus.relative)%>%
  distinct() #remove the duplicate names

#make the tree
rel <- bind.relative(sp.list=sp_file, tree=GBOTB.extended, nodes=nodes.info.1)
whole_tree <- phylo.maker(rel$species.list, tree=rel$phylo, nodes=rel$nodes.info,
                          scenarios=c("S1","S2","S3"))
sce1_tree <- whole_tree$scenario.1
par(mfrow = c(1, 1))
plot.phylo(as.phylo(sce1_tree), cex = 0.6, main = "scenario.1",
           type="fan")

#Calculate abundance-weighted mean phylogenetic distance (awMPD) 
#for each plant species in each community (plot)
plot_df_list <- list()
table <- abun_traits
plot_list <- unique(table$plot)
for(i in 1:length(plot_list)) {
  table_plot <- table%>%
    filter(plot== plot_list[i])
  
  sp_list <- unique(table_plot$species)
  species_list0 <- unique(table_plot$species)
  species_list <- unique(table_plot$specieslowbar)
  sp_df_list <- list()
  
  cut_tree<-drop.tip(sce1_tree,
                     setdiff(sce1_tree$tip.label,species_list))
  phylo_dist <- cophenetic.phylo(cut_tree)%>%
    as.data.frame()%>%
    rownames_to_column()
  
  for(j in 1:length(sp_list)){
    
    table_plot_sp <- table_plot%>%
      filter(species== sp_list[j])
    
    table_plot_rest <- table_plot%>%
      filter(species!= sp_list[j])
    
    dist_x <- phylo_dist%>%
      filter(rowname==species_list[j])%>%
      select(-rowname)%>%
      t()%>%
      data.frame()%>%
      rownames_to_column()%>%
      filter(rowname!=species_list[j])%>%
      mutate(species=gsub("_", " ", rowname))
    
    names(dist_x)[2] <- "dist"
    
    dist_x2 <- table_plot_rest%>%
      left_join(dist_x, by="species")%>%
      select(-rowname)
    
    dist_x2%>%
      select(plot,
             species,
             dist,
             abundance)
    
    phylo_w_dist <- weighted.mean(dist_x2$dist,
                                  dist_x2$abundance)
    
    table_plot_sp <- table_plot_sp%>%
      mutate(phylo_w_dist=phylo_w_dist)
    
    sp_df_list[[j]] <- table_plot_sp
    
  }
  
  table_plot2 <- bind_rows(sp_df_list)
  plot_df_list[[i]] <- table_plot2
}

phylo_dist_all <- bind_rows(plot_df_list)
#The relevant awMPD value for the focal (planted) species in each plot was included  
#in a column to the 'survival_data.csv' file (the 'awMPD' column)

#Now make the tree for focal 16 species and test if effect sizes have 
#phylogenetic signal

effectsize <- read.csv("phylo_effect.csv", header=T)
str(effectsize)

sp_file <- effectsize%>%
  select(species, genus, family, species.relative, genus.relative)%>%
  distinct() #remove the duplicate names

#make the tree
rel <- bind.relative(sp.list=sp_file, tree=GBOTB.extended, nodes=nodes.info.1)
whole_tree <- phylo.maker(rel$species.list, tree=rel$phylo, nodes=rel$nodes.info,
                          scenarios=c("S1","S2","S3"))
sce1_tree <- whole_tree$scenario.1
write.tree(sce1_tree, file="species.tre")

par(mfrow = c(1, 1))
plot.phylo(as.phylo(sce1_tree), cex = 0.6, main = "scenario.1",
           type="fan")

#Define the tip data

trait_file <- effectsize%>%
  select(specieslowbar, effectsizeI, effectsizeIF)
rownames(trait_file) <- trait_file[ ,1]
trait_file[ ,1] <- NULL
write.csv(trait_file, file="traits.csv")

#Link tree with trait data

phylofile <- read.p4d("species.tre", "traits.csv", phylo.format = "newick", 
                      data.format="csv")
#Run the test

phyloSignal(phylofile)

#Now repeat while removing the two 'outliers' (P. virgatum and M. fistulosa)

effectsize2 <- filter(effectsize, specieslowbar != "Panicum_virgatum") %>%
  filter(specieslowbar != "Monarda_fistulosa")
sp_file2 <- effectsize2%>%
  select(species, genus, family, species.relative, genus.relative)%>%
  distinct() #remove the duplicate names

#make the tree
rel2 <- bind.relative(sp.list=sp_file2, tree=GBOTB.extended, nodes=nodes.info.1)
whole_tree2 <- phylo.maker(rel2$species.list, tree=rel$phylo, nodes=rel$nodes.info,
                           scenarios=c("S1","S2","S3"))
sce1_tree2 <- whole_tree2$scenario.1
write.tree(sce1_tree2, file="species2.tre")
par(mfrow = c(1, 1))
plot.phylo(as.phylo(sce1_tree2), cex = 0.6, main = "scenario.1",
           type="fan")

#Define the tip data
trait_file2 <- effectsize2%>%
  select(specieslowbar, effectsizeI, effectsizeIF)
rownames(trait_file2) <- trait_file2[ ,1]
trait_file2[ ,1] <- NULL
write.csv(trait_file2, file="traits2.csv")

#Link tree with trait data
phylofile <- read.p4d("species2.tre", "traits2.csv", phylo.format = "newick", 
                      data.format="csv")

#Run the test
phyloSignal(phylofile)

###########################################################################################

#Part 4: Traits

#Read in and prepare data
traitspace <- read.csv("species_traits_all.csv", header=T)
traitspace2 <- traitspace[,-1]
rownames(traitspace2) <- traitspace[,1]

#Aboveground traits
#limit to SLA, LDMC, leaf area, leaf N, leaf P
#other traits are available in dataset (leaf dry mass, leaf fresh mass, height, seed mass - 
#feel free to explore!)

#Run PCAs
traits.pca <- prcomp(traitspace2[ ,c(6:10)], scale. = T)
summary(traits.pca)
pcaresult <- traits.pca$x
write.csv(pcaresult, file="abovegroundpcaresults.csv")

#Belowground traits

#Run PCAs
root.pca <- prcomp(traitspace2[ ,c(2:5)], scale. = T)
summary(root.pca)
rootpcaresult <- root.pca$x
write.csv(rootpcaresult, file="rootpcaresults.csv")

#Modified versions of 'pcaresult' and 'rootpcaresult' were then combined with
#effect sizes to produce 'pca_results_modified.csv', which we will now read in

pcaplotting <- read.csv("pcaresultsmodified.csv", header=T, stringsAsFactors = T)

#See if there is a relationship between aboveground traits and effect sizes
#Include random effect of functional group (grass or forb)

#Effect of I vs. PC1aboveground
mod1 <- lmer(logeffectsizeI ~ PC1aboveground + (1|funcgroup), data=pcaplotting)
Anova(mod1, type="III")
summary(mod1)

#Effect of I vs. PC2aboveground
mod2 <- lmer(logeffectsizeI ~ PC2aboveground + (1|funcgroup), data=pcaplotting)
Anova(mod2, type="III")
summary(mod2)

#Effect of IF vs. PC1aboveground
mod3 <- lmer(logeffectsizeIF ~ PC1aboveground + (1|funcgroup), data=pcaplotting)
Anova(mod3, type="III")
summary(mod3)

#Effect of IF vs. PC2aboveground
mod4 <- lmer(logeffectsizeIF ~ PC2aboveground + (1|funcgroup), data=pcaplotting)
Anova(mod4, type="III")
summary(mod4)

#See if there is a relationship between root traits and effect sizes

#Effect of I vs. PC1root
mod5 <- lmer(logeffectsizeI ~ PC1root + (1|funcgroup), data=pcaplotting)
Anova(mod5, type="III")
summary(mod5)

#Effect of I vs. PC2root
mod6 <- lmer(logeffectsizeI ~ PC2root + (1|funcgroup), data=pcaplotting)
Anova(mod6, type="III")
summary(mod6)

#Effect of IF vs. PC1root
mod7 <- lmer(logeffectsizeIF ~ PC1root + (1|funcgroup), data=pcaplotting)
Anova(mod7, type="III")
summary(mod7)

#Effect of IF vs. PC2root
mod8 <- lmer(logeffectsizeIF ~ PC2root + (1|funcgroup), data=pcaplotting)
Anova(mod8, type="III")
summary(mod8)

#See if there is a relationship for SLA or leaf N separately 
mod9 <- lmer(logeffectsizeI ~ SLA + (1|funcgroup), data=pcaplotting)
Anova(mod9, type="III")
summary(mod9)

mod10 <- lmer(logeffectsizeI ~ N + (1|funcgroup), data=pcaplotting)
Anova(mod10, type="III")
summary(mod10)

mod11 <- lmer(logeffectsizeIF ~ SLA + (1|funcgroup), data=pcaplotting)
Anova(mod11, type="III")
summary(mod11)

mod12 <- lmer(logeffectsizeIF ~ N + (1|funcgroup), data=pcaplotting)
Anova(mod12, type="III")
summary(mod12)

###########################################################################################

#Part 5: Survival analysis 

survival <- read.csv("survival_data.csv", header=T, stringsAsFactors=T)
survival$plotnumber <- as.factor(survival$plotnumber)
survival$date <- as.factor(survival$date)

#scale continuous variables
survival <- survival %>% 
  mutate(across(c(awMPD, lightavailable, moisture, heightinitial), scale)) %>%
  mutate(across(c(awMPD, lightavailable, moisture, heightinitial), as.numeric))

#add unique ID for each plant
survival$group <- paste(survival$plot, survival$plantnumber, sep="_")
str(survival)
survival$group <- as.factor(survival$group)

#Split by sample date
time <- split(survival, survival$date)
June22nd <- time$`1`
June30th <- time$`2`
July6th <- time$`3`
July18th <- time$`4`
July25th <- time$`5`
July31st <- time$`6`
August8th <- time$`7`

#Split by species
allspecies <- split(survival, survival$species)
Achmi <- allspecies$Achmi
Agafo <- allspecies$Agafo
Andge <- allspecies$Andge
Ascsy <- allspecies$Ascsy
Asctu <- allspecies$Asctu
Astaz <- allspecies$Astaz
Boucu <- allspecies$Boucu
Bougr <- allspecies$Bougr
Liaas <- allspecies$Liaas
Monfi <- allspecies$Monfi
Muhra <- allspecies$Muhra
Panvi <- allspecies$Panvi
Rudhi <- allspecies$Rudhi
Schsc <- allspecies$Schsc
Sornu <- allspecies$Sornu
Spohe <- allspecies$Spohe

#Analyse proportions data
#Model set 1a

#Make new dataframe with the proportion of each plot that survived as the 
#response variable (column = 'mean')
propdata <- list()
for(i in 1:288){
  testJn22 <- filter(survival, plotnumber==i, date=="1")
  testJn30 <- filter(survival, plotnumber==i, date=="2")
  testJl6 <- filter(survival, plotnumber==i, date=="3")
  testJl18 <- filter(survival, plotnumber==i, date=="4")
  testJl25 <- filter(survival, plotnumber==i, date=="5")
  testJl31 <- filter(survival, plotnumber==i, date=="6")
  testAu8 <- filter(survival, plotnumber==i, date=="7")
  meansurvivedJn22 <- mean(testJn22$survived)
  lineJn22 <- mutate(testJn22, mean=meansurvivedJn22) %>%
    subset(select=-c(survived, plantnumber, heightinitial, group)) %>%
    distinct() 
  propdata[[i]] <- lineJn22
  meansurvivedJn30 <- mean(testJn30$survived)
  lineJn30 <- mutate(testJn30, mean=meansurvivedJn30) %>%
    subset(select=-c(survived, plantnumber, heightinitial, group)) %>%
    distinct() 
  propdata[[i+288]] <- lineJn30
  meansurvivedJl6 <- mean(testJl6$survived)
  lineJl6 <- mutate(testJl6, mean=meansurvivedJl6) %>%
    subset(select=-c(survived, plantnumber, heightinitial, group)) %>%
    distinct() 
  propdata[[i+576]] <- lineJl6
  meansurvivedJl18 <- mean(testJl18$survived)
  lineJl18 <- mutate(testJl18, mean=meansurvivedJl18) %>%
    subset(select=-c(survived, plantnumber, heightinitial, group)) %>%
    distinct() 
  propdata[[i+864]] <- lineJl18
  meansurvivedJl25 <- mean(testJl25$survived)
  lineJl25 <- mutate(testJl25, mean=meansurvivedJl25) %>%
    subset(select=-c(survived, plantnumber, heightinitial, group)) %>%
    distinct() 
  propdata[[i+1152]] <- lineJl25
  meansurvivedJl31 <- mean(testJl31$survived)
  lineJl31 <- mutate(testJl31, mean=meansurvivedJl31) %>%
    subset(select=-c(survived, plantnumber, heightinitial, group)) %>%
    distinct() 
  propdata[[i+1440]] <- lineJl31
  meansurvivedAu8 <- mean(testAu8$survived)
  lineAu8 <- mutate(testAu8, mean=meansurvivedAu8) %>%
    subset(select=-c(survived, plantnumber, heightinitial, group)) %>%
    distinct() 
  propdata[[i+1728]] <- lineAu8
}
propdataall <- bind_rows(propdata)

#add plot richness as a variable
richness <- read.csv("plot_richness.csv", header=T)
richness$plotnumber <- as.factor(richness$plotnumber)
str(richness)
repeats <- rep(richness$richness, 7)

propdataall2 <- bind_cols(propdataall, repeats)
propdataall2 <- dplyr::rename(propdataall2, richness = ...11)
propdataall2$richness <- scale(propdataall2$richness)
propdataall2$richness <- as.numeric(propdataall2$richness)

proptime2 <- split(propdataall2, propdataall2$date)
June22ndprop2 <- proptime2$`1`
June30thprop2 <- proptime2$`2`
July6thprop2 <- proptime2$`3`
July18thprop2 <- proptime2$`4`
July25thprop2 <- proptime2$`5`
July31stprop2 <- proptime2$`6`
August8thprop2 <- proptime2$`7`

#Analyse time period by time period
#(Interactions tested and non-significant)

Ju22prop2 <- glm(mean ~ treatment + species + awMPD +
                   stage + lightavailable + moisture + richness,
                 family=quasibinomial(link="logit"), 
                 data=June22ndprop2)
summary(Ju22prop2)
Anova(Ju22prop2, type="III")
r.squaredGLMM(Ju22prop2)
confint(Ju22prop2)

Ju30prop2 <- glm(mean ~ treatment + species + awMPD +
                   stage + lightavailable + moisture + richness,
                 family=quasibinomial(link="logit"), 
                 data=June30thprop2)
summary(Ju30prop2)
Anova(Ju30prop2, type="III")
r.squaredGLMM(Ju30prop2)
confint(Ju30prop2)

Jl6prop2 <- glm(mean ~ treatment + species + awMPD +
                  stage + lightavailable + moisture + richness,
                family=quasibinomial(link="logit"), 
                data=July6thprop2)
summary(Jl6prop2)
Anova(Jl6prop2, type="III")
r.squaredGLMM(Jl6prop2)
confint(Jl6prop2)

Jl18prop2 <- glm(mean ~ treatment + species + awMPD +
                   stage + lightavailable + moisture + richness,
                 family=quasibinomial(link="logit"), 
                 data=July18thprop2)
summary(Jl18prop2)
Anova(Jl18prop2, type="III")
r.squaredGLMM(Jl18prop2)
confint(Jl18prop2)

Jl25prop2 <- glm(mean ~ treatment + species + awMPD +
                   stage + lightavailable + moisture + richness,
                 family=quasibinomial(link="logit"), 
                 data=July25thprop2)
summary(Jl25prop2)
Anova(Jl25prop2, type="III")
r.squaredGLMM(Jl25prop2)
confint(Jl25prop2)

Jl31prop2 <- glm(mean ~ treatment + species + awMPD +
                   stage + lightavailable + moisture + richness,
                 family=quasibinomial(link="logit"), 
                 data=July31stprop2)
summary(Jl31prop2)
Anova(Jl31prop2, type="III")
r.squaredGLMM(Jl31prop2)
confint(Jl31prop2)

Au8prop2 <- glm(mean ~ treatment + species + awMPD +
                  stage + lightavailable + moisture + richness,
                family=quasibinomial(link="logit"), 
                data=August8thprop2)
summary(Au8prop2)
Anova(Au8prop2, type="III")
r.squaredGLMM(Au8prop2)
confint(Au8prop2)

#Poisson model of new dead
#Model set 1b

#Read in data containing new dead
#Note that this is a modified version of the proportions data file, all variables
#have already been scaled as per modifications to the 'survival_data' file
dead <- read.csv("new_dead.csv", header=T, stringsAsFactors = T)

newdeadmodel <- glm(newdead ~ treatment*date + species + awMPD +
                      stage + lightavailable + moisture + richness,
                    family=poisson, 
                    data=dead)
summary(newdeadmodel)
Anova(newdeadmodel, type="III")

1 - pchisq(summary(newdeadmodel)$deviance, summary(newdeadmodel)$df.residual)
#model fits the data well

#Analyse individual survival at each time period
#Model set 2.
#(Interactions tested and all non-significant, inclusion of interactions 
#leads to less power to detect main effects of treatment and stage)

#June 22nd
fullsurvivalmodelJn22 <- glmmTMB(survived ~ treatment + stage + heightinitial +
                                  (1|species/plotnumber), 
                                family=binomial(link="logit"), 
                                data=June22nd)
summary(fullsurvivalmodelJn22)
Anova(fullsurvivalmodelJn22, type="III")
confint(fullsurvivalmodelJn22)

#diagnostics confirmed to be suitable for this and all models
simulationOutput <- simulateResiduals(fittedModel = fullsurvivalmodelJn22, plot=F)
plot(simulationOutput)

#June 30th
fullsurvivalmodelJn30 <- glmmTMB(survived ~ treatment + stage + heightinitial +
                                   (1|species/plotnumber), 
                                 family=binomial(link="logit"), 
                                 data=June30th)
summary(fullsurvivalmodelJn30)
Anova(fullsurvivalmodelJn30, type="III")
confint(fullsurvivalmodelJn30)

#July 6th
fullsurvivalmodelJl6 <- glmmTMB(survived ~ treatment + stage + heightinitial +
                                   (1|species/plotnumber), 
                                 family=binomial(link="logit"), 
                                 data=July6th)
summary(fullsurvivalmodelJl6)
Anova(fullsurvivalmodelJl6, type="III")
confint(fullsurvivalmodelJl6)

#July 18th
fullsurvivalmodelJl18 <- glmmTMB(survived ~ treatment + stage + heightinitial +
                                  (1|species/plotnumber), 
                                family=binomial(link="logit"), 
                                data=July18th)
summary(fullsurvivalmodelJl18)
Anova(fullsurvivalmodelJl18, type="III")
confint(fullsurvivalmodelJl18)

#July 25th
fullsurvivalmodelJl25 <- glmmTMB(survived ~ treatment + stage + heightinitial +
                                   (1|species/plotnumber), 
                                 family=binomial(link="logit"), 
                                 data=July25th)
summary(fullsurvivalmodelJl25)
Anova(fullsurvivalmodelJl25, type="III")
confint(fullsurvivalmodelJl25)

#July 31st
fullsurvivalmodelJl31 <- glmmTMB(survived ~ treatment + stage + heightinitial +
                                   (1|species/plotnumber), 
                                 family=binomial(link="logit"), 
                                 data=July31st)
summary(fullsurvivalmodelJl31)
Anova(fullsurvivalmodelJl31, type="III")
confint(fullsurvivalmodelJl31)

#August 8th
fullsurvivalmodelAu8 <- glmmTMB(survived ~ treatment + stage +
                                  heightinitial +
                                  (1|species/plotnumber), 
                                family=binomial(link="logit"), 
                                data=August8th)
summary(fullsurvivalmodelAu8)
Anova(fullsurvivalmodelAu8, type="III")
confint(fullsurvivalmodelAu8)

#Analyse individual survival for each species at final time period
#Model set 3

finaltimespecies <- split(August8th, August8th$species)
Achmifinaltime <- finaltimespecies$Achmi
Agafofinaltime <- finaltimespecies$Agafo
Andgefinaltime <- finaltimespecies$Andge
Ascsyfinaltime <- finaltimespecies$Ascsy
Asctufinaltime <- finaltimespecies$Asctu
Astazfinaltime <- finaltimespecies$Astaz
Boucufinaltime <- finaltimespecies$Boucu
Bougrfinaltime <- finaltimespecies$Bougr
Liaasfinaltime <- finaltimespecies$Liaas
Monfifinaltime <- finaltimespecies$Monfi
Muhrafinaltime <- finaltimespecies$Muhra
Panvifinaltime <- finaltimespecies$Panvi
Rudhifinaltime <- finaltimespecies$Rudhi
Schscfinaltime <- finaltimespecies$Schsc
Sornufinaltime <- finaltimespecies$Sornu
Spohefinaltime <- finaltimespecies$Spohe

#Achmi 
Achmimodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Achmifinaltime)
summary(Achmimodel)
Anova(Achmimodel, type="III")
confint(Achmimodel)

#Agafo
Agafomodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Agafofinaltime)
summary(Agafomodel)
Anova(Agafomodel, type="III")
confint(Agafomodel)

#Andge
Andgemodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Andgefinaltime)
summary(Andgemodel)
Anova(Andgemodel, type="III")
confint(Andgemodel)

#Ascsy
Ascsymodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Ascsyfinaltime)
summary(Ascsymodel)
Anova(Ascsymodel, type="III")
confint(Ascsymodel)

#Asctu
Asctumodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Asctufinaltime)
summary(Asctumodel)
Anova(Asctumodel, type="III")
confint(Asctumodel)

#Astaz
Astazmodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Astazfinaltime)
summary(Astazmodel)
Anova(Astazmodel, type="III")
confint(Astazmodel)

#Boucu
Boucumodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Boucufinaltime)
summary(Boucumodel)
Anova(Boucumodel, type="III")
confint(Boucumodel)

#Bougr
Bougrmodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Bougrfinaltime)
summary(Bougrmodel)
Anova(Bougrmodel, type="III")
confint(Bougrmodel)

#Liaas
Liaasmodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Liaasfinaltime)
summary(Liaasmodel)
Anova(Liaasmodel, type="III")
confint(Liaasmodel)

#Monfi 
Monfimodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Monfifinaltime)
summary(Monfimodel)
Anova(Monfimodel, type="III")
confint(Monfimodel)

#Muhra 
Muhramodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Muhrafinaltime)
summary(Muhramodel)
Anova(Muhramodel, type="III")
confint(Muhramodel)

#Panvi 
Panvimodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Panvifinaltime)
summary(Panvimodel)
Anova(Panvimodel, type="III")
confint(Panvimodel)

#Rudhi 
Rudhimodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Rudhifinaltime)
summary(Rudhimodel)
Anova(Rudhimodel, type="III")
confint(Rudhimodel)

#Schsc 
Schscmodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Schscfinaltime)
summary(Schscmodel)
Anova(Schscmodel, type="III")
confint(Schscmodel)

#Sornu 
Sornumodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Sornufinaltime)
summary(Sornumodel)
Anova(Sornumodel, type="III")
confint(Sornumodel)

#Spohe
Spohemodel <- glmmTMB(survived ~ treatment + stage + heightinitial +
                        (1|plotnumber),
                      family=binomial(link="logit"), 
                      data=Spohefinaltime)
summary(Spohemodel)
Anova(Spohemodel, type="III")
confint(Spohemodel)

###########################################################################################

#Part 6: Main Figures
#Note: this code provides the basis for all figures; additional modifications to enhance 
#visualisation only were made in the program Inkscape.

#Figure 1

#Data frames grouped by type of damage (from descriptive statistics above)
insectdamage <- data.frame(Treatment = c("C", "I", "IF", "C", "I", "IF"),
                           Stage = c("Early", "Early", "Early", "Mid", "Mid", "Mid"),
                           mean = c(10.3, 2.65, 2.91, 7.29, 1.53, 2.13), 
                           se = c(0.7397884, 0.2438133, 0.3442206, 0.5346492, 0.1684607, 0.2529999))
fungaldamage <- data.frame(Treatment = c("C", "I", "IF", "C", "I", "IF"),
                           Stage = c("Early", "Early", "Early", "Mid", "Mid", "Mid"),
                           mean = c(2.74, 2.25, 1.00, 2.98, 2.28, 0.71), 
                           se = c(0.1354408, 0.1217978, 0.1137342, 0.25354379, 0.18256005, 0.03940073))

fig1a <- ggplot(data=insectdamage, aes(x=Stage, y=mean, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylab("Mean percent of leaf area damaged (±SE)") + xlab("Community successional stage") +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,
                position=position_dodge(.9)) +
  theme(legend.position=c(0.3, 0.7)) + ylim(0, 12)
fig1a

fig1b <- ggplot(data=fungaldamage, aes(x=Stage, y=mean, fill=Treatment)) +
  geom_bar(stat="identity", position=position_dodge()) +
  ylab("Mean percent of leaf area damaged (±SE)") + xlab("Community successional stage") +
  scale_fill_brewer(palette="Blues") + theme_classic() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se), width=.1,
                position=position_dodge(.9)) +
  theme(legend.position=c(0.3, 0.7)) + ylim(0, 12)
fig1b

fig1 <- ggarrange(fig1a, fig1b, nrow=2, labels=c("A", "B"))
fig1

#Figure 2

#Figure 2a
propdataallplotting <- propdataall
#Add extra data for 100% survival at week 0 (i.e. when planted)
wk0 <- propdataallplotting[propdataallplotting$date == '1', ]
wk0['mean'] <- 1
wk0['date'] <- 0
wk0$date <- as.factor(wk0$date)
propdataallplotting <- bind_rows(propdataallplotting, wk0)
propdataallplotting$date <- recode_factor(propdataallplotting$date, `0`="0",
                                          `1`="1", `2`="2.14", `3`="3", `4`="4.71",
                                          `5`="5.71", `6`="6.57", `7`="7.71")
propdataallplotting$date <- as.numeric(as.character(propdataallplotting$date))

fig2a <- ggplot(data=propdataallplotting, aes(x=date, y=mean, color=treatment)) +
  geom_smooth(aes(fill=treatment, linetype=stage), alpha=0.1) +
  scale_colour_manual(values=c('#000000', '#999999','#E69F00')) +
  scale_fill_manual(values=c('#000000', '#999999','#E69F00')) +
  xlab("Weeks since transplanting") + ylab("Mean proportion of seedlings survived") +
  theme_classic() +
  labs(color="Treatment", linetype="Stage") +
  scale_y_continuous(breaks=seq(0.4, 1, 0.2)) +
  scale_x_continuous(breaks=seq(0, 8, 2))
fig2a

#Figure 2b
#Effect sizes of treatment through time - plot level (outcomes of model set 1a)
#read in file containing backtransformed (exponentiated) effect sizes
#note time since sampling is recorded as continuous variable to fairly reflect
#intervals of sampling
effect <- read.csv("effect_sizes_plotlevel.csv", header=T)
fig2b <- ggplot(data=effect, aes(x=time, y=ef, group=treatment, color=treatment)) +
  geom_point(position = position_dodge(0.5), size=2.5) + 
  geom_errorbar(aes(ymin=cilower, ymax=ciupper), width=0.5,
                position=position_dodge(0.5)) +
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00')) +
  xlab("Weeks after transplanting") + ylab("Odds of survival relative to control") + labs(color="Treatment") +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(0.5, 3) + theme(legend.position = c(0.22, 0.82))
fig2b

#Figure 3

#Figure 3a (light)
visreg(Au8prop2, "lightavailable")

it_vr <- visreg(Au8prop2, "lightavailable", ylab="Relative survival", xlab="Light availability")
it_vr_df <- as.data.frame(cbind(it_vr$fit$light, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$light,
                                it_vr$res$visregRes,
                                August8thprop2$treatment))
colnames(it_vr_pr) <- c("it", "resids", "treatment")
it_vr_pr$treatment <- as.factor(it_vr_pr$treatment)
it_vr_pr$treatment <- recode_factor(it_vr_pr$treatment, `1` = "C", `2` = "I", `3` = "IF")

fig3a <- ggplot() +
  geom_point(aes(y=resids, x= it, color=treatment), alpha = 0.8, data = it_vr_pr) +
  geom_ribbon(aes(y= fit, x= it, ymin=lwr, ymax=upr),
              data=it_vr_df, linetype=0, alpha=0.3) +
  geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2, linetype="dashed") +
  scale_colour_manual(values=c('#000000', '#999999','#E69F00')) +
  theme_classic() +
  ylab("Relative seedling survival") +
  xlab("Percent light availability (scaled)") +
  labs(color="Treatment")
fig3a
ggsave("light survival.svg")

#Figure 3b (richness)
visreg(Au8prop2, "richness")

it_vr <- visreg(Au8prop2, "richness", ylab="Relative survival", xlab="Light reduction")
it_vr_df <- as.data.frame(cbind(it_vr$fit$richness, it_vr$fit$visregFit,
                                it_vr$fit$visregLwr, it_vr$fit$visregUpr))
colnames(it_vr_df) <- c("it", "fit", "lwr", "upr")
it_vr_pr <- as.data.frame(cbind(it_vr$res$richness,
                                it_vr$res$visregRes,
                                August8thprop2$treatment))
colnames(it_vr_pr) <- c("it", "resids", "treatment")
it_vr_pr$treatment <- as.factor(it_vr_pr$treatment)
it_vr_pr$treatment <- recode_factor(it_vr_pr$treatment, `1` = "C", `2` = "I", `3` = "IF")

fig3b <- ggplot() +
  geom_point(aes(y=resids, x= it, color=treatment), alpha = 0.8, 
             position=position_jitter(width=0.08), data = it_vr_pr) +
  geom_ribbon(aes(y= fit, x= it, ymin=lwr, ymax=upr),
              data=it_vr_df, linetype=0, alpha=0.3) +
  geom_line(aes(x=it, y=fit), data=it_vr_df, size = 1.2) +
  scale_colour_manual(values=c('#000000', '#999999','#E69F00')) +
  theme_classic() +
  ylab("Relative seedling survival") +
  xlab("Plot richness (scaled)")
fig3b

fig3 <- ggarrange(fig3a, fig3b, nrow=2, labels=c("A", "B"))
fig3

#Figure 4

speffect <- read.csv("effect_size_species.csv", header=T)

fig4 <- ggplot(data=speffect, aes(x=species, y=ef, group=treatment, color=treatment)) +
  geom_point(position = position_dodge(0.5), size=2.1) + 
  geom_errorbar(aes(ymin=cilower, ymax=ciupper), width=0.5,
                position=position_dodge(0.5)) +
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00')) +
  xlab("Species") + ylab("Odds of survival \nrelative to control") + labs(color="Treatment") +
  geom_hline(yintercept=1, linetype="dashed") +
  scale_x_discrete(limits=c("A. azureus", "A. foeniculum", "L. aspera",
                            "A. millefolium", "A. syriaca", "A. tuberosa", "R. hirta", 
                            "M. fistulosa", "S. scoparium", "S. nutans", 
                            "S. heterolepsis", "B. curtipendula", "B. gracilis",  
                            "M. racemosa", "A. gerardium", "P. virgatum")) +
  theme(axis.text.x = element_text(angle=45)) +
  scale_y_continuous(trans='log10')
fig4

#Figure 5

fig5 <- barplot.phylo4d(phylofile)

#Figure 6

#Figure 6a
fig6a <- ggbiplot(traits.pca, groups=traitspace2$funcgroup, labels=NULL, 
                   scale=0) +
  geom_point(aes(color=traitspace2$funcgroup), size=2) + 
  scale_colour_manual(values=c('#999999','#E69F00')) + 
  xlim(-4, 2) + theme_classic() + theme(legend.position='none') +
  xlab("PC1 Leaf Traits (49.0% explained var.)") + ylab("PC2 Leaf Traits (18.5% explained variation)")
fig6a

#Figure 6b
fig6b <- ggplot(pcaplotting, aes(x=PC1aboveground, y=logeffectsizeI, 
                                       color=funcgroup, fill=funcgroup)) +
  geom_point(size=2) + geom_smooth(method=lm) +
  scale_colour_manual(values=c('#999999','#E69F00')) + 
  scale_fill_manual(values=c('#999999','#E69F00')) +
  xlab("PC1 Leaf Traits") + ylab("Insecticide effect size (log)") +
  labs(color="Functional \ngroup", fill="Functional \ngroup") + 
  xlim(-4, 2) + theme_classic() +
  theme(legend.position=c(0.82, 0.87)) 
fig6b
ggsave("effectIplot.svg", width=95, height=95, units="mm")

fig6 <- ggarrange(fig6a, fig6b, nrow=2, labels=c("A", "B"))
fig6

###########################################################################################

#Part 7: Supplementary figures

#Figure S2

figS2a <- ggplot(pcaplotting, aes(x=PC1aboveground, y=PC2aboveground, 
                                          color=funcgroup)) +
  geom_point(size=2) + geom_text(label=pcaplotting$species, size=3.2) + 
  scale_colour_manual(values=c('#999999','#E69F00')) + 
  xlab("PC1 Leaf Traits (49.0% explained var.)") + 
  ylab("PC2 Leaf Traits (18.5% explained variation)") +
  labs(color="Functional \ngroup") + 
  theme_classic() + xlim(-4, 3) + theme(legend.position = c(0.7, 0.87))
figS2a

figS2b <- ggplot(pcaplotting, aes(x=PC1root, y=PC2root, 
                                          color=funcgroup)) +
  geom_point(size=2) + geom_text(label=pcaplotting$species, size=3.2) + 
  scale_colour_manual(values=c('#999999','#E69F00')) + 
  xlab("PC1 Root Traits (51.2% explained var.)") + 
  ylab("PC2 Root Traits (34.0% explained variation)") +
  theme_classic() + theme(legend.position = "none") + xlim(-3, 3)
figS2b

#Figure S4

#Effect sizes of treatment through time - individual level (outcomes of model set 2)
#read in file containing backtransformed (exponentiated) effect sizes
#note time since sampling is recorded as continuous variable to fairly reflect
#intervals of sampling
effectind <- read.csv("effect_sizes_individuallevel.csv", header=T)
figS4 <- ggplot(data=effectind, aes(x=time, y=ef, group=treatment, color=treatment)) +
  geom_point(position = position_dodge(0.5), size=2.5) + 
  geom_errorbar(aes(ymin=cilower, ymax=ciupper), width=0.5,
                position=position_dodge(0.5)) +
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00')) +
  xlab("Weeks after planting") + ylab("Odds of survival \nrelative to control") + labs(color="Treatment") +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(0.5, 3)
figS4

#Figure S5

figS5 <- ggplot(data=dead, aes(x=date, y=newdead, color=treatment)) +
  geom_smooth(aes(fill=treatment), alpha=0.1) +
  scale_colour_manual(values=c('#000000', '#999999','#E69F00')) +
  scale_fill_manual(values=c('#000000', '#999999','#E69F00')) +
  xlab("Weeks since planting") + ylab("Mean number of new \ndead seedlings per plot ") +
  labs(color="Treatment", fill="Treatment") + 
  theme_classic() +
  theme(legend.position = c(0.7, 0.7)) +
  scale_x_continuous(breaks=seq(1, 8, 1))
figS5

#Figure S6
#Plot the effect sizes of stage, light and richness through time

effects <- read.csv("effect_sizes_other.csv", header=T)

stage <- ggplot(data=effects, aes(x=time, y=efstage)) +
  geom_point(position = position_dodge(0.5), size=2.5, color='#E69F00') + 
  geom_errorbar(aes(ymin=cilowerstage, ymax=ciupperstage), width=0.5,
                position=position_dodge(0.5), color='#E69F00') +
  theme_classic() +
  xlab("Weeks after planting") + 
  ylab("Odds of survival in mid-succesion \nrelative to early-succesion community") +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(0, 2) +
  theme(axis.title.x = element_blank())
stage

light <- ggplot(data=effects, aes(x=time, y=eflight)) +
  geom_point(position = position_dodge(0.5), size=2.5, color='#E69F00') + 
  geom_errorbar(aes(ymin=cilowerlight, ymax=ciupperlight), width=0.5,
                position=position_dodge(0.5), color='#E69F00') +
  theme_classic() +
  xlab("Weeks since planting") + 
  ylab("Odds of survival for each percent \nincrease in light availability") +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(0, 2) +
  theme(axis.title.x = element_blank())
light

richness <- ggplot(data=effects, aes(x=time, y=efrichness)) +
  geom_point(position = position_dodge(0.5), size=2.5, color='#E69F00') + 
  geom_errorbar(aes(ymin=cilowerrichness, ymax=ciupperrichness), width=0.5,
                position=position_dodge(0.5), color='#E69F00') +
  theme_classic() +
  xlab("Weeks since planting") + 
  ylab("Odds of survival for each additional \nspecies in the plot") +
  geom_hline(yintercept=1, linetype="dashed") +
  ylim(0, 2)
richness

figS6 <- ggarrange(plot4, plot5, plot6, nrow=3, 
                             labels=c("A", "B", "C"))
figS6

#Figure S7

figS7 <- ggplot(data=August8thprop2, aes(x=species, y=mean, 
                                                 color=stage, group=stage)) +
  geom_jitter(position=position_jitterdodge(jitter.width=0.5, jitter.height = 0.02, 
                                            dodge.width=0.6), alpha=0.5) +
  stat_summary(fun.data="mean_sdl",  fun.args = list(mult=1), 
               size=0.5, shape=17, position=position_dodge(0.8)) +
  theme_classic() +
  scale_x_discrete(limits=c("Achmi", "Agafo", "Ascsy", "Asctu", 
                            "Astaz", "Liaas", "Monfi", "Rudhi", 
                            "Andge", "Boucu", "Bougr","Muhra",
                            "Panvi" ,"Spohe", "Schsc", "Sornu")) +
  theme(axis.text.x = element_text(angle=45)) +
  xlab("Species") + ylab("Proportion of seedlings \nsurviving per plot") + labs(fill="Stage") +
  scale_color_manual(values=c('#999999','#E69F00'))
figS7

#Figure S8

figS8a <- ggplot(insecticide, aes(x=meaninsectdamageC, y=ef, color=growth)) +
  geom_point(size=2.5) + scale_colour_manual(values=c('#999999','#E69F00')) + 
  xlab("Mean insect damage on control plants") + 
  ylab("Insecticide \ntreatment effect size") +
  labs(color="Functional \ngroup") + 
  theme_classic() 
figS8a

figS8b <- ggplot(fungicide, aes(x=meaninsectfungaldamageC, y=ef, color=growth)) +
  geom_point(size=2.5) + scale_colour_manual(values=c('#999999','#E69F00')) + 
  xlab("Mean insect plus fungal damage on control plants") + 
  ylab("Insecticide plus fungicide \ntreatment effect size") +
  labs(color="Functional \ngroup") + 
  theme_classic() 
figS8b

figS8 <- ggarrange(figS8a, figS8b, nrow=2, labels=c("A", "B"))
figS8

###########################################################################################

#A joke to finish

#Q: Where do bad desserts get taken?

#A: Into custardy
