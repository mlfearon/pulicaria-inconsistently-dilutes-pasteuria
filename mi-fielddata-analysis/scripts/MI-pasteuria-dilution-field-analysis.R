# Analysis of field Pasteuria data and dilution effect for "Inconsistent dilution: Experimental but not field evidence for a dilution effect in Daphnia–bacteria interactions"


# Submitted to: Oecologia


# Code written by Michelle Fearon
# Last updated: June 29, 2023

## This code produces four models to test how host density, species richness, and year correlate with
## either maximal Pasteuria prevalence (models A & B) or area under the prevalence curve (models C & D) in D. dentifera, 
## with two models for each response variable. Models A & C included the host densities and species richness 
## from the date of max prevalence, while models B & D included the mean host densities and total species richness 
## for each lake and year combination. We then used model selection to determine the best version of each model A-D,
## and complemented that approach with model averaging.


#libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(DescTools)
library(lme4)
library(lmerTest)
library(glmmTMB)
library(car)
library(DHARMa)
library(performance)
library(MuMIn)
library(ggeffects)
library(ggpubr)
library(vegan)
library(pscl)
library(grid)
library(gridExtra)
library(reshape)
library(ggiraph)
library(ggiraphExtra)
library(moonBook)
library(RColorBrewer)
library(here)



# set the path to the script relative to the project root directory
here::i_am("mi-fielddata-analysis/scripts/MI-pasteuria-dilution-field-analysis.R")



# function to make a grid plot with a shared legend
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position="none"))
  gl <- c(gl, ncol = ncol, nrow = nrow)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  
  grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}



## updated overdispersion function from Ben Bolker
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}


#==========================================================================================================

# load cleaned parasite and host density data 2014 - 2017
data_all <- read.csv(here("mi-fielddata-analysis/data/Clean-Data-2014-2021_All-Host-Densities.csv"), stringsAsFactors = F)
# Total = total number of daphnia hosts for that species/site/date
# .inf = the total number of infections for each parasite as a column
# .prev = the calculated prevalence of each parasite (only when the total host density for that species/site/date > 20)
# host.density = total area density of the host species for that site/date
# .den = density of each host species for that site/date
# total.density = summed density of all host species for that site/date


# Determine our minimum detection threshold for host density (consistent across all host species)
min_detection <- filter(data_all, Host.Density > 0, Host.Density < 100, Year < 2018)
min(min_detection$Host.Density)  # 53.8 animals per m squared is the minimum density that we can detect


# get data to only include dentifera for main analyses
data <- data_all %>% 
  filter(Year != 2018 & Year != 2020 & Year != 2021) %>% 
  filter(Host.Species == "dentifera") %>% # get infection prevalence data for dentifera only
  filter(!is.na(Host.Density)) %>%  # remove sample dates with missing density data
  mutate(past.density.inf = pasteuria.prev * Host.Density, # add pasteuria infected density to data set
         Diluter.density = Pulicaria.density + Retrocurva.density, # diluter host densities
         OtherHost.density = Pulicaria.density + Retrocurva.density + Ceriodaphnia.density + Dubia.density +
           Parvula.density + Ambigua.density + Mendotae.density) %>%
  select(Host.Species, Lake, Year, Julian.Day, Julian, Total, starts_with("past"), Host.Density:OtherHost.density) # Select only the variables needed in cleaned data set
dim(data)

# diversity metrics per sampling date
data$species.richness <- rowSums(select(data, Dentifera.density:Mendotae.density) > 0)
data$shannon <- diversity(select(data, Dentifera.density:Mendotae.density), "shannon")

# set all zero densities to half of the minimum detection threshold (lowest observed areal density)
# to reduce zeros in the data for later log transformation and correct for rare hosts that may have been present but below our detection threshold
min_density <- 53.8/2

# initial check of zero densities
data %>% select(Host.Density:OtherHost.density)

# replace zeros with 'min_density" value for each host density column [Note: an error is produced but it still works when you check below]
data$Host.Density <- replace(data$Host.Density, data$Host.Density == 0, min_density)
data$Total.Density <- replace(data$Total.Density, data$Total.Density == 0, min_density)
data$Dentifera.density <- replace(data$Dentifera.density, data$Dentifera.density == 0, min_density)
data$Pulicaria.density <- replace(data$Pulicaria.density, data$Pulicaria.density == 0, min_density)
data$Retrocurva.density <- replace(data$Retrocurva.density, data$Retrocurva.density == 0, min_density)
data$Ceriodaphnia.density <- replace(data$Ceriodaphnia.density, data$Ceriodaphnia.density == 0, min_density)
data$Dubia.density <- replace(data$Dubia.density, data$Dubia.density == 0, min_density)
data$Parvula.density <- replace(data$Parvula.density, data$Parvula.density == 0, min_density)
data$Ambigua.density <- replace(data$Ambigua.density, data$Ambigua.density == 0, min_density)
data$Mendotae.density <- replace(data$Mendotae.density, data$Mendotae.density == 0, min_density)
data$Diluter.density <- replace(data$Diluter.density, data$Diluter.density == 0, min_density)
data$OtherHost.density <- replace(data$OtherHost.density, data$OtherHost.density == 0, min_density)

# check that it worked
data %>% select(Host.Density:OtherHost.density)
sum(data %>% select(Host.Density:OtherHost.density)==0) # should be zero if all have been changed

# check data set
str(data)
head(data)
dim(data)
#View(data)

range(data$Julian)


######################################################################
### Calculate the maximal and mean infection prevalence per host/lake/year
######################################################################


# calculate the max, mean, sd, and se of pasteuria prevalence for each host, lake and year
sum.data <- data %>%
  group_by(Host.Species, Lake, Year) %>%  # group by host species, lake and year
  filter(Julian < 305) %>% # remove November dates, only include July - October for max prev and mean densities
  summarize(past.max.prev = max(pasteuria.prev, na.rm = T),  # calculate max pasteuria infection prevalence
            past.mean.prev = mean(pasteuria.prev, na.rm = T),# calculate mean pasteuria infection prevalence
            past.sd.prev = sd(pasteuria.prev, na.rm = T),    # calculate sd pasteuria infection prevalence
            past.se.prev = sd(pasteuria.prev, na.rm = T)/sqrt(length(pasteuria.prev)),  # calculate standard error pasteuria infection prevalence
            past.max.inf = max(past.density.inf, na.rm = T),  # maximum pf pasteuria infection density
            mean.total.count = mean(Total, na.rm = T),   # mean of the total dentifera hosts counted
            mean.density = mean(Host.Density, na.rm = T),    # mean host density of each species
            mean.tot.density = mean(Total.Density, na.rm = T),# mean total density of all hosts
            mean.dentifera.density = mean(Dentifera.density, na.rm = T),# mean dentifera density
            mean.pulicaria.density = mean(Pulicaria.density, na.rm = T),# mean pulicaria density
            mean.retrocurva.density = mean(Retrocurva.density, na.rm = T), # mean retrocurva density
            mean.ceriodaphnia.density = mean(Ceriodaphnia.density, na.rm = T), # mean ceriodaphnia density
            mean.dubia.density = mean(Dubia.density, na.rm = T), # mean Dubia density
            mean.parvula.density = mean(Parvula.density, na.rm = T), # mean parvula density
            mean.ambigua.density = mean(Ambigua.density, na.rm = T), # mean Ambigua density
            mean.mendotae.density = mean(Mendotae.density, na.rm = T), # mean Mendotae density
            mean.diluter.density = mean(Diluter.density, na.rm = T),  # mean of diluter hosts (pulicaria + retrocurva)
            mean.other.density = mean(OtherHost.density, na.rm = T), # mean of all other hotsts (non-dentifera)
            mean.richness = mean(species.richness, na.rm = T),       # mean species richness
            mean.shannon = mean(shannon, na.rm = T)) %>%             # mean shannon diversity index
  ungroup()
 # note that the error produced is for one lake year (2015 Whitemore) that does not have any records of pasteuria in dentifera that were greater than 20 individuals counted, therefore prevalence could not be calculated for those dates.
#View(sum.data)
dim(sum.data)


### Calculate mean August densities
Aug.densities <- data %>% 
  filter(Julian > 212, Julian < 244) %>%  # restrict data to just August dates for each year
  group_by(Host.Species, Lake, Year) %>%
  summarize(mean.Aug.dentifera = mean(Dentifera.density, na.rm = T), # mean dentifera density
            mean.Aug.pulicaria = mean(Pulicaria.density, na.rm = T),# mean pulicaria density
            mean.Aug.retrocurva = mean(Retrocurva.density, na.rm = T), # mean retrocurva density
            mean.Aug.richness = mean(species.richness, na.rm = T))       # mean species richness

dim(Aug.densities)

# add mean August densities to sum.data
sum.data  <- left_join(sum.data, Aug.densities)

sum.data$Year <- as.character(sum.data$Year)
sum.data$mean.total.count[ is.nan(sum.data$mean.total.count)] <- 0



# add species richness and simpson diversity index to data set
sum.data$species.richness <- rowSums(select(sum.data, mean.dentifera.density:mean.mendotae.density) > min_density)


unique(sum.data$Lake)

#########################
# Area Under the Curve
##########################

# load data with area under the curve calculated, and host densities on the date of the max prevalence
auc_data <- read.csv(here("mi-fielddata-analysis/data/auc_14to21_prev.csv"), stringsAsFactors = F)
auc_data$Year <- as.character(auc_data$Year)


# select only Pasteuria data and prep to join with the sum.data 
auc_data_past <- auc_data %>%
  filter(Parasite.Species == "pasteuria.prev" & Host.Species == "dentifera"  & Year != 2020 & Year != 2021) %>%
  arrange(Host.Species, Lake) %>%
  dplyr::rename(pasteuria.auc = AUC.prev2) %>%  # Use AUC.prev2 instead of AUC.prev because it is more accurately calculated
  select(Host.Species, Lake, Year, pasteuria.auc, Max.Prevalence, Count.At.Max:shannon.at.max)
#View(auc_data_past)

str(auc_data_past)
str(sum.data)

sum.data  <- left_join(sum.data, auc_data_past)

# replace zero densities at max with the "min_density"
# replace zeros with 'min_density" value for each host density column [Note: an error is produced but it still works when you check below]
sum.data$dentifera.at.max <- replace(sum.data$dentifera.at.max, sum.data$dentifera.at.max == 0, min_density)
sum.data$pulicaria.at.max <- replace(sum.data$pulicaria.at.max, sum.data$pulicaria.at.max == 0, min_density)
sum.data$retrocurva.at.max <- replace(sum.data$retrocurva.at.max, sum.data$retrocurva.at.max == 0, min_density)
sum.data$diluter.at.max <- replace(sum.data$diluter.at.max, sum.data$diluter.at.max == 0, min_density)
sum.data$otherhost.at.max <- replace(sum.data$otherhost.at.max, sum.data$otherhost.at.max == 0, min_density)

str(sum.data)


# make year a factor to work in the figures below
sum.data$Year <- as.character(sum.data$Year)
sum.data$Host.Species <- as.factor(sum.data$Host.Species)

View(sum.data)

# filter out the lake year that doesn't have a valid max prevalence because there was not enough dentifera detected to calculate prevalence
sum.data <- filter(sum.data, !is.na(Max.Prevalence)) # removed Whitmore in 2015 row.


# Create a data set where lake-years without an epidemic are removed (Reviewer requested: compare analyses of this data set to the original full data set)
sum.data_epidemics <- filter(sum.data, Max.Prevalence > 0)  # removes 10 lake-years that did not have Pasteuria epidemics in D. dentifera

#convert to a data frame
sum.data_dentifera <- as.data.frame(sum.data)
sum.data_epi_dentifera <- as.data.frame(sum.data_epidemics)
dim(sum.data_dentifera)
dim(sum.data_epi_dentifera)

range(sum.data_dentifera$richness.at.max)
range(sum.data_dentifera$species.richness)

range(sum.data_epi_dentifera$richness.at.max)
range(sum.data_epi_dentifera$species.richness)

# Center and scale all continuous variables, densities get log + 1 transformation
sum.data_dentifera$mean.tot.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.tot.density)))
sum.data_dentifera$mean.dentifera.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.dentifera.density)))
sum.data_dentifera$mean.pulicaria.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.pulicaria.density)))
sum.data_dentifera$mean.retrocurva.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.retrocurva.density)))
sum.data_dentifera$mean.other.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.other.density)))
sum.data_dentifera$species.richness_z <- as.numeric(scale(sum.data_dentifera$species.richness))
sum.data_dentifera$mean.richness_z <- as.numeric(scale(sum.data_dentifera$mean.richness))
sum.data_dentifera$mean.shannon_z <- as.numeric(scale(sum.data_dentifera$mean.shannon))
sum.data_dentifera$dentifera.at.max_z <- as.numeric(scale(log(sum.data_dentifera$dentifera.at.max)))
sum.data_dentifera$pulicaria.at.max_z <- as.numeric(scale(log(sum.data_dentifera$pulicaria.at.max)))
sum.data_dentifera$retrocurva.at.max_z <- as.numeric(scale(log(sum.data_dentifera$retrocurva.at.max)))
sum.data_dentifera$diluter.at.max_z <- as.numeric(scale(log(sum.data_dentifera$diluter.at.max)))
sum.data_dentifera$otherhost.at.max_z <- as.numeric(scale(log(sum.data_dentifera$otherhost.at.max)))
sum.data_dentifera$richness.at.max_z <- as.numeric(scale(sum.data_dentifera$richness.at.max))
sum.data_dentifera$shannon.at.max_z <- as.numeric(scale(sum.data_dentifera$shannon.at.max))
sum.data_dentifera$pasteuria.auc_log <- as.numeric(log(sum.data_dentifera$pasteuria.auc+1))
sum.data_dentifera$mean.Aug.dentifera.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.Aug.dentifera)))
sum.data_dentifera$mean.Aug.pulicaria.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.Aug.pulicaria)))
sum.data_dentifera$mean.Aug.retrocurva.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.Aug.retrocurva)))
sum.data_dentifera$mean.Aug.richness_z <- as.numeric(scale(log(sum.data_dentifera$mean.Aug.richness)))

sum.data_epi_dentifera$mean.tot.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.tot.density)))
sum.data_epi_dentifera$mean.dentifera.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.dentifera.density)))
sum.data_epi_dentifera$mean.pulicaria.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.pulicaria.density)))
sum.data_epi_dentifera$mean.retrocurva.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.retrocurva.density)))
sum.data_epi_dentifera$mean.other.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.other.density)))
sum.data_epi_dentifera$species.richness_z <- as.numeric(scale(sum.data_epi_dentifera$species.richness))
sum.data_epi_dentifera$mean.richness_z <- as.numeric(scale(sum.data_epi_dentifera$mean.richness))
sum.data_epi_dentifera$mean.shannon_z <- as.numeric(scale(sum.data_epi_dentifera$mean.shannon))
sum.data_epi_dentifera$dentifera.at.max_z <- as.numeric(scale(log(sum.data_epi_dentifera$dentifera.at.max)))
sum.data_epi_dentifera$pulicaria.at.max_z <- as.numeric(scale(log(sum.data_epi_dentifera$pulicaria.at.max)))
sum.data_epi_dentifera$retrocurva.at.max_z <- as.numeric(scale(log(sum.data_epi_dentifera$retrocurva.at.max)))
sum.data_epi_dentifera$diluter.at.max_z <- as.numeric(scale(log(sum.data_epi_dentifera$diluter.at.max)))
sum.data_epi_dentifera$otherhost.at.max_z <- as.numeric(scale(log(sum.data_epi_dentifera$otherhost.at.max)))
sum.data_epi_dentifera$richness.at.max_z <- as.numeric(scale(sum.data_epi_dentifera$richness.at.max))
sum.data_epi_dentifera$shannon.at.max_z <- as.numeric(scale(sum.data_epi_dentifera$shannon.at.max))
sum.data_epi_dentifera$pasteuria.auc_log <- as.numeric(log(sum.data_epi_dentifera$pasteuria.auc+1))
sum.data_epi_dentifera$mean.Aug.dentifera.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.Aug.dentifera)))
sum.data_epi_dentifera$mean.Aug.pulicaria.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.Aug.pulicaria)))
sum.data_epi_dentifera$mean.Aug.retrocurva.density_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.Aug.retrocurva)))
sum.data_epi_dentifera$mean.Aug.richness_z <- as.numeric(scale(log(sum.data_epi_dentifera$mean.Aug.richness)))

#include observation level random effect
sum.data_dentifera$OLRE <- seq_len(nrow(sum.data_dentifera))

sum.data_epi_dentifera$OLRE <- seq_len(nrow(sum.data_epi_dentifera))

# Test the correlation between Pasteuria Maximum prevalence and Area Under the prevalence curve
cor.test(sum.data_dentifera$past.max.prev, sum.data_dentifera$pasteuria.auc_log)

cor.test(sum.data_epi_dentifera$past.max.prev, sum.data_epi_dentifera$pasteuria.auc_log)


# species richness is correlated with total density (r = 0.35, p = 0.006), but not individual species densities
cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.tot.density)
cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.dentifera.density)
cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.retrocurva.density)
cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.pulicaria.density)

cor.test(sum.data_epi_dentifera$species.richness, sum.data_epi_dentifera$mean.tot.density)
cor.test(sum.data_epi_dentifera$species.richness, sum.data_epi_dentifera$mean.dentifera.density)
cor.test(sum.data_epi_dentifera$species.richness, sum.data_epi_dentifera$mean.retrocurva.density)
cor.test(sum.data_epi_dentifera$species.richness, sum.data_epi_dentifera$mean.pulicaria.density)


####################################################
## MODEL A: Max Pasteuria prevalence in D. dentifera vs host densities and richness at the max prevalence date

###  Initial Model of MAX Pasteuria prev in dentifera, with densities at max prev, used beta-binomial to control overdispersion 
modA_initial <- glmmTMB(past.max.prev ~ pulicaria.at.max_z + dentifera.at.max_z + retrocurva.at.max_z + richness.at.max_z + (1|Year) + (1|Lake), 
               family = betabinomial(), weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modA_initial)
check_collinearity(modA_initial)
testDispersion(modA_initial)
testZeroInflation(modA_initial)
modA_simResid <- simulateResiduals(fittedModel = modA_initial)
plot(modA_simResid)  # model looks good

# Reviewer requested: model with data set where max prev is > zero (no infections in dentifera) [MANUSCRIPT APPENDIX ONLY]
modA_initial2 <- glmmTMB(past.max.prev ~ pulicaria.at.max_z + dentifera.at.max_z + retrocurva.at.max_z + richness.at.max_z + (1|Year) + (1|Lake), 
               family = betabinomial(), weights = Count.At.Max, data = sum.data_epi_dentifera)
summary(modA_initial2)
check_collinearity(modA_initial2)
testDispersion(modA_initial2)
testZeroInflation(modA_initial2)  # zero inflated
modA2_simResid <- simulateResiduals(fittedModel = modA_initial2)
plot(modA2_simResid)  # model looks good


#### Reviewer suggested nested model approach [CODE ONLY, NOT USED IN MANUSCRIPT]
# M1: Richness only
modA1 <- glmmTMB(past.max.prev ~ richness.at.max_z + (1|Year) + (1|Lake), 
               family = betabinomial(), weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modA1)  # richness is not significant

# M2: Richness + dentifera density
modA2 <- glmmTMB(past.max.prev ~ richness.at.max_z + dentifera.at.max_z + (1|Year) + (1|Lake), 
               family = betabinomial(), weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modA2) # none of the factors are significant

# M3: Richness + dentifera density + pulicaria density
modA3 <- glmmTMB(past.max.prev ~ richness.at.max_z + dentifera.at.max_z + pulicaria.at.max_z + (1|Year) + (1|Lake), 
               family = betabinomial(), weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modA3) # none of the factors are significant

# M4: Richness + dentifera density + pulicaria density + retrocurva density
modA4 <- glmmTMB(past.max.prev ~ richness.at.max_z + dentifera.at.max_z + pulicaria.at.max_z + retrocurva.at.max_z + (1|Year) + (1|Lake), 
               family = betabinomial(), weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modA4) # none of the factors are significant

# M5: Richness + dentifera density + all other host density
modA5 <- glmmTMB(past.max.prev ~ richness.at.max_z + dentifera.at.max_z + otherhost.at.max_z + (1|Year) + (1|Lake), 
               family = betabinomial(), weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modA5) # none of the factors are significant

# compare AIC among the 3 models
AIC(modA1, modA2, modA3, modA4, modA5)  # modA1 has the lowest AIC




##### Model Selection and Averaging [APPROACH USED IN THE MANUSCRIPT]
# model section ranking by AICc using ML
options(na.action = "na.fail")
msc <- dredge(modA_initial, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc)
(attr(msc, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList <- get.models(msc, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:

#### Table 1
summary(model.avg(fmList))



## Appendix S1 Table S5
# Export AIC table for Model A (delta AIC < 4)
msc_A <- filter(msc, delta < 4)
write.csv(msc_A, here("mi-fielddata-analysis/results/ModelA_AIC_table.csv"), quote = F, row.names = F)


## Appendix S1 Table S9
# Top Model A
modA_top <- glmmTMB(past.max.prev ~ pulicaria.at.max_z + (1|Year) + (1|Lake), 
              family = betabinomial(), weights = Count.At.Max, data = sum.data_dentifera)
summary(modA_top)
testDispersion(modA_intial)
testZeroInflation(modA_intial)
modA_simResid <- simulateResiduals(fittedModel = modA_intial)
plot(modA_simResid)  # model looks good



# Marginal effects from top model
me1_a <- ggpredict(modA_top, c("pulicaria.at.max_z [all]"))
plot(me1_a, add.data = T)


# recalculate original Pulicaria density at Max 
me1_a$x 
pulic_mean <- mean(log(sum.data_dentifera$pulicaria.at.max+1)) # mean of original richness from disease data set
pulic_sd <- sd(log(sum.data_dentifera$pulicaria.at.max+1)) # sd of original richness from disease data set
me1_a$PulicMax_log <- t((t(me1_a$x) * pulic_sd) + pulic_mean)
me1_a$PulicMax <- exp(me1_a$PulicMax_log)

display.brewer.pal(3, "Dark2")
brewer.pal(3, "Dark2")
range(me1_a$PulicMax)

# Figure 2A: Model A
# plot of Max Pasteuria prevalence vs pulicaria density at max prevalence
MaxPrev_Pulic <- ggplot(me1_a, aes(x = PulicMax, y = predicted)) +
  geom_line(linewidth = 1,  color = "#D95F02") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#D95F02") +
  geom_jitter(data = sum.data_dentifera, aes(x = pulicaria.at.max, y = past.max.prev), size = 2, alpha = 0.4, width = 0.05, color = "#D95F02", shape = 17) +
  labs(x = bquote(italic("D. pulicaria") ~ "Density at Max Prevalence (no." ~ "m"^-2*")"), y = bquote(atop("Maximal " ~ italic("Pasteuria") ~ "Prevalence", "in" ~italic("D. dentifera")))) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(x = 100, y= 0.1), label = "N.S.") +
  scale_x_log10(labels = scales::comma, lim = c(20, 150000)) +
  theme_classic()+
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(MaxPrev_Pulic)
ggsave(here("mi-fielddata-analysis/figures/MaxPastPrev_PulicariaAtMax_predict_color.tiff"), plot = MaxPrev_Pulic, dpi = 300, width = 10, height = 10, units = "cm", compression="lzw")


# second top model is the NULL model


# third top model
modA_third <- glmmTMB(past.max.prev ~ dentifera.at.max_z + (1|Year) + (1|Lake), 
                    family = betabinomial(), weights = Count.At.Max, data = sum.data_dentifera)
summary(modA_3)





############################################
## MODEL B: Max Pasteuria prevalence in D. dentifera vs mean host densities and richness from throughout the season

### Initial Model of MAX Pasteuria prev in dentifera, with mean densities, used beta-binomial to control overdispersion 
modB_initial <- glmmTMB(past.max.prev ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + (1|Year) + (1|Lake), 
              family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modB_initial)
check_collinearity(modB_initial)
testDispersion(modB_initial)
testZeroInflation(modB_initial)
modB_simResid <- simulateResiduals(fittedModel = modB_initial)
plot(modB_simResid)

# Reviewer requested: model with data set where max prev is > zero (no infections in dentifera) [MANUSCRIPT APPENDIX ONLY]
modB_initial2 <- glmmTMB(past.max.prev ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + (1|Year) + (1|Lake), 
                        family = betabinomial, weights = round(Count.At.Max), data = sum.data_epi_dentifera)
summary(modB_initial2)
check_collinearity(modB_initial2)
testDispersion(modB_initial2)
testZeroInflation(modB_initial2)
modB2_simResid <- simulateResiduals(fittedModel = modB_initial2)
plot(modB2_simResid)


# Reviewer suggested nested models
# M1: Richness only
modB1 <- glmmTMB(past.max.prev ~ species.richness_z + (1|Year) + (1|Lake), 
                        family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modB1)  # richness not significant

# M2: Richness + dentifera density
modB2 <- glmmTMB(past.max.prev ~ species.richness_z + mean.dentifera.density_z + (1|Year) + (1|Lake), 
               family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modB2) # none of the factors are significant

# M3: Richness + dentifera density + pulicaria density
modB3 <- glmmTMB(past.max.prev ~ species.richness_z + mean.dentifera.density_z + mean.pulicaria.density_z + (1|Year) + (1|Lake), 
               family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modB3) # none of the factors are significant

# M4: Richness + dentifera density + pulicaria density + retrocurva density
modB4 <- glmmTMB(past.max.prev ~ species.richness_z + mean.dentifera.density_z + mean.pulicaria.density_z + mean.retrocurva.density_z + (1|Year) + (1|Lake), 
               family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modB4) # none of the factors are significant

# M5: Richness + dentifera density + all other host density
modB5 <- glmmTMB(past.max.prev ~ species.richness_z + mean.dentifera.density_z + mean.other.density_z + (1|Year) + (1|Lake), 
               family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modB5) # none of the factors are significant

# compare AIC among the 3 models
AIC(modB1, modB2, modB3, modB4, modB5)  # modB2 has the lowest AIC (Richness and Denitfera density)



##### Model Selection and Averaging [APPROACH USED IN THE MANUSCRIPT]
# model section ranking by AICc using ML
options(na.action = "na.fail")
msc2 <- dredge(modB_initial, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc2)
(attr(msc2, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList2 <- get.models(msc2, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:

#### Table 1
summary(model.avg(fmList2))



## Appendix S1 Table S6
# Export AIC table for Model B (delta AIC < 4)
msc_B <- filter(msc2, delta < 4)
write.csv(msc_B, here("mi-fielddata-analysis/results/ModelB_AIC_table.csv"), quote = F, row.names = F)


## Appendix S1 Table S9
# Top Model B
modB_top <- glmmTMB(past.max.prev ~ mean.dentifera.density_z + (1|Year) + (1|Lake), 
             family = betabinomial(), weights = Count.At.Max, data = sum.data_dentifera)
summary(modB_top)


# second top model
modB_second <- glmmTMB(past.max.prev ~ mean.dentifera.density_z + mean.retrocurva.density_z + (1|Year) + (1|Lake), 
                    family = betabinomial(), weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modB_second)



# Marginal effects from top model
me2_b <- ggpredict(modB_top, c("mean.dentifera.density_z [all]"))
plot(me2_b, add.data = T)



# recalculate original mean Dentifera density
me2_b$x 
dent_mean <- mean(log(sum.data_dentifera$mean.dentifera.density+1)) # mean of original richness from disease data set
dent_sd <- sd(log(sum.data_dentifera$mean.dentifera.density+1)) # sd of original richness from disease data set
me2_b$DentMean_log <- t((t(me2_b$x) * dent_sd) + dent_mean)
me2_b$DentMean <- exp(me2_b$DentMean_log)


## Figure 2B: Model B
# plot of Max Pasteuria prevalence vs mean dentifera density
MaxPrev_Dent <- ggplot(me2_b, aes(x = DentMean, y = predicted)) +
  geom_line(linewidth = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera, aes(x = mean.dentifera.density, y = past.max.prev), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote("Mean" ~italic("D. dentifera") ~ "Density (no." ~ "m"^-2*")"), y = bquote(atop("Maximal " ~ italic("Pasteuria") ~ "Prevalence", "in" ~italic("D. dentifera")))) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(x = 1000, y= 0.1), label = "p = 0.017") +
  scale_x_log10(labels = scales::comma) +
  theme_classic()+
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(MaxPrev_Dent)
ggsave(here("mi-fielddata-analysis/figures/MaxPastPrev_MeanDentifera_predict_color.tiff"), plot = MaxPrev_Dent, dpi = 300, width = 10, height = 10, units = "cm", compression="lzw")




############################################
## MODEL Testing if early host densities impact peak maximum Pasteuria prevalence 
# Max Pasteuria prevalence in D. dentifera vs mean August host densities and richness

### Initial Model of MAX Pasteuria prev in dentifera, with mean AUGUST densities, used beta-binomial to control overdispersion 
modAUG_initial <- glmmTMB(past.max.prev ~ mean.Aug.pulicaria.density_z + mean.Aug.dentifera.density_z + mean.Aug.retrocurva.density_z + mean.Aug.richness_z + (1|Year) + (1|Lake), 
                        family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modAUG_initial)
check_collinearity(modAUG_initial)
testDispersion(modAUG_initial)
testZeroInflation(modAUG_initial)
modAUG_simResid <- simulateResiduals(fittedModel = modAUG_initial)
plot(modAUG_simResid)


##### Model Selection and Averaging [APPROACH USED IN THE MANUSCRIPT]
# model section ranking by AICc using ML
options(na.action = "na.fail")
mscAUG <- dredge(modAUG_initial, rank = "AICc", trace = TRUE, REML = FALSE)
print(mscAUG)
(attr(mscAUG, "rank.call"))
# Get the models (fitted by REML, as in the global model)
mscAUG_update <- filter(mscAUG, delta < 2)
fmListAUG <- get.models(mscAUG_update, delta < 2) # subsetting isn't working, did it manually below
fmListAUG[-c(3:16)] # (need to remove the last model that doesn't converge)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:

#### MODEL AVERAGE 
summary(model.avg(fmListAUG[-c(3:16)], delta < 2))


## Appendix S1 Table SXX
# Export AIC table for Model XX (delta AIC < 4)
msc_AUG <- filter(mscAUG, delta < 4)
write.csv(msc_AUG, here("mi-fielddata-analysis/results/Model_MaxPrev_AugDensities_AIC_table.csv"), quote = F, row.names = F)



# Top model
modAUG_top <- glmmTMB(past.max.prev ~ mean.Aug.dentifera.density_z + (1|Year) + (1|Lake), 
                          family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modAUG_top)

# second top model
modAUG_second <- glmmTMB(past.max.prev ~ mean.Aug.dentifera.density_z + mean.Aug.pulicaria.density_z + (1|Year) + (1|Lake), 
                      family = betabinomial, weights = round(Count.At.Max), data = sum.data_dentifera)
summary(modAUG_second)




##################################
## MODEL C: AUC (Integrated) Pasteuria prevalence in D. dentifera vs host densities and richness at the max prevalence

# remove one instance of NA where there was only a single date where enough dentifera were counted during the season, so a max prev could be caluclated but not AUC (Whitmore 2016 removed)
sum.data_dentifera2 <- filter(sum.data_dentifera, !is.na(pasteuria.auc))  

### Initial Model of AUC Pasteuria prev in dentifera, with densities at the max prevalence
modC_initial <- lmer(pasteuria.auc_log ~ pulicaria.at.max_z + dentifera.at.max_z + retrocurva.at.max_z + richness.at.max_z + (1|Year) + (1|Lake), 
             data = sum.data_dentifera2) 
summary(modC_initial)
AIC(modC_initial)
vif(modC_initial)
plot(modC_initial)
qqnorm(resid(modC_initial))
qqline(resid(modC_initial)) # somewhat heavy lower tail that falls off of the QQ normal line (looks better when using data with max prev > 0)

# Reviewer requested: model with data set where max prev is > zero (no infections in dentifera) [MANUSCRIPT APPENDIX ONLY]
modC_initial2 <- lmer(pasteuria.auc_log ~ pulicaria.at.max_z + dentifera.at.max_z + retrocurva.at.max_z + richness.at.max_z + (1|Year) + (1|Lake), 
             data = sum.data_epi_dentifera) 
summary(modC_initial2)
AIC(modC_initial2)
vif(modC_initial2)
plot(modC_initial2)
qqnorm(resid(modC_initial2))
qqline(resid(modC_initial2)) # looks better when using data with max prev > 0



# Reviewer suggested nested models
# M1: Richness only
modC1 <- lmer(pasteuria.auc_log ~ richness.at.max_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modC1)  # richness is not sig
AIC(modC1)

# M2: Richness + dentifera density
modC2 <- lmer(pasteuria.auc_log ~ richness.at.max_z + dentifera.at.max_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modC2) # none of the factors are significant

# M3: Richness + dentifera density + pulicaria density
modC3 <- lmer(pasteuria.auc_log ~ richness.at.max_z + dentifera.at.max_z + pulicaria.at.max_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modC3) # none of the factors are significant

# M4: Richness + dentifera density + pulicaria density + retrocurva density
modC4 <- lmer(pasteuria.auc_log ~ richness.at.max_z + dentifera.at.max_z + pulicaria.at.max_z + retrocurva.at.max_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modC4) # none of the factors are significant

# M5: Richness + dentifera density + all other host density
modC5 <- lmer(pasteuria.auc_log ~ richness.at.max_z + dentifera.at.max_z + otherhost.at.max_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modC5) # none of the factors are significant

# compare AIC among the 3 models
AIC(modC1, modC2, modC3, modC4, modC5)  # modC1 has the lowest AIC




# model section ranking by AICc using ML
options(na.action = "na.fail")
msc3 <- dredge(modC_initial, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc3)
(attr(msc3, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList3 <- get.models(msc3, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:

#### Table 1
summary(model.avg(fmList3))


## Appendix S1 Table S7
# Export AIC table for Model C (delta AIC < 4)
msc_C <- filter(msc3, delta < 4)
write.csv(msc_C, here("mi-fielddata-analysis/results/ModelC_AIC_table.csv"), quote = F, row.names = F)


## Appendix S1 Table S9
# Top Model C
modC_top <- lmer(pasteuria.auc_log ~ dentifera.at.max_z + (1|Year)+ (1|Lake), data = sum.data_dentifera2)
summary(modC_top)


# second top model is the Null model


# third top model
modC_third <- lmer(pasteuria.auc_log ~ dentifera.at.max_z + retrocurva.at.max_z + (1|Lake), data = sum.data_dentifera2)
summary(modC_third)




# Marginal effects from top model
me3_c <- ggpredict(modC_top, c("dentifera.at.max_z"))
plot(me3_c, add.data = T)



# recalculate original Dentifera density at max prevalence
me3_c$x 
dent_mean <- mean(log(sum.data_dentifera$dentifera.at.max+1)) # mean of original richness from disease data set
dent_sd <- sd(log(sum.data_dentifera$dentifera.at.max+1)) # sd of original richness from disease data set
me3_c$DentMax_log <- t((t(me3_c$x) * dent_sd) + dent_mean)
me3_c$DentMax <- exp(me3_c$DentMax_log)
me3_c$predicted_backtransformed <- exp(me3_c$predicted)-1  #back transform the predicted and conf interval values from log+1 
me3_c$conf.low_backtransformed <- exp(me3_c$conf.low)-1
me3_c$conf.high_backtransformed <- exp(me3_c$conf.high)-1

range(me3_c$predicted)
range(me3_c$predicted_backtransformed)

# Figure 2C: Model C
# plot of Max Pasteuria prevalence vs dentifera density at max prevalence
AUC_Dent <- ggplot(me3_c, aes(x = DentMax, y = predicted_backtransformed)) +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low_backtransformed, ymax = conf.high_backtransformed), alpha = 0.2, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera2, aes(x = dentifera.at.max, y = pasteuria.auc), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote(italic("D. dentifera") ~ "Density at Max Prevalence (no." ~ "m"^-2*")"), y = bquote(atop("Integrated" ~ italic(" Pasteuria") ~ "Prevalence",  "in" ~italic("D. dentifera")~ "(prevalence x day)"))) +
  geom_text(aes(x = 2000, y= 4), label = "N.S.") +
  scale_x_log10(labels = scales::comma) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(AUC_Dent)
ggsave(here("mi-fielddata-analysis/figures/AUCPastPrev_DentiferaAtMax_predict_color.tiff"), plot = AUC_Dent, dpi = 300, width = 10, height = 10, units = "cm", compression="lzw")



################################
## MODEL D: AUC (Integrated) Pasteuria prevalence in D. dentifera vs host densities and richness at the max prevalence

### Initial Model of AUC Pasteuria prev in dentifera, with mean densities
modD_initial <- lmer(pasteuria.auc_log ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + (1|Year) + (1|Lake), 
                     data = sum.data_dentifera2) 
summary(modD_initial)
AIC(modD_initial)
vif(modD_initial)
plot(modD_initial)
qqnorm(resid(modD_initial))
qqline(resid(modD_initial)) # somewhat heavy lower tail that falls off of the QQ normal line (looks better when using data with max prev > 0)

# Reviewer requested: model with data set where max prev is > zero (no infections in dentifera) [MANUSCRIPT APPENDIX ONLY]
modD_initial2 <- lmer(pasteuria.auc_log ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + (1|Year) + (1|Lake), 
                      data = sum.data_epi_dentifera) 
summary(modD_initial2)
AIC(modD_initial2)
vif(modD_initial2)
plot(modD_initial2)
qqnorm(resid(modD_initial2))
qqline(resid(modD_initial2)) # looks better when using data with max prev > 0


# Reviewer suggested nested models
# M1: Richness only
modD1 <- lmer(pasteuria.auc_log ~ species.richness_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modD1)  # richness is not sig
AIC(modD1)

# M2: Richness + dentifera density
modD2 <- lmer(pasteuria.auc_log ~ species.richness_z + mean.dentifera.density_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modD2) # none of the factors are significant

# M3: Richness + dentifera density + pulicaria density
modD3 <- lmer(pasteuria.auc_log ~ species.richness_z + mean.dentifera.density_z + mean.pulicaria.density_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modD3) # none of the factors are significant

# M4: Richness + dentifera density + pulicaria density + retrocurva density
modD4 <- lmer(pasteuria.auc_log ~ species.richness_z + mean.dentifera.density_z + mean.pulicaria.density_z + mean.retrocurva.density_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modD4) # none of the factors are significant

# M5: Richness + dentifera density + all other host density
modD5 <- lmer(pasteuria.auc_log ~ species.richness_z + mean.dentifera.density_z + mean.other.density_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modD5) # none of the factors are significant

# compare AIC among the 3 models
AIC(modD1, modD2, modD3, modD4, modD5)  # modD1 and modD2 are equivalent, with the lowest AIC






# model section ranking by AICc using ML
options(na.action = "na.fail")
msc4 <- dredge(modD_initial, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc4)
(attr(msc4, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList4 <- get.models(msc4, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:

#### Table 1
summary(model.avg(fmList4))  # top model is only model within delta < 2

## Appendix S1 Table S8
# Export AIC table for Model D (delta AIC < 4)
msc_D <- filter(msc4, delta < 4)
write.csv(msc_D, here("mi-fielddata-analysis/results/ModelD_AIC_table.csv"), quote = F, row.names = F)


## Appendix S1 Table S9
# Top Model D
modD_top <- lmer(pasteuria.auc_log ~ mean.dentifera.density_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modD_top)


# second top model
modD_second <- lmer(pasteuria.auc_log ~ mean.dentifera.density_z + mean.retrocurva.density_z + (1|Lake), data = sum.data_dentifera2)
summary(modD_second)


# Marginal effects from top model
me4_d <- ggpredict(modD_top, c("mean.dentifera.density_z"))
plot(me4_d, add.data = T)



# recalculate original mean Dentifera density 
me4_d$x 
dent_mean <- mean(log(sum.data_dentifera$mean.dentifera.density+1)) # mean of original mean Dentifera density 
dent_sd <- sd(log(sum.data_dentifera$mean.dentifera.density+1)) # sd of original mean Dentifera density 
me4_d$DentMean_log <- t((t(me4_d$x) * dent_sd) + dent_mean)
me4_d$DentMean <- exp(me4_d$DentMean_log)
me4_d$predicted_backtransformed <- exp(me4_d$predicted)-1  #back transform the predicted and conf interval values from log+1 
me4_d$conf.low_backtransformed <- exp(me4_d$conf.low)-1
me4_d$conf.high_backtransformed <- exp(me4_d$conf.high)-1

# Figure 2D: Model D
# plot of AUC Pasteuria prevalence vs mean dentifera density
AUC_Dent2 <- ggplot(me4_d, aes(x = DentMean, y = predicted_backtransformed)) +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low_backtransformed, ymax = conf.high_backtransformed), alpha = 0.2, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera2, aes(x = mean.dentifera.density, y = pasteuria.auc), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote("Mean" ~italic("D. dentifera") ~ "Density (no." ~ "m"^-2*")"), y = bquote(atop("Integrated" ~ italic("Pasteuria") ~ "Prevalence", "in" ~italic("D. dentifera")~ "(prevalence x day)"))) +
  geom_text(aes(x = 2000, y= 4), label = "p = 0.014") +
  scale_x_log10(labels = scales::comma) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(AUC_Dent2)
ggsave(here("mi-fielddata-analysis/figures/AUCPastPrev_MeanDentifera_predict_color.tiff"), plot = AUC_Dent2, dpi = 300, width = 10, height = 10, units = "cm", compression="lzw")


# Figure 2, panels A-D
### four panel plot of MI field analyses
MI_full <- ggarrange(MaxPrev_Pulic, MaxPrev_Dent, AUC_Dent, AUC_Dent2, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(here("mi-fielddata-analysis/figures/Fig2_MI_Field_Analysis_predict_color.tiff"), plot = MI_full, dpi = 300, width = 8.1, height = 8, units = "in", compression="lzw")




############################################
## MODEL Testing if early host densities impact integrated Pasteuria prevalence (whole season epidemic size)
# Integrated Pasteuria prevalence (AUC) in D. dentifera vs mean August host densities and richness

### Initial Model of Integrated Pasteuria prevalence in dentifera, with mean AUGUST densities
modAUG2_initial <- lmer(pasteuria.auc_log ~ mean.Aug.pulicaria.density_z + mean.Aug.dentifera.density_z + mean.Aug.retrocurva.density_z + mean.Aug.richness_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)

summary(modAUG2_initial)
check_collinearity(modAUG2_initial)
plot(modAUG2_initial)
qqnorm(resid(modAUG2_initial))
qqline(resid(modAUG2_initial))


##### Model Selection and Averaging [APPROACH USED IN THE MANUSCRIPT]
# model section ranking by AICc using ML
options(na.action = "na.fail")
mscAUG2 <- dredge(modAUG2_initial, rank = "AICc", trace = TRUE, REML = FALSE)
print(mscAUG2)
(attr(mscAUG2, "rank.call"))
# Get the models (fitted by REML, as in the global model)
mscAUG_update <- filter(mscAUG2, delta < 2)
fmListAUG2 <- get.models(mscAUG2, delta < 2) # subsetting isn't working, did it manually below
fmListAUG2[-c(3:16)] # (need to remove the last model that doesn't converge)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:

#### MODEL AVERAGE 
summary(model.avg(fmListAUG2, delta < 2))  # only the top model has a delta < 2


## Appendix S1 Table SXX
# Export AIC table for Model XX (delta AIC < 4)
msc_AUG2 <- filter(mscAUG2, delta < 4)
write.csv(msc_AUG2, here("mi-fielddata-analysis/results/Model_IntegratedPrev_AugDensities_AIC_table.csv"), quote = F, row.names = F)



# Top model
modAUG2_top <- lmer(pasteuria.auc_log ~ mean.Aug.dentifera.density_z + (1|Year) + (1|Lake), data = sum.data_dentifera2)
summary(modAUG2_top)









#######################################
### Other versions of the analyses using infection density in dentifera that are NOT included in the main text
#######################################


##### Models of Pasteuria infection density in dentifera, with densities at max prevalence * Year (did not include D. dentifera density since that is included in the calculation for the infection density response variable already)
mod5 <- lmer(log(past.max.inf+1) ~ pulicaria.at.max_z + retrocurva.at.max_z + richness.at.max_z + Year + 
                pulicaria.at.max_z:Year + retrocurva.at.max_z:Year + (1|Lake), data = sum.data_dentifera) 
summary(mod5)
vif(mod5)
plot(mod5)   # pretty good fit
qqnorm(resid(mod5))
qqline(resid(mod5))
overdisp_fun(mod5)
Anova(mod5)


# model section ranking by AICc using ML
options(na.action = "na.fail")
msc5 <- dredge(mod5, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc5)
(attr(msc5, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList5 <- get.models(msc5, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList5))


mod5a <- lmer(log(past.max.inf+1) ~ pulicaria.at.max_z + (1|Lake), data = sum.data_dentifera) 
summary(mod5a)
plot(mod5a)
qqnorm(resid(mod5a))
qqline(resid(mod5a))


mod5b <- lmer(log(past.max.inf+1) ~ pulicaria.at.max_z + richness.at.max_z + (1|Lake), data = sum.data_dentifera) 
summary(mod5b)


me1 <- ggpredict(mod5a, c("pulicaria.at.max_z [all]"))
plot(me1, add.data = T)



##### Models of Pasteuria infection density in dentifera, with mean densities * Year
mod6 <- lmer(log(past.max.inf+1) ~ mean.pulicaria.density_z +  mean.retrocurva.density_z + species.richness_z + Year +
              mean.pulicaria.density_z:Year + mean.retrocurva.density_z:Year + (1|Lake), data = sum.data_dentifera) 
summary(mod6)
vif(mod6)
plot(mod6)   # data could fit the model better...it's not great
qqnorm(resid(mod6))
qqline(resid(mod6))
overdisp_fun(mod6)
Anova(mod6)


# model section ranking by AICc using ML
options(na.action = "na.fail")
msc6 <- dredge(mod6, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc6)
(attr(msc6, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList6 <- get.models(msc6, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList6))


# top model is the null model

#  second top model
mod6a <- lmer(log(past.max.inf+1) ~ mean.retrocurva.density_z + (1|Lake), data = sum.data_dentifera) 
summary(mod6a)
plot(mod6a)


me1 <- ggpredict(mod6a, c("mean.pulicaria.density_z [all]"))
plot(me1, add.data = T)

me2 <- ggpredict(mod6a, c("mean.retrocurva.density_z [all]"))
plot(me2, add.data = T)

me3 <- ggpredict(mod6a, c("species.richness_z [all]"))
plot(me3, add.data = T)



#######################################
# Assessing Maximal Prevalence and AUC for Pasteuria in Pulicaria in the field
######################################

# load cleaned parasite and host density data 2014 - 2017
data2 <- read.csv(here("mi-fielddata-analysis/data/Clean-Data-2014-2021_All-Host-Densities.csv"), stringsAsFactors = F)


# get data with both dentifera and pulicaria hosts
data_pulic <- data2 %>% 
  filter(Year != 2018 & Year != 2020 & Year != 2021) %>% 
  filter(Host.Species == "dentifera" | Host.Species == "pulicaria") %>% # get infection prevalence data for dentifera only
  filter(!is.na(Host.Density)) %>%  # remove sample dates with missing density data
  mutate(past.density.inf = pasteuria.prev * Host.Density, # add pasteuria infected density to data set
         Diluter.density = Pulicaria.density + Retrocurva.density, # diluter host densities
         OtherHost.density = Pulicaria.density + Retrocurva.density + Ceriodaphnia.density + Dubia.density +
           Parvula.density + Ambigua.density + Mendotae.density) %>%
  select(Host.Species, Lake, Year, Julian.Day, Julian, Total, starts_with("past"), Host.Density:OtherHost.density) # Select only the variables needed in cleaned data set
dim(data_pulic)
#View(data_pulic)

# diversity metrics per sampling date
data_pulic$species.richness <- rowSums(select(data_pulic, Dentifera.density:Mendotae.density) > 0)
data_pulic$shannon <- diversity(select(data_pulic, Dentifera.density:Mendotae.density), "shannon")

# calculate the max, mean, sd, and se of pasteuria prevalence for each host, lake and year
sum.data_pulic <- data_pulic %>%
  group_by(Host.Species, Lake, Year) %>%  # group by host species, lake and year
  filter(Julian < 305) %>% # remove November dates, only include July - October for max prev and mean densities
  summarize(past.max.prev = max(pasteuria.prev, na.rm = T),  # calculate max pasteuria infection prevalence
            past.mean.prev = mean(pasteuria.prev, na.rm = T),# calculate mean pasteuria infection prevalence
            past.sd.prev = sd(pasteuria.prev, na.rm = T),    # calculate sd pasteuria infection prevalence
            past.se.prev = sd(pasteuria.prev, na.rm = T)/sqrt(length(pasteuria.prev)),  # calculate standard error pasteuria infection prevalence
            past.max.inf = max(past.density.inf, na.rm = T),  # maximum pf pasteuria infection density
            mean.total.count = mean(Total, na.rm = T),   # mean of the total dentifera hosts counted
            mean.density = mean(Host.Density, na.rm = T),    # mean host density of each species
            mean.tot.density = mean(Total.Density, na.rm = T),# mean total density of all hosts
            mean.dentifera.density = mean(Dentifera.density, na.rm = T),# mean dentifera density
            mean.pulicaria.density = mean(Pulicaria.density, na.rm = T),# mean pulicaria density
            mean.retrocurva.density = mean(Retrocurva.density, na.rm = T), # mean retrocurva density
            mean.ceriodaphnia.density = mean(Ceriodaphnia.density, na.rm = T), # mean ceriodaphnia density
            mean.dubia.density = mean(Dubia.density, na.rm = T), # mean Dubia density
            mean.parvula.density = mean(Parvula.density, na.rm = T), # mean parvula density
            mean.ambigua.density = mean(Ambigua.density, na.rm = T), # mean Ambigua density
            mean.mendotae.density = mean(Mendotae.density, na.rm = T), # mean Mendotae density
            mean.diluter.density = mean(Diluter.density, na.rm = T),  # mean of diluter hosts (pulicaria + retrocurva)
            mean.other.density = mean(OtherHost.density, na.rm = T), # mean of all other hotsts (non-dentifera)
            mean.richness = mean(species.richness, na.rm = T),       # mean species richness
            mean.shannon = mean(shannon, na.rm = T)) %>%             # mean shannon diversity index
  ungroup()
# note that the error produced is for one lake year (2015 Whitemore) that does not have any records of pasteuria in dentifera that were greater than 20 individuals counted, therefore prevalence could not be calculated for those dates.
#View(sum.data_pulic)
dim(sum.data_pulic)

sum.data_pulic$Year <- as.character(sum.data_pulic$Year)
sum.data_pulic$mean.total.count[ is.nan(sum.data_pulic$mean.total.count)] <- 0


# add species richness and simpson diversity index to data set
sum.data_pulic$species.richness <- rowSums(select(sum.data_pulic, mean.dentifera.density:mean.mendotae.density) > 0)


unique(sum.data_pulic$Lake)


# load data with area under the curve calculated, and host densities on the date of the max prevalence
auc_data2 <- read.csv(here("mi-fielddata-analysis/data/auc_14to21_prev.csv"), stringsAsFactors = F)
auc_data2$Year <- as.character(auc_data$Year)


# select only Pasteuria data and prep to join with the sum.data 
auc_data_past_pulic <- auc_data2 %>%
  filter(Parasite.Species == "pasteuria.prev" & Year != 2020 & Year != 2021) %>%
  filter(Host.Species == "dentifera" | Host.Species == "pulicaria" ) %>%
  arrange(Host.Species, Lake) %>%
  dplyr::rename(pasteuria.auc = AUC.prev2) %>%  # Use AUC.prev2 instead of AUC.prev because it is more accurately calculated
  select(Host.Species, Lake, Year, pasteuria.auc, Max.Prevalence, Count.At.Max:shannon.at.max)
#View(auc_data_past_pulic)

str(auc_data_past_pulic)
str(sum.data_pulic)

sum.data_pulic  <- left_join(sum.data_pulic, auc_data_past_pulic)

sum.data_pulic$past.max.prev[sum.data_pulic$past.max.prev == "-Inf"] <- NA

# pivot wider to get max prev and AUC of pasteuria in dentifera and pulicaria in separate columns
data_dent.pulic <- sum.data_pulic %>%
  select(Lake, Year, Host.Species, past.max.prev, pasteuria.auc) %>%
  tidyr::pivot_wider(names_from = "Host.Species", values_from = c(past.max.prev, pasteuria.auc), names_sep = ".") %>%
  mutate(missing.prev = if_else(is.na(past.max.prev.dentifera),"Pulicaria only", if_else(is.na(past.max.prev.pulicaria), "Dentifera only", "Both present")),
         missing.auc = if_else(is.na(pasteuria.auc.dentifera),"Pulicaria only", if_else(is.na(pasteuria.auc.pulicaria), "Dentifera only", "Both present")))  # make a categorical variable to indicate if there were no dentifera or pulicaria to calculate a max prev or AUC
View(data_dent.pulic)

# replace NAs with Zeros so that points will show up for both host species, but I plan to use the missing categorical variables I just made
# to indicate which points are missing values because none of that host species was detected.
data_dent.pulic[is.na(data_dent.pulic)] <- 0

# make year a factor to work in the figures below
data_dent.pulic$Year <- as.character(data_dent.pulic$Year)

data_dent.pulic <- as.data.frame(data_dent.pulic)


slope <- lm(past.max.prev.pulicaria ~ past.max.prev.dentifera, data = data_dent.pulic)
summary(slope)

## Appendix S1: Figure S1
# Make a figure that shows a 1:1 line and paired pasteuria prevalence in pulicaria and dentifera
prev_dent_pulic <- ggplot(data = data_dent.pulic, aes(x = past.max.prev.dentifera, y = past.max.prev.pulicaria)) +
  geom_jitter(aes(shape = missing.prev), alpha = 0.5, size = 2, height = 0.0005) +
  geom_smooth(method = "lm", se = FALSE, col = "black") +
  geom_abline(slope = 1, intercept = 0, size = 1, linetype = "dashed", ylim(0,0.07)) +
  labs(x = bquote("Maximal " ~ italic("Pasteuria") ~ "Prevalence in" ~ italic("D. dentifera")), y = bquote("Maximal " ~ italic("Pasteuria") ~ "Prevalence in" ~ italic("D. pulicaria")),
       shape = "Host Species \nDetected") +
  scale_y_continuous(labels = scales::percent, limits = c(-0.001,0.07)) +
  scale_x_continuous(labels = scales::percent) +
  annotate("text", x = 0.04, y = 0.065, label = "1:1") +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
prev_dent_pulic
ggsave(here("mi-fielddata-analysis/figures/MaxPastPrev_PulicariaVsDentifera.tiff"), plot = prev_dent_pulic, dpi = 300, width = 5, height = 4, units = "in", compression="lzw")





##############################################################
# Make figures of prevalence and host densities over time
##############################################################

# full time series data for all hosts 
head(data_all)
summary(data_all)

# data for just dentifera
head(data)
summary(data)

data$Year_factor <- as.factor(data$Year)

# Figure of Pasteuria prevalence in dentifera for all lakes and years (2014-2017)
data_dent_prev_all_years <- filter(data, Host.Species == "dentifera", Julian < 305, Julian > 190)

pastprev_dent_time_allyears <- ggplot(data_dent_prev_all_years, aes(x = Julian, y = pasteuria.prev, color = Year_factor)) +
  labs(title = "Pasteuria prevalence in D. dentifera through time", y = "Pasteuria Prevalence in D. dentifera", x = "Julian Day") +
  geom_line() +
  geom_point(size = 1) +
  geom_line(aes(y = 0.01), linetype = "dashed", color = "gray") +
  facet_wrap( ~ Lake, scales = "free_y") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
print(pastprev_dent_time_allyears)
ggsave("mi-fielddata-analysis/figures/PastPrev_Dentifera_JulianDay_AllYears.tiff", plot = pastprev_dent_time_allyears, dpi = 300, width = 7.5, height = 6, units = "in", compression="lzw")





# Figures of Pasteuria prevalence in dentifera, retrocurva, and pulicaria for all lakes and separated by each year
data_allhost_prev_all_years <- filter(data_all, Host.Species == "dentifera" | Host.Species == "retrocurva" | Host.Species == "pulicaria", Year < 2018, Julian < 305, Julian > 190)
unique(data_allhost_prev_all_years$Host.Species)
range(data_allhost_prev_all_years$Julian)

pastprev_allhost_time_allyears <- ggplot(data_allhost_prev_all_years, aes(x = Julian, y = pasteuria.prev, color = Host.Species)) +
  labs(y = "Pasteuria Prevalence", x = "Julian Day") +
  geom_line() +
  geom_point(size = 1) +
  geom_line(aes(y = 0.01), linetype = "dashed", color = "black") +
  facet_grid(Year ~ Lake, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")
print(pastprev_allhost_time_allyears)
ggsave("mi-fielddata-analysis/figures/PastPrev_AllHost_JulianDay_AllYears.tiff", plot = pastprev_allhost_time_allyears, dpi = 300, width = 11, height = 6.5, units = "in", compression="lzw")





# Figures of host density for dentifera, retrocurva, and pulicaria for all lakes and separated by each year

density_allhost_time_allyears <- ggplot(data_allhost_prev_all_years, aes(x = Julian, y = Host.Density+1, color = Host.Species)) +
  labs(y = "log Host Density + 1", x = "Julian Day") +
  geom_line() +
  geom_point(size = 1) +
  geom_line(aes(y = 1), linetype = "dashed", color = "black") +
  facet_grid(Year ~ Lake, scales = "free_y") +
  scale_y_log10() +
  scale_color_brewer(palette = "Dark2") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "bottom")
print(density_allhost_time_allyears)
ggsave("mi-fielddata-analysis/figures/Density_AllHost_JulianDay_AllYears.tiff", plot = density_allhost_time_allyears, dpi = 300, width = 11, height = 6.5, units = "in", compression="lzw")







