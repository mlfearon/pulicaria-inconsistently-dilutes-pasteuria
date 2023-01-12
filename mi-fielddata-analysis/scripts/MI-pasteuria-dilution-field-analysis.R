# Analysis of field Pasteuria data and dilution effect for "Mixed evidence for a dilution effect in Daphnia communities infected by a bacterial parasite"


# Submitted to: Oecologia


# Code written by Michelle Fearon
# Last updated: Dec 17, 2022

## This code produces four models to test how host density, species richness, and year correlate with
## either maximal Pasteuria prevalence (models A & B) or area under the prevalence curve (models C & D) in D. dentifera, 
## with two models for each response variable. Models A & C included the host densities and species richness 
## from the date of max prevalence, while models B & D included the mean host densities and total species richness 
## for each lake and year combination. We then used model selection to determine the best version of each model A-D.


#libraries
library(ggplot2)
library(dplyr)
library(vegan)
library(DescTools)
library(lme4)
library(lmerTest)
library(car)
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
data <- read.csv(here("mi-fielddata-analysis/data/Clean-Data-2014-2021_All-Host-Densities.csv"), stringsAsFactors = F)
# Total = total number of daphnia hosts for that species/site/date
# .inf = the total number of infections for each parasite as a column
# .prev = the calculated prevalence of each parasite (only when the total host density for that species/site/date > 20)
# host.density = total area density of the host species for that site/date
# .den = density of each host species for that site/date
# total.density = summed density of all host species for that site/date


# get data to only include dentifera for main analyeses
data <- data %>% 
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


sum.data$Year <- as.character(sum.data$Year)
sum.data$mean.total.count[ is.nan(sum.data$mean.total.count)] <- 0



# add species richness and simpson diversity index to data set
sum.data$species.richness <- rowSums(select(sum.data, mean.dentifera.density:mean.mendotae.density) > 0)


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



# make year a factor to work in the figures below
sum.data$Year <- as.character(sum.data$Year)
sum.data$Host.Species <- as.factor(sum.data$Host.Species)


# filter out the lake year that doesn't have a valid max prevalence because there was not enough dentifera detected to calculate prevalence
sum.data <- filter(sum.data, !is.na(Max.Prevalence)) # removed Whitmore in 2015 row.
sum.data_dentifera <- as.data.frame(sum.data)


# Center and scale all continuous variables, densities get log + 1 transformation
sum.data_dentifera$mean.tot.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.tot.density)))
sum.data_dentifera$mean.dentifera.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.dentifera.density+1)))
sum.data_dentifera$mean.pulicaria.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.pulicaria.density+1)))
sum.data_dentifera$mean.retrocurva.density_z <- as.numeric(scale(log(sum.data_dentifera$mean.retrocurva.density+1)))
sum.data_dentifera$species.richness_z <- as.numeric(scale(sum.data_dentifera$species.richness))
sum.data_dentifera$mean.richness_z <- as.numeric(scale(sum.data_dentifera$mean.richness))
sum.data_dentifera$mean.shannon_z <- as.numeric(scale(sum.data_dentifera$mean.shannon))
sum.data_dentifera$dentifera.at.max_z <- as.numeric(scale(log(sum.data_dentifera$dentifera.at.max+1)))
sum.data_dentifera$pulicaria.at.max_z <- as.numeric(scale(log(sum.data_dentifera$pulicaria.at.max+1)))
sum.data_dentifera$retrocurva.at.max_z <- as.numeric(scale(log(sum.data_dentifera$retrocurva.at.max+1)))
sum.data_dentifera$diluter.at.max_z <- as.numeric(scale(log(sum.data_dentifera$diluter.at.max+1)))
sum.data_dentifera$otherhost.at.max_z <- as.numeric(scale(log(sum.data_dentifera$otherhost.at.max+1)))
sum.data_dentifera$richness.at.max_z <- as.numeric(scale(sum.data_dentifera$richness.at.max))
sum.data_dentifera$shannon.at.max_z <- as.numeric(scale(sum.data_dentifera$shannon.at.max))
sum.data_dentifera$pasteuria.auc_log <- as.numeric(log(sum.data_dentifera$pasteuria.auc+1))

#include observation level random effect
sum.data_dentifera$OLRE <- seq_len(nrow(sum.data_dentifera))


# Test the correlation between Pasteuria Maximum prevalence and Area Under the prevalence curve
cor.test(sum.data_dentifera$past.max.prev, sum.data_dentifera$pasteuria.auc)
cor.test(sum.data_dentifera$past.max.prev, sum.data_dentifera$pasteuria.auc_log)


cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.tot.density)
cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.dentifera.density)
cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.retrocurva.density)
cor.test(sum.data_dentifera$species.richness, sum.data_dentifera$mean.pulicaria.density)



## MODEL A
###  Initial Model of MAX Pasteuria prev in dentifera, with densities at max prev * Year 
mod <- glmer(past.max.prev ~ pulicaria.at.max_z + dentifera.at.max_z + retrocurva.at.max_z + richness.at.max_z +
               Year + pulicaria.at.max_z:Year + dentifera.at.max_z:Year + retrocurva.at.max_z:Year + (1|Lake)+ (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod)
vif(mod)
overdisp_fun(mod)
Anova(mod)




# model section ranking by AICc using ML
options(na.action = "na.fail")
msc <- dredge(mod, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc)
(attr(msc, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList <- get.models(msc, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList))

# Export AIC table for Model A (delta AIC < 4)
msc_A <- filter(msc, delta < 4)
write.csv(msc_A, here("mi-fielddata-analysis/results/ModelA_AIC_table.csv"), quote = F, row.names = F)


# top model A
modA <- glmer(past.max.prev ~ pulicaria.at.max_z + (1|Lake) + (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(modA)
overdisp_fun(modA)

# Marginal effects from top model
me1_a <- ggpredict(modA, c("pulicaria.at.max_z [all]"))
plot(me1_a, add.data = T)


# recalculate original Pulicaria density at Max 
me1_a$x 
pulic_mean <- mean(log(sum.data_dentifera$pulicaria.at.max+1)) # mean of original richness from disease data set
pulic_sd <- sd(log(sum.data_dentifera$pulicaria.at.max+1)) # sd of original richness from disease data set
me1_a$PulicMax_log <- t((t(me1_a$x) * pulic_sd) + pulic_mean)
me1_a$PulicMax <- exp(me1_a$PulicMax_log)

display.brewer.pal(3, "Dark2")
brewer.pal(3, "Dark2")

# Figure 2A: Model A
# plot of Max Pasteuria prevalence vs pulicaria density at max prevalence
MaxPrev_Pulic <- ggplot(me1_a, aes(x = PulicMax, y = predicted)) +
  geom_line(size = 1,  color = "#D95F02") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#D95F02") +
  geom_jitter(data = sum.data_dentifera, aes(x = pulicaria.at.max+1, y = past.max.prev), size = 2, alpha = 0.4, width = 0.05, color = "#D95F02", shape = 17) +
  labs(x = bquote(italic("D. pulicaria") ~ "Density at Max Prevalence (no./" ~ "M"^2*")"), y = bquote(atop( "","Maximal " ~ italic("Pasteuria") ~ "Prevalence"))) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(x = 10, y= 0.1), label = "p = 0.097") +
  scale_x_log10(labels = scales::comma) +
  theme_classic()+
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(MaxPrev_Pulic)
ggsave(here("mi-fielddata-analysis/figures/MaxPastPrev_PulicariaAtMax_predict_color.tiff"), plot = MaxPrev_Pulic, dpi = 300, width = 10, height = 10, units = "cm", compression="lzw")


# second top model is the NULL model


# third top model
modA_3 <- glmer(past.max.prev ~ dentifera.at.max_z + (1|Lake) + (1|OLRE), 
              family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(modA_3)




## MODEL B
### Model of MAX Pasteuria prev in dentifera, with mean densities * Year 
mod2 <- glmer(past.max.prev ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + Year +
               mean.pulicaria.density_z:Year + mean.retrocurva.density_z:Year + mean.dentifera.density_z:Year + (1|Lake) + (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2)
vif(mod2)
overdisp_fun(mod2)

# model section ranking by AICc using ML
options(na.action = "na.fail")
msc2 <- dredge(mod2, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc2)
(attr(msc2, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList2 <- get.models(msc2, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList2))



# Export AIC table for Model B (delta AIC < 4)
msc_B <- filter(msc2, delta < 4)
write.csv(msc_B, here("mi-fielddata-analysis/results/ModelB_AIC_table.csv"), quote = F, row.names = F)



# top model B
mod2B <- glmer(past.max.prev ~ mean.dentifera.density_z + (1|Lake) + (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2B)
overdisp_fun(mod2B)

# second top model
mod2B_2 <- glmer(past.max.prev ~ mean.dentifera.density_z + Year + (1|Lake) + (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2B_2)
vif(mod2B_2)
overdisp_fun(mod2B_2)
Anova(mod2B_2)



# Marginal effects from top model
me2_b <- ggpredict(mod2B, c("mean.dentifera.density_z [all]"))
plot(me2_b, add.data = T)



# recalculate original mean Dentifera density
me2_b$x 
dent_mean <- mean(log(sum.data_dentifera$mean.dentifera.density+1)) # mean of original richness from disease data set
dent_sd <- sd(log(sum.data_dentifera$mean.dentifera.density+1)) # sd of original richness from disease data set
me2_b$DentMean_log <- t((t(me2_b$x) * dent_sd) + dent_mean)
me2_b$DentMean <- exp(me2_b$DentMean_log)



# plot of Max Pasteuria prevalence vs mean dentifera density
MaxPrev_Dent <- ggplot(me2_b, aes(x = DentMean, y = predicted)) +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera, aes(x = mean.dentifera.density+1, y = past.max.prev), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote("Mean" ~italic("D. dentifera") ~ "Density (no./" ~ "M"^2*")"), y = bquote(atop( "","Maximal " ~ italic("Pasteuria") ~ "Prevalence"))) +
  scale_y_continuous(labels = scales::percent) +
  geom_text(aes(x = 1000, y= 0.1), label = "p = 0.018") +
  scale_x_log10(labels = scales::comma) +
  theme_classic()+
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(MaxPrev_Dent)
ggsave(here("mi-fielddata-analysis/figures/MaxPastPrev_MeanDentifera_predict_color.tiff"), plot = MaxPrev_Dent, dpi = 300, width = 10, height = 10, units = "cm", compression="lzw")




## MODEL C
##### Models of AUC Pasteuria prev in dentifera, with densities at max prev * Year 
sum.data_dentifera <- filter(sum.data_dentifera, !is.na(pasteuria.auc))  # remove one instance of NA where there was only a single date where enough dentifera were counted during the season, so a max prev could be caluclated but not AUC

mod3 <- lmer(pasteuria.auc_log ~ pulicaria.at.max_z + dentifera.at.max_z + retrocurva.at.max_z + richness.at.max_z + Year + 
               pulicaria.at.max_z:Year + dentifera.at.max_z:Year + retrocurva.at.max_z:Year + (1|Lake), data = sum.data_dentifera) 
summary(mod3)
vif(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqline(resid(mod3)) # somewhat heavy lower tail that falls off of the QQ normal line
overdisp_fun(mod3)
Anova(mod3)


# model section ranking by AICc using ML
options(na.action = "na.fail")
msc3 <- dredge(mod3, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc3)
(attr(msc3, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList3 <- get.models(msc3, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList3))


# Export AIC table for Model C (delta AIC < 4)
msc_C <- filter(msc3, delta < 4)
write.csv(msc_C, here("mi-fielddata-analysis/results/ModelC_AIC_table.csv"), quote = F, row.names = F)



# top model C
mod3C <- lmer(pasteuria.auc_log ~ dentifera.at.max_z + (1|Lake), data = sum.data_dentifera)
summary(mod3C)


# second top model is the Null model



# third top model
mod3C_2 <- lmer(pasteuria.auc_log ~ dentifera.at.max_z + retrocurva.at.max_z + (1|Lake), data = sum.data_dentifera)
summary(mod3C_2)
vif(mod3C_2)
overdisp_fun(mod3C_2)
Anova(mod3C_2)



# Marginal effects from top model
me3_c <- ggpredict(mod3C, c("dentifera.at.max_z"))
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

# plot of Max Pasteuria prevalence vs dentifera density at max prevalence
AUC_Dent <- ggplot(me3_c, aes(x = DentMax, y = predicted_backtransformed)) +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low_backtransformed, ymax = conf.high_backtransformed), alpha = 0.2, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera, aes(x = dentifera.at.max+1, y = pasteuria.auc), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote(italic("D. dentifera") ~ "Density at Max Prevalence (no./" ~ "M"^2*")"), y = bquote(atop("Integrated" ~ italic(" Pasteuria") ~ "Prevalence", "(prevalence x day)"))) +
  geom_text(aes(x = 2000, y= 4), label = "p = 0.06") +
  scale_x_log10(labels = scales::comma) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(AUC_Dent)
ggsave(here("mi-fielddata-analysis/figures/AUCPastPrev_DentiferaAtMax_predict_color.tiff"), plot = AUC_Dent, dpi = 300, width = 12, height = 10, units = "cm", compression="lzw")




## MODEL D
##### Models of AUC Pasteuria prev in dentifera, with mean densities * Year 
mod4 <- lmer(pasteuria.auc_log ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + Year +
               mean.pulicaria.density_z:Year + mean.dentifera.density_z:Year + mean.retrocurva.density_z:Year + (1|Lake), data = sum.data_dentifera)
summary(mod4)
vif(mod4)
plot(mod4)
qqnorm(resid(mod4))
qqline(resid(mod4))
overdisp_fun(mod4)
Anova(mod4)




# model section ranking by AICc using ML
options(na.action = "na.fail")
msc4 <- dredge(mod4, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc4)
(attr(msc4, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList4 <- get.models(msc4, delta < 2)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList4))


# Export AIC table for Model D (delta AIC < 4)
msc_D <- filter(msc4, delta < 4)
write.csv(msc_D, here("mi-fielddata-analysis/results/ModelD_AIC_table.csv"), quote = F, row.names = F)



# top model D
mod4D <- lmer(pasteuria.auc_log ~ mean.dentifera.density_z + (1|Lake), data = sum.data_dentifera)
summary(mod4D)




# second top model
mod4D_2 <- lmer(pasteuria.auc_log ~ mean.dentifera.density_z + mean.retrocurva.density_z + (1|Lake), data = sum.data_dentifera)
summary(mod4D_2)



# Marginal effects from top model
me4_d <- ggpredict(mod4D, c("mean.dentifera.density_z"))
plot(me4_d, add.data = T)



# recalculate original mean Dentifera density 
me4_d$x 
dent_mean <- mean(log(sum.data_dentifera$mean.dentifera.density+1)) # mean of original richness from disease data set
dent_sd <- sd(log(sum.data_dentifera$mean.dentifera.density+1)) # sd of original richness from disease data set
me4_d$DentMean_log <- t((t(me4_d$x) * dent_sd) + dent_mean)
me4_d$DentMean <- exp(me4_d$DentMean_log)
me4_d$predicted_backtransformed <- exp(me4_d$predicted)-1  #back transform the predicted and conf interval values from log+1 
me4_d$conf.low_backtransformed <- exp(me4_d$conf.low)-1
me4_d$conf.high_backtransformed <- exp(me4_d$conf.high)-1


# plot of AUC Pasteuria prevalence vs mean dentifera density
AUC_Dent2 <- ggplot(me4_d, aes(x = DentMean, y = predicted_backtransformed)) +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low_backtransformed, ymax = conf.high_backtransformed), alpha = 0.2, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera, aes(x = mean.dentifera.density+1, y = pasteuria.auc), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote("Mean" ~italic("D. dentifera") ~ "Density (no./" ~ "M"^2*")"), y = bquote(atop("Integrated" ~ italic("Pasteuria") ~ "Prevalence", "(prevalence x day)"))) +
  geom_text(aes(x = 2000, y= 4), label = "p = 0.014") +
  scale_x_log10(labels = scales::comma) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(AUC_Dent2)
ggsave(here("mi-fielddata-analysis/figures/AUCPastPrev_MeanDentifera_predict_color.tiff"), plot = AUC_Dent2, dpi = 300, width = 10, height = 10, units = "cm", compression="lzw")



### four panel plot of MI field analyses
MI_full <- ggarrange(MaxPrev_Pulic, MaxPrev_Dent, AUC_Dent, AUC_Dent2, ncol = 2, nrow = 2, labels = c("A", "B", "C", "D"))
ggsave(here("mi-fielddata-analysis/figures/MI_Field_Analysis_predict_color.tiff"), plot = MI_full, dpi = 300, width = 8, height = 8, units = "in", compression="lzw")











### TO CUT?????


##### Models of Pasteuria infection density in dentifera, with densities at max prevalence * Year
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
















##### Models of MEAN Pasteuria prev in dentifera

mod10 <- glmer(past.mean.prev ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + Year +
               mean.pulicaria.density_z:Year + mean.retrocurva.density_z:Year + mean.dentifera.density_z:Year + (1|OLRE), 
             family = "binomial", weights = mean.total.count, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5))) 
summary(mod10)
vif(mod10)
overdisp_fun(mod10)
Anova(mod10)



# model section ranking by AICc using ML
options(na.action = "na.fail")
msc10 <- dredge(mod10, rank = "AICc", trace = TRUE, REML = FALSE)
print(msc10)
(attr(msc10, "rank.call"))
# Get the models (fitted by REML, as in the global model)
fmList10 <- get.models(msc10, delta < 4)
# Because the models originate from 'dredge(..., rank = AICc, REML = FALSE)',
# the default weights in 'model.avg' are ML based:
summary(model.avg(fmList10))


# top model
mod10a <- glmer(past.mean.prev ~ mean.dentifera.density_z + (1|OLRE), 
               family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod10a)


# second top model
mod10b <- glmer(past.max.prev ~ mean.dentifera.density_z + mean.retrocurva.density_z + (1|OLRE), 
               family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod10b)
vif(mod10b)
overdisp_fun(mod10b)
Anova(mod10b)

me1 <- ggpredict(mod10, c("mean.pulicaria.density_z [all]", "Year"))
plot(me1, add.data = T)

me2 <- ggpredict(mod10, c("mean.retrocurva.density_z [all]"))
plot(me2, add.data = T)

me3 <- ggpredict(mod10, c("mean.dentifera.density_z [all]"))
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
