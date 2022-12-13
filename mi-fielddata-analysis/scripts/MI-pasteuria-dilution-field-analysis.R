# Analysis of field Pasteuria data and dilution effect for "Mixed evidence for a dilution effect in Daphnia communities infected by a bacterial parasite"


# Submitted to: Oecologia


# Code written by Michelle Fearon
# Last updated: Dec 13, 2022

## This code produces four models to test how host density, species richness, and year correlate with
## either maximal Pasteuria prevalence (models A & B) or area under the prevalence curve (models C & D) in D. dentifera, 
## with two models for each response variable. Models A & C included the host densities and species richness 
## from the date of max prevalence, while models B & D included the mean host densities and total species richness 
## for each lake and year combination. We then used model selection to determine the best version of each model A-D.


# set the path to the script relative to the project root directory
here::i_am("mi-fielddata-analysis/scripts/MI-pasteuria-dilution-field-analysis.R")


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
library(vegan)
library(pscl)
library(grid)
library(gridExtra)
library(reshape)
library(ggiraph)
library(ggiraphExtra)
library(moonBook)
library(RColorBrewer)

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


data <- data %>% 
  filter(Year != 2018 & Year != 2020 & Year != 2021) %>% 
  filter(Host.Species == "dentifera") %>% # get infection prevalence data for dentifera only
  filter(!is.na(Host.Density)) %>%  # remove sample dates with missing density data
  mutate(past.density.inf = pasteuria.prev * Host.Density, # add pasteuria infected density to data set
         Diluter.density = Pulicaria.density + Retrocurva.density, # diluter host densities
         OtherHost.density = Pulicaria.density + Retrocurva.density + Ceriodaphnia.density + Dubia.density +
           Parvula.density + Ambigua.density + Mendotae.density) %>%
  select(Host.Species, Lake, Year, Julian.Day, Julian, Total, starts_with("past"), Host.Density:OtherHost.density) # Select only the variables needed in cleaned data set


# diversity metrics per sampling date
data$species.richness <- rowSums(select(data, Dentifera.density:Mendotae.density) > 0)
data$shannon <- diversity(select(data, Dentifera.density:Mendotae.density), "shannon")

# check data set
str(data)
head(data)
dim(data)
View(data)

range(data$Julian)



######################################################################
### Calculate the maximal and mean infection prevalence per host/lake/year
######################################################################


# calculate the max, mean, sd, and se of pasteuria prevalence for each host, lake and year
sum.data <- data %>%
  group_by(Host.Species, Lake, Year) %>%  # group by host species, lake and year
  filter(Julian < 305) %>% # remove November dates, only include July - October for max prev and mean densities
  summarize(past.max.prev = max(pasteuria.prev),  # calculate max pasteuria infection prevalence
            past.mean.prev = mean(pasteuria.prev),# calculate mean pasteuria infection prevalence
            past.sd.prev = sd(pasteuria.prev),    # calculate sd pasteuria infection prevalence
            past.se.prev = sd(pasteuria.prev)/sqrt(length(pasteuria.prev)),  # calculate standard error pasteuria infection prevalence
            past.max.inf = max(past.density.inf),  # maximum pf pasteuria infection density
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

View(sum.data)
dim(sum.data)


sum.data$Year <- as.character(sum.data$Year)
sum.data$mean.total.count[ is.nan(sum.data$mean.total.count)] <- 0



# add species richness and simpson diversity index to data set
sum.data$species.richness <- rowSums(select(sum.data, mean.dentifera.density:mean.mendotae.density) > 0)


unique(sum.data$Lake)

#########################
# Area Under the Curve
##########################

# load data with area under the curve calculated
auc_data <- read.csv(here("mi-fielddata-analysis/data/auc_14to21_prev.csv"), stringsAsFactors = F)
auc_data$Year <- as.character(auc_data$Year)


# select only Pasteuria data and prep to join with 
auc_data_past <- auc_data %>%
  filter(Parasite.Species == "pasteuria.prev") %>%
  arrange(Host.Species, Lake) %>%
  dplyr::rename(pasteuria.auc = AUC.prev2) %>%
  select(Host.Species, Lake, Year, pasteuria.auc, Max.Prevalence, Count.At.Max:shannon.at.max)


str(auc_data_past)
str(sum.data)

sum.data  <- left_join(sum.data, auc_data_past)



# make year a factor to work in the figures below
sum.data$Year <- as.character(sum.data$Year)
sum.data$Host.Species <- as.factor(sum.data$Host.Species)


#filter to get only pasteuria prevalence in dentifera
sum.data_dentifera <- sum.data %>%
  filter(Host.Species == "dentifera")
View(sum.data_dentifera)
str(sum.data_dentifera)
sum.data_dentifera <- as.data.frame(sum.data_dentifera)

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



hist(sum.data_dentifera$richness.at.max)
sum(sum.data_dentifera$species.richness == 6)

#include observation level random effect
sum.data_dentifera$OLRE <- seq_len(nrow(sum.data_dentifera))



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


# top model A
moda <- glmer(past.max.prev ~ pulicaria.at.max_z + (1|Lake) + (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(moda)


# Marginal effects from top model
me1a <- ggpredict(moda, c("pulicaria.at.max_z [all]"))
plot(me1a, add.data = T)


# recalculate original Pulicaria density at Max 
me1a$x 
pulic_mean <- mean(log(sum.data_dentifera$pulicaria.at.max+1)) # mean of original richness from disease data set
pulic_sd <- sd(log(sum.data_dentifera$pulicaria.at.max+1)) # sd of original richness from disease data set
me1a$PulicMax_log <- t((t(me1a$x) * pulic_sd) + pulic_mean)
me1a$PulicMax <- exp(me1a$PulicMax_log)

display.brewer.pal(3, "Dark2")
brewer.pal(3, "Dark2")

# plot of Max Pasteuria prevalence vs pulicaria density at max prevalence
MaxPrev_Pulic <- ggplot(me1a, aes(x = PulicMax, y = predicted)) +
  labs(tag = "A") +
  geom_line(size = 1,  color = "#D95F02") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, outline.type = NULL, fill = "#D95F02") +
  geom_jitter(data = sum.data_dentifera, aes(x = pulicaria.at.max+1, y = past.max.prev), size = 2, alpha = 0.4, width = 0.05, color = "#D95F02", shape = 17) +
  labs(x = bquote(italic("D. pulicaria") ~ "Density at Max Prevalence (no./" ~ "M"^2 ~")"), y = bquote(atop( "","Maximal " ~ italic("Pasteuria") ~ "Prevalence (%)"))) +
  #ylim(0,1) +
  geom_text(aes(x = 50, y= 0.1), label = "p = 0.06") +
  scale_x_log10(labels = scales::comma) +
  theme_classic()+
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(MaxPrev_Pulic)
ggsave("MaxPastPrev_PulicariaAtMax_predict_color.png", plot = MaxPrev_Pulic, dpi = 300, width = 10, height = 10, units = "cm")





# third top model
modb <- glmer(past.max.prev ~ pulicaria.at.max_z + retrocurva.at.max_z + (1|Lake) + (1|OLRE), 
              family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(modb)





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


# top model
mod2a <- glmer(past.max.prev ~ mean.dentifera.density_z + (1|Lake) + (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2a)


# second top model
mod2b <- glmer(past.max.prev ~ mean.dentifera.density_z + Year + (1|Lake) + (1|OLRE), 
             family = "binomial", weights = Count.At.Max, data = sum.data_dentifera, control=glmerControl(optimizer="bobyqa", optCtrl=list(maxfun=2e5)))
summary(mod2b)
vif(mod2b)
overdisp_fun(mod2b)
Anova(mod2b)



# Marginal effects from top model
me2a <- ggpredict(mod2a, c("mean.dentifera.density_z [all]"))
plot(me2a, add.data = T)



# recalculate original mean Dentifera density
me2a$x 
dent_mean <- mean(log(sum.data_dentifera$mean.dentifera.density+1)) # mean of original richness from disease data set
dent_sd <- sd(log(sum.data_dentifera$mean.dentifera.density+1)) # sd of original richness from disease data set
me2a$DentMean_log <- t((t(me2a$x) * dent_sd) + dent_mean)
me2a$DentMean <- exp(me2a$DentMean_log)



# plot of Max Pasteuria prevalence vs mean dentifera density
MaxPrev_Dent <- ggplot(me2a, aes(x = DentMean, y = predicted)) +
  labs(tag = "B") +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, outline.type = NULL, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera, aes(x = mean.dentifera.density+1, y = past.max.prev), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote("Mean" ~italic("D. dentifera") ~ "Density (no./" ~ "M"^2 ~")"), y = bquote(atop( "","Maximal " ~ italic("Pasteuria") ~ "Prevalence (%)"))) +
  #ylim(0,1) +
  geom_text(aes(x = 1000, y= 0.1), label = "p = 0.039") +
  scale_x_log10(labels = scales::comma) +
  theme_classic()+
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(MaxPrev_Dent)
ggsave("MaxPastPrev_MeanDentifera_predict_color.png", plot = MaxPrev_Dent, dpi = 300, width = 10, height = 10, units = "cm")









##### Models of AUC Pasteuria prev in dentifera, with densities at max prev * Year 
mod3 <- lmer(log(pasteuria.auc+1) ~ pulicaria.at.max_z + dentifera.at.max_z + retrocurva.at.max_z + richness.at.max_z + Year + 
               pulicaria.at.max_z:Year + dentifera.at.max_z:Year + retrocurva.at.max_z:Year + (1|Lake), data = sum.data_dentifera) 
summary(mod3)
vif(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqline(resid(mod3))
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


# top model
mod3a <- lmer(log(pasteuria.auc+1) ~ dentifera.at.max_z + (1|Lake), data = sum.data_dentifera)
summary(mod3a)


# second top model
mod3b <- lmer(log(pasteuria.auc+1) ~ dentifera.at.max_z + retrocurva.at.max_z + (1|Lake), data = sum.data_dentifera)
summary(mod3b)
vif(mod3b)
overdisp_fun(mod3b)
Anova(mod3b)



# Marginal effects from top model
me3a <- ggpredict(mod3a, c("dentifera.at.max_z [all]"))
plot(me3a, add.data = T)



# recalculate original Dentifera density at max prevalence
me3a$x 
dent_mean <- mean(log(sum.data_dentifera$dentifera.at.max+1)) # mean of original richness from disease data set
dent_sd <- sd(log(sum.data_dentifera$dentifera.at.max+1)) # sd of original richness from disease data set
me3a$DentMax_log <- t((t(me3a$x) * dent_sd) + dent_mean)
me3a$DentMax <- exp(me3a$DentMax_log)



# plot of Max Pasteuria prevalence vs dentifera density at max prevalence
AUC_Dent <- ggplot(me3a, aes(x = DentMax, y = predicted)) +
  labs(tag = "C") +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, outline.type = NULL, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera, aes(x = dentifera.at.max+1, y = pasteuria.auc), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote(italic("D. dentifera") ~ "Density at Max Prevalence (no./" ~ "M"^2 ~")"), y = bquote(atop("Area Under the" ~ italic("Pasteuria"), "Prevalence Curve (prevalence x day)"))) +
  geom_text(aes(x = 1000, y= 4), label = "p = 0.02") +
  scale_x_log10(labels = scales::comma) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(AUC_Dent)
ggsave("AUCPastPrev_DentiferaAtMax_predict_color.png", plot = AUC_Dent, dpi = 300, width = 12, height = 10, units = "cm")





##### Models of AUC Pasteuria prev in dentifera, with mean densities * Year 
mod4 <- lmer(log(pasteuria.auc+1) ~ mean.pulicaria.density_z + mean.dentifera.density_z + mean.retrocurva.density_z + species.richness_z + Year +
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



# top model
mod4a <- lmer(log(pasteuria.auc+1) ~ mean.dentifera.density_z + (1|Lake), data = sum.data_dentifera)
summary(mod4a)




# second top model
mod4b <- lmer(log(pasteuria.auc+1) ~ mean.dentifera.density_z + mean.retrocurva.density_z + (1|Lake), data = sum.data_dentifera)
summary(mod4b)



# Marginal effects from top model
me4a <- ggpredict(mod4a, c("mean.dentifera.density_z [all]"))
plot(me4a, add.data = T)



# recalculate original mean Dentifera density 
me4a$x 
dent_mean <- mean(log(sum.data_dentifera$mean.dentifera.density+1)) # mean of original richness from disease data set
dent_sd <- sd(log(sum.data_dentifera$mean.dentifera.density+1)) # sd of original richness from disease data set
me4a$DentMean_log <- t((t(me4a$x) * dent_sd) + dent_mean)
me4a$DentMean <- exp(me4a$DentMean_log)



# plot of AUC Pasteuria prevalence vs mean dentifera density
AUC_Dent2 <- ggplot(me4a, aes(x = DentMean, y = predicted)) +
  labs(tag = "D") +
  geom_line(size = 1, color = "#1B9E77") +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2, outline.type = NULL, fill = "#1B9E77") +
  geom_point(data = sum.data_dentifera, aes(x = mean.dentifera.density+1, y = pasteuria.auc), size = 2, alpha = 0.4, color = "#1B9E77") +
  labs(x = bquote("Mean" ~italic("D. dentifera") ~ "Density (no./" ~ "M"^2 ~")"), y = bquote(atop("Area Under the" ~ italic("Pasteuria"), "Prevalence Curve (prevalence x day)"))) +
  geom_text(aes(x = 1000, y= 4), label = "p = 0.0066") +
  scale_x_log10(labels = scales::comma) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"))
print(AUC_Dent2)
ggsave("AUCPastPrev_MeanDentifera_predict_color.png", plot = AUC_Dent2, dpi = 300, width = 10, height = 10, units = "cm")



### four panel plot of MI field analyses
MI_full <- grid.arrange(MaxPrev_Pulic, MaxPrev_Dent, AUC_Dent, AUC_Dent2, nrow = 2)
ggsave("MI_Field_Analysis_predict_color.png", plot = MI_full, dpi = 300, width = 21, height = 20, units = "cm")


grid_arrange_shared_legend(MaxPrev_Pulic, MaxPrev_Dent, AUC_Dent, AUC_Dent2)











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









