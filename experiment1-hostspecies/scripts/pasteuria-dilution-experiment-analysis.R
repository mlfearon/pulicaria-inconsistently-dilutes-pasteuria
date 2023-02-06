## Experiment 1 Code for "Mixed evidence for a dilution effect in Daphnia communities infected by a bacterial parasite"


# Submitted to: Oecologia


# Code written by Michelle Fearon
# Last updated: Jan 12, 2023

## This code analyzes Camden's 2015 Pasteuria Dilution Experiment 
## Testing whether resistant Daphnia dentifera, D. retrocurva, or D. pulicaria host species 
## can dilute pasteuria using similar experimental conditions to previous tests showing 
## D. pulicaria diluting Metschnikowia (Hall et al 2009)


# libraries
library(dplyr)
library(lme4)
library(ggeffects)
library(ggplot2)
library(car)
library(emmeans)
library(RColorBrewer)
library(here)


# set the path to the script relative to the project root directory
here::i_am("experiment1-hostspecies/scripts/pasteuria-dilution-experiment-analysis.R")


## updated overdispersion function from Ben Bolker
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}




##### load experimental data

experiment <- read.csv(here("experiment1-hostspecies", "data", "pasteuria-dilution-data.csv"), stringsAsFactors = F, header = T)
head(experiment)
tail(experiment)
dim(experiment)

experiment <- experiment %>%
  mutate(Total.Tested = Total.Infected + Total.Uninfected) %>%
  arrange(Diluter.Species) %>%
  mutate(ID = 1:119)

experiment$Diluter.Density_factor <- as.factor(experiment$Diluter.Density)

table(experiment$Diluter.Species)
table(experiment$Diluter.Density)

# model to test for dilution in pasteuria infected dentifera

mod <- glm(Prevalence ~ Diluter.Species * Diluter.Density, family = "binomial", weights = Total.Tested, data = experiment)
summary(mod)
vif(mod)
Anova(mod, test.statistic = "Wald")
overdisp_fun(mod)   # model is overdispersed
plot(mod)


# updated model to test for dilution in pasteria infected dentifera (with ID random effect to control for overdispersion)
mod2 <- glmer(Prevalence ~ Diluter.Species * Diluter.Density + (1|ID), family = "binomial", weights = Total.Tested, data = experiment)
summary(mod2)
vif(mod2)
plot(mod2)
Anova(mod2)
overdisp_fun(mod2)  # overdispersion is controlled for now.



# plot of predicted values of prevalence by diluter density and diluter host species
me <- ggpredict(mod2, c("Diluter.Density", "Diluter.Species"))
plot(me, add.data = T)


plot_predict <- ggplot(data = me, aes(x = x, y = predicted)) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2) +
  labs(x = "Density of Diluter Host (# per beaker)", y = bquote("Proportion of Susceptible" ~ italic("D. dentifera") ~ "Infected"), title = NULL, color = "Diluter Species") +
  scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  scale_fill_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  theme_classic()
print(plot_predict)
ggsave(here("experiment1-hostspecies", "figures", "DilutionExperiment_predicted.tiff"), plot = plot_predict, dpi = 300, width = 12, height = 10, units = "cm", compression="lzw")



# add points to the plot above.
plot_predict2 <- plot_predict + 
  geom_jitter(data = experiment, aes(x = Diluter.Density, y = Prevalence, color = Diluter.Species, shape = Diluter.Species), size = 2, width = 0.4, height = 0.02, alpha = 0.7) +
  scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  scale_shape(name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 10, color = "black"), axis.title.y = element_text(size = 9, color = "black"), 
        legend.position = "top", legend.title = element_text(size = 8, color = "black"), legend.text = element_text(size = 7, color = "black"), legend.justification = "left", legend.margin = margin(0,0,0,-10), legend.box.margin = margin(0,0,0,-20))
print(plot_predict2)
ggsave(here("experiment1-hostspecies", "figures", "Figure1.tiff"), plot = plot_predict2, dpi = 300, width = 3.5, height = 4, units = "in", compression="lzw")



# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each host spp for a given host density = 4)
a <- emmeans(mod2, specs = pairwise ~ Diluter.Species | Diluter.Density, type = "response")
a

# estimates of the slopes of each diluter species based on changes in density
ab <- emtrends(mod2, pairwise ~ Diluter.Species, var="Diluter.Density", type = "response", infer = c(TRUE, TRUE), at=list(Diluter.Density=c(4)))
ab


mean_prev <- experiment %>%
  group_by(Diluter.Species, Diluter.Density_factor) %>%
  summarize(Prev_mean = mean(Prevalence))


## CUT FOR FINAL VERSION OF CODE??


# model to test for dilution in pasteria infected dentifera (density as a factor rather than continuous)

mod3 <- glm(Prevalence ~ Diluter.Species * Diluter.Density_factor, family = "binomial", weights = Total.Tested, data = experiment)
summary(mod3)
vif(mod3)
plot(mod3)
Anova(mod3)



me3 <- ggpredict(mod3, c("Diluter.Density_factor", "Diluter.Species"))
plot(me3)


# calculates the pairwise tests for each genus within each density category 
b <- emmeans(mod3, specs = pairwise ~ Diluter.Species | Diluter.Density_factor, type = "response")
b

# pulicaria is consistently different from dentifera at all densities except for 0
# pulicaria is different from retrocurva at density of 2, 6, and 8
# retrocurva is different from dentifera at density of 4 and 6





# figure of prevalence across diluter density by host species 
# based on raw data

exp <- experiment %>% 
  group_by(Diluter.Species, Diluter.Density) %>%
  summarize(prevalence = mean(Prevalence), sd = sd(Prevalence), se = sd(Prevalence)/sqrt(length(Prevalence))) %>%
  arrange(Diluter.Species)

# similar to Camden's original plot
plot <- ggplot(data = exp, aes(x = Diluter.Density, y = prevalence, color = Diluter.Species, shape = Diluter.Species)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), width = 0.3) +
  geom_line() +
  labs(x = "Density of Diluter Host (# per beaker)", y = bquote("Proportion of Susceptible" ~ italic("D. dentifera") ~ "Infected")) +
  scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  scale_shape_discrete(name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  theme_classic()
print(plot)
ggsave(here("experiment1-hostspecies", "figures","DilutionExperiment.png"), plot = plot, dpi = 300, width = 12, height = 10, units = "cm")
