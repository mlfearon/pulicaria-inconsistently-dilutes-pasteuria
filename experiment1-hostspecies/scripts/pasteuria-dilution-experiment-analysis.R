## This code analyzes Camden's 2015 Dilution Experiment

# Code written by Michelle Fearon

setwd("C:/Users/mlfea/OneDrive/Documents/PROJECTS/MHMP Daphnia Duffy/Generality of dilution (MHMP)/Puliaria dilution paper/Pulicaria dilution experiment 1 (Camden)")



# libraries
library(dplyr)
library(lme4)
library(ggeffects)
library(ggplot2)
library(car)
library(MuMIn)
library(emmeans)
library(RColorBrewer)


##### load experimental data

experiment <- read.csv("AllDilutionData.csv", stringsAsFactors = F, header = T)
head(experiment)
tail(experiment)

experiment <- experiment %>%
  mutate(Total.Tested = Total.Infected + Total.Uninfected) %>%
  arrange(Diluter.Species)

experiment$Diluter.Density_factor <- as.factor(experiment$Diluter.Density)



# model to test for dilution in pasteuria infected dentifera

mod <- glm(Prevalence ~ Diluter.Species * Diluter.Density, family = "binomial", weights = Total.Tested, data = experiment)
summary(mod)
vif(mod)
plot(mod)
Anova(mod, test.statistic = "Wald")



# plot of predicted values of prevalence by diluter density and diluter host species
me <- ggpredict(mod, c("Diluter.Density", "Diluter.Species"))
plot(me, add.data = T)


plot_predict <- ggplot(data = me, aes(x = x, y = predicted)) +
  geom_line(aes(color = group), size = 1) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, outline.type = NULL) +
  labs(x = "Density of Diluter Host (# per beaker)", y = bquote("Proportion of Susceptible" ~ italic("D. dentifera") ~ "Infected"), title = NULL, color = "Diluter Species") +
  scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  scale_fill_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  theme_classic()
print(plot_predict)
ggsave("DilutionExperiment_predicted.tiff", plot = plot_predict, dpi = 300, width = 12, height = 10, units = "cm", compression="lzw")



# attempt to add points to the plot above.
plot_predict2 <- plot_predict + 
  geom_jitter(data = experiment, aes(x = Diluter.Density, y = Prevalence, color = Diluter.Species, shape = Diluter.Species), size = 2, width = 0.4, height = 0.02, alpha = 0.7) +
  scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  scale_shape(name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva"))))
print(plot_predict2)
ggsave("DilutionExperiment_predicted2.tiff", plot = plot_predict2, dpi = 300, width = 12, height = 10, units = "cm", compression="lzw")



# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each host spp for a given host density = 4)
a <- emmeans(mod, specs = pairwise ~ Diluter.Species | Diluter.Density, type = "response")
a

# estimates of the slopes of each diluter species based on changes in density
ab <- emtrends(mod, pairwise ~ Diluter.Species, var="Diluter.Density", type = "response")
ab


# model to test for dilution in pasteria infected dentifera (with random effects)

mod2 <- glmer(Prevalence ~ Diluter.Species * Diluter.Density + (1|Replicate), family = "binomial", weights = Total.Tested, data = experiment)
summary(mod2)
vif(mod2)
plot(mod2)
Anova(mod2)

anova(mod, mod2)


me2 <- ggpredict(mod2, c("Diluter.Density", "Diluter.Species"))
plot(me2)






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


plot <- ggplot(data = exp, aes(x = Diluter.Density, y = prevalence, color = Diluter.Species, shape = Diluter.Species)) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = prevalence - se, ymax = prevalence + se), width = 0.3) +
  geom_line() +
  labs(x = "Density of Diluter Host (# per beaker)", y = bquote("Proportion of Susceptible" ~ italic("D. dentifera") ~ "Infected")) +
  scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  scale_shape_discrete(name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
  theme_classic()
print(plot)
ggsave("DilutionExperiment.png", plot = plot, dpi = 300, width = 12, height = 10, units = "cm")
