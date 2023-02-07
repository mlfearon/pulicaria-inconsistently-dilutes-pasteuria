## Experiment 2 Code for "Inconsistent dilution: Experimental but not field evidence for a dilution effect in Daphnia–bacteria interactions"


# Submitted to: Oecologia


# Code written by Michelle Fearon
# Last updated: Feb 6, 2023

## This code analyzes Michelle's 2021 Dilution Follow-up Experiment 
## Testing whether different D. pulicaria genotypes have a different capacity
## to dilute Metschnikowia and Pasteuria parasites


# libraries
library(tidyverse)
library(dplyr)
library(lme4)
library(glmmTMB)
library(logistf)
library(ggeffects)
library(ggplot2)
library(ggpubr)
library(car)
library(MuMIn)
library(emmeans)
library(olsrr)
library(RColorBrewer)
library(brms)
library(here)

# set the path to the script relative to the project root directory
here::i_am("experiment2-pulicariagenotypes/scripts/pulicaria-genotype-dilution-expt-analysis.R")


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

###################################
##### load experimental data

# Pasteuria and Metschnikowia prevalence data in Daphnia dentifera
dentifera_data <- read.csv(here("experiment2-pulicariagenotypes", "data", "DilutionDentiferaInfectionPrevalence.csv"), stringsAsFactors = F, header = T)
head(dentifera_data)
dim(dentifera_data)

# Pasteuria and Metschnikowia prevalence data in Daphnia pulicaria
pulicaria_data <- read.csv(here("experiment2-pulicariagenotypes", "data", "DilutionPulicariaInfectionPrevalence.csv"), stringsAsFactors = F, header = T)
head(pulicaria_data)

# check for any infections in the diluter pulicaria
metsch_prev_pulic <- pulicaria_data %>%
  filter(Parasite == "Metschnikowia") %>%
  summarize(Total.Metsch.Prev = sum(Total_Infected)/sum(Total_N))
metsch_prev_pulic  # no Metsch infections in pulicaria

past_prev_pulic <- pulicaria_data %>%
  filter(Parasite == "Pasteuria") %>%
  summarize(Total.Past.Prev = sum(Total_Infected)/sum(Total_N))
past_prev_pulic  # no Pasteuria infections in pulicaria




#############################
## Body Size data 
pulic_bodysize <- read.csv(here("experiment2-pulicariagenotypes", "data", "DilutionPulicariaBodySize.csv"), stringsAsFactors = F, header = T)
head(pulic_bodysize)
pulic_bodysize <- arrange(pulic_bodysize, PulicariaLine)



# calculate mean and sd of pulicaria body sizes
pulic_bodysize2 <- pulic_bodysize %>% 
  group_by(PulicariaLine, Parasite) %>%
  summarize(BodySize_um_mean = mean(BodySize_um), BodySize_um_sd = sd(BodySize_um), BodySize_mm_mean = mean(BodySize_mm), BodySize_mm_sd = sd(BodySize_mm))

View(pulic_bodysize2)

# visualize body size among different pulicaria genotypes
bodysizeplot <- ggplot(data = pulic_bodysize, aes(x = PulicariaLine, y = BodySize_um)) +
  geom_boxplot() +
  geom_jitter(alpha = 0.5, width = 0.1)
bodysizeplot

# check if there are significant differences in bodysize among pulicaria genotypes
size_mod <- lm(BodySize_um ~ PulicariaLine, data = pulic_bodysize)
summary(size_mod)
Anova(size_mod)
plot(size_mod) # residuals are fairly normally distributed
ols_plot_resid_qq(size_mod)
ols_test_normality(size_mod)  # Shapiro-Wilk test is not significant, indicates reasonably normal distribution
ols_test_breusch_pagan(size_mod) # test of homogenicity of variance is not significant, variance is constant


# check sig differences among each pair of genotypes
size_dif <- emmeans(size_mod, specs = pairwise ~ PulicariaLine, type = "response")
size_dif   # Yes, body size varies among different pulicaria genotypes

### Appendix S1: Figure S2
# visualize differences in bodysize among different pulicaria genotypes
me_size <- ggpredict(size_mod, c("PulicariaLine"))
plot(me_size, add.data = T) +
  labs(y= bquote("Body Size ("*mu*"m)"), x = "Pulicaria Genotype", title = NULL) +
  geom_text(aes(y = 2100), label = c("a","a","b","ab","b", "ab"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 11, color = "black"))
ggsave(here("experiment2-pulicariagenotypes", "figures", "Pulicaria_BodySize.tiff"), dpi = 600, width = 3, height = 3.5, units = "in", compression="lzw")


## pulicaria genotypes have significantly different bodysizes




#############################################
# add body size to dentifera and pulicaria prevalence data
dentifera_data2 <- full_join(dentifera_data, pulic_bodysize)
dentifera_data2 <- full_join(dentifera_data2, pulic_bodysize2)
View(dentifera_data2)


pulicaria_data2 <- full_join(pulicaria_data, pulic_bodysize)
pulicaria_data2 <- full_join(pulicaria_data2, pulic_bodysize2)
View(pulicaria_data2)


# split dentifera data by parasite
dentifera_metsch <- dentifera_data2 %>%
  filter(Parasite == "Metschnikowia")
dentifera_metsch <- arrange(dentifera_metsch, group_by = Treatment, PulicariaLine)
dentifera_metsch <- mutate(dentifera_metsch, Treatment2 = if_else(Treatment == "control", 1, 2))
str(dentifera_metsch)
dentifera_metsch$PulicariaLine <- factor(dentifera_metsch$PulicariaLine, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control1", "Control2"))
dentifera_metsch$PulicariaLine2 <- factor(dentifera_metsch$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


dentifera_metsch_lines <- filter(dentifera_metsch, PulicariaLine2 != "Control")
dentifera_metsch_lines$PulicariaLine2 <- factor(dentifera_metsch_lines$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W"))
dentifera_metsch_lines2 <- filter(dentifera_metsch_lines, PulicariaLine2 != "Pine") # remove genotype (Pine) that did not get infected
str(dentifera_metsch_lines)
dim(dentifera_metsch)
dim(dentifera_metsch_lines)
dim(dentifera_metsch_lines2)


dentifera_past <- dentifera_data2 %>%
  filter(Parasite == "Pasteuria")
dentifera_past <- arrange(dentifera_past, group_by = Treatment, PulicariaLine)
dentifera_past <- mutate(dentifera_past, Treatment2 = if_else(Treatment == "control", 1, 2))
str(dentifera_past)
dentifera_past$PulicariaLine <- factor(dentifera_past$PulicariaLine, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control1", "Control2"))
dentifera_past$PulicariaLine2 <- factor(dentifera_past$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


dentifera_past_lines <- filter(dentifera_past, PulicariaLine2 != "Control")
dentifera_past_lines$PulicariaLine2 <- factor(dentifera_past_lines$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W"))
dentifera_past_lines2 <- filter(dentifera_past_lines, PulicariaLine2 != "Pine", PulicariaLine2 != "Clear5") # remove genotype (Pine and Clear5) that did not get infected
dim(dentifera_past)
dim(dentifera_past_lines)
dim(dentifera_past_lines2)




#####################################
### Dentifera Metsch prevalence models

# first test comparing diluter treatments to controls
dent_metsch_mod <- glm(Prevalence ~ Treatment, family = "binomial", weights = N, data = dentifera_metsch)
summary(dent_metsch_mod)
Anova(dent_metsch_mod, test.statistic = "Wald")  # not significant interaction between genotype and bodysize on prevalence
overdisp_fun(dent_metsch_mod)

# pairwise comparison of control vs diluter treatments
metsch_treatment <- emmeans(dent_metsch_mod, specs = pairwise ~ Treatment, type = "response")
metsch_treatment

metsch_treatment_prob <- as.data.frame(metsch_treatment$emmeans)

## Figure 4C
# plot of predicted values of prevalence by diluter treatments vs controls
me_metsch_treatment <- ggpredict(dent_metsch_mod, c("Treatment"))
me_metsch_treatment$Treatment2 <- c("Control", "Diluter")
pulic_colors1 <- c("#666666", "#d95f02")
metsch_predict_treatment <- ggplot() + 
  scale_x_discrete(labels = c("Control", "Diluter")) +
  geom_jitter(data = dentifera_metsch, aes(x = Treatment, y = Prevalence, color = Treatment), size = 2, width = 0.2, height = 0.01, alpha = 0.6) +
  geom_pointrange(data = metsch_treatment_prob, aes(x = Treatment, y = prob, ymax = asymp.UCL, ymin = asymp.LCL), color = "black", linewidth = 1) +
  scale_color_manual(values = pulic_colors1) +
  scale_y_continuous(labels = scales::percent, limit = c(-0.05,1.02)) +
  labs(x = "Treatment", y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0.2), legend.position = "none")
metsch_predict_treatment
ggsave(here("experiment2-pulicariagenotypes", "figures", "Metsch_Diluter-v-Control.tiff"), plot = metsch_predict_treatment, dpi = 600, width = 3, height = 4, units = "in", compression="lzw")



#  Next test for effects of pulicaria genotype and body size, and their interaction (removed the controls and the Pine genotype that did not have any infection)
dent_metsch_mod2 <- glm(Prevalence ~ PulicariaLine2 * BodySize_mm, family = "binomial", weights = N, data = dentifera_metsch_lines2)
summary(dent_metsch_mod2)
vif(dent_metsch_mod2, type = "predictor")
overdisp_fun(dent_metsch_mod2)

### Appendix S1: Table S10
Anova(dent_metsch_mod2, test.statistic = "Wald")  # no significant interaction between genotype and bodysize on prevalence


# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in metsch prev between each pulicaria genotpye for a given body size = 1.89 mm)
a <- emmeans(dent_metsch_mod2, specs = pairwise ~ PulicariaLine2 | BodySize_mm, type = "response")
a


# plot of predicted values of prevalence by diluter density and diluter host species
me <- ggpredict(dent_metsch_mod2, c("PulicariaLine2", "BodySize_mm"))
plot(me, add.data = F)

#There is no significant interaction between bodysize and pulicaria genotype. But genotype and bodysize are also quite correlated together
# shown by the high VIFs for the model. Try testing genotype and bodysize effects on disease prevalence independently.

# plot_predict <- ggplot(data = me, aes(x = x, y = predicted)) +
#   geom_line(aes(color = x), size = 1) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, outline.type = NULL) +
#   labs(x = "Density of Diluter Host (# per beaker)", y = bquote("Proportion of Susceptible" ~ italic("D. dentifera") ~ "Infected"), title = NULL, color = "Diluter Species") +
#   scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
#   scale_fill_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
#   theme_classic()
# print(plot_predict)
# ggsave("DilutionExperiment_predicted.png", plot = plot_predict, dpi = 300, width = 12, height = 10, units = "cm")


# Since bodysize was not an important factor in the model, we removed it and focused on genotype effects only (model results presented in main text)
# Test of effect of pulicaria genotype only on Metsch prevalence (Pine genotype removed from analysis)
dent_metsch_mod3 <- glm(Prevalence ~ PulicariaLine2, family = "binomial", weights = N, data = dentifera_metsch_lines2)
summary(dent_metsch_mod3)
Anova(dent_metsch_mod3, test.statistic = "Wald")
overdisp_fun(dent_metsch_mod3)

predict(dent_metsch_mod3, type = "response")
dent_metsch_mod3$fitted.values
dent_metsch_mod3$ci.upper

# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotype)
b<- emmeans(dent_metsch_mod3, specs = pairwise ~ PulicariaLine2, type = "response")
b


#extract predicted infection probabilities and confidence intervals
metsch_prob <- as.data.frame(b$emmeans)

# plot of predicted values of prevalence by diluter genotype
me2 <- ggpredict(dent_metsch_mod3, c("PulicariaLine2"))
pulic_colors2 <- c("#d95f02")

## Figure 4D
metsch_genotype_predict <- ggplot() +
  geom_jitter(data = dentifera_metsch_lines, aes(x = PulicariaLine2, y = Prevalence), color = pulic_colors2, size = 2, width = 0.2, height = 0.01, alpha = 0.6) +
  geom_pointrange(data = metsch_prob, aes(x = PulicariaLine2, y = prob, ymax = asymp.UCL, ymin = asymp.LCL), color = "black", linewidth = 1) +
  scale_y_continuous(labels = scales::percent, limit = c(-0.05,1.02)) +
  labs(x = bquote(italic("D. pulicaria")~ "Diluter Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0.2))
metsch_genotype_predict

ggsave(here("experiment2-pulicariagenotypes", "figures", "Metsch_PulicariaGenotype_update.tiff"), plot = metsch_genotype_predict, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")



#### Other versions of the final Metschnikowia model presented in the text that we tried for dealing with separation in the data (i.e. treatments with zero infection prevalence) 
###  that ultimately did not produce successful models

# try with zero-inflated model
dent_metsch_mod3a <- glmmTMB(Prevalence ~ PulicariaLine2, ziformula = ~PulicariaLine2, family = "binomial", weights = N, data = dentifera_metsch_lines2, na.action=na.omit)
summary(dent_metsch_mod3a)
Anova(dent_metsch_mod3a, test.statistic = "Chisq")
#plot(dent_metsch_mod2)
#overdisp_fun(dent_metsch_mod2a)

# Zero-inflated model is 10 AIC higher than the mod3 above.
AIC(dent_metsch_mod3, dent_metsch_mod3a)


# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotype)
c <- emmeans(dent_metsch_mod3a, specs = pairwise ~ PulicariaLine2, type = "response")
c


# plot of predicted values of prevalence by diluter genotype
me3a <- ggpredict(dent_metsch_mod3a, c("PulicariaLine2"), type = "zi.prob")
me3a$x <- factor(me3a$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))
View(me3a)

pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_a <- plot(me3a, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
metsch_predict_a
ggsave(here("experiment2-pulicariagenotypes", "figures", "Metsch_PulicariaGenotype_zeroinf.tiff"), dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# try the model above with a logistic regression
dent_metsch_mod3b <- logistf(Prevalence ~ PulicariaLine2, weights = N, data = dentifera_metsch)
summary(dent_metsch_mod3b)
anova(dent_metsch_mod3b)

predict(dent_metsch_mod3b, type = "response")
coef(dent_metsch_mod3b, type = "response")
exp(confint(dent_metsch_mod3b, type = "response"))
dent_metsch_mod3b$ci.lower
dent_metsch_mod3b$ci.upper

me3b <- ggpredict(dent_metsch_mod3b, c("PulicariaLine2"))
me3b$x <- factor(me3b$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_b <- plot(me3b, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", linewidth = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
metsch_predict_b
ggsave(here("experiment2-pulicariagenotypes", "figures", "Metsch_PulicariaGenotype_logistic.tiff"), dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# try the model above with a rare events logistic regression
library(Zelig)
dent_metsch_mod3c <- relogit(Prevalence ~ PulicariaLine2, weights = N, data = dentifera_metsch)
summary(dent_metsch_mod3c)
anova(dent_metsch_mod3c)


me3c <- ggpredict(dent_metsch_mod3c, c("PulicariaLine2"))
me3c$x <- factor(me3c$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_c <- plot(me3c, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", linewidth = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
metsch_predict_c
ggsave(here("experiment2-pulicariagenotypes", "figures", "Metsch_PulicariaGenotype_logistic2.tiff"), dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


library(logistf)
# try the model above with a logistic regression with added covariate (to help model predictions!)
dent_metsch_mod3d <- flac(lfobject = dent_metsch_mod2b, data = dentifera_metsch)
summary(dent_metsch_mod3d)

predict(dent_metsch_mod3d, type = "response")


# try the model above with a logistic regression with adjusted intercept (to help model predictions!)
dent_metsch_mod3e <- flic(lfobject = dent_metsch_mod2b)
dent_metsch_mod3e <- logistf(Prevalence ~ PulicariaLine2, weights = N, flic = T, data = dentifera_metsch)
summary(dent_metsch_mod3e)


predict(dent_metsch_mod3e, type = "response")
coef(dent_metsch_mod3e, type = "response")
exp(confint(dent_metsch_mod3e, type = "response"))
dent_metsch_mod3e$ci.lower
dent_metsch_mod3e$ci.upper
get_sigma(dent_metsch_mod3e)

me3e <- ggpredict(dent_metsch_mod3e, c("PulicariaLine2"))
me3e$x <- factor(me3e$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_e <- plot(me3e, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", linewidth = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
metsch_predict_e
ggsave(here("experiment2-pulicariagenotypes", "figures", "Metsch_PulicariaGenotype_logistic3.tiff"), dpi = 600, width = 5, height = 4, units = "in", compression="lzw")

# ultimately none of these other attempts produced a better model than the one presented in the main text




### There is a significant difference in Metsch prevalence between pulicaria diluters and control treatments, but there is not a significant 
### difference in Metsch prevalence among the pulicaria lines, nor with bodysize

#  comparison of pulicaria bodysize against Metsch prevalence in dentifera
dent_metsch_mod4 <- glm(Prevalence ~ BodySize_mm_mean, family = "binomial", weights = N, data = dentifera_metsch_lines2)
summary(dent_metsch_mod4)
Anova(dent_metsch_mod4, test.statistic = "Wald")
overdisp_fun(dent_metsch_mod4)


# plot of predicted values of prevalence by pulicaria bodysize
me4 <- ggpredict(dent_metsch_mod4, c("BodySize_mm_mean [all]"))
plot(me4, add.data = F) + 
  labs(x = "Body Size (mm)", y = "Fungal Prevalence (%)", title = NULL)

## No significant effect of body size on Metschnikowia prevalence






#############################################
### Dentifera Pasteuria prevalence models

# first test comparing diluter treatments to controls
# Comparison of diluter vs control treatments for pasteuria
#  comparison of pulicaria genotype and controls only
dent_past_mod <- glm(Prevalence ~ Treatment, family = "binomial", weights = N, data = dentifera_past)
summary(dent_past_mod)
Anova(dent_past_mod, test.statistic = "Wald")
overdisp_fun(dent_past_mod)


past_treatment <- emmeans(dent_past_mod, specs = pairwise ~ Treatment, type = "response")
past_treatment

past_treatment_prob <- as.data.frame(past_treatment$emmeans)


# Figure 4A
# plot of predicted values of prevalence by diluter treatments vs controls
me_past_treatment <- ggpredict(dent_past_mod, c("Treatment"))
me_past_treatment$Treatment2 <- c("Control", "Diluter")
past_predict_treatment <- ggplot() +
  scale_x_discrete(labels = c("Control", "Diluter")) +
  geom_jitter(data = dentifera_past, aes(x = Treatment, y = Prevalence, color = Treatment), size = 2, width = 0.2, height = 0.01, alpha = 0.6) +
  geom_pointrange(data = past_treatment_prob, aes(x = Treatment, y = prob, ymax = asymp.UCL, ymin = asymp.LCL), color = "black", linewidth = 1) +
  scale_color_manual(values = pulic_colors1) +
  scale_y_continuous(labels = scales::percent, limit = c(-0.05,1.02)) +
  labs(x = "Treatment", y = bquote(italic("Pasteuria ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"), legend.position = "none")
past_predict_treatment
ggsave(here("experiment2-pulicariagenotypes", "figures", "Past_Diluter-v-Control.tiff"), plot = past_predict_treatment, dpi = 600, width = 3, height = 4, units = "in", compression="lzw")




#  Next test for effects of pulicaria genotype and body size, and their interaction (removed the controls and the Pine genotype that did not have any infection)
dent_past_mod2 <- glm(Prevalence ~ PulicariaLine2 * BodySize_mm, family = "binomial", weights = N, data = dentifera_past_lines2)
summary(dent_past_mod2)  #something isn't quite right about this model, weird predictions
vif(dent_past_mod2)
overdisp_fun(dent_past_mod2)

## Appendix S1: Table S10
Anova(dent_past_mod2, test.statistic = "Wald")

# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotype for a given body size = 1.89 mm)
d <- emmeans(dent_past_mod2, specs = pairwise ~ PulicariaLine2 | BodySize_mm, type = "response")
d


# plot of predicted values of prevalence by diluter density and diluter host species
me4 <- ggpredict(dent_past_mod2, c("PulicariaLine2", "BodySize_mm"))
plot(me4, add.data = F)


# plot_predict <- ggplot(data = me, aes(x = x, y = predicted)) +
#   geom_line(aes(color = x), size = 1) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, outline.type = NULL) +
#   labs(x = "Density of Diluter Host (# per beaker)", y = bquote("Proportion of Susceptible" ~ italic("D. dentifera") ~ "Infected"), title = NULL, color = "Diluter Species") +
#   scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
#   scale_fill_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
#   theme_classic()
# print(plot_predict)
# ggsave("DilutionExperiment_predicted.png", plot = plot_predict, dpi = 300, width = 12, height = 10, units = "cm")



# Since bodysize was not an important factor in the model, we removed it and focused on genotype effects only (model results presented in main text)
# Test of effect of pulicaria genotype only on Pasteuria prevalence (Pine and Clear5 genotypes removed from analysis due to no infection in those genotypes)
#  comparison of pulicaria genotype and controls only
dent_past_mod3 <- glm(Prevalence ~ PulicariaLine2, family = "binomial", weights = N, data = dentifera_past_lines2)
summary(dent_past_mod3)
Anova(dent_past_mod3, test.statistic = "Wald")
overdisp_fun(dent_past_mod3)




# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotype)
e<- emmeans(dent_past_mod3, specs = pairwise ~ PulicariaLine2, type = "response")
e


#extract predicted infection probabilities and confidence intervals
past_prob <- as.data.frame(e$emmeans)

# plot of predicted values of prevalence by diluter genotype
me5 <- ggpredict(dent_past_mod3, c("PulicariaLine2"))
pulic_colors2 <- c("#d95f02")

## Figure 4D
past_genotype_predict <- ggplot() +
  geom_jitter(data = dentifera_past_lines, aes(x = PulicariaLine2, y = Prevalence), color = pulic_colors2, size = 2, width = 0.2, height = 0.01, alpha = 0.6) +
  geom_pointrange(data = past_prob, aes(x = PulicariaLine2, y = prob, ymax = asymp.UCL, ymin = asymp.LCL), color = "black", linewidth = 1) +
  scale_y_continuous(labels = scales::percent, limit = c(-0.05,1.02)) +
  labs(x = bquote(italic("D. pulicaria")~ "Diluter Genotype"), y = bquote(italic("Pasteuria ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
past_genotype_predict

ggsave(here("experiment2-pulicariagenotypes", "figures", "Past_PulicariaGenotype_update.tiff"), plot = past_genotype_predict, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")



## Figure 4, panels A-D
# joint figure of pasteuria and metsch prevalence control vs diluters and prevalence vs pulicaria genotype
Figure4 <- ggarrange(past_predict_treatment, past_genotype_predict, metsch_predict_treatment, metsch_genotype_predict, labels = c("A", "B", "C", "D"), ncol = 2, nrow = 2, widths = c(3,5))
ggsave(here("experiment2-pulicariagenotypes", "figures", "Figure4.tiff"), plot = Figure3, dpi = 600, width = 7, height = 6, units = "in", compression="lzw")




### There is a significant difference in Pasteuria prevalence between pulicaria diluters and control treatments, but there is not a significant 
### difference in Pasteuria prevalence among the pulicaria lines, nor with bodysize
# prevalence vs body size
dent_past_mod4 <- glm(Prevalence ~ BodySize_mm_mean, family = "binomial", weights = N, data = dentifera_past_lines2)
summary(dent_past_mod4)  #not significant
Anova(dent_past_mod4, test.statistic = "Wald")
overdisp_fun((dent_past_mod4))

# Pasteuria prevalence in dentifera is not affected by pulicaria body size



#### Other versions of the final Pasteuria model presented in the text that we tried for dealing with separation in the data (i.e. treatments with zero infection prevalence) 
###  that ultimately did not produce successful models

# try the model above with a logistic regression
dent_past_mod2b <- logistf(Prevalence ~ PulicariaLine2, weights = N, data = dentifera_past)
summary(dent_past_mod2b)


predict(dent_past_mod2b, type = "response")
coef(dent_past_mod2b, type = "response")
exp(confint(dent_past_mod2b, type = "response"))
dent_past_mod2b$ci.lower
dent_past_mod2b$ci.upper

me2b <- ggpredict(dent_past_mod2b, c("PulicariaLine2"))
me2b$x <- factor(me2b$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

past_predict_b <- plot(me2b, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Pasteuria ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", linewidth = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
past_predict_b
ggsave(here("experiment2-pulicariagenotypes", "figures", "past_PulicariaGenotype_logistic.tiff"), plot = past_predict_b, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")





# try the model above with a rare events logistic regression
library(Zelig)
dent_past_mod2c <- relogit(Prevalence ~ PulicariaLine2, weights = N, data = dentifera_past)
summary(dent_past_mod2c)
anova(dent_past_mod2c)



me2c <- ggpredict(dent_past_mod2c, c("PulicariaLine2"))
me2c$x <- factor(me2c$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

past_predict_c <- plot(me2c, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Pasteuria ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", linewidth = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
past_predict_c
ggsave(here("experiment2-pulicariagenotypes", "figures", "past_PulicariaGenotype_logistic2.tiff"), plot = past_predict_c, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# try the model above with a logistic regression with added covariate (to help model predictions!)
dent_past_mod2d <- flac(lfobject = dent_past_mod2b, data = dentifera_past)
summary(dent_past_mod2d)

predict(dent_past_mod2d, type = "response")


# try the model above with a logistic regression with adjusted intercept (to help model predictions!)
dent_past_mod2e <- flic(lfobject = dent_past_mod2b)
dent_past_mod2e <- logistf(Prevalence ~ PulicariaLine2, weights = N, flic = T, data = dentifera_past)
summary(dent_past_mod2e)


predict(dent_past_mod2e, type = "response", se.fit = TRUE)
coef(dent_past_mod2e, type = "response")
exp(confint(dent_past_mod2e, type = "response"))
dent_past_mod2e$ci.lower
dent_past_mod2e$ci.upper
get_sigma(dent_past_mod2e)

me2e <- ggpredict(dent_past_mod2e, c("PulicariaLine2"))
me2e$x <- factor(me2e$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

past_predict_e <- plot(me2e, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Pasteuria ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", linewidth = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black"))
past_predict_e
ggsave(here("experiment2-pulicariagenotypes", "figures", "past_PulicariaGenotype_logistic3.tiff"), plot = past_predict_e, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")
