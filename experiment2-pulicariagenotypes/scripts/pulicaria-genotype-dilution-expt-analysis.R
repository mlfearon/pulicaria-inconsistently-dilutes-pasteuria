## This code analyzes Michelle's 2021 Dilution Follow-up Experiment 

# Code written by Michelle Fearon


# libraries
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
library(RColorBrewer)
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



##### load experimental data


# Pasteuria and Metschnikowia prevalence data in Daphnia dentifera
dentifera_data <- read.csv(here("experiment2-pulicariagenotypes", "data", "DilutionDentiferaInfectionPrevalence.csv"), stringsAsFactors = F, header = T)
head(dentifera_data)


# Pasteuria and Metschnikowia prevalence data in Daphnia pulicaria
pulicaria_data <- read.csv(here("experiment2-pulicariagenotypes", "data", "DilutionPulicariaInfectionPrevalence.csv"), stringsAsFactors = F, header = T)
head(pulicaria_data)


## Body Size data 
pulic_bodysize <- read.csv(here("experiment2-pulicariagenotypes", "data", "DilutionPulicariaBodySize.csv"), stringsAsFactors = F, header = T)
head(pulic_bodysize)
pulic_bodysize <- arrange(pulic_bodysize, PulicariaLine)

?order()

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
plot(size_mod) # residuals are fairly normally distributed
Anova(size_mod)

# check sig differences among each pair of genotypes
size_dif <- emmeans(size_mod, specs = pairwise ~ PulicariaLine, type = "response")
size_dif

# visualize differences in bodysize among different pulicaria genotypes
me_size <- ggpredict(size_mod, c("PulicariaLine"))
plot(me_size, add.data = T) +
  labs(y= bquote("Body Size ("~ mu~"m)"), x = "Pulicaria Genotype", title = NULL) +
  geom_text(aes(y = 2100), label = c("a","a","b","ab","b", "ab"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title = element_text(size = 11, color = "black"))
ggsave("Pulicaria_BodySize.tiff", dpi = 600, width = 3, height = 3.5, units = "in", compression="lzw")





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
str(dentifera_metsch)
dentifera_metsch$PulicariaLine <- factor(dentifera_metsch$PulicariaLine, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control1", "Control2"))
dentifera_metsch$PulicariaLine2 <- factor(dentifera_metsch$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


dentifera_metsch_lines <- filter(dentifera_metsch, PulicariaLine2 != "Control")
dentifera_metsch_lines$PulicariaLine2 <- factor(dentifera_metsch_lines$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W"))
str(dentifera_metsch_lines)



dentifera_past <- dentifera_data2 %>%
  filter(Parasite == "Pasteuria")
dentifera_past <- arrange(dentifera_past, group_by = Treatment, PulicariaLine)
str(dentifera_past)
dentifera_past$PulicariaLine <- factor(dentifera_past$PulicariaLine, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control1", "Control2"))
dentifera_past$PulicariaLine2 <- factor(dentifera_past$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


dentifera_past_lines <- filter(dentifera_past, PulicariaLine2 != "Control")
dentifera_past_lines$PulicariaLine2 <- factor(dentifera_past_lines$PulicariaLine2, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W"))







### Dentifera Metsch prevalence models


#  comparison of pulicaria genotype and body size
dent_metsch_mod <- glm(Prevalence ~ PulicariaLine * BodySize_mm, family = "binomial", weights = N, data = dentifera_metsch)
summary(dent_metsch_mod)
vif(dent_metsch_mod)
plot(dent_metsch_mod)
Anova(dent_metsch_mod, test.statistic = "Wald")  # not significant interaction between genotype and bodysize on prevalence
overdisp_fun(dent_metsch_mod)

# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotpye for a given body size = 1.89 mm)
a <- emmeans(dent_metsch_mod, specs = pairwise ~ PulicariaLine | BodySize_mm, type = "response")
a


# plot of predicted values of prevalence by diluter density and diluter host species
me <- ggpredict(dent_metsch_mod, c("PulicariaLine", "BodySize_mm"))
plot(me, add.data = F)


# plot_predict <- ggplot(data = me, aes(x = x, y = predicted)) +
#   geom_line(aes(color = x), size = 1) +
#   geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, outline.type = NULL) +
#   labs(x = "Density of Diluter Host (# per beaker)", y = bquote("Proportion of Susceptible" ~ italic("D. dentifera") ~ "Infected"), title = NULL, color = "Diluter Species") +
#   scale_color_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
#   scale_fill_brewer(palette = "Dark2", name = "Diluter Species", labels = c(bquote(italic("D. dentifera")), bquote(italic("D. pulicaria")), bquote(italic("D. retrocurva")))) +
#   theme_classic()
# print(plot_predict)
# ggsave("DilutionExperiment_predicted.png", plot = plot_predict, dpi = 300, width = 12, height = 10, units = "cm")



#  comparison of pulicaria genotype and controls only
dent_metsch_mod2 <- glm(Prevalence ~ PulicariaLine2, family = "binomial", weights = N, data = dentifera_metsch)
summary(dent_metsch_mod2)
Anova(dent_metsch_mod2, test.statistic = "Wald")
#plot(dent_metsch_mod2)
overdisp_fun(dent_metsch_mod2)

predict(dent_metsch_mod2, type = "response")
dent_metsch_mod2$fitted.values
dent_metsch_mod2$ci.upper

# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotype)
b<- emmeans(dent_metsch_mod2, specs = pairwise ~ PulicariaLine2, type = "response")
b


# plot of predicted values of prevalence by diluter genotype
me2 <- ggpredict(dent_metsch_mod2, c("PulicariaLine2"), type = "zi.prob")
me2$x <- factor(me2$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))
pulic_colors2 <- c(rep("#d95f02", 6))

metsch_predict <- plot(me2, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
metsch_predict
ggsave("Metsch_PulicariaGenotype.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")




# try with zero-inflated model
dent_metsch_mod2a <- glmmTMB(Prevalence ~ PulicariaLine2, ziformula = ~PulicariaLine2, family = "binomial", weights = N, data = dentifera_metsch)
summary(dent_metsch_mod2a)
Anova(dent_metsch_mod2a, test.statistic = "Wald")
#plot(dent_metsch_mod2)
overdisp_fun(dent_metsch_mod2a)



# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotype)
b<- emmeans(dent_metsch_mod2a, specs = pairwise ~ PulicariaLine2, type = "response")
b


# plot of predicted values of prevalence by diluter genotype
me2a <- ggpredict(dent_metsch_mod2a, c("PulicariaLine2"), type = "zi.prob")
me2a$x <- factor(me2a$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_a <- plot(me2a, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
metsch_predict_a
ggsave("Metsch_PulicariaGenotype_zeroinf.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")






# try the model above with a logistic regression
dent_metsch_mod2b <- logistf(Prevalence ~ PulicariaLine2, weights = N, data = dentifera_metsch)
summary(dent_metsch_mod2b)
anova(dent_metsch_mod2b)

predict(dent_metsch_mod2b, type = "response")
coef(dent_metsch_mod2b, type = "response")
exp(confint(dent_metsch_mod2b, type = "response"))
dent_metsch_mod2b$ci.lower
dent_metsch_mod2b$ci.upper

me2b <- ggpredict(dent_metsch_mod2b, c("PulicariaLine2"))
me2b$x <- factor(me2b$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_b <- plot(me2b, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
metsch_predict_b
ggsave("Metsch_PulicariaGenotype_logistic.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")





# try the model above with a rare events logistic regression
library(Zelig)
dent_metsch_mod2c <- relogit(Prevalence ~ PulicariaLine2, weights = N, data = dentifera_metsch)
summary(dent_metsch_mod2c)
anova(dent_metsch_mod2c)



me2c <- ggpredict(dent_metsch_mod2c, c("PulicariaLine2"))
me2c$x <- factor(me2c$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_c <- plot(me2c, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
metsch_predict_c
ggsave("Metsch_PulicariaGenotype_logistic2.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


library(logistf)
# try the model above with a logistic regression with added covariate (to help model predictions!)
dent_metsch_mod2d <- flac(lfobject = dent_metsch_mod2b)
summary(dent_metsch_mod2d)

predict(dent_metsch_mod2d, type = "response")


# try the model above with a logistic regression with adjusted intercept (to help model predictions!)
dent_metsch_mod2e <- flic(lfobject = dent_metsch_mod2b)
dent_metsch_mod2e <- logistf(Prevalence ~ PulicariaLine2, weights = N, flic = T, data = dentifera_metsch)
summary(dent_metsch_mod2e)


predict(dent_metsch_mod2e, type = "response")
coef(dent_metsch_mod2e, type = "response")
exp(confint(dent_metsch_mod2e, type = "response"))
dent_metsch_mod2e$ci.lower
dent_metsch_mod2e$ci.upper
get_sigma(dent_metsch_mod2e)

me2e <- ggpredict(dent_metsch_mod2e, c("PulicariaLine2"))
me2e$x <- factor(me2e$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

metsch_predict_e <- plot(me2e, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
metsch_predict_e
ggsave("Metsch_PulicariaGenotype_logistic3.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")









### There is a significant difference between pulicaria diluters and control treatments, but there is not a significant 
### difference among the pulicaria lines, nor with boysize


#  comparison of pulicaria bodysize against Metsch prevalence in dentifera
dent_metsch_mod3 <- glm(Prevalence ~ BodySize_mm_mean, family = "binomial", weights = N, data = dentifera_metsch)
summary(dent_metsch_mod3)
Anova(dent_metsch_mod3, test.statistic = "Wald")
plot(dent_metsch_mod3)


# plot of predicted values of prevalence by pulicaria bodysize
me3 <- ggpredict(dent_metsch_mod3, c("BodySize_mm_mean [all]"))
plot(me3, add.data = F) + 
  labs(x = "Body Size (mm)", y = "Fungal Prevalence (%)", title = NULL)

## No significant effect of body size


# Comparison of diluter vs control treatments for pasteuria
#  comparison of pulicaria genotype and controls only
dent_metsch_mod4 <- glm(Prevalence ~ Treatment, family = "binomial", weights = N, data = dentifera_metsch)
summary(dent_metsch_mod4)
Anova(dent_metsch_mod4, test.statistic = "Wald")
plot(dent_metsch_mod4)



# plot of predicted values of prevalence by diluter treatments vs controls
me4 <- ggpredict(dent_metsch_mod4, c("Treatment"))
plot(me4, add.data = F) + 
  labs(x = "Treatment", y = "Fungal Prevalence (%)", title = NULL)





### Dentifera Pasteuria prevalence models


#  comparison of pulicaria genotype and body size
dent_past_mod <- glm(Prevalence ~ PulicariaLine * BodySize_mm, family = "binomial", weights = N, data = dentifera_past)
summary(dent_past_mod)  #something isn't quite right about this model, weird predictions
vif(dent_past_mod)
Anova(dent_past_mod, test.statistic = "Wald")
plot(dent_past_mod)

# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotpye for a given body size = 1.89 mm)
c <- emmeans(dent_past_mod, specs = pairwise ~ PulicariaLine | BodySize_mm, type = "response")
c


# plot of predicted values of prevalence by diluter density and diluter host species
me4 <- ggpredict(dent_past_mod, c("PulicariaLine", "BodySize_mm"))
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



#  comparison of pulicaria genotype and controls only
dent_past_mod2 <- glm(Prevalence ~ PulicariaLine2, family = "binomial", weights = N, data = dentifera_past)
summary(dent_past_mod2)
Anova(dent_past_mod2, test.statistic = "Wald")
plot(dent_past_mod2)
overdisp_fun(dent_past_mod2)


# calculates the pairwise tests for each genus within each site 
#(determines if there are sig differences in Pasteuria prev between each pulicaria genotype)
d<- emmeans(dent_past_mod2, specs = pairwise ~ PulicariaLine2, type = "response")
d


# plot of predicted values of prevalence by diluter density and diluter host species
me5 <- ggpredict(dent_past_mod2, c("PulicariaLine2"))
View(me5)
str(me5)
me5$x <- factor(me5$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))




pulic_colors <- c(rep("#d95f02", 6), rep("black",1))


past_predict <- plot(me5, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Pasteuria ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("ab","ab","ab","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 9, color = "black"))
past_predict
ggsave("Pasteuria_PulicariaGenotype.tiff", plot = past_predict, dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


## Figure 3
# joint figure of pasteuria and metsch prevalence vs pulicaria genotype
Figure <- ggarrange(past_predict, metsch_predict, labels = c("A", "B"),ncol = 1, nrow = 2)
ggsave("Prevalence_PulicariaGenotype.tiff", plot = Figure, dpi = 600, width = 4, height = 6, units = "in", compression="lzw")


# past_predict <- ggplot(data = dentifera_past, aes(x = PulicariaLine, y = Prevalence)) +
#   geom_jitter(width = 0.2, height = 0, alpha = 0.5, color = "#d95f02") +
#   geom_vline(xintercept = 6.5, color = "black") +
#   labs(x = "Pulicaria Genotype", y = "Bacterial Prevalence (%)") +
#   theme_classic()
# past_predict
# 
# past_predict_update <- past_predict +
#   geom_point(data = me5, aes(x= x, y = predicted), size = 0.6)
# #geom_errorbar(data = me5, aes(x= x,  ymin = conf.low, ymax = conf.high))
# past_predict_update
# ggsave("Pasteuria_PulicariaGenotype.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")
# 



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
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("pastnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
past_predict_b
ggsave("past_PulicariaGenotype_logistic.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")





# try the model above with a rare events logistic regression
library(Zelig)
dent_past_mod2c <- relogit(Prevalence ~ PulicariaLine2, weights = N, data = dentifera_past)
summary(dent_past_mod2c)
anova(dent_past_mod2c)



me2c <- ggpredict(dent_past_mod2c, c("PulicariaLine2"))
me2c$x <- factor(me2c$x, levels = c("BA", "Clear5", "Clover", "Mid67", "Pine", "W", "Control"))


pulic_colors <- c(rep("#d95f02", 6), rep("gray",1))

past_predict_c <- plot(me2c, add.data = T, jitter = c(0.5,0.01), colors = pulic_colors) + 
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
past_predict_c
ggsave("past_PulicariaGenotype_logistic2.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")


# try the model above with a logistic regression with added covariate (to help model predictions!)
dent_past_mod2d <- flac(lfobject = dent_past_mod2b)
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
  labs(x = bquote(italic("D. pulicaria")~ "Genotype"), y = bquote(italic("Metschnikowia ")~"Prevalence in" ~ italic("D. dentifera")), title = NULL) + 
  geom_vline(xintercept = 6.5, color = "black", size = 1.5) +
  geom_text(aes(y = 1.05), label = c("b","b","b","b","ab", "b", "a"), position = position_dodge(width = 0.4), show.legend = F, size = 10/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 8, color = "black"), axis.title.x = element_text(size = 11, color = "black"), axis.title.y = element_text(size = 8.5, color = "black", hjust = 0))
past_predict_e
ggsave("past_PulicariaGenotype_logistic3.tiff", dpi = 600, width = 5, height = 4, units = "in", compression="lzw")





?vegan::adonis




# prevalence vs body size
dent_past_mod3 <- glm(Prevalence ~ BodySize_mm_mean, family = "binomial", weights = N, data = dentifera_past)
summary(dent_past_mod3)  #not significant
Anova(dent_past_mod3, test.statistic = "Wald")



# Comparison of diluter vs control treatments for pasteuria
#  comparison of pulicaria genotype and controls only
dent_past_mod4 <- glm(Prevalence ~ Treatment, family = "binomial", weights = N, data = dentifera_past)
summary(dent_past_mod4)
Anova(dent_past_mod4, test.statistic = "Wald")
plot(dent_past_mod4)



# plot of predicted values of prevalence by diluter treatments vs controls
me6 <- ggpredict(dent_past_mod4, c("Treatment"))
plot(me6, add.data = F) + 
  labs(x = "Treatment", y = "Bacteria Prevalence (%)", title = NULL)


# comparison of prevalence in dentifera based on diluter genotype
dentifera_past_dilute <- dentifera_past %>%
  filter(Treatment == "diluter")

dent_past_mod7 <- glm(Prevalence ~ PulicariaLine, family = "binomial", weights = N, data = dentifera_past_dilute)
summary(dent_past_mod7)
e<- emmeans(dent_past_mod7, specs = pairwise ~ PulicariaLine, type = "response")
e


dim(dentifera_past)



### attempt to analyze the data in a different way...ALL THE SAME OUTCOME
data2 <- read.csv("DilutionDentiferaInfectionPrevalence_Combined.csv", stringsAsFactors = F, header = T)
head(data2)
pdata <- filter(data, Parasite == "Pasteuria")

model <- glm(Prevalence ~ PulicariaLine,family = "binomial", weights = N, data = pdata)
summary(model)
Anova(model)
f <- emmeans(model, specs = pairwise ~ PulicariaLine, type = "response")


data <- read.csv("DilutionDentiferaInfectionPrevalence_Expanded.csv", stringsAsFactors = F, header = T)
head(data)
pdata <- filter(data, Parasite == "Pasteuria")

model2 <- glm(Infected ~ PulicariaLine,family = "binomial", data = data)
summary(model2)
Anova(model2)
g <- emmeans(model2, specs = pairwise ~ PulicariaLine, type = "response")


model3 <- logistf(Infected ~ PulicariaLine, data = data)
summary(model3)
Anova(model3)

model4 <- flac(Infected ~ PulicariaLine, data = data)
summary(model4)
model4$linear.predictors

me_4 <- ggpredict(model4, "PulicariaLine")
plot(me4)
class(model4)





#### NPMANOVA example


## Load iris dataset:
data(iris)
## Scatterplot for two variables:
library(lattice)
xyplot(Sepal.Length ~ Sepal.Width, groups = Species, data = iris,
       pch = 16, auto.key = TRUE)

library(vegan)
Y <- iris[, c("Sepal.Length", "Sepal.Width")]
## Perform the NPMANOVA:
adonis2(Y ~ iris$Species, method = "euclidean")


# post-hoc test for adonis (need to update R for this to work)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)
pairwise.adonis(Y, iris$Species, sim.method = "euclidean",
                p.adjust.m = "bonferroni")

?pairwise.adonis

# try with dentifera past data
xyplot(Prevalence ~ BodySize_mm, groups = PulicariaLine2, data = dentifera_past_lines,
       pch = 16, auto.key = TRUE)

?vegan::betadisper()
?vegan::vegdist()


A <- dentifera_past_lines[ , c("Prevalence", "BodySize_mm")]
adonis(A ~ dentifera_past_lines$PulicariaLine2, method = "euclidean", permutations = 999)
pairwise.adonis(A, dentifera_past_lines$PulicariaLine2, sim.method = "bray", p.adjust.m = "bonferroni")


# this doesn't work
X <- dentifera_past_lines[ , "Prevalence"]
adonis(X ~ dentifera_past_lines$PulicariaLine2, method = "bray", permutations = 999)



# Try PERMANOVA with both past and metsch prevalence as response variables



dentifera_past_lines <- dplyr::rename(dentifera_past_lines, Past.Prev = Prevalence)
dentifera_metsch_lines <- dplyr::rename(dentifera_metsch_lines, Metsch.Prev = Prevalence)
dentifera_metsch_lines2 <- select(dentifera_metsch_lines, PulicariaLine2, Rep, Metsch.Prev)


new_data <- full_join(dentifera_past_lines, dentifera_metsch_lines2)

?xyplot
xyplot(Past.Prev ~ Metsch.Prev, groups = PulicariaLine2, data = new_data,
       pch = 16, auto.key = TRUE)

B <- new_data[ , c("Past.Prev", "Metsch.Prev")]
adonis(B ~ new_data$PulicariaLine2, method = "euclidean", permutations = 999)
pairwise.adonis(A, new_data$PulicariaLine2, sim.method = "euclidean", p.adjust.m = "bonferroni")


adonis(B ~ new_data$PulicariaLine2 * new_data$BodySize_mm, method = "euclidean", permutations = 999)
