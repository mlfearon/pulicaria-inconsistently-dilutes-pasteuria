# Varun's Age Experiment Analysis

# checks what your working directory is set to
getwd()

# set your working directory to where your data is
setwd("/Users/varunravichandran/Documents/Daphnia/Age Experiment Data")


#Michelle's working directory
setwd("C:/Users/mlfea/OneDrive/Documents/PROJECTS/MHMP Daphnia Duffy/Generality of dilution (MHMP)/Puliaria dilution paper/Varun Age Pulicaria Experiment")



# load libraries
library(dplyr)      # useful for manipulating data
library(ggplot2)    # visualizing data
library(ggeffects)  # visualizing model results
library(lme4)       # running models!!
library(emmeans)    # testing for differences among categories after you run a linear model
library(car)
library(MuMIn)      # useful for model selection (comparing multiple models)
library(colorBlindness)
library(gridExtra)
library(grid)

## updated overdispersion function from Ben Bolker
overdisp_fun <- function(model) {
  rdf <- df.residual(model)
  rp <- residuals(model,type="pearson")
  Pearson.chisq <- sum(rp^2)
  prat <- Pearson.chisq/rdf
  pval <- pchisq(Pearson.chisq, df=rdf, lower.tail=FALSE)
  c(chisq=Pearson.chisq,ratio=prat,rdf=rdf,p=pval)
}
# this function is an important test for glms to ensure model fit
# if the p-value is significant then the model is overdispersed and we will need to fix that.


# this is a colorblind friendly color pallette 
# color blind palette (from https://zenodo.org/record/3381072#.X8VKS2hKjZt)
Tol_muted <- c('#882255','#44AA99', '#117733', '#332288', '#DDCC77','#CC6677', '#88CCEE', '#AA4499', '#DDDDDD', '#999933')
# Goedhart, Joachim. (2019, August 29). Material related to the blog "Dataviz with Flying Colors". Zenodo. http://doi.org/10.5281/zenodo.3381072
Okabe_Ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000")


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


####################################################
# Read in our data and finish cleaning our data
####################################################

# Read in the  Gut Spores and bodysize data and name it
gutspores <- read.csv("GutSporesData.csv", stringsAsFactors = F, header = T)

# look at first 6 rows
head(gutspores)

# look at last 6 rows
tail(gutspores)


# There are a bunch of NAs in the data set that we want to remove before working with it
gutspores <- select(gutspores, line:punctured)  # select function allows us to select the columns that we want by their names

# check that it worked
head(gutspores)




# look at structure of data set (check that all of your variables came in with the right class)
str(gutspores)

# the variables all look to be the right class, but if we needed to fix them then we would 
# use the as.factor(), as.numeric(), or as.character() functions to convert the class


# This is an example of what that would look like
#gutspores$bodysize <- as.numeric(gutspores$bodysize)



# next we need to check that all of our categorical variables are all correct (no typos or extra spaces in any row)

unique(gutspores$age)  # looks like we have a typo in young as "yough" in our data set that we will have to fix

unique(gutspores$line) # looks good

unique(gutspores$dose) # looks good


# now we need to fix the one spot that says "yough" instead of "young"
gutspores$age[gutspores$age == "yough"] <- "young"

# check that it was fixed
unique(gutspores$age) # it looks good now :)



# convert our categories to factors
gutspores$line <- factor(gutspores$line, levels = c("BA", "Clover", "Mid67", "Pine"))
gutspores$age <- factor(gutspores$age, levels = c("young", "middle", "old"))
gutspores$dose <- factor(gutspores$dose, levels = c("high", "low", "control"))

# calculate the total number of spores
gutspores <- mutate(gutspores, total_spores = embedded + punctured)

head(gutspores)



# split data based on dose treatment (control vs exposed)
gutspores_exposed <- filter(gutspores, dose == "high")
head(gutspores_exposed)



## a simple model for embedded spores based on age and line
model <- glm(embedded ~ age * line, family = "poisson", data = gutspores_exposed)
summary(model)
Anova(model)  # the interaction is not significant


# simplify model to remove the interaction
model2 <- glm(embedded ~ age + line, family = "poisson", data = gutspores_exposed)
summary(model2)
Anova(model2)

# this is a post-hoc test that allows us to determine if there are significant differences among our categories
a <- emmeans(model2, specs = pairwise ~ line | age, type = "response")
a

# visualize the model predictions
me_embed <- ggpredict(model2, c("age", "line"))
plot(me_embed)

# same plot with the age and pulicaria line flipped
me_embed2 <- ggpredict(model2, c("line", "age"))
plot(me_embed2)



# more complex model with age, line, bodysize, and gut cell width
model3 <- glm(embedded ~ gutcellwidth + bodysize + age + line, family = "poisson", data = gutspores_exposed)
summary(model3)
Anova(model3)

model20 <- glm(embedded ~ gutcellwidth + bodysize + line, family = "poisson", data = gutspores_exposed)
summary(model20)
Anova(model20)

model21 <- glm(embedded ~ gutcellwidth + age + line, family = "poisson", data = gutspores_exposed)
summary(model21)
Anova(model21)


# model selection
out.put_subset <- model.sel(model, model2, model20, model21)
out.put_subset
write.csv(out.put_subset, "Model Selection Table_subset.csv")








# embedded spore data
model4 <- glm(punctured ~ age * line, family = "poisson", data = gutspores_exposed)
summary(model4)

model5 <- glm(punctured ~ age * line, family = "poisson", data = gutspores_exposed)
summary(model5)
Anova(model5)

b <- emmeans(model5, specs = pairwise ~ age | line, type = "response")
b

me_punct <- ggpredict(model5, c("age", "line"))
plot(me_punct)

me_punct2 <- ggpredict(model5, c("line", "age"))
plot(me_punct2)

model6 <- glm(punctured ~ gutcellwidth + bodysize + age + line, family = "poisson", data = gutspores_exposed)
summary(model6)

model7 <- glm(punctured ~ gutcellwidth * bodysize, family = "poisson", data = gutspores_exposed)
summary(model7)

model8 <- glm(punctured ~ bodysize * gutcellwidth, family = "poisson", data = gutspores_exposed)
summary(model8)

model9 <- glm(embedded ~ bodysize * gutcellwidth, family = "poisson", data = gutspores_exposed)
summary(model9)

model10 <- glm(embedded ~ gutcellwidth * bodysize, family = "poisson", data = gutspores_exposed)
summary(model10)

model11 <- glm(total_spores ~ age * line, family = "poisson", data = gutspores_exposed)
summary(model11)

me_total <- ggpredict(model11, c("age", "line"))
plot(me_total)

model12 <- glm(total_spores ~ age, family = "poisson", data = gutspores_exposed)
summary(model12)

model22 <- glm(total_spores ~ age * gutcellwidth, family = "poisson", data = gutspores_exposed)
summary(model22)

model23 <- glm(total_spores ~ gutcellwidth + age + line, family = "poisson", data = gutspores_exposed)
summary(model23)
Anova(model23)





# check if bodysize varies with age or line
model24 <- lm(bodysize ~ age * line, data = gutspores)
summary(model24)
Anova(model24)
plot(model24) # residual plots don't look great, need to fix!!!!

X <- emmeans(model24, specs = pairwise ~ line | age, type = "response")
X

# create a data frame with the significant differences among lines within each age group
sig_labels <- data.frame(age = c(rep("young", 4), rep("middle", 4), rep("old", 4)), line = rep(c("BA", "Clover", "Mid67", "Pine"),3), 
                         height = c(2000,1800,1800,1700,2050,2050,1950,1775,2150,2150,2050,1900), 
                         label = c("a","b","b","c","a", "ab", "b", "c", "a", "a", "b", "c"))




hist(sqrt(gutspores$bodysize))
gutspores <- as.data.frame(gutspores)

# make a figure that shows how body size varies among age groups and genotypes
ggplot(gutspores, aes(x = age, y= bodysize, color = line)) +
  geom_boxplot() +
  geom_jitter(aes(color = line), position = position_jitterdodge(0.2), alpha = 0.5) +
  labs(x = "Age", y = bquote("Body Size ("~ mu~"m)"), title = NULL) +
  ylim(1400,2200) +
  geom_text(data = sig_labels, aes(x = age, y = height, color = line), label = c("a","b","b","c","a", "ab", "b", "c", "a", "a", "b", "c"), position = position_dodge(width = 0.75), show.legend = F, size = 10/.pt) +
  scale_color_manual(values = Tol_muted, name = "Genotype", labels = c("BA", "Clover", "Mid67", "Pine")) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))
ggsave("Bodysize vs age and line.tiff", dpi = 600, width = 5, height = 3.5, units = "in", compression="lzw") # save the figure



# this is another way to show that same data based on the predictions from the model above
    # it shows pretty much the same thing, but I'm having trouble getting the points to line up correctly for some reason
    # so we may want to go with the one above for now.
me_bodysize <- ggpredict(model24, c("age", "line"))
View(me_bodysize)
me_bodysize$x <- factor(me_bodysize$x, levels = c("young", "middle", "old"))
me_bodysize$group <- factor(me_bodysize$group, levels = c("BA", "Clover", "Mid67", "Pine"))
size <- plot(me_bodysize, add.data = T, jitter = F, dodge = T) +
  labs(x = "Age", y = bquote("Body Size ("~ mu~"m)"), title = NULL) +
  ylim(1400,2200) +
  scale_color_manual(values = Tol_muted, name = "Genotype", labels = c("BA", "Clover", "Mid67", "Pine")) +
  #scale_shape_discrete(name = "Species", labels = c(bquote(italic("Apis mellifera")), bquote(italic("Bombus impatiens")), 
  #                                                  bquote(italic("Lasioglossum")~ "spp."), bquote(italic("Eucera pruinosa")))) +
  #geom_text(aes(y = predicted + 0.13), label = c("a","b","c","c","a", "b", "c", "c", "a", "a", "b", "b"), position = position_dodge(width = 0.4), show.legend = F, size = 9/.pt) +
  theme_classic() +
  theme(axis.text = element_text(size = 10, color = "black"), axis.title = element_text(size = 11, color = "black"), 
        legend.text = element_text(size = 10), legend.title = element_text(size = 10))

print(size)
ggsave("Bodysize vs age and line2.tiff", plot = size, dpi = 600, width = 5, height = 3.5, units = "in", compression="lzw")







# color blind check
cvdPlot(prev_host)
  



# check if gut cell width varies with age or line
model25 <- lm(gutcellwidth ~ age * line, data = gutspores)
summary(model25)
plot(model25) 

Anova(model25)

hist(gutspores$gutcellwidth)



me_gut <- ggpredict(model25, c("age", "line"))
plot(me_gut, add.data = TRUE)

y <- emmeans(model25, specs = pairwise ~ line | age, type = "response")
y


# now you can use the example above to make a similar figure but instead of bodysize on the y axis
# we should have gut cell width in micrometers













# Quick test if bodysize correlates with gut cell width
cor.test(gutspores$bodysize, gutspores$gutcellwidth)
# no they are not significantly correlated,-0.14







# Now we need to read in the infection prevalence data and name it
InfectionBinary <- read.csv("InfectionBinary.csv", stringsAsFactors = F, header = T)

head(InfectionBinary)
tail(InfectionBinary)


str(InfectionBinary)

# need to remove the NA's from the data set
InfectionBinary <- InfectionBinary %>%
  select(line:status) %>%
  filter(!is.na(rep))



# check that each category is correct
unique(InfectionBinary$line)
unique(InfectionBinary$dose) # looks like there are some high treatments that have an extra space that needs to be fixed
unique(InfectionBinary$age)


InfectionBinary$dose[ InfectionBinary$dose == "high " ] <- "high"
unique(InfectionBinary$dose)  # now dose is fixed

# convert our categories to factors
InfectionBinary$line <- factor(InfectionBinary$line, levels = c("BA", "Clover", "Mid67", "Pine"))
InfectionBinary$age <- factor(InfectionBinary$age, levels = c("young", "middle", "old"))
InfectionBinary$dose <- factor(InfectionBinary$dose, levels = c("high", "low", "control"))



## add body size data and gut cell width to the data set




inf_mod <- glm(Metsch ~ age + line + dose, family = "binomial", data = InfectionBinary)
summary(inf_mod)
Anova(inf_mod)


X <- emmeans(inf_mod, specs = pairwise ~ line | age, type = "response")
X

me_inf <- ggpredict(inf_mod, c("age", "line"))
plot(me_inf, add.data = F)


inf_mod <- glm(Metsch ~ bodysize * gutcellwidth, family = "binomial", data = InfectionBinary)
summary(inf_mod)
Anova(inf_mod)






# read in the Spore Count data and name it
SporeCount <- read.csv("SporeCountData.csv", stringsAsFactors = F, header = T)

head(SporeCount)
tail(SporeCount)

str(SporeCount)

# need to remove the NA's from the data set
SporeCount <- SporeCount %>%
  select(line:mature_avg) %>%
  filter(!is.na(beaker))



# check that each category is correct
unique(SporeCount$line)
unique(SporeCount$dose)
unique(SporeCount$age)


# convert our categories to factors
SporeCount$line <- factor(SporeCount$line, levels = c("BA", "Clover", "Mid67", "Pine"))
SporeCount$age <- factor(SporeCount$age, levels = c("young", "middle", "old"))
SporeCount$dose <- factor(SporeCount$dose, levels = c("high", "low", "control"))


# remove samples without spores
SporeCount_subset <- SporeCount %>%
  filter(spores == "Y") %>%
  filter(bud_avg < 1000) %>%  # remove a couple of outliers in the data that I don't fully trust
  mutate(total_spores = bud_avg + immature_avg + mature_avg, bud_per_animal = (bud_avg * 10000*0.01), immature_per_animal = (immature_avg * 10000*0.01), 
         mature_per_animal = (mature_avg * 10000*0.01), total_per_animal = bud_per_animal + immature_per_animal + mature_per_animal)

head(SporeCount_subset)

# plot of budding spores per animal by age and genotype
bud <- ggplot(data = SporeCount_subset, aes(x = age, y = bud_per_animal, color = line)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Budding Spores per Animal") +
  scale_color_manual(values = Tol_muted, name = "Genotype", labels = c("BA", "Clover", "Mid67", "Pine")) +
  theme_classic()
print(bud)
ggsave("BuddingSpores vs age and genotype.tiff", plot = bud, dpi = 600, width =4.5, height = 3.5, units = "in", compression="lzw")

# plot of immature spores per animal by age and genotype
immature <- ggplot(data = SporeCount_subset, aes(x = age, y = immature_per_animal, color = line)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Immature Spores per Animal") +
  scale_color_manual(values = Tol_muted, name = "Genotype", labels = c("BA", "Clover", "Mid67", "Pine")) +
  theme_classic()
print(immature)
ggsave("ImmatureSpores vs age and genotype.tiff", plot = immature, dpi = 600, width = 4.5, height = 3.5, units = "in", compression="lzw")


# plot of mature spores per animal by age and genotype
mature <- ggplot(data = SporeCount_subset, aes(x = age, y = mature_per_animal, color = line)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Mature Spores per animal") +
  scale_color_manual(values = Tol_muted, name = "Genotype", labels = c("BA", "Clover", "Mid67", "Pine")) +
  theme_classic()
print(mature)
ggsave("MatureSpores vs age and genotype.tiff", plot = mature, dpi = 600, width = 5, height = 3.5, units = "in", compression="lzw")


# plot of total spores per animal by age and genotype
total <- ggplot(data = SporeCount_subset, aes(x = age, y = total_per_animal, color = line)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Total Spores per animal") +
  scale_color_manual(values = Tol_muted, name = "Genotype", labels = c("BA", "Clover", "Mid67", "Pine")) +
  theme_classic()
print(total)
ggsave("TotalSpores vs age and genotype.tiff", plot = total, dpi = 600, width = 5, height = 3.5, units = "in", compression="lzw")


# combine all four plots into one
all_spores <- grid_arrange_shared_legend(bud, immature, mature, total, nrow = 2, ncol = 2, position = "right")
ggsave("Spores vs age and genotype_four plots.tiff", plot = all_spores, dpi = 600, width = 7, height = 6.5, units = "in", compression="lzw")







# look at differences in only Mid67
SporeCount_Mid67 <- SporeCount_subset %>%
  filter(line == "Mid67")


# plot of budding spores per animal by age and dose in Mid67 only
bud_mid67 <- ggplot(data = SporeCount_Mid67, aes(x = age, y = bud_per_animal, color = dose)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Budding Spores per Animal") +
  scale_color_manual(values = Okabe_Ito[5:6], name = "Dose") +
  ylim(0,30000) +
  theme_classic()
print(bud_mid67)
ggsave("Mid67_BuddingSpores vs age and dose.tiff", plot= bud_mid67, dpi = 600, width =4.5, height = 3.5, units = "in", compression="lzw")

# plot of Immature spores per animal by age and dose in Mid67 only
immature_mid67 <- ggplot(data = SporeCount_Mid67, aes(x = age, y = immature_per_animal, color = dose)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Immature Spores per Animal") +
  scale_color_manual(values = Okabe_Ito[5:6], name = "Dose") +
  theme_classic()
print(immature_mid67)
ggsave("Mid67_ImmatureSpores vs age and dose.tiff", plot = immature_mid67, dpi = 600, width = 4.5, height = 3.5, units = "in", compression="lzw")


# plot of mature spores per animal by age and dose in Mid67 only
mature_mid67 <- ggplot(data = SporeCount_Mid67, aes(x = age, y = mature_per_animal, color = dose)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Mature Spores per animal") +
  scale_color_manual(values = Okabe_Ito[5:6], name = "Dose") +
  theme_classic()
print(mature_mid67)
ggsave("Mid67_MatureSpores vs age and dose.tiff", plot = mature_mid67, dpi = 600, width = 5, height = 3.5, units = "in", compression="lzw")


# plot of total spores per animal by age and dose in Mid67 only
total_mid67 <- ggplot(data = SporeCount_Mid67, aes(x = age, y = total_per_animal, color = dose)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(), alpha = 0.5) +
  labs(x = "Age", y = "Number of Total Spores per animal") +
  scale_color_manual(values = Okabe_Ito[5:6], name = "Dose") +
  theme_classic()
print(total_mid67)
ggsave("Mid67_MatureSpores vs age and dose.tiff", plot = total_mid67, dpi = 600, width = 5, height = 3.5, units = "in", compression="lzw")


# combine all four plots into one
Mid67_spores <- grid_arrange_shared_legend(bud_mid67, immature_mid67, mature_mid67, total_mid67, nrow = 2, ncol = 2, position = "right")
ggsave("Mid67_Spores vs age and genotype_four plots.tiff", plot = Mid67_spores, dpi = 600, width = 7, height = 6.5, units = "in", compression="lzw")
