# Area under the prevalence curve calculations for 2014 - 2018, 2020 - 2021

# This code generates the calculations for area under the prevalence curve for each year, lake, host species, and parasite combination from July to November.

# Written by: Michelle Fearon
# Last updated: 17 Dec 2022


# AUC OF PREVALENCE  ============================================================

library(tidyverse)
library(DescTools)
library(vegan)
library(here)

# set the path to the script relative to the project root directory
here::i_am("mi-fielddata-analysis/scripts/Area Under the Prevalence Curve 2014-2021.R")


loop_data <- read.csv(here("mi-fielddata-analysis/data/Clean-Data-2014-2021_All-Host-Densities.csv"), header = TRUE)
str(loop_data)


# diluter host densities (for Michelle's pastueria dilution project)
loop_data$Diluter.density <- loop_data$Pulicaria.density + loop_data$Retrocurva.density

# all non-dentifera densities
loop_data$OtherHost.density <- loop_data$Pulicaria.density + loop_data$Retrocurva.density + loop_data$Ceriodaphnia.density + loop_data$Dubia.density + 
  loop_data$Parvula.density + loop_data$Ambigua.density + loop_data$Mendotae.density

# diversity metrics per sampling date
loop_data$species.richness <- rowSums(select(loop_data, Dentifera.density:Mendotae.density) > 0)
loop_data$shannon <- diversity(select(loop_data, Dentifera.density:Mendotae.density), "shannon")


# gather into tall format; select appropriate columns; filter out early sampling dates so that all years start in July
loop_data <- loop_data %>%
  gather(Parasite.Species, Prevalence, pasteuria.prev:coke.prev) %>%
  select(Host.Species, Year, Lake, Parasite.Species, Julian.Day, Julian, Prevalence, Total, Total.Density:shannon) %>%
  filter(Julian > 181)

# arrange data 
loop_data <- loop_data %>%
  group_by(Host.Species, Year, Lake, Parasite.Species, Julian.Day) %>%
  dplyr::arrange( .by_group = T)


# set year and julian day as correct class
loop_data$Year <- as.factor(loop_data$Year)
loop_data$Julian.Day <- as.integer(loop_data$Julian.Day)

# order everything by julian day
#loop_data <- loop_data[order(loop_data$Julian.Day), ] #make sure julian days are in order

# make sure it's in data frame format
loop_data <- as.data.frame(loop_data)


##### REMOVE 2017 and 2018 sampling dates with missing density data
loop_data <- loop_data %>%
  filter(!is.na(Total.Density))


str(loop_data)
View(loop_data)

# write.csv(loop_data, "loop_data.csv", row.names = FALSE)
# AUC loop -------------------------------------------------------------------------

YEAR<-unique(loop_data$Year)

OUT <- NULL
for (i in YEAR) {
  thisyear <- loop_data[loop_data$Year == i, ] #look at 1 year at a time
  LAKE <- unique(thisyear$Lake)
  
  for (j in LAKE) {
    thislake <- thisyear[thisyear$Lake == j, ] #within a year, look at 1 lake at a time
    HOST <- unique(thislake$Host.Species)
    
    for (h in HOST) {
      thishost <- thislake[thislake$Host.Species == h, ] # within a lake, look at 1 host at a time
      PARA <- unique(thishost$Parasite.Species)
      
      for (p in PARA) {
        thisstuff <- thishost[thishost$Parasite.Species == p, ] #within that host, look at 1 parasite at a time.
        
        if(sum(!is.na(thisstuff$Prevalence)) > 0 ){  # if there is at least one Prevalence value that is not an NA, then calculate the max prev, etc
          maxdata<-filter(thisstuff, Prevalence == max(Prevalence, na.rm = T)) # creates a data frame with a single row if there is a max prevalence (mean of a single value is that value), or if prevalence is all zeros then code below calculates the mean across the season
          max.prev<-maxdata$Prevalence[1] #find the maximum prevalence of that parasite in the lake.
          count.at.max<-mean(maxdata$Total, na.rm = T) #If max.prev is 0, we will use the mean daphnia that were counted across the season.
          dentifera.at.max <- mean(maxdata$Dentifera.density, na.rm = T)
          pulicaria.at.max <- mean(maxdata$Pulicaria.density, na.rm = T)
          retrocurva.at.max <- mean(maxdata$Retrocurva.density, na.rm = T)
          diluter.at.max <- mean(maxdata$Diluter.density, na.rm = T)
          otherhost.at.max <- mean(maxdata$OtherHost.density, na.rm = T)
          richness.at.max <- mean(maxdata$species.richness, na.rm = T)
          shannon.at.max <- mean(maxdata$shannon, na.rm = T)
          
          if(sum(!is.na(thisstuff$Prevalence)) > 1 ){  # need at least 2 prevalence values to calculate the Area Under the Curve
          TOTS <- NULL #make an empty matrix to include sequential prevalences. This will be used to calculate integrated areas.
          
          for (m in 1:length(thisstuff$Prevalence)) {
            total <- thisstuff[m, "Prevalence"] #prevalence
            jul <- thisstuff[m, "Julian"] #julian day
            tots <- c(m, jul, total) #new row in dataframe with m, julian day, and total
            TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each lake/year/parasite)
            colnames(TOTS) <- c("m", "date", "total")
          }
          
          rownames(TOTS) <- NULL
          TOTS <- as.data.frame(TOTS)
          CURVE <- NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
          
          for (n in 1:(length(TOTS$total) - 1)) {
            curve <- 0.5*((TOTS[n, 3] + TOTS[n+1, 3])*(TOTS[n+1, 2] - TOTS[n, 2])) #trapezoid rule; this calculates area of each chunk
            CURVE <- c(CURVE, curve)
          } #write down each chunk
          
          area <- sum(CURVE, na.rm = T) #take the sum
          AREA <- AUC(TOTS[ , "date"], TOTS[ , "total"], method = "trapezoid", na.rm = T)
          }else{
            area <- NA
            AREA <- NA
          }
        } else{ # if all Prevalence values are NA (b/c no data collected), then put in an NA for max.prev and still calculate the mean densities
          max.prev<-NA #find the maximum prevalence of that parasite in the lake.
          count.at.max<-mean(thisstuff$Total, na.rm = T) #If max.prev is 0, we will use the mean daphnia that were counted across the season.
          dentifera.at.max <- mean(thisstuff$Dentifera.density, na.rm = T)
          pulicaria.at.max <- mean(thisstuff$Pulicaria.density, na.rm = T)
          retrocurva.at.max <- mean(thisstuff$Retrocurva.density, na.rm = T)
          diluter.at.max <- mean(thisstuff$Diluter.density, na.rm = T)
          otherhost.at.max <- mean(thisstuff$OtherHost.density, na.rm = T)
          richness.at.max <- mean(thisstuff$species.richness, na.rm = T)
          shannon.at.max <- mean(thisstuff$shannon, na.rm = T)
          
          #AUC under the curve is NA if there are no prevalence values for this lake, year, host, parasite combination
          area <- NA
          AREA <- NA
        }
        
        
        output <- c(i, j, h, p, max.prev, area, AREA, count.at.max, dentifera.at.max, pulicaria.at.max, 
                    retrocurva.at.max, diluter.at.max, otherhost.at.max, richness.at.max, shannon.at.max)
        OUT <- rbind(OUT, output)
      }
    }
  }
}



temp <- filter(loop_data, Year == 2015, Host.Species == "dentifera", Lake == "Whitmore", Parasite.Species == "pasteuria.prev")

if(sum(!is.na(temp$Prevalence)) > 0 ){  # if there is at least one Prevalence value that is not an NA, then calculate the max prev, etc
  maxdata<-filter(temp, Prevalence == max(Prevalence, na.rm = T))
  max.prev<-maxdata$Prevalence[1] #find the maximum prevalence of that parasite in the lake.
  count.at.max<-mean(maxdata$Total, na.rm = T) #If max.prev is 0, we will use the mean daphnia that were counted across the season.
  dentifera.at.max <- mean(maxdata$Dentifera.density, na.rm = T)
  pulicaria.at.max <- mean(maxdata$Pulicaria.density, na.rm = T)
  retrocurva.at.max <- mean(maxdata$Retrocurva.density, na.rm = T)
  diluter.at.max <- mean(maxdata$Diluter.density, na.rm = T)
  otherhost.at.max <- mean(maxdata$OtherHost.density, na.rm = T)
  richness.at.max <- mean(maxdata$species.richness, na.rm = T)
  shannon.at.max <- mean(maxdata$shannon, na.rm = T)
  
  if(sum(!is.na(temp$Prevalence)) > 1 ){  # need at least 2 prevlance values to calculate the Area Under the Curve
    TOTS <- NULL #make an empty matrix to include sequential prevalences. This will be used to calculate integrated areas.
    
    for (m in 1:length(temp$Prevalence)) {
      total <- temp[m, "Prevalence"] #prevalence
      jul <- temp[m, "Julian"] #julian day
      tots <- c(m, jul, total) #new row in dataframe with m, julian day, and total
      TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each lake/year/parasite)
      colnames(TOTS) <- c("m", "date", "total")
    }
    
    rownames(TOTS) <- NULL
    TOTS <- as.data.frame(TOTS)
    CURVE <- NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
    
    for (n in 1:(length(TOTS$total) - 1)) {
      curve <- 0.5*((TOTS[n, 3] + TOTS[n+1, 3])*(TOTS[n+1, 2] - TOTS[n, 2])) #trapezoid rule; this calculates area of each chunk
      CURVE <- c(CURVE, curve)
    } #write down each chunk
    
    area <- sum(CURVE, na.rm = T) #take the sum
    AREA <- AUC(TOTS[ , "date"], TOTS[ , "total"], method = "trapezoid", na.rm = T)
  }else{
    area <- NA
    AREA <- NA
  }
}else{ # if all Prevalence values are NA (b/c no data collected), then put in an NA for max.prev and still calculate the mean densities
  max.prev<-NA #find the maximum prevalence of that parasite in the lake.
  count.at.max<-mean(temp$Total, na.rm = T) #If max.prev is 0, we will use the mean daphnia that were counted across the season.
  dentifera.at.max <- mean(temp$Dentifera.density, na.rm = T)
  pulicaria.at.max <- mean(temp$Pulicaria.density, na.rm = T)
  retrocurva.at.max <- mean(temp$Retrocurva.density, na.rm = T)
  diluter.at.max <- mean(temp$Diluter.density, na.rm = T)
  otherhost.at.max <- mean(temp$OtherHost.density, na.rm = T)
  richness.at.max <- mean(temp$species.richness, na.rm = T)
  shannon.at.max <- mean(temp$shannon, na.rm = T)
  
  TOTS <- NULL #make an empty matrix to include sequential prevalences. This will be used to calculate integrated areas.
  
  for (m in 1:length(temp$Prevalence)) {
    total <- temp[m, "Prevalence"] #prevalence
    jul <- temp[m, "Julian"] #julian day
    tots <- c(m, jul, total) #new row in dataframe with m, julian day, and total
    TOTS <- rbind(TOTS, tots) #bind all rows together (different one for each lake/year/parasite)
    colnames(TOTS) <- c("m", "date", "total")
  }
  
  rownames(TOTS) <- NULL
  TOTS <- as.data.frame(TOTS)
  CURVE <- NULL #make an empty matrix for calulating each chunk of the integrated area by trapezoid rule.
  
  for (n in 1:(length(TOTS$total) - 1)) {
    curve <- 0.5*((TOTS[n, 3] + TOTS[n+1, 3])*(TOTS[n+1, 2] - TOTS[n, 2])) #trapezoid rule; this calculates area of each chunk
    CURVE <- c(CURVE, curve)
  } #write down each chunk
  
  area <- sum(CURVE) #take the sum
  AREA <- NA
}





colnames(OUT) <- c("Year", "Lake", "Host.Species", "Parasite.Species", "Max.Prevalence", "AUC.prev", "AUC.prev2", 
                   "Count.At.Max", 'dentifera.at.max', 'pulicaria.at.max', 'retrocurva.at.max', 'diluter.at.max', 
                   'otherhost.at.max', 'richness.at.max', 'shannon.at.max')
#write.table(OUT, "auc_14to17_prev.txt", col.names=TRUE, row.names=FALSE, sep="\t", append=FALSE) 


auc_data <- as.data.frame(OUT)
row.names(auc_data) <- c()


sapply(auc_data, class)
c <- c("Max.Prevalence", "AUC.prev", "AUC.prev2", "Count.At.Max", 'dentifera.at.max', 'pulicaria.at.max', 
       'retrocurva.at.max', 'diluter.at.max', 'otherhost.at.max', 'richness.at.max', 'shannon.at.max')
auc_data[c] <- sapply(auc_data[c], as.numeric)
#sapply(auc_data, class)


str(auc_data)
tail(auc_data)
View(auc_data)

unique(auc_data$Lake)


#check AUC calculations
auc_data_check <- auc_data %>%
  mutate(AUC_check = if_else(AUC.prev == AUC.prev2, T, F)) %>%
  filter(AUC_check == F)
      # it appears that Clara's calculation (AUC.prev) is different from the AUC function (AUC.prev2) when there is an
      # interruption in our detections of parasite prevalence. As of 12/17/2022, Meg recommends to go with the function
      # that produces a higher value due to calculation around the NAs with na.rm = T.


write.csv(auc_data, here("mi-fielddata-analysis/data/auc_14to21_prev.csv"), quote = F)

