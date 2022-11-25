library(ggplot2)
library("readxl")
library(tidyverse)
library(dplyr)
library(ggpubr)
library(cowplot)
library(ggthemes)
library(viridis)
library(ggrepel)
library('gridExtra')
library('data.table')
library('scales')
library('lubridate')
library("patchwork")
library("zoo")
library(sp)
library(rworldmap)


world_data<-getMap(resolution='low')@data

country2continent = function(x)
{  
  country_data = subset(world_data, NAME==x)
  
  return (as.character(country_data[['REGION']]))   # returns the continent (7 continent model)
}


country2continent_region = function(x)
{  
  country_data = subset(world_data, NAME==x)
  
  return (as.character(country_data[['IMAGE24']]))  
}


country2lat = function(x)
{  
  country_data = subset(world_data, NAME==x)
  
  return (as.numeric(country_data[['LAT']]))  
}


country2long = function(x)
{  
  country_data = subset(world_data, NAME==x)
  
  return (as.numeric(country_data[['LON']]))  
}




export2type = function(x,y)
{  
  if (x==y) {
    
    return (as.character("Regional"))
  }
  
  else {
    
    return (as.character("Global"))
  }
}

stat_box_data <- function(y, upper_limit = max(iris$Sepal.Length) * 1.15) {
  return( 
    data.frame(
      y = 350,
      label = paste('n=',length(y),' ',
                    'mean delay=',round(mean(y)),'days'
      )
    )
  )
}


############## Reading Data #####################


#Setting up VOC TMRCAs from published literature - median and lower/upper bound of 95% confidence interval reported

alpha_first_date<-as.Date("2020-08-28", format="%Y-%m-%d") #Alpha (Hill et al. 2022) (2020-08-15 - 2020-09-09)
beta_first_date<-as.Date("2020-08-05", format="%Y-%m-%d") #(Tegally et al. 2021) (2020-07-15 - 2020-08-30)
delta_first_date<-as.Date("2020-10-19", format="%Y-%m-%d") #(McCrone et al. 2022) (2020-09-06 - 2020-11-29)
gamma_first_date<-as.Date("2020-11-15", format="%Y-%m-%d") #(Faria et al. 2021) (2020-10-06 – 2020-11-24)
ba1_first_date<-as.Date("2021-10-09", format="%Y-%m-%d") #(Viana et al. 2022) (2021-09-30 – 2021-10-20)
ba2_first_date<-as.Date("2021-11-05", format="%Y-%m-%d") #(Tegally et al. 2022 - BA4/5 paper) (2021-10-09 – 2021-11-29)
ba4_5_first_date<-as.Date("2021-12-15", format="%Y-%m-%d") #(Tegally et al. 2022 - BA4/5 paper) (2021-11-25 – 2022-02-06)

alpha_first_date_lo<-as.Date("2020-08-15", format="%Y-%m-%d") #Alpha (Hill et al. 2022) (2020-08-15 - 2020-09-09)
alpha_first_date_hi<-as.Date("2020-09-09", format="%Y-%m-%d") #Alpha (Hill et al. 2022) (2020-08-15 - 2020-09-09)
beta_first_date_lo<-as.Date("2020-07-15", format="%Y-%m-%d") #(Tegally et al. 2021) (2020-07-15 - 2020-08-30)
beta_first_date_hi<-as.Date("2020-08-30", format="%Y-%m-%d") #(Tegally et al. 2021) (2020-07-15 - 2020-08-30)
delta_first_date_lo<-as.Date("2020-09-06", format="%Y-%m-%d") #(McCrone et al. 2022) (2020-09-06 - 2020-11-29)
delta_first_date_hi<-as.Date("2020-11-29", format="%Y-%m-%d") #(McCrone et al. 2022) (2020-09-06 - 2020-11-29)
gamma_first_date_lo<-as.Date("2020-10-06", format="%Y-%m-%d") #(Faria et al. 2021) (2020-10-06 – 2020-11-24)
gamma_first_date_hi<-as.Date("2020-11-24", format="%Y-%m-%d") #(Faria et al. 2021) (2020-10-06 – 2020-11-24)
ba1_first_date_lo<-as.Date("2021-09-30", format="%Y-%m-%d") #(Viana et al. 2022) (2021-09-30 – 2021-10-20)
ba1_first_date_hi<-as.Date("2021-10-20", format="%Y-%m-%d") #(Viana et al. 2022) (2021-09-30 – 2021-10-20)
ba2_first_date_lo<-as.Date("2021-10-09", format="%Y-%m-%d") #(Tegally et al. 2022 - BA4/5 paper) (2021-10-09 – 2021-11-29)
ba2_first_date_hi<-as.Date("2021-11-29", format="%Y-%m-%d") #(Tegally et al. 2022 - BA4/5 paper) (2021-10-09 – 2021-11-29)
ba4_5_first_date_lo<-as.Date("2021-11-25", format="%Y-%m-%d") #(Tegally et al. 2022 - BA4/5 paper) (2021-11-25 – 2022-02-06)
ba4_5_first_date_hi<-as.Date("2022-02-06", format="%Y-%m-%d") #(Tegally et al. 2022 - BA4/5 paper) (2021-11-25 – 2022-02-06)




#Reading in all data for import/export analysis
#This step can take some time (reading and annotating large data files) - Between 10 and 30 mins depending on your resources

# Alpha

replicate1<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed1234.csv', sep = '\t', header = TRUE)
replicate1$replicate <- "1"
replicate2<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed1898.csv', sep = '\t', header = TRUE)
replicate2$replicate <- "2"
replicate3<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed2007.csv', sep = '\t', header = TRUE)
replicate3$replicate <- "3"
replicate4<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed2498.csv', sep = '\t', header = TRUE)
replicate4$replicate <- "4"
replicate5<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed3107.csv', sep = '\t', header = TRUE)
replicate5$replicate <- "5"
replicate6<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed4321.csv', sep = '\t', header = TRUE)
replicate6$replicate <- "6"
replicate7<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed4891.csv', sep = '\t', header = TRUE)
replicate7$replicate <- "7"
replicate8<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed5327.csv', sep = '\t', header = TRUE)
replicate8$replicate <- "8"
replicate9<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed6404.csv', sep = '\t', header = TRUE)
replicate9$replicate <- "9"
replicate10<-read.table(file='../ViralTransitions/Alpha/annottated_tree_events_seed7102.csv', sep = '\t', header = TRUE)
replicate10$replicate <- "10"

alpha_all_replicates<-rbind(replicate1,replicate2,replicate3,replicate4,
                            replicate5,replicate6,replicate7,replicate8,
                            replicate9,replicate10)
alpha_all_replicates[alpha_all_replicates == "USA"] <- "United States"
alpha_all_replicates[alpha_all_replicates == "Republic of the Congo"] <- "Congo (Brazzaville)"
alpha_all_replicates[alpha_all_replicates == "Democratic Republic of the Congo"] <- "Congo (Kinshasa)"

alpha_all_replicates[alpha_all_replicates == "Coted Ivoire"] <- "Ivory Coast"

alpha_all_replicates<-subset(alpha_all_replicates,Origin!='Reunion')
alpha_all_replicates<-subset(alpha_all_replicates,Origin!='Mayotte')

alpha_all_replicates$date<-date_decimal(alpha_all_replicates$EventTime)

alpha_all_replicates$Origin_Continent<-lapply(alpha_all_replicates$Origin,country2continent)
alpha_all_replicates$Destination_Continent<-lapply(alpha_all_replicates$Destination,country2continent)

alpha_all_replicates$Origin_Continent_Region<-lapply(alpha_all_replicates$Origin,country2continent_region)
alpha_all_replicates$Destination_Continent_Region<-lapply(alpha_all_replicates$Destination,country2continent_region)

alpha_all_replicates$Variant <- "Alpha"

alpha_all_replicates$days<-as.Date(cut(alpha_all_replicates$date,breaks = "day",start.on.monday = FALSE))
alpha_all_replicates$date<-as.Date(cut(alpha_all_replicates$date,breaks = "week",start.on.monday = FALSE))
alpha_all_replicates$date2<-as.Date(cut(alpha_all_replicates$date,breaks = "2 week",start.on.monday = FALSE))
alpha_all_replicates$date4<-as.Date(cut(alpha_all_replicates$date,breaks = "1 month",start.on.monday = FALSE))




# Beta

replicate1<-read.table(file='../ViralTransitions/Beta/seed1344_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate1$replicate <- "1"
replicate2<-read.table(file='../ViralTransitions/Beta/seed2007_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate2$replicate <- "2"
replicate3<-read.table(file='../ViralTransitions/Beta/seed2099_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate3$replicate <- "3"
replicate4<-read.table(file='../ViralTransitions/Beta/seed3278_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate4$replicate <- "4"
replicate5<-read.table(file='../ViralTransitions/Beta/seed3500_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate5$replicate <- "5"
replicate6<-read.table(file='../ViralTransitions/Beta/seed3989_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate6$replicate <- "6"
replicate7<-read.table(file='../ViralTransitions/Beta/seed4034_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate7$replicate <- "7"
replicate8<-read.table(file='../ViralTransitions/Beta/seed4704_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate8$replicate <- "8"
replicate9<-read.table(file='../ViralTransitions/Beta/seed5342_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate9$replicate <- "9"
replicate10<-read.table(file='../ViralTransitions/Beta/seed6767_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate10$replicate <- "10"



beta_all_replicates<-rbind(replicate1,replicate2,replicate3,replicate4,
                           replicate5,replicate6,replicate7,replicate8,
                           replicate9,replicate10)
beta_all_replicates[beta_all_replicates == "USA"] <- "United States"
beta_all_replicates[beta_all_replicates == "Republic of the Congo"] <- "Congo (Brazzaville)"
beta_all_replicates[beta_all_replicates == "Democratic Republic of the Congo"] <- "Congo (Kinshasa)"

beta_all_replicates[beta_all_replicates == "Coted Ivoire"] <- "Ivory Coast"

beta_all_replicates<-subset(beta_all_replicates,Origin!='Reunion')
beta_all_replicates<-subset(beta_all_replicates,Origin!='Mayotte')

beta_all_replicates$date<-date_decimal(beta_all_replicates$EventTime)

beta_all_replicates$Origin_Continent<-lapply(beta_all_replicates$Origin,country2continent)
beta_all_replicates$Destination_Continent<-lapply(beta_all_replicates$Destination,country2continent)

beta_all_replicates$Origin_Continent_Region<-lapply(beta_all_replicates$Origin,country2continent_region)
beta_all_replicates$Destination_Continent_Region<-lapply(beta_all_replicates$Destination,country2continent_region)

beta_all_replicates$Variant <- "Beta"

beta_all_replicates$days<-as.Date(cut(beta_all_replicates$date,breaks = "day",start.on.monday = FALSE))
beta_all_replicates$date<-as.Date(cut(beta_all_replicates$date,breaks = "week",start.on.monday = FALSE))
beta_all_replicates$date2<-as.Date(cut(beta_all_replicates$date,breaks = "2 week",start.on.monday = FALSE))
beta_all_replicates$date4<-as.Date(cut(beta_all_replicates$date,breaks = "1 month",start.on.monday = FALSE))



# Gamma

replicate1<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed1298.csv', sep = '\t', header = TRUE)
replicate1$replicate <- "1"
replicate2<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed2007.csv', sep = '\t', header = TRUE)
replicate2$replicate <- "2"
replicate3<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed3007.csv', sep = '\t', header = TRUE)
replicate3$replicate <- "3"
replicate4<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed3334.csv', sep = '\t', header = TRUE)
replicate4$replicate <- "4"
replicate5<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed4010.csv', sep = '\t', header = TRUE)
replicate5$replicate <- "5"
replicate6<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed4507.csv', sep = '\t', header = TRUE)
replicate6$replicate <- "6"
replicate7<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed4987.csv', sep = '\t', header = TRUE)
replicate7$replicate <- "7"
replicate8<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed5346.csv', sep = '\t', header = TRUE)
replicate8$replicate <- "8"
replicate9<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed6129.csv', sep = '\t', header = TRUE)
replicate9$replicate <- "9"
replicate10<-read.table(file='../ViralTransitions/Gamma/annottated_tree_events_seed6667.csv', sep = '\t', header = TRUE)
replicate10$replicate <- "10"

gamma_all_replicates<-rbind(replicate1,replicate2,replicate3,replicate4,
                            replicate5,replicate6,replicate7,replicate8,
                            replicate9,replicate10)
gamma_all_replicates[gamma_all_replicates == "USA"] <- "United States"
gamma_all_replicates[gamma_all_replicates == "Republic of the Congo"] <- "Congo (Brazzaville)"
gamma_all_replicates[gamma_all_replicates == "Democratic Republic of the Congo"] <- "Congo (Kinshasa)"

gamma_all_replicates[gamma_all_replicates == "Coted Ivoire"] <- "Ivory Coast"

gamma_all_replicates<-subset(gamma_all_replicates,Origin!='Reunion')
gamma_all_replicates<-subset(gamma_all_replicates,Origin!='Mayotte')

gamma_all_replicates$date<-date_decimal(gamma_all_replicates$EventTime)

gamma_all_replicates$Origin_Continent<-lapply(gamma_all_replicates$Origin,country2continent)
gamma_all_replicates$Destination_Continent<-lapply(gamma_all_replicates$Destination,country2continent)

gamma_all_replicates$Origin_Continent_Region<-lapply(gamma_all_replicates$Origin,country2continent_region)
gamma_all_replicates$Destination_Continent_Region<-lapply(gamma_all_replicates$Destination,country2continent_region)

gamma_all_replicates$Variant <- "Gamma"

gamma_all_replicates$days<-as.Date(cut(gamma_all_replicates$date,breaks = "day",start.on.monday = FALSE))
gamma_all_replicates$date<-as.Date(cut(gamma_all_replicates$date,breaks = "week",start.on.monday = FALSE))
gamma_all_replicates$date2<-as.Date(cut(gamma_all_replicates$date,breaks = "2 week",start.on.monday = FALSE))
gamma_all_replicates$date4<-as.Date(cut(gamma_all_replicates$date,breaks = "1 month",start.on.monday = FALSE))


## Delta



replicate1<-read.table(file='../ViralTransitions/Delta/delta0546_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate1$replicate <- "1"
replicate2<-read.table(file='../ViralTransitions/Delta/delta1321_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate2$replicate <- "2"
replicate3<-read.table(file='../ViralTransitions/Delta/voc/delta2765_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate3$replicate <- "3"
replicate4<-read.table(file='../ViralTransitions/Delta/delta3876_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate4$replicate <- "4"
replicate5<-read.table(file='../ViralTransitions/Delta/delta4012_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate5$replicate <- "5"
replicate6<-read.table(file='../ViralTransitions/Delta/delta5555_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate6$replicate <- "6"
replicate7<-read.table(file='../ViralTransitions/Delta/delta6667_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate7$replicate <- "7"
replicate8<-read.table(file='../ViralTransitions/Delta/delta7776_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate8$replicate <- "8"
replicate9<-read.table(file='../ViralTransitions/Delta/delta8878_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate9$replicate <- "9"
replicate10<-read.table(file='../ViralTransitions/Delta/delta9898_voc_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate10$replicate <- "10"


delta_all_replicates<-rbind(replicate1,replicate2,replicate3,replicate4,
                            replicate5,replicate6,replicate7,replicate8,
                            replicate9,replicate10)
delta_all_replicates[delta_all_replicates == "USA"] <- "United States"
delta_all_replicates[delta_all_replicates == "Republic of the Congo"] <- "Congo (Brazzaville)"
delta_all_replicates[delta_all_replicates == "Democratic Republic of the Congo"] <- "Congo (Kinshasa)"

delta_all_replicates[delta_all_replicates == "Coted Ivoire"] <- "Ivory Coast"

delta_all_replicates<-subset(delta_all_replicates,Origin!='Reunion')
delta_all_replicates<-subset(delta_all_replicates,Origin!='Mayotte')

delta_all_replicates$date<-date_decimal(delta_all_replicates$EventTime)

delta_all_replicates$Origin_Continent<-lapply(delta_all_replicates$Origin,country2continent)
delta_all_replicates$Destination_Continent<-lapply(delta_all_replicates$Destination,country2continent)

delta_all_replicates$Origin_Continent_Region<-lapply(delta_all_replicates$Origin,country2continent_region)
delta_all_replicates$Destination_Continent_Region<-lapply(delta_all_replicates$Destination,country2continent_region)

delta_all_replicates$Variant <- "Delta"

delta_all_replicates$days<-as.Date(cut(delta_all_replicates$date,breaks = "day",start.on.monday = FALSE))
delta_all_replicates$date<-as.Date(cut(delta_all_replicates$date,breaks = "week",start.on.monday = FALSE))
delta_all_replicates$date2<-as.Date(cut(delta_all_replicates$date,breaks = "2 week",start.on.monday = FALSE))
delta_all_replicates$date4<-as.Date(cut(delta_all_replicates$date,breaks = "1 month",start.on.monday = FALSE))

## BA.1

replicate1<-read.table(file='../ViralTransitions/BA1/voc_seed0987_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate1$replicate <- "1"
replicate2<-read.table(file='../ViralTransitions/BA1/voc_seed1098_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate2$replicate <- "2"
replicate3<-read.table(file='../ViralTransitions/BA1/voc_seed2109_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate3$replicate <- "3"
replicate4<-read.table(file='../ViralTransitions/BA1/voc_seed3210_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate4$replicate <- "4"
replicate5<-read.table(file='../ViralTransitions/BA1/voc_seed4321_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate5$replicate <- "5"
replicate6<-read.table(file='../ViralTransitions/BA1/voc_seed5432_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate6$replicate <- "6"
replicate7<-read.table(file='../ViralTransitions/BA1/voc_seed6543_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate7$replicate <- "7"
replicate8<-read.table(file='../ViralTransitions/BA1/voc_seed7654_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate8$replicate <- "8"
replicate9<-read.table(file='../ViralTransitions/BA1/voc_seed8765_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate9$replicate <- "9"
replicate10<-read.table(file='../ViralTransitions/BA1/voc_seed9876_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate10$replicate <- "10"


BA1_all_replicates<-rbind(replicate1,replicate2,replicate3,replicate4,
                          replicate5,replicate6,replicate7,replicate8,
                          replicate9,replicate10)

BA1_all_replicates[BA1_all_replicates == "USA"] <- "United States"
BA1_all_replicates[BA1_all_replicates == "Republic of the Congo"] <- "Congo (Brazzaville)"
BA1_all_replicates[BA1_all_replicates == "Democratic Republic of the Congo"] <- "Congo (Kinshasa)"

BA1_all_replicates[BA1_all_replicates == "Coted Ivoire"] <- "Ivory Coast"

BA1_all_replicates<-subset(BA1_all_replicates,Origin!='Reunion')
BA1_all_replicates<-subset(BA1_all_replicates,Origin!='Mayotte')

BA1_all_replicates$date<-date_decimal(BA1_all_replicates$EventTime)

BA1_all_replicates$Origin_Continent<-lapply(BA1_all_replicates$Origin,country2continent)
BA1_all_replicates$Destination_Continent<-lapply(BA1_all_replicates$Destination,country2continent)

BA1_all_replicates$Origin_Continent_Region<-lapply(BA1_all_replicates$Origin,country2continent_region)
BA1_all_replicates$Destination_Continent_Region<-lapply(BA1_all_replicates$Destination,country2continent_region)

BA1_all_replicates$Variant <- "BA1"

BA1_all_replicates$days<-as.Date(cut(BA1_all_replicates$date,breaks = "day",start.on.monday = FALSE))
BA1_all_replicates$date<-as.Date(cut(BA1_all_replicates$date,breaks = "week",start.on.monday = FALSE))
BA1_all_replicates$date2<-as.Date(cut(BA1_all_replicates$date,breaks = "2 week",start.on.monday = FALSE))
BA1_all_replicates$date4<-as.Date(cut(BA1_all_replicates$date,breaks = "1 month",start.on.monday = FALSE))


## BA.2

replicate1<-read.table(file='../ViralTransitions/BA2/voc_seed0987_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate1$replicate <- "1"
replicate2<-read.table(file='../ViralTransitions/BA2/voc_seed1234_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate2$replicate <- "2"
replicate3<-read.table(file='../ViralTransitions/BA2/voc_seed2345_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate3$replicate <- "3"
replicate4<-read.table(file='../ViralTransitions/BA2/voc_seed3456_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate4$replicate <- "4"
replicate5<-read.table(file='../ViralTransitions/BA2/voc_seed4567_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate5$replicate <- "5"
replicate6<-read.table(file='../ViralTransitions/BA2/voc_seed5678_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate6$replicate <- "6"
replicate7<-read.table(file='../ViralTransitions/BA2/voc_seed6789_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate7$replicate <- "7"
replicate8<-read.table(file='../ViralTransitions/BA2/voc_seed7890_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate8$replicate <- "8"
replicate9<-read.table(file='../ViralTransitions/BA2/voc_seed8901_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate9$replicate <- "9"
replicate10<-read.table(file='../ViralTransitions/BA2/voc_seed9012_annottated_tree_events.csv', sep = '\t', header = TRUE)
replicate10$replicate <- "10"

BA2_all_replicates<-rbind(replicate1,replicate2,replicate3,replicate4,
                          replicate5,replicate6,replicate7,replicate8,
                          replicate9,replicate10)
BA2_all_replicates[BA2_all_replicates == "USA"] <- "United States"
BA2_all_replicates[BA2_all_replicates == "Republic of the Congo"] <- "Congo (Brazzaville)"
BA2_all_replicates[BA2_all_replicates == "Democratic Republic of the Congo"] <- "Congo (Kinshasa)"

BA2_all_replicates[BA2_all_replicates == "Coted Ivoire"] <- "Ivory Coast"

BA2_all_replicates<-subset(BA2_all_replicates,Origin!='Reunion')
BA2_all_replicates<-subset(BA2_all_replicates,Origin!='Mayotte')

BA2_all_replicates$date<-date_decimal(BA2_all_replicates$EventTime)

BA2_all_replicates$Origin_Continent<-lapply(BA2_all_replicates$Origin,country2continent)
BA2_all_replicates$Destination_Continent<-lapply(BA2_all_replicates$Destination,country2continent)

BA2_all_replicates$Origin_Continent_Region<-lapply(BA2_all_replicates$Origin,country2continent_region)
BA2_all_replicates$Destination_Continent_Region<-lapply(BA2_all_replicates$Destination,country2continent_region)

BA2_all_replicates$Variant <- "BA2"

BA2_all_replicates$days<-as.Date(cut(BA2_all_replicates$date,breaks = "day",start.on.monday = FALSE))
BA2_all_replicates$date<-as.Date(cut(BA2_all_replicates$date,breaks = "week",start.on.monday = FALSE))
BA2_all_replicates$date2<-as.Date(cut(BA2_all_replicates$date,breaks = "2 week",start.on.monday = FALSE))
BA2_all_replicates$date4<-as.Date(cut(BA2_all_replicates$date,breaks = "1 month",start.on.monday = FALSE))


########## FIGURE 1 ###############


#Figure 1A:
########### Global Dissemination maps

worldmap <- getMap()
world_map_data<-worldmap@data
world_map_data_short<-world_map_data  %>% dplyr::select("NAME","POP_EST","REGION")


#Alpha


alpha_all_replicates$destination_lat<- lapply(alpha_all_replicates$Destination,country2lat)
alpha_all_replicates$destination_long<- lapply(alpha_all_replicates$Destination,country2long)

alpha_all_replicates$origin_lat<- lapply(alpha_all_replicates$Origin,country2lat)
alpha_all_replicates$origin_long<- lapply(alpha_all_replicates$Origin,country2long)


world1 <- map_data("world")

world3 <- map_data("world")

world3<- world3 %>% 
  left_join(alpha_all_replicates, by = c("region" = "Destination"))

alpha_all_replicates$Origin_Continent_Region<-as.character(alpha_all_replicates$Origin_Continent_Region)
alpha_all_replicates$Origin_Continent<-as.character(alpha_all_replicates$Origin_Continent)

alpha_all_replicates<-subset(alpha_all_replicates,Destination_Continent!='character(0)')
alpha_all_replicates$Destination_Continent<-unlist(alpha_all_replicates$Destination_Continent)

alpha_all_replicates<-alpha_all_replicates %>% 
  group_by(Origin_Continent_Region) %>% 
  mutate(mean_origin_lat=mean(unlist(origin_lat)), mean_origin_long=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent_Region) %>% 
  mutate(mean_destination_lat=mean(unlist(destination_lat)), mean_destination_long=mean(unlist(destination_long))) %>% 
  ungroup()
alpha_all_replicates<-subset(subset(subset(alpha_all_replicates, Origin_Continent!="character(0)"),mean_origin_lat!='numeric(0)'),mean_origin_long!='numeric(0)')


alpha_all_replicates<-alpha_all_replicates %>% 
  group_by(Origin_Continent_Region, Destination_Continent_Region) %>% 
  mutate(mean_date=mean(days)) %>% 
  ungroup() 


alpha_all_replicates<-alpha_all_replicates %>% 
  group_by(Origin_Continent) %>% 
  mutate(mean_origin_lat2=mean(unlist(origin_lat)), mean_origin_long2=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent) %>% 
  mutate(mean_destination_lat2=mean(unlist(destination_lat)), mean_destination_long2=mean(unlist(destination_long))) %>% 
  ungroup()

d_alpha_links<-alpha_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,Destination_Continent_Region,mean_origin_lat,mean_origin_long,mean_destination_lat,mean_destination_long, mean_date)  %>%
  count()
d_alpha_links<-subset(subset(subset(d_alpha_links,n>10),Origin_Continent_Region!='character(0)'),Destination_Continent_Region!='character(0)')

d_alpha_origins<-alpha_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,mean_origin_lat,mean_origin_long)  %>%
  count()


d_alpha_origins2<-alpha_all_replicates %>%
  dplyr::group_by(replicate, Origin_Continent,mean_origin_lat2,mean_origin_long2)  %>%
  dplyr::count()

d_alpha_origins3<-alpha_all_replicates %>%
  group_by(replicate,Origin,mean_origin_lat,mean_origin_long)  %>%
  count() %>%
  ungroup() %>%
  group_by(Origin)  %>%
  mutate(n_mean=mean(n))
d_alpha_origins3



d_alpha_links2<-alpha_all_replicates %>%
  group_by(replicate, Origin_Continent,Destination_Continent,mean_origin_lat2,mean_origin_long2,mean_destination_lat2,mean_destination_long2)  %>%
  count()

alpha_origin_continental<-alpha_all_replicates %>%
  group_by(replicate)  %>%
  count(Origin_Continent)
colnames(alpha_origin_continental)<-c("replicate","Continent","Exports")


alpha_destination_continental<-alpha_all_replicates %>%
  group_by(replicate)  %>%
  count(Destination_Continent)
colnames(alpha_destination_continental)<-c("replicate","Continent","Imports")

alpha_source_sink<-left_join(subset(alpha_origin_continental,Continent!='character(0)'),subset(alpha_destination_continental,Continent!='character(0)'))
alpha_source_sink$net_movements<-alpha_source_sink$Exports-alpha_source_sink$Imports
colnames(alpha_source_sink)<-c("replicate","Origin_Continent","Exports","Imports",'net_movements')


alpha_net<-ggplot(alpha_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('red3','dodgerblue3'),name="",labels=c("Sink","Source"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=net_movements,fill=..y..>0))+
  # ggtitle("Alpha")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('Net Exports')+
  theme(axis.text=element_text(size=7))+
  
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
alpha_net


alpha_absolutes<-ggplot(alpha_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('dodgerblue4','red4'),name="",labels=c("Exports","Imports"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=Exports, fill='Exports'))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=-Imports, fill='Imports'))+
  
  # ggtitle("alpha")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('No. Viral Movements')+
  theme(axis.text=element_text(size=7))+
  scale_y_continuous(limits=c(-5000,5000),breaks=seq(-5000,5000,1000),labels=abs(seq(-5000,5000,1000)))+
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
alpha_absolutes


regional_cols<-c("Europe" = "#1b9e77", "Africa" = "#d95f02", "Asia" = "#7570b3", "North America" = "#e7298a",'South America'='gold2','Australia'='dodgerblue3')

alpha_all_replicates_global_regional<-subset(subset(alpha_all_replicates, Origin_Continent!='character(0)'),Destination_Continent!='character(0)')
alpha_all_replicates_global_regional$Export_Type<- mapply(export2type,alpha_all_replicates_global_regional$Origin_Continent,alpha_all_replicates_global_regional$Destination_Continent)


alpha_all_replicates_global_regional_count<-alpha_all_replicates_global_regional %>%
  group_by(Origin,Export_Type)  %>%
  count(Export_Type)  %>%
  ungroup()

colnames(alpha_all_replicates_global_regional_count)<-c("Origin","Export_Type","Counts")

alpha_all_replicates_global_regional_count$Variant<-"Alpha"


alpha_all_replicates_regional_count<-subset(alpha_all_replicates_global_regional_count,Export_Type=='Regional')
total_regional_counts<-sum(alpha_all_replicates_regional_count$Counts)
alpha_all_replicates_regional_count$ProportionsOfRegional<-alpha_all_replicates_regional_count$Counts/total_regional_counts
alpha_all_replicates_regional_count<-alpha_all_replicates_regional_count%>%select("Origin","ProportionsOfRegional")


alpha_all_replicates_global_count<-subset(alpha_all_replicates_global_regional_count,Export_Type=='Global')
total_global_counts<-sum(alpha_all_replicates_global_count$Counts)
alpha_all_replicates_global_count$ProportionsOfGlobal<-alpha_all_replicates_global_count$Counts/total_global_counts
alpha_all_replicates_global_count<-alpha_all_replicates_global_count%>%select("Origin","ProportionsOfGlobal")

alpha_all_replicates_regional_global_count<-left_join(alpha_all_replicates_regional_count,alpha_all_replicates_global_count)
alpha_all_replicates_regional_global_count$Variant<-"Alpha"


alpha_destination_continental<-alpha_all_replicates %>%
  group_by(replicate)  %>%
  count(Destination_Continent)
colnames(alpha_destination_continental)<-c("replicate","Continent","Imports")

d_alpha_links$mean_date<-as.Date(d_alpha_links$mean_date)
lab_dates <- pretty(d_alpha_links$mean_date)

alpha_map<-ggplot() +
  theme_void()+
  geom_map(data=world1,map=world1, aes(long, lat,map_id=region), color="white", fill='snow3',size=0.2)+
  
  geom_curve(data = d_alpha_links,
             aes(x = as.double(mean_origin_long), 
                 y = as.double(mean_origin_lat), 
                 xend = as.double(mean_destination_long), 
                 yend = as.double(mean_destination_lat), colour=as.Date(mean_date)),size=2,
             alpha=0.8)+
  
  geom_point(data = d_alpha_origins,
             aes(x = as.double(mean_origin_long), y = as.double(mean_origin_lat), size=n),
             color='black', shape=21)+
  scale_color_gradientn(colours = c('#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221'),breaks = lab_dates, 
                        labels = lab_dates, name='Inferred Dispersal Dates (Mean)')+
  scale_fill_manual(values=regional_cols, name='Regions of Origin')+
  theme(legend.position = 'top')+
  theme(legend.direction = 'horizontal')+
   coord_fixed()+
  guides(colour = guide_colourbar(barwidth = 20, barheight = 0.5,title.position = 'top',title.hjust = 1,ticks.colour = "white",
                                  ticks.linewidth = 2),
         size='none')+
  theme(plot.title = element_text(family="Helvetica"))+
  scale_y_continuous(limits=c(-55,90))+
  scale_x_continuous(limits=c(-170,170))

#alpha_map

alpha_fig2<-alpha_map +
  inset_element(alpha_net, 0.0, 0.0, 0.3, 0.5)
alpha_fig2

#Beta

beta_all_replicates$destination_lat<- lapply(beta_all_replicates$Destination,country2lat)
beta_all_replicates$destination_long<- lapply(beta_all_replicates$Destination,country2long)

beta_all_replicates$origin_lat<- lapply(beta_all_replicates$Origin,country2lat)
beta_all_replicates$origin_long<- lapply(beta_all_replicates$Origin,country2long)


world3 <- map_data("world")

world3<- world3 %>% 
  left_join(beta_all_replicates, by = c("region" = "Destination"))

beta_all_replicates$Origin_Continent_Region<-as.character(beta_all_replicates$Origin_Continent_Region)
beta_all_replicates$Origin_Continent<-as.character(beta_all_replicates$Origin_Continent)


beta_all_replicates<-subset(beta_all_replicates,Destination_Continent!='character(0)')
beta_all_replicates$Destination_Continent<-unlist(beta_all_replicates$Destination_Continent)


beta_all_replicates<-beta_all_replicates %>% 
  group_by(Origin_Continent_Region) %>% 
  mutate(mean_origin_lat=mean(unlist(origin_lat)), mean_origin_long=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent_Region) %>% 
  mutate(mean_destination_lat=mean(unlist(destination_lat)), mean_destination_long=mean(unlist(destination_long))) %>% 
  ungroup()


beta_all_replicates<-beta_all_replicates %>% 
  group_by(Origin_Continent) %>% 
  mutate(mean_origin_lat2=mean(unlist(origin_lat)), mean_origin_long2=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent) %>% 
  mutate(mean_destination_lat2=mean(unlist(destination_lat)), mean_destination_long2=mean(unlist(destination_long))) %>% 
  ungroup()

beta_all_replicates<-beta_all_replicates %>% 
  group_by(Origin_Continent_Region, Destination_Continent_Region) %>% 
  mutate(mean_date=mean(days)) %>% 
  ungroup() 


d_beta_links<-beta_all_replicates %>%
  dplyr::group_by(Origin_Continent,Origin_Continent_Region,Destination_Continent_Region,mean_origin_lat,mean_origin_long,mean_destination_lat,mean_destination_long,mean_date)  %>%
  dplyr::count()
d_beta_links<-subset(subset(subset(d_beta_links,n>10),Origin_Continent_Region!='character(0)'),Destination_Continent_Region!='character(0)')


d_beta_origins<-beta_all_replicates %>%
  dplyr::group_by(Origin_Continent,Origin_Continent_Region,mean_origin_lat,mean_origin_long)  %>%
  dplyr::count()



d_beta_origins3<-beta_all_replicates %>%
  group_by(replicate,Origin,mean_origin_lat,mean_origin_long)  %>%
  count() %>%
  ungroup() %>%
  group_by(Origin)  %>%
  mutate(n_mean=mean(n))
d_beta_origins3


d_beta_links2<-beta_all_replicates %>%
  dplyr::group_by(replicate, Origin_Continent,Destination_Continent,mean_origin_lat2,mean_origin_long2,mean_destination_lat2,mean_destination_long2)  %>%
  dplyr::count()


beta_origin_continental<-beta_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Origin_Continent)
colnames(beta_origin_continental)<-c("replicate","Continent","Exports")


beta_destination_continental<-beta_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Destination_Continent)
colnames(beta_destination_continental)<-c("replicate","Continent","Imports")

beta_source_sink<-left_join(subset(beta_origin_continental,Continent!='character(0)'),subset(beta_destination_continental,Continent!='character(0)'))
beta_source_sink$net_movements<-beta_source_sink$Exports-beta_source_sink$Imports
colnames(beta_source_sink)<-c("replicate","Origin_Continent","Exports","Imports",'net_movements')


beta_net<-ggplot(beta_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('red3','dodgerblue3'),name="",labels=c("Sink","Source"))+
  geom_bar(stat='summary',fun.y='mean',aes(x=Origin_Continent,y=net_movements,fill=..y..>0))+
  #ggtitle("Beta")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('Net Exports')+
  theme(axis.text=element_text(size=7))+
  
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
beta_net
plot_grid(alpha_net,beta_net)


beta_absolutes<-ggplot(beta_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('dodgerblue4','red4'),name="",labels=c("Exports","Imports"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=Exports, fill='Exports'))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=-Imports, fill='Imports'))+
  
  # ggtitle("beta")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('No. Viral Movements')+
  theme(axis.text=element_text(size=7))+
  scale_y_continuous(limits=c(-5000,5000),breaks=seq(-5000,5000,1000),labels=abs(seq(-5000,5000,1000)))+
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
beta_absolutes


beta_all_replicates_global_regional<-subset(subset(beta_all_replicates, Origin_Continent!='character(0)'),Destination_Continent!='character(0)')
beta_all_replicates_global_regional$Export_Type<- mapply(export2type,beta_all_replicates_global_regional$Origin_Continent,beta_all_replicates_global_regional$Destination_Continent)

#prop.table(table(data2$days, data2$Nextstrain_variants))


beta_all_replicates_global_regional_count<-beta_all_replicates_global_regional %>%
  group_by(Origin,Export_Type)  %>%
  count(Export_Type)  %>%
  ungroup()

colnames(beta_all_replicates_global_regional_count)<-c("Origin","Export_Type","Counts")
beta_all_replicates_global_regional_count$Variant<-"Beta"



beta_all_replicates_regional_count<-subset(beta_all_replicates_global_regional_count,Export_Type=='Regional')
total_regional_counts<-sum(beta_all_replicates_regional_count$Counts)
beta_all_replicates_regional_count$ProportionsOfRegional<-beta_all_replicates_regional_count$Counts/total_regional_counts
beta_all_replicates_regional_count<-beta_all_replicates_regional_count%>%select("Origin","ProportionsOfRegional")


beta_all_replicates_global_count<-subset(beta_all_replicates_global_regional_count,Export_Type=='Global')
total_global_counts<-sum(beta_all_replicates_global_count$Counts)
beta_all_replicates_global_count$ProportionsOfGlobal<-beta_all_replicates_global_count$Counts/total_global_counts
beta_all_replicates_global_count<-beta_all_replicates_global_count%>%select("Origin","ProportionsOfGlobal")

beta_all_replicates_regional_global_count<-left_join(beta_all_replicates_regional_count,beta_all_replicates_global_count)
beta_all_replicates_regional_global_count$Variant<-"Beta"


regional_cols<-c("Europe" = "#1b9e77", "Africa" = "#d95f02", "Asia" = "#7570b3", "North America" = "#e7298a",'South America'='gold2','Australia'='dodgerblue3')


d_beta_links$mean_date<-as.Date(d_beta_links$mean_date)
lab_dates <- pretty(d_beta_links$mean_date)

beta_map<-ggplot() +
  theme_void()+
  geom_map(data=world1,map=world1, aes(long, lat,map_id=region), color="white", fill='snow3',size=0.2)+
  
  geom_curve(data = d_beta_links,
             aes(x = as.double(mean_origin_long), 
                 y = as.double(mean_origin_lat), 
                 xend = as.double(mean_destination_long), 
                 yend = as.double(mean_destination_lat), colour=as.Date(mean_date)),size=2,
             alpha=0.8)+
  
  geom_point(data = d_beta_origins,
             aes(x = as.double(mean_origin_long), y = as.double(mean_origin_lat), size=n),
             color='black', shape=21)+
  scale_color_gradientn(colours = c('#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221'),breaks = lab_dates, 
                        labels = lab_dates, name='Inferred Dispersal Dates (Mean)')+
  scale_fill_manual(values=regional_cols, name='Regions of Origin')+
  theme(legend.position = 'top')+
  theme(legend.direction = 'horizontal')+
  coord_fixed()+
  guides(colour = guide_colourbar(barwidth = 20, barheight = 0.5,title.position = 'top',title.hjust = 1,ticks.colour = "white",
                                  ticks.linewidth = 2),
         size='none')+
  theme(plot.title = element_text(family="Helvetica"))+
  scale_y_continuous(limits=c(-55,90))+
  scale_x_continuous(limits=c(-170,170))


beta_map


beta_fig2<-beta_map +
  inset_element(beta_net, 0.0, 0.0, 0.3, 0.5)
beta_fig2

#Gamma


gamma_all_replicates$destination_lat<- lapply(gamma_all_replicates$Destination,country2lat)
gamma_all_replicates$destination_long<- lapply(gamma_all_replicates$Destination,country2long)

gamma_all_replicates$origin_lat<- lapply(gamma_all_replicates$Origin,country2lat)
gamma_all_replicates$origin_long<- lapply(gamma_all_replicates$Origin,country2long)

world3 <- map_data("world")

world3<- world3 %>% 
  left_join(gamma_all_replicates, by = c("region" = "Destination"))

gamma_all_replicates$Origin_Continent_Region<-as.character(gamma_all_replicates$Origin_Continent_Region)
gamma_all_replicates$Origin_Continent<-as.character(gamma_all_replicates$Origin_Continent)



gamma_all_replicates<-subset(gamma_all_replicates,Destination_Continent!='character(0)')
gamma_all_replicates$Destination_Continent<-unlist(gamma_all_replicates$Destination_Continent)

gamma_all_replicates<-gamma_all_replicates %>% 
  group_by(Origin_Continent_Region) %>% 
  mutate(mean_origin_lat=mean(unlist(origin_lat)), mean_origin_long=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent_Region) %>% 
  mutate(mean_destination_lat=mean(unlist(destination_lat)), mean_destination_long=mean(unlist(destination_long))) %>% 
  ungroup()


gamma_all_replicates<-gamma_all_replicates %>% 
  group_by(Origin_Continent) %>% 
  mutate(mean_origin_lat2=mean(unlist(origin_lat)), mean_origin_long2=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent) %>% 
  mutate(mean_destination_lat2=mean(unlist(destination_lat)), mean_destination_long2=mean(unlist(destination_long))) %>% 
  ungroup()


gamma_all_replicates<-gamma_all_replicates %>% 
  group_by(Origin_Continent_Region, Destination_Continent_Region) %>% 
  mutate(mean_date=mean(days)) %>% 
  ungroup() 


d_gamma_links<-gamma_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,Destination_Continent_Region,mean_origin_lat,mean_origin_long,mean_destination_lat,mean_destination_long,mean_date)  %>%
  count()
d_gamma_links<-subset(subset(subset(d_gamma_links,n>10),Origin_Continent_Region!='character(0)'),Destination_Continent_Region!='character(0)')


d_gamma_origins<-gamma_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,mean_origin_lat,mean_origin_long)  %>%
  count()



d_gamma_origins3<-gamma_all_replicates %>%
  group_by(replicate,Origin,mean_origin_lat,mean_origin_long)  %>%
  count() %>%
  ungroup() %>%
  group_by(Origin)  %>%
  mutate(n_mean=mean(n))
d_gamma_origins3


d_gamma_links2<-gamma_all_replicates %>%
  group_by(replicate, Origin_Continent,Destination_Continent,mean_origin_lat2,mean_origin_long2,mean_destination_lat2,mean_destination_long2)  %>%
  count()


gamma_origin_continental<-gamma_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Origin_Continent)
colnames(gamma_origin_continental)<-c("replicate","Continent","Exports")


gamma_destination_continental<-gamma_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Destination_Continent)
colnames(gamma_destination_continental)<-c("replicate","Continent","Imports")

gamma_source_sink<-left_join(subset(gamma_origin_continental,Continent!='character(0)'),subset(gamma_destination_continental,Continent!='character(0)'))
gamma_source_sink$net_movements<-gamma_source_sink$Exports-gamma_source_sink$Imports
colnames(gamma_source_sink)<-c("replicate","Origin_Continent","Exports","Imports",'net_movements')



gamma_net<-ggplot(gamma_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('red3','dodgerblue3'),name="",labels=c("Sink","Source"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=net_movements,fill=..y..>0))+
  # ggtitle("Gamma")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('Net Exports')+
  theme(axis.text=element_text(size=7))+
  
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
gamma_net


gamma_absolutes<-ggplot(gamma_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('dodgerblue4','red4'),name="",labels=c("Exports","Imports"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=Exports, fill='Exports'))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=-Imports, fill='Imports'))+
  
  # ggtitle("Gamma")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('No. Viral Movements')+
  theme(axis.text=element_text(size=7))+
  scale_y_continuous(limits=c(-5000,5000),breaks=seq(-5000,5000,1000),labels=abs(seq(-5000,5000,1000)))+
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
gamma_absolutes


gamma_all_replicates_global_regional<-subset(subset(gamma_all_replicates, Origin_Continent!='character(0)'),Destination_Continent!='character(0)')
gamma_all_replicates_global_regional$Export_Type<- mapply(export2type,gamma_all_replicates_global_regional$Origin_Continent,gamma_all_replicates_global_regional$Destination_Continent)

gamma_all_replicates_global_regional_count<-gamma_all_replicates_global_regional %>%
  group_by(Origin,Export_Type)  %>%
  count(Export_Type)  %>%
  ungroup()

colnames(gamma_all_replicates_global_regional_count)<-c("Origin","Export_Type","Counts")
gamma_all_replicates_global_regional_count$Variant<-"Gamma"


gamma_all_replicates_regional_count<-subset(gamma_all_replicates_global_regional_count,Export_Type=='Regional')
total_regional_counts<-sum(gamma_all_replicates_regional_count$Counts)
gamma_all_replicates_regional_count$ProportionsOfRegional<-gamma_all_replicates_regional_count$Counts/total_regional_counts
gamma_all_replicates_regional_count<-gamma_all_replicates_regional_count%>%select("Origin","ProportionsOfRegional")


gamma_all_replicates_global_count<-subset(gamma_all_replicates_global_regional_count,Export_Type=='Global')
total_global_counts<-sum(gamma_all_replicates_global_count$Counts)
gamma_all_replicates_global_count$ProportionsOfGlobal<-gamma_all_replicates_global_count$Counts/total_global_counts
gamma_all_replicates_global_count<-gamma_all_replicates_global_count%>%select("Origin","ProportionsOfGlobal")

gamma_all_replicates_regional_global_count<-left_join(gamma_all_replicates_regional_count,gamma_all_replicates_global_count)
gamma_all_replicates_regional_global_count$Variant<-"Gamma"


regional_cols<-c("Europe" = "#1b9e77", "Africa" = "#d95f02", "Asia" = "#7570b3", "North America" = "#e7298a",'South America'='gold2','Australia'='dodgerblue3')

d_gamma_links$mean_date<-as.Date(d_gamma_links$mean_date)
lab_dates <- pretty(d_gamma_links$mean_date)

gamma_map<-ggplot() +
  theme_void()+
  geom_map(data=world1,map=world1, aes(long, lat,map_id=region), color="white", fill='snow3',size=0.2)+
  
  geom_curve(data = d_gamma_links,
             aes(x = as.double(mean_origin_long), 
                 y = as.double(mean_origin_lat), 
                 xend = as.double(mean_destination_long), 
                 yend = as.double(mean_destination_lat), colour=as.Date(mean_date)),size=2,
             alpha=0.8)+
  
  geom_point(data = d_gamma_origins,
             aes(x = as.double(mean_origin_long), y = as.double(mean_origin_lat), size=n),
             color='black', shape=21)+
  scale_color_gradientn(colours = c('#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221'),breaks = lab_dates, 
                        labels = lab_dates, name='Inferred Dispersal Dates (Mean)')+
  scale_fill_manual(values=regional_cols, name='Regions of Origin')+
  theme(legend.position = 'top')+
  theme(legend.direction = 'horizontal')+
  coord_fixed()+
  guides(colour = guide_colourbar(barwidth = 20, barheight = 0.5,title.position = 'top',title.hjust = 1,ticks.colour = "white",
                                  ticks.linewidth = 2),
         size='none')+
  theme(plot.title = element_text(family="Helvetica"))+
  scale_y_continuous(limits=c(-55,90))+
  scale_x_continuous(limits=c(-170,170))
gamma_map

gamma_fig2<-gamma_map +
  inset_element(gamma_net, 0.0, 0.0, 0.3, 0.5)
gamma_fig2


#Delta


delta_all_replicates$destination_lat<- lapply(delta_all_replicates$Destination,country2lat)
delta_all_replicates$destination_long<- lapply(delta_all_replicates$Destination,country2long)

delta_all_replicates$origin_lat<- lapply(delta_all_replicates$Origin,country2lat)
delta_all_replicates$origin_long<- lapply(delta_all_replicates$Origin,country2long)

world3 <- map_data("world")

world3<- world3 %>% 
  left_join(delta_all_replicates, by = c("region" = "Destination"))

delta_all_replicates$Origin_Continent_Region<-as.character(delta_all_replicates$Origin_Continent_Region)
delta_all_replicates$Origin_Continent<-as.character(delta_all_replicates$Origin_Continent)


delta_all_replicates<-subset(delta_all_replicates,Destination_Continent!='character(0)')
delta_all_replicates$Destination_Continent<-unlist(delta_all_replicates$Destination_Continent)

delta_all_replicates<-delta_all_replicates %>% 
  group_by(Origin_Continent_Region) %>% 
  mutate(mean_origin_lat=mean(unlist(origin_lat)), mean_origin_long=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent_Region) %>% 
  mutate(mean_destination_lat=mean(unlist(destination_lat)), mean_destination_long=mean(unlist(destination_long))) %>% 
  ungroup()


delta_all_replicates<-delta_all_replicates %>% 
  group_by(Origin_Continent) %>% 
  mutate(mean_origin_lat2=mean(unlist(origin_lat)), mean_origin_long2=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent) %>% 
  mutate(mean_destination_lat2=mean(unlist(destination_lat)), mean_destination_long2=mean(unlist(destination_long))) %>% 
  ungroup()

delta_all_replicates<-delta_all_replicates %>% 
  group_by(Origin_Continent_Region, Destination_Continent_Region) %>% 
  mutate(mean_date=mean(days)) %>% 
  ungroup() 


d_delta_links<-delta_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,Destination_Continent_Region,mean_origin_lat,mean_origin_long,mean_destination_lat,mean_destination_long,mean_date)  %>%
  count()
d_delta_links<-subset(subset(subset(d_delta_links,n>10),Origin_Continent_Region!='character(0)'),Destination_Continent_Region!='character(0)')


d_delta_origins<-delta_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,mean_origin_lat,mean_origin_long)  %>%
  count()


d_delta_origins3<-delta_all_replicates %>%
  group_by(replicate,Origin)  %>%
  count() %>%
  ungroup() %>%
  group_by(Origin)  %>%
  mutate(n_mean=mean(n))
d_delta_origins3



d_delta_links2<-delta_all_replicates %>%
  group_by(replicate, Origin_Continent,Destination_Continent,mean_origin_lat2,mean_origin_long2,mean_destination_lat2,mean_destination_long2)  %>%
  count()



delta_origin_continental<-delta_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Origin_Continent)
colnames(delta_origin_continental)<-c("replicate","Continent","Exports")


delta_destination_continental<-delta_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Destination_Continent)
colnames(delta_destination_continental)<-c("replicate","Continent","Imports")

delta_source_sink<-left_join(subset(delta_origin_continental,Continent!='character(0)'),subset(delta_destination_continental,Continent!='character(0)'))
delta_source_sink$net_movements<-delta_source_sink$Exports-delta_source_sink$Imports
colnames(delta_source_sink)<-c("replicate","Origin_Continent","Exports","Imports",'net_movements')



delta_net<-ggplot(delta_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('red3','dodgerblue3'),name="",labels=c("Sink","Source"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=net_movements,fill=..y..>0))+
  #ggtitle("Delta")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('Net Exports')+
  theme(axis.text=element_text(size=7))+
  
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
delta_net


delta_absolutes<-ggplot(delta_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('dodgerblue4','red4'),name="",labels=c("Exports","Imports"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=Exports, fill='Exports'))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=-Imports, fill='Imports'))+
  
  # ggtitle("delta")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('No. Viral Movements')+
  theme(axis.text=element_text(size=7))+
  scale_y_continuous(limits=c(-5000,5000),breaks=seq(-5000,5000,1000),labels=abs(seq(-5000,5000,1000)))+
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
delta_absolutes


delta_all_replicates_global_regional<-subset(subset(delta_all_replicates, Origin_Continent!='character(0)'),Destination_Continent!='character(0)')
delta_all_replicates_global_regional$Export_Type<- mapply(export2type,delta_all_replicates_global_regional$Origin_Continent,delta_all_replicates_global_regional$Destination_Continent)

delta_all_replicates_global_regional_count<-delta_all_replicates_global_regional %>%
  group_by(Origin,Export_Type)  %>%
  count(Export_Type)  %>%
  ungroup()

colnames(delta_all_replicates_global_regional_count)<-c("Origin","Export_Type","Counts")
delta_all_replicates_global_regional_count$Variant<-"Delta"


delta_all_replicates_regional_count<-subset(delta_all_replicates_global_regional_count,Export_Type=='Regional')
total_regional_counts<-sum(delta_all_replicates_regional_count$Counts)
delta_all_replicates_regional_count$ProportionsOfRegional<-delta_all_replicates_regional_count$Counts/total_regional_counts
delta_all_replicates_regional_count<-delta_all_replicates_regional_count%>%select("Origin","ProportionsOfRegional")


delta_all_replicates_global_count<-subset(delta_all_replicates_global_regional_count,Export_Type=='Global')
total_global_counts<-sum(delta_all_replicates_global_count$Counts)
delta_all_replicates_global_count$ProportionsOfGlobal<-delta_all_replicates_global_count$Counts/total_global_counts
delta_all_replicates_global_count<-delta_all_replicates_global_count%>%select("Origin","ProportionsOfGlobal")

delta_all_replicates_regional_global_count<-left_join(delta_all_replicates_regional_count,delta_all_replicates_global_count)
delta_all_replicates_regional_global_count$Variant<-"Delta"

regional_cols<-c("Europe" = "#1b9e77", "Africa" = "#d95f02", "Asia" = "#7570b3", "North America" = "#e7298a",'South America'='gold2','Australia'='dodgerblue3')

d_delta_links$mean_date<-as.Date(d_delta_links$mean_date)
lab_dates <- pretty(d_delta_links$mean_date)

delta_map<-ggplot() +
  theme_void()+
 geom_map(data=world1,map=world1, aes(long, lat,map_id=region), color="white", fill='snow3',size=0.2)+
  
  geom_curve(data = d_delta_links,
             aes(x = as.double(mean_origin_long), 
                 y = as.double(mean_origin_lat), 
                 xend = as.double(mean_destination_long), 
                 yend = as.double(mean_destination_lat), colour=as.Date(mean_date)),size=2,
             alpha=0.8)+
  
  geom_point(data = d_delta_origins,
             aes(x = as.double(mean_origin_long), y = as.double(mean_origin_lat), size=n),
             color='black', shape=21)+
  scale_color_gradientn(colours = c('#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221'),breaks = lab_dates, 
                        labels = lab_dates, name='Inferred Dispersal Dates (Mean)')+
  scale_fill_manual(values=regional_cols, name='Regions of Origin')+
  theme(legend.position = 'top')+
  theme(legend.direction = 'horizontal')+
  coord_fixed()+
  guides(colour = guide_colourbar(barwidth = 25, barheight = 0.5,title.position = 'top',title.hjust = 1,ticks.colour = "white",
                                  ticks.linewidth = 2),
         size='none')+
  theme(plot.title = element_text(family="Helvetica"))+
  scale_y_continuous(limits=c(-55,90))+
  scale_x_continuous(limits=c(-170,170))


delta_map

delta_fig2<-delta_map +
  inset_element(delta_net, 0.0, 0.0, 0.3, 0.5)
delta_fig2



#BA1


BA1_all_replicates$destination_lat<- lapply(BA1_all_replicates$Destination,country2lat)
BA1_all_replicates$destination_long<- lapply(BA1_all_replicates$Destination,country2long)

BA1_all_replicates$origin_lat<- lapply(BA1_all_replicates$Origin,country2lat)
BA1_all_replicates$origin_long<- lapply(BA1_all_replicates$Origin,country2long)

world3 <- map_data("world")

world3<- world3 %>% 
  left_join(BA1_all_replicates, by = c("region" = "Destination"))

BA1_all_replicates$Origin_Continent_Region<-as.character(BA1_all_replicates$Origin_Continent_Region)
BA1_all_replicates$Origin_Continent<-as.character(BA1_all_replicates$Origin_Continent)


BA1_all_replicates<-subset(BA1_all_replicates,Destination_Continent!='character(0)')
BA1_all_replicates$Destination_Continent<-unlist(BA1_all_replicates$Destination_Continent)

BA1_all_replicates<-BA1_all_replicates %>% 
  group_by(Origin_Continent_Region) %>% 
  mutate(mean_origin_lat=mean(unlist(origin_lat)), mean_origin_long=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent_Region) %>% 
  mutate(mean_destination_lat=mean(unlist(destination_lat)), mean_destination_long=mean(unlist(destination_long))) %>% 
  ungroup()


BA1_all_replicates<-BA1_all_replicates %>% 
  group_by(Origin_Continent) %>% 
  mutate(mean_origin_lat2=mean(unlist(origin_lat)), mean_origin_long2=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent) %>% 
  mutate(mean_destination_lat2=mean(unlist(destination_lat)), mean_destination_long2=mean(unlist(destination_long))) %>% 
  ungroup()

BA1_all_replicates<-BA1_all_replicates %>% 
  group_by(Origin_Continent_Region, Destination_Continent_Region) %>% 
  mutate(mean_date=mean(days)) %>% 
  ungroup() 


d_BA1_links<-BA1_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,Destination_Continent_Region,mean_origin_lat,mean_origin_long,mean_destination_lat,mean_destination_long,mean_date)  %>%
  count()
d_BA1_links<-subset(subset(subset(d_BA1_links,n>10),Origin_Continent_Region!='character(0)'),Destination_Continent_Region!='character(0)')


d_BA1_origins<-BA1_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,mean_origin_lat,mean_origin_long)  %>%
  count()



d_BA1_origins3<-BA1_all_replicates %>%
  group_by(replicate,Origin)  %>%
  count() %>%
  ungroup() %>%
  group_by(Origin)  %>%
  mutate(n_mean=mean(n))
d_BA1_origins3

d_BA1_links2<-BA1_all_replicates %>%
  group_by(replicate, Origin_Continent,Destination_Continent,mean_origin_lat2,mean_origin_long2,mean_destination_lat2,mean_destination_long2)  %>%
  count()



BA1_origin_continental<-BA1_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Origin_Continent)
colnames(BA1_origin_continental)<-c("replicate","Continent","Exports")


BA1_destination_continental<-BA1_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Destination_Continent)
colnames(BA1_destination_continental)<-c("replicate","Continent","Imports")

BA1_source_sink<-left_join(subset(BA1_origin_continental,Continent!='character(0)'),subset(BA1_destination_continental,Continent!='character(0)'))
BA1_source_sink$net_movements<-BA1_source_sink$Exports-BA1_source_sink$Imports
colnames(BA1_source_sink)<-c("replicate","Origin_Continent","Exports","Imports",'net_movements')



BA1_net<-ggplot(BA1_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('red3','dodgerblue3'),name="",labels=c("Sink","Source"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=net_movements,fill=..y..>0))+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('Net Exports')+
  theme(axis.text=element_text(size=7))+
  
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
BA1_net


BA1_absolutes<-ggplot(BA1_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('dodgerblue4','red4'),name="",labels=c("Exports","Imports"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=Exports, fill='Exports'))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=-Imports, fill='Imports'))+
    theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('No. Viral Movements')+
  theme(axis.text=element_text(size=7))+
  scale_y_continuous(limits=c(-5000,5000),breaks=seq(-5000,5000,1000),labels=abs(seq(-5000,5000,1000)))+
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
BA1_absolutes


BA1_all_replicates_global_regional<-subset(subset(BA1_all_replicates, Origin_Continent!='character(0)'),Destination_Continent!='character(0)')
BA1_all_replicates_global_regional$Export_Type<- mapply(export2type,BA1_all_replicates_global_regional$Origin_Continent,BA1_all_replicates_global_regional$Destination_Continent)

BA1_all_replicates_global_regional_count<-BA1_all_replicates_global_regional %>%
  group_by(Origin,Export_Type)  %>%
  count(Export_Type)  %>%
  ungroup()

colnames(BA1_all_replicates_global_regional_count)<-c("Origin","Export_Type","Counts")
BA1_all_replicates_global_regional_count$Variant<-"Omicron BA1"

BA1_all_replicates_regional_count<-subset(BA1_all_replicates_global_regional_count,Export_Type=='Regional')
total_regional_counts<-sum(BA1_all_replicates_regional_count$Counts)
BA1_all_replicates_regional_count$ProportionsOfRegional<-BA1_all_replicates_regional_count$Counts/total_regional_counts
BA1_all_replicates_regional_count<-BA1_all_replicates_regional_count%>%select("Origin","ProportionsOfRegional")


BA1_all_replicates_global_count<-subset(BA1_all_replicates_global_regional_count,Export_Type=='Global')
total_global_counts<-sum(BA1_all_replicates_global_count$Counts)
BA1_all_replicates_global_count$ProportionsOfGlobal<-BA1_all_replicates_global_count$Counts/total_global_counts
BA1_all_replicates_global_count<-BA1_all_replicates_global_count%>%select("Origin","ProportionsOfGlobal")

BA1_all_replicates_regional_global_count<-left_join(BA1_all_replicates_regional_count,BA1_all_replicates_global_count)
BA1_all_replicates_regional_global_count$Variant<-"Omicron BA1"


regional_cols<-c("Europe" = "#1b9e77", "Africa" = "#d95f02", "Asia" = "#7570b3", "North America" = "#e7298a",'South America'='gold2','Australia'='dodgerblue3')

d_BA1_links$mean_date<-as.Date(d_BA1_links$mean_date)
lab_dates <- pretty(d_BA1_links$mean_date)

BA1_map<-ggplot() +
  theme_void()+
  geom_map(data=world1,map=world1, aes(long, lat,map_id=region), color="white", fill='snow3',size=0.2)+
  
  geom_curve(data = d_BA1_links,
             aes(x = as.double(mean_origin_long), 
                 y = as.double(mean_origin_lat), 
                 xend = as.double(mean_destination_long), 
                 yend = as.double(mean_destination_lat), colour=as.Date(mean_date)),size=2,
             alpha=0.8)+
  
  geom_point(data = d_BA1_origins,
             aes(x = as.double(mean_origin_long), y = as.double(mean_origin_lat), size=n),
             color='black', shape=21)+
  scale_color_gradientn(colours = c('#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221'),breaks = lab_dates, 
                        labels = lab_dates, name='Inferred Dispersal Dates (Mean)')+
  scale_fill_manual(values=regional_cols, name='Regions of Origin')+
  theme(legend.position = 'top')+
  theme(legend.direction = 'horizontal')+
  coord_fixed()+
  #ggtitle("Alpha")+
  guides(colour = guide_colourbar(barwidth = 20, barheight = 0.5,title.position = 'top',title.hjust = 1,ticks.colour = "white",
                                  ticks.linewidth = 2),
         size='none')+
  #  scale_size(guide = guide_legend(direction = "vertical",position=c(0.2,0.5)))+
  theme(plot.title = element_text(family="Helvetica"))+
  scale_y_continuous(limits=c(-55,90))+
  scale_x_continuous(limits=c(-170,170))


BA1_map

BA1_fig2<-BA1_map +
  inset_element(BA1_net, 0.0, 0.0, 0.3, 0.5)
BA1_fig2


#BA2


BA2_all_replicates$destination_lat<- lapply(BA2_all_replicates$Destination,country2lat)
BA2_all_replicates$destination_long<- lapply(BA2_all_replicates$Destination,country2long)

BA2_all_replicates$origin_lat<- lapply(BA2_all_replicates$Origin,country2lat)
BA2_all_replicates$origin_long<- lapply(BA2_all_replicates$Origin,country2long)

world3 <- map_data("world")

world3<- world3 %>% 
  left_join(BA2_all_replicates, by = c("region" = "Destination"))

BA2_all_replicates$Origin_Continent_Region<-as.character(BA2_all_replicates$Origin_Continent_Region)
BA2_all_replicates$Origin_Continent<-as.character(BA2_all_replicates$Origin_Continent)


BA2_all_replicates<-subset(BA2_all_replicates,Destination_Continent!='character(0)')
BA2_all_replicates$Destination_Continent<-unlist(BA2_all_replicates$Destination_Continent)

BA2_all_replicates<-BA2_all_replicates %>% 
  group_by(Origin_Continent_Region) %>% 
  mutate(mean_origin_lat=mean(unlist(origin_lat)), mean_origin_long=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent_Region) %>% 
  mutate(mean_destination_lat=mean(unlist(destination_lat)), mean_destination_long=mean(unlist(destination_long))) %>% 
  ungroup()


BA2_all_replicates<-BA2_all_replicates %>% 
  group_by(Origin_Continent) %>% 
  mutate(mean_origin_lat2=mean(unlist(origin_lat)), mean_origin_long2=mean(unlist(origin_long))) %>% 
  ungroup() %>% 
  group_by(Destination_Continent) %>% 
  mutate(mean_destination_lat2=mean(unlist(destination_lat)), mean_destination_long2=mean(unlist(destination_long))) %>% 
  ungroup()

BA2_all_replicates<-BA2_all_replicates %>% 
  group_by(Origin_Continent_Region, Destination_Continent_Region) %>% 
  mutate(mean_date=mean(days)) %>% 
  ungroup() 


d_BA2_links<-BA2_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,Destination_Continent_Region,mean_origin_lat,mean_origin_long,mean_destination_lat,mean_destination_long,mean_date)  %>%
  count()
d_BA2_links<-subset(subset(subset(d_BA2_links,n>10),Origin_Continent_Region!='character(0)'),Destination_Continent_Region!='character(0)')


d_BA2_origins<-BA2_all_replicates %>%
  group_by(Origin_Continent,Origin_Continent_Region,mean_origin_lat,mean_origin_long)  %>%
  count()


d_BA2_origins3<-BA2_all_replicates %>%
  group_by(replicate,Origin)  %>%
  count() %>%
  ungroup() %>%
  group_by(Origin)  %>%
  mutate(n_mean=mean(n))
d_BA2_origins3

d_BA2_links2<-BA2_all_replicates %>%
  group_by(replicate, Origin_Continent,Destination_Continent,mean_origin_lat2,mean_origin_long2,mean_destination_lat2,mean_destination_long2)  %>%
  count()



BA2_origin_continental<-BA2_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Origin_Continent)
colnames(BA2_origin_continental)<-c("replicate","Continent","Exports")


BA2_destination_continental<-BA2_all_replicates %>%
  dplyr::group_by(replicate)  %>%
  dplyr::count(Destination_Continent)
colnames(BA2_destination_continental)<-c("replicate","Continent","Imports")

BA2_source_sink<-left_join(subset(BA2_origin_continental,Continent!='character(0)'),subset(BA2_destination_continental,Continent!='character(0)'))
BA2_source_sink$net_movements<-BA2_source_sink$Exports-BA2_source_sink$Imports
colnames(BA2_source_sink)<-c("replicate","Origin_Continent","Exports","Imports",'net_movements')



BA2_net<-ggplot(BA2_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('red3','dodgerblue3'),name="",labels=c("Sink","Source"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=net_movements,fill=..y..>0))+
  #ggtitle("BA2")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('Net Exports')+
  theme(axis.text=element_text(size=7))+
  
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
BA2_net



BA2_absolutes<-ggplot(BA2_source_sink)+
  theme_minimal()+
  scale_fill_manual(values=c('dodgerblue4','red4'),name="",labels=c("Exports","Imports"))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=Exports, fill='Exports'))+
  geom_bar(stat='summary',fun.y='mean',aes(x = Origin_Continent,y=-Imports, fill='Imports'))+
  
  # ggtitle("BA2")+
  theme(legend.position = "none")+
  theme(legend.justification = "left")+
  ylab('No. Viral Movements')+
  theme(axis.text=element_text(size=7))+
  scale_y_continuous(limits=c(-5000,5000),breaks=seq(-5000,5000,1000),labels=abs(seq(-5000,5000,1000)))+
  theme(axis.title.x=element_text(size=8))+
  theme(axis.title.y=element_blank())+
  coord_flip()+
  theme(legend.key.size = unit(0.2, 'cm'), #change legend key size
        legend.key.height = unit(0.2, 'cm'), #change legend key height
        legend.key.width = unit(0.2, 'cm'))
BA2_absolutes


BA2_all_replicates_global_regional<-subset(subset(BA2_all_replicates, Origin_Continent!='character(0)'),Destination_Continent!='character(0)')
BA2_all_replicates_global_regional$Export_Type<- mapply(export2type,BA2_all_replicates_global_regional$Origin_Continent,BA2_all_replicates_global_regional$Destination_Continent)

BA2_all_replicates_global_regional_count<-BA2_all_replicates_global_regional %>%
  group_by(Origin,Export_Type)  %>%
  count(Export_Type)  %>%
  ungroup()

colnames(BA2_all_replicates_global_regional_count)<-c("Origin","Export_Type","Counts")
BA2_all_replicates_global_regional_count$Variant<-"Omicron BA2"


BA2_all_replicates_regional_per_region_count<-subset(BA2_all_replicates_global_regional,Export_Type=='Regional')
BA2_all_replicates_regional_per_region_count<-BA2_all_replicates_regional_per_region_count%>%
  group_by(Origin_Continent,Export_Type)%>%
  count()
colnames(BA2_all_replicates_regional_per_region_count)<-c("Origin_Continent","Export_Type","RegionCounts")

BA2_all_replicates_regional_count<-subset(BA2_all_replicates_global_regional_count,Export_Type=='Regional')
BA2_all_replicates_regional_count<-left_join(BA2_all_replicates_regional_count,BA2_all_replicates_global_regional%>%dplyr::select("Origin","Origin_Continent"))
BA2_all_replicates_regional_count<-left_join(BA2_all_replicates_regional_count,BA2_all_replicates_regional_per_region_count)

BA2_all_replicates_regional_count$ProportionsOfRegional<-BA2_all_replicates_regional_count$Counts/BA2_all_replicates_regional_count$RegionCounts
BA2_all_replicates_regional_count<-BA2_all_replicates_regional_count%>%dplyr::select("Origin","ProportionsOfRegional")


BA2_all_replicates_global_count<-subset(BA2_all_replicates_global_regional_count,Export_Type=='Global')
total_global_counts<-sum(BA2_all_replicates_global_count$Counts)
BA2_all_replicates_global_count$ProportionsOfGlobal<-BA2_all_replicates_global_count$Counts/total_global_counts
BA2_all_replicates_global_count<-BA2_all_replicates_global_count%>%dplyr::select("Origin","ProportionsOfGlobal")

BA2_all_replicates_regional_global_count<-left_join(BA2_all_replicates_regional_count,BA2_all_replicates_global_count)
BA2_all_replicates_regional_global_count$Variant<-"Omicron BA2"
BA2_all_replicates_regional_global_count<-unique(BA2_all_replicates_regional_global_count)

regional_cols<-c("Europe" = "#1b9e77", "Africa" = "#d95f02", "Asia" = "#7570b3", "North America" = "#e7298a",'South America'='gold2','Australia'='dodgerblue3')

d_BA2_links$mean_date<-as.Date(d_BA2_links$mean_date)
lab_dates <- pretty(d_BA2_links$mean_date)

BA2_map<-ggplot() +
  theme_void()+
  geom_map(data=world1,map=world1, aes(long, lat,map_id=region), color="white", fill='snow3',size=0.2)+
  
  geom_curve(data = d_BA2_links,
             aes(x = as.double(mean_origin_long), 
                 y = as.double(mean_origin_lat), 
                 xend = as.double(mean_destination_long), 
                 yend = as.double(mean_destination_lat), colour=as.Date(mean_date)),size=2,
             alpha=0.8)+
  
  geom_point(data = d_BA2_origins,
             aes(x = as.double(mean_origin_long), y = as.double(mean_origin_lat), size=n),
             color='black', shape=21)+
  scale_color_gradientn(colours = c('#c51b7d','#e9a3c9','#fde0ef','#e6f5d0','#a1d76a','#4d9221'),breaks = lab_dates, 
                        labels = lab_dates, name='Inferred Dispersal Dates (Mean)')+
  scale_fill_manual(values=regional_cols, name='Regions of Origin')+
  theme(legend.position = 'top')+
  theme(legend.direction = 'horizontal')+
  coord_fixed()+
  guides(colour = guide_colourbar(barwidth = 20, barheight = 0.5,title.position = 'top',title.hjust = 1,ticks.colour = "white",
                                  ticks.linewidth = 2),
         size='none')+
  theme(plot.title = element_text(family="Helvetica"))+
  scale_y_continuous(limits=c(-55,90))+
  scale_x_continuous(limits=c(-170,170))


BA2_map

BA2_fig2<-BA2_map +
  inset_element(BA2_net, 0.0, 0.0, 0.3, 0.5)
BA2_fig2

plot_grid(alpha_fig2,beta_fig2,gamma_fig2,delta_fig2,BA1_fig2,BA2_fig2,ncol=2,labels=c("Alpha","Beta","Gamma","Delta","Omicron BA.1","Omicron BA.2"))


#Figure 1B:

all_global_plot<-ggplot()+
  theme_minimal()+
  scale_fill_manual(values=c('#f3e8d6',"#ffd166","#d74e09","#2081c3","#99c5b5","#335145"))+
  geom_bar(data=subset(na.omit(all_both),total_global>0.005),color='black',size=0.1,stat='identity',position='stack',aes(fill=Variant,y = ProportionsOfGlobal,x=reorder(Origin,-total_global, na.last = TRUE)))+
  theme(plot.title = element_text(family="Helvetica"))+
  theme(legend.position = "top")+
  theme(legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))+
  theme(axis.title=element_blank()) +
  coord_trans(y="sqrt")+
  theme(rect = element_rect(fill = 'grey50'),
        line = element_blank(),
        panel.background = element_rect(fill = 'grey50'))
all_global_plot



africa_regional_plot<-ggplot()+
  theme_minimal()+
  scale_fill_manual(values=c('#f3e8d6',"#ffd166","#d74e09","#99c5b5","#335145"))+
  geom_bar(data=subset(na.omit(subset(all_both,Origin_Continent=='Africa')),total_global>0.001),color='black',size=0.1,stat='identity',position='stack',aes(fill=Variant,y = ProportionsOfRegional,x=reorder(Origin,-total_regional, na.last = TRUE)))+
  theme(plot.title = element_text(family="Helvetica"))+
  theme(legend.position = "none")+
  theme(legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8,colour='white'))+
  theme(axis.text.y=element_text(colour='white',size=8))+
  theme(axis.title=element_blank()) +
  coord_trans(y="sqrt")+
  theme(
    panel.background = element_rect(fill = 'white'))
africa_regional_plot


europe_regional_plot<-ggplot()+
  theme_minimal()+
  scale_fill_manual(values=c('#f3e8d6',"#ffd166","#d74e09","#2081c3","#99c5b5","#335145"))+
  geom_bar(data=subset(na.omit(subset(all_both,Origin_Continent=='Europe')),total_global>0.005),color='black',size=0.1,stat='identity',position='stack',aes(fill=Variant,y = ProportionsOfRegional,x=reorder(Origin,-total_regional, na.last = TRUE)))+
  theme(plot.title = element_text(family="Helvetica"))+
  theme(legend.position = "none")+
  theme(legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,color='white',size=8))+
  theme(axis.text.y=element_text(color='white',size=8))+
  
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10))+
  theme(axis.title=element_blank()) +
  coord_trans(y="sqrt")+
  theme(
    panel.background = element_rect(fill = 'white'))
europe_regional_plot


asia_regional_plot<-ggplot()+
  theme_minimal()+
  scale_fill_manual(values=c('#f3e8d6',"#ffd166","#d74e09","#99c5b5","#335145"))+
  geom_bar(data=subset(na.omit(subset(all_both,Origin_Continent=='Asia')),total_global>0.001),color='black',size=0.1,stat='identity',position='stack',aes(fill=Variant,y = ProportionsOfRegional,x=reorder(Origin,-total_regional, na.last = TRUE)))+
  theme(plot.title = element_text(family="Helvetica"))+
  theme(legend.position = "none")+
  theme(legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8,color='white'))+
  theme(axis.text.y=element_text(size=8,colour='white'))+
  theme(axis.title=element_blank()) +
  coord_trans(y="sqrt")+
  theme(
    panel.background = element_rect(fill = 'white'))
asia_regional_plot


samerica_regional_plot<-ggplot()+
  theme_minimal()+
  scale_fill_manual(values=c('#f3e8d6',"#d74e09","#2081c3","#99c5b5","#335145"))+
  geom_bar(data=subset(na.omit(subset(all_both,Origin_Continent=='South America')),total_global>0.001),color='black',size=0.1,stat='identity',position='stack',aes(fill=Variant,y = ProportionsOfRegional,x=reorder(Origin,-total_regional, na.last = TRUE)))+
  theme(plot.title = element_text(family="Helvetica"))+
  theme(legend.position = "none")+
  theme(legend.direction = "horizontal")+
  guides(fill = guide_legend(nrow = 1))+
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size=8,colour='white'))+
  theme(axis.text.y=element_text(size=8,colour='white'))+
  theme(axis.title=element_blank()) +
  coord_trans(y="sqrt")+
  theme(
    panel.background = element_rect(fill = 'white'))
samerica_regional_plot

contributions<-all_global_plot +
  inset_element(europe_regional_plot, 0.55, 0.1, 1, 0.6) +
  inset_element(asia_regional_plot, 0.1, 0.4, 0.4, 1) +
  inset_element(africa_regional_plot, 0.4, 0.48, 0.7, 1) +
  inset_element(samerica_regional_plot, 0.7, 0.55, 1, 1)
contributions



########## FIGURE 2 ###############


####### Figure 2A: Investigating exports of VOC from "source" country

#Alpha
alpha_all_replicates_table <- subset(subset(alpha_all_replicates,Origin!='UNKNOWN'),EventTime>2020.7) %>% count(date2, replicate,Origin=='United Kingdom')
colnames(alpha_all_replicates_table)<-c('date2', 'replicate', 'OriginUK', 'n')
alpha_all_replicates_table_summarize <- alpha_all_replicates_table %>% group_by(date2,OriginUK)  %>% 
  summarise(mean = mean(n), sd = sd(n))

#Calculating error bars for histogram
alpha_all_replicates_table_summarize$y_pos = NA

for (date in unique(alpha_all_replicates_table_summarize$date2)) {
  x=nrow(subset(subset(alpha_all_replicates_table_summarize, date2==date),OriginUK==TRUE))
  if (x==0) {
    alpha_all_replicates_table_summarize$y_pos[alpha_all_replicates_table_summarize$OriginUK == FALSE & alpha_all_replicates_table_summarize$date2 == date] = alpha_all_replicates_table_summarize$mean[alpha_all_replicates_table_summarize$OriginUK == FALSE & alpha_all_replicates_table_summarize$date2 == date]
  } else {
    alpha_all_replicates_table_summarize$y_pos[alpha_all_replicates_table_summarize$OriginUK == TRUE & alpha_all_replicates_table_summarize$date2 == date] = alpha_all_replicates_table_summarize$mean[alpha_all_replicates_table_summarize$OriginUK == TRUE & alpha_all_replicates_table_summarize$date2 == date]
    
    alpha_all_replicates_table_summarize$y_pos[alpha_all_replicates_table_summarize$OriginUK == FALSE & alpha_all_replicates_table_summarize$date2 == date] = alpha_all_replicates_table_summarize$mean[alpha_all_replicates_table_summarize$OriginUK == TRUE & alpha_all_replicates_table_summarize$date2 == date] + 
      alpha_all_replicates_table_summarize$mean[alpha_all_replicates_table_summarize$OriginUK == FALSE & alpha_all_replicates_table_summarize$date2 == date]
  }
}



alpha_all_replicates_table_summarize2<-alpha_all_replicates_table %>% group_by(date2,OriginUK) %>%
  summarise(n = n()) %>%
  mutate(freq = n / sum(n))

d_alpha<-subset(subset(alpha_all_replicates,Origin!='UNKNOWN'),EventTime>2020.7) %>%
  group_by(date2,replicate)  %>%
  mutate(count = n_distinct(Origin)) %>%
  group_by(date2)%>% 
  summarise(mean = mean(count), sd = sd(count))

p_alpha2<-ggplot()+
  theme_minimal()+
  scale_color_manual(values=c('purple4'),name='')+
  scale_fill_manual(values=c('darkorange1','dodgerblue2'), name='Inferred Origin',labels=c('Other Countries','United Kingdom'))+
  geom_bar(stat='identity',data=alpha_all_replicates_table_summarize,aes(x=date2,y=mean,fill=OriginUK),alpha=0.6,color='black', size=0.1)+
  
  geom_errorbar(data=alpha_all_replicates_table_summarize,aes(ymin=y_pos-sd, ymax=y_pos+sd,x=date2), width=.5,
                position='identity')+
  geom_ribbon(data=d_alpha,aes(x=date2, y=mean*8, ymin=(mean-sd)*8,ymax=(mean+sd)*8),fill='purple4',alpha=0.2)+
  geom_line(data=d_alpha,aes(x=date2,y=mean*8,colour='No. of\nSources'),show.legend=FALSE,shape=21,size=0.7)+
  ylab('')+
  xlab('All Inferred Introductions')+
  theme(legend.position='bottom')+
  theme(legend.direction='vertical')+
  
  ggtitle('Alpha')+
  theme(plot.title = element_text(hjust=0.5, family="Helvetica"))+
  theme(legend.key.size = unit(0.2, 'cm'),legend.title = element_text(size=8), legend.text = element_text(size=8))+ #change legend key size
  theme(axis.title.x=element_text(hjust=0.5,family='Helvetica',size=8), axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=8))+
  theme(axis.text.y = element_text(size=8))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 months", limits=as.Date(c("2020/09/01","2021/05/05")))+
  scale_y_continuous(
    name='',
    sec.axis = sec_axis(~ . / 8,name=''))
#p_alpha2

#Beta
beta_all_replicates_table <- beta_all_replicates %>% count(date2, replicate,Origin=='South Africa')
colnames(beta_all_replicates_table)<-c('date2', 'replicate', 'OriginSA', 'n')
beta_all_replicates_table_summarize <- beta_all_replicates_table %>% group_by(date2,OriginSA)  %>% 
  summarise(mean = mean(n), sd = sd(n))

beta_all_replicates_table_summarize$y_pos = NA



for (date in unique(beta_all_replicates_table_summarize$date2)) {
  x=nrow(subset(subset(beta_all_replicates_table_summarize, date2==date),OriginSA==TRUE))
  if (x==0) {
    #print('true is zero')
    #delta_all_replicates_table_summarize$y_pos[delta_all_replicates_table_summarize$OriginIndia == TRUE & delta_all_replicates_table_summarize$date2 == date] = 0
    
    beta_all_replicates_table_summarize$y_pos[beta_all_replicates_table_summarize$OriginSA == FALSE & beta_all_replicates_table_summarize$date2 == date] = beta_all_replicates_table_summarize$mean[beta_all_replicates_table_summarize$OriginSA == FALSE & beta_all_replicates_table_summarize$date2 == date]
    
  } else {
    beta_all_replicates_table_summarize$y_pos[beta_all_replicates_table_summarize$OriginSA == TRUE & beta_all_replicates_table_summarize$date2 == date] = beta_all_replicates_table_summarize$mean[beta_all_replicates_table_summarize$OriginSA == TRUE & beta_all_replicates_table_summarize$date2 == date]
    
    beta_all_replicates_table_summarize$y_pos[beta_all_replicates_table_summarize$OriginSA == FALSE & beta_all_replicates_table_summarize$date2 == date] = beta_all_replicates_table_summarize$mean[beta_all_replicates_table_summarize$OriginSA == TRUE & beta_all_replicates_table_summarize$date2 == date] + 
      beta_all_replicates_table_summarize$mean[beta_all_replicates_table_summarize$OriginSA == FALSE & beta_all_replicates_table_summarize$date2 == date]
  }
}





d_beta<-beta_all_replicates %>%
  group_by(date2,replicate)  %>%
  mutate(count = n_distinct(Origin)) %>%
  group_by(date2)%>% 
  summarise(mean = mean(count), sd = sd(count))

p_beta2<-ggplot()+
  theme_minimal()+
  scale_color_manual(values=c('purple4'),name='')+
  scale_fill_manual(values=c('darkorange1','dodgerblue2'), name='Inferred Origin',labels=c('Other Countries','South Africa'))+
  geom_bar(stat='identity',data=beta_all_replicates_table_summarize,aes(x=date2,y=mean,fill=OriginSA),alpha=0.6,color='black', size=0.1)+
  
  geom_errorbar(data=beta_all_replicates_table_summarize,aes(ymin=y_pos-sd, ymax=y_pos+sd,x=date2), width=.5,
                position='identity')+
  geom_ribbon(data=d_beta,aes(x=date2, y=mean*8, ymin=(mean-sd)*8,ymax=(mean+sd)*8),fill='purple4',alpha=0.2)+
  geom_line(data=d_beta,aes(x=date2,y=mean*8,colour='Number of\nSources'),show.legend = FALSE,shape=21,size=0.7)+
  ylab('')+
  xlab('All Inferred Introductions')+
  theme(legend.position='bottom')+
  theme(legend.direction='vertical')+
  
  ggtitle('Beta')+
  theme(plot.title = element_text(hjust=0.5, family="Helvetica"))+
  theme(legend.key.size = unit(0.2, 'cm'),legend.title = element_text(size=8), legend.text = element_text(size=8))+ #change legend key size
  theme(axis.title.x=element_text(hjust=0.5,family='Helvetica',size=8), axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=8))+
  theme(axis.text.y = element_text(size=8))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "3 months", limits=as.Date(c("2020/08/20","2021/06/15")))+
  scale_y_continuous(
    name='',
    sec.axis = sec_axis(~ . / 8,name=''))



#p_beta2

#Gamma
gamma_all_replicates_table <- subset(gamma_all_replicates,Origin!='UNKNOWN') %>% count(date2, replicate,Origin=='Brazil')
colnames(gamma_all_replicates_table)<-c('date2', 'replicate', 'OriginBrazil', 'n')
gamma_all_replicates_table_summarize <- gamma_all_replicates_table %>% group_by(date2,OriginBrazil)  %>% 
  summarise(mean = mean(n), sd = sd(n))

gamma_all_replicates_table_summarize$y_pos = NA



for (date in unique(gamma_all_replicates_table_summarize$date2)) {
  x=nrow(subset(subset(gamma_all_replicates_table_summarize, date2==date),OriginBrazil==TRUE))
  if (x==0) {
    gamma_all_replicates_table_summarize$y_pos[gamma_all_replicates_table_summarize$OriginBrazil == FALSE & gamma_all_replicates_table_summarize$date2 == date] = gamma_all_replicates_table_summarize$mean[gamma_all_replicates_table_summarize$OriginBrazil == FALSE & gamma_all_replicates_table_summarize$date2 == date]
    
  } else {
    gamma_all_replicates_table_summarize$y_pos[gamma_all_replicates_table_summarize$OriginBrazil == TRUE & gamma_all_replicates_table_summarize$date2 == date] = gamma_all_replicates_table_summarize$mean[gamma_all_replicates_table_summarize$OriginBrazil == TRUE & gamma_all_replicates_table_summarize$date2 == date]
    
    gamma_all_replicates_table_summarize$y_pos[gamma_all_replicates_table_summarize$OriginBrazil == FALSE & gamma_all_replicates_table_summarize$date2 == date] = gamma_all_replicates_table_summarize$mean[gamma_all_replicates_table_summarize$OriginBrazil == TRUE & gamma_all_replicates_table_summarize$date2 == date] + 
      gamma_all_replicates_table_summarize$mean[gamma_all_replicates_table_summarize$OriginBrazil == FALSE & gamma_all_replicates_table_summarize$date2 == date]
  }
}


d_gamma<-subset(gamma_all_replicates,Origin!='UNKNOWN') %>%
  group_by(date2,replicate)  %>%
  mutate(count = n_distinct(Origin)) %>%
  group_by(date2)%>% 
  summarise(mean = mean(count), sd = sd(count))



p_gamma2<-ggplot()+
  theme_minimal()+
  scale_color_manual(values=c('purple4'),name='')+
  scale_fill_manual(values=c('darkorange1','dodgerblue2'), name='Inferred Origin',labels=c('Other Countries','Brazil'))+
  geom_bar(stat='identity',data=gamma_all_replicates_table_summarize,aes(x=date2,y=mean,fill=OriginBrazil),alpha=0.6,color='black', size=0.1)+
  
  geom_errorbar(data=gamma_all_replicates_table_summarize,aes(ymin=y_pos-sd, ymax=y_pos+sd,x=date2), width=.5,
                position='identity')+
  geom_ribbon(data=d_gamma,aes(x=date2, y=mean*8, ymin=(mean-sd)*8,ymax=(mean+sd)*8),fill='purple4',alpha=0.2)+
  geom_line(data=d_gamma,aes(x=date2,y=mean*8,colour='Number of\nSources'), show.legend=FALSE,shape=21,size=0.7)+
  ylab('')+
  xlab('All Inferred Introductions')+
  theme(legend.position='bottom')+
  theme(legend.direction='vertical')+
  
  ggtitle('Gamma')+
  theme(plot.title = element_text(hjust=0.5, family="Helvetica"))+
  theme(legend.key.size = unit(0.2, 'cm'),legend.title = element_text(size=8), legend.text = element_text(size=8))+ #change legend key size
  theme(axis.title.x=element_text(hjust=0.5,family='Helvetica',size=8), axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=8))+
  theme(axis.text.y = element_text(size=8))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "2 months", limits=as.Date(c("2020/11/01","2021/04/30")))+
  scale_y_continuous(
    name='',
    sec.axis = sec_axis(~ . / 8,name=''))



#p_gamma2


#Delta
delta_all_replicates_table <- subset(delta_all_replicates,Origin!='UNKNOWN') %>% count(date2, replicate,Origin=='India')
colnames(delta_all_replicates_table)<-c('date2', 'replicate', 'OriginIndia', 'n')
delta_all_replicates_table_summarize <- delta_all_replicates_table %>% group_by(date2,OriginIndia)  %>% 
  summarise(mean = mean(n), sd = sd(n))

delta_all_replicates_table_summarize$y_pos = NA



for (date in unique(delta_all_replicates_table_summarize$date2)) {
  x=nrow(subset(subset(delta_all_replicates_table_summarize, date2==date),OriginIndia==TRUE))
  if (x==0) {
    delta_all_replicates_table_summarize$y_pos[delta_all_replicates_table_summarize$OriginIndia == FALSE & delta_all_replicates_table_summarize$date2 == date] = delta_all_replicates_table_summarize$mean[delta_all_replicates_table_summarize$OriginIndia == FALSE & delta_all_replicates_table_summarize$date2 == date]
    
  } else {
    delta_all_replicates_table_summarize$y_pos[delta_all_replicates_table_summarize$OriginIndia == TRUE & delta_all_replicates_table_summarize$date2 == date] = delta_all_replicates_table_summarize$mean[delta_all_replicates_table_summarize$OriginIndia == TRUE & delta_all_replicates_table_summarize$date2 == date]
    
    delta_all_replicates_table_summarize$y_pos[delta_all_replicates_table_summarize$OriginIndia == FALSE & delta_all_replicates_table_summarize$date2 == date] = delta_all_replicates_table_summarize$mean[delta_all_replicates_table_summarize$OriginIndia == TRUE & delta_all_replicates_table_summarize$date2 == date] + 
      delta_all_replicates_table_summarize$mean[delta_all_replicates_table_summarize$OriginIndia == FALSE & delta_all_replicates_table_summarize$date2 == date]
  }
}


d_delta<-subset(delta_all_replicates,Origin!='UNKNOWN') %>%
  group_by(date2,replicate)  %>%
  mutate(count = n_distinct(Origin)) %>%
  group_by(date2)%>% 
  summarise(mean = mean(count), sd = sd(count))


p_delta2<-ggplot()+
  theme_minimal()+
  scale_color_manual(values=c('purple4'),name='')+
  scale_fill_manual(values=c('darkorange1','dodgerblue2'), name='Inferred Origin',labels=c('Other Countries','India'))+
  geom_bar(stat='identity',data=delta_all_replicates_table_summarize,aes(x=date2,y=mean,fill=OriginIndia),alpha=0.6,color='black', size=0.1)+
  
  geom_errorbar(data=delta_all_replicates_table_summarize,aes(ymin=y_pos-sd, ymax=y_pos+sd,x=date2), width=.5,
                position='identity')+
  geom_ribbon(data=d_delta,aes(x=date2, y=mean*8, ymin=(mean-sd)*8,ymax=(mean+sd)*8),fill='purple4',alpha=0.2)+
  geom_line(data=d_delta,aes(x=date2,y=mean*8,colour='Number of\nSources'),show.legend=FALSE,shape=21,size=0.7)+
  ylab('')+
  xlab('All Inferred Introductions')+
  theme(legend.position='bottom')+
  theme(legend.direction='vertical')+
  
  ggtitle('Delta')+
  theme(plot.title = element_text(hjust=0.5, family="Helvetica"))+
  theme(legend.key.size = unit(0.2, 'cm'),legend.title = element_text(size=8), legend.text = element_text(size=8))+ #change legend key size
  theme(axis.title.x=element_text(hjust=0.5,family='Helvetica',size=8), axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=8))+
  theme(axis.text.y = element_text(size=8))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "4 months", limits=as.Date(c("2020/08/10","2021/11/25")))+
  scale_y_continuous(
    name='',
    sec.axis = sec_axis(~ . / 8,name=''))



#p_delta2


#BA1

BA1_all_replicates %>% group_by(replicate,Origin=='South Africa') %>%  summarise(n = n()) %>% mutate(perc = (n / sum(n))*100) ### percentage from SA per replicate

BA1_all_replicates_table <- BA1_all_replicates %>% count(date, replicate,Origin=='South Africa')
colnames(BA1_all_replicates_table)<-c('date', 'replicate', 'OriginSA', 'n')
BA1_all_replicates_table_summarize <- BA1_all_replicates_table %>% group_by(date,OriginSA)  %>% 
  summarise(mean = mean(n), sd = sd(n))

BA1_all_replicates_table_summarize$y_pos = NA



for (date in unique(BA1_all_replicates_table_summarize$date)) {
  x=nrow(subset(subset(BA1_all_replicates_table_summarize, date==date),OriginSA==TRUE))
  if (x==0) {
    
    BA1_all_replicates_table_summarize$y_pos[BA1_all_replicates_table_summarize$OriginSA == FALSE & BA1_all_replicates_table_summarize$date == date] = BA1_all_replicates_table_summarize$mean[BA1_all_replicates_table_summarize$OriginSA == FALSE & BA1_all_replicates_table_summarize$date == date]
    
  } else {
    BA1_all_replicates_table_summarize$y_pos[BA1_all_replicates_table_summarize$OriginSA == TRUE & BA1_all_replicates_table_summarize$date == date] = BA1_all_replicates_table_summarize$mean[BA1_all_replicates_table_summarize$OriginSA == TRUE & BA1_all_replicates_table_summarize$date == date]
    
    BA1_all_replicates_table_summarize$y_pos[BA1_all_replicates_table_summarize$OriginSA == FALSE & BA1_all_replicates_table_summarize$date == date] = BA1_all_replicates_table_summarize$mean[BA1_all_replicates_table_summarize$OriginSA == TRUE & BA1_all_replicates_table_summarize$date == date] + 
      BA1_all_replicates_table_summarize$mean[BA1_all_replicates_table_summarize$OriginSA == FALSE & BA1_all_replicates_table_summarize$date == date]
  }
}


d_BA1<-BA1_all_replicates %>%
  group_by(date,replicate)  %>%
  mutate(count = n_distinct(Origin)) %>%
  group_by(date)%>% 
  summarise(mean = mean(count), sd = sd(count))


p_BA12<-ggplot()+
  theme_minimal()+
  scale_color_manual(values=c('purple4'),name='')+
  scale_fill_manual(values=c('darkorange1','dodgerblue2'), name='Inferred Origin',labels=c('Other Countries','South Africa'))+
  geom_bar(stat='identity',data=BA1_all_replicates_table_summarize,aes(x=date,y=mean,fill=OriginSA),alpha=0.6,color='black', size=0.1)+
  
  geom_errorbar(data=BA1_all_replicates_table_summarize,aes(ymin=y_pos-sd, ymax=y_pos+sd,x=date), width=.5,
                position='identity')+
  geom_ribbon(data=d_BA1,aes(x=date, y=mean*20, ymin=(mean-sd)*20,ymax=(mean+sd)*20),fill='purple4',alpha=0.2)+
  geom_line(data=d_BA1,aes(x=date,y=mean*20,colour='Number of\nSources'),show.legend=FALSE,shape=21,size=0.7)+
  ylab('')+
  xlab('All Inferred Introductions')+
  theme(legend.position='bottom')+
  theme(legend.direction='vertical')+
  
  ggtitle('Omicron BA.1')+
  theme(plot.title = element_text(hjust=0.5, family="Helvetica"))+
  theme(legend.key.size = unit(0.2, 'cm'),legend.title = element_text(size=8), legend.text = element_text(size=8))+ #change legend key size
  theme(axis.title.x=element_text(hjust=0.5,family='Helvetica',size=8), axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=8))+
  theme(axis.text.y = element_text(size=8))+
  
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "1 months", limits=as.Date(c("2021/11/05","2022/02/10")))+
  scale_y_continuous(
    name='',
    sec.axis = sec_axis(~ . / 20,name=''))

#p_BA12

#BA2

BA2_all_replicates_table <- BA2_all_replicates %>% count(date, replicate,Origin=='South Africa')
colnames(BA2_all_replicates_table)<-c('date', 'replicate', 'OriginSA', 'n')
BA2_all_replicates_table_summarize <- BA2_all_replicates_table %>% group_by(date,OriginSA)  %>% 
  summarise(mean = mean(n), sd = sd(n))

BA2_all_replicates_table_summarize$y_pos = NA



for (date in unique(BA2_all_replicates_table_summarize$date)) {
  x=nrow(subset(subset(BA2_all_replicates_table_summarize, date==date),OriginSA==TRUE))
  if (x==0) {
    BA2_all_replicates_table_summarize$y_pos[BA2_all_replicates_table_summarize$OriginSA == FALSE & BA2_all_replicates_table_summarize$date == date] = BA2_all_replicates_table_summarize$mean[BA2_all_replicates_table_summarize$OriginSA == FALSE & BA2_all_replicates_table_summarize$date == date]
    
  } else {
    BA2_all_replicates_table_summarize$y_pos[BA2_all_replicates_table_summarize$OriginSA == TRUE & BA2_all_replicates_table_summarize$date == date] = BA2_all_replicates_table_summarize$mean[BA2_all_replicates_table_summarize$OriginSA == TRUE & BA2_all_replicates_table_summarize$date == date]
    
    BA2_all_replicates_table_summarize$y_pos[BA2_all_replicates_table_summarize$OriginSA == FALSE & BA2_all_replicates_table_summarize$date == date] = BA2_all_replicates_table_summarize$mean[BA2_all_replicates_table_summarize$OriginSA == TRUE & BA2_all_replicates_table_summarize$date == date] + 
      BA2_all_replicates_table_summarize$mean[BA2_all_replicates_table_summarize$OriginSA == FALSE & BA2_all_replicates_table_summarize$date == date]
  }
}


d_BA2<-BA2_all_replicates %>%
  group_by(date,replicate)  %>%
  mutate(count = n_distinct(Origin)) %>%
  group_by(date)%>% 
  summarise(mean = mean(count), sd = sd(count))



p_BA22<-ggplot()+
  theme_minimal()+
  scale_color_manual(values=c('purple4'),name='')+
  scale_fill_manual(values=c('darkorange1','dodgerblue2'), name='Inferred Origin',labels=c('Other Countries','South Africa'))+
  geom_bar(stat='identity',data=BA2_all_replicates_table_summarize,aes(x=date,y=mean,fill=OriginSA),alpha=0.6,color='black', size=0.1)+
  
  geom_errorbar(data=BA2_all_replicates_table_summarize,aes(ymin=y_pos-sd, ymax=y_pos+sd,x=date), width=.5,
                position='identity')+
  geom_ribbon(data=d_BA2,aes(x=date, y=mean*20, ymin=(mean-sd)*20,ymax=(mean+sd)*20),fill='purple4',alpha=0.2)+
  geom_line(data=d_BA2,aes(x=date,y=mean*20,colour='Number of\nSources'),show.legend=FALSE,shape=21,size=0.7)+
  ylab('')+
  xlab('All Inferred Introductions')+
  theme(legend.position='bottom')+
  theme(legend.direction='vertical')+
  
  ggtitle('Omicron BA.2')+
  theme(plot.title = element_text(hjust=0.5, family="Helvetica"))+
  theme(legend.key.size = unit(0.2, 'cm'),legend.title = element_text(size=8), legend.text = element_text(size=8))+ #change legend key size
  theme(axis.title.x=element_text(hjust=0.5,family='Helvetica',size=8), axis.title.y=element_blank())+
  theme(axis.text.x = element_text(size=8))+
  theme(axis.text.y = element_text(size=8))+
  
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "1 months", limits=as.Date(c("2021/11/25","2022/03/05")))+
  scale_y_continuous(
    name='',
    sec.axis = sec_axis(~ . / 20,name=''))

#p_BA22


###Figure 2A:
plot_grid(p_alpha2,p_beta2,p_gamma2,p_delta2,p_BA12,p_BA22)


### Figure 2B: Investigating where the earliest VOC introduction per country came from (original source country or not)

#get dataset of first VOC introduction per country (From inferred import/export analysis)
#then count source as being fromoriginal country vs others
#plot over time to see as time goes by, if source country less likely to be original country

#Alpha

d_alpha_earliest<-subset(alpha_all_replicates,Origin!='UNKNOWN') %>%
  group_by(replicate, Destination)  %>%
  arrange(date) %>% 
  slice(1L)

hundred_after_tmrca_alpha<-alpha_first_date+100

p_alpha_earliest_intros<-ggplot(data=subset(d_alpha_earliest,Destination!='United Kingdom'),aes(x=days, y=reorder(Destination,days,median)))+
  theme_bw()+ 
  geom_smooth(method='loess',se=FALSE,size=0.7,aes(x=days, y=reorder(Destination,days,median),group=1),color='grey60')+
  geom_vline(xintercept = hundred_after_tmrca_alpha,size=0.7,linetype='dashed',colour='red3')+
  
  geom_point(alpha=0.3,size=2,aes(colour=Origin=='United Kingdom'), position=position_jitter(h=0, w=0))+
  
  scale_colour_manual(values=c('darkorange2','dodgerblue3'),name='Inferred Origin',labels=c('Other Places','UK'))+
  theme(axis.text.y = element_text(size=7))+
  ylab('')+
  xlab('')+
  
  ggtitle('First Inferred Introductions')+
  theme(plot.title=element_text(hjust=0.5,family='Helvetica',size=8))+
  theme(axis.title=element_blank(), axis.text.x = element_text(size=8))+
  theme(legend.position="none",legend.background = element_rect(colour = NA, fill = NA),legend.key = element_rect(colour = NA, fill = NA))+
  scale_y_discrete(labels=rev(1:length(unique(subset(d_alpha_earliest,Destination!='United Kingdom')$Destination))),country_numbers,limits=rev)+
  
  scale_x_date(position='top',date_labels = "%b\n%Y",date_breaks = "2 months", limits=as.Date(c("2020/09/01","2021/05/05")))

p_alpha_earliest_intros

Alpha_movement<-plot_grid(p_alpha2,p_alpha_earliest_intros,ncol=1,align='v',rel_heights = c(0.2,0.8))
#Alpha_movement

#Beta

d_beta_earliest<-beta_all_replicates %>%
  group_by(replicate, Destination)  %>%
  arrange(date) %>% 
  slice(1L)

hundred_after_tmrca_beta<-beta_first_date+100
p_beta_earliest_intros<-ggplot(data=subset(d_beta_earliest,Destination!='South Africa'),aes(x=days, y=reorder(Destination,days,median)))+
  theme_bw()+ 
  geom_smooth(method='loess',se=FALSE,size=0.7,aes(x=days, y=reorder(Destination,days,median),group=1),color='grey60')+
  geom_vline(xintercept = hundred_after_tmrca_beta,size=0.7,linetype='dashed',colour='red3')+
  
  geom_point(alpha=0.3,size=2,aes(colour=Origin=='South Africa'), position=position_jitter(h=0, w=0))+
  scale_colour_manual(values=c('darkorange2','dodgerblue3'),name='Inferred Origin',labels=c('Other Places','South Africa'))+
  theme(axis.text.y = element_text(size=7))+
  ylab('')+
  xlab('')+
  
  ggtitle('First Inferred Introductions')+
  theme(plot.title=element_text(hjust=0.5,family='Helvetica',size=8))+
  theme(axis.title=element_blank(), axis.text.x = element_text(size=8))+
  theme(legend.position="none",legend.background = element_rect(colour = NA, fill = NA),legend.key = element_rect(colour = NA, fill = NA))+
  scale_y_discrete(labels=rev(1:length(unique(subset(d_beta_earliest,Destination!='South Africa')$Destination))),limits=rev)+
  
  scale_x_date(position='top',date_labels = "%b\n%Y",date_breaks = "3 months", limits=as.Date(c("2020/08/20","2021/06/15")))
p_beta_earliest_intros

Beta_movement<-plot_grid(p_beta2,p_beta_earliest_intros,ncol=1,align='v',rel_heights = c(0.2,0.8))
#Beta_movement


#Gamma
d_gamma_earliest<-subset(gamma_all_replicates,Origin!='UNKNOWN') %>%
  group_by(replicate, Destination)  %>%
  arrange(date) %>% 
  slice(1L)

d_gamma_earliest$arrival<-as.numeric(d_gamma_earliest$days-gamma_first_date)

label_minor<-pretty(subset(d_gamma_earliest,arrival>0)$arrival)

hundred_after_tmrca_gamma<-gamma_first_date+100
p_gamma_earliest_intros<-ggplot(data=subset(d_gamma_earliest,Destination!='Brazil'),aes(x=days, y=reorder(Destination,days,median)))+
  theme_bw()+ 
  geom_smooth(method='loess',se=FALSE,size=0.7,aes(x=days, y=reorder(Destination,days,median),group=1),color='grey60')+
  geom_vline(xintercept = hundred_after_tmrca_gamma,size=0.7,linetype='dashed',colour='red3')+
  
  geom_point(alpha=0.3,size=2,aes(colour=Origin=='Brazil'), position=position_jitter(h=0, w=0))+
  scale_colour_manual(values=c('darkorange2','dodgerblue3'),name='Inferred Origin',labels=c('Other Places','Brazil'))+
  theme(axis.text.y = element_text(size=7))+
  ylab('')+
  xlab('')+
  
  ggtitle('First Inferred Introductions')+
  theme(plot.title=element_text(hjust=0.5,family='Helvetica',size=8))+
  theme(axis.title=element_blank(), axis.text.x = element_text(size=8))+
  theme(legend.position="none",legend.background = element_rect(colour = NA, fill = NA),legend.key = element_rect(colour = NA, fill = NA))+
  scale_y_discrete(labels=rev(1:length(unique(subset(d_gamma_earliest,Destination!='Brazil')$Destination))),limits=rev)+
  
  scale_x_date(position='top',date_labels = "%b\n%Y",date_breaks = "2 months", limits=as.Date(c("2020/11/01","2021/04/30")))
p_gamma_earliest_intros

Gamma_movement<-plot_grid(p_gamma2,p_gamma_earliest_intros,ncol=1,align='v',rel_heights = c(0.2,0.8))
#Gamma_movement


#Delta
d_delta_earliest<-subset(delta_all_replicates,Origin!='UNKNOWN') %>%
  group_by(replicate, Destination)  %>%
  arrange(date) %>% 
  slice(1L)
hundred_after_tmrca_delta<-delta_first_date+100

p_delta_earliest_intros<-ggplot(data=subset(d_delta_earliest,Destination!='India'),aes(x=days, y=reorder(Destination,days)))+
  theme_bw()+ 
  geom_smooth(method='loess',se=FALSE,size=0.7,aes(x=days, y=reorder(Destination,days,median),group=1),color='grey60')+
  geom_vline(xintercept = hundred_after_tmrca_delta,size=0.7,linetype='dashed',colour='red3')+
  
  geom_point(alpha=0.3,size=2,aes(colour=Origin=='India'), position=position_jitter(h=0, w=0))+
  scale_colour_manual(values=c('darkorange2','dodgerblue3'),name='Inferred Origin',labels=c('Other Places','India'))+
  theme(axis.text.y = element_text(size=7))+
  ylab('')+
  xlab('')+
  ggtitle('First Inferred Introductions')+
  theme(plot.title=element_text(hjust=0.5,family='Helvetica',size=8))+
  theme(axis.title=element_blank(), axis.text.x = element_text(size=8))+
  theme(legend.position="none",legend.background = element_rect(colour = NA, fill = NA),legend.key = element_rect(colour = NA, fill = NA))+
  scale_y_discrete(labels=rev(1:length(unique(subset(d_delta_earliest,Destination!='India')$Destination))),limits=rev)+
  
  scale_x_date(position='top',date_labels = "%b\n%Y",date_breaks = "4 months", limits=as.Date(c("2020/08/10","2021/11/25")))

p_delta_earliest_intros

Delta_movement<-plot_grid(p_delta2,p_delta_earliest_intros,ncol=1,align='v',rel_heights = c(0.2,0.8))
#Delta_movement

#BA1

d_BA1_earliest<-BA1_all_replicates %>%
  group_by(replicate, Destination)  %>%
  arrange(days) %>% 
  slice(1L)

d_BA1_earliest_median<-d_BA1_earliest%>%
  group_by(Destination)  %>%
  summarise(median_days=median(days))

hundred_after_tmrca_ba1<-ba1_first_date+100

p_BA1_earliest_intros<-ggplot(data=subset(d_BA1_earliest,Destination!='South Africa'),aes(x=days, y=reorder(Destination,days)))+
  theme_bw()+ 
  geom_smooth(method='loess',se=FALSE,size=0.7,aes(x=days, y=reorder(Destination,days,median),group=1),color='grey60')+
  geom_vline(xintercept = hundred_after_tmrca_ba1,size=0.7,linetype='dashed',colour='red3')+
  
  geom_point(alpha=0.3,size=2,aes(colour=Origin=='South Africa'), position=position_jitter(h=0, w=0))+
  scale_colour_manual(values=c('darkorange2','dodgerblue3'),name='Inferred Origin',labels=c('Other Places','India'))+
  theme(axis.text.y = element_text(size=7))+
  ylab('')+
  xlab('')+
  ggtitle('First Inferred Introductions')+
  theme(plot.title=element_text(hjust=0.5,family='Helvetica',size=8))+
  theme(axis.title=element_blank(), axis.text.x = element_text(size=8))+
  theme(legend.position="none",legend.background = element_rect(colour = NA, fill = NA),legend.key = element_rect(colour = NA, fill = NA))+
  scale_y_discrete(labels=rev(1:length(unique(subset(d_BA1_earliest,Destination!='South Africa')$Destination))),limits=rev)+
  
  scale_x_date(position='top',date_labels = "%b\n%Y",date_breaks = "1 months", limits=as.Date(c("2021/11/05","2022/02/10")))

p_BA1_earliest_intros

BA1_movement<-plot_grid(p_BA12,p_BA1_earliest_intros,ncol=1,align='v',rel_heights = c(0.2,0.8))
#BA1_movement


#BA2

d_BA2_earliest<-BA2_all_replicates %>%
  group_by(replicate, Destination)  %>%
  arrange(date) %>% 
  slice(1L)

country_numbers<-1:100
hundred_after_tmrca_ba2<-ba2_first_date+100

p_BA2_earliest_intros<-ggplot(data=subset(d_BA2_earliest,Destination!='South Africa'),aes(x=days, y=reorder(Destination,days,median)))+
  theme_bw()+ 
  geom_smooth(method='loess',se=FALSE,size=0.7,aes(x=days, y=reorder(Destination,days,median),group=1),color='grey60')+
  geom_vline(xintercept = hundred_after_tmrca_ba2,size=0.7,linetype='dashed',colour='red3')+
  
  geom_point(alpha=0.3,size=2,aes(colour=Origin=='South Africa'), position=position_jitter(h=0, w=0))+
  scale_colour_manual(values=c('darkorange2','dodgerblue3'),name='Inferred Origin',labels=c('Other Places','South Africa'))+
  theme(axis.text.y = element_text(size=7))+
  ylab('')+
  xlab('')+
  ggtitle('First Inferred Introductions')+
  theme(plot.title=element_text(hjust=0.5,family='Helvetica',size=8))+
  theme(axis.title=element_blank(), axis.text.x = element_text(size=8))+
  theme(legend.position="none",legend.background = element_rect(colour = NA, fill = NA),legend.key = element_rect(colour = NA, fill = NA))+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "1 months", limits=as.Date(c("2021/11/25","2022/03/05")),position='top')+
  scale_y_discrete(labels=rev(1:length(unique(subset(d_BA2_earliest,Destination!='South Africa')$Destination))),limits=rev)
p_BA2_earliest_intros

BA2_movement<-plot_grid(p_BA22,p_BA2_earliest_intros,ncol=1,align='v',rel_heights = c(0.2,0.8))
#BA2_movement

fig2<-plot_grid(Alpha_movement,Beta_movement,Gamma_movement,Delta_movement,BA1_movement,BA2_movement,ncol=6)

fig2
#ggsave(file="FIG2_panel.pdf",plot=fig3,width = 12.35, height = 10, units = "in", limitsize = FALSE)


########## FIGURE 3 ###############


#Reading genomic metadata (Takes some time to annotate)
gisaid_data_original<-read.csv('metadata_2022-09-17_00-48.tsv',sep='\t')

gisaid_data<-gisaid_data_original
gisaid_data$date<-as.Date(gisaid_data$date, format="%Y-%m-%d")
gisaid_data$date2<-as.Date(cut(gisaid_data$date,breaks = "1 weeks",start.on.monday = FALSE))
gisaid_data$days<-as.Date(cut(gisaid_data$date,breaks = "day",start.on.monday = FALSE))
gisaid_data$date4<-as.Date(cut(gisaid_data$date,breaks = "1 month",start.on.monday = FALSE))


VOCs_df<-gisaid_data[grepl("B.1.1.7/.",gisaid_data$pango_lineage) | gisaid_data$pango_lineage=='B.1.1.7' | grepl("B.1.617.2",gisaid_data$pango_lineage) | grepl("AY.",gisaid_data$pango_lineage) |
                       grepl("B.1.351",gisaid_data$pango_lineage) | grepl("P.1.",gisaid_data$pango_lineage) | gisaid_data$pango_lineage=='P.1' |
                       grepl("BA.",gisaid_data$pango_lineage) | grepl("B.1.1.529",gisaid_data$pango_lineage),]

VOC_list<-unique(as.vector(VOCs_df$pango_lineage))

shouldBecomeOther<-!(gisaid_data$pango_lineage %in% VOC_list)

gisaid_data$Nextstrain_clade[shouldBecomeOther]<- "Others"

alpha<-grepl("B.1.1.7/.",gisaid_data$pango_lineage) | gisaid_data$pango_lineage=='B.1.1.7'
delta<-grepl("B.1.617.2",gisaid_data$pango_lineage) | grepl("AY.",gisaid_data$pango_lineage) 
BA1<-grepl("BA.1",gisaid_data$pango_lineage) | grepl("B.1.1.529",gisaid_data$pango_lineage)
BA2<-grepl("BA.2",gisaid_data$pango_lineage)
BA4_5<-grepl("BA.4",gisaid_data$pango_lineage) | grepl("BA.5",gisaid_data$pango_lineage)

beta<-grepl("B.1.351",gisaid_data$pango_lineage)
gamma<-grepl("P.1.",gisaid_data$pango_lineage) | gisaid_data$pango_lineage=='P.1'

gisaid_data$Nextstrain_clade[alpha]<- "Alpha"
gisaid_data$Nextstrain_clade[delta]<- "Delta"
gisaid_data$Nextstrain_clade[BA1]<- "BA1"

gisaid_data$Nextstrain_clade[BA2]<- "BA2"
gisaid_data$Nextstrain_clade[BA4_5]<- "BA4_5"

gisaid_data$Nextstrain_clade[beta]<- "Beta"
gisaid_data$Nextstrain_clade[gamma]<- "Gamma"


early_Alphas<-subset(gisaid_data, Nextstrain_clade=='Alpha' & date<alpha_first_date)$strain
early_Betas<-subset(gisaid_data, Nextstrain_clade=='Beta' & date<beta_first_date)$strain
early_Gammas<-subset(gisaid_data, Nextstrain_clade=='Gamma' & date<gamma_first_date)$strain
early_Deltas<-subset(gisaid_data, Nextstrain_clade=='Delta' & date<delta_first_date)$strain
early_BA1s<-subset(gisaid_data, Nextstrain_clade=='BA1' & date<ba1_first_date)$strain
early_BA2s<-subset(gisaid_data, Nextstrain_clade=='BA2' & date<ba2_first_date)$strain
early_BA4_5s<-subset(gisaid_data, Nextstrain_clade=='BA4_5' & date<ba4_5_first_date)$strain


gisaid_data <- gisaid_data[ !(gisaid_data$strain %in% c(early_Alphas)), ]
gisaid_data <- gisaid_data[ !(gisaid_data$strain %in% c(early_Betas)), ]
gisaid_data <- gisaid_data[ !(gisaid_data$strain %in% c(early_Gammas)), ]
gisaid_data <- gisaid_data[ !(gisaid_data$strain %in% c(early_Deltas)), ]
gisaid_data <- gisaid_data[ !(gisaid_data$strain %in% c(early_BA1s)), ]
gisaid_data <- gisaid_data[ !(gisaid_data$strain %in% c(early_BA2s)), ]
gisaid_data <- gisaid_data[ !(gisaid_data$strain %in% c(early_BA4_5s)), ]


gisaid_data$date_submitted<-as.Date(gisaid_data$date_submitted)
gisaid_data$days_submitted<-as.Date(cut(gisaid_data$date_submitted,breaks = "day",start.on.monday = FALSE))
gisaid_data$submission_lag=gisaid_data$days_submitted-gisaid_data$days

first_VOC_country<- subset(subset(gisaid_data,Nextstrain_clade!='Others'),submission_lag<366) %>% group_by(country, Nextstrain_clade) %>% slice(which.min(date))
first_VOC_country<-first_VOC_country %>% dplyr::select("country","Nextstrain_clade","date")
first_VOC_country_table<-spread(first_VOC_country, Nextstrain_clade, date)

first_VOC_country_detection<- subset(subset(gisaid_data,Nextstrain_clade!='Others'),submission_lag<366) %>% group_by(country, Nextstrain_clade) %>% slice(which.min(date_submitted))
first_VOC_country_detection<-first_VOC_country_detection %>% dplyr::select("country","Nextstrain_clade","date_submitted")
first_VOC_country_detection_table<-spread(first_VOC_country_detection, Nextstrain_clade, date_submitted)


###### Travel data

###Total global

#total number of international passengers aggregated by month from 2020 to 2022
aggregate_total <- read.csv('../Data/Global_passengers_monthly.csv')
head(aggregate_total)

#total number of international passengers aggregated by country from 2020 to 2022
country_total <- read.csv('../Data/Global_passengers_per_country_2020_2022.csv')

###Travel from presumed origin countries for each VOC during their relevant period of global dissemination
#B117-UK September 2020 - March 2021
b117_UK_travel_data<-read_csv("../Data/UK_passengers_b117_period.csv")

#B1351-SA September 2020 - March 2021
b1351_SA_travel_data<-read_csv("../Data/SA_passengers_b351_period.csv")

#P1_Brazil November 2020 - May 2021
P1_BRA_travel_data<-read_csv("../Data/Brazil_passengers_p1_period.csv")

#Delta_India September 2020 - September 2021
Delta_IND_travel_data<-read_csv("../Data/India_passengers_delta_period.csv")


#Omicron_SA November 2021 - March 2022

Omicron_SA_travel_data<-read_csv("../Data/SA_passengers_omicron_period.csv")



# Global VOC introduction, arrival delay vs travel

first_VOC_country_table$Alpha_arrival=first_VOC_country_table$Alpha-alpha_first_date
first_VOC_country_table$Beta_arrival=first_VOC_country_table$Beta-beta_first_date
first_VOC_country_table$Delta_arrival=first_VOC_country_table$Delta-delta_first_date
first_VOC_country_table$Gamma_arrival=first_VOC_country_table$Gamma-gamma_first_date
first_VOC_country_table$BA1_arrival=first_VOC_country_table$BA1-ba1_first_date
first_VOC_country_table$BA2_arrival=first_VOC_country_table$BA2-ba2_first_date
first_VOC_country_table$BA4_5_arrival=first_VOC_country_table$BA2-ba4_5_first_date


alpha_inferred_arrival<-d_alpha_earliest ##### FROM FIG2 script
alpha_inferred_arrival<-alpha_inferred_arrival %>% 
  dplyr::select("Origin",'Destination','replicate','date') %>% 
  group_by(Destination) %>%
  summarize(alpha_inferred_arrival_mean_date=mean(date))
colnames(alpha_inferred_arrival)<-c('country','alpha_inferred_arrival_mean_date')
alpha_inferred_arrival


beta_inferred_arrival<-d_beta_earliest ##### FROM FIG2 script
beta_inferred_arrival<-beta_inferred_arrival %>% 
  dplyr::select("Origin",'Destination','replicate','date') %>% 
  group_by(Destination) %>%
  summarize(beta_inferred_arrival_mean_date=mean(date))
colnames(beta_inferred_arrival)<-c('country','beta_inferred_arrival_mean_date')
beta_inferred_arrival

gamma_inferred_arrival<-d_gamma_earliest ##### FROM FIG2 script
gamma_inferred_arrival<-gamma_inferred_arrival %>% 
  dplyr::select("Origin",'Destination','replicate','date') %>% 
  group_by(Destination) %>%
  summarize(gamma_inferred_arrival_mean_date=mean(date))
colnames(gamma_inferred_arrival)<-c('country','gamma_inferred_arrival_mean_date')
gamma_inferred_arrival

delta_inferred_arrival<-d_delta_earliest ##### FROM FIG2 script
delta_inferred_arrival<-delta_inferred_arrival %>% 
  dplyr::select("Origin",'Destination','replicate','date') %>% 
  group_by(Destination) %>%
  summarize(delta_inferred_arrival_mean_date=mean(date))
colnames(delta_inferred_arrival)<-c('country','delta_inferred_arrival_mean_date')
delta_inferred_arrival


BA1_inferred_arrival<-d_BA1_earliest ##### FROM FIG2 script
BA1_inferred_arrival<-BA1_inferred_arrival %>% 
  dplyr::select("Origin",'Destination','replicate','date') %>% 
  group_by(Destination) %>%
  summarize(BA1_inferred_arrival_mean_date=mean(date))
colnames(BA1_inferred_arrival)<-c('country','BA1_inferred_arrival_mean_date')
BA1_inferred_arrival


BA2_inferred_arrival<-d_BA2_earliest ##### FROM FIG2 script
BA2_inferred_arrival<-BA2_inferred_arrival %>% 
  dplyr::select("Origin",'Destination','replicate','date') %>% 
  group_by(Destination) %>%
  summarize(BA2_inferred_arrival_mean_date=mean(date))
colnames(BA2_inferred_arrival)<-c('country','BA2_inferred_arrival_mean_date')
BA2_inferred_arrival


BA4_5_inferred_arrival<-d_BA4_5_earliest ##### FROM FIG2 script
BA4_5_inferred_arrival<-BA4_5_inferred_arrival %>% 
  dplyr::select("Origin",'Destination','replicate','date') %>% 
  group_by(Destination) %>%
  summarize(BA4_5_inferred_arrival_mean_date=mean(date))
colnames(BA4_5_inferred_arrival)<-c('country','BA4_5_inferred_arrival_mean_date')
BA4_5_inferred_arrival

all_inferred_arrivals<-plyr::join_all(list(alpha_inferred_arrival,beta_inferred_arrival,gamma_inferred_arrival,delta_inferred_arrival,BA1_inferred_arrival,BA2_inferred_arrival))


#Global summary of VOC arrival vs global travel

#V1: arrival=first detected sequenced case
p_VOC_arrival<-ggplot()+
  theme_bw()+
  
  geom_segment(aes(x=alpha_first_date_lo,xend=alpha_first_date_hi,y=-20,yend=-20),color='#d5bdaf',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=alpha_first_date,xend=alpha_first_date,y=-30,yend=-10),color='#d5bdaf',size=1,show.legend=FALSE)+
  geom_segment(aes(x=beta_first_date_lo,xend=beta_first_date_hi,y=-10,yend=-10),color='#ffd166',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=beta_first_date,xend=beta_first_date,y=-20,yend=0),color='#ffd166',size=1,show.legend=FALSE)+
  geom_segment(aes(x=delta_first_date_lo,xend=delta_first_date_hi,y=-10,yend=-10),color='#d74e09',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=delta_first_date,xend=delta_first_date,y=-20,yend=0),color='#d74e09',size=1,show.legend=FALSE)+
  geom_segment(aes(x=gamma_first_date_lo,xend=gamma_first_date_hi,y=-20,yend=-20),color='#2081c3',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=gamma_first_date,xend=gamma_first_date,y=-30,yend=-10),color='#2081c3',size=1,show.legend=FALSE)+
  geom_segment(aes(x=ba1_first_date_lo,xend=ba1_first_date_hi,y=-10,yend=-10),color='#a3b18a',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=ba1_first_date,xend=ba1_first_date,y=-20,yend=0),color='#a3b18a',size=1,show.legend=FALSE)+
  geom_segment(aes(x=ba2_first_date_lo,xend=ba2_first_date_hi,y=-20,yend=-20),color='#588157',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=ba2_first_date,xend=ba2_first_date,y=-30,yend=-10),color='#588157',size=1,show.legend=FALSE)+
  geom_segment(aes(x=ba4_5_first_date_lo,xend=ba4_5_first_date_hi,y=-10,yend=-10),color='#3a5a40',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=ba4_5_first_date,xend=ba4_5_first_date,y=-20,yend=0),color='#3a5a40',size=1,show.legend=FALSE)+
  
  geom_line(data = subset(aggregate_total,Date>"2020-07-01"), aes(x=Date, y=passengers/150000, color="Global Travel\n(Monthly Passengers)"),size=2)+
  
  geom_violin(width=20,data=subset(first_VOC_country_table, Alpha_arrival>0), aes(x=alpha_first_date,y=Alpha_arrival, fill='Alpha'), alpha=0.9)+
  geom_violin(width=20,data=subset(first_VOC_country_table, Beta_arrival>0), aes(x=beta_first_date,y=Beta_arrival,fill='Beta'), alpha=0.9)+
  geom_violin(width=20,data=subset(first_VOC_country_table, Delta_arrival>0), aes(x=delta_first_date,y=Delta_arrival,fill='Delta'), alpha=0.9)+
  geom_violin(width=20,data=subset(first_VOC_country_table, Gamma_arrival>0), aes(x=gamma_first_date,y=Gamma_arrival,fill='Gamma'), alpha=0.9)+
  geom_violin(width=20,data=subset(first_VOC_country_table, BA1_arrival>0), aes(x=ba1_first_date,y=BA1_arrival,fill='Omicron BA.1'), alpha=0.9)+
  geom_violin(width=20,data=subset(first_VOC_country_table, BA2_arrival>0), aes(x=ba2_first_date,y=BA2_arrival,fill='Omicron BA.2'), alpha=0.9)+
  geom_violin(width=20,data=subset(first_VOC_country_table, BA4_5_arrival>0), aes(x=ba4_5_first_date,y=BA4_5_arrival,fill='Omicron BA.4/5'), alpha=0.9)+
  
  scale_fill_manual(values=c('#d5bdaf',"#ffd166","#d74e09","#2081c3","#a3b18a","#588157","#3a5a40"),name='Variant')+
  scale_color_manual(values=c("grey50"), name='')+
  stat_summary(data=subset(first_VOC_country_table, Alpha_arrival>0), aes(x=alpha_first_date,y=Alpha_arrival),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(first_VOC_country_table, Alpha_arrival>0), aes(x=alpha_first_date,y=Alpha_arrival),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.4,colour='#d5bdaf')+
  
  stat_summary(data=subset(first_VOC_country_table, Beta_arrival>0), aes(x=beta_first_date,y=Beta_arrival),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(first_VOC_country_table, Beta_arrival>0), aes(x=beta_first_date,y=Beta_arrival),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.6,colour='#ffd166')+
  stat_summary(data=subset(first_VOC_country_table, Delta_arrival>0), aes(x=delta_first_date,y=Delta_arrival),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(first_VOC_country_table, Delta_arrival>0), aes(x=delta_first_date,y=Delta_arrival),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-1,colour='#d74e09')+
  stat_summary(data=subset(first_VOC_country_table, Gamma_arrival>0), aes(x=gamma_first_date,y=Gamma_arrival),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(first_VOC_country_table, Gamma_arrival>0), aes(x=gamma_first_date,y=Gamma_arrival),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.4,colour='#2081c3')+
  stat_summary(data=subset(first_VOC_country_table, BA1_arrival>0), aes(x=ba1_first_date,y=BA1_arrival),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(first_VOC_country_table, BA1_arrival>0), aes(x=ba1_first_date,y=BA1_arrival),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.4,hjust=1.3,colour='#a3b18a')+
  stat_summary(data=subset(first_VOC_country_table, BA2_arrival>0), aes(x=ba2_first_date,y=BA2_arrival),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(first_VOC_country_table, BA2_arrival>0), aes(x=ba2_first_date,y=BA2_arrival),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.4,hjust=1.2,colour='#588157')+
  stat_summary(data=subset(first_VOC_country_table, BA4_5_arrival>0), aes(x=ba4_5_first_date,y=BA4_5_arrival),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(first_VOC_country_table, BA4_5_arrival>0), aes(x=ba4_5_first_date,y=BA4_5_arrival),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.4,hjust=1.2,colour='#3a5a40')+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "4 months")+
  xlab('Date of VOC origin')+
  ylab('Delay in VOC arrival in countries worldwide\n(days between VOC origin and first sampling in n countries)')+
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    name = "Delay between VOC TMRCA and\nfirst sampling in n countries (No. of Days)",
    sec.axis = sec_axis( trans=~.*150000, name="Number of Passengers",  labels=label_number_si())) +
  ggtitle("Delay in VOC arrival in countries worldwide")+
  theme(legend.position='bottom',plot.title = element_text(hjust=0.5, family="Helvetica"))

p_VOC_arrival

#V2: arrival=first inferred introduction
p_VOC_arrival_inferred<-ggplot()+
  theme_bw()+
  geom_segment(aes(x=alpha_first_date_lo,xend=alpha_first_date_hi,y=-20,yend=-20),color='#d5bdaf',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=alpha_first_date,xend=alpha_first_date,y=-30,yend=-10),color='#d5bdaf',size=1,show.legend=FALSE)+
  geom_segment(aes(x=beta_first_date_lo,xend=beta_first_date_hi,y=-10,yend=-10),color='#ffd166',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=beta_first_date,xend=beta_first_date,y=-20,yend=0),color='#ffd166',size=1,show.legend=FALSE)+
  geom_segment(aes(x=delta_first_date_lo,xend=delta_first_date_hi,y=-10,yend=-10),color='#d74e09',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=delta_first_date,xend=delta_first_date,y=-20,yend=0),color='#d74e09',size=1,show.legend=FALSE)+
  geom_segment(aes(x=gamma_first_date_lo,xend=gamma_first_date_hi,y=-20,yend=-20),color='#2081c3',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=gamma_first_date,xend=gamma_first_date,y=-30,yend=-10),color='#2081c3',size=1,show.legend=FALSE)+
  geom_segment(aes(x=ba1_first_date_lo,xend=ba1_first_date_hi,y=-10,yend=-10),color='#a3b18a',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=ba1_first_date,xend=ba1_first_date,y=-20,yend=0),color='#a3b18a',size=1,show.legend=FALSE)+
  geom_segment(aes(x=ba2_first_date_lo,xend=ba2_first_date_hi,y=-20,yend=-20),color='#588157',size=1.5,show.legend=FALSE)+
  geom_segment(aes(x=ba2_first_date,xend=ba2_first_date,y=-30,yend=-10),color='#588157',size=1,show.legend=FALSE)+
  geom_line(data = subset(aggregate_total,Date>"2020-07-01"), aes(x=Date, y=passengers/150000, color="Global Travel\n(Monthly Passengers)"),size=2)+
  geom_violin(width=20,data=subset(VOC_arrival_sampling_inferred_country_table, alpha_inferred_arrival_mean_date-alpha_first_date>0), aes(x=alpha_first_date,y=alpha_inferred_arrival_mean_date-alpha_first_date, fill='Alpha'),alpha=0.9)+
  geom_violin(width=20,data=subset(VOC_arrival_sampling_inferred_country_table, beta_inferred_arrival_mean_date-beta_first_date>0), aes(x=beta_first_date,y=beta_inferred_arrival_mean_date-beta_first_date,fill='Beta'),alpha=0.9)+
  geom_violin(width=20,data=subset(VOC_arrival_sampling_inferred_country_table, gamma_inferred_arrival_mean_date-gamma_first_date>0), aes(x=gamma_first_date,y=gamma_inferred_arrival_mean_date-gamma_first_date,fill='Gamma'),alpha=0.9)+
  geom_violin(width=20,data=subset(VOC_arrival_sampling_inferred_country_table, delta_inferred_arrival_mean_date-delta_first_date>0), aes(x=delta_first_date,y=delta_inferred_arrival_mean_date-delta_first_date,fill='Delta'),alpha=0.9)+
  geom_violin(width=20,data=subset(VOC_arrival_sampling_inferred_country_table, BA1_inferred_arrival_mean_date-ba1_first_date>0), aes(x=ba1_first_date,y=BA1_inferred_arrival_mean_date-ba1_first_date,fill='Omicron BA.1'),alpha=0.9)+
  geom_violin(width=20,data=subset(VOC_arrival_sampling_inferred_country_table, BA2_inferred_arrival_mean_date-ba2_first_date>0), aes(x=ba2_first_date,y=BA2_inferred_arrival_mean_date-ba2_first_date,fill='Omicron BA.2'),alpha=0.9)+
  
  scale_fill_manual(values=c('#d5bdaf',"#ffd166","#d74e09","#2081c3","#a3b18a","#588157","#3a5a40"),name='Variant')+
  scale_color_manual(values=c("grey50"), name='')+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, alpha_inferred_arrival_mean_date-alpha_first_date>0), aes(x=alpha_first_date,y=alpha_inferred_arrival_mean_date-alpha_first_date),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, alpha_inferred_arrival_mean_date-alpha_first_date>0), aes(x=alpha_first_date,y=alpha_inferred_arrival_mean_date-alpha_first_date),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,hjust=0.8,colour='#d5bdaf')+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, beta_inferred_arrival_mean_date-beta_first_date>0), aes(x=beta_first_date,y=beta_inferred_arrival_mean_date-beta_first_date),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, beta_inferred_arrival_mean_date-beta_first_date>0), aes(x=beta_first_date,y=beta_inferred_arrival_mean_date-beta_first_date),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.6,colour='#ffd166')+
  
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, delta_inferred_arrival_mean_date-delta_first_date>0), aes(x=delta_first_date,y=delta_inferred_arrival_mean_date-delta_first_date),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, delta_inferred_arrival_mean_date-delta_first_date>0), aes(x=delta_first_date,y=delta_inferred_arrival_mean_date-delta_first_date),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=-0.6,colour='#d74e09')+
  
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, gamma_inferred_arrival_mean_date-gamma_first_date>0), aes(x=gamma_first_date,y=gamma_inferred_arrival_mean_date-gamma_first_date),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, gamma_inferred_arrival_mean_date-gamma_first_date>0), aes(x=gamma_first_date,y=gamma_inferred_arrival_mean_date-gamma_first_date),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=0,hjust=1.4,colour='#2081c3')+
  
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, BA1_inferred_arrival_mean_date-ba1_first_date>0), aes(x=ba1_first_date,y=BA1_inferred_arrival_mean_date-ba1_first_date),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, BA1_inferred_arrival_mean_date-ba1_first_date>0), aes(x=ba1_first_date,y=BA1_inferred_arrival_mean_date-ba1_first_date),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=0,hjust=1.8,colour='#a3b18a')+
  
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, BA2_inferred_arrival_mean_date-ba2_first_date>0), aes(x=ba2_first_date,y=BA2_inferred_arrival_mean_date-ba2_first_date),fun.data=median_hilow, fun.args = (conf.int=.5),
               geom="pointrange", color="black", size = 0.75,
               position = position_dodge(width = 0.9))+
  stat_summary(data=subset(VOC_arrival_sampling_inferred_country_table, BA2_inferred_arrival_mean_date-ba2_first_date>0), aes(x=ba2_first_date,y=BA2_inferred_arrival_mean_date-ba2_first_date),fun.data=stat_box_data,
               geom = "text",angle=90,
               size=3,vjust=0,hjust=1.8,colour='#588157')+
  xlab('Date of VOC origin')+
  ylab('Delay in VOC arrival in countries worldwide\n(days between VOC origin and inferred introduction)')+
  scale_x_date(date_labels = "%b\n%Y",date_breaks = "4 months")+
  scale_y_continuous(
    labels = label_number(accuracy = 1),
    name = "Days between VOC TMRCA\nand inferred introduction",
    sec.axis = sec_axis( trans=~.*150000, name="Number of Passengers",  labels=label_number_si())) +
  ggtitle("Delay in VOC arrival in countries worldwide")+
  theme(legend.position='bottom',plot.title = element_text(hjust=0.5, family="Helvetica"))



p_VOC_arrival_inferred


plot_grid(p_VOC_arrival,p_VOC_arrival_inferred)



#Alpha arrival vs travel

# Sampling dates vs global travel
first_VOC_country_table_alphag<-unique(subset(first_VOC_country_table,Alpha_arrival>0) %>% dplyr::select("country","Alpha",'Alpha_arrival'))

global_flights_specified_time<-subset(subset(total,year==2020), month>8) %>% #Global travel from September to December 2020
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_alphag<- first_VOC_country_table_alphag %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_alphag<-unique(subset(subset(first_VOC_country_table_alphag,!is.na(total)),Alpha_arrival>0))

first_VOC_country_table_alphag<- first_VOC_country_table_alphag %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Inferred dates vs global travel

first_VOC_country_table_alphag2<-unique(subset(VOC_arrival_sampling_inferred_country_table, alpha_inferred_arrival_mean_date-alpha_first_date>0) %>% dplyr::select("country","alpha_inferred_arrival_mean_date"))

global_flights_specified_time<-subset(subset(total,year==2020), month>8) %>% #Global travel from September to December 2020
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_alphag2<- first_VOC_country_table_alphag2 %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_alphag2<-unique(subset(subset(first_VOC_country_table_alphag2,!is.na(total)),alpha_inferred_arrival_mean_date-alpha_first_date>0))

first_VOC_country_table_alphag2<- first_VOC_country_table_alphag2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Sampling dates vs UK travel
first_VOC_country_table_alpha<-unique(subset(first_VOC_country_table,Alpha_arrival>0) %>% dplyr::select("country","Alpha",'Alpha_arrival'))


first_VOC_country_table_alpha<- first_VOC_country_table_alpha %>% 
  left_join(b117_UK_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_alpha<-unique(subset(subset(first_VOC_country_table_alpha,!is.na(TotalUKarrival)),Alpha_arrival>0))

first_VOC_country_table_alpha<- first_VOC_country_table_alpha %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Inferred dates vs UK travel

first_VOC_country_table_alpha2<-unique(subset(VOC_arrival_sampling_inferred_country_table, alpha_inferred_arrival_mean_date-alpha_first_date>0) %>% dplyr::select("country","alpha_inferred_arrival_mean_date"))


first_VOC_country_table_alpha2<- first_VOC_country_table_alpha2 %>% 
  left_join(b117_UK_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_alpha2<-unique(subset(subset(first_VOC_country_table_alpha2,!is.na(TotalUKarrival)),alpha_inferred_arrival_mean_date-alpha_first_date>0))

first_VOC_country_table_alpha2<- first_VOC_country_table_alpha2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Alpha correlation plot - sampling dates
corr_alpha_sampling_plot=ggplot()+
  theme_bw()+
  geom_smooth(colour='#d5bdaf',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_alpha,aes(as.numeric(Alpha_arrival),TotalUKarrival))+
  geom_point(colour='#d5bdaf',size=1,data=first_VOC_country_table_alpha,aes(as.numeric(Alpha_arrival),TotalUKarrival))+
  stat_cor(data=first_VOC_country_table_alpha,aes(as.numeric(Alpha_arrival),TotalUKarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d5bdaf',fontface=1,label.x = 1, label.y = 2.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_alpha,aes(as.numeric(Alpha_arrival),TotalUKarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d5bdaf',fontface=1,label.x = 1, label.y = 2.5,method = "spearman")+
  geom_smooth(colour='#d5bdaf',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_alphag,aes(as.numeric(Alpha_arrival),total))+
  geom_point(colour='#d5bdaf',size=1,shape=21,data=first_VOC_country_table_alphag,aes(as.numeric(Alpha_arrival),total))+
  stat_cor(data=first_VOC_country_table_alphag,aes(as.numeric(Alpha_arrival),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                        breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                        labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d5bdaf',fontface=1,label.x = 2, label.y = 7.5,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Alpha")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_alpha_sampling_plot

# Alpha correlation plot - inferred dates
corr_alpha_sampling_plot2=ggplot()+
  theme_bw()+
  geom_smooth(colour='#d5bdaf',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_alpha2,aes(as.numeric(alpha_inferred_arrival_mean_date-alpha_first_date),TotalUKarrival))+
  geom_point(colour='#d5bdaf',size=1,data=first_VOC_country_table_alpha2,aes(as.numeric(alpha_inferred_arrival_mean_date-alpha_first_date),TotalUKarrival))+
  stat_cor(data=first_VOC_country_table_alpha2,aes(as.numeric(alpha_inferred_arrival_mean_date-alpha_first_date),TotalUKarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                     breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                     labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d5bdaf',fontface=1,label.x = 1.7, label.y = 2.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_alpha2,aes(as.numeric(alpha_inferred_arrival_mean_date-alpha_first_date),TotalUKarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                     breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                     labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d5bdaf',fontface=1,label.x = 1.7, label.y = 2.5,method = "spearman")+
  geom_smooth(colour='#d5bdaf',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_alphag2,aes(as.numeric(alpha_inferred_arrival_mean_date-alpha_first_date),total))+
  geom_point(colour='#d5bdaf',size=1,shape=21,data=first_VOC_country_table_alphag2,aes(as.numeric(alpha_inferred_arrival_mean_date-alpha_first_date),total))+
  stat_cor(data=first_VOC_country_table_alphag2,aes(as.numeric(alpha_inferred_arrival_mean_date-alpha_first_date),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                             breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                             labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d5bdaf',fontface=1,label.x = 2, label.y = 7.2,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Alpha")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_alpha_sampling_plot2

plot_grid(corr_alpha_sampling_plot,corr_alpha_sampling_plot2)

#Beta


#Beta arrival vs travel

# Sampling dates vs global travel
first_VOC_country_table_betag<-unique(subset(first_VOC_country_table,Beta_arrival>0) %>% dplyr::select("country","Beta",'Beta_arrival'))

global_flights_specified_time<-rbind(subset(total,year==2020 & month>10), subset(total,year==2021 & month<2)) %>% #Global travel from November 2020 to January 2021
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_betag<- first_VOC_country_table_betag %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_betag<-unique(subset(subset(first_VOC_country_table_betag,!is.na(total)),Beta_arrival>0))

first_VOC_country_table_betag<- first_VOC_country_table_betag %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

corr_beta_sampling_global_p=cortest2df(cor.test(as.numeric(first_VOC_country_table_beta$Beta_arrival), first_VOC_country_table_beta$total, method = "pearson", use = "complete.obs"),"Beta","Sampled Arrival vs Global Travel","pearson")

# Inferred dates vs global travel

first_VOC_country_table_betag2<-unique(subset(VOC_arrival_sampling_inferred_country_table, beta_inferred_arrival_mean_date-beta_first_date>0) %>% dplyr::select("country","beta_inferred_arrival_mean_date"))

global_flights_specified_time<-rbind(subset(total,year==2020 & month>10), subset(total,year==2021 & month<2)) %>% #Global travel from November 2020 to January 2021
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_betag2<- first_VOC_country_table_betag2 %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_betag2<-unique(subset(subset(first_VOC_country_table_betag2,!is.na(total)),beta_inferred_arrival_mean_date-beta_first_date>0))

first_VOC_country_table_betag2<- first_VOC_country_table_betag2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Sampling dates vs SA travel

first_VOC_country_table_beta<-unique(subset(first_VOC_country_table,Beta_arrival>0) %>% dplyr::select("country","Beta",'Beta_arrival'))


first_VOC_country_table_beta<- first_VOC_country_table_beta %>% 
  left_join(b1351_SA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_beta<-unique(subset(subset(first_VOC_country_table_beta,!is.na(TotalSAarrival)),Beta_arrival>0))

first_VOC_country_table_beta<- first_VOC_country_table_beta %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

corr_beta_sampling_sa_p=cortest2df(cor.test(as.numeric(first_VOC_country_table_beta$Beta_arrival), first_VOC_country_table_beta$TotalSAarrival, method = "pearson", use = "complete.obs"),"Beta","Sampled Arrival vs Travel from First Reporting Country","pearson")

# Inferred dates vs SA travel

first_VOC_country_table_beta2<-unique(subset(VOC_arrival_sampling_inferred_country_table, beta_inferred_arrival_mean_date-beta_first_date>0) %>% dplyr::select("country","beta_inferred_arrival_mean_date"))


first_VOC_country_table_beta2<- first_VOC_country_table_beta2 %>% 
  left_join(b1351_SA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_beta2<-unique(subset(subset(first_VOC_country_table_beta2,!is.na(TotalSAarrival)),beta_inferred_arrival_mean_date-beta_first_date>0))

first_VOC_country_table_beta2<- first_VOC_country_table_beta2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

#Beta correlation plot (sampling dates)
corr_beta_sampling_plot=ggplot()+
  theme_bw()+
  geom_smooth(colour='#ffd166',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_beta,aes(as.numeric(Beta_arrival),TotalSAarrival))+
  geom_point(colour='#ffd166',size=1,data=first_VOC_country_table_beta,aes(as.numeric(Beta_arrival),TotalSAarrival))+
  stat_cor(data=first_VOC_country_table_beta,aes(as.numeric(Beta_arrival),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                              breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                              labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#ffd166',fontface=1,label.x = 1.5, label.y = 2.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_beta,aes(as.numeric(Beta_arrival),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                              breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                              labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#ffd166',fontface=1,label.x = 1.5, label.y = 2.5,method = "spearman")+
  geom_smooth(colour='#ffd166',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_betag,aes(as.numeric(Beta_arrival),total))+
  geom_point(colour='#ffd166',size=1,shape=21,data=first_VOC_country_table_betag,aes(as.numeric(Beta_arrival),total))+
  stat_cor(data=first_VOC_country_table_betag,aes(as.numeric(Beta_arrival),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                      labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#ffd166',fontface=1,label.x = 1.5, label.y = 7.5,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Beta")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_beta_sampling_plot

#Beta correlation plot (inferred dates)
corr_beta_sampling_plot2=ggplot()+
  theme_bw()+
  geom_smooth(colour='#ffd166',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_beta2,aes(as.numeric(beta_inferred_arrival_mean_date-beta_first_date),TotalSAarrival))+
  geom_point(colour='#ffd166',size=1,data=first_VOC_country_table_beta2,aes(as.numeric(beta_inferred_arrival_mean_date-beta_first_date),TotalSAarrival))+
  stat_cor(data=first_VOC_country_table_beta2,aes(as.numeric(beta_inferred_arrival_mean_date-beta_first_date),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                  labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#ffd166',fontface=1,label.x = 2, label.y = 2.2,method = "spearman")+
  stat_cor(data=first_VOC_country_table_beta2,aes(as.numeric(beta_inferred_arrival_mean_date-beta_first_date),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                  breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                  labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#ffd166',fontface=1,label.x = 2, label.y = 2.2,method = "spearman")+
  geom_smooth(colour='#ffd166',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_betag2,aes(as.numeric(beta_inferred_arrival_mean_date-beta_first_date),total))+
  geom_point(colour='#ffd166',size=1,shape=21,data=first_VOC_country_table_betag2,aes(as.numeric(beta_inferred_arrival_mean_date-beta_first_date),total))+
  stat_cor(data=first_VOC_country_table_betag2,aes(as.numeric(beta_inferred_arrival_mean_date-beta_first_date),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                          breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                          labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#ffd166',fontface=1,label.x = 2, label.y = 7.5,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Beta")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_beta_sampling_plot2

plot_grid(corr_beta_sampling_plot,corr_beta_sampling_plot2)


#Gamma arrival vs travel

# Sampling dates vs global travel
first_VOC_country_table_gammag<-unique(subset(first_VOC_country_table,Gamma_arrival>0) %>% dplyr::select("country","Gamma",'Gamma_arrival'))

global_flights_specified_time<-rbind(subset(total,year==2020 & month>11), subset(total,year==2021 & month<4)) %>% #Global travel from December 2020 to March 2021
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_gammag<- first_VOC_country_table_gammag %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_gammag<-unique(subset(subset(first_VOC_country_table_gammag,!is.na(total)),Gamma_arrival>0))

first_VOC_country_table_gammag<- first_VOC_country_table_gammag %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))


# Inferred dates vs global travel

first_VOC_country_table_gammag2<-unique(subset(VOC_arrival_sampling_inferred_country_table, gamma_inferred_arrival_mean_date-gamma_first_date>0) %>% dplyr::select("country","gamma_inferred_arrival_mean_date"))

global_flights_specified_time<-rbind(subset(total,year==2020 & month>11), subset(total,year==2021 & month<4)) %>% #Global travel from December 2020 to March 2021
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_gammag2<- first_VOC_country_table_gammag2 %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_gammag2<-unique(subset(subset(first_VOC_country_table_gammag2,!is.na(total)),gamma_inferred_arrival_mean_date-gamma_first_date>0))

first_VOC_country_table_gammag2<- first_VOC_country_table_gammag2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))


# Sampling dates vs Brazil travel

first_VOC_country_table_gamma<-unique(subset(first_VOC_country_table,Gamma_arrival>0) %>% dplyr::select("country","Gamma",'Gamma_arrival'))


first_VOC_country_table_gamma<- first_VOC_country_table_gamma %>% 
  left_join(P1_BRA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_gamma<-unique(subset(subset(first_VOC_country_table_gamma,!is.na(TotalBRAarrival)),Gamma_arrival>0))

first_VOC_country_table_gamma<- first_VOC_country_table_gamma %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))


# Inferred dates vs SA travel


first_VOC_country_table_gamma2<-unique(subset(VOC_arrival_sampling_inferred_country_table, gamma_inferred_arrival_mean_date-gamma_first_date>0) %>% dplyr::select("country","gamma_inferred_arrival_mean_date"))


first_VOC_country_table_gamma2<- first_VOC_country_table_gamma2 %>% 
  left_join(P1_BRA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_gamma2<-unique(subset(subset(first_VOC_country_table_gamma2,!is.na(TotalBRAarrival)),gamma_inferred_arrival_mean_date-gamma_first_date>0))

first_VOC_country_table_gamma2<- first_VOC_country_table_gamma2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))


#Gamma correlation plot (sampling dates)
corr_gamma_sampling_plot=ggplot()+
  theme_bw()+
  geom_smooth(colour='#2081c3',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_gamma,aes(as.numeric(Gamma_arrival),TotalBRAarrival))+
  geom_point(colour='#2081c3',size=1,data=first_VOC_country_table_gamma,aes(as.numeric(Gamma_arrival),TotalBRAarrival))+
  stat_cor(data=first_VOC_country_table_gamma,aes(as.numeric(Gamma_arrival),TotalBRAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                 breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                 labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#2081c3',fontface=1,label.x = 1, label.y = 2.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_gamma,aes(as.numeric(Gamma_arrival),TotalBRAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                 breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                 labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#2081c3',fontface=1,label.x = 1, label.y = 2.5,method = "spearman")+
  geom_smooth(colour='#2081c3',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_gammag,aes(as.numeric(Gamma_arrival),total))+
  geom_point(colour='#2081c3',size=1,shape=21,data=first_VOC_country_table_gammag,aes(as.numeric(Gamma_arrival),total))+
  stat_cor(data=first_VOC_country_table_gammag,aes(as.numeric(Gamma_arrival),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                        breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                        labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#2081c3',fontface=1,label.x = 2, label.y = 7.35,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Gamma")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_gamma_sampling_plot


#Gamma correlation plot (inferred dates)
corr_gamma_sampling_plot2=ggplot()+
  theme_bw()+
  geom_smooth(colour='#2081c3',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_gamma2,aes(as.numeric(gamma_inferred_arrival_mean_date-gamma_first_date),TotalBRAarrival))+
  geom_point(colour='#2081c3',size=1,data=first_VOC_country_table_gamma2,aes(as.numeric(gamma_inferred_arrival_mean_date-gamma_first_date),TotalBRAarrival))+
  stat_cor(data=first_VOC_country_table_gamma2,aes(as.numeric(gamma_inferred_arrival_mean_date-gamma_first_date),TotalBRAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                      labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#2081c3',fontface=1,label.x = 1.5, label.y = 3,method = "spearman")+
  stat_cor(data=first_VOC_country_table_gamma2,aes(as.numeric(gamma_inferred_arrival_mean_date-gamma_first_date),TotalBRAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                      labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#2081c3',fontface=1,label.x = 1.5, label.y = 3,method = "spearman")+
  geom_smooth(colour='#2081c3',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_gammag2,aes(as.numeric(gamma_inferred_arrival_mean_date-gamma_first_date),total))+
  geom_point(colour='#2081c3',size=1,shape=21,data=first_VOC_country_table_gammag2,aes(as.numeric(gamma_inferred_arrival_mean_date-gamma_first_date),total))+
  stat_cor(data=first_VOC_country_table_gammag2,aes(as.numeric(gamma_inferred_arrival_mean_date-gamma_first_date),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                             breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                             labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#2081c3',fontface=1,label.x = 1.8, label.y = 7.35,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Gamma")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_gamma_sampling_plot2

plot_grid(corr_gamma_sampling_plot,corr_gamma_sampling_plot2)

#Delta arrival vs travel

# Sampling dates vs global travel
first_VOC_country_table_deltag<-unique(subset(first_VOC_country_table,Delta_arrival>0) %>% dplyr::select("country","Delta",'Delta_arrival'))

global_flights_specified_time<-rbind(subset(total,year==2020 & month>11), subset(total,year==2021 & month<4)) %>% #Global travel from December 2020 to March 2021
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_deltag<- first_VOC_country_table_deltag %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_deltag<-unique(subset(subset(first_VOC_country_table_deltag,!is.na(total)),Delta_arrival>0))

first_VOC_country_table_deltag<- first_VOC_country_table_deltag %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))



# Inferred dates vs global travel

first_VOC_country_table_deltag2<-unique(subset(VOC_arrival_sampling_inferred_country_table, delta_inferred_arrival_mean_date-delta_first_date>0) %>% dplyr::select("country","delta_inferred_arrival_mean_date"))

global_flights_specified_time<-rbind(subset(total,year==2020 & month>11), subset(total,year==2021 & month<4)) %>% #Global travel from December 2020 to March 2021
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_deltag2<- first_VOC_country_table_deltag2 %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_deltag2<-unique(subset(subset(first_VOC_country_table_deltag2,!is.na(total)),delta_inferred_arrival_mean_date-delta_first_date>0))

first_VOC_country_table_deltag2<- first_VOC_country_table_deltag2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Sampling dates vs India travel


first_VOC_country_table_delta<-unique(subset(first_VOC_country_table,Delta_arrival>0) %>% dplyr::select("country","Delta",'Delta_arrival'))


first_VOC_country_table_delta<- first_VOC_country_table_delta %>% 
  left_join(Delta_IND_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_delta<-unique(subset(subset(first_VOC_country_table_delta,!is.na(TotalINDarrival)),Delta_arrival>0))

first_VOC_country_table_delta<- first_VOC_country_table_delta %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

corr_delta_sampling_india_p=cortest2df(cor.test(as.numeric(first_VOC_country_table_delta$Delta_arrival), first_VOC_country_table_delta$TotalINDarrival, method = "pearson", use = "complete.obs"),"Delta","Sampled Arrival vs Travel from First Reporting Country","pearson")


# Inferred dates vs India travel

first_VOC_country_table_delta2<-unique(subset(VOC_arrival_sampling_inferred_country_table, delta_inferred_arrival_mean_date-delta_first_date>0) %>% dplyr::select("country","delta_inferred_arrival_mean_date"))


first_VOC_country_table_delta2<- first_VOC_country_table_delta2 %>% 
  left_join(Delta_IND_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_delta2<-unique(subset(subset(first_VOC_country_table_delta2,!is.na(TotalINDarrival)),delta_inferred_arrival_mean_date-delta_first_date>0))

first_VOC_country_table_delta2<- first_VOC_country_table_delta2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

#Delta correlation plot (sampling dates)
corr_delta_sampling_plot=ggplot()+
  theme_bw()+
  geom_smooth(colour='#d74e09',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_delta,aes(as.numeric(Delta_arrival),TotalINDarrival))+
  geom_point(colour='#d74e09',size=1,data=first_VOC_country_table_delta,aes(as.numeric(Delta_arrival),TotalINDarrival))+
  stat_cor(data=first_VOC_country_table_delta,aes(as.numeric(Delta_arrival),TotalINDarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                 breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                 labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d74e09',fontface=1,label.x = 1, label.y = 2.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_delta,aes(as.numeric(Delta_arrival),TotalINDarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                 breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                 labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d74e09',fontface=1,label.x = 1, label.y = 2.5,method = "spearman")+
  
  geom_smooth(colour='#d74e09',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_deltag,aes(as.numeric(Delta_arrival),total))+
  geom_point(colour='#d74e09',size=1,shape=21,data=first_VOC_country_table_deltag,aes(as.numeric(Delta_arrival),total))+
  stat_cor(data=first_VOC_country_table_deltag,aes(as.numeric(Delta_arrival),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                        breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                        labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d74e09',label.x = 2, label.y = 7.5,fontface=1,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Delta")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_delta_sampling_plot


#Delta correlation plot (sampling dates)
corr_delta_sampling_plot2=ggplot()+
  theme_bw()+
  geom_smooth(colour='#d74e09',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_delta2,aes(as.numeric(delta_inferred_arrival_mean_date-delta_first_date),TotalINDarrival))+
  geom_point(colour='#d74e09',size=1,data=first_VOC_country_table_delta2,aes(as.numeric(delta_inferred_arrival_mean_date-delta_first_date),TotalINDarrival))+
  stat_cor(data=first_VOC_country_table_delta2,aes(as.numeric(delta_inferred_arrival_mean_date-delta_first_date),TotalINDarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                      labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d74e09',fontface=1,label.x = 1.5, label.y = 2.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_delta2,aes(as.numeric(delta_inferred_arrival_mean_date-delta_first_date),TotalINDarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                                      breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                                      labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d74e09',fontface=1,label.x = 1.5, label.y = 2.5,method = "spearman")+
  
  geom_smooth(colour='#d74e09',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_deltag2,aes(as.numeric(delta_inferred_arrival_mean_date-delta_first_date),total))+
  geom_point(colour='#d74e09',size=1,shape=21,data=first_VOC_country_table_deltag2,aes(as.numeric(delta_inferred_arrival_mean_date-delta_first_date),total))+
  stat_cor(data=first_VOC_country_table_deltag2,aes(as.numeric(delta_inferred_arrival_mean_date-delta_first_date),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                             breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                             labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#d74e09',label.x = 2, label.y = 7.5,fontface=1,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Delta")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_delta_sampling_plot2

plot_grid(corr_delta_sampling_plot,corr_delta_sampling_plot2)


#BA.1 arrival vs travel

# Sampling dates vs global travel
first_VOC_country_table_BA1g<-unique(subset(first_VOC_country_table,BA1_arrival>0) %>% dplyr::select("country","BA1",'BA1_arrival'))

global_flights_specified_time<-rbind(subset(total,year==2021 & month>10), subset(total,year==2022 & month<3)) %>% #Global travel from November 2021 to February 2022
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_BA1g<- first_VOC_country_table_BA1g %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_BA1g<-unique(subset(subset(first_VOC_country_table_BA1g,!is.na(total)),BA1_arrival>0))

first_VOC_country_table_BA1g<- first_VOC_country_table_BA1g %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Inferred dates vs global travel

first_VOC_country_table_BA1g2<-unique(subset(VOC_arrival_sampling_inferred_country_table, BA1_inferred_arrival_mean_date-ba1_first_date>0) %>% dplyr::select("country","BA1_inferred_arrival_mean_date"))

global_flights_specified_time<-rbind(subset(total,year==2021 & month>10), subset(total,year==2022 & month<3)) %>% 
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_BA1g2<- first_VOC_country_table_BA1g2 %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_BA1g2<-unique(subset(subset(first_VOC_country_table_BA1g2,!is.na(total)),BA1_inferred_arrival_mean_date-ba1_first_date>0))

first_VOC_country_table_BA1g2<- first_VOC_country_table_BA1g2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))


# Sampling dates vs SA travel

first_VOC_country_table_BA1<-unique(subset(VOC_arrival_sampling_inferred_country_table,BA1_arrival>0) %>% dplyr::select("country","BA1",'BA1_arrival'))


first_VOC_country_table_BA1<- first_VOC_country_table_BA1 %>% 
  left_join(Omicron_SA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_BA1<-unique(subset(subset(first_VOC_country_table_BA1,!is.na(TotalSAarrival)),BA1_arrival>0))

first_VOC_country_table_BA1<- first_VOC_country_table_BA1 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

corr_ba1_sampling_sa_p=cortest2df(cor.test(as.numeric(first_VOC_country_table_BA1$BA1_arrival), first_VOC_country_table_BA1$TotalSAarrival, method = "pearson", use = "complete.obs"),"Omicron BA.1","Sampled Arrival vs Travel from First Reporting Country","pearson")


# Inferred dates vs SA travel

first_VOC_country_table_BA12<-unique(subset(VOC_arrival_sampling_inferred_country_table, BA1_inferred_arrival_mean_date-ba1_first_date>0) %>% dplyr::select("country","BA1_inferred_arrival_mean_date"))


first_VOC_country_table_BA12<- first_VOC_country_table_BA12 %>% 
  left_join(Omicron_SA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_BA12<-unique(subset(subset(first_VOC_country_table_BA12,!is.na(TotalSAarrival)),BA1_inferred_arrival_mean_date-ba1_first_date>0))

first_VOC_country_table_BA12<- first_VOC_country_table_BA12 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

#Omicron BA.1 correlation plot (sampling dates)
corr_ba1_sampling_plot=ggplot()+
  theme_bw()+
  geom_smooth(colour='#a3b18a',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_BA1,aes(as.numeric(BA1_arrival),TotalSAarrival))+
  geom_point(colour='#a3b18a',size=1,data=first_VOC_country_table_BA1,aes(as.numeric(BA1_arrival),TotalSAarrival))+
  stat_cor(data=first_VOC_country_table_BA1,aes(as.numeric(BA1_arrival),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                            breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#a3b18a',fontface=1,label.x = 0.5, label.y = 1.6,method = "spearman")+
  stat_cor(data=first_VOC_country_table_BA1,aes(as.numeric(BA1_arrival),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                            breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#a3b18a',fontface=1,label.x = 0.5, label.y = 1.6,method = "spearman")+
  
  geom_smooth(colour='#a3b18a',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_BA1g,aes(as.numeric(BA1_arrival),total))+
  geom_point(colour='#a3b18a',size=1,shape=21,data=first_VOC_country_table_BA1g,aes(as.numeric(BA1_arrival),total))+
  stat_cor(data=first_VOC_country_table_BA1g,aes(as.numeric(BA1_arrival),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                    labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#a3b18a',fontface=1,label.x = 1.8, label.y = 7.1,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Omicron BA.1")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_ba1_sampling_plot


#Omicron BA.1 correlation plot (inferred dates)
corr_ba1_sampling_plot2=ggplot()+
  theme_bw()+
  geom_smooth(colour='#a3b18a',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_BA12,aes(as.numeric(BA1_inferred_arrival_mean_date-ba1_first_date),TotalSAarrival))+
  geom_point(colour='#a3b18a',size=1,data=first_VOC_country_table_BA12,aes(as.numeric(BA1_inferred_arrival_mean_date-ba1_first_date),TotalSAarrival))+
  stat_cor(data=first_VOC_country_table_BA12,aes(as.numeric(BA1_inferred_arrival_mean_date-ba1_first_date),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                               breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                               labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#a3b18a',fontface=1,label.x = 1.5, label.y = 1.6,method = "spearman")+
  stat_cor(data=first_VOC_country_table_BA12,aes(as.numeric(BA1_inferred_arrival_mean_date-ba1_first_date),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                               breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                               labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#a3b18a',fontface=1,label.x = 1.5, label.y = 1.6,method = "spearman")+
  
  geom_smooth(colour='#a3b18a',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_BA1g2,aes(as.numeric(BA1_inferred_arrival_mean_date-ba1_first_date),total))+
  geom_point(colour='#a3b18a',size=1,shape=21,data=first_VOC_country_table_BA1g2,aes(as.numeric(BA1_inferred_arrival_mean_date-ba1_first_date),total))+
  stat_cor(data=first_VOC_country_table_BA1g2,aes(as.numeric(BA1_inferred_arrival_mean_date-ba1_first_date),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                       breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                       labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#a3b18a',fontface=1,label.x = 1.8, label.y = 7.1,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Omicron BA.1")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_ba1_sampling_plot2

plot_grid(corr_ba1_sampling_plot,corr_ba1_sampling_plot2)

#BA.2 arrival vs travel

# Sampling dates vs global travel
first_VOC_country_table_BA2g<-unique(subset(first_VOC_country_table,BA2_arrival>0) %>% dplyr::select("country","BA2",'BA2_arrival'))

global_flights_specified_time<-rbind(subset(total,year==2021 & month>10), subset(total,year==2022 & month<4)) %>% #Global travel from November 2021 to March 2022
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_BA2g<- first_VOC_country_table_BA2g %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_BA2g<-unique(subset(subset(first_VOC_country_table_BA2g,!is.na(total)),BA2_arrival>0))

first_VOC_country_table_BA2g<- first_VOC_country_table_BA2g %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Inferred dates vs global travel

first_VOC_country_table_BA2g2<-unique(subset(VOC_arrival_sampling_inferred_country_table, BA2_inferred_arrival_mean_date-ba2_first_date>0) %>% dplyr::select("country","BA2_inferred_arrival_mean_date"))

global_flights_specified_time<-rbind(subset(total,year==2021 & month>10), subset(total,year==2022 & month<4)) %>% 
  group_by(origCtryName) %>%
  summarize(total=sum(totalVol))

first_VOC_country_table_BA2g2<- first_VOC_country_table_BA2g2 %>% 
  left_join(global_flights_specified_time, by = c("country" = "origCtryName"))

first_VOC_country_table_BA2g2<-unique(subset(subset(first_VOC_country_table_BA2g2,!is.na(total)),BA2_inferred_arrival_mean_date-ba2_first_date>0))

first_VOC_country_table_BA2g2<- first_VOC_country_table_BA2g2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

# Sampling dates vs SA travel

first_VOC_country_table_BA2<-unique(subset(VOC_arrival_sampling_inferred_country_table,BA2_arrival>0) %>% dplyr::select("country","BA2",'BA2_arrival'))


first_VOC_country_table_BA2<- first_VOC_country_table_BA2 %>% 
  left_join(Omicron_SA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_BA2<-unique(subset(subset(first_VOC_country_table_BA2,!is.na(TotalSAarrival)),BA2_arrival>0))

first_VOC_country_table_BA2<- first_VOC_country_table_BA2 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))

corr_ba2_sampling_sa_p=cortest2df(cor.test(as.numeric(first_VOC_country_table_BA2$BA2_arrival), first_VOC_country_table_BA2$TotalSAarrival,  method = "pearson", use = "complete.obs"),"Omicron BA.2","Sampled Arrival vs Travel from First Reporting Country","pearson")


# Inferred dates vs SA travel

first_VOC_country_table_BA22<-unique(subset(VOC_arrival_sampling_inferred_country_table, BA2_inferred_arrival_mean_date-ba2_first_date>0) %>% dplyr::select("country","BA2_inferred_arrival_mean_date"))


first_VOC_country_table_BA22<- first_VOC_country_table_BA22 %>% 
  left_join(Omicron_SA_travel_data, by = c("country" = "destCtryName"))

first_VOC_country_table_BA22<-unique(subset(subset(first_VOC_country_table_BA22,!is.na(TotalSAarrival)),BA2_inferred_arrival_mean_date-ba2_first_date>0))

first_VOC_country_table_BA22<- first_VOC_country_table_BA22 %>% 
  left_join(world_map_data_short, by = c("country" = "NAME"))


#Omicron BA.2 correlation plot (sampling dates)
corr_ba2_sampling_plot=ggplot()+
  theme_bw()+
  geom_smooth(colour='#588157',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_BA2,aes(as.numeric(BA2_arrival),TotalSAarrival))+
  geom_point(colour='#588157',size=1,data=first_VOC_country_table_BA2,aes(as.numeric(BA2_arrival),TotalSAarrival))+
  stat_cor(data=first_VOC_country_table_BA2,aes(as.numeric(BA2_arrival),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                            breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#588157',fontface=1,label.x = 1, label.y = 1.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_BA2,aes(as.numeric(BA2_arrival),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                            breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                            labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#588157',fontface=1,label.x = 1, label.y = 1.5,method = "spearman")+
  
  geom_smooth(colour='#588157',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_BA2g,aes(as.numeric(BA2_arrival),total))+
  geom_point(colour='#588157',size=1,shape=21,data=first_VOC_country_table_BA2g,aes(as.numeric(BA2_arrival),total))+
  stat_cor(data=first_VOC_country_table_BA2g,aes(as.numeric(BA2_arrival),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                    breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                    labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#588157',fontface=1,label.x = 1.9, label.y = 7.1,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Omicron BA.2")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_ba2_sampling_plot


#Omicron BA.2 correlation plot (inferred dates)
corr_ba2_sampling_plot2=ggplot()+
  theme_bw()+
  geom_smooth(colour='#588157',method='lm',size=0.5,se=FALSE,data=first_VOC_country_table_BA22,aes(as.numeric(BA2_inferred_arrival_mean_date-ba2_first_date),TotalSAarrival))+
  geom_point(colour='#588157',size=1,data=first_VOC_country_table_BA22,aes(as.numeric(BA2_inferred_arrival_mean_date-ba2_first_date),TotalSAarrival))+
  stat_cor(data=first_VOC_country_table_BA22,aes(as.numeric(BA2_inferred_arrival_mean_date-ba2_first_date),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                               breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                               labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#588157',fontface=1,label.x = 1.5, label.y = 1.5,method = "spearman")+
  stat_cor(data=first_VOC_country_table_BA22,aes(as.numeric(BA2_inferred_arrival_mean_date-ba2_first_date),TotalSAarrival,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                               breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                               labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#588157',fontface=1,label.x = 1.5, label.y = 1.5,method = "spearman")+
  
  geom_smooth(colour='#588157',linetype='dashed',size=0.5,se=FALSE,method='lm',data=first_VOC_country_table_BA2g2,aes(as.numeric(BA2_inferred_arrival_mean_date-ba2_first_date),total))+
  geom_point(colour='#588157',size=1,shape=21,data=first_VOC_country_table_BA2g2,aes(as.numeric(BA2_inferred_arrival_mean_date-ba2_first_date),total))+
  stat_cor(data=first_VOC_country_table_BA2g2,aes(as.numeric(BA2_inferred_arrival_mean_date-ba2_first_date),total,label = paste("rho", "'='",..r..,cut(..p.., 
                                                                                                                                                       breaks = c(-Inf, 0.0001, 0.001, 0.01, 0.05, Inf),
                                                                                                                                                       labels = c("'****'", "'***'", "'**'", "'*'", "'ns'")), sep = "~")),size=3,colour='#588157',fontface=1,label.x = 1.7, label.y = 7.1,method = "spearman")+
  scale_y_continuous(trans='log10',label=label_number_si())+
  scale_x_continuous(trans='log10')+
  theme(axis.title = element_blank(),axis.text = element_text(size=8))+
  ggtitle("Omicron BA.2")+
  theme(plot.title = element_text(size=10,hjust=0.5))
corr_ba2_sampling_plot2

plot_grid(corr_ba2_sampling_plot,corr_ba2_sampling_plot2)



fig3b_new<-plot_grid(corr_alpha_sampling_plot,corr_beta_sampling_plot,corr_gamma_sampling_plot,corr_delta_sampling_plot,corr_ba1_sampling_plot,corr_ba2_sampling_plot,corr_ba4_5_sampling_plot,ncol=2)
plot_grid(p_VOC_arrival,fig3b_new,ncol=2,rel_widths = c(0.5,0.5),labels=c("A","B"))

suppfigb_new<-plot_grid(corr_alpha_sampling_plot2,corr_beta_sampling_plot2,corr_gamma_sampling_plot2,corr_delta_sampling_plot2,corr_ba1_sampling_plot2,corr_ba2_sampling_plot2,ncol=2)
plot_grid(p_VOC_arrival_inferred,suppfigb_new,ncol=2,rel_widths = c(0.5,0.5),labels=c("A","B"))

