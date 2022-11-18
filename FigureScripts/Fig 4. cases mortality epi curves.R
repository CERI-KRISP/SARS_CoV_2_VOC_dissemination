library(ggplot2)
library(tidyverse)
library(gridExtra)
library(lubridate)
library(data.table)
library(readxl)
library(plyr)
library(dplyr)
library(writexl)
library(cowplot)
library(ggpubr)
library(zoo)

####import GISAID dataset and cleaning
options(scipen = 100)
data <-fread('metadata.tsv', data.table=FALSE)
names(data)[names(data) == 'Collection date'] <- 'Cdate'
data <- subset(data, select = c('Accession ID', 'Cdate', 'Location', 'Pango lineage', 'Variant'))
data <- filter(data, data$Cdate != 2020)
data <- filter(data, data$Cdate != 2021)
data <- filter(data, data$Cdate != 2022)
data[c('Continent', 'Country', 'State')] <- str_split_fixed(data$Location, ' /', 3)
data$Country <- as.factor(data$Country)
data <- filter(data, data$Country != "")

data$Variants2 <- with(data, factor(Variant, 
                                    levels = c('VOC Alpha GRY (B.1.1.7+Q.*) first detected in the UK', 'VOC Beta GH/501Y.V2 (B.1.351+B.1.351.2+B.1.351.3) first detected in South Africa', 'VOC Gamma GR/501Y.V3 (P.1+P.1.*) first detected in Brazil/Japan','VOC Delta GK (B.1.617.2+AY.*) first detected in India', 'VOC Omicron GRA (B.1.1.529+BA.*) first detected in Botswana/Hong Kong/South Africa'), 
                                    labels = c("Alpha","Beta", "Gamma", "Delta", "Omicron")))

data$Variants2[is.na(data$Variants2)] = "Other"
names(data)[names(data) == 'Pango lineage'] <- 'Pango'

data[data$Pango %like% 'BA.1' & data$Variants2 == "Omicron", "Variants2"] <- "BA.1"
data[data$Pango %like% 'BA.2' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BA.4' & data$Variants2 == "Omicron", "Variants2"] <- "BA.4"
data[data$Pango %like% 'BA.5' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BA.3' & data$Variants2 == "Omicron", "Variants2"] <- "Omicron"
data[data$Pango %like% 'BC.' & data$Variants2 == "Omicron", "Variants2"] <- "BA.1"
data[data$Pango %like% 'BD' & data$Variants2 == "Omicron", "Variants2"] <- "BA.1"
data[data$Pango %like% 'BE' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BF' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BG' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BH' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BJ' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BK' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BL' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BM' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BN' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BP' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BQ' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BR' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BS' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BT' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BU' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BV' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BW' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'BY' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'BZ' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'CA' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'CB' & data$Variants2 == "Omicron", "Variants2"] <- "BA.2"
data[data$Pango %like% 'CC' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'CD' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'CE' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
data[data$Pango %like% 'CF' & data$Variants2 == "Omicron", "Variants2"] <- "BA.5"
levels(data$Variants2)
table(data$Variants2)
data$Variants2[data$Variants2=="Omicron"]<-"Other Omicron"

data$days<-as.Date(cut(data$Cdate,breaks = "day",start.on.monday = FALSE))

africa <- subset(data, Continent == 'Africa', select = c("Accession ID", "Cdate", "Location", "Variant", "Variants2","days", "Continent", "Country", "State"))

asia <- subset(data, Continent == 'Asia', select = c("Accession ID", "Cdate", "Location", "Variant", "Variants2","days", "Continent", "Country", "State"))
asia = subset(asia, asia$Cdate > '2020-01-01')

europe <- subset(data, Continent == 'Europe', select = c("Accession ID", "Cdate", "Location", "Variant", "Variants2","days", "Continent", "Country", "State"))

oceania <- subset(data, Continent == 'Oceania', select = c("Accession ID", "Cdate", "Location", "Variant", "Variants2","days", "Continent", "Country", "State"))

Namerica <- subset(data, Continent == 'North America', select = c("Accession ID", "Cdate", "Location", "Variant","days", "Variants2", "Continent", "Country", "State"))

Samerica <- subset(data, Continent == 'South America', select = c("Accession ID", "Cdate", "Location", "Variant","days", "Variants2", "Continent", "Country", "State"))


###OWID data - subset to continent

owid<-read_excel('owid.xlsx')

asia_owid <- subset(owid, location == 'Asia')
asia_owid$days<-as.Date(cut(asia_owid$date,breaks = "day",start.on.monday = FALSE))

africa_owid <- subset(owid, location == 'Africa')
africa_owid$days<-as.Date(cut(africa_owid$date,breaks = "day",start.on.monday = FALSE))

europe_owid <- subset(owid, location == 'Europe')
europe_owid$days<-as.Date(cut(europe_owid$date,breaks = "day",start.on.monday = FALSE))

oceania_owid <- subset(owid, location == 'Oceania')
oceania_owid$days<-as.Date(cut(oceania_owid$date,breaks = "day",start.on.monday = FALSE))

Namerica_owid <- subset(owid, location == 'North America')
Namerica_owid$days<-as.Date(cut(Namerica_owid$date,breaks = "day",start.on.monday = FALSE))

Samerica_owid <- subset(owid, location == 'South America')
Samerica_owid$days<-as.Date(cut(Samerica_owid$date,breaks = "day",start.on.monday = FALSE))


####Testing data 

africa_test <- subset (owid, continent=='Africa', select = c("date", 'new_tests_smoothed_per_thousand'))
africa_test$days<-as.Date(cut(africa_test$date,breaks = "day",start.on.monday = FALSE))
africa_new_ave <- aggregate(new_tests_smoothed_per_thousand ~ days, africa_test, mean)

asia_test <- subset (owid, continent=='Asia', select = c("date", 'new_tests_smoothed_per_thousand'))
asia_test$days<-as.Date(cut(asia_test$date,breaks = "day",start.on.monday = FALSE))
asia_new_ave <- aggregate(new_tests_smoothed_per_thousand ~ days, asia_test, mean)

europe_test <- subset (owid, continent=='Europe', select = c("date", 'new_tests_smoothed_per_thousand'))
europe_test$days<-as.Date(cut(europe_test$date,breaks = "day",start.on.monday = FALSE))
europe_new_ave <- aggregate(new_tests_smoothed_per_thousand ~ days, europe_test, mean)

oceania_test <- subset (owid, continent=='Oceania', select = c("date", 'new_tests_smoothed_per_thousand'))
oceania_test$days<-as.Date(cut(oceania_test$date,breaks = "day",start.on.monday = FALSE))
oceania_new_ave <- aggregate(new_tests_smoothed_per_thousand ~ days, oceania_test, mean)

Samerica_test <- subset (owid, continent=='South America', select = c("date", 'new_tests_smoothed_per_thousand'))
Samerica_test$days<-as.Date(cut(Samerica_test$date,breaks = "day",start.on.monday = FALSE))
Samerica_new_ave <- aggregate(new_tests_smoothed_per_thousand ~ days, Samerica_test, mean)

Namerica_test <- subset (owid, continent=='North America', select = c("date", 'new_tests_smoothed_per_thousand'))
Namerica_test$days<-as.Date(cut(Namerica_test$date,breaks = "day",start.on.monday = FALSE))
Namerica_new_ave <- aggregate(new_tests_smoothed_per_thousand ~ days, Namerica_test, mean)

###case per variant plots
###Asia
prop.table(table(asia$days, asia$Variants2))

P_asia <- prop.table(table(asia$days, asia$Variants2), margin=1)

temp_asia<-as.data.frame(P_asia)
names(temp_asia)[1] <- 'days'

head(temp_asia)
temp_asia$days<-as.Date(cut(as.Date(temp_asia$days,format='%Y-%m-%d'),
                            breaks = "day",
                            start.on.monday = FALSE))

temp2_asia<-asia_owid[c("days","new_cases_per_million")]

head(temp2_asia)

temp3_asia<-join(temp_asia, temp2_asia,
                 type = "left")

tail(temp3_asia)

temp3_asia$days<-as.Date(cut(as.Date(temp3_asia$days,format='%Y-%m-%d'),
                             breaks = "day",
                             start.on.monday = FALSE))

temp3_asia$cases_per_variant=temp3_asia$new_cases_per_million*temp3_asia$Freq

###smoothing
asia_new_ave2 <- asia_new_ave %>%
  dplyr::mutate(new_tests_smoothed = zoo::rollmean(new_tests_smoothed_per_thousand, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(asia_new_ave2)

temp3_asia <- temp3_asia %>%
  group_by(Var2) %>%
  dplyr::mutate(cases_per_variant_7day = zoo::rollmean(cases_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_asia)

dateVec <- seq(from = as.Date("2020-02-01"), to = as.Date("2022-10-01"), by = "days")

p_Epi_asia<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_asia, aes(x = days, y = cases_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+ 
  geom_line(data=asia_new_ave2, aes(x=days, y=new_tests_smoothed*120),linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 2000),breaks = c(20, 100, 200, 500, 1000, 1500, 2000),labels = c(20, 100,200, 500, 1000, 1500, 2000),
    sec.axis = sec_axis(~./120, name=" ",breaks = c(1, 2.5, 5, 10, 15),labels = c(1, 2.5, 5, 10, 15))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Asia')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

p_Epi_asia 


###Africa
prop.table(table(africa$days, africa$Variants2))

P_africa <- prop.table(table(africa$days, africa$Variants2), margin=1)

temp_africa<-as.data.frame(P_africa)
names(temp_africa)[1] <- 'days'

###smoothing for genomic data
temp_africa <- temp_africa %>%
  group_by(Var2) %>%
  dplyr::mutate(Freq_7day = zoo::rollmean(Freq, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp_africa)

temp_africa$days<-as.Date(cut(as.Date(temp_africa$days,format='%Y-%m-%d'),
                              breaks = "day",
                              start.on.monday = FALSE))

temp2_africa<-africa_owid[c("days","new_cases_per_million")]
head(temp2_africa)

temp3_africa<-join(temp_africa, temp2_africa,
                   type = "left")

tail(temp3_africa)

temp3_africa$days<-as.Date(cut(as.Date(temp3_africa$days,format='%Y-%m-%d'),
                               breaks = "day",
                               start.on.monday = FALSE))

temp3_africa$cases_per_variant=temp3_africa$new_cases_per_million*temp3_africa$Freq_7day

###smoothing for cases
temp3_africa <- temp3_africa %>%
  group_by(Var2) %>%
  dplyr::mutate(cases_per_variant_7day = zoo::rollmean(cases_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_africa)

africa_new_ave2 <- africa_new_ave %>%
  dplyr::mutate(new_tests_smoothed = zoo::rollmean(new_tests_smoothed_per_thousand, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(africa_new_ave2)

p_Epi_africa<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_africa, aes(x = days, y = cases_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+ 
  geom_line(data=africa_new_ave2, aes(x=days, y=new_tests_smoothed*120),linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 2000),breaks = c(20, 100, 200, 500, 1000, 1500, 2000),labels = c(20, 100,200, 500, 1000, 1500, 2000),
    sec.axis = sec_axis(~./120, name=" ",breaks = c(1, 2.5, 5, 10, 15),labels = c(1, 2.5, 5, 10, 15))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Africa')+
  theme(axis.title.y = element_text(hjust = 0.5, size=11, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

p_Epi_africa

###Europe
P_europe <- prop.table(table(europe$days, europe$Variants2), margin=1)

temp_europe<-as.data.frame(P_europe)
names(temp_europe)[1] <- 'days'

head(temp_europe)
temp_europe$days<-as.Date(cut(as.Date(temp_europe$days,format='%Y-%m-%d'),
                              breaks = "day",
                              start.on.monday = FALSE))

temp2_europe<-europe_owid[c("days","new_cases_per_million")]
head(temp2_europe)


temp3_europe<-join(temp_europe, temp2_europe,
                   type = "left")

tail(temp3_europe)

temp3_europe$days<-as.Date(cut(as.Date(temp3_europe$days,format='%Y-%m-%d'),
                               breaks = "day",
                               start.on.monday = FALSE))

temp3_europe$cases_per_variant=temp3_europe$new_cases_per_million*temp3_europe$Freq

###smoothing
temp3_europe <- temp3_europe %>%
  group_by(Var2) %>%
  dplyr::mutate(cases_per_variant_7day = zoo::rollmean(cases_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_europe)

europe_new_ave2 <- europe_new_ave %>%
  dplyr::mutate(new_tests_smoothed = zoo::rollmean(new_tests_smoothed_per_thousand, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(europe_new_ave2)

p_Epi_europe<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_europe, aes(x = days, y = cases_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+ 
  geom_line(data=europe_new_ave2, aes(x=days, y=new_tests_smoothed*120),linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 2000),breaks = c(20, 100, 200, 500, 1000, 1500, 2000),labels = c(20, 100,200, 500, 1000, 1500, 2000),
    sec.axis = sec_axis(~./120, name=" ",breaks = c(1, 2.5, 5, 10, 15),labels = c(1, 2.5, 5, 10, 15))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(axis.title.y = element_text(size = 8.5))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Europe')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

p_Epi_europe

###Oceania
P_oceania <- prop.table(table(oceania$days, oceania$Variants2), margin=1)

temp_oceania<-as.data.frame(P_oceania)
names(temp_oceania)[1] <- 'days'

head(temp_oceania)
temp_oceania$days<-as.Date(cut(as.Date(temp_oceania$days,format='%Y-%m-%d'),
                               breaks = "day",
                               start.on.monday = FALSE))

temp2_oceania<-oceania_owid[c("days","new_cases_per_million")]
head(temp2_oceania)

temp3_oceania<-join(temp_oceania, temp2_oceania,
                    type = "left")

tail(temp3_oceania)

temp3_oceania$days<-as.Date(cut(as.Date(temp3_oceania$days,format='%Y-%m-%d'),
                                breaks = "day",
                                start.on.monday = FALSE))

temp3_oceania$cases_per_variant=temp3_oceania$new_cases_per_million*temp3_oceania$Freq

###smoothing
temp3_oceania <- temp3_oceania %>%
  group_by(Var2) %>%
  dplyr::mutate(cases_per_variant_7day = zoo::rollmean(cases_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_oceania)

oceania_new_ave2 <- oceania_new_ave %>%
  dplyr::mutate(new_tests_smoothed = zoo::rollmean(new_tests_smoothed_per_thousand, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(oceania_new_ave2)

p_Epi_oceania<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_oceania, aes(x = days, y = cases_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+ 
  geom_line(data=oceania_new_ave2, aes(x=days, y=new_tests_smoothed*120),linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 2000),breaks = c(20, 100, 200, 500, 1000, 1500, 2000),labels = c(20, 100,200, 500, 1000, 1500, 2000),
    sec.axis = sec_axis(~./120, name=" ",breaks = c(1, 2.5, 5, 10, 15),labels = c(1, 2.5, 5, 10, 15))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Oceania')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

p_Epi_oceania

###Namerica
P_Namerica <- prop.table(table(Namerica$days, Namerica$Variants2), margin=1)

temp_Namerica<-as.data.frame(P_Namerica)
names(temp_Namerica)[1] <- 'days'

head(temp_Namerica)
temp_Namerica$days<-as.Date(cut(as.Date(temp_Namerica$days,format='%Y-%m-%d'),
                                breaks = "day",
                                start.on.monday = FALSE))

temp2_Namerica<-Namerica_owid[c("days","new_cases_per_million")]
head(temp2_Namerica)

temp3_Namerica<-join(temp_Namerica, temp2_Namerica,
                     type = "left")

tail(temp3_Namerica)

temp3_Namerica$days<-as.Date(cut(as.Date(temp3_Namerica$days,format='%Y-%m-%d'),
                                 breaks = "day",
                                 start.on.monday = FALSE))

temp3_Namerica$cases_per_variant=temp3_Namerica$new_cases_per_million*temp3_Namerica$Freq

###smoothing
temp3_Namerica <- temp3_Namerica %>%
  group_by(Var2) %>%
  dplyr::mutate(cases_per_variant_7day = zoo::rollmean(cases_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_Namerica)

Namerica_new_ave2 <- Namerica_new_ave %>%
  dplyr::mutate(new_tests_smoothed = zoo::rollmean(new_tests_smoothed_per_thousand, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(Namerica_new_ave2)

p_Epi_Namerica<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_Namerica, aes(x = days, y = cases_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+ 
  geom_line(data=Namerica_new_ave2, aes(x=days, y=new_tests_smoothed*120),linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = "", trans = "sqrt",limits=c(0, 2000),breaks = c(20, 100, 200, 500, 1000, 1500, 2000),labels = c(20, 100,200, 500, 1000, 1500, 2000),
    sec.axis = sec_axis(~./120, name=" ",breaks = c(1, 2.5, 5, 10, 15),labels = c(1, 2.5, 5, 10, 15))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('North America')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

p_Epi_Namerica

###Samerica
P_Samerica <- prop.table(table(Samerica$days, Samerica$Variants2), margin=1)

temp_Samerica<-as.data.frame(P_Samerica)
names(temp_Samerica)[1] <- 'days'

###smoothing for genomic data
temp_Samerica <- temp_Samerica %>%
  group_by(Var2) %>%
  dplyr::mutate(Freq_7day = zoo::rollmean(Freq, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp_Samerica)

temp_Samerica$days<-as.Date(cut(as.Date(temp_Samerica$days,format='%Y-%m-%d'),
                                breaks = "day",
                                start.on.monday = FALSE))

temp2_Samerica<-Samerica_owid[c("days","new_cases_per_million")]
head(temp2_Samerica)

temp3_Samerica<-join(temp_Samerica, temp2_Samerica,
                     type = "left")

tail(temp3_Samerica)

temp3_Samerica$days<-as.Date(cut(as.Date(temp3_Samerica$days,format='%Y-%m-%d'),
                                 breaks = "day",
                                 start.on.monday = FALSE))

temp3_Samerica$cases_per_variant=temp3_Samerica$new_cases_per_million*temp3_Samerica$Freq_7day

###smoothing
temp3_Samerica <- temp3_Samerica %>%
  group_by(Var2) %>%
  dplyr::mutate(cases_per_variant_7day = zoo::rollmean(cases_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_Samerica)

Samerica_new_ave2 <- Samerica_new_ave %>%
  dplyr::mutate(new_tests_smoothed = zoo::rollmean(new_tests_smoothed_per_thousand, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(Samerica_new_ave2)

p_Epi_Samerica<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_Samerica, aes(x = days, y = cases_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+ 
  geom_line(data=Samerica_new_ave2, aes(x=days, y=new_tests_smoothed*120),linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = "", trans = "sqrt",limits=c(0, 2000),breaks = c(20, 100, 200, 500, 1000, 1500, 2000),labels = c(20, 100,200, 500, 1000, 1500, 2000),
    sec.axis = sec_axis(~./120, name=" ",breaks = c(1, 2.5, 5, 10, 15),labels = c(1, 2.5, 5, 10, 15)))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('South America')+
  theme(axis.title.y = element_text(hjust = 0.5, size=11, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

p_Epi_Samerica

####Deaths per variant plots
###GISAID data
africa$days20 <- as.Date(africa$days)+20
asia$days20 <- as.Date(asia$days)+20
europe$days20 <- as.Date(europe$days)+20
oceania$days20 <- as.Date(oceania$days)+20
Namerica$days20 <- as.Date(Namerica$days)+20
Samerica$days20 <- as.Date(Samerica$days)+20

###OWID data
owid<-read_excel('owid.xlsx')
asia_owid <- subset(owid, location == 'Asia')
asia_owid$days<-as.Date(cut(asia_owid$date,breaks = "day",start.on.monday = FALSE))
asia_owid$people_fully_vaccinated[is.na(asia_owid$people_fully_vaccinated)] <- 0
asia_owid$proportion_people_fully_vaccinated <- (asia_owid$people_fully_vaccinated/asia_owid$population)

africa_owid <- subset(owid, location == 'Africa')
africa_owid$days<-as.Date(cut(africa_owid$date,breaks = "day",start.on.monday = FALSE))
africa_owid$proportion_people_fully_vaccinated <- (africa_owid$people_fully_vaccinated/africa_owid$population)

europe_owid <- subset(owid, location == 'Europe')
europe_owid$days<-as.Date(cut(europe_owid$date,breaks = "day",start.on.monday = FALSE))
europe_owid$proportion_people_fully_vaccinated <- (europe_owid$people_fully_vaccinated/europe_owid$population)

oceania_owid <- subset(owid, location == 'Oceania')
oceania_owid$days<-as.Date(cut(oceania_owid$date,breaks = "day",start.on.monday = FALSE))
oceania_owid$proportion_people_fully_vaccinated <- (oceania_owid$people_fully_vaccinated/oceania_owid$population)

Namerica_owid <- subset(owid, location == 'North America')
Namerica_owid$days<-as.Date(cut(Namerica_owid$date,breaks = "day",start.on.monday = FALSE))
Namerica_owid$proportion_people_fully_vaccinated <- (Namerica_owid$people_fully_vaccinated/Namerica_owid$population)

Samerica_owid <- subset(owid, location == 'South America')
Samerica_owid$days<-as.Date(cut(Samerica_owid$date,breaks = "day",start.on.monday = FALSE))
Samerica_owid$proportion_people_fully_vaccinated <- (Samerica_owid$people_fully_vaccinated/Samerica_owid$population)

###Asia 20 DAYS
P_asia20 <- prop.table(table(asia$days20, asia$Variants2), margin=1)

temp_asia20<-as.data.frame(P_asia20)
names(temp_asia20)[1] <- 'days'

head(temp_asia20)
temp_asia20$days<-as.Date(cut(as.Date(temp_asia20$days,format='%Y-%m-%d'),
                              breaks = "day",
                              start.on.monday = FALSE))

temp2_asia20<-asia_owid[c("days","new_deaths_per_million")]
head(temp2_asia20)

temp3_asia20<-join(temp_asia20, temp2_asia20,
                   type = "left")

tail(temp3_asia20)

temp3_asia20$days<-as.Date(cut(as.Date(temp3_asia20$days,format='%Y-%m-%d'),
                               breaks = "day",
                               start.on.monday = FALSE))

temp3_asia20$deaths_per_variant=temp3_asia20$new_deaths_per_million*temp3_asia20$Freq

###smoothing
temp3_asia20 <- temp3_asia20 %>%
  group_by(Var2) %>%
  dplyr::mutate(deaths_per_variant_7day = zoo::rollmean(deaths_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_asia20)

deaths_asia20<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_asia20, aes(x = days, y = deaths_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+
  geom_line(data=asia_owid, aes(x=days, y=proportion_people_fully_vaccinated*8),linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 8),breaks = c(0.2, 1, 2, 4, 6, 8),labels = c(0.2, 1, 2, 4, 6, 8),
    sec.axis = sec_axis(~./8, name=" ", breaks = c(0.1, 0.25, 0.5, 0.75, 1.00),labels = c(0.1, 0.25, 0.5, 0.75, 1.00))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Asia')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

deaths_asia20

###africa 20 DAYS
P_africa20 <- prop.table(table(africa$days20, africa$Variants2), margin=1)

temp_africa20<-as.data.frame(P_africa20)
names(temp_africa20)[1] <- 'days'

head(temp_africa20)
temp_africa20$days<-as.Date(cut(as.Date(temp_africa20$days,format='%Y-%m-%d'),
                                breaks = "day",
                                start.on.monday = FALSE))

temp2_africa20<-africa_owid[c("days","new_deaths_per_million")]
head(temp2_africa20)

temp3_africa20<-join(temp_africa20, temp2_africa20,
                     type = "left")

tail(temp3_africa20)

temp3_africa20$days<-as.Date(cut(as.Date(temp3_africa20$days,format='%Y-%m-%d'),
                                 breaks = "day",
                                 start.on.monday = FALSE))

temp3_africa20$deaths_per_variant=temp3_africa20$new_deaths_per_million*temp3_africa20$Freq

###smoothing
temp3_africa20 <- temp3_africa20 %>%
  group_by(Var2) %>%
  dplyr::mutate(deaths_per_variant_7day = zoo::rollmean(deaths_per_variant, k = 20, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_africa20)

deaths_africa20<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_africa20, aes(x = days, y = deaths_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+
  geom_line(data=africa_owid, aes(x=days, y=proportion_people_fully_vaccinated*8), linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 8),breaks = c(0.2, 1, 2, 4, 6, 8),labels = c(0.2, 1, 2, 4, 6, 8),
    sec.axis = sec_axis(~./8, name=" ", breaks = c(0.1, 0.25, 0.5, 0.75, 1.00),labels = c(0.1, 0.25, 0.5, 0.75, 1.00))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Africa')+
  theme(axis.title.y = element_text(hjust = 0.5, size=11, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

deaths_africa20

###europe 20 DAYS
P_europe20 <- prop.table(table(europe$days20, europe$Variants2), margin=1)

temp_europe20<-as.data.frame(P_europe20)
names(temp_europe20)[1] <- 'days'

head(temp_europe20)
temp_europe20$days<-as.Date(cut(as.Date(temp_europe20$days,format='%Y-%m-%d'),
                                breaks = "day",
                                start.on.monday = FALSE))

temp2_europe20<-europe_owid[c("days","new_deaths_per_million")]
head(temp2_europe20)

temp3_europe20<-join(temp_europe20, temp2_europe20,
                     type = "left")

tail(temp3_europe20)

temp3_europe20$days<-as.Date(cut(as.Date(temp3_europe20$days,format='%Y-%m-%d'),
                                 breaks = "day",
                                 start.on.monday = FALSE))

temp3_europe20$deaths_per_variant=temp3_europe20$new_deaths_per_million*temp3_europe20$Freq

###smoothing
temp3_europe20 <- temp3_europe20 %>%
  group_by(Var2) %>%
  dplyr::mutate(deaths_per_variant_7day = zoo::rollmean(deaths_per_variant, k = 14, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_europe20)

deaths_europe20<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_europe20, aes(x = days, y = deaths_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+
  geom_line(data=europe_owid, aes(x=days, y=proportion_people_fully_vaccinated*8), linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 8),breaks = c(0.2, 1, 2, 4, 6, 8),labels = c(0.2, 1, 2, 4, 6, 8),
    sec.axis = sec_axis(~./8, name=" ", breaks = c(0.1, 0.25, 0.5, 0.75, 1.00),labels = c(0.1, 0.25, 0.5, 0.75, 1.00))
   )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Europe')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))+
  theme(axis.title.y = element_text(hjust = 0.5, size=9.5, face="bold"))

deaths_europe20

###oceania 20 DAYS
P_oceania20 <- prop.table(table(oceania$days20, oceania$Variants2), margin=1)

temp_oceania20<-as.data.frame(P_oceania20)
names(temp_oceania20)[1] <- 'days'

head(temp_oceania20)
temp_oceania20$days<-as.Date(cut(as.Date(temp_oceania20$days,format='%Y-%m-%d'),
                                 breaks = "day",
                                 start.on.monday = FALSE))

temp2_oceania20<-oceania_owid[c("days","new_deaths_per_million")]
head(temp2_oceania20)

temp3_oceania20<-join(temp_oceania20, temp2_oceania20,
                      type = "left")

tail(temp3_oceania20)

temp3_oceania20$days<-as.Date(cut(as.Date(temp3_oceania20$days,format='%Y-%m-%d'),
                                  breaks = "day",
                                  start.on.monday = FALSE))

temp3_oceania20$deaths_per_variant=temp3_oceania20$new_deaths_per_million*temp3_oceania20$Freq

###smoothing
temp3_oceania20 <- temp3_oceania20 %>%
  group_by(Var2) %>%
  dplyr::mutate(deaths_per_variant_7day = zoo::rollmean(deaths_per_variant, k = 22, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_oceania20)

deaths_oceania20<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_oceania20, aes(x = days, y = deaths_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+
  geom_line(data=oceania_owid, aes(x=days, y=proportion_people_fully_vaccinated*8), linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 8),breaks = c(0.2, 1, 2, 4, 6, 8),labels = c(0.2, 1, 2, 4, 6, 8),
    sec.axis = sec_axis(~./8, name=" ", breaks = c(0.1, 0.25, 0.5, 0.75, 1.00),labels = c(0.1, 0.25, 0.5, 0.75, 1.00))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('Oceania')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

deaths_oceania20

###Namerica 20 DAYS
P_Namerica20 <- prop.table(table(Namerica$days20, Namerica$Variants2), margin=1)

temp_Namerica20<-as.data.frame(P_Namerica20)
names(temp_Namerica20)[1] <- 'days'

head(temp_Namerica20)
temp_Namerica20$days<-as.Date(cut(as.Date(temp_Namerica20$days,format='%Y-%m-%d'),
                                  breaks = "day",
                                  start.on.monday = FALSE))

temp2_Namerica20<-Namerica_owid[c("days","new_deaths_per_million")]
head(temp2_Namerica20)

temp3_Namerica20<-join(temp_Namerica20, temp2_Namerica20,
                       type = "left")

tail(temp3_Namerica20)

temp3_Namerica20$days<-as.Date(cut(as.Date(temp3_Namerica20$days,format='%Y-%m-%d'),
                                   breaks = "day",
                                   start.on.monday = FALSE))

temp3_Namerica20$deaths_per_variant=temp3_Namerica20$new_deaths_per_million*temp3_Namerica20$Freq

###smoothing
temp3_Namerica20 <- temp3_Namerica20 %>%
  group_by(Var2) %>%
  dplyr::mutate(deaths_per_variant_7day = zoo::rollmean(deaths_per_variant, k = 14, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_Namerica20)

deaths_Namerica20<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_Namerica20, aes(x = days, y = deaths_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+
  geom_line(data=Namerica_owid, aes(x=days, y=proportion_people_fully_vaccinated*8), linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 8),breaks = c(0.2, 1, 2, 4, 6, 8),labels = c(0.2, 1, 2, 4, 6, 8),
    sec.axis = sec_axis(~./8, name=" ", breaks = c(0.1, 0.25, 0.5, 0.75, 1.00),labels = c(0.1, 0.25, 0.5, 0.75, 1.00)))+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('North America')+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

deaths_Namerica20

###Samerica 20 DAYS
P_Samerica20 <- prop.table(table(Samerica$days20, Samerica$Variants2), margin=1)

temp_Samerica20<-as.data.frame(P_Samerica20)
names(temp_Samerica20)[1] <- 'days'

head(temp_Samerica20)
temp_Samerica20$days<-as.Date(cut(as.Date(temp_Samerica20$days,format='%Y-%m-%d'),
                                  breaks = "day",
                                  start.on.monday = FALSE))

temp2_Samerica20<-Samerica_owid[c("days","new_deaths_per_million")]
head(temp2_Samerica20)

temp3_Samerica20<-join(temp_Samerica20, temp2_Samerica20,
                       type = "left")

tail(temp3_Samerica20)

temp3_Samerica20$days<-as.Date(cut(as.Date(temp3_Samerica20$days,format='%Y-%m-%d'),
                                   breaks = "day",
                                   start.on.monday = FALSE))

temp3_Samerica20$deaths_per_variant=temp3_Samerica20$new_deaths_per_million*temp3_Samerica20$Freq

###smoothing
temp3_Samerica20 <- temp3_Samerica20 %>%
  group_by(Var2) %>%
  dplyr::mutate(deaths_per_variant_7day = zoo::rollmean(deaths_per_variant, k = 22, fill='extend')) %>%
  dplyr::ungroup()
head(temp3_Samerica20)

deaths_Samerica20<-
  ggplot() + 
  theme_classic()+
  scale_x_date(limits = c(min(dateVec), max=max(dateVec)), date_labels = "%b\n%y",date_breaks = "6 month")+
  scale_colour_manual(values=c('#a7c957',"#354f52","#52796f","#84a98c","#2081c3","#d74e09","#ffd166",'#ddb892',"grey"), labels=c("Omicron", "BA.4/5","BA.2","BA.1","Gamma","Delta",'Beta', "Alpha","Other"))+
  geom_density(data=temp3_Samerica20, aes(x = days, y = deaths_per_variant_7day, color = Var2),stat="identity",size=1, alpha = 0.55)+
  geom_line(data=Samerica_owid, aes(x=days, y=proportion_people_fully_vaccinated*8), linetype = "twodash", size=1, color = "black")+
  xlab('')+
  scale_y_continuous(
    name = " ", trans = "sqrt",limits=c(0, 8),breaks = c(0.2, 1, 2, 4, 6, 8),labels = c(0.2, 1, 2, 4, 6, 8),
    sec.axis = sec_axis(~./8, name=" ", breaks = c(0.1, 0.25, 0.5, 0.75, 1.00),labels = c(0.1, 0.25, 0.5, 0.75, 1.00))
  )+
  theme(plot.title = element_text(hjust = 0.5))+
  theme(legend.title=element_text(size=0))+
  theme(legend.position = "top", legend.text= element_text(size=10),legend.title=element_text(size=0), legend.spacing.x = unit(0.2, 'cm')) +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))+
  ggtitle('South America')+
  theme(axis.title.y = element_text(hjust = 0.5, size=11, face="bold"))+
  theme(plot.title = element_text(hjust = 0.5, size=10, face="bold"))

deaths_Samerica20
