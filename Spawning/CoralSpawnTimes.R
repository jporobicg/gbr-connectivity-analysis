# coral spawn times for M&DS
# Doropoulos Nov 2021
# data source: Baird et al 2021
# central GBR only
# filter to Acropora and Merulinidae only

rm(list=ls())

library(here)
library(tidyverse)

data<-read.csv("SpawningDatabase.csv", header=T)

str(data)

data$Genus<-as.factor(data$Genus)
levels(data$Genus)

# Acroporidae only include Acropora (i.e., not Alveopora, Anacropora, Astreopora, Montipora)

# Merulindae =  Astrea, Caulastraea, Coelastrea, Cyphastrea, Dipsastraea, Echinopora, Favites, Goniastrea, Leptoria, Merulina,
#               Oulophyllia, Paragoniastrea, Pectinia, Platygyra

data<-data %>% filter(Genus=="Acropora"|
                        Genus=="Astrea"|
                        Genus=="Caulastraea"|
                        Genus=="Coelastrea"|
                        Genus=="Cyphastrea"|
                        Genus=="Dipsastraea"|
                        Genus=="Echinopora"|
                        Genus=="Favites"|
                        Genus=="Goniastrea"|
                        Genus=="Leptoria"|
                        Genus=="Merulina"|
                        Genus=="Oulophyllia"|
                        Genus=="Paragoniastrea"|
                        Genus=="Pectinia"|
                        Genus=="Platygyra")

family<-function(x) {
   if(x=="Acropora") y = "Acropora"
  if(x=="Astrea") y = "Merulinidae"
  if(x=="Caulastraea") y = "Merulinidae"
  if(x=="Coelastrea") y = "Merulinidae"
  if(x=="Cyphastrea") y = "Merulinidae"
  if(x=="Dipsastraea") y = "Merulinidae"
  if(x=="Echinopora") y = "Merulinidae"
  if(x=="Favites") y = "Merulinidae"
   if(x=="Goniastrea") y = "Merulinidae"
  if(x=="Leptoria") y = "Merulinidae"
  if(x=="Merulina") y = "Merulinidae"
  if(x=="Oulophyllia") y = "Merulinidae"
  if(x=="Paragoniastrea") y = "Merulinidae"
  if(x=="Pectinia") y = "Merulinidae"
   if(x=="Platygyra") y = "Merulinidae"
  return(y)
}

data[,"Family"]<-NA;
data$Family<-sapply(data$Genus,family)

# change any negative min values for Night After Full moon to 0 and summarise release time
names(data)

data$Date2 <-as.Date(data$Date, format = "%d-%b-%y")


str(data$Date2)

#data <- data %>% separate(Date, sep="-", into = c("day", "month", "year"))

data <- data %>% separate(Date2, sep="-", into = c("year", "month", "day"))
str(data$year)
data$Year_Fact<-data$year; data$Year_Fact<-as.factor(data$Year_Fact); levels(data$Year_Fact)

data$DoSRtNFM<-as.numeric(data$DoSRtNFM)

# Date of Spawning Relative to Nearest Full Moon - ALL YEARS #
data$lat_text=as.character(data$Latitude)
Spawning_NightsAfterFullMoon <- data %>%
  mutate(DoSRtNFM = if_else(DoSRtNFM < 0, 0, DoSRtNFM)) %>%   #making every negative value 0; turn off if want to see negative values
  #group_by(Family) %>%
  group_by(Family, year,lat_text) %>%                                 #turn on if want summary table by year
  summarize(min = min(DoSRtNFM),
            max = max(DoSRtNFM),
            #quantile = quantile(DoSRtNFM),
            median = median(DoSRtNFM),
            mean = mean(DoSRtNFM),
            sd = sd(DoSRtNFM),
            se = sd(DoSRtNFM)/sqrt(n()),
            n = n())
Spawning_NightsAfterFullMoon
length(data$DoSRtNFM)
length(L
ggplot(data, aes(y=Latitude, x=DoSRtNFM,colour=Family))+geom_point() +facet_wrap(.~year)
data$Latclass<- ifelse(data$Latitude > -19, '>-19-North','<-19-South')
Acropora <- data %>% filter(Family=="Acropora")
ggplot(data, aes(y=Latitude, x=DoSRtNFM, colour=Family))+geom_point()

Merulinidae<- data %>% filter(Family=="Merulinidae")

data2<-data %>% filter(year >= 2010)
cbPalette <- c('#f98b65', '#cae0fe',  '#6d798b',  '#894c38',  '#f1b435', "#56B4E9", '#516091', '#F67280', '#D6F8B8')
#png('Spawning_215.png', height=900,width=1100)
ggplot(data, aes(x=DoSRtNFM, fill=Family)) +
        geom_bar(position='dodge') + scale_fill_manual(values=cbPalette)+
        ggtitle("All years (1981-2016)")+xlab("Date of Spawning Relative to Nearest Full Moon")+
        scale_x_continuous(breaks = c(0:12), limits = c(0,12))  +
        theme_bw()+theme(text = element_text(size=20))+
        geom_vline(aes(xintercept = median(Acropora$DoSRtNFM, na.rm=TRUE)), colour=cbPalette[1],lwd=2)+
        geom_vline(aes(xintercept = median(Merulinidae$DoSRtNFM, na.rm=TRUE)), colour=cbPalette[3],lwd=2)
#dev.off()
# start times of spawning ALL YEARS #
str(data)
data$Start_decimal

Spawning_StartTimeAfterSunset <- data %>%
  drop_na(Start_decimal) %>%
  group_by(Family) %>%
  #group_by(Family, year) %>%                                 #turn on if want summary table by year
  summarize(min = min(Start_decimal),
            max = max(Start_decimal),
            median = median(Start_decimal),
            mean = mean(Start_decimal),
            sd = sd(Start_decimal),
            se = sd(Start_decimal)/sqrt(n()),
            n = n())
Spawning_StartTimeAfterSunset

Acropora <- data %>% filter(Family=="Acropora")
Merulinidae<- data %>% filter(Family=="Merulinidae")

library(gridExtra)
library(ggpubr)
p1<-ggplot(data %>% drop_na(Start_decimal), aes(x=Start_decimal, fill=Family)) +
  geom_density(alpha=0.5) + scale_fill_manual(values=cbPalette)+
  ggtitle("All years (1981-2016)")+xlab("")+
  scale_x_continuous(breaks = c(17:25), limits = c(17,25))  +
  theme_bw()+theme(text = element_text(size=20),legend.position = c(0.1, 0.85))+
  geom_vline(aes(xintercept = median(Acropora$Start_decimal, na.rm=TRUE)), colour=cbPalette[1],lwd=2)+
  geom_vline(aes(xintercept = median(Merulinidae$Start_decimal, na.rm=TRUE)), colour=cbPalette[2],lwd=2)

# end times of spawning ALL YEARS #
str(data)
data$End_decimal

Spawning_EndTimeAfterSunset <- data %>%
  drop_na(End_decimal) %>%
  group_by(Family) %>%
  #group_by(Family, year) %>%                                 #turn on if want summary table by year
  summarize(min = min(End_decimal),
            max = max(End_decimal),
            median = median(End_decimal),
            mean = mean(End_decimal),
            sd = sd(End_decimal),
            se = sd(End_decimal)/sqrt(n()),
            n = n())
Spawning_EndTimeAfterSunset

Acropora <- data %>% filter(Family=="Acropora")
Merulinidae<- data %>% filter(Family=="Merulinidae")

p2<-ggplot(data %>% drop_na(End_decimal), aes(x=End_decimal, fill=Family)) +
  geom_density(alpha=0.5) + scale_fill_manual(values=cbPalette)+
  ggtitle("All years (1981-2016)")+xlab("Hour of Spawning Time")+
  scale_x_continuous(breaks = c(17:25), limits = c(17,25))  +
  theme_bw()+theme(text = element_text(size=20),legend.position = c(0.1, 0.85))+
  geom_vline(aes(xintercept = median(Acropora$End_decimal, na.rm=TRUE)), colour=cbPalette[1],lwd=2)+
  geom_vline(aes(xintercept = median(Merulinidae$End_decimal, na.rm=TRUE)), colour=cbPalette[2],lwd=2)




p3<-ggplot(data %>% drop_na(End_decimal), aes(x=End_decimal, fill=Family)) +
  geom_bar(aes(y = after_stat(prop), group = Family), alpha = 0.5,position="dodge") + scale_fill_manual(values=cbPalette)+
  ggtitle("All years (1981-2016)")+xlab("Hour of Spawning Time")+
  scale_x_continuous(breaks = c(17:25), limits = c(17,25))  +
  theme_bw()+theme(text = element_text(size=20),legend.position = c(0.1, 0.85))+
  geom_vline(aes(xintercept = median(Acropora$End_decimal, na.rm=TRUE)), colour=cbPalette[1],lwd=2)+
  geom_vline(aes(xintercept = median(Merulinidae$End_decimal, na.rm=TRUE)), colour=cbPalette[2],lwd=2)

data_round <- data %>% mutate(Round_start = round(Start_decimal), Round_end = round(End_decimal))

p4<-ggplot(data_round %>% drop_na(Round_start), aes(x=Round_start, fill=Family)) +
  geom_bar(aes(y = after_stat(prop), group = Family), alpha = 1,position="dodge") + scale_fill_manual(values=cbPalette)+
  ggtitle("Start Spawning")+xlab("Hour of Spawning Time")+
  scale_x_continuous(breaks = c(17:25), limits = c(17,25))  +
  theme_bw()+theme(text = element_text(size=20),legend.position = c(0.1, 0.85))+
  geom_vline(aes(xintercept = median(Acropora$End_decimal, na.rm=TRUE)), colour=cbPalette[1],lwd=2)+
    geom_vline(aes(xintercept = median(Merulinidae$End_decimal, na.rm=TRUE)), colour=cbPalette[2],lwd=2)

p5<-ggplot(data_round %>% drop_na(Round_end), aes(x=Round_end, fill=Family)) +
  geom_bar(aes(y = after_stat(prop), group = Family), alpha = 1.0,position="dodge") + scale_fill_manual(values=cbPalette)+
  ggtitle("End Spawnin")+xlab("Hour of Spawning Time")+
  scale_x_continuous(breaks = c(17:25), limits = c(17,25))  +
  theme_bw()+theme(text = element_text(size=20),legend.position = c(0.1, 0.85))+
  geom_vline(aes(xintercept = median(Acropora$End_decimal, na.rm=TRUE)), colour=cbPalette[1],lwd=2)+
  geom_vline(aes(xintercept = median(Merulinidae$End_decimal, na.rm=TRUE)), colour=cbPalette[2],lwd=2)


library(patchwork)

p4 + p5          # side by side
p4 / p5          # one above the other
(p4 | p5) / p2   # combine multiple

prop_table <- data_round %>%
  drop_na(Round_start) %>%
  count(Family, Round_start) %>%
  group_by(Family) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

prop_table



summary(data)


aes(y = after_stat(prop), group = Family), alpha = 0.5



png('Spawning_time.png',width=1200,height=900)
ggarrange(p1, p2, labels = c("A", "B"),
          ncol = 1, nrow = 2, align = "v")
dev.off()


ggplot(data, aes(x = DoSRtNFM, fill = Family)) +
  geom_bar(aes(y = ..prop.., group = Family), position = "dodge") +
  scale_fill_manual(values = cbPalette) +
  ggtitle("All years (1981-2016)") +
  xlab("Date of Spawning Relative to Nearest Full Moon") +
  ylab("Probability") +
  scale_x_continuous(breaks = c(0:12), limits = c(0, 12)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  geom_vline(aes(xintercept = median(Acropora$DoSRtNFM, na.rm = TRUE)),
             colour = cbPalette[1], lwd = 2) +
  geom_vline(aes(xintercept = median(Merulinidae$DoSRtNFM, na.rm = TRUE)),
             colour = cbPalette[3], lwd = 2)

library(dplyr)

data_prop <- data %>%
  count(Family, DoSRtNFM) %>%
  mutate(prob = n / sum(n))

ggplot(data_prop, aes(x = DoSRtNFM, y = prob, fill = Family)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = cbPalette) +
  ggtitle("All years (1981-2016)") +
  xlab("Date of Spawning Relative to Nearest Full Moon") +
  ylab("Probability") +
  scale_x_continuous(breaks = 0:12, limits = c(0, 12)) +
  theme_bw() +
  theme(text = element_text(size = 20)) +
  geom_vline(aes(xintercept = median(Acropora$DoSRtNFM, na.rm = TRUE)),
             colour = cbPalette[1], lwd = 2) +
  geom_vline(aes(xintercept = median(Merulinidae$DoSRtNFM, na.rm = TRUE)),
             colour = cbPalette[3], lwd = 2)
