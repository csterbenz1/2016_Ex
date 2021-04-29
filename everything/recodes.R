library(tidyverse)


path_data= "/Users/Ciara/Dropbox/kpop/data/"

#setwd("~/Dropbox/Documents/2020__2021/work/kpop/data/")

abc_wapo <- read_csv(paste0(path_data, "abc_wapo_oct20_nov06/31115521.csv"))

names(abc_wapo)[which(str_detect(names(abc_wapo), "party"))]
table(abc_wapo$partlean)

## Age: q910 has age (with q910 followup for refused), and agebreak has (5 way age: 18-29, 30-39, 40-49, 50-64, 65+)
## Gender: qsex
## race: racenet
## Region: stcode (state)
## PID: q901 + q904 or partlean, depending on how we want leaners
## Educ: educnew (4 way: HS or less, some college, college degree, post-grad) (q909 more fine grained)
## Income: income and income2 (income has up to 100k, income2 has over 100k)
## Relig: relnet (also in q911 for more detailed info)
## Bornagain: q911b
## Church attendance: none
## horserace: q2a (already voted folks) + q3 + q4 (for leaners)
## intend to vote: q2

cbs_nytimes <- read_csv("./cbs_nytimes_oct28_nov1/31114163.csv")

names(cbs_nytimes)[which(str_detect(names(cbs_nytimes), "educ"))]
table(cbs_nytimes$HRCTrump16)


## Age: age, with agea follow-up for refused
## Gender: sex
## race: raceeth
## Region: state and state2
## PID: prty (lean has leaners)
## Educ: educ, educ4
## Income: income
## Relig: religion
## Bornagain: evan
## Church attendance: none
## horserace: HRCTrump16 + LnHRCTrump16 + Vt16AlrdyVtd
## intend to vote: LklyVote16

library(foreign)
cnn_orc <- read.spss("./cnn_orc_oct20_oct23/usorccnn2016-1023.por") %>% data.frame() %>%
  as_tibble()

names(cnn_orc)[which(str_detect(names(cnn_orc), "REL"))]
table(cnn_orc$V1)


## Age: AAGE with AAGE2 as follow-up
## Gender: SEX
## race: RACE
## Region: STATE
## PID: PARTYID
## Educ: EEDUC
## Income: IINCOME (but also INCOME and INCEXIT for slightly different breaks)
## Relig: none
## Bornagain: bornagn
## Church attendance: none
## horserace: P1 + P1A
## intend to vote: V1

fox <- read_csv("./fox_news_nov1_nov3/31116949.csv")


names(fox)[which(str_detect(names(fox), "age"))]
table(fox$q0c)

## Age: age (5 way bucket: 18-29, 30-45, 45-54, 55-64, 65+)
## Gender: gender
## race: raceq
## Region: N/A (there is a region, but no mapping for what it is)
## PID: party
## Educ: educ
## Income: income
## Relig: N/A
## Bornagain: brnagain
## Church attendance: attend
## horserace: q4
## intend to vote: q0c