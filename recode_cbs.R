#### Recoding CBS/NYT Pre-Wave (Cleaned March 2021)
rm(list = ls())


### region definitions
northeast_states = c("Massachusetts", "Pennsylvania", "Connecticut",
                     "New Hampshire", "New York", "New Jersey", "Maine",
                     "Rhode Island", "Vermont")
midwest_states = c("Minnesota", "Iowa", "Nebraska", "Ohio", "Indiana",
                   "Wisconsin", "Kansas", "Missouri", "Illinois", "Michigan",
                   "South Dakota", "North Dakota") 
south_states = c("Delaware", "Texas", "North Carolina", "Georgia", "Alabama",
                 "West Virginia", "Maryland", "District of Columbia",
                 "Kentucky", "South Carolina", "Tennessee", "Louisiana",
                 "Oklahoma", "Florida", "Arkansas", "Virginia", "Mississippi") 
west_states = c("California", "Oregon", "Washington", "Hawaii", "New Mexico",
                "Alaska", "Colorado", "Utah", "Wyoming", "Arizona", "Nevada",
                "Montana", "Idaho")


path_data= "/Users/Ciara/Dropbox/kpop/data/"

cbs_nytimes <- read_csv(paste0(path_data, "cbs_nytimes_oct28_nov1/31114163.csv"))

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

# HRCTrump16 asks who you plan to vote for; if depends or dn/no answer then LnHRCTrumpGJ16
# we can see that this isn't totally how the data look tho, we get some lean answers even after answering R, D
View(cbs_nytimes %>% group_by(HRCTrump16, LnHRCTrumpGJ16 ) %>% summarise(n()))


#there are 17 ppl who voted, who not answer who for on any questions (dn ask lean)
#to be consistent w pew, we will code refused voters as othter
cbs_nytimes %>% filter(Vt16AlrdyVtd == "Don't know/No answer") %>% group_by( HRCTrump16, LnHRCTrump16 ) %>% summarise(n())

#to match waht was done in pre-pew, don't know leaners will be considered other
cbs_nytimes <- cbs_nytimes %>% 
    mutate(recode_vote_2016 = 
               case_when( (str_detect(Vt16AlrdyVtd, "Trump") |
                               str_detect(HRCTrump16, "Trump") |
                               str_detect(LnHRCTrump16, "Trump")) ~ "Republican",
                          (str_detect(Vt16AlrdyVtd, "Clinton") |
                               str_detect(HRCTrump16, "Clinton") |
                               str_detect(LnHRCTrump16, "Clinton") ) ~ "Democrat",
                          (str_detect(Vt16AlrdyVtd, "Johnson|Stein|Other|No answer") |
                               str_detect(HRCTrump16, "Other") |
                               str_detect(LnHRCTrump16, "Johnson|Stein|Other|Don't"))
                          ~ "Other",
                         #to be explicit about ppl who won't vote
                          (str_detect(HRCTrump16, "Won't vote")|
                               str_detect(LnHRCTrump16, "Won't vote")) ~ NA_character_,
                         #remaining NAs
                          TRUE ~ NA_character_) 
               )
#double checking NAs
#NA is only those that won't vote or we have NA for even lean, or NA overall
cbs_nytimes %>% filter(is.na(recode_vote_2016)) %>% group_by(Vt16AlrdyVtd, HRCTrump16, LnHRCTrump16 ) %>% summarise(n())
#none w low vote probably answered R,D,O
cbs_nytimes %>% filter(is.na(recode_vote_2016)) %>% group_by(LklyVote16, HRCTrump16, LnHRCTrump16) %>% summarise(n())


###### how bias are we?

pres <- readRDS( "/Users/Ciara/Dropbox/kpop/application/data/election.rds"))

natl_margin <- pres %>%
    summarise(margin = (sum(demtotal) - sum(reptotal)) /
                  (sum(demtotal) + sum(reptotal))) %>%
    as.numeric()

natl_diff <- pres %>%
    summarise(margin = (sum(demtotal) - sum(reptotal)) /
                  (sum(totalvotes))) %>%
    as.numeric()*100


vote_contrast <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /
                           (recode_vote_2016Democrat + recode_vote_2016Republican))

vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 


pew <- readRDS("/Users/Ciara/Dropbox/kpop/application/data/pew_new.rds")
pew_srs <- svydesign(ids = ~1, data = pew)


#it's really not bad at least compared to pew
cbs_srs <- svydesign(ids = ~1, data = cbs_nytimes)
svycontrast(svymean(~recode_vote_2016, 
                    cbs_srs, na.rm = TRUE),
            vote_contrast)*100
svycontrast(svymean(~recode_vote_2016, 
                    pew_srs, na.rm = TRUE),
            vote_contrast)*100
natl_margin*100


svycontrast(svymean(~recode_vote_2016, 
                    cbs_srs, na.rm = TRUE),
            vote_diff)*100
svycontrast(svymean(~recode_vote_2016, 
                    pew_srs, na.rm = TRUE),
            vote_diff)*100
natl_diff



