#### Recoding ABC Pre-Wave (Cleaned March 2021)
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

#for comparability with pew: ask who you vote for -> if dnk -> ask how lean -> if lean or vote for not R, D -> ask between R,D who would you choose


#outcome
#q2a = already voted, top4 or other
#q3 who vote for if eleciton today among top 4
#q4 = if q3 is none of these, dn/no opion, or na/refused, push to lean amnog top 4;  (556 total)
#since this did not force a horserace it's not directly comparable with pew
#following what we did with the pre pew, DK and refused are grouped as other


abc_wapo <- abc_wapo %>% 
    mutate(recode_vote_2016 =
           case_when( (str_detect(q2a, "Trump") |
                         str_detect(q3, "Trump") | 
                           str_detect(q4, "Trump") ) ~ "Republican",
                      (str_detect(q2a, "Clinton") |
                           str_detect(q3, "Clinton") | 
                           str_detect(q4, "Clinton")) ~ "Democrat",
                       (str_detect(q2a, "Johnson|Stein|Other|Refused") |
                            str_detect(q3, "Johnson|Stein|Other") | 
                            str_detect(q4, "Johnson|Stein|Other|None|DK|Refused")) ~ "Other",
                        #to be explicit
                     (str_detect(q3, "Would not vote")|
                          str_detect(q4, "Would not vote")) ~ NA_character_,
                     TRUE ~ NA_character_) )
#double checking NAs
#NA is would not vote (q3 OR Q4) OR NA on all three

abc_wapo %>% filter(is.na(recode_vote_2016)) %>% group_by( q2a, q3, q4) %>% summarise(n())
#note that would not vote dn match when asked separately lieklihood of voting (q2)
abc_wapo %>% filter(is.na(recode_vote_2016)) %>% group_by(q2, q2a, q3, q4) %>% summarise(n())




############ How off does this look? Compare w/ pew?

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


#it's really not bad
abc_srs <- svydesign(ids = ~1, data = abc_wapo)
svycontrast(svymean(~recode_vote_2016, 
                    abc_srs, na.rm = TRUE),
            vote_contrast)*100
svycontrast(svymean(~recode_vote_2016, 
                 pew_srs, na.rm = TRUE),
            vote_contrast)*100
natl_margin*100

svycontrast(svymean(~recode_vote_2016, 
                    abc_srs, na.rm = TRUE),
            vote_diff)*100
svycontrast(svymean(~recode_vote_2016, 
                    pew_srs, na.rm = TRUE),
            vote_diff)*100
natl_diff






#drop people who say they will not vote ( note that this still includes people who report less that 50% chance of voting) (62)
abc_wapo %>% filter(q2 == "(VOL) Don't think will vote") %>% group_by(recode_vote_2016) %>% summarise(n())

abc_wapo <- abc_wapo %>% filter(q2 != "(VOL) Don't think will vote")
