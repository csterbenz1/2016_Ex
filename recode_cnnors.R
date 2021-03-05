#### Recoding CNN Pre-Wave (Cleaned March 2021)
rm(list = ls())
library(foreign)

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


cnn_orc <- read.spss(paste0(path_data,"cnn_orc_oct20_oct23/usorccnn2016-1023.por")) %>% data.frame() %>%
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
## preference is P5
## intend to vote: V1

#outcome
#p1 asks horse race allowing other or neitther or dn/refused/undecided (9) if not R,D
#if p1 = niether or dn/refused asked to lean (followed in data hurrah)
#to be consistent w pre_pew, if dnk/refused/undecided then other
cnn_orc <- cnn_orc %>% 
    mutate(recode_vote_2016 = 
        case_when( (str_detect(P1, "TRUMP") |
                        str_detect(P1A, "TRUMP")) ~ "Republican",
                   (str_detect(P1, "CLINTON") |
                        str_detect(P1A, "CLINTON") ) ~ "Democrat",
                   (str_detect(P1, "OTHER|9") |
                        str_detect(P1A, "OTHER|NEITHER"))
                   ~ "Other",
                   #remaining NAs
                   TRUE ~ NA_character_) 
    )

#check NAs
#only NAs on everything
cnn_orc %>% filter(is.na(recode_vote_2016)) %>% group_by(P1, P1A) %>% summarise(n())






############ How off does this look? Compare w/ pew?

pres <- readRDS( "/Users/Ciara/Dropbox/kpop/application/data/election.rds")

natl_margin <- pres %>%
    summarise(margin = (sum(demtotal) - sum(reptotal)) /
                  (sum(demtotal) + sum(reptotal))) %>%
    as.numeric()*100

natl_diff <- pres %>%
    summarise(margin = (sum(demtotal) - sum(reptotal)) /
                  (sum(totalvotes))) %>%
    as.numeric()*100


vote_contrast <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /
                           (recode_vote_2016Democrat + recode_vote_2016Republican))

vote_diff <- quote((recode_vote_2016Democrat - recode_vote_2016Republican) /(recode_vote_2016Democrat + recode_vote_2016Republican + recode_vote_2016Other)) 


pew <- readRDS("/Users/Ciara/Dropbox/kpop/application/data/pew_new.rds")
pew_srs <- svydesign(ids = ~1, data = pew)


#We end up exactly split R and D:
cnn_orc %>% group_by(recode_vote_2016) %>% summarise(n())


cnn_srs <- svydesign(ids = ~1, data = cnn_orc)
svycontrast(svymean(~recode_vote_2016, 
                    cnn_srs, na.rm = TRUE),
            vote_contrast)*100
svycontrast(svymean(~recode_vote_2016, 
                    pew_srs, na.rm = TRUE),
            vote_contrast)*100
natl_margin

svycontrast(svymean(~recode_vote_2016, 
                    cnn_srs, na.rm = TRUE),
            vote_diff)*100
svycontrast(svymean(~recode_vote_2016, 
                    pew_srs, na.rm = TRUE),
            vote_diff)*100
natl_diff


