rm(list = ls())

library(MASS)
library(tidyverse)
library(survey)
library(srvyr)
library(randomForest)
library(foreign)

### functions
load_to_env <- function (rdata, env = new.env()) {
    load(rdata, env)
    return(env)
}

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

################################################################################
#### Pew #######################################################################
################################################################################

### load original Pew .sav file
#late october pew
pew <- read.spss("./data/Pew/Oct16 public.sav", to.data.frame=TRUE)

#post election pew
pew <- read.spss("./data/Pew/Nov16-Post-Election-public/Nov16 Post-Election public.sav", to.data.frame=TRUE)
    

#check pre and post 3rd party voting
# pew_pre %>% group_by(q10) %>% summarise(n())
# pew_pre %>% group_by(q10a) %>% summarise(n())
# #124 for Johnson/lean Johnson (108+16 lean), 
# #57 lean for Stein/ lean Stein (51 + 6 lean), 
# #31 for other/refused to lean, 
# nrow(pew_pre)
# #so not excluding non response this is
# nas <- 463
# (124/(nrow(pew_pre)-nas))*100
# (48/nrow(pew))*100
# 
# (57/(nrow(pew_pre)-nas))*100
# (26/nrow(pew))*100
# 
# (31/(nrow(pew_pre)-nas))*100
# (29/nrow(pew))*100
# pew %>% group_by(q4) %>% summarise(n())
# 
# (124+57+31)/(nrow(pew_pre))
# (48+26+29)/nrow(pew)
# 
# #gary johnson 48 jill stein 26, oher 29


#what about excluding NAs?

## recode forced major party vote (Q11)
#was pew originally loaded with load_to_env() so this worked?

#in pre-election poll Q11 forced a horse race betwen dem and rep (if you were forced to choose among) with q11a asking which way do you lean with other options
#q10 asked if you voted today who would you vote for among trump, clinton, stein, and johnson, 10a asked if they lean for these candidates if responded dnk
#i'm not actually sure how they made the q11horse2 varibale 
#it is not q10 + q10a
#q11 is only asked if they ansewr a third candidate for q10 or dnk
#q11horse2 is not q11 + q11a + q10 + q10a? that totals 1065 for hilary but its only 1052
#ok so that drops 13 which is probably the q11a
#is it the q11 + q10 + q10a? trump: 816+42 + 70 = 928 yes and hilary: 944 +28 + 80 = 1052 yes
#what is other then: 53 = q11, don't know 87 = q11
#so that seems like? not what we want????

#in post election poll we have: orignQ11HORSE2, origq10horse (same as origq10horseGP), origq13HORSE (original q13 in pre election poll was is your vote more against hilary or pro trump etc so the answers of lean here make no sense, though it also has origq13a and b that are these answers so basically juts no idea what this horse variable is)
#also a origQ14HOURSE2 question

#FOR POST:
cleannames <- colnames(pew)
colnames(pew)[-22] <- gsub( "orig","", cleannames)[-22]

pew <- #pew_env$pew_201610 %>%
    pew %>%
    mutate(recode_vote_2016 = 
             case_when(str_detect(q4, "Trump") ~ "Republican",
                       str_detect(q4, "Clinton") ~ "Democrat",
                       is.na(q4) ~ NA_character_,
                       TRUE ~ "Other" 
                       ))


#FOR PRE:
# pew <- #pew_env$pew_201610 %>%
#     pew %>%
#     mutate(recode_vote_2016 = 
#                case_when(str_detect(Q11HORSE2, "Trump") ~ "Republican",
#                          str_detect(Q11HORSE2, "Clinton") ~ "Democrat",
#                          is.na(Q11HORSE2) ~ NA_character_,
#                          TRUE ~ "Other" 
#                ))

## Start by looking at missingness (fraction):
lapply(pew[, c("age", "sex", "racethn", "state", "party", "educ2")], 
       function(x) (sum(is.na(x)) + sum(grepl("Don't know/Refused", x))) / length(x))

#age and party and educ2 have no NAs, only don't know/refused
#racethn does have NAs and no don't know/refused

## First recodes
pew <- pew %>% mutate(
  # age
  recode_age = ifelse(age == "Don't know/Refused (VOL.)", 
                      NA, 
                      as.numeric(as.character(age))),

  # gender
  recode_female = case_when(sex == "Female" ~ "Female",
                            TRUE ~ "Male"),
  # race/ethnicity
  # Note: combines missing with Other
  recode_race = case_when(racethn == "White, non-Hisp" ~ "White",
                          racethn == "Black, non-Hisp" ~ "Black",
                          racethn == "Hispanic" ~ "Hispanic",
                          TRUE ~ "Other"),
  # region (note there are no cases that do not fall into these groups such that we fill TRUE ~ South)
  recode_region = case_when(state %in% northeast_states ~ "Northeast",
                            state %in% west_states ~ "West",
                            state %in% midwest_states ~ "Midwest",
                            state %in% south_states ~ "South",
                            TRUE ~ "South"),
  
  # party -- combines refused + no preference with Independents
  recode_pid_3way = case_when( party == "Democrat" ~ "Dem",
                               party == "Republican" ~ "Rep",
                               TRUE ~ "Ind"),
  
  # state
  recode_inputstate = state,
  
  # education -- leaving missing for now (14 Don't now)
  recode_educ = factor(case_when( 
    educ2 == "Less than high school (Grades 1-8 or no formal schooling) " ~ "No HS",
    educ2 == "High school incomplete (Grades 9-11 or Grade 12 with NO diploma)" ~ "No HS",
    educ2 == "High school graduate (Grade 12 with diploma or GED certificate)" ~ "High school graduate",
    educ2 == "Some college, no degree (includes some community college)" ~ "Some college",
    educ2 == "Two year associate degree from a college or university" ~ "2-year",
    educ2 == "Four year college or university degree/Bachelor's degree (e.g., BS, BA, AB)" ~ "4-year",
    educ2 == "Some postgraduate or professional schooling, no postgraduate degree" ~ "Post-grad",
    educ2 == "Postgraduate or professional degree, including master's, doctorate, medical or law degree" ~ "Post-grad",
    TRUE ~ NA_character_), 
    levels = c("No HS", "High school graduate", "Some college", "2-year", "4-year", "Post-grad"))
)


#XXX check above wraning

#### adding indicator for any missingness to use in simulations 
# XXXXXXX ADD IN MISSING ON NEW VARS
#note that region and inputstate and pid and female should not ever have any missing
pew <- pew %>% 
  mutate(missing = ifelse(is.na(recode_age) | is.na(recode_educ), 1, 0))

## modeling age
## listwise delete missing on age/education
#45 observations dropped= 41 missingon age + 4 missing on educ
age_model <- glm(recode_age ~ recode_female + recode_race + recode_region + 
                   recode_pid_3way + recode_educ + child,
                 family = Gamma(),
                 data = pew)


#reports 41 observations dropped = those missing on age
age_model_no_educ <- glm(recode_age ~ recode_female + recode_race + 
                           recode_region + recode_pid_3way + child,
                 family = Gamma(),
                 data = pew)

## recode age and age buckets
pew <- pew %>%
  mutate(recode_age = case_when(
    ## if age is not missing, keep current age
    !is.na(recode_age) ~ recode_age,
    ## if age is missing but not eduction, use age_model
    is.na(recode_age) & !is.na(recode_educ) ~ round(predict(age_model, ., type = "response"), 0),
    ## if age and education are missing, use age_model_no_educ
    is.na(recode_educ) ~ round(predict(age_model_no_educ, ., type = "response"), 0)
  ),
  
  ## four way age bucket; should be no na's at this point
  recode_age_bucket = factor(case_when( recode_age <= 35 ~ "18 to 35",
                                        recode_age <= 50 ~ "36 to 50",
                                        recode_age < 65 ~ "51 to 64",
                                        !is.na(recode_age) ~ "65+"),
                             levels = c("18 to 35", "36 to 50", "51 to 64", "65+")),
  
  ## three way age bucket; should be no NAs at this point
  recode_age_3way = case_when( recode_age <= 50 ~ "a_18to50",
                               recode_age < 65 ~ "b_51to64",
                               !is.na(recode_age) ~ "c_65")
  )

## education model (ordered logit)
#drops 14 missing on educ
educ_model <- polr(recode_educ ~ recode_female + recode_race + recode_region + 
                     recode_pid_3way + recode_age + child,
                   data = pew)

pew <- pew %>%
  ## impute missing educ bucket
  mutate(recode_educ = case_when(
    ## if not missing, use current education
    !is.na(recode_educ) ~ recode_educ,
    ## if missing, use educ_model prediction
    is.na(recode_educ) ~ predict(educ_model, newdata =., type = "class")
  ),
  
  ## education 3 way bucket among white voters, no split among non-white
  recode_educ_wh_3way = factor(case_when(
    recode_race != "White" ~ "No Split",
    recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
    recode_educ %in% c("2-year", "4-year") ~ "College",
    TRUE ~ "Post-grad"), levels = c("No Split", "No College", "College", "Post-grad")),
  
  recode_educ_3way = factor(case_when(
      recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
      recode_educ %in% c("2-year", "4-year") ~ "College",
      TRUE ~ "Post-grad"), 
  levels = c("No College", "College", "Post-grad"))
)  
  
## Make some interactions 
pew <- pew %>%
  mutate(
  # interactions
    recode_race_educ_reg = as.factor(paste(recode_region, case_when(
          recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
          TRUE ~ recode_race), sep = ", "))
)




####### New Variables: income, bornagain, church attendance,
### Pew
pew <- pew %>% mutate(
    recode_relig = factor(case_when(
        #say unitarians are protestant (only 4 of these), Christians are protestant (269)
        relig == "Protestant (Baptist, Methodist, Non-denominational, Lutheran...)" ~ "Protestant",
        relig == "Christian (VOL.)" ~ "Protestant",
        relig == "Unitarian (Universalist) (VOL.)" ~ "Protestant",
        relig == "Jewish (Judaism)" ~ "Jewish",
        relig == "Roman Catholic (Catholic)" ~ "Catholic",
        relig == "Nothing in particular" ~ "Nothing in particular",
        relig == "Agnostic (not sure if there is a God)" ~ "Agnostic",
        relig == "Atheist (do not believe in God)" ~ "Atheist",
        relig == "Buddhist" ~ "Buddhist",
        relig == "Muslim (Islam)" ~ "Muslim",
        relig == "Mormon (Church of Jesus Christ of Latter-day Saints/LDS)" ~ "Mormon",
        relig == "Hindu" ~ "Hindu",
        relig == "Orthodox (Greek, Russian, or some other orthodox church)" ~ "Orthodox",
        #these people were asked a follow up if considered christain in $chr
        #out of 41 (31) of these: 11 (8) said yes christ, 29 (22) no, 1 refused
        relig == "Something else (SPECIFY:______)" ~ "Something else",
        #those who respond don't know/refuse =NA (36) (NB: also asked if consider christian in chr)
        #out of these 50 (36): 30 (24) said yes christian, 11 (5) said no, 9 (7) refused
        relig == "Don't Know/Refused (VOL.)" & chr == "Yes" ~ "Protestant",
        #leaves only 20 (12) missing
        TRUE ~ NA_character_),
        levels = c("Protestant", "Jewish", "Catholic", "Nothing in particular", 
                   "Agnostic", "Atheist", "Buddhist", "Muslim", "Mormon", "Hindu",
                   "Orthodox", "Something else")),
    
    #supposedly only asked:if chr = yes or if relig = 1-4,13
    #(Protestant, Catholic, Mormon, Orthodox, Christian)
    #though looking at the NAs that's not totally true some in these groups still have NAs
    recode_born = case_when(born == "Yes, would" ~ "Yes",
                            born == "No, would not" ~ "No",
                            #46 refused but we think it's safe to say these people are not
                            born == "Don't know/Refused (VOL.)" ~ "No",
                            #NA's (589) these are mainly non-christian denominations
                            #so I'd say safe to say they are not born again
                            TRUE ~ "No"),
    
    #the names are actually ok here
    #leaving the 35 (28) that response don't know as they are 
    recode_attndch = attend, 
    
    recode_income = factor(case_when(
        income == "Less than $10,000" ~ "<10k",
        income == "10 to under $20,000" ~ "10-20k",
        income == "20 to under $30,000" ~ "20-30k",
        income == "30 to under $40,000" ~ "30-40k",
        income == "40 to under $50,000" ~ "40-50k",
        income == "50 to under $75,000" ~ "50-100k",
        income == "75 to under $100,000" ~ "50-100k",
        income == "100 to under $150,000 [OR]" ~ "100-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs, but as precaution lump them with prefer not to say
        TRUE ~ NA_character_),
        levels = c("<10k", "10-20k", "20-30k", "30-40k", "40-50k", "50-100k",
                   "100-150k",">150k", "Prefer not to say")),
    recode_income_5way = factor(case_when(
        income == "Less than $10,000" ~ "<20k",
        income == "10 to under $20,000" ~ "<20k",
        income == "20 to under $30,000" ~ "20-50k",
        income == "30 to under $40,000" ~ "20-50k",
        income == "40 to under $50,000" ~ "20-50k",
        income == "50 to under $75,000" ~ "50-100k",
        income == "75 to under $100,000" ~ "50-100k",
        income == "100 to under $150,000 [OR]" ~ "100-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs
        TRUE ~ NA_character_),
        levels = c("<20k", "20-50k", "50-100k", "100-150k",">150k", "Prefer not to say")),
    
    recode_income_4way = factor(case_when(
        income == "Less than $10,000" ~ "<50k",
        income == "10 to under $20,000" ~ "<50k",
        income == "20 to under $30,000" ~ "<50k",
        income == "30 to under $40,000" ~ "<50k",
        income == "40 to under $50,000" ~ "<50k",
        income == "50 to under $75,000" ~ "50-100k",
        income == "75 to under $100,000" ~ "50-100k",
        income == "100 to under $150,000 [OR]" ~ "100-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs
        TRUE ~ NA_character_),
        levels = c("<50k", "50-100k", "100-150k",">150k", "Prefer not to say")),
    
    recode_income_3way = factor(case_when(
        income == "Less than $10,000" ~ "<50k",
        income == "10 to under $20,000" ~ "<50k",
        income == "20 to under $30,000" ~ "<50k",
        income == "30 to under $40,000" ~ "<50k",
        income == "40 to under $50,000" ~ "<50k",
        income == "50 to under $75,000" ~ "50-150k",
        income == "75 to under $100,000" ~ "50-150k",
        income == "100 to under $150,000 [OR]" ~ "50-150k",
        income == "$150,000 or more" ~ ">150k",
        #222 refused
        income == "[VOL. DO NOT READ] Don't know/Refused" ~ "Prefer not to say",
        #there should be no NAs
        TRUE ~ NA_character_),
        levels = c("<50k", "50-150k",">150k","Prefer not to say"))
)

## modeling missing income (uses filled in age and educ)
#this doesnt perform well so leaving refused as a category
# income_model <- polr(recode_income ~ recode_female + recode_race + recode_region + 
#                          recode_pid_3way + recode_educ + child + recode_age,
#                      data = pew)
# 
# #performs pretty poorly
# sum(predict(income_model,, type = "class") == na.omit(pew$recode_income))/sum(!is.na(pew$recode_income))
# unique(predict(income_model, type = "class"))
# unique( predict(income_model,newdata = pew[is.na(pew$recode_income),], type = "class"))


pew <- pew %>% mutate(
    #NOTE THAT THIS IS ACTUALLY 5WAY I CANT COUNT
    recode_relig_6way = factor(case_when(
        #clean up by saying those that are christian and not catholic is protestant
        recode_relig == "Protestant" ~ "Protestant",
        recode_relig == "Jewish" ~ "Jewish",
        recode_relig == "Catholic" ~ "Catholic",
        recode_relig == "Nothing in particular" ~ "Not religious",
        recode_relig == "Agnostic" ~ "Not religious",
        recode_relig == "Atheist" ~ "Not religious",
        recode_relig == "Buddhist" ~ "Other",
        recode_relig == "Muslim" ~ "Other",
        recode_relig == "Hindu" ~ "Other",
        recode_relig == "Mormon" ~ "Other",
        recode_relig == "Orthodox" ~ "Other",
        recode_relig == "Something else" & chr == "Yes"  ~ "Protestant",
        recode_relig == "Something else" & chr == "No"  ~ "Other",
        recode_relig == "Something else" & chr == "Don't know/Refused (VOL.)"  ~ "Not religious",
        #20 remaining NAs as Not religious
        TRUE ~ "Not religious"),
        levels = c("Protestant", "Jewish", "Catholic","Not religious", "Other")),

    recode_attndch_4way = factor(case_when(
        recode_attndch == "Once a week" ~ "Weekly",
        recode_attndch == "More than once a week" ~ "Weekly",
        recode_attndch == "Once or twice a month" ~ "Monthly",
        recode_attndch == "A few times a year" ~ "Yearly",
        recode_attndch == "Seldom" ~ "Never",
        recode_attndch == "Never" ~ "Never",
        recode_attndch == " Don't know/Refused (VOL.)" ~ "Never",
        #the 35 who refused say never go
        TRUE ~ "Never"),
        levels = c("Weekly","Monthly", "Yearly", "Never"))
)




# FOR PRE POLL ONLY:
#post-poll q1 is did you vote and if thy didnt they end the interview so we already subsetted to actual voters and now do not need to filter for plan to vote

## Remove those who say they definitely will not vote 
#463 NAs (all also have NA vote choice) and 46 don't plan to vote (some have vote pref but if they refuse then fine to drop)
#(5 of those that are NA here are missing also on age/educ so we get from 45 missing to 40 missing total)

pew <- pew%>%
  filter(plan1 %in% c("Plan to vote", "Already voted", "Don't know/Refused (VOL.)"))





#updated August 2020 for pre-pew
#post-pew done oct 2020
saveRDS(pew, "./data/pew_post_CS_q4.rds")





























################################################################################
#### CCES ######################################################################
################################################################################

cces_env <- load_to_env("./data/CCES/CCES16_Common_Content.RData")


# Common Content Weight -- commonweight_vv_post (post wave with vote validation)
cces_2016 <- cces_env$cces_2016 %>% 
    filter((CC16_401 == "I definitely voted in the General Election.") &
               !is.na(commonweight_vv_post))

## make outcome

# Horserace (post) question 410a
#this makes those that said they were not sure or did not vote to NA
cces_2016 <- cces_2016 %>%
    mutate( recode_vote_2016 = case_when(
                str_detect(CC16_410a, "Democrat") ~ "Democrat",
                str_detect(CC16_410a, "Republican") ~ "Republican",
                (as.numeric(CC16_410a) < 6 | as.numeric(CC16_410a) == 8) ~ "Other",
                #XXX this still includes 81 people who claim to have both definitly voted and on CC16_401 and didn't vote on CC16_410a; perhaps didnt vote for pres only
                #360: 229 = I'm not sure, 50 = NA, 81 = I didn't vote in this election
                TRUE ~ NA_character_) )


#what's missing in CCES? just a few party id's
lapply(cces_2016[, c("birthyr", "inputstate_post", "gender", "race", "pid3", "educ")], 
       function(x) (sum(is.na(x)) + sum(grepl("Skipped", x)) + sum(grepl("Not Asked", x))))


## demographic recodes
cces_2016$inputstate_post = droplevels(cces_2016$inputstate_post)

cces_2016 <- cces_2016 %>% mutate(
  recode_id = 1:nrow(cces_2016),
  
  # age
  recode_age = 2016 - as.numeric(as.character(birthyr)),
  recode_age_bucket = factor(case_when( recode_age <= 35 ~ "18 to 35",
                          recode_age <= 50 ~ "36 to 50",
                          recode_age < 65 ~ "51 to 64",
                          !is.na(recode_age) ~ "65+"), 
                          levels = c("18 to 35", "36 to 50", "51 to 64", "65+")),
  
  recode_age_3way = case_when( recode_age <= 50 ~ "a_18to50",
                          recode_age < 65 ~ "b_51to64",
                          !is.na(recode_age) ~ "c_65"),
  
  # gender
  recode_female = case_when(gender == "Female" ~ "Female",
                            TRUE ~ "Male"),
  
  # race: this means Mixed, Asian, Other, Native American, Middle Eastern all in Other
  recode_race = ifelse(race %in% c("White", "Black", "Hispanic"), as.character(race), "Other"),
  
  # region: there should be no states that do not fall in one of these categories (no NAs)
  recode_region = case_when(inputstate_post %in% northeast_states ~ "Northeast",
                            inputstate_post %in% west_states ~ "West",
                            inputstate_post %in% midwest_states ~ "Midwest",
                            inputstate_post %in% south_states ~ "South",
                            TRUE ~ "South"),
  
  # party: note 13 missing on pid are categorized as Indep
  recode_pid_3way = case_when( str_detect(pid3, "Democrat") ~ "Dem",
                        str_detect(pid3, "Republican") ~ "Rep",
                        TRUE ~ "Ind"),

  # educ
  recode_educ = factor(educ, levels = c("No HS", "High school graduate", 
                                        "Some college", "2-year", "4-year", 
                                        "Post-grad")),
  
  ## education 3 way bucket among white voters, no split among non-white
  recode_educ_wh_3way = factor(case_when(
    recode_race != "White" ~ "No Split",
    recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
    recode_educ %in% c("2-year", "4-year") ~ "College",
    TRUE ~ "Post-grad"), 
    levels = c("No Split", "No College", "College", "Post-grad")),
  
  recode_educ_3way = factor(case_when(
      recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
      recode_educ %in% c("2-year", "4-year") ~ "College",
      TRUE ~ "Post-grad"), 
      levels = c("No College", "College", "Post-grad")),
  
  # state
  recode_inputstate = inputstate_post,
  
  # interactions
  recode_race_educ_reg = as.factor(paste(recode_region, case_when(
      recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
      TRUE ~ recode_race), sep = ", "))
)

cces_2016$recode_inputstate <- droplevels(cces_2016$recode_inputstate)
cces_2016$recode_educ <- droplevels(cces_2016$recode_educ)



## save
#08/05 updated to fix the na's on vote choice ONLY
#saveRDS(cces_2016, "data/cces_CS.rds")



#### Inputting New Variables
cces_2016 <- cces_2016 %>% mutate(
    #if responded don't know/refused 
    #unitarians are protestant (only 4 of these), Christains are protestant (269)
    recode_relig = factor(case_when(
        religpew == "Protestant" ~ "Protestant",
        religpew == "Jewish" ~ "Jewish",
        religpew == "Roman Catholic" ~ "Catholic",
        religpew == "Nothing in particular" ~ "Nothing in particular",
        religpew == "Agnostic" ~ "Agnostic",
        religpew == "Atheist" ~ "Atheist",
        religpew == "Buddhist" ~ "Buddhist",
        religpew == "Muslim" ~ "Muslim",
        religpew == "Mormon" ~ "Mormon",
        religpew == "Hindu" ~ "Hindu",
        religpew == "Eastern or Greek Orthodox" ~ "Orthodox",
        religpew == "Something else" ~ "Something else",
        #NAs remain NAs 38
        TRUE ~ NA_character_),
        levels = c("Protestant", "Jewish", "Catholic", "Nothing in particular", 
                   "Agnostic", "Atheist", "Buddhist", "Muslim", "Mormon", "Hindu",
                   "Orthodox", "Something else")),
    
    #22 NAs
    recode_born = case_when(pew_bornagain == "No" ~ "No",
                            pew_bornagain == "Yes" ~ "Yes", 
                            #NAs (22) assume to be not
                            TRUE ~ "No"),
    
    #names ok; NAs (23) and don't know (319) 
    recode_attndch = pew_churatd,
    
    recode_income = factor(case_when(
        faminc == "Less than $10,000" ~ "<10k",
        faminc == "$10,000 - $19,999" ~ "10-20k",
        faminc == "$20,000 - $29,999" ~ "20-30k",
        faminc == "$30,000 - $39,999" ~ "30-40k",
        faminc == "$40,000 - $49,999" ~ "40-50k",
        faminc == "$50,000 - $59,999" ~ "50-100k",
        faminc == "$60,000 - $69,999" ~ "50-100k",
        faminc == "$70,000 - $79,999" ~ "50-100k",
        faminc == "$80,000 - $99,999" ~ "50-100k",
        faminc == "$100,000 - $119,999" ~ "100-150k",
        faminc == "$120,000 - $149,999" ~ "100-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs say they are also prefer not to say
        TRUE ~ "Prefer not to say"),
        levels = c("<10k", "10-20k", "20-30k", "30-40k", "40-50k", "50-100k",
                   "100-150k",">150k","Prefer not to say")),
    
    recode_income_5way = factor(case_when(
        faminc == "Less than $10,000" ~ "<20k",
        faminc == "$10,000 - $19,999" ~ "<20k",
        faminc == "$20,000 - $29,999" ~ "20-50k",
        faminc == "$30,000 - $39,999" ~ "20-50k",
        faminc == "$40,000 - $49,999" ~ "20-50k",
        faminc == "$50,000 - $59,999" ~ "50-100k",
        faminc == "$60,000 - $69,999" ~ "50-100k",
        faminc == "$70,000 - $79,999" ~ "50-100k",
        faminc == "$80,000 - $99,999" ~ "50-100k",
        faminc == "$100,000 - $119,999" ~ "100-150k",
        faminc == "$120,000 - $149,999" ~ "100-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs
        TRUE ~ "Prefer not to say"),
        levels = c("<20k", "20-50k", "50-100k", "100-150k",">150k","Prefer not to say")),
    recode_income_4way = factor(case_when(
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "$10,000 - $19,999" ~ "<50k",
        faminc == "$20,000 - $29,999" ~ "<50k",
        faminc == "$30,000 - $39,999" ~ "<50k",
        faminc == "$40,000 - $49,999" ~ "<50k",
        faminc == "$50,000 - $59,999" ~ "50-100k",
        faminc == "$60,000 - $69,999" ~ "50-100k",
        faminc == "$70,000 - $79,999" ~ "50-100k",
        faminc == "$80,000 - $99,999" ~ "50-100k",
        faminc == "$100,000 - $119,999" ~ "100-150k",
        faminc == "$120,000 - $149,999" ~ "100-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs
        TRUE ~ "Prefer not to say"),
        levels = c("<50k", "50-100k", "100-150k",">150k", "Prefer not to say")),
    
    recode_income_3way = factor(case_when(
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "Less than $10,000" ~ "<50k",
        faminc == "$10,000 - $19,999" ~ "<50k",
        faminc == "$20,000 - $29,999" ~ "<50k",
        faminc == "$30,000 - $39,999" ~ "<50k",
        faminc == "$40,000 - $49,999" ~ "<50k",
        faminc == "$50,000 - $59,999" ~ "50-150k",
        faminc == "$60,000 - $69,999" ~ "50-150k",
        faminc == "$70,000 - $79,999" ~ "50-150k",
        faminc == "$80,000 - $99,999" ~ "50-150k",
        faminc == "$100,000 - $119,999" ~ "50-150k",
        faminc == "$120,000 - $149,999" ~ "50-150k",
        faminc == "$150,000 - $199,999" ~ ">150k",
        faminc == "$200,000 - $249,999" ~ ">150k",
        faminc == "$250,000 - $349,999" ~ ">150k",
        faminc == "$350,000 - $499,999" ~ ">150k",
        faminc == "$500,000 or more" ~ ">150k",
        faminc == "$150,000 or more" ~ ">150k",
        #4788 do not say
        faminc == "Prefer not to say" ~ "Prefer not to say",
        #9 NAs
        TRUE ~ "Prefer not to say"),
        levels = c("<50k", "50-150k",">150k","Prefer not to say"))
)

#doesn't work so well so allow refused to be its own cat
# sum(is.na(cces_2016$recode_income))
# ## modeling missing income (uses filled in age and educ)
# # also 85 NAs on child18 
# income_model_cces <- polr(recode_income ~ recode_female + recode_race + recode_region + 
#                          recode_pid_3way + recode_educ + child18 + recode_age,
#                      data = cces_2016)
# 
# length(predict(income_model_cces,, type = "class"))
# 
# sum(!is.na(cces_2016$recode_income))
# #why are there so many nas
# sum(is.na(predict(income_model_cces,cces_2016[!is.na(cces_2016$recode_income),], type = "class")))
# 
# sum(predict(income_model_cces,, type = "class") == cces_2016[!is.na(cces_2016$recode_income) & !is.na(cces_2016$child18)  ,"recode_income"] %>% pull())/sum(!is.na(cces_2016$recode_income) & !is.na(cces_2016$child18) )



#CCES: cleaning up new vars
cces_2016 <- cces_2016 %>% mutate(
    recode_relig_6way = factor(case_when(
        recode_relig == "Protestant" ~ "Protestant",
        recode_relig == "Jewish" ~ "Jewish",
        recode_relig == "Catholic" ~ "Catholic",
        recode_relig == "Nothing in particular" ~ "Not religious",
        recode_relig == "Agnostic" ~ "Not religious",
        recode_relig == "Atheist" ~ "Not religious",
        recode_relig == "Buddhist" ~ "Other",
        recode_relig == "Muslim" ~ "Other",
        recode_relig == "Hindu" ~ "Other",
        recode_relig == "Mormon" ~ "Other",
        recode_relig == "Orthodox" ~ "Other",
        #there's a lot of these (2k) 
        #quite of few of them are v religios going to church weekly (300)
        #when asked what type of protestant many identified so call them protestant
        recode_relig == "Something else" & !is.na(religpew_protestant) ~ "Protestant",
        recode_relig == "Something else" & !is.na(religpew_protestant) ~ "Other",
        #38 NAs to not religious
        TRUE ~ "Not religious"),
        levels = c("Protestant", "Jewish", "Catholic","Not religious", "Other")),
    
    recode_attndch_4way = factor(case_when(
        recode_attndch == "Once a week" ~ "Weekly",
        recode_attndch == "More than once a week" ~ "Weekly",
        recode_attndch == "Once or twice a month" ~ "Monthly",
        recode_attndch == "A few times a year" ~ "Yearly",
        recode_attndch == "Seldom" ~ "Never",
        recode_attndch == "Never" ~ "Never",
        recode_attndch == "Don't know" ~ "Never",
        #the 35 who refused say never go
        TRUE ~ "Never"),
        levels = c("Weekly","Monthly", "Yearly", "Never"))
)




# save
#9/12 updated to add new variables and fix name of educ_3way to wh and add real 3way
saveRDS(cces_2016, "data/cces_new.rds")
saveRDS(pew, "data/pew_new.rds")



################################################################################
#### ELECTION RESULTS ##########################################################
################################################################################

pres_env <- load_to_env("data/1976-2016-president.RData")

election_2016 <- pres_env$x %>%
    filter(year == 2016) %>%
    select(state, party, candidatevotes, totalvotes) %>% 
    group_by(state) %>%
    summarize(demtotal = sum(candidatevotes[party == "democrat"], na.rm = TRUE),
              reptotal = sum(candidatevotes[party == "republican"],
                             na.rm = TRUE),
              totalvotes = unique(totalvotes)) %>%
    ungroup()

#saveRDS(election_2016, "data/election.rds")


############ graveyard (temp)

# income_model2 <- polr(recode_income_5way ~ recode_female + recode_race + recode_region + 
#                           recode_pid_3way + recode_educ + child + recode_age,
#                       data = pew)
# 
# unique( predict(income_model2,, type = "class"))
# sum(predict(income_model2,, type = "class") == na.omit(pew$recode_income_5way))/sum(!is.na(pew$recode_income_5way))
# 
# income_model3 <- polr(recode_income_4way ~ recode_female + recode_race + recode_region + 
#                           recode_pid_3way + recode_educ + child + recode_age,
#                       data = pew)
# 
# unique( predict(income_model3,, type = "class"))
# sum(predict(income_model3,, type = "class") == na.omit(pew$recode_income_4way))/sum(!is.na(pew$recode_income_4way))
# 
# 
# income_model4 <- polr(recode_income_3way ~ recode_female + recode_race + recode_region + 
#                           recode_pid_3way + recode_educ + child + recode_age,
#                       data = pew)
# 
# unique( predict(income_model4,, type = "class"))
# sum(predict(income_model4,, type = "class") == na.omit(pew$recode_income_3way))/sum(!is.na(pew$recode_income_3way))
# 
# 
# no_match_3way <- predict(income_model4,, type = "class") != na.omit(pew$recode_income_3way)
# 
# 
# comp<- cbind(as.character(predict(income_model4,, type = "class")[no_match_3way]), as.character(na.omit(pew$recode_income_3way)[no_match_3way]))
# 
# 
# #alright so then how does education do?
# educ_model
# sum(predict(educ_model,, type = "class") == na.omit(pew$recode_educ))/sum(!is.na(pew$recode_educ))
# 
# 
# pew <- pew %>% mutate(
#     recode_income == is.na(recode_income) ~ predict(income_model,  newdata =., type = "class")
# )

#### cces 
##### 5 way
# income_model2_cces <- polr(recode_income_5way ~ recode_female + recode_race +
#                                recode_region + recode_pid_3way + recode_educ +
#                                child18 + recode_age,
#                            data = cces_2016)
# 
# unique( predict(income_model2_cces,, type = "class"))
# sum(predict(income_model2_cces,, type = "class") == cces_2016[!is.na(cces_2016$recode_income_5way) & !is.na(cces_2016$child18)  ,"recode_income_5way"] %>% pull())/sum(!is.na(cces_2016$recode_income_5way) & !is.na(cces_2016$child18) )
# 
# 
# ### 4 way
# income_model3_cces <- polr(recode_income_4way ~ recode_female + recode_race + 
#                                recode_region + recode_pid_3way + recode_educ +
#                                child18 + recode_age,
#                            data = cces_2016)
# 
# unique( predict(income_model3_cces,, type = "class"))
# sum(predict(income_model3_cces,, type = "class") == cces_2016[!is.na(cces_2016$recode_income_4way) & !is.na(cces_2016$child18)  ,"recode_income_4way"] %>% pull())/sum(!is.na(cces_2016$recode_income_4way) & !is.na(cces_2016$child18) )
# 
# 
# #### 3 way
# income_model4_cces <- polr(recode_income_3way ~ recode_female + recode_race +
#                                recode_region +  recode_pid_3way + recode_educ + 
#                                child18 + recode_age,
#                            data = cces_2016)
# 
# unique( predict(income_model4_cces,, type = "class"))
# sum(predict(income_model4_cces,, type = "class") == cces_2016[!is.na(cces_2016$recode_income_4way) & !is.na(cces_2016$child18)  ,"recode_income_3way"] %>% pull())/sum(!is.na(cces_2016$recode_income_4way) & !is.na(cces_2016$child18) )
# 
# 
# 
# 
# 
# #### tables
# pew_new <- pew %>% dplyr::select(recode_relig, recode_born, recode_attndch, recode_income)
# cces_new <- cces_2016 %>% dplyr::select(recode_relig, recode_born, recode_attndch, recode_income)
# 
# relig <- cbind(pew_new %>% group_by(recode_relig) %>% summarise(n()),
#                cces_new  %>% group_by(recode_relig) %>% summarise(n()))
# 
# born <- cbind(pew_new %>% group_by(recode_born) %>% summarise(n()),
#               cces_new  %>% group_by(recode_born) %>% summarise(n()))
# 
# 
# church <- cbind(rbind(pew_new %>% group_by(recode_attndch) %>% summarise(n()), c(NA, 0)),
#                 cces_new  %>% group_by(recode_attndch) %>% summarise(n()))
# 
# income <- cbind(pew_new %>% group_by(recode_income) %>% summarise(n()),
#                 cces_new  %>% group_by(recode_income) %>% summarise(n()))
# 
# 
# kable(relig[,c(1,2,4)],
#       format = "latex",
#       caption = "Religion",
#       col.names = c("Level", "Pew",
#                     "CCES"),
#       booktabs = T, digits= 3) %>%
#     kable_styling(position = "center", latex_options = "hold_position")
# 
# kable(born[,c(1,2,4)],
#       format = "latex",
#       caption = "Born-Again",
#       col.names = c("Level", "Pew",
#                     "CCES"),
#       booktabs = T, digits= 3) %>%
#     kable_styling(position = "center", latex_options = "hold_position")
# 
# kable(church[,c(1,2,4)],
#       format = "latex",
#       caption = "Church Attendance",
#       col.names = c("Level", "Pew",
#                     "CCES"),
#       booktabs = T, digits= 3) %>%
#     kable_styling(position = "center", latex_options = "hold_position")
# 
# 
# kable(income[,c(1,2,4)],
#       format = "latex",
#       caption = "Income",
#       col.names = c("Level", "Pew",
#                     "CCES"),
#       booktabs = T, digits= 3) %>%
#     kable_styling(position = "center", latex_options = "hold_position")
# 
# 
# 
# #uncomment if want don't knows to become NAs
# # levels(cces_2016$recode_attndch)[levels(cces_2016$recode_attndch) == "Don't know"] <- NA_character_
# 
# 
# ## Digging into pew$child a bit
# 
# #Children: pew$child (may need HH1 and HH3), cces_2016$child18 and cces_2016$child18num
# #hh1 = how many people including yourself live in your household
# #hh3 = how many including yourself are 18+
# unique(pew$child)
# hh1 <- as.character(pew$hh1)
# hh1 <- ifelse(grepl("Don't",hh1),"-99", hh1)
# hh1 <- ifelse(grepl("8 or",hh1),"8", hh1) #8 or more to just 8
# hh1 <- as.numeric(hh1)
# 
# hh3 <- as.character(pew$hh3)
# hh3 <- ifelse(grepl("Don't",hh3),"-99", hh3)
# hh3 <- as.numeric(hh3)
# 
# child <- pew$child 
# child <- ifelse(grepl("No", child), 0, child)
# #oops not looking grat
# child == hh1-hh3
# #hmm appears to not be the same thing?
# child[1:20]
# (hh1-hh3)[1:20]
# #it could could be that someone else's kids live in the same house such that hh1-hh3 > child
# sum(na.omit(child < hh1-hh3))
# #generally appears not to be true so probably child = how many children total and not limiting to under 18s :(


#NB: some additional info in pewrelig_t, but kind of a weird group who answers so not sure who its a follow up for (not all protestants just some, not all something else's just some)
# cces_2016 %>% filter(religpew_t != "__NA__") %>% group_by(religpew) %>% summarise(n())


