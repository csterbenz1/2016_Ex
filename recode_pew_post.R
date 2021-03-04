#### Pew POST-Wave Recoding (Cleaned March 2021)
rm(list = ls())

library(MASS)
library(tidyverse)
library(survey)
library(srvyr)
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

################################################################################
#### Pew #######################################################################
################################################################################

### load original Pew .sav file
#late october pew
path_dat = "/Users/Ciara/Dropbox/kpop/application/data/"

#post election pew
pew <- read.spss(paste0(path_dat, "Pew/Nov16-Post-Election-public/Nov16 Post-Election public.sav"), to.data.frame=TRUE)

#FOR POST: to use the same cleaning code, we need to drop orig from the front of the var names
cleannames <- colnames(pew)
colnames(pew)[-22] <- gsub( "orig","", cleannames)[-22]

#q4 is reported vote choice
pew <- pew %>%
    mutate(recode_vote_2016 = 
               case_when(str_detect(q4, "Trump") ~ "Republican",
                         str_detect(q4, "Clinton") ~ "Democrat",
                         is.na(q4) ~ NA_character_,
                         TRUE ~ "Other" 
               ))

## Start by looking at missingness (fraction):
lapply(pew[, c("age", "sex", "racethn", "state", "party", "educ2")], 
       function(x) (sum(is.na(x)) + sum(grepl("Don't know/Refused", x))) / length(x))

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

pew <- pew %>% 
    mutate(missing = ifelse(is.na(recode_age) | is.na(recode_educ), 1, 0))

## modeling age to fill in missing
age_model <- glm(recode_age ~ recode_female + recode_race + recode_region + 
                     recode_pid_3way + recode_educ + child,
                 family = Gamma(),
                 data = pew)

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

#March 2021: updated names of income_Xway to count pref not to say in X, relig name to be 5 way (counted wrong before and had 6)
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
        #out of 41 of these: 11 said yes christ, 29 no, 1 refused
        relig == "Something else (SPECIFY:______)" ~ "Something else",
        #those who respond don't know/refuse  (50) 
        #out of these 50: 30 said yes christian, 11 said no, 9 refused
        relig == "Don't Know/Refused (VOL.)" & chr == "Yes" ~ "Protestant",
        #leaves only 20 (9+11) we'll call these NA
        TRUE ~ NA_character_),
        levels = c("Protestant", "Jewish", "Catholic", "Nothing in particular", 
                   "Agnostic", "Atheist", "Buddhist", "Muslim", "Mormon", "Hindu",
                   "Orthodox", "Something else")),
    
    #supposedly only asked born:if chr = yes or if relig = 1-4,13
    #(Protestant, Catholic, Mormon, Orthodox, Christian)
    #though looking at the NAs that's not totally true some in these groups still have NAs
    recode_born = case_when(born == "Yes, would" ~ "Yes",
                            born == "No, would not" ~ "No",
                            #59 refused but we think it's safe to say these people are not
                            born == "Don't know/Refused (VOL.)" ~ "No",
                            #NA's (769) these are mainly non-christian denominations
                            #so say safe to say they are not born again
                            TRUE ~ "No"),
    
    #the names are actually ok here
    #leaving the 35 that response don't know as they are 
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
    recode_income_6way = factor(case_when(
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
    
    recode_income_5way = factor(case_when(
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
    
    recode_income_4way = factor(case_when(
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

pew <- pew %>% mutate(
    recode_relig_5way = factor(case_when(
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
        #20 remaining NAs as Not religious (don't knows + refused to answer/no on chr)
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
        #assuming the 35 who refused say never go
        TRUE ~ "Never"),
        levels = c("Weekly","Monthly", "Yearly", "Never"))
)

## Make some interactions 
#edited in March 2021 to add more interactions: XXX should update to final interactions used
pew <- pew %>%
    mutate(
        # interactions
        recode_race_educ_reg = as.factor(paste(recode_region, case_when(
            recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
            TRUE ~ recode_race), sep = ", ")),
        #36 levels
        recode_educ_pid_race = as.factor(paste(recode_educ_3way,
                                               recode_pid_3way, 
                                               recode_race, sep = ", ")),
        #12 levels
        recode_pid_race = as.factor(paste(recode_pid_3way, 
                                          recode_race, sep = ", ")),
        #9 levels
        recode_educ_pid = as.factor(paste(recode_educ_3way,
                                          recode_pid_3way, sep = ", ")),
        recode_agesq = recode_age*recode_age,
        recode_agecubed = recode_age*recode_age*recode_age,
        #13 levles
        recode_midwest_edu_race = as.factor(
            case_when(recode_region == "Midwest" ~ paste(recode_region,
                                                         recode_race,
                                                         recode_educ_3way, 
                                                         sep = ", "),
                      TRUE ~ "No Split")), 
        #4 levels ( low educated whites in midwest only)
        recode_midwest_wh_edu = factor(case_when(
            (recode_race != "White" | recode_region != "Midwest") ~ "No Split", 
            TRUE ~ as.character(recode_educ_3way)), 
            levels = c("No Split", "No College", "College", "Post-grad")), 
        recode_age_factor = factor(case_when(recode_age <92 ~ as.character(recode_age), 
                                             TRUE ~ "92+"))
    )





#post-poll q1 is did you vote and if thy didnt they end the interview so we already subsetted to actual voters and now do not need to filter for did you vote or anythign 
nrow(pew)
#pew_post_CS_q4 done oct 2020
#pew_post_clean dont march 2021 (just added interactiosn and fixed Xway names)
saveRDS(pew, paste0(path_dat, "pew_post_clean.rds"))
