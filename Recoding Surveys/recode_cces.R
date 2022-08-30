#### CCES Recoding (Cleaned March 2021)
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
#### CCES #######################################################################
################################################################################

load_to_env <- function (rdata, env = new.env()) {
    load(rdata, env)
    return(env)
}
path_dat = "/Users/Ciara/Dropbox/kpop/application/data/"
cces_env <- load_to_env(paste0(path_dat, "CCES/CCES16_Common_Content.RData"))


# Common Content Weight -- commonweight_vv_post (post wave with vote validation)
#keeping only those weighted
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
        #this still includes 81 people who claim to have both definitly voted and on CC16_401 and didn't vote on CC16_410a; interp as didnt vote for pres only
        #360 NAs created: 229 = I'm not sure, 50 = NA, 81 = I didn't vote in this election
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
                                 #Other, Not Sure, Independent , NA
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
    recode_inputstate = inputstate_post
)

cces_2016$recode_inputstate <- droplevels(cces_2016$recode_inputstate)
cces_2016$recode_educ <- droplevels(cces_2016$recode_educ)


#### New Variables
cces_2016 <- cces_2016 %>% mutate(
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

    recode_born = case_when(pew_bornagain == "No" ~ "No",
                            pew_bornagain == "Yes" ~ "Yes", 
                            #NAs (43) assume to be not
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
    
    recode_income_6way = factor(case_when(
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
    
    recode_income_5way = factor(case_when(
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
    
    recode_income_4way = factor(case_when(
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

# cleaning up new vars
cces_2016 <- cces_2016 %>% mutate(
    recode_relig_5way = factor(case_when(
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
        #there's a lot of something else (2315) 
        #quite of few of them are v religios going to church weekly (300)
        #when asked what type of protestant many identified so call them protestant
        #if answered, call prot (includes 71 who answered somethign else to religpew_protestant)
        recode_relig == "Something else" & !is.na(religpew_protestant) ~ "Protestant",
        #if no answer other (1785)
        recode_relig == "Something else" & is.na(religpew_protestant) ~ "Other",
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

## Make some interactions 
#edited in March 2021 to add more interactions: XXX should update to final interactions used
cces_2016 <- cces_2016 %>%
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







# save
#9/12 updated to add new variables and fix name of educ_3way to wh and add real 3way in cces_nwe
#march 2021 updated to include more interactions, fix Xway names, and importantly caught typo in relig_5way something else & is.na = other (prev was soemthing else $!is.na so 530 protestants were called other)
nrow(cces_2016)
saveRDS(cces_2016, paste0(path_dat, "cces_clean.rds"))
