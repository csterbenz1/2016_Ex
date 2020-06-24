library(tidyverse)
library(survey)
library(srvyr)

### functions
load_to_env <- function (rdata, env = new.env()) {
    load(rdata, env)
    return(env)
}

### definitions
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

### load
pew_env <- load_to_env("data/Oct16 public.RData")

## recode and filter
pew <- pew_env$pew_201610 %>%
    mutate(recode_vote_2016 = case_when(str_detect(Q11HORSE2, "Trump") ~ "Republican",
                                                  str_detect(Q11HORSE2, "Clinton") ~ "Democrat",
                                                  is.na(Q11HORSE2) ~ NA_character_,
                                                  TRUE ~ "Other"
                                                  ))


pew <- pew %>% mutate(
  # age
  recode_age = ifelse(age == "Don't know/Refused (VOL.)", NA, age),
  recode_age_bucket = factor(case_when( recode_age <= 35 ~ "18 to 35",
                          recode_age <= 50 ~ "36 to 50",
                          recode_age < 65 ~ "51 to 64",
                          !is.na(recode_age) ~ "65+",
                          TRUE ~ "51 to 64"), levels = c("18 to 35", "36 to 50", "51 to 64", "65+")),
  
  recode_age_3way = case_when( recode_age <= 50 ~ "a_18to50",
                          recode_age < 65 ~ "b_51to64",
                          !is.na(recode_age) ~ "c_65",
                          TRUE ~ "b_51to64"),
  
  # gender
  recode_female = case_when(sex == "Female" ~ "Female",
                            TRUE ~ "Male"),
  
  # race
  recode_race = case_when(racethn == "White, non-Hisp" ~ "White",
                          racethn == "Black, non-Hisp" ~ "Black",
                          racethn == "Hispanic" ~ "Hispanic",
                          TRUE ~ "Other"),
  
  # region
  recode_region = case_when(state %in% northeast_states ~ "Northeast",
                            state %in% west_states ~ "West",
                            state %in% midwest_states ~ "Midwest",
                            state %in% south_states ~ "South",
                            TRUE ~ "South"),
  
  # party
  recode_pid_3way = case_when( party == "Democrat" ~ "Dem",
                        party == "Republican" ~ "Rep",
                        TRUE ~ "Ind"),
  
  # education
  recode_educ = factor(case_when(educ2 == "Less than high school (Grades 1-8 or no formal schooling)" ~ "No HS",
                          educ2 == "High school incomplete (Grades 9-11 or Grade 12 with NO diploma)" ~ "No HS",
                          educ2 == "High school graduate (Grade 12 with diploma or GED certificate)" ~ "High school graduate",
                          educ2 == "Some college, no degree (includes some community college)" ~ "Some college",
                          educ2 == "Two year associate degree from a college or university" ~ "2-year",
                          educ2 == "Four year college or university degree/Bachelor's degree (e.g., BS, BA, AB)" ~ "4-year",
                          educ2 == "Some postgraduate or professional schooling, no postgraduate degree" ~ "Post-grad",
                          educ2 == "Postgraduate or professional degree, including master's, doctorate, medical or law degree" ~ "Post-grad",
                          TRUE ~ "Some college"), levels = c("No HS", "High school graduate", "Some college", "2-year", "4-year", "Post-grad")),
  
  recode_educ_3way = factor(case_when(recode_race != "White" ~ "No Split",
                                      recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
                                      recode_educ %in% c("2-year", "4-year") ~ "College",
                                      TRUE ~ "Post-grad"), levels = c("No Split", "No College", "College", "Post-grad")),
  
  # state
  recode_inputstate = state,
  
  # interactions
  recode_race_weduc = as.factor(case_when(
    recode_race == "White" ~ paste(recode_race, recode_educ),
    TRUE ~ recode_race)),
  recode_race_educ = as.factor(paste(recode_race, recode_educ)),
  recode_race_educ_reg = as.factor(paste(recode_region, case_when(
    recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
    TRUE ~ recode_race), sep = ", ")),
  recode_race_educ_reg_female = as.factor(paste(recode_female, paste(recode_region, case_when(
    recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
    TRUE ~ recode_race), sep = ", "), sep = ", ")),
  recode_race_college_state = as.factor(paste(recode_inputstate, case_when(
    recode_race == "White" ~ paste(recode_race, case_when(recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
                                                          TRUE ~ "College")),
    TRUE ~ recode_race))),
  recode_age_reg = paste(recode_age_bucket, recode_region, sep = ", "),
  recode_female_reg = paste(recode_female, recode_region, sep = ", ")
)

## Save
saveRDS(pew, "data/pew_full.rds")

pew <- pew %>%
  filter(plan1 %in% c("Plan to vote", "Already voted", "Don't know/Refused (VOL.)"))
saveRDS(pew, "data/pew.rds")

################################################################################
#### CCES ######################################################################
################################################################################

cces_env <- load_to_env("./data/CCES16_Common_Content.RData")

## make outcome

# Horserace (post)
cces_2016 <- cces_env$cces_2016 %>%
    mutate( recode_vote_2016 = case_when(
                str_detect(CC16_410a, "Democrat") ~ "Democrat",
                str_detect(CC16_410a, "Republican") ~ "Republican",
                as.numeric(CC16_410a) < 6 ~ "Other",
                TRUE ~ NA_character_) )

# Common Content Weight -- commonweight_vv_post (post wave with vote validation)

## demographic recodes
cces_2016$inputstate_post = droplevels(cces_2016$inputstate_post)

cces_2016 <- cces_2016 %>% mutate(
  recode_id = 1:nrow(cces_2016),
  
  # age
  recode_age = 2012 - as.numeric(as.character(birthyr)),
  recode_age_bucket = factor(case_when( recode_age <= 35 ~ "18 to 35",
                          recode_age <= 50 ~ "36 to 50",
                          recode_age < 65 ~ "51 to 64",
                          !is.na(recode_age) ~ "65+",
                          TRUE ~ "51 to 64"), levels = c("18 to 35", "36 to 50", "51 to 64", "65+")),
  
  recode_age_3way = case_when( recode_age <= 50 ~ "a_18to50",
                          recode_age < 65 ~ "b_51to64",
                          !is.na(recode_age) ~ "c_65",
                          TRUE ~ "b_51to64"),
  
  # gender
  recode_female = case_when(gender == "Female" ~ "Female",
                            TRUE ~ "Male"),
  
  # race
  recode_race = ifelse(race %in% c("White", "Black", "Hispanic"), as.character(race), "Other"),
  
  # region
  recode_region = case_when(inputstate_post %in% northeast_states ~ "Northeast",
                            inputstate_post %in% west_states ~ "West",
                            inputstate_post %in% midwest_states ~ "Midwest",
                            inputstate_post %in% south_states ~ "South",
                            TRUE ~ "South"),
  
  # party
  recode_pid_3way = case_when( str_detect(pid3, "Democrat") ~ "Dem",
                        str_detect(pid3, "Republican") ~ "Rep",
                        TRUE ~ "Ind"),
  
  

  # educ
  recode_educ = factor(educ, levels = c("No HS", "High school graduate", "Some college", "2-year", "4-year", "Post-grad")),
  
  recode_educ_3way = factor(case_when(recode_race != "White" ~ "No Split",
                                      recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
                                      recode_educ %in% c("2-year", "4-year") ~ "College",
                                      TRUE ~ "Post-grad"), levels = c("No Split", "No College", "College", "Post-grad")),
  
  # state
  recode_inputstate = inputstate_post,
  
  # interactions
  recode_race_weduc = as.factor(case_when(
    recode_race == "White" ~ paste(recode_race, recode_educ),
    TRUE ~ recode_race)),
  recode_race_educ = as.factor(paste(recode_race, recode_educ)),
  recode_race_educ_reg = as.factor(paste(recode_region, case_when(
    recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
    TRUE ~ recode_race), sep = ", ")),
  recode_race_educ_reg_female = as.factor(paste(recode_female, paste(recode_region, case_when(
    recode_race == "White" ~ paste(recode_race, recode_educ, sep = ", "),
    TRUE ~ recode_race), sep = ", "), sep = ", ")),
  recode_race_college_state = as.factor(paste(recode_inputstate, case_when(
    recode_race == "White" ~ paste(recode_race, case_when(recode_educ %in% c("No HS", "High school graduate", "Some college") ~ "No College",
                                                          TRUE ~ "College")),
    TRUE ~ recode_race))),
  recode_age_reg = paste(recode_age_bucket, recode_region, sep = ", "),
  recode_female_reg = paste(recode_female, recode_region, sep = ", ")
)

cces_2016$recode_inputstate <- droplevels(cces_2016$recode_inputstate)
cces_2016$recode_educ <- droplevels(cces_2016$recode_educ)

## save

saveRDS(cces_2016, "data/cces.rds")

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

saveRDS(election_2016, "data/election.rds")
