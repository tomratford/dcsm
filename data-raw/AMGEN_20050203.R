## code to prepare `AMGEN_20050203` dataset goes here

library(dplyr)

adsl <- haven::read_sas("../../Lancaster/MATH590/AllProvidedFiles_309/PDS_DSA_20050203/adsl_pds2019.sas7bdat")
adrs <- haven::read_sas("../../Lancaster/MATH590/AllProvidedFiles_309/PDS_DSA_20050203/adrsp_pds2019.sas7bdat") %>%
  filter(RSREADER == "Oncologist") # use oncologist only for simplicity, in reality should impute with BICR rules

# find all who have known progression
find_delta1 <- adrs %>%
  group_by(SUBJID) %>%
  filter(any(RSRESP == "Progressive disease")) %>%
  mutate(
    LASTDY = lag(VISITDY),
    LASTVISIT = lag(VISIT)
  ) %>%
  arrange(desc(VISITDY)) %>%
  slice(1) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    LASTDY = max(0, LASTDY), # Truncate negative values
    delta1 = 1
  ) %>%
  ungroup() %>%
  select(SUBJID, L=LASTDY, R=VISITDY, delta1)

# Find all who have not progressed, and get last visit prior to their censoring/death
find_delta0 <- adrs %>%
  left_join(select(adsl, SUBJID, DTHDY), by = "SUBJID") %>%
  group_by(SUBJID) %>%
  filter(all(RSRESP != "Progressive disease")) %>%
  mutate(assessment_post_censoring = max(as.numeric(VISITDY >= DTHDY &
                                                      RSRESP != "Unable to evaluate"))) %>% # if they're known to have not progressed (what about subject 7?)
  filter(VISITDY <= DTHDY) %>%
  mutate(LASTDY = lag(VISITDY), LASTVISIT = lag(VISIT)) %>%
  arrange(desc(VISITDY)) %>%
  slice(1) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    VISITDY = if_else(is.na(LASTDY), DTHDY, VISITDY),
    # If only screening visit use right censoring time
    LASTDY = max(0, LASTDY),
    # Truncate negative values (i.e screening time)
    LASTDY = if_else(is.na(LASTDY), 0, LASTDY),
    delta0 = assessment_post_censoring # how can we be sure of no known progression?
  ) %>%
  ungroup() %>%
  select(SUBJID, L = LASTDY, R = VISITDY, delta0)

find_ints <- bind_rows(find_delta1, find_delta0) %>%
  mutate(
    delta1 = if_else(is.na(delta1), 0, delta1),
    delta0 = if_else(is.na(delta0), 0, delta0),
  )

# dual censored data
addc <- adsl %>%
  select(SUBJID, TRT, ATRT, V = DTHDY, delta2 = DTH) %>%
  left_join(find_ints, "SUBJID") %>%
  mutate(
    R = if_else(R > V, V, R),
    R = if_else(is.na(R), V, R), # if no right value, set to censoring
    L = if_else(is.na(L), 0, L), # if no left value, set to 0
    delta1 = if_else(is.na(delta1), 0, delta1),
    delta0 = if_else(is.na(delta0), 0, delta0),
    ATRTN = case_match(ATRT, "FOLFOX alone" ~ 0, .default = 1)
  )

# change to year and remove problematic subject
AMGEN_20050203 <- addc %>%
  mutate(L = L / 365.25,
         R = R / 365.25,
         V = V / 365.25) %>%
  filter(SUBJID != "000205")

usethis::use_data(AMGEN_20050203, overwrite = TRUE)
