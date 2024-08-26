library(dplyr)

adsl <- haven::read_sas("../../Lancaster/MATH590/AllProvidedFiles_310/PDS_DSA_20020408/adsl_pds2019.sas7bdat")
adrs <- haven::read_sas("../../Lancaster/MATH590/AllProvidedFiles_310/PDS_DSA_20020408/adrsp_pds2019.sas7bdat") %>%
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
  left_join(select(adsl, SUBJID, DTHDYX), by = "SUBJID") %>%
  group_by(SUBJID) %>%
  filter(all(RSRESP != "Progressive disease")) %>%
  mutate(assessment_post_censoring = max(as.numeric(VISITDY >= DTHDYX &
                                                      (RSRESP != "Unable to evaluate" | RSRESP != "Unknown")))) %>% # if they're known to have not progressed (what about subject 7?)
  filter(VISITDY <= DTHDYX) %>%
  mutate(LASTDY = lag(VISITDY), LASTVISIT = lag(VISIT)) %>%
  arrange(desc(VISITDY)) %>%
  slice(1) %>%
  ungroup() %>%
  rowwise() %>%
  mutate(
    VISITDY = if_else(is.na(LASTDY), DTHDYX, VISITDY),
    # If only screening visit use right censoring time
    LASTDY = max(0, LASTDY),
    # Truncate negative values (i.e screening time)
    LASTDY = if_else(is.na(LASTDY), 0, LASTDY),
    delta0 = assessment_post_censoring # how can we be sure of no known progression?
  ) %>%
  ungroup() %>%
  select(SUBJID, L = LASTDY, R = VISITDY, delta0)

# Alternative approach, if last response was stable/partial and left within 3 months then known to not?
find_delta0_2 <- adrs %>%
  left_join(select(adsl, SUBJID, DTHDYX, DTHX), by = "SUBJID") %>%
  group_by(SUBJID) %>%
  filter(all(RSRESP != "Progressive disease")) %>%
  mutate(LASTDY = lag(VISITDY), LASTVISIT = lag(VISIT)) %>%
  ungroup() %>%
  filter(RSRESP == "Stable disease") %>%
  arrange(SUBJID, desc(VISITDY)) %>%
  distinct(SUBJID, .keep_all = T) %>%
  mutate(dy_diff = DTHDYX - VISITDY) %>%
  filter(dy_diff < 92) %>%
  mutate(delta0_tmp = 1) %>%
  select(SUBJID, delta0_tmp)

find_ints <- bind_rows(find_delta1, find_delta0) %>%
  mutate(
    delta1 = if_else(is.na(delta1), 0, delta1),
    delta0 = if_else(is.na(delta0), 0, delta0),
  )

# dual censored data
addc <- adsl %>%
  select(SUBJID, TRT, ATRT, V = DTHDYX, delta2 = DTHX) %>%
  left_join(find_ints, "SUBJID") %>%
  left_join(find_delta0_2, "SUBJID") %>%
  mutate(
    R = if_else(R > V, V, R),
    R = if_else(is.na(R), V, R), # if no right value, set to censoring
    L = if_else(is.na(L), 0, L), # if no left value, set to 0
  ) %>%
  mutate(
    delta1 = if_else(is.na(delta1), 0, delta1),
    delta0 = if_else(!is.na(delta0_tmp), 1, if_else(!is.na(delta0), delta0, 0)),
    ATRTN = case_match(ATRT, "Best supportive care" ~ 0, .default = 1)
  ) %>%
  select(-delta0_tmp)

AMGEN_20020408 <- addc %>%
  mutate(L = L / 365.25,
         R = R / 365.25,
         V = V / 365.25)

usethis::use_data(AMGEN_20020408, overwrite = TRUE)
