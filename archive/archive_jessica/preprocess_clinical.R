library(tidyverse)
library(lubridate)

# parse all date columns as date datatype
clinical <- 
  read_csv("../clinical.csv") %>%
  mutate(date_dx=date(ymd(`Dx Date`)),
         date_last_fu=date(ymd(`Last F/U`)),
         date_death=date(ymd(`Date of Death`)),
         date_lf=date(ymd(`Date of LF`)),
         date_rf=date(ymd(`Date RF`)),
         date_dm=date(ymd(`Date DM`)))

# fill missing last follow-up dates with date of death
clinical <- clinical %>%
  mutate(date_last_fu=coalesce(date_last_fu, date_death))
  
clinical <- clinical %>%
  rowwise() %>%
  mutate(date_lrf=min(date_lf,
                      date_rf,
                      date_last_fu,
                      na.rm=T),
         date_dfs=min(date_lrf,
                      date_death,
                      date_last_fu,
                      na.rm=T)) %>%
  ungroup()

# fill missing last follow-up dates with date of death
clinical <- clinical %>%
  mutate(date_last_fu=coalesce(date_last_fu, date_death))

# give event columns consistent names
clinical <- clinical %>% 
  mutate(event_os=`Vital Status (0=alive, 1=dead)`,
         event_lf=`Local Failure`,
         event_rf=`RF (0=np, 1=yes)`,
         event_lrf=LRF,
         event_dm=`DM (0=no, 1=yes)`,
         event_dfs=DFS,
         patient_id=`Study Number`)

# compute time to death/failure
clinical <- clinical %>% 
  mutate_at(vars(contains("date_"), -date_dx, -date_last_fu),
            .funs=list(time=~ifelse(is.na(.),
                                    date_last_fu - date_dx,
                                    . - date_dx))) %>%
  rename_at(vars(contains("time")), list(~gsub("date_", "", .)))

clinical <- clinical %>%
  mutate(`Concurrent chemo (1=yes, 0=no)`=ifelse(`Concurrent chemo (1=yes, 0=no)`=="Yes", 1, 0),
         SEX_F=ifelse(SEX=="Female", 1, 0),
         cT=as.double(na_if(cT, "x")),
         `HIV (0=no, 1=yes)`=ifelse(`HIV (0=no, 1=yes)`=="luk", 0, `HIV (0=no, 1=yes)`))

clinical %>%
  filter(`Study Number` != 1018) %>%
  write_csv("clinical_clean.csv")
