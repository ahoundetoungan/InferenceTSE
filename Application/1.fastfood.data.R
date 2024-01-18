#' This file prepares the data set to be used.
#' The data set will be saved under the name "DataFastFood.2Stage.rda" and will be used in other files.
#' It also computes descriptive statistics

rm(list = ls())
library(haven)
library(dplyr)
library(glue)
library(PartialNetwork)

## Add your working directory
proot <- c("~/AHdata/FastFood")
root  <- sapply(proot, dir.exists)
root  <- proot[root][1]
setwd(root)

filname    <- "DataFastFood.2Stage.rda"  
if (file.exists(filname)) {
  load(file = filname)
} else {
  schinfo  <- read_xpt("School Files/School Information Data/Schinfo.xpt")
  # SAT_SCHL All students at this school were selected for the In-Home Interview1 1 = yes
  # SCID school ID
  wave1    <- read_xpt("In Home Interview Files/allwave1.xpt")
  wave2    <- read_xpt("In Home Interview Files/wave2.xpt")
  hfriend2 <- read_xpt("Friend Files/Wave II In-Home Nominations/hfriend2.xpt")
  
  # All ID with school and sister school
  allID    <- wave2 %>% select(c("AID", "SCID2", "SSCID2"))
  
  # Saturated schools
  schid    <- unlist(schinfo %>% filter(SAT_SCHL == 1) %>% select("SCID"))
  datafull <- wave2 %>% left_join(wave1, by = "AID") %>% left_join(hfriend2, by = "AID") %>%
    filter(SCID2 %in% schid)
  
  # School ID, Gender, Age, Hispanic, White, Black, Asian, Grade, Fast food consumption, 
  # Weekly allowance, lives with father, lives with mother
  datafull <- datafull %>% mutate(SCID = factor(SCID2), Male = ifelse(BIO_SEX2 == 1, 1, 0), Female = 1 - Male, Age = CALCAGE2, 
                                  Hispanic = ifelse(is.na(H1GI4), 0, ifelse(H1GI4 == 1, 1, 0)), White = ifelse(is.na(H1GI6A), 0, ifelse(H1GI6A == 1, 1, 0)),
                                  Black = ifelse(is.na(H1GI6B), 0, ifelse(H1GI6B == 1 & White == 0, 1, 0)), 
                                  Asian = ifelse(is.na(H1GI6D), 0, ifelse(H1GI6D == 1 & White == 0 & Black == 0, 1, 0)),
                                  Otherrace = 1 - White - Black - Asian, Grade = H2GI9, Grade7_8 = ifelse(H2GI9 %in% 7:8, 1, 0), Grade9_10 = ifelse(H2GI9 %in% 9:10, 1, 0),
                                  Grade11_12 = ifelse(H2GI9 %in% 11:12, 1, 0), Fastfood = H2NU77, Wkyallowance = ifelse(H2EE8 < 96, H2EE8, NA),
                                  WithFather = apply(datafull %>% select(starts_with("H2HR4")) == 11, 1, any),
                                  WithMother = apply(datafull %>% select(starts_with("H2HR4")) == 14, 1, any),
                                  WithBoth = ifelse((WithFather == 1) & (WithMother == 1), 1, 0),
                                  MJobProf = ifelse(is.na(H2RM4), 0, ifelse(H2RM4 %in% 1:2, 1, 0)), MJobNone = ifelse(is.na(H2RM4), 0, ifelse(H2RM4 == 16, 1, 0)), 
                                  MJobMiss = ifelse(is.na(H2RM4), 1, ifelse(H2RM4 > 16, 1, 0)), MJobOhter = 1 - MJobProf - MJobNone - MJobMiss, 
                                  FJobProf = ifelse(is.na(H2RF4 ), 0, ifelse(H2RF4  %in% 1:2, 1, 0)), FJobNone = ifelse(is.na(H2RF4 ), 0, ifelse(H2RF4  == 16, 1, 0)), 
                                  FJobMiss = ifelse(is.na(H2RF4 ), 1, ifelse(H2RF4  > 16, 1, 0)), FJobOhter = 1 - FJobProf - FJobNone - FJobMiss)
  
  # Mother education
  # lower than High 1, 2, 10
  # High 3, 4, 5
  # Some college 6, 7
  # Graduate from college or university 8, 9
  # Missing >= 11
  # H1NM4 non resident biological mother is Wave 1
  # H1RM1 resident mother wave 1 
  # H2RM1 resident mother in Wave 2
  # If H2HR7A--H2HR7R is 7 then lives with biological mother
  # # If H1HR6A--H1HR6K is 7 then Mother education is H1RM1
  # # Otherwise Mother education is H1NM4
  # Otherwise Mother education is H2RM1
  # Do the same for father education
  datafull <- datafull %>% mutate(WithBioMo1 = apply(datafull %>% select(starts_with("H1HR6")) == 7, 1, any),
                                  WithBioMo1 = ifelse(is.na(WithBioMo1), 0, WithBioMo1),
                                  WithBioMo2 = apply(datafull %>% select(starts_with("H2HR7")) == 7, 1, any),
                                  WithBioMo2 = ifelse(is.na(WithBioMo2), 0, WithBioMo2),
                                  WithBioFa1 = apply(datafull %>% select(starts_with("H1HR6")) == 1, 1, any),
                                  WithBioFa1 = ifelse(is.na(WithBioFa1), 0, WithBioFa1),
                                  WithBioFa2 = apply(datafull %>% select(starts_with("H2HR7")) == 1, 1, any),
                                  WithBioFa2 = ifelse(is.na(WithBioFa2), 0, WithBioFa2),
                                  MoEduc     = ifelse(WithBioMo2, ifelse(WithBioMo1, H1RM1, H1NM4), H2RM1),
                                  FaEduc     = ifelse(WithBioFa2, ifelse(WithBioFa1, H1RF1, H1NF4), H2RF1),
                                  MoEduLHigh = ifelse(is.na(MoEduc), 0, ifelse(MoEduc %in% c(1, 2, 10), 1, 0)),
                                  MoEduHigh  = ifelse(is.na(MoEduc), 0, ifelse(MoEduc %in% 3:5, 1, 0)),
                                  MoEduSomCO = ifelse(is.na(MoEduc), 0, ifelse(MoEduc %in% 6:7, 1, 0)),
                                  MoEduUniv  = ifelse(is.na(MoEduc), 0, ifelse(MoEduc %in% 8:9, 1, 0)),
                                  MoEduMiss  = ifelse(is.na(MoEduc), 1, ifelse(MoEduc >= 11, 1, 0)),
                                  FaEduLHigh = ifelse(is.na(FaEduc), 0, ifelse(FaEduc %in% c(1, 2, 10), 1, 0)),
                                  FaEduHigh  = ifelse(is.na(FaEduc), 0, ifelse(FaEduc %in% 3:5, 1, 0)),
                                  FaEduSomCO = ifelse(is.na(FaEduc), 0, ifelse(FaEduc %in% 6:7, 1, 0)),
                                  FaEduUniv  = ifelse(is.na(FaEduc), 0, ifelse(FaEduc %in% 8:9, 1, 0)),
                                  FaEduMiss = ifelse(is.na(FaEduc), 1, ifelse(FaEduc >= 11, 1, 0))) %>% 
    select(!c("MoEduc", "FaEduc")) 
  
  # Network
  # 55555555 nominated friend was also nominated as one of the partners, this is a missing value
  # 77777777 nominated friend doesn’t go to sister or sample school, this is not a missing value
  # 88888888 nominated friend goes to sister school- not on list, this is not a missing value
  # 99999999 nominated friend goes to sample school- not on list, this is not a missing value
  datafull <- datafull %>%
    mutate(NMF    = colSums(apply(datafull %>% select(starts_with("MF_AID")), 1, nchar) == 8), #number of male friends (declared)
           NFF    = colSums(apply(datafull %>% select(starts_with("FF_AID")), 1, nchar) == 8),
           MMISS5 = rowSums((datafull %>% select(starts_with("MF_AID"))) == "55555555"),
           MMISS7 = rowSums((datafull %>% select(starts_with("MF_AID"))) == "77777777"),
           MMISS8 = rowSums((datafull %>% select(starts_with("MF_AID"))) == "88888888"),
           MMISS9 = rowSums((datafull %>% select(starts_with("MF_AID"))) == "99999999"),
           FMISS5 = rowSums((datafull %>% select(starts_with("FF_AID"))) == "55555555"),
           FMISS7 = rowSums((datafull %>% select(starts_with("FF_AID"))) == "77777777"),
           FMISS8 = rowSums((datafull %>% select(starts_with("FF_AID"))) == "88888888"),
           FMISS9 = rowSums((datafull %>% select(starts_with("FF_AID"))) == "99999999"))  %>%                               # Replacing values
    mutate(across(starts_with(c("MF_AID", "FF_AID")), .fns = function(u) ifelse(u == AID, "", u)),
           across(starts_with(c("MF_AID", "FF_AID")), .fns = function(u) ifelse(u %in% c("55555555"), "", u)),
           across(starts_with(c("MF_AID", "FF_AID")), .fns = function(u) ifelse(u %in% c("77777777"), "", u)),
           across(starts_with(c("MF_AID", "FF_AID")), .fns = function(u) ifelse(u %in% c("88888888"), "", u)),
           across(starts_with(c("MF_AID", "FF_AID")), .fns = function(u) ifelse(u %in% c("99999999"), "", u))) %>% 
    arrange(SCID, AID) %>% filter(!is.na(Wkyallowance), Fastfood < 8) 

  # number of school
  SCH      <- datafull %>% filter(!duplicated(paste0(SCID2, "_", SSCID2))) %>% select(c("SCID2", "SSCID2"))
  nsc      <- nrow(SCH)
  
  # Construct the adjacency matrix
  GM       <- list()
  GF       <- list()
  for (k in 1:nsc) {
    sc         <- SCH$SCID2[k]
    ssc        <- SCH$SSCID2[k]
    datIDk     <- datafull %>% filter(SCID2 == sc) %>% select(c("AID", "Male", starts_with(c("MF_AID", "FF_AID"))))
    datk       <- datIDk %>% select(!c("AID", "Male"))
    Nk         <- nrow(datk)
    GMmisk     <- matrix(rep(datIDk$Male, Nk), Nk, Nk, byrow = TRUE)
    GFmisk     <- matrix(rep(1 - datIDk$Male, Nk), Nk, Nk, byrow = TRUE)
    GMk        <- t(apply(datk, 1, function(u) unlist(sapply(datIDk$AID, function(v) v %in% u))))*1
    GF[[k]]    <- GMk*GFmisk
    GM[[k]]    <- GMk*GMmisk
  }
  
  datafull     <- datafull %>% mutate(MFr = unlist(lapply(GM, function(x) rowSums(x))),
                                      FFr = unlist(lapply(GF, function(x) rowSums(x)))) %>%
    group_by(SCID) %>% mutate(scsize = n()) %>% ungroup()
  scsize       <- unlist(datafull %>% filter(!duplicated(SCID)) %>% select(scsize))
  exp.var      <- c("Female", "Age", "Hispanic", "Black", "Asian", "Otherrace", "Grade7_8", "Grade9_10", 
                 "Grade11_12", "Wkyallowance", "WithBoth", "MJobProf", "MJobMiss", "MJobOhter", 
                 "FJobProf", "FJobMiss", "FJobOhter", "MoEduLHigh", "MoEduSomCO", "MoEduUniv", 
                 "MoEduMiss", "FaEduLHigh", "FaEduSomCO", "FaEduUniv", "FaEduMiss")
  all.var      <- c("Fastfood", "Female", "Age", "Hispanic", "White", "Black", "Asian", "Otherrace", "Grade7_8", "Grade9_10", 
                    "Grade11_12", "Wkyallowance", "WithBoth", "MJobNone", "MJobProf", "MJobMiss", "MJobOhter", 
                    "FJobNone", "FJobProf", "FJobMiss", "FJobOhter", "MoEduHigh", "MoEduLHigh", "MoEduSomCO", "MoEduUniv", 
                    "MoEduMiss", "FaEduHigh", "FaEduLHigh", "FaEduSomCO", "FaEduUniv", "FaEduMiss")
  
  save(list = c("datafull", "GF", "all.var", "exp.var", "schinfo", "GM", "SCH", "nsc", "scsize", "schid"), 
       file = filname)
}

############################ Descriptive stat ##################################
datafull     <- datafull %>% mutate(Fr = MFr + FFr)

my.stat.des  <- function(x) {
  out <- c(mean(x), sd(x), min(x), max(x))
  names(out) <- c("Mean", "St. Dev.",   "Min",   "Max")
  out
} 
(sdesc       <- apply(datafull %>% select(!!all.var, MFr, FFr, Fr), 2, my.stat.des))
write.csv(t(sdesc), file = "sdesc.csv")
