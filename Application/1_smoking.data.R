#' This file prepares the data set to be used.
#' The data set will be saved under the name "AHdata.rda" and will be used in other files.
#' It also computes descriptive statistics

##################### Libraries
rm(list = ls())
library(dplyr)
library(tidyr)
library(labelled)
library(readstata13)
library(ggplot2)
library(DescTools)
library(PartialNetwork)

##################### Data
filname        <- "AHdata.rda"  
if (file.exists(filname)) {
  load(file = filname)
} else {
  # Data
  mydata       <- read.dta13("w1s.dta")  # data from Stata
  mydata       <- mydata[order(mydata$sschlcde),]
  #In school there are only 3 girls but some students declare more than 3 female friends
  #In schools 93 and 193 either girls are declared as male friends or boys are declared as female friends
  mydata       <- mydata %>% filter(!(sschlcde %in% c(31, 93, 193))) %>% 
    group_by(sschlcde) %>%
    mutate(male        = replace_na(male, Mode(male, na.rm = TRUE)[1]),
           female      = 1 - male,
           hispanic    = replace_na(hispanic, Mode(hispanic, na.rm = TRUE)[1]),
           age         = replace_na(age, Mode(age, na.rm = TRUE)[1]),
           yearinschl  = replace_na(yearinschl, mean(yearinschl, na.rm = TRUE)),
           withbothpar = replace_na(withbothpar, Mode(withbothpar, na.rm = TRUE)[1]),
           smoke       = ifelse(cigarettes > 0, 1, 0),
           smoke       = replace_na(smoke, Mode(smoke, na.rm = TRUE)[1]),
           club        = ifelse(nclubs > 0, 1, 0),
           ID          = as.factor(aid)) %>% ungroup()
  
  mislist      <- c(55555555, 77777777, 88888888, 99999999, 99959995)
  mf_coln      <- paste0("mf", 1:5, "aid")
  ff_coln      <- paste0("ff", 1:5, "aid")
  f_coln       <- c(mf_coln, ff_coln)
  
  # mislist is an ID
  if (sum(mydata$aid %in% mislist) > 0) {
    stop("mislist is an ID")
  } else {
    cat("mislist is not an ID: OK", "\n")
  }
  
  # list of variable (excluding reference variables)
  va.all.names   <- c("male", "female", "age", "hispanic", "racewhite", "raceblack", "raceasian", "raceother", 
                      "withbothpar", "yearinschl", "mehigh", "melhigh", "memhigh", "memiss", "mjprof", 
                      "mjhome", "mjother", "mjmiss", "smoke")
  
  # list of variable (excluding reference variables for identification)
  va.names      <- c("female", "age", "hispanic", "raceblack", "raceasian", "raceother", 
                     "withbothpar", "yearinschl", "melhigh", "memhigh", "memiss", "mjprof", 
                     "mjother", "mjmiss", "smoke")
  
  # Are there NAs?
  apply(mydata[,va.names], 2, function(w) sum(is.na(w)))
  
  # remove friend from different groups
  # remove self friendship
  # remove friend non found
  N       <- nrow(mydata)
  dscf    <- rep(0,N)
  sfre    <- rep(0,N)
  nffr    <- rep(0,N)
  for (i in 1:N) {
    cat("Student ", i, "/", N, sep = "", "\n")
    for (j in f_coln) {
      k   <- which(mydata$aid == mydata[i, j, drop = TRUE])
      if (length(k) != 0) {
        # If the friend is found,
        # remove if different school
        if(mydata[i, "sschlcde", drop = TRUE] != mydata[k, "sschlcde", drop = TRUE]) {
          mydata[i, j, drop = TRUE]   <- -1
          dscf[i]                     <- dscf[i] + 1
        }
        # remove if self friendship
        if(mydata[i, "aid", drop = TRUE] == mydata[k, "aid", drop = TRUE]) {
          mydata[i, j, drop = TRUE]   <- -2
          sfre[i]                     <- sfre[i] + 1
        }
      }
      else {
        # If the friend is not found,
        if (!((mydata[i, j, drop = TRUE] %in% mislist) | is.na(mydata[i, j, drop = TRUE]))) {
          mydata[i, j]   <- -3
          nffr[i]        <- nffr[i] + 1
        }
      }
    }
  }
  
  cat("remove", sum(dscf), "link(s) because students from different schools: their code are recode as -1", "\n")
  cat("remove", sum(sfre), "self-friendship(s): their code are recoded as -2", "\n")
  cat("remove", sum(nffr), "non-found friends: their code are recoded as -3", "\n")
  rm(list = c("i", "j", "k"))
  
  # Keep schools < School size max
  schtable      <- table(mydata$sschlcde)
  school        <- as.numeric(names(schtable))
  sch.size      <- as.numeric(schtable)
  nsch          <- length(school)
  
  # dependent variable
  mislistmis   <- c(55555555, 99999999, 99959995) #error code
  
  # This function prepares the data and the network 
  gen.data  <- function(db) {
    G       <- vector("list", nsch) #the true network
    Gobtmis <- vector("list", nsch) #is 1 if the entry of G is observed and o otherwise
    nmatchm <- vector("list", nsch) #Number of unmatched male friends
    nmatchf <- vector("list", nsch) #Number of unmatched female friends
    matchm  <- vector("list", nsch) #Number of matched male friends
    matchf  <- vector("list", nsch) #Number of matched female friends
    friendm <- vector("list", nsch) #Number of declared male friends
    friendf <- vector("list", nsch) #Number of declared male friends
    
    for (i in 1:nsch) {
      cat("School: ", i, "/", nsch, "\n", sep = "")
      schi         <- school[i]
      dbi          <- db[db$sschlcde == schi,]
      Ni           <- nrow(dbi)
      Gi           <- matrix(0, Ni, Ni)
      Giom         <- matrix(1, Ni, Ni)
      Giot         <- matrix(1, Ni, Ni)
      Giotm        <- matrix(1, Ni, Ni)
      
      nmatchim     <- numeric() #will contain missing male links
      nmatchif     <- numeric() #will contain missing female links
      matchim      <- numeric() #will contain male links
      matchif      <- numeric() #will contain female links
      for (j in 1:Ni) {
        idxm       <- which(dbi$aid %in% unlist(dbi[j, mf_coln]))
        idxf       <- which(dbi$aid %in% unlist(dbi[j, ff_coln]))
        idx        <- c(idxm, idxf)
        matchim[j] <- length(idxm) # observed male links
        matchif[j] <- length(idxf) # observed female links
        Gi[j,idx]  <- 1
        
        # missing links
        idxm.miss  <- which(unlist(dbi[j, mf_coln]) %in% c(mislistmis, -3))
        idxf.miss  <- which(unlist(dbi[j, ff_coln]) %in% c(mislistmis, -3))
        nmatchim[j]<- length(idxm.miss) # Number of missing male links
        nmatchif[j]<- length(idxf.miss) # Number of missing female links
        
        
        # If there are messing links, then we doubt the values
        if(nmatchim[j] > 0) {
          # Giom[j, (Gi[j,] == 0)] <- 0
          Giom[j, dbi$female == 0] <- 0
        }
        if(nmatchif[j] > 0) {
          # Giom[j, (Gi[j,] == 0)] <- 0
          Giom[j, dbi$female == 1] <- 0
        }        
        
        # Top coding
        # If top coding for any sex, then we doubt all zeros in Gi[j,] for the same sex
        # male
        if ((length(idxm) + nmatchim[j]) == 5){# There are 5 male friends declared
          Giot[j, (Gi[j,] == 0) & (dbi$female == 0)] <- 0
        }
        # female
        if ((length(idxf) + nmatchif[j]) == 5){# There are 5 female friends declared
          Giot[j, (Gi[j,] == 0) & (dbi$female == 1)] <- 0
        }
        
        # Missing and top coding
        Giotm[j,]  <- Giom[j,]*Giot[j,]
      }
      
      #  store G
      diag(Gi)     <- 0
      diag(Giotm)  <- 0
      G[[i]]       <- Gi
      Gobtmis[[i]] <- Giotm
      
      # unmatched
      matchm[[i]]  <- matchim
      matchf[[i]]  <- matchif
      
      # unmatched
      nmatchm[[i]] <- nmatchim
      nmatchf[[i]] <- nmatchif
      
      # friends 
      friendm[[i]] <- matchim + nmatchim
      friendf[[i]] <- matchif + nmatchif
    }
    
    
    list(G       = G,
         Gobtmis = Gobtmis,
         matchm  = matchm,
         matchf  = matchf,
         nmatchm = nmatchm,
         nmatchf = nmatchf,
         friendm = friendm,
         friendf = friendf)
  }
  
  # Use the function to prepare the data
  tmp     <- gen.data(mydata)
  G       <- tmp$G
  Gobtmis <- tmp$Gobtmis
  matchm  <- tmp$matchm
  matchf  <- tmp$matchf
  nmatchm <- tmp$nmatchm
  nmatchf <- tmp$nmatchf
  friendm <- tmp$friendm
  friendf <- tmp$friendf
  
  # Add degree in mydata
  mydata  <- mydata %>% mutate(
    # number of observed friends
    matchmall  = unlist(matchm),
    matchfall  = unlist(matchf),
    matchall   = matchmall + matchfall, 
    # number of unmatched friends
    nmatchmall = unlist(nmatchm),
    nmatchfall = unlist(nmatchf),
    nmatchall  = nmatchmall + nmatchfall, 
    # declared friends
    friendmall = unlist(friendm),
    friendfall = unlist(friendf),
    friendall  = friendmall + friendfall)
  
  # dataset for the logistic model
  # variable to include 1 (indicator == 1 if same value)
  va.log1       <- c("female", "hispanic", "racewhite", "raceblack", "raceasian", "withbothpar",
                     "melhigh", "memhigh", "mjprof")
  # variable to include  2 (absolute value of the difference)
  va.log2       <- c("age", "yearinschl")
  # distance
  dist1         <- function(x, y) as.numeric(x == y)
  dist2         <- function(x, y) abs(x - y)
  dist3         <- function(x, y) x
  dist4         <- function(x, y) y
  
  # data
  tmp           <- c(0, cumsum(sch.size))
  # for va.log1
  X1tmp         <- do.call("cbind", lapply(va.log1, function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx       <- mydata[c(tmp[x] + 1):tmp[x+1], z, drop = TRUE]   
      matrix(kronecker(va.zx, va.zx, FUN = dist1), sch.size[x])}))
  }))
  # for va.log2
  X2tmp         <- do.call("cbind", lapply(va.log2, function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx     <- mydata[c(tmp[x] + 1):tmp[x+1], z, drop = TRUE]   
      matrix(kronecker(va.zx, va.zx, FUN = dist2), sch.size[x])}))
  }))
  # for va.log1 and 2 with dist 3
  X3tmp         <- do.call("cbind", lapply(c("ID", va.log1, va.log2), function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx     <- mydata[c(tmp[x] + 1):tmp[x+1], z, drop = TRUE]   
      matrix(kronecker(va.zx, va.zx, FUN = dist3), sch.size[x])}))
  }))
  # for va.log1 and 2 with dist 4
  X4tmp         <- do.call("cbind", lapply(c("ID", "sschlcde", va.log1, va.log2), function(z) {
    mat.to.vec(lapply(1:nsch, function(x) {
      va.zx     <- mydata[c(tmp[x] + 1):tmp[x+1], z, drop = TRUE]   
      matrix(kronecker(va.zx, va.zx, FUN = dist4), sch.size[x])}))
  }))
  Xlogit        <- cbind(X1tmp, X2tmp, X3tmp, X4tmp)  
  colnames(Xlogit)  <- c("same.sex", va.log1[-1], "diff.age", "diff.yearinschl",
                         paste0(c("ID", va.log1, va.log2), ".i"), 
                         paste0(c("ID", "sschlcde", va.log1, va.log2), ".j"))
  
  rm(list = c("X1tmp", "X2tmp", "X3tmp", "X4tmp", "tmp"))
  gc()
  save(list = ls(all = TRUE), file = filname)
}

############################ Descriptive stat ##################################
# Descriptive statistic function 
my.stat.des <- function(x) {
  out <- c(mean(x), sd(x), min(x), max(x))
  names(out) <- c("Mean", "St. Dev.", "Min", "Max")
  out
}

# all the variables 
allVar        <- mydata[,va.all.names]

# Descriptive stat
sdes           <- round(t(apply(allVar, 2, my.stat.des)), 3)
print(sdes)

write.csv(sdes, file = "des.stat.csv")

#Proportion of missing links
sum(mydata$nmatchall)/sum(mydata$friendall) 

# Proportion of missing to be inferred due to error code and top coding
sum(sapply(Gobtmis, function(x) sum(x == 0) - nrow(x)))/sum(sch.size*(sch.size - 1))

# graph missing links
ggplot(data = data.frame(mm = mydata$nmatchall), aes(x = mm)) +
  geom_bar(color = "black", fill = "#eeeeee") + 
  theme_bw() + xlab("") + ylab("Frequency") + 
  scale_x_discrete(limits = 0:10) 
# size 7 × 4 inch