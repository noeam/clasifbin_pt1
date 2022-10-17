#############################################################
#
#   Program: datasim.R
#  
#   Project: Simulate data for reproducible research
#
#   Biostatisticians: Meridith Blevins, MS
#                     Bryan E Shepherd, PhD
#
#   Purpose: This program simulates CCASAnet data for use 
#            in reproducible research.  By simulating data
#            from the actual data, we provide reviwers and
#            readers a valid platform to rerun and test our
#            analysis scripts, without sharing the actual
#            datasets which are not openly available for
#            one reason or another as discussed in our IWHOD
#            abstract.
#
#   Notes: 
#          
#   Created: 16 November 2015
#	  
#   Revisions: 
#   
#
#############################################################

rm(list=ls()) # SHOULD BE FIRST LINE IN ALL R PROGRAMS - CLEARS NAMESPACE

# if(.Platform$OS.type == "windows") {
#   baseDir <- file.path("B:", "Projects/") } else {
#     baseDir <- file.path("~", "Projects/") }
# setwd(baseDir)
# setwd("CCASAnet/datasim"); 

warning("IN RStudio, go to Session -> Set Working Directory -> To Source File Location")
setwd(gsub("/code","",getwd()))
set.seed(201612)

## load rms for rcs() function
library(rms)
## load nnet for multinom() function
library(nnet)

## READ IN DATA FILES
readtables <- paste0(c("basic","art","follow","lab_cd4","lab_rna","visit"),"_sim")

## previously simulated data is on dropbox in /output
## CHOOSE FIRST SELECTS THE TEXT STRING OCCURING BEFORE THE SPECIFIED SEPARATER
choosefirst <- function(var,sep=".") unlist(lapply(strsplit(var,sep,fixed=TRUE),function(x) x[1]))
## DETERMINE WHICH TABLES EXIST IN '/output'
existingtables <- choosefirst(list.files("output"))
if(!all(readtables %in% existingtables)) stop(paste0("Please place the expected tables in /output folder named as specified: (",paste(readtables,collapse=", "),")"))
## READ IN ALL EXISTING TABLES
for(i in 1:length(readtables)){
  if(!is.na(readtables[i])){
    readcsv <- read.csv(paste("output/",readtables[i],".csv",sep=""),header=TRUE,stringsAsFactors = FALSE,na.strings=c(NA,""))
    names(readcsv) <- tolower(names(readcsv))
    assign(gsub("_sim","",readtables[i]),readcsv)
  }
}

## convert and then remove baseline dates that were impercise to more than just the day
basic <- cleanup.import(basic,datevars=c("birth_d","baseline_d","hivdiagnosis_d","aids_d","aids_cl_d"),dateformat="%Y-%m-%d")
basic$baseline_d[!(is.na(basic$baseline_d_a) | basic$baseline_d_a=="D")] <- NA
basic$age <- with(basic,as.numeric(baseline_d - birth_d))/365.25 ## in YEARS
basic$aids.miss <- ifelse(is.na(basic$aids_y) | basic$aids_y==9,1,0)

## subset data to adults
adultpatients <- basic$patient[basic$age >= 18]
basic   <- basic[basic$patient %in% adultpatients,]
art     <- art[art$patient %in% adultpatients,]
follow  <- follow[follow$patient %in% adultpatients,]
lab_cd4 <- lab_cd4[lab_cd4$patient %in% adultpatients,]
lab_rna <- lab_rna[lab_rna$patient %in% adultpatients,]
visit   <- visit[visit$patient %in% adultpatients,]

## to call the function basicsim()
sites <- unique(basic$site)

## for testing the function
country <- sites[1]

basicsim <- function(country){
  j <- which(sites==country)
  
  ## simulate data by site, so rather than contstantly subset, we will do this inside of a function
  basicsub   <- basic[basic$site==country,]
  artsub     <- art[art$site==country,]
  followsub  <- follow[follow$site==country,]
  lab_cd4sub <- lab_cd4[lab_cd4$site==country,]
  lab_rnasub <- lab_rna[lab_rna$site==country,]
  visitsub   <- visit[visit$site==country,]
  
  ## we start simulation by selecting the baseline date using sampling with replacement 
  ## had tested the distribution, but decided it wasn't really necessary...
  # repeat{
    datesamp <- with(basicsub,sample(baseline_d[!is.na(baseline_d)],nrow(basicsub),replace=TRUE))
#     ## Kolmogorov-Smirnov Tests to ensure our sample matches the distribution of the original data (as it should!)
#     ksp <- ks.test(as.numeric(datesamp),as.numeric(basicsub$baseline_d))$p.value
#     if(ksp <= 0.10) warning(paste("Sampling distribution was too different for",country))
#     if(ksp > 0.10) break
#     }
    basic_sim <- data.frame(patient=paste(substr(country,1,2),1:length(datesamp),sep="."),site=country,baseline_d=datesamp)

  ## need numeric dates for modeling and prediction
  basicsub$baseline_d_num <- as.numeric(basicsub$baseline_d) - median(as.numeric(basicsub$baseline_d),na.rm=TRUE)
  basic_sim$baseline_d_num <- as.numeric(basic_sim$baseline_d) - median(as.numeric(basicsub$baseline_d),na.rm=TRUE)
  
  ## next we'll generate sex from the distribution of sex across baseline_d
  m1 <- glm(male ~ rcs(baseline_d_num,5),data=basicsub,family="binomial")
  basic_sim$male <- mapply(rbinom,n=1,size=1,prob=predict(m1,basic_sim,type="response"))

  ## next we'll generate age and then birth_d from the distribution of age across sex and baseline_d
  m2 <- lm(age ~ male * rcs(baseline_d_num,5),data=basicsub)
  basic_sim$age <- mapply(rnorm,n=1,mean=predict(m2,basic_sim,type="response"),sd=summary(m2)$sigma)
  basic_sim$age[basic_sim$age<18] <- 18 ## impute as minimum 18 years old
  par(mfrow=c(2,1)); hist(basicsub$age); hist(basic_sim$age)
  basic_sim$birth_d <- basic_sim$baseline_d - round(basic_sim$age*365.25)   ###  Bryan fixed (because dates thought age was measured in days)
  
  ## next we'll generate hiv diagnosis date
  basicsub$hivdiagnosis_d_num <- as.numeric(basicsub$hivdiagnosis_d) - median(as.numeric(basicsub$hivdiagnosis_d),na.rm=TRUE)
  m3 <- lm(hivdiagnosis_d_num ~ rcs(age,4) * male + rcs(baseline_d_num,5),data=basicsub)
  basic_sim$hivdiagnosis_d_num <- mapply(rnorm,n=1,mean=predict(m3,basic_sim,type="response"),sd=summary(m3)$sigma)
  basic_sim$hivdiagnosis_d <- as.Date(basic_sim$hivdiagnosis_d_num + median(as.numeric(basicsub$hivdiagnosis_d),na.rm=TRUE),origin=c("1970-01-01"))
  par(mfrow=c(2,1)); hist(basicsub$hivdiagnosis_d_num); hist(basic_sim$hivdiagnosis_d_num)
  ## in case some dates are generated following date of birth, replace those with date of birth (safeguard as this doesn't appear to actually happen)
  basic_sim$hivdiagnosis_d[basic_sim$hivdiagnosis_d<basic_sim$birth_d] <- basic_sim$birth_d[basic_sim$hivdiagnosis_d<basic_sim$birth_d]

  ## aids_y indicator (skipping aids_d because only brazil and honduras complete this variable)
  if(!all(basicsub$aids.miss==1)){
    basicsub$aids_y[basicsub$aids_y==9] <- NA
    m4 <- glm(aids_y ~ rcs(age,4) * male + rcs(baseline_d_num,5) + rcs(hivdiagnosis_d_num,5),data=basicsub,family="binomial")
    basic_sim$aids_y <- mapply(rbinom,n=1,size=1,prob=predict(m4,basic_sim,type="response"))
    mean(basic_sim$aids_y); mean(basicsub$aids_y,na.rm=TRUE)
    m4.miss <- glm(aids.miss ~ rcs(age,4) * male + rcs(baseline_d_num,5) + rcs(hivdiagnosis_d_num,5),data=basicsub,family="binomial")
    basic_sim$aids.miss <- mapply(rbinom,n=1,size=1,prob=predict(m4,basic_sim,type="response"))
    basic_sim$aids_y[basic_sim$aids.miss] <- 9
  }
  if(all(basicsub$aids.miss==1)){
    basic_sim$aids_y <- 9
    basic_sim$aids.miss <- 1
  }
  ## mode of transmission
  if(country!="haiti"){
    basicsub$mode <- as.character(basicsub$mode)
    m5 <- multinom(mode ~ rcs(age,4) * male + rcs(baseline_d_num,5) + rcs(hivdiagnosis_d_num,5),data=basicsub)
    getmode <- function(prob){
      x <- rmultinom(n=1,size=1,prob=prob)
      return(row.names(x)[x==1])
    }
    basic_sim$mode <- apply(predict(m5,basic_sim,type="probs"),1,getmode)
    warning("these distributions didn't look great for some reason -- can bryan help resolve this?")
    table(basicsub$mode)/nrow(basicsub)
    table(basic_sim$mode)/nrow(basic_sim)
  }
  if(country=="haiti"){
    basic_sim$mode <- "Unknown"
  }

  with(basic_sim,data.frame(patient,birth_d,site,male,mode,hivdiagnosis_d,aids_y))  
  ## recart_y WILL BE GENERATED FROM art data
  ## center, mode_oth, all *_d_a will be ignored
  
  return(basic_sim)
}

basic_sim <- do.call( rbind, lapply(sites,basicsim))
table(basic_sim$site)
table(basic$site)
if(anyDuplicated(basic_sim$patient)>0) stop("There are duplicates in the patient ID.")

####  New things added on 2016-03-08
p.clin.trial<-table(basic$clinicaltrial_y,basic$site)[2,]/table(basic$site)
basic_sim$clinicaltrial_y<-c(rbinom(sum(basic_sim$site=="argentina"),1,p.clin.trial["argentina"]),
                           rbinom(sum(basic_sim$site=="brazil"),1,p.clin.trial["brazil"]),
                           rbinom(sum(basic_sim$site=="chile"),1,p.clin.trial["chile"]),
                           rbinom(sum(basic_sim$site=="haiti"),1,p.clin.trial["haiti"]),
                           rbinom(sum(basic_sim$site=="honduras"),1,p.clin.trial["honduras"]),
                           rbinom(sum(basic_sim$site=="mexico"),1,p.clin.trial["mexico"]),
                           rbinom(sum(basic_sim$site=="peru"),1,p.clin.trial["peru"]))
table(basic_sim$site,basic_sim$clinicaltrial_y)

basic_sim$mode_oth<-NA

u1<-runif(sum(basic_sim$site=="argentina"))
rec1<-ifelse(u1<table(basic$recart_y,basic$site)[1,"argentina"]/sum(basic$site=="argentina"),0,ifelse(u1>1-table(basic$recart_y,basic$site)[3,"argentina"]/sum(basic$site=="argentina"),9,1))
u2<-runif(sum(basic_sim$site=="brazil"))
rec2<-ifelse(u2<table(basic$recart_y,basic$site)[1,"brazil"]/sum(basic$site=="brazil"),0,ifelse(u2>1-table(basic$recart_y,basic$site)[3,"brazil"]/sum(basic$site=="brazil"),9,1))
u3<-runif(sum(basic_sim$site=="chile"))
rec3<-ifelse(u3<table(basic$recart_y,basic$site)[1,"chile"]/sum(basic$site=="chile"),0,ifelse(u3>1-table(basic$recart_y,basic$site)[3,"chile"]/sum(basic$site=="chile"),9,1))
u4<-runif(sum(basic_sim$site=="haiti"))
rec4<-ifelse(u4<table(basic$recart_y,basic$site)[1,"haiti"]/sum(basic$site=="haiti"),0,ifelse(u4>1-table(basic$recart_y,basic$site)[3,"haiti"]/sum(basic$site=="haiti"),9,1))
u5<-runif(sum(basic_sim$site=="honduras"))
rec5<-ifelse(u5<table(basic$recart_y,basic$site)[1,"honduras"]/sum(basic$site=="honduras"),0,ifelse(u5>1-table(basic$recart_y,basic$site)[3,"honduras"]/sum(basic$site=="honduras"),9,1))
u6<-runif(sum(basic_sim$site=="mexico"))
rec6<-ifelse(u6<table(basic$recart_y,basic$site)[1,"mexico"]/sum(basic$site=="mexico"),0,ifelse(u6>1-table(basic$recart_y,basic$site)[3,"mexico"]/sum(basic$site=="mexico"),9,1))
u7<-runif(sum(basic_sim$site=="peru"))
rec7<-ifelse(u7<table(basic$recart_y,basic$site)[1,"peru"]/sum(basic$site=="peru"),0,ifelse(u7>1-table(basic$recart_y,basic$site)[3,"peru"]/sum(basic$site=="peru"),9,1))
basic_sim$recart_y<-c(rec1,rec2,rec3,rec4,rec5,rec6,rec7)

u<-runif(length(basic_sim$site))
basic_sim$aids_y<-with(basic_sim, ifelse(site=="haiti"&u<.25,1,
                                  ifelse(site=="haiti"&u>.25,0,
                                  ifelse(site=="mexico"&u<.5,1,
                                  ifelse(site=="mexico"&u<.9,0,
                                  ifelse(site=="mexico"&u>.9,9,
                                  ifelse(site=="peru"&u<.4,1,
                                  ifelse(site=="peru"&u<.9,0,
                                  ifelse(site=="peru"&u>.9,9,aids_y)))))))))



###### NOW THAT PATIENT DEMOGRAPHICS and OTHER tblBASIC VARIABLES WERE SIMULATED STRATIFIED BY SITE, 
###### WE GO AHEAD AND POOL ALL THE DATA FOR CLINICAL VARIABLES AND INCLUDE SITE AS A COVARIATE

######  FIRST CD4
# basic$age<-floor(basic$age/365.25)
# basic_sim$age<-floor(basic_sim$age/365.25)
basic.sim<-basic_sim
#basic.sim$baseline_d_num<-basic.sim$hivdiagnosis_d_num<-NULL

## create a data frame of first CD4 among patients
lab_cd4$cd4_d<-as.Date(lab_cd4$cd4_d)
lab_cd4<-lab_cd4[order(lab_cd4$patient,lab_cd4$cd4_d),]
first.cd4<-lab_cd4[!duplicated(lab_cd4$patient),]
## merge with basic so we can calculate difference between baseline CD4 and baseline date
cd4.basic<-merge(basic,first.cd4,by=c("site","patient"),all.x=TRUE)
cd4.basic$diff.time<-as.numeric(cd4.basic$cd4_d-cd4.basic$baseline_d)
cd4.basic$baseline_d_n<-with(cd4.basic,as.numeric(baseline_d))
## create indicator for missing aids_y
# cd4.basic$aids.miss<-with(cd4.basic, as.character(ifelse(is.na(aids_y)|aids_y==9,"missing",
#                                                          ifelse(aids_y==1,"AIDS","not AIDS"))))
cd4.basic$diff.time<-with(cd4.basic, ifelse(diff.time< -500, -500,
                                     ifelse(diff.time>  500,  500, diff.time)))
mod.diff.time<-lm(diff.time~rcs(baseline_d_n,3)+aids.miss+site+mode+male+rcs(age,3), data=cd4.basic)

cd4.basic.sim<-basic.sim
# cd4.basic.sim$aids.miss<-with(cd4.basic.sim, as.character(ifelse(is.na(aids_y)|aids_y==9,"missing",
#                                                              ifelse(aids_y==1,"AIDS","not AIDS"))))
cd4.basic.sim$baseline_d_n<-with(cd4.basic.sim,as.numeric(baseline_d))
## simultation of time from baseline to first CD4
cd4.basic.sim$diff.time<-predict(mod.diff.time, newdata=cd4.basic.sim)+sample(mod.diff.time$residuals,length(cd4.basic.sim$patient), replace=TRUE)

cd4.basic$sqrt.cd4<-with(cd4.basic, sqrt(cd4_v))
mod.cd4<-lm(sqrt.cd4~rcs(baseline_d_n,3)+aids.miss+site+mode+male+rcs(age,3)+rcs(diff.time,3), data=cd4.basic)
## simulation of first CD4 count
cd4.basic.sim$sqrt.cd4<-predict(mod.cd4, newdata=cd4.basic.sim)+sample(mod.cd4$residuals,length(cd4.basic.sim$patient), replace=TRUE)
cd4.basic.sim$cd4_v.all<-cd4.basic.sim$sqrt.cd4^2

## impute missing values for first CD4 count
cd4.basic$cd4.miss<-with(cd4.basic, ifelse(is.na(cd4_v),1,0))
mod.cd4.miss<-glm(cd4.miss ~ rcs(baseline_d_n,3)+aids.miss+site+mode+male+rcs(age,3), data=cd4.basic, family="binomial")
p.cd4.miss<-predict(mod.cd4.miss, newdata=cd4.basic.sim, type="response")
cd4.basic.sim$cd4.miss<-with(cd4.basic.sim, rbinom(length(p.cd4.miss),1,p.cd4.miss))
cd4.basic.sim$cd4_v<-with(cd4.basic.sim, ifelse(cd4.miss==1,NA,cd4_v.all))
cd4.basic.sim$cd4_d<-with(cd4.basic.sim, as.Date(baseline_d+diff.time))
cd4.basic.sim$cd4_d<-with(cd4.basic.sim, cd4_d+ifelse(cd4.miss==1,NA,0))

######  ART INITIATION
art$art_sd<-as.Date(art$art_sd)
art<-art[order(art$site,art$patient,art$art_sd),]
## create data frame of first instance of ART drugs
art1<-art[!duplicated(art$patient),]

d<-merge(cd4.basic,art1, by=c("site","patient"),all.x=TRUE)
d$no.art<-with(d, ifelse(is.na(art_id),1,0))

d.s<-cd4.basic.sim
d.s$cd4.yes<-1-d.s$cd4.miss
## Fill in missing CD4 counts with zero, will also include dummy variable for missing CD4
## this is important because in many settings, patients are not prescribed ART until CD4 drops below a given threshold
d.s$cd4.v.yes<-with(d.s, ifelse(is.na(cd4_v),0,cd4_v))
d.s$cd4.v.yes<-round(d.s$cd4.v.yes)

d$cd4.yes<-1-d$cd4.miss
d$cd4.v.yes<-with(d, ifelse(is.na(cd4_v),0,cd4_v))
mod.no.art<-glm(no.art~rcs(baseline_d_n,3)+aids.miss+site+mode+male+rcs(age,3)+cd4.miss+rcs(sqrt(cd4.v.yes),3),  family="binomial", data=d)
p.no.art<-predict(mod.no.art, newdata=d.s, type="response")
d.s$no.art<-rbinom(length(p.no.art), 1, p.no.art)

##  Next I need to simulate the date of ART initiation, given that I have started ART.  
## It is possible to start ART before the baseline date.
d$diff.time<-with(d, as.numeric(art_sd-baseline_d))
d$diff.time<-with(d, ifelse(diff.time< -300, -300, diff.time))
mod.art.time<-lm(diff.time~rcs(baseline_d_n,3)+aids.miss+site+mode+male+rcs(age,3)+cd4.miss+rcs(sqrt(cd4.v.yes),3), data=d)
d.s$diff.time<-predict(mod.art.time, newdata=d.s)+sample(mod.art.time$residuals,length(d.s$patient), replace=TRUE)
d.s$art_sd <- d.s$baseline_d + d.s$diff.time

#########  Commented out on 2016-03-08
# #####  I can then simulate the specific regimen.
# drop.sparse <- with(as.data.frame.table(table(d$art_id)),Var1[Freq/nrow(d)<0.02])
# table(d$art_id %in% drop.sparse)
# d$art_id[d$art_id %in% drop.sparse] <- NA
# reg1 <- multinom(art_id ~ rcs(age,4) * male + rcs(baseline_d_n,5)+mode+cd4.miss+rcs(sqrt(cd4.v.yes),3)+site+aids.miss,data=d)
# getmode <- function(prob){
#   x <- rmultinom(n=1,size=1,prob=prob)
#   return(row.names(x)[x==1])
# }
# d.s$art_id <- apply(predict(reg1,d,type="probs"),1,getmode)
# table(d$art_id)/nrow(d)
# table(d.s$art_id)/nrow(d.s)
# ## distributions are okay, not identical, but not way off either.

#########  Added on 2016-03-08
art.unique<-art[duplicated(art$patient)==FALSE,]

art.unique1<-art.unique[art.unique$site=="argentina",]
samp1<-sample(1:length(art.unique1$site),sum(d.s$site=="argentina"),replace=TRUE)
art1<-art.unique1[samp1,]
art1a<-art1[,c("art_id","pi","nnrti1","nnrti2","nnrti","nrti","t20","ccr5","ii1","ii2","rtv_drug","numdrugs","art_class")]
art.unique2<-art.unique[art.unique$site=="brazil",]
samp2<-sample(1:length(art.unique2$site),sum(d.s$site=="brazil"),replace=TRUE)
art2<-art.unique2[samp2,]
art2a<-art2[,c("art_id","pi","nnrti1","nnrti2","nnrti","nrti","t20","ccr5","ii1","ii2","rtv_drug","numdrugs","art_class")]
art.unique3<-art.unique[art.unique$site=="chile",]
samp3<-sample(1:length(art.unique3$site),sum(d.s$site=="chile"),replace=TRUE)
art3<-art.unique3[samp3,]
art3a<-art3[,c("art_id","pi","nnrti1","nnrti2","nnrti","nrti","t20","ccr5","ii1","ii2","rtv_drug","numdrugs","art_class")]
art.unique4<-art.unique[art.unique$site=="haiti",]
samp4<-sample(1:length(art.unique4$site),sum(d.s$site=="haiti"),replace=TRUE)
art4<-art.unique4[samp4,]
art4a<-art4[,c("art_id","pi","nnrti1","nnrti2","nnrti","nrti","t20","ccr5","ii1","ii2","rtv_drug","numdrugs","art_class")]
art.unique5<-art.unique[art.unique$site=="honduras",]
samp5<-sample(1:length(art.unique5$site),sum(d.s$site=="honduras"),replace=TRUE)
art5<-art.unique5[samp5,]
art5a<-art5[,c("art_id","pi","nnrti1","nnrti2","nnrti","nrti","t20","ccr5","ii1","ii2","rtv_drug","numdrugs","art_class")]
art.unique6<-art.unique[art.unique$site=="mexico",]
samp6<-sample(1:length(art.unique6$site),sum(d.s$site=="mexico"),replace=TRUE)
art6<-art.unique6[samp6,]
art6a<-art6[,c("art_id","pi","nnrti1","nnrti2","nnrti","nrti","t20","ccr5","ii1","ii2","rtv_drug","numdrugs","art_class")]
art.unique7<-art.unique[art.unique$site=="peru",]
samp7<-sample(1:length(art.unique7$site),sum(d.s$site=="peru"),replace=TRUE)
art7<-art.unique7[samp7,]
art7a<-art7[,c("art_id","pi","nnrti1","nnrti2","nnrti","nrti","t20","ccr5","ii1","ii2","rtv_drug","numdrugs","art_class")]
art.tot<-rbind(art1a,art2a,art3a,art4a,art5a,art6a,art7a)

d.s$art_id<-art.tot$art_id
d.s$pi<-art.tot$pi
d.s$nnrti1<-art.tot$nnrti1
d.s$nnrti2<-art.tot$nnrti2
d.s$nnrti<-art.tot$nnrti
d.s$nrti<-art.tot$nrti
d.s$t20<-art.tot$t20
d.s$ccr5<-art.tot$ccr5
d.s$ii1<-art.tot$ii1
d.s$ii2<-art.tot$ii2
d.s$rtv_drug<-art.tot$rtv_drug
d.s$numdrugs<-art.tot$numdrugs
d.s$art_class<-art.tot$art_class


## FOR NOW WE ARE HAPPY JUST HAVING ART INITIATION INFO
art_sim <- with(d.s[!is.na(d.s$art_sd),],data.frame(patient,site,art_id,art_sd,pi,nnrti1,nnrti2,nnrti,nrti,t20,ccr5,ii1,ii2,rtv_drug,numdrugs,art_class))
art_sim$art_ed<-NA
art_sim$art_rs<-""

#####  Then I'll simulate CD4 trajectories 


d2<-merge(d,lab_cd4,by=c("site","patient"),all=TRUE)
d2<-d2[,c("site","patient","baseline_d","male","age","aids_y","mode","cd4_v.x","cd4_d.x","cd4_v.y","cd4_d.y","art_sd")]
# d2<-d2[!is.na(d2$cd4_d.y),]
d2$on.ART<-with(d2, ifelse(is.na(art_sd), 0,
                    ifelse(art_sd<cd4_d.y, 1, 0)))
d2$t<-with(d2, as.numeric(cd4_d.y-cd4_d.x))
d2$t.ART<-with(d2, as.numeric(cd4_d.y-art_sd))
d2$t.ART1<-with(d2, ifelse(is.na(t.ART),0,t.ART))
d2$log.t.ART<-with(d2, ifelse(t.ART1<1,0,log(t.ART1)))
#d2$log.t.ART<-with(d2, ifelse(is.na(t.ART),0,
#                              ifelse(t.ART<1,0,log(t.ART))))
d2$sqrt.cd4<-sqrt(d2$cd4_v.y)
d2$sqrt.cd4<-ifelse(d2$sqrt.cd4>50,50,d2$sqrt.cd4)
d2$sqrt.cd4.init<-sqrt(d2$cd4_v.x)


library("lme4")

d2a<-d2[d2$t>0,]
d2a$log.t<-log(d2a$t)

mod<-lmer(sqrt.cd4 ~ rcs(log.t,5) + rcs(sqrt.cd4.init,4) + male + rcs(age,4)  + on.ART + rcs(log.t.ART,5) + (1|patient), data=d2a)
summary(mod)
stuff<-predict(mod)

###  Next I need to predict from this model

d.s1<-d.s[order(d.s$cd4_d),]

id.samp<-sample(unique(d2$patient),sum(d.s1$cd4.yes),replace=FALSE)

d3<-d2[which(d2$patient %in% id.samp),]
dim(d2)
dim(d3)
d3<-d3[order(d3$patient),]
d3.unique<-d3[!duplicated(d3$patient),]
d3.unique<-d3.unique[order(d3.unique$cd4_d.x),]

d.s2<-d.s1[d.s1$cd4.yes==1,]
d.s2$match.id<-d3.unique$patient

d4<-d3[,c("patient","t")]

d.s3<-merge(d.s2,d4,by.x="match.id",by.y="patient")
d.s3<-d.s3[order(d.s3$patient,d.s3$t),]
d.s3$sqrt.cd4.init<-sqrt(d.s3$cd4.v.yes)
d.s3$on.ART<-with(d.s3, ifelse(cd4_d+t>art_sd,1,0))
d.s3$t.ART<-with(d.s3, ifelse(on.ART==1, cd4_d+t-art_sd, 0))
d.s3$log.t.ART<-with(d.s3, ifelse(t.ART<1,0,log(t.ART)))
d.s3$log.t<-log(d.s3$t)

d.s4<-d.s3[d.s3$t>0,c("patient","male","log.t","age","t","on.ART","sqrt.cd4.init","t.ART","log.t.ART")]
d.s4$sqrt.cd4 <-predict(mod, newdata=d.s4, re.form=NA) + rnorm(n=nrow(d.s4),mean=0,sd=attr(summary(mod)$varcor$patient,"stddev"))  #  random effects set at 0
## merge with CD4 baseline data to get full lab_cd4 dataset!
lab_cd4_sim <- merge(with(cd4.basic.sim,data.frame(patient,site,cd4_d)),
                     rbind(with(d.s4,data.frame(patient,cd4_v=sqrt.cd4^2,time=t)),
                           with(cd4.basic.sim,data.frame(patient,cd4_v,time=0))),all.y=TRUE)
lab_cd4_sim$cd4_d <- with(lab_cd4_sim,cd4_d + time)
## for now, we are not determining cd4_per since this is an adult only simulation

lab_cd4_sim<-lab_cd4_sim[order(lab_cd4_sim$patient,lab_cd4_sim$cd4_d),]
lab_cd4_sim$cd4_v<-round(lab_cd4_sim$cd4_v)
summary(lab_cd4_sim$cd4_v)    ##  There are a lot more missing values than in the original data

#####  and VL trajectories.  

#####  Then mortality 

#####  My strategy for simulating mortality is to use table basic and to use the last recorded CD4 value

### Creating a dataset that only has the last CD4 per patient
d.cd4.j<-lab_cd4[order(lab_cd4$patient,lab_cd4$cd4_d),]
d.cd4.j$stuff<-duplicated(d.cd4.j$patient)
d.cd4.j$stuff1<-c(d.cd4.j$stuff[-1],FALSE)
d.cd4.j<-d.cd4.j[d.cd4.j$stuff1==FALSE,]
d.cd4.j<-d.cd4.j[,c("patient","site","cd4_d","cd4_v")]      

### Merging the last CD4 per patient with basic and mortality data
ds1<-merge(basic,d.cd4.j,by=c('patient','site'),all.x=TRUE)
ds2<-merge(ds1,follow,by=c("patient","site"),all.x=TRUE)

###  Fitting model for whether a patient dies
ds2$birth_d<-as.Date(ds2$birth_d)
ds2$baseline_d<-as.Date(ds2$baseline_d)
ds2$cd4_d<-as.Date(ds2$cd4_d)
ds2$cd4.miss<-with(ds2, ifelse(is.na(cd4_v),1,0))
ds2$sqrt.cd4<-with(ds2, ifelse(cd4.miss==1,0,sqrt(cd4_v)))
ds2$aids.cat<-with(ds2, as.character(ifelse(is.na(aids_cl_y), "missing",
                                     ifelse(aids_cl_y==9,"missing",
                                     ifelse(aids_cl_y==0,"not AIDS",
                                     ifelse(aids_cl_y==1,"AIDS","Bogey"))))))
ds2$cd4.d.n<-with(ds2, ifelse(cd4.miss==1,0,as.numeric(cd4_d)))
ds2$birth.d.n<-with(ds2, ifelse(birth_d=="1900-01-01",NA,as.numeric(birth_d)))
ds2$baseline.d.n<-as.numeric(ds2$baseline_d)
fit.death<-glm(death_y~site+rcs(birth.d.n,4)+rcs(baseline.d.n,3)+mode+male+cd4.miss+rcs(sqrt.cd4,4) + rcs(cd4.d.n,4) + aids.cat, family=binomial, data=ds2)

###  Given a patient dies, date of death based on date of last CD4 measurement (or baseline date if don't have a CD4 measurement)
ds2$death_d<-as.Date(ds2$death_d)
ds2$diff.death.cd4<-as.numeric(ds2$death_d-ds2$cd4_d)
ds2$diff.death.cd4.log<-log(ds2$diff.death.cd4+1)
fit.time1<-lm(diff.death.cd4.log~site + rcs(cd4.d.n,4) + rcs(sqrt.cd4,4), subset=cd4.miss==0&death_y==1, data=ds2)  ## model if have a CD4 measurement
ds2$diff.death.base<-as.numeric(ds2$death_d-ds2$baseline_d)
ds2$diff.death.base.log<-log(ds2$diff.death.base+1)

fit.time2<-lm(diff.death.base.log ~ site + rcs(baseline.d.n,6), subset=cd4.miss==1&death_y==1, data=ds2)  ## model if don't have a CD4 measurement

############## SIMULATE DEATH DATA ##############
############## SIMULATE DEATH DATA ##############
############## SIMULATE DEATH DATA ##############
############## SIMULATE DEATH DATA ##############
############## SIMULATE DEATH DATA ##############
ds2$diff.lalive.base<-as.numeric(as.Date(ds2$l_alive_d)-ds2$baseline_d)
# ds2$diff.lalive.cd4<-as.numeric(as.Date(ds2$l_alive_d)-ds2$cd4_d)

cox.death<-cph(Surv(diff.lalive.base,death_y)~site+rcs(birth.d.n,4)+rcs(baseline.d.n,3)+mode+male+cd4.miss+sqrt.cd4 + rcs(cd4.d.n,4) + aids.cat, data=ds2)
## removing splines so that the design matrix is much easier -- could come back and redo this part
cox.death1<-cph(Surv(diff.lalive.base,death_y)~site+birth.d.n+baseline.d.n+mode+male+cd4.miss+sqrt.cd4 + cd4.d.n + aids.cat, data=ds2,x=TRUE,y=TRUE)

###  Now predicting with the simulated data

### Creating a dataset that only has the last CD4 per patient
d.cd4.k<-lab_cd4_sim[order(lab_cd4_sim$patient,lab_cd4_sim$cd4_d),]
d.cd4.k$stuff<-duplicated(d.cd4.k$patient)
d.cd4.k$stuff1<-c(d.cd4.k$stuff[-1],FALSE)
d.cd4.k<-d.cd4.k[d.cd4.k$stuff1==FALSE,]
d.cd4.k<-d.cd4.k[,c("patient","site","cd4_d","cd4_v")]      

ds2s<-merge(basic_sim,d.cd4.k,by=c('patient','site'),all.x=TRUE)

###  Predicting whether a patient dies
ds2s$birth_d<-as.Date(ds2s$birth_d)
ds2s$baseline_d<-as.Date(ds2s$baseline_d)
ds2s$cd4_d<-as.Date(ds2s$cd4_d)
ds2s$cd4.miss<-with(ds2s, ifelse(is.na(cd4_v),1,0))
ds2s$sqrt.cd4<-with(ds2s, ifelse(cd4.miss==1,0,sqrt(cd4_v)))
ds2s$aids.cat<-with(ds2s, as.character(ifelse(is.na(aids_y), "missing",
                                            ifelse(aids_y==9,"missing",
                                                   ifelse(aids_y==0,"not AIDS",
                                                          ifelse(aids_y==1,"AIDS","Bogey"))))))
ds2s$cd4.d.n<-with(ds2s, ifelse(cd4.miss==1,0,as.numeric(cd4_d)))
ds2s$birth.d.n<-with(ds2s, ifelse(birth_d=="1900-01-01",NA,as.numeric(birth_d)))
ds2s$baseline.d.n<-as.numeric(ds2s$baseline_d)

lp<-predict(fit.death, newdata=ds2s)
prob.death<-exp(lp)/(1+exp(lp))
ds2s$death_y<-rbinom(length(lp),1,prob.death)

###  Predicting date of death

## predicting follow-up time
x1 <- with(ds2s,data.frame(site=="brazil",site=="chile",site=="haiti",site=="honduras",site=="mexico",site=="peru",
                           birth.d.n,baseline.d.n,
                           mode=="Generic Sexual",mode=="Hemophiliac",mode=="Heterosexual contact",
                           mode=="Heterosexual contact and Injecting drug user",mode=="Homo/Bisexual and Injecting drug user",
                           mode=="Homosexual contact",mode=="Injecting drug user",mode=="Other (specify in mode_oth)",
                           mode=="Perinatal",mode=="Transfusion nonhemophilia related",mode=="Unknown",
                           male,cd4.miss,sqrt.cd4,cd4.d.n,aids.cat=="missing",aids.cat=="not AIDS"))
x <- apply(x1,2,as.numeric)
# the three lines beneath produced too many times near
# myrates <- exp(x%*%coefficients(cox.death1)+1)  # the risk exp(beta*x), parameters for exp r.v.
# ds2s$diff.lalive.base <- mapply(rexp,n=1,rate=myrates) # generates the r.v.
# ds2s$diff.lalive.base[ds2s$diff.lalive.base>max(ds2$diff.lalive.base)] <- sample(ds2$diff.lalive.base,size=sum(ds2s$diff.lalive.base>max(ds2$diff.lalive.base)),replace=TRUE)

# Weibull latent event times
lambda <- 0.1
rho <- 1
v <- runif(n=nrow(ds2))
ds2s$diff.lalive.base <- (- log(v) / (lambda * exp(x%*%coefficients(cox.death1))))^(1 / rho)
sum(ds2s$diff.lalive.base[ds2s$death_y==1]>max(ds2$diff.lalive.base[ds2$death_y==1]),na.rm=TRUE)
ds2s$diff.lalive.base[ds2s$diff.lalive.base>max(ds2$diff.lalive.base)] <- sample(ds2$diff.lalive.base,size=sum(ds2s$diff.lalive.base>max(ds2$diff.lalive.base)),replace=TRUE)
# plot(ds2s$diff.lalive.base,ds2$diff.lalive.base)
# hist(ds2s$diff.lalive.base)
# hist(ds2$diff.lalive.base)

ds2s$death_d <- ds2s$baseline_d + ds2s$diff.lalive.base
ds2s$death_d[ds2s$death_y!=1] <- NA
# sum(ds2s$death_d>as.Date("2016-12-31"),na.rm=TRUE)

# ds3s<-ds2s[ds2s$cd4.miss==0&ds2s$death_y==1,]
# death.log.time<-predict(fit.time1, newdata=ds3s)+sample(fit.time1$residuals,length(ds3s$patient),replace=TRUE)
# death.time<-round(exp(death.log.time)-1)
# death.time<-ifelse(death.time<0,0,death.time)
# ds3s$death_d<-ds3s$baseline_d+death.time
# 
# ds4s<-ds2s[ds2s$cd4.miss==1&ds2s$death_y==1,]
# death.log.time2<-predict(fit.time2, newdata=ds4s)+sample(fit.time2$residuals,length(ds4s$patient),replace=TRUE)
# death.time2<-round(exp(death.log.time2)-1)
# death.time2<-ifelse(death.time2<0,0,death.time2)
# ds4s$death_d<-ds4s$baseline_d+death.time2
# 
# ds5s<-ds2s[ds2s$death_y==0,]
# ds5s$death_d<-NA
# 
# ds6s<-rbind(ds3s,ds4s,ds5s)
# 
safe.ifelse <- function(cond, yes, no){
  class.y <- class(yes)
  X <- ifelse(cond,yes,no)
  class(X) <-class.y; return(X)
}
## writing over about 30 that were past the database close
ds2s$death_d<-with(ds2s, safe.ifelse(death_d>as.Date("2016-12-31"),as.Date("2016-12-31"),death_d))
ds7s<-merge(ds2s,art_sim[,c("patient","site","art_sd","art_ed")],by=c("patient","site"),all.x=TRUE)

######  l_alive_d set as date of death if died, date of last CD4 if have a non-missing CD4, and otherwise baseline_d.  
######  I think this can be improved to incorporate VL and ART data.  Also, some cd4_d are very late -- we should fix.
ds7s$l_alive_d<-with(ds7s, safe.ifelse(!is.na(death_d),death_d,
                           safe.ifelse(!is.na(cd4_d)&!is.na(art_sd),pmax(cd4_d,art_sd,baseline_d),
                           safe.ifelse(!is.na(cd4_d),pmax(cd4_d,baseline_d),
                           safe.ifelse(!is.na(art_sd),pmax(art_sd,baseline_d),
                                       baseline_d)))))

follow_sim<-ds7s[,c("patient","site","l_alive_d","death_y","death_d")]
follow_sim<-follow_sim[order(follow_sim$patient),]

######  Now getting rid of all data after the date of death (l_alive_d), and fixing up dates.
art_sim<-merge(art_sim,follow_sim[,c("patient","site","l_alive_d")],by=c("patient","site"),all.x=TRUE)
art_sim<-art_sim[art_sim$l_alive_d>=art_sim$art_sd,]
art_sim$l_alive_d<-NULL
ord<-with(art_sim, order(site,patient,art_sd))
art_sim<-art_sim[ord,]

lab_cd4_sim<-merge(lab_cd4_sim,follow_sim[,c("patient","site","l_alive_d")],by=c("patient","site"),all.x=TRUE)
lab_cd4_sim<-lab_cd4_sim[lab_cd4_sim$l_alive_d>=lab_cd4_sim$cd4_d,]
lab_cd4_sim$l_alive_d<-NULL
ord<-with(lab_cd4_sim, order(site,patient,cd4_d))
lab_cd4_sim<-lab_cd4_sim[ord,]

basic_sim$hivdiagnosis_d<-with(basic_sim, safe.ifelse(hivdiagnosis_d>baseline_d,baseline_d,hivdiagnosis_d))

junk<-merge(basic_sim,follow_sim,by=c("patient","site"),all.x=TRUE)


#####  next we'll simulate "lab_rna"

rna<-lab_cd4_sim
rna$rna_d<-rna$cd4_d
rna$cd4_d<-NULL
rna1<-merge(rna,art_sim[,c("patient","site","art_sd")],by=c("patient","site"),all.x=TRUE)
rna1$t.from.art<-round(as.numeric(with(rna1, rna_d-art_sd)))
rna2<-rna1[-which(rna1$site=="haiti"),]

lab_rna1<-merge(lab_rna,art.unique,by=c("patient","site"))
lab_rna1$detect<-ifelse(lab_rna1$rna_v<=0,0,1)
lab_rna1$t.from.art<-round(as.numeric(with(lab_rna1, as.Date(rna_d)-as.Date(art_sd))))

rna2$v1<-rbinom(length(rna2$patient),1,mean(lab_rna1$detect[lab_rna1$t.from.art<=0],na.rm=TRUE))
rna2$v2<-rbinom(length(rna2$patient),1,mean(lab_rna1$detect[lab_rna1$t.from.art>0],na.rm=TRUE))
rna2$v3<-sample(lab_rna1$rna_v[lab_rna1$t.from.art<=0&lab_rna1$detect==1],length(rna2$v1),replace=TRUE)
rna2$v4<-sample(lab_rna1$rna_v[lab_rna1$t.from.art>0&lab_rna1$detect==1],length(rna2$v1),replace=TRUE)
rna2$v5<-sample(lab_rna1$rna_v[lab_rna1$detect==0],length(rna2$v1),replace=TRUE)

rna2$rna_v<-with(rna2, ifelse(t.from.art<=0 &v1==1, v3,
                       ifelse(t.from.art<=0 &v1==0, v5,
                       ifelse(t.from.art>0 & v2==1, v4,
                       ifelse(t.from.art>0 & v2==0, v5, -9999)))))

lab_rna_sim<-rna2[,c("patient","site","rna_d","rna_v")]

#####  then "visit" data for last 

visit_sim<- basic_sim[,c("patient","site","baseline_d")]
visit_sim$visit_d<-visit_sim$baseline_d
visit_sim$baseline_d<-NULL
visit_sim$whostage<-visit_sim$cdcstage<-""
  
#####  save regimen changes (ie, more records for "art") for last
#####  ever need "center" data?  not for dataviz.

basic_sim$aids_d<-basic_sim$baseline_d
basic_sim$aids_d<-with(basic_sim, as.Date(ifelse(aids_y==1, as.character(aids_d), NA)))
basic_sim$birth_d_a<-"D"
basic_sim$aids_cl_y<-basic_sim$aids_y
basic_sim$aids_cl_d<-basic_sim$aids_d

follow_sim$death_d_a<-"D"
follow_sim$drop_rs<-" "
u<-rbinom(length(follow_sim$drop_rs),1,0.002)
follow_sim$drop_rs <- safe.ifelse(u==1,"Patient transfered",follow_sim$drop_rs)
follow_sim$drop_rs_oth <- follow_sim$drop_rs


## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## create a dataset with ART data and date of ART start (called haart_d)
art_change <- merge(art,with(art,aggregate(list(haart_d=art_sd),by=list(patient=patient),min)))
art_change$art_start2sd <- as.numeric(with(art_change,art_sd - haart_d))
art_change$art_start2ed <- as.numeric(with(art_change,as.Date(art_ed) - haart_d))
art_change0 <- art_change[art_change$art_start2sd==0,]
names(art_change0)[names(art_change0)=="patient"] <- "actual_patient"
art_sim$art_start2sd <- 0


## assign random patient trajectory based on same baseline regimen
x <- with(art_sim,split(patient,art_id))
temp <- function(i){
  if(length(art_change0$actual_patient[art_change0$art_id %in% names(x)[[i]]])>0){
    sample(art_change0$actual_patient[art_change0$art_id %in% names(x)[[i]]],size=length(x[[i]]),replace=TRUE)
  }
  else{
    return(rep(NA,length(x[[i]])))
  }
}
x1 <- lapply(1:length(x),temp)
art_sim$actual_patient <- with(art_sim,unlist(x1,art_id))
## jitter time between regimens
myjitter5 <- sample(0:5,nrow(art_change),replace=TRUE)
art_change$art_start2sd <- art_change$art_start2sd + myjitter5
art_change$art_start2ed <- art_change$art_start2ed + myjitter5
names(art_change)[names(art_change)=="patient"] <- "actual_patient"
art_change$haart_d <- art_change$art_rs_oth <- art_change$art_ed_a <- art_change$art_sd_a <- art_change$art_ed <- art_change$art_sd <- NULL
# art_sim$art_start2sd <- NULL
intersect(names(art_change),names(art_sim))
setdiff(names(art_change),names(art_sim))
setdiff(names(art_sim),names(art_change))
art_change_nozero <- art_change[art_change$art_start2sd!=0,]

removesparse <- function(x) names(x[x<10])
## not sampling patient trajectories for baseline regimens that are very sparse (<10)
art_sim$actual_patient[art_sim$art_id %in% removesparse(table(art_sim$art_id))] <- NA
art_sim$haart_d <- art_sim$art_sd
art_change_nozero1 <- merge(art_sim[c("patient","actual_patient","haart_d")],art_change_nozero)
## this merge essentially sets them on top of eachother
art_sim1 <- merge(art_change_nozero1,art_sim,all=TRUE)
## need to add start and end times to haart_d
art_sim1$art_sd[is.na(art_sim1$art_sd)] <- with(art_sim1[is.na(art_sim1$art_sd),],haart_d + art_start2sd)
art_sim1$art_ed <- with(art_sim1,haart_d + art_start2ed)
## need to impute first art_ed as art_sd of first switch (if any)
art_sim1$art_ed2 <- with(art_sim1,unsplit(lapply(split(art_sd,patient),function(x) ifelse(length(x)==1,NA,as.character(sort(x)[2]))),patient))
art_sim1$art_ed[is.na(art_sim1$art_ed) & art_sim1$art_start2sd==0] <- as.Date(art_sim1$art_ed2[is.na(art_sim1$art_ed) & art_sim1$art_start2sd==0])
## table(is.na(art_sim1$art_ed)) <- these are single regimens (no change) so end date should be blank

art_sim2 <- merge(art_sim1,follow_sim[c("patient","l_alive_d")])
art_sim3 <- art_sim2[!(art_sim2$l_alive_d < art_sim2$art_sd),]
art_sim3$art_ed[art_sim3$l_alive_d < art_sim3$art_ed & !is.na(art_sim3$art_ed)] <- art_sim3$l_alive_d[art_sim3$l_alive_d < art_sim3$art_ed & !is.na(art_sim3$art_ed)]
## THIS IS THE NEW DATASET
## NUMBER OF SWITCHES: table(table(art_sim3$patient))
## NUMBER OF PATIENTS: length(unique(art_sim3$patient))
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)
## ART TRAJECTORY (SWITCHES IN REGIMENS)

write.csv(basic_sim,"output/basic_sim2.csv",na="",row.names=FALSE)
write.csv(art_sim3[c('patient','site','art_id','art_sd','pi',
                     'nnrti1','nnrti2','nnrti','nrti','t20','ccr5','ii1','ii2',
                     'rtv_drug','numdrugs','art_class','art_ed','art_rs')],"output/art_sim2.csv",na="",row.names=FALSE)
write.csv(follow_sim,"output/follow_sim2.csv",na="",row.names=FALSE)
write.csv(lab_cd4_sim,"output/lab_cd4_sim2.csv",na="",row.names=FALSE)
write.csv(lab_rna_sim,"output/lab_rna_sim2.csv",na="",row.names=FALSE)
write.csv(visit_sim,"output/visit_sim2.csv",na="",row.names=FALSE)


