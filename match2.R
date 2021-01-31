library(foreign)
library(tidyverse)
setwd("~/Documents/Year1/causalinf/siblingcash")
casheffects <- read.dta("data/publicdata.dta")
varbaseline = c("s_utilities", "s_durables", "s_infraest_hh", "s_age_sorteo", 
                "s_age_sorteo2", "s_years_back", "s_sexo", "s_single", "s_edadhead", 
                "s_yrshead", "s_tpersona", "s_num18", "s_puntaje", "s_ingtotal", "suba", "s_over_age",
                "factor(grade)", "factor(s_teneviv)", "factor(s_estcivil)", "factor(s_estrato)")


#Subsetting data as in paper "Improving the Design of Conditional Transfer Programs" Section V Sibling Effects
sibs <- casheffects %>%
  filter(fu_observed==1) %>%
  filter(grade!=11) %>% 
  group_by(fu_nim_hogar) %>%
  mutate(num_rsib=n()) %>%
  filter(num_rsib==2) %>%
  arrange(fu_nim_hogar) %>%
  mutate(tsib = if_else(row_number()==1, dplyr::lead(treatment, n = 1), dplyr::lag(treatment, n = 1))) %>%
  mutate(num_tsib=sum(treatment)) %>%
  filter(suba == 0 | grade >= 9) %>%
  ungroup()

#########Effects in untreated############################
sibs_ut <- sibs %>%
  filter(control==T)

#Regression from paper
untreatedreg_at <- lm(formula(paste("at_msamean ~ tsib + suba + factor(school_code) + ", 
                                    paste(varbaseline, collapse=" + "))), 
                      data = sibs_ut)
summary(untreatedreg_at)  

#Fitting propensity scores
propfit <- glm(formula(paste("tsib ~ suba + factor(school_code) + ", 
                              paste(varbaseline, collapse=" + "))),
                family=binomial(), data=sibs_ut)

sibs_ut <- sibs_ut %>%
  mutate(propscores=propfit$fitted.values)
boxplot(propscores~tsib, data=sibs_ut)


###First match
library(DOS2)
library(optmatch)
covariates<- c("s_utilities", "s_durables", "s_infraest_hh", "s_age_sorteo", 
               "s_age_sorteo2", "s_years_back", "s_sexo", "s_single", "s_edadhead", 
               "s_yrshead", "s_tpersona", "s_num18", "s_puntaje", "s_ingtotal", "suba", "s_over_age") #Just the continuous variables from varbaseline


dmat <- smahal(sibs_ut$tsib, select(sibs_ut, covariates))
dmat.cal <- DOS2::addcaliper(dmat, sibs_ut$tsib, sibs_ut$propscores, caliper=.2)
cashmatch = pairmatch(t(dmat.cal), data=sibs_ut)
sibs_ut_mpop <- sibs_ut %>% 
  mutate(cashmatch=cashmatch) %>%
  filter(!is.na(cashmatch))

sibs_ut_pairs <- sibs_ut_mpop %>%
  pivot_wider(id_cols=c(cashmatch), names_from = tsib, values_from=c(at_msamean, m_enrolled)) %>%
  mutate(attendancedif=at_msamean_1-at_msamean_0)

wilcox.test(sibs_ut_pairs$attendancedif, conf.int = T, conf.level = .90)

control_enroll=filter(sibs_ut_pairs, !is.na(m_enrolled_0), !is.na(m_enrolled_1))$m_enrolled_0
treatment_enroll=filter(sibs_ut_pairs, !is.na(m_enrolled_0), !is.na(m_enrolled_1))$m_enrolled_1
mcnemar.test(control_enroll, treatment_enroll)


###Exact match on gender
library(DiPs)

#Females
sibs_ut_female <- sibs_ut %>%
  filter(s_sexo==0)

dmat.female <- maha_sparse(1-sibs_ut_female$tsib, 
                       select(sibs_ut_female, covariates))
dmat.female <- addcaliper(dmat.female, 1-sibs_ut_female$tsib, 1-sibs_ut_female$propscores, c(-.5,.5), stdev=T)
cashmatch.female <- match(1-sibs_ut_female$tsib, dmat.female, sibs_ut_female)
sibs_ut_femalepairs <- cashmatch.female$data %>%
  pivot_wider(id_cols=c(mset), names_from = tsib, values_from=c(at_msamean, m_enrolled)) %>%
  mutate(attendancedif=at_msamean_1-at_msamean_0)
  
wilcox.test(sibs_ut_femalepairs$attendancedif, conf.int = T, conf.level = .90)

control_enroll_f=filter(sibs_ut_femalepairs, !is.na(m_enrolled_0), !is.na(m_enrolled_1))$m_enrolled_0
treatment_enroll_f=filter(sibs_ut_femalepairs, !is.na(m_enrolled_0), !is.na(m_enrolled_1))$m_enrolled_1
mcnemar.test(control_enroll_f, treatment_enroll_f)

#Males
sibs_ut_male <- sibs_ut %>%
  filter(s_sexo==1)

dmat.male <- maha_sparse(1-sibs_ut_male$tsib, 
                           select(sibs_ut_male, covariates))
dmat.male <- addcaliper(dmat.male, 1-sibs_ut_male$tsib, 1-sibs_ut_male$propscores, c(-.5,.5), stdev=T)
cashmatch.male <- match(1-sibs_ut_male$tsib, dmat.male, sibs_ut_male)

sibs_ut_malepairs <- cashmatch.male$data %>%
  pivot_wider(id_cols=c(mset), names_from = tsib, values_from=c(at_msamean, m_enrolled)) %>%
  mutate(attendancedif=at_msamean_1-at_msamean_0)

wilcox.test(sibs_ut_malepairs$attendancedif, conf.int = T, conf.level = .90)

###Checking covariate balance for our matches

#Function for computing standard differences
library(lazyeval)
standdifs <- function(covariate, prematchdata, postmatchdata) {
  presummary <- prematchdata %>% 
    group_by(tsib) %>%
    summarise_(mn=interp(~mean(var), var=as.name(covariate)), sd=interp(~sd(var), var=as.name(covariate)))
  premeancontrol=presummary$mn[1]
  premeantreat=presummary$mn[2]
  sdcontrol=presummary$sd[1]
  sdtreat=presummary$sd[2]
  postsummary <- postmatchdata %>% 
    group_by(tsib) %>%
    summarise_(mn=interp(~mean(var), var=as.name(covariate)))
  postmeancontrol=postsummary$mn[1]
  postmeantreat=postsummary$mn[2]  
  predif=abs(premeantreat-premeancontrol)/sqrt((sdtreat^2+sdcontrol^2)/2)
  postdif=abs(postmeantreat-postmeancontrol)/sqrt((sdtreat^2+sdcontrol^2)/2)
  return(c(predif, postdif))
} 

#Checking propensity score balance
paste("Standard difference of propensity scores before match:", standdifs("propscores", sibs_ut, sibs_ut)[2])
paste("Standard difference of propensity scores for first match:", standdifs("propscores", sibs_ut, sibs_ut_mpop)[2])
sibs_ut_mpop2 <- rbind(cashmatch.female$data, cashmatch.male$data)
paste("Standard difference of propensity scores for second match:", standdifs("propscores", sibs_ut, sibs_ut_mpop2)[2])


#Other covariates
unmatcheddifs <- c()
for (covariate in covariates) {
  difs=standdifs(covariate, sibs_ut, sibs_ut)
  unmatcheddifs <- c(unmatcheddifs, difs[1])
}

matcheddifs <- c()
for (covariate in covariates) {
  difs <- standdifs(covariate, sibs_ut, sibs_ut_mpop)
  matcheddifs <- c(matcheddifs, difs[2])
}

matcheddifs2 <- c()
for (covariate in covariates) {
  difs <- standdifs(covariate, sibs_ut, sibs_ut_mpop2)
  matcheddifs2 <- c(matcheddifs2, difs[2])
}

#Boxplots of standard differences
standards<-list(unmatcheddifs, matcheddifs, matcheddifs2)
names(standards) <- c("Pre match", "Post first match", "Post second match")
par(mgp=c(3,2,0), cex.axis=.75)
boxplot(standards, ylab="standard differences")







