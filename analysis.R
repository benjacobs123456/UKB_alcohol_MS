# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
setwd("/data/Wolfson-UKBB-Dobson/UKB_alcohol")

# read data
df = read_tsv("ukb_pheno_21_01_22.tsv")

# read in source of report data
source_of_report_data = read_tsv("../ukb_pheno_17_03_21/ukb_pheno_nocognitive_17032021.tsv.gz",col_types=cols_only(
  `Source of report of G35 (multiple sclerosis).0.0` = col_character(),
  `Date G35 first reported (multiple sclerosis).0.0` = col_character(),
  EID = col_double()
  ))

# define MS status
source_of_report_data = source_of_report_data %>% mutate(MS_status = ifelse(!is.na(`Source of report of G35 (multiple sclerosis).0.0`),1,0))

# combine
df = df %>% dplyr::select(-MS_status) %>% left_join(source_of_report_data,by="EID")

count_ms = function(){
  df %>%
    filter(!is.na(MS_status)) %>%
    dplyr::count(MS_status) %>%
    mutate(total = sum(n),percent = n/sum(n)*100)
}

count_ms()

# exclude participants who have withdrawn
withdrawn = read_tsv("../helper_progs_and_key/excluded_indivs",col_names=FALSE)
df = df %>% filter(!EID %in% withdrawn$X1)
count_ms()


# restrict to genetic cohort
# Remove non-European participants
df = df %>% filter(`Genetic ethnic grouping.0.0` =="Caucasian")
count_ms()

# rename cols
df$alcohol_freq = df$`Alcohol intake frequency..0.0.y`
df$alcohol_status = df$`Alcohol drinker status.0.0.y`

# exclude NAs
df = df %>% filter(!(`Alcohol drinker status.0.0.y` == "Prefer not to answer"))
df$never_drinker = df$alcohol_status=="Never"
df$drb_15_carrier = df$DRB_15>0

count_ms()

# define age at ms diagnosis 
df = df %>%
  mutate(`Date G35 first reported (multiple sclerosis).0.0` = ifelse(`Date G35 first reported (multiple sclerosis).0.0` == "1902-02-02",NA,`Date G35 first reported (multiple sclerosis).0.0`)) %>%
  mutate(ms_dx_date = as.Date(`Date G35 first reported (multiple sclerosis).0.0`,
                              format = "%Y-%m-%d")) %>%
  mutate(year_of_recruitment = `Age at recruitment.0.0` + `Year of birth.0.0`) %>%
  mutate(approx_date_recruitment = as.Date(paste0(year_of_recruitment,"-01-01"),
                                           format = "%Y-%m-%d") ) %>%
  mutate(delta_recruitment_ms_date = as.numeric(approx_date_recruitment - ms_dx_date)/365.25) %>%
  mutate(age_ms_dx = as.numeric(ms_dx_date - as.Date(paste0(`Year of birth.0.0`,"-01-01"),
                                        format = "%Y-%m-%d") ) / 365.25) 

# one person was diagnosed at 0.7 years - implausible. Remove. 
df = df %>% mutate(age_ms_dx = ifelse(age_ms_dx <1,NA,age_ms_dx) )
              

# descriptive stats
df %>% group_by(MS_status) %>%
  summarise_at(.vars = c("Age at recruitment.0.0","age_ms_dx","delta_recruitment_ms_date"),
                .funs = c("median","IQR"),na.rm=T)
df %>% group_by(MS_status) %>% dplyr::count(Sex.0.0) %>% mutate(prop = n/sum(n))


df %>% group_by(MS_status) %>% dplyr::count(alcohol_status) %>% mutate(prop = n/sum(n))

table(df$MS_status,df$alcohol_status)/rowSums(table(df$MS_status,df$alcohol_status))*100
table(df$MS_status,df$DRB_15)
table(df$MS_status,df$DRB_15)/rowSums(table(df$MS_status,df$DRB_15))

table(df$MS_status,df$never_drinker,df$drb_15_carrier)
table(df$MS_status,df$never_drinker,df$`Alcohol intake versus 10 years previously.0.0`)

# look at who drinks same or less now
table(df$MS_status,df$`Alcohol intake versus 10 years previously.0.0`)/rowSums(table(df$MS_status,df$`Alcohol intake versus 10 years previously.0.0`))*100

# plot
drb_labs = c("DRB1*15 Negative","DRB1*15 Positive")
names(drb_labs) = levels(factor(df$drb_15_carrier))
df$alcohol_status = factor(df$alcohol_status,levels=c(
"Current",
"Previous",
"Never"
))

p0=ggplot(df,aes(factor(MS_status),fill=alcohol_status))+facet_wrap(~drb_15_carrier,labeller = labeller(drb_15_carrier = drb_labs))+
  geom_bar(position="fill",col="black")+
scale_x_discrete(labels=c("Controls","Cases"))+
theme_bw()+
labs(x="MS status",y="Proportion",fill="Alcohol drinking status")+
  scale_fill_brewer(palette="Paired")


# initialise df for model results 
res_df = data.frame()

# models
print_ci = function(z,name){
  x=summary(z)$coef
  beta = x[nrow(x),1]
  se = x[nrow(x),2]
  p=se = x[nrow(x),4]
  lower_ci = exp(beta - 1.96*se)
  upper_ci = exp(beta + 1.96*se)
  or = exp(beta)
  df = data.frame(name,beta,se,p,lower_ci,upper_ci,or)
  res_df <<- bind_rows(res_df,df)
  paste0("OR=",round(or,2),"(95% CI ",round(lower_ci,2)," - ",round(upper_ci,2),"), p=",round(p,3))
}

# match cases to 10 controls for age and sex
matched_df = df
cases = matched_df %>% filter(MS_status==1)
controls = matched_df %>% filter(MS_status==0)
overall_control_df = data.frame()
for(i in c(1:nrow(cases))){
  message("Doing case ",i)
  this_case = cases[i,]
  these_controls = controls %>%
    filter(!EID %in% overall_control_df$EID) %>%
    filter(abs(`Age at recruitment.0.0` - this_case$`Age at recruitment.0.0`)<1) %>%
    filter(`Sex.0.0`== this_case$`Sex.0.0`) %>%
    sample_n(size=10,replace=F)
  overall_control_df <<- bind_rows(overall_control_df,these_controls)
}
matched_df = bind_rows(overall_control_df,cases)

table(matched_df$MS_status)
multi_model_matched = glm(data=matched_df,MS_status ~  `Age at recruitment.0.0` + Sex.0.0 + never_drinker,family=binomial(link="logit"))
print_ci(multi_model_matched,name="matched")

# adjusted models in full dataset
multi_model = glm(data=df,
                  MS_status ~  `Age at recruitment.0.0` + 
                    Sex.0.0 + 
                    never_drinker,
                  family=binomial(link="logit"))
print_ci(multi_model,name="whole_cohort_age_sex")


multi_model_townsend = glm(data=df,
                           MS_status ~  `Age at recruitment.0.0` + 
                             Sex.0.0 + 
                             `Townsend deprivation index at recruitment.0.0` + 
                             never_drinker,
                           family=binomial(link="logit"))
print_ci(multi_model_townsend,name = "whole_cohort_age_sex_townsend")

# associations with confounders
## townsend
ggplot(df,
       aes(never_drinker,
           `Townsend deprivation index at recruitment.0.0`))+
  geom_boxplot()

a = df[df$never_drinker==T,][['Townsend deprivation index at recruitment.0.0']]
b = df[df$never_drinker==F,][['Townsend deprivation index at recruitment.0.0']]
t.test(a,b)


# cbmi 
ggplot(df,
       aes(`Comparative body size at age 10.0.0`,
           fill=never_drinker))+
  geom_bar(position = "fill")

df = df %>% 
  mutate(cbmi = `Comparative body size at age 10.0.0`) %>%
  mutate(cbmi = ifelse(!cbmi %in% c("About average","Thinner","Plumper"),NA,cbmi))
table_chisq = table(df$never_drinker,df$cbmi)
chisq.test(table_chisq)

# smoking
df = df %>% 
  mutate(smok = `Smoking status.0.0.x`) %>%
  mutate(smok = ifelse(smok == "Prefer not to answer",NA,smok))
table_chisq = table(df$never_drinker,df$smok)
chisq.test(table_chisq)
df %>% group_by(never_drinker) %>%
  dplyr::count(smok) %>%
  filter(!is.na(smok)) %>%
  mutate(prop = n/sum(n)*100)

ggplot(df,
       aes(smok,
           fill=never_drinker))+
  geom_bar(position = "fill")


# further models 

multi_model_smok = glm(data=df,
                           MS_status ~  `Age at recruitment.0.0` + 
                             Sex.0.0 + 
                             smok + 
                             never_drinker,
                           family=binomial(link="logit"))
print_ci(multi_model_smok,name = "whole_cohort_age_sex_smoking")
df = df %>% mutate(smok_binary = ifelse(smok == "Never","never","ever"))
multi_model_smok = glm(data=df,
                       MS_status ~  `Age at recruitment.0.0` + 
                         Sex.0.0 + 
                         smok_binary + 
                         never_drinker,
                       family=binomial(link="logit"))
print_ci(multi_model_smok,name = "whole_cohort_age_sex_smoking_binary")

multi_model_smok_townsend = glm(data=df,
                       MS_status ~  `Age at recruitment.0.0` + 
                         Sex.0.0 + 
                         smok + `Townsend deprivation index at recruitment.0.0`+
                         never_drinker,
                       family=binomial(link="logit"))
print_ci(multi_model_smok_townsend,name = "whole_cohort_age_sex_smoking_townsend")

multi_model_smok_townsend = glm(data=df,
                                MS_status ~  `Age at recruitment.0.0` + 
                                  Sex.0.0 + 
                                  smok_binary + `Townsend deprivation index at recruitment.0.0`+
                                  never_drinker,
                                family=binomial(link="logit"))
print_ci(multi_model_smok_townsend,name = "whole_cohort_age_sex_smoking_binary_townsend")

# define age at smoking
df = df %>%
  mutate(age_smok = ifelse(
    `Smoking status.0.0.x` == "Current",`Age started smoking in current smokers.0.0.x`,NA
  ))  %>%
  mutate(age_smok = ifelse(
    `Smoking status.0.0.x` == "Previous",`Age started smoking in former smokers.0.0.x`,NA
  ))

df %>% filter(age_ms_dx<age_smok)

multi_model_smok = glm(data=df %>% filter(! (MS_status == 1 & age_ms_dx<age_smok)),
                       MS_status ~  `Age at recruitment.0.0` + 
                         Sex.0.0 + 
                         smok + 
                         never_drinker,
                       family=binomial(link="logit"))
summary(multi_model_smok)

multi_model_smok = glm(data=df %>% filter(! (MS_status == 1 & age_ms_dx<age_smok)),
                       MS_status ~  `Age at recruitment.0.0` + 
                         Sex.0.0 + 
                         smok_binary + 
                         never_drinker,
                       family=binomial(link="logit"))
summary(multi_model_smok)

df %>% dplyr::count(smok,never_drinker,MS_status)

# HLA
df %>%
  group_by(MS_status) %>%
  dplyr::count(drb_15_carrier) %>%
  mutate(prop = n/sum(n))
hla_multi_model = glm(data=df,MS_status ~ `Age at recruitment.0.0` + Sex.0.0 + `Genetic principal components.0.1` + `Genetic principal components.0.2` + `Genetic principal components.0.3` + `Genetic principal components.0.4` +drb_15_carrier,family=binomial(link="logit"))
summary(hla_multi_model)
print_ci(hla_multi_model,name="hla")


hla_alcohol_interaction_multi_model = glm(data=df,MS_status ~ `Age at recruitment.0.0` + Sex.0.0 + `Genetic principal components.0.1` + `Genetic principal components.0.2` + `Genetic principal components.0.3` + `Genetic principal components.0.4` + never_drinker * drb_15_carrier,family=binomial(link="logit"))
summary(hla_alcohol_interaction_multi_model)
print_ci(hla_alcohol_interaction_multi_model,name="hla_interaction")

# fewer covars
hla_alcohol_interaction_multi_model = glm(data=df,MS_status ~ never_drinker * drb_15_carrier,family=binomial(link="logit"))
summary(hla_alcohol_interaction_multi_model)
print_ci(hla_alcohol_interaction_multi_model,name="hla_interaction_no_covar")

# quant
hla_alcohol_interaction_multi_model = glm(data=df,MS_status ~ never_drinker * DRB_15,family=binomial(link="logit"))
summary(hla_alcohol_interaction_multi_model)
print_ci(hla_alcohol_interaction_multi_model,name="hla_interaction_additive")



# additive interaction

library(doParallel)
registerDoParallel(cores=20)
trials = 1000
r = foreach(icount(trials), .combine=cbind) %dopar% {
  df = df %>% dplyr::select(MS_status,`Age at recruitment.0.0`,`Sex.0.0`,never_drinker,drb_15_carrier,contains("enetic"))
  df = sample_n(df, size=nrow(df), replace=TRUE)
  model = glm(data=df,
              MS_status~`Age at recruitment.0.0`+
                `Sex.0.0`+
                `Genetic principal components.0.1` + `Genetic principal components.0.2` + `Genetic principal components.0.3` + `Genetic principal components.0.4` +
                never_drinker*drb_15_carrier,
              family=binomial(link="logit"))
  coeffs=as.numeric(coef(model))
  RERI=exp(coeffs[8]+coeffs[9]+coeffs[10])-exp(coeffs[8])-exp(coeffs[9])+1
  AP=RERI/(exp(coeffs[8]+coeffs[9]+coeffs[10]))
}


# actual estimate
df = df %>% dplyr::select(MS_status,`Age at recruitment.0.0`,`Sex.0.0`,never_drinker,drb_15_carrier,contains("enetic"))
model = glm(data=df,
            MS_status~`Age at recruitment.0.0`+
              `Sex.0.0`+
              `Genetic principal components.0.1` + `Genetic principal components.0.2` + `Genetic principal components.0.3` + `Genetic principal components.0.4` +
              never_drinker*drb_15_carrier,
            family=binomial(link="logit"))
coeffs=as.numeric(coef(model))
RERI=exp(coeffs[8]+coeffs[9]+coeffs[10])-exp(coeffs[8])-exp(coeffs[9])+1
AP=RERI/(exp(coeffs[8]+coeffs[9]+coeffs[10]))

median_ap = median(r)
cis = quantile(r[1,],c(0.025,0.975))
one_tailed_p = 1-((sum(r>0)+1)/(length(r)+1))
res = data.frame(AP,median_ap,cis[1],cis[2],one_tailed_p,row.names = NULL)
write_csv(res,"alcohol_results.csv")

# stratified models 
hla_neg = df %>% filter(drb_15_carrier==F)
hla_pos = df %>% filter(drb_15_carrier==T)

# adjusted models in full dataset
hla_strat_neg_model = glm(data=hla_neg,
                  MS_status ~  `Age at recruitment.0.0` + 
                    Sex.0.0 + 
                    never_drinker,
                  family=binomial(link="logit"))
print_ci(hla_strat_neg_model,name="hla_negative_model")

hla_strat_pos_model = glm(data=hla_pos,
                          MS_status ~  `Age at recruitment.0.0` + 
                            Sex.0.0 + 
                            never_drinker,
                          family=binomial(link="logit"))
print_ci(hla_strat_pos_model,name="hla_pos_model")

# plot
res_df = res_df[-c(8:11),]
res_df$name = 
  c("Matched cohort (10:1)",
    "Whole cohort (age + sex)",
    "Whole cohort (age + sex + SES)",
    "Whole cohort (age + sex + smoking)",
    "Whole cohort (age + sex + smoking [binary])",
    "Whole cohort (age + sex + smoking + SES)",
    "Whole cohort (age + sex + smoking [binary] + SES)",
    "HLA-DRB1*15:01 negative",
    "HLA-DRB1*15:01 positive")
res_df = res_df %>% filter(!grepl("binary",name))
p=ggplot(res_df,
       aes(or,name))+
  geom_point(shape=15,size=3)+
  geom_errorbarh(mapping = aes(xmin = lower_ci,xmax=upper_ci,y=name),height=0.1)+
  scale_x_log10()+
  theme_minimal()+
  labs(x="Odds Ratio for MS",y="Model")+
  geom_vline(xintercept = 1,alpha=0.1)


png("fig1.png",res=300,units="in",width=8,height=3)
grid.arrange(p,p0)
dev.off()

# power calcs
n = nrow(df)
n_ms = nrow(df %>% filter(MS_status==1))
ms_prevalence = n_ms/n
n_non_drink = nrow(df %>% filter(never_drinker==T))
exposure_prevalence = (n_non_drink)/n
rr = 1/0.75

pvals=list()
for(i in c(1:1000)){
# simulate phenos
message("Doing ",i)
df = data.frame(exposure = rbinom(n=n,size=1,prob=exposure_prevalence))
non_exposed = df %>% filter(exposure==0)
exposed = df %>% filter(exposure==1)
non_exposed$case_status = rbinom(n=nrow(non_exposed),size=1,prob=ms_prevalence)
exposed$case_status = rbinom(n=nrow(exposed),size=1,prob=ms_prevalence*rr)

df = bind_rows(exposed,non_exposed)
model = glm(data = df,case_status ~ exposure,family=binomial(link="logit"))
pval = summary(model)$coefficients[2,4]
pvals[i] = pval
}

power = table(unlist(pvals)<0.05)[2]/sum(table(unlist(pvals)<0.05))
message("Power:")
print(power)
