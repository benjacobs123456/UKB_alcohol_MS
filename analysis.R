# load libraries
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(TwoSampleMR)
setwd("/data/Wolfson-UKBB-Dobson/UKB_alcohol")

# read data
df = read_tsv("ukb_pheno_21_01_22.tsv")

# read in source of report data
source_of_report_data = read_tsv("../ukb_pheno_17_03_21/ukb_pheno_nocognitive_17032021.tsv.gz",col_types=cols_only(
  `Source of report of G35 (multiple sclerosis).0.0` = col_character(),
  EID = col_double()
  ))

# define MS status
source_of_report_data = source_of_report_data %>% mutate(MS_status = ifelse(!is.na(`Source of report of G35 (multiple sclerosis).0.0`),1,0))


# combine
df = df %>% select(-MS_status) %>% left_join(source_of_report_data,by="EID")

# exclude participants who have withdrawn
withdrawn = read_tsv("../helper_progs_and_key/excluded_indivs",col_names=FALSE)
df = df %>% filter(!EID %in% withdrawn$X1)

# restrict to genetic cohort
# Remove non-European participants
df = df %>% filter(`Genetic ethnic grouping.0.0` =="Caucasian")
table(df$MS_status)
# filter relatedness
# nb this is taking care to exclude the non-ms control from each pair to boost case numbers
kin = read_table2("/data/Wolfson-UKBB-Dobson/helper_progs_and_key/ukb43101_rel_s488282.dat")
highly_related = kin %>% filter(Kinship>0.0884) %>% filter(ID1 %in% df$EID) %>% filter(ID2 %in% df$EID)
highly_related = highly_related %>% left_join((df %>% filter(EID %in% highly_related$ID1) %>% select(EID,MS_status)  %>% rename(ID1 = EID,MS_status_ID1 = MS_status)),by="ID1") %>%
left_join((df %>% filter(EID %in% highly_related$ID2) %>% select(EID,MS_status) %>% rename(ID2 = EID, MS_status_ID2 = MS_status)),by="ID2")

exclusion = bind_rows(highly_related %>% filter(MS_status_ID1==1 & MS_status_ID2==0) %>% select(ID2) %>% rename(EID = ID2),
highly_related %>% filter(MS_status_ID1==0 & MS_status_ID2==1) %>% select(ID1) %>% rename(EID = ID1),
highly_related %>% filter(MS_status_ID1==0 & MS_status_ID2==0) %>% select(ID1) %>% rename(EID = ID1))
df = df %>% filter(!EID %in% exclusion$EID)
table(df$MS_status)


# rename cols
df$alcohol_freq = df$`Alcohol intake frequency..0.0.y`
df$alcohol_status = df$`Alcohol drinker status.0.0.y`

# exclude NAs
df = df %>% filter(!(`Alcohol drinker status.0.0.y` == "Prefer not to answer"))
df$never_drinker = df$alcohol_status=="Never"
df$drb_15_carrier = df$DRB_15>0

# descriptive
table(df$MS_status)
table(df$MS_status) %>% sum
table(df$MS_status) / sum(table(df$MS_status))*100

table(df$alcohol_status)
table(df$DRB_15)

table(df$MS_status,df$alcohol_status)
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


png("prop_plot.png",res=300,units="in",width=6,height=4)
ggplot(df,aes(factor(MS_status),fill=alcohol_status))+facet_wrap(~drb_15_carrier,labeller = labeller(drb_15_carrier = drb_labs))+geom_bar(position="fill")+
scale_x_discrete(labels=c("Controls","Cases"))+
theme_bw()+
labs(x="MS status",y="Proportion",fill="Alcohol drinking status")
dev.off()

# models
print_ci = function(z){
  x=summary(z)$coef
  beta = x[nrow(x),1]
  se = x[nrow(x),2]
  p=se = x[nrow(x),4]
  lower_ci = exp(beta - 1.96*se)
  upper_ci = exp(beta + 1.96*se)
  or = exp(beta)
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
uni_model_matched = glm(data=matched_df,MS_status ~ never_drinker,family=binomial(link="logit"))
print_ci(uni_model_matched)
multi_model_matched = glm(data=matched_df,MS_status ~  `Age at recruitment.0.0` + Sex.0.0 + never_drinker,family=binomial(link="logit"))
print_ci(multi_model_matched)

# adjusted models in full dataset
multi_model = glm(data=df,MS_status ~  `Age at recruitment.0.0` + Sex.0.0 + never_drinker,family=binomial(link="logit"))
print_ci(multi_model)

multi_model_townsend = glm(data=df,MS_status ~  `Age at recruitment.0.0` + Sex.0.0 + `Townsend deprivation index at recruitment.0.0` + never_drinker,family=binomial(link="logit"))
print_ci(multi_model_townsend)

hla_multi_model = glm(data=df,MS_status ~ `Age at recruitment.0.0` + Sex.0.0 + `Genetic principal components.0.1` + `Genetic principal components.0.2` + `Genetic principal components.0.3` + `Genetic principal components.0.4` +drb_15_carrier,family=binomial(link="logit"))
summary(hla_multi_model)

hla_alcohol_model = glm(data=df,MS_status ~ never_drinker + drb_15_carrier,family=binomial(link="logit"))
summary(hla_alcohol_model)
hla_alcohol_interaction_model = glm(data=df,MS_status ~ never_drinker * drb_15_carrier,family=binomial(link="logit"))
summary(hla_alcohol_interaction_model)
hla_alcohol_interaction_multi_model = glm(data=df,MS_status ~ `Age at recruitment.0.0` + Sex.0.0 + `Genetic principal components.0.1` + `Genetic principal components.0.2` + `Genetic principal components.0.3` + `Genetic principal components.0.4` + never_drinker * drb_15_carrier,family=binomial(link="logit"))
summary(hla_alcohol_interaction_multi_model)



# additive interaction

library(doParallel)
registerDoParallel(cores=20)
trials = 10000
r = foreach(icount(trials), .combine=cbind) %dopar% {
  df = df %>% select(MS_status,`Age at recruitment.0.0`,`Sex.0.0`,never_drinker,drb_15_carrier,contains("enetic"))
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
df = df %>% select(MS_status,`Age at recruitment.0.0`,`Sex.0.0`,never_drinker,drb_15_carrier,contains("enetic"))
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


# power calcs
n = 378353
n_ms = 2100
ms_prevalence = n_ms/n
exposure_prevalence = (80+11641)/n
rr = 1.25

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
