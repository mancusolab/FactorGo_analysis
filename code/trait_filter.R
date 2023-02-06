### check traits in PanUKBB
library(tidyverse)
library(data.table)

## phenocode = FieldID
lookup <- read_tsv("./data/ukbb/PanUKBB/PanUKB_phe_04112022.tsv")

lookup %>% filter(!is.na(n_cases_EUR)) %>% 
  filter(!is.na(filename))

summary(lookup$n_cases_EUR) # 63 NAs (not available to EUR pop)

####### clean Sample size and trait_type ######
## Five trait_types: "biomarkers","categorical","continuous", "icd10","phecode","prescriptions"
## continuous traits are "biomarkers" and "continuous"
## Prescription data was extracted from UKB data field 42039
## "biomarkers","categorical","continuous" has coding field from UKB

# total: 7228
lookup <- lookup %>% 
  mutate(n_controls_EUR=as.numeric(n_controls_EUR), 
         n_cases_EUR=as.numeric(n_cases_EUR)) %>% 
  mutate(BIN_QT = ifelse(trait_type %in% c("biomarkers", "continuous"), "QT", "BIN"),
         N = ifelse(BIN_QT == "QT", n_cases_EUR, n_cases_EUR+n_controls_EUR))

#lookup %>% #distinct(phenotype_qc_EUR)
table(lookup$trait_type)
summary(lookup$n_cases_EUR)
summary(lookup$N) # 68 NAs; don't know why
#################################################

####### Examine SAIGE h2 ######
#################################################


####### Examine identifier and grouping of phe ######
## The first 7 column define a distinct trait
## phenocode: 
# for continuous, biomarkers, and categorical traits, this corresponds to the field ID as described by UKB, e.g. 21001 for BMI)
# for phecode, icd10, prescription this is unique
## coding: (only) For categorical variables, this corresponds to the coding that was used (e.g. coding 2 for field 1747). For all other trait_types, this field is blank
## modifer: eg. irnt for QT, version, etc
lookup %>% 
  distinct(trait_type, phenocode, pheno_sex, coding, modifier, description, description_more, coding_description)

table(lookup$modifier)

#################################################
######### Filter by sample size ################
## filter by sample size and h2
## note: #cases for QT are total sample size
## SAIGE h2 tends to underestimate h2
h2_pass <- 0 #0.01
N_pass <- 1000
trait_pass <- "PASS"

## 3813
phe_keep <- lookup %>% 
  filter(n_cases_EUR > N_pass)

dim(phe_keep)
summary(phe_keep$N)
table(phe_keep$phenotype_qc_EUR)

#################################################

######### Filter traits ################
## 831
pheQT <- phe_keep %>% 
  filter(BIN_QT == "QT") 

dim(pheQT)
table(pheQT$trait_type) 

#### biomarkers: 30 ["Biological samples > Assay results > Blood assays > Blood biochemistry"]
pheQT_biomarkers <- phe_keep %>% filter(trait_type == "biomarkers")
dim(pheQT_biomarkers)
pheQT_biomarkers %>% View

#### continuous: 801
pheQT_continuous <- phe_keep %>% filter(trait_type == "continuous")
dim(pheQT_continuous)
pheQT_continuous %>% View
##### --------------------------------------------------#####

## 2982
pheBIN <- phe_keep %>% 
  filter(BIN_QT == "BIN") 
dim(pheBIN)
table(pheBIN$trait_type)

##### icd10: 375, [chapters of ICD10]
pheBIN_icd10 <- pheBIN %>% filter(trait_type == "icd10")
dim(pheBIN_icd10)
pheBIN_icd10 %>% View

##### phecode: 662 [has category in disease classification]
pheBIN_phecode <- pheBIN %>% filter(trait_type == "phecode")
dim(pheBIN_phecode)
pheBIN_phecode %>% View

##### prescriptions: 373
pheBIN_prescriptions <- pheBIN %>% filter(trait_type == "prescriptions")
dim(pheBIN_prescriptions)
pheBIN_prescriptions %>% View

##### categorical: 1572
pheBIN_categorical <- pheBIN %>% filter(trait_type == "categorical")
dim(pheBIN_categorical)
pheBIN_categorical %>% View

# keep field: 20001, 20002, 20107, 20110, 20003 (267);
pheBIN_categorical <- pheBIN_categorical %>% 
  filter(phenocode %in% c("20001", "20002", "20107", "20110", "20003", "COVID19"))
dim(pheBIN_categorical)
pheBIN_categorical %>% View

##### --------------------------------------------------#####

#### Refine QT 
## QT: prefer outcome is irnt transformed: 831 --> 806
## remove Phenocode == "random"
pheQT_keep <- bind_rows(pheQT_biomarkers, pheQT_continuous) %>% 
  group_by(trait_type, phenocode, coding, description) %>% 
  add_count() %>% 
  mutate(to_rm = ifelse((n > 1 & !grepl("irnt", modifier)), 1, 0)) %>% 
  filter(to_rm == 0) %>% 
  ungroup() %>% 
  select(-c(to_rm, n)) %>% 
  filter(phenocode != "random") %>% 
  filter(!grepl("Genomics", category))

dim(pheQT_keep)

pheQT_keep %>% View
  # filter(grepl("Genomics", category))
  filter(phenocode==50) %>% View
  
pheQT_keep %>% filter(grepl("Believed safe to perform", description, ignore.case = T)) %>% View

##### --------------------------------------------------#####

#### Refine BIN 1677 --> 1677

pheBIN_keep <- bind_rows(pheBIN_icd10, pheBIN_phecode, pheBIN_prescriptions, pheBIN_categorical)
dim(pheBIN_keep)

## remove Genomics: this simply reflect comparison of self-report grouping with PCA grouping 2562 --> 2561
## To remove: Genomics --> 1677 (not selected in above step, so no reduction)
pheBIN_keep <- pheBIN_keep %>% filter(!grepl("Genomics", category))
dim(pheBIN_keep)

############ Write out phe list ############
## Note: put non-NA value on top, otherwise read_tsv will force all values to be NA
## 1677 BIN + 806 QT = 2483
allphe <- bind_rows(pheBIN_keep, pheQT_keep)
dim(allphe)
allphe %>% distinct(aws_link)

## final refine
allphe %>% filter(grepl("Population characteristics", category)) %>% pull(description_more)
allphe %>% filter(grepl("Recruitment > Reception", category, ignore.case = T)) %>% View
allphe %>% filter(grepl("Home location", category, ignore.case = T)) %>% View
allphe %>% filter(grepl("Health-related outcomes > Hospital inpatient", category, ignore.case = T)) %>% View

table(allphe$trait_type)
table(allphe$BIN_QT)

allphe %>% filter(BIN_QT == "QT") %>% View
allphe %>% filter(phenocode == "20003") %>% View

allphe <- read_tsv("./data/ukbb/PanUKBB/phe2483.PanUKBB.tsv")

## create grouping:
allphe <- allphe %>% 
  mutate(phe_group = ifelse(grepl("lifestyle | Diet | Touchscreen | environment", category, ignore.case = T) | phenocode %in% c("Smoking") , "LIFE", "NA")) %>% 
  mutate(phe_group = ifelse(grepl("physical measures | imaging", category, ignore.case = T) | phenocode %in% c("DBP","SBP", "MAP","MCP","PP", "whr", "eGFR", "eGFRcreacys","eGFRcys","FEV1FVC"), "PM", phe_group)) %>% 
  mutate(phe_group = ifelse(grepl("Cognitive | Mental", category, ignore.case = T) & trait_type=="continuous", "MENT", phe_group)) %>% 
  mutate(phe_group = ifelse(phenocode %in% c("20107", "20110"), "FH", phe_group)) %>% 
  mutate(phe_group = ifelse(trait_type %in% c("icd10", "phecode") | phenocode %in% c("20002", "COVID19"), "D", phe_group)) %>%
  mutate(phe_group = ifelse(trait_type != "continuous" & grepl("cancer", description, ignore.case = T) & !grepl("Non-cancer", description, ignore.case = T) & phenocode!="FH", "Cancer", phe_group)) %>% 
  mutate(phe_group = ifelse(phenocode %in% c("20003") | trait_type == "prescriptions", "MED", phe_group)) %>%
  mutate(phe_group = ifelse(grepl("assay", category, ignore.case = T)  | trait_type == "biomarkers" | phenocode %in% c("AG", "IBil", "LDLC","NAP"), "Assay", phe_group))


allphe %>% filter(phe_group == "PM") %>% 
  arrange(category) %>% view
  pull(category) %>% table()

allphe %>% select(N) %>%  
  add_count(N) %>% 
  mutate(prop = n/nrow(allphe)) %>% 
  distinct(N,.keep_all = T) %>% 
  select(-n) %>%
  write_tsv("./VB_add_N/panukbb_N.tsv")

## total 2483 phe
allphe %>% write_tsv("./data/ukbb/PanUKBB/phe2483.PanUKBB.tsv")
