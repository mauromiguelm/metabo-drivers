### IMPORT PACKAGES ###
setwd('C:\\Users\\masierom\\polybox\\Data')
path_fig = path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures_mean'
library(ggplot2)
library(tidyr)
library(readxl)
library(KEGGREST)
library(dplyr)

### IMPORT DATA ###

drug_data <- list()

metab.data <- list()

drugs = read.csv('CANCER60GI50_updatedExperCpds.lst', stringsAsFactors = F) 

drugs$CONCUNIT[drugs$CONCUNIT == "log10(M)"] <- "M"

drugs$NLOGGI50 <- as.numeric(as.character(drugs$NLOGGI50))

drug_to_target <- read.csv('drugbank_all_target_polypeptide_ids.csv', stringsAsFactors = F)

db_to_pubchem_drug <-  read.csv('DB_Pubchem_drug links.csv', stringsAsFactors = F)

colnames(db_to_pubchem_drug)[3] <- 'CAS'

NSC_to_CAS  <-  read.csv('NSC_CAS_Sept2013.csv', header = F, stringsAsFactors = F)

#extract all sheets from excel file for metabolic data

read_excel_allsheets <- function(filename) {
  sheets <- readxl::excel_sheets(filename)
  x <-    lapply(sheets, function(X) readxl::read_excel(filename, sheet = X))
  names(x) <- sheets
  x
}

AZ_metab = read_excel_allsheets('AstraScreenPipe_metabolome_Improved.xlsx')


AZ_metab.ions = data.frame(AZ_metab$ions)
AZ_metab.intensities1 = data.frame(AZ_metab$intensities1)

rm(read_excel_allsheets)

# ### ORGANIZE DATA ### ---------------------------------------------------

# remove duplicated CAS Numbers in NSC_to_CAS, keep 1 NSC to 1 CAS, as it should be

colnames(NSC_to_CAS)[1] = 'NSC' 

colnames(NSC_to_CAS)[2] = 'CAS'

NSC_to_CAS <- NSC_to_CAS[NSC_to_CAS$NSC %in% unique(drugs$NSC),] #use this trick to avoid keeping unique values from duplicated in which one NSC is paired and the other is not to the drugs$NSC

NSC_to_CAS <- NSC_to_CAS[!(duplicated(NSC_to_CAS[2]) | duplicated(NSC_to_CAS[2], fromLast = F)), ] #FIXME since I dont know which is the corect NSC, first parse the NCI_60 DSD and then add these things

#including everolymus (733504) cas in NSC_to_CAS


NSC_to_CAS[NSC_to_CAS$CAS == '152459-95-5', 'CAS'] = '220127-57-1' #imatinib

NSC_to_CAS[NSC_to_CAS$CAS == '50-18-0', 'CAS'] = '6055-19-2' #cyclophosphamide

NSC_to_CAS[NSC_to_CAS$CAS == '37076-68-9', 'CAS'] = '17902-23-7' #tengafur

NSC_to_CAS[NSC_to_CAS$CAS == '444731-52-6', 'CAS'] = '635702-64-6' #Pazopanib

NSC_to_CAS[NSC_to_CAS$CAS == '849217-68-1', 'CAS'] = '1140909-48-3' #Cabozantinib

Cabozantinib <- data.frame(NSC = 761068, CAS = '1140909-48-3')

BPTES <- data.frame(NSC = 798303, CAS = '314045-39-1') 

CB839 <- data.frame(NSC = 783415, CAS = '1439399-58-2')

Metformin <- data.frame(NSC = 91485, CAS = '1115-70-4') 

Panzem <- data.frame(NSC = 659853, CAS = '362-07-2')

Geldanamycin <- data.frame(NSC = 330507, CAS = '75747-14-7')

YC_1 <- data.frame(NSC = 728165, CAS = '170632-47-0')

Ponatinib <- data.frame(NSC = 758487, CAS = '1114544-31-8')

Ixazomib <- data.frame(NSC = 758254, CAS = '1239908-20-3')

Osimertinib <- data.frame(NSC = 779217, CAS = '1421373-66-1')

Alectinib <- data.frame(NSC = 764040, CAS = '1256580-46-7')

Lenvatinib <- data.frame(NSC = 755980, CAS = '417716-92-8')

Ceritinib <- data.frame(NSC = 776422, CAS = '1032900-25-6')

Belinostat <- data.frame(NSC = 758774, CAS = '414864-00-9')

Ibrutinib <- data.frame(NSC = 761910, CAS = '936563-96-1')

Capecitabine <- data.frame(NSC = 712807 , CAS = '154361-50-9')

Gemcitabine <- data.frame(NSC = 613327, CAS = '122111-03-9')

pemetrexed_disodium <- data.frame(NSC = 698037, CAS = '150399-23-8')

lapatinib <- data.frame(NSC = 745750, CAS = '388082-78-8')

erlotinib <- data.frame(NSC = 718781, CAS = '183319-69-9')

everolymus <- data.frame(NSC = 733504, CAS = '159351-69-6')

temsirolimus <- data.frame(NSC = 683864, CAS = '162635-04-3')

Bosutinib <- data.frame(NSC = 765694, CAS = '380843-75-4')

Regorafenib <- data.frame(NSC = 763932, CAS = '755037-03-7')

Panobinostat <- data.frame(NSC = 761190, CAS = '404950-80-7')

NSC_to_CAS <- rbind(NSC_to_CAS, Panobinostat, Regorafenib, Bosutinib, everolymus, temsirolimus, erlotinib, lapatinib, pemetrexed_disodium, Gemcitabine, Capecitabine,
                    Ibrutinib, Belinostat, Ceritinib, Lenvatinib, Alectinib, Osimertinib, Ixazomib, Ponatinib, Cabozantinib, BPTES, CB839, Metformin,
                    Panzem, Geldanamycin, YC_1)

rm(Bosutinib, everolymus, temsirolimus, erlotinib, lapatinib, pemetrexed_disodium, Gemcitabine, Capecitabine,
   Ibrutinib, Belinostat, Ceritinib, Lenvatinib, Alectinib, Osimertinib, Ixazomib, Ponatinib, Cabozantinib, Regorafenib, Panobinostat, BPTES, CB839, Metformin,
   Panzem, Geldanamycin, YC_1)

# remove speciec !human in drug_to_target

drug_to_target = drug_to_target[drug_to_target$Species == 'Human',]

# merging drug IDs and targets prior to merging with datasets 

colnames(drug_to_target)[13] = 'DrugBank.ID'

colnames(drug_to_target)[2] = 'Target_name'

drug_to_target = drug_to_target %>% #split the target into columns
  mutate(DrugBank.ID = strsplit(DrugBank.ID,';')) %>%
  unnest(DrugBank.ID)

drug_to_target$DrugBank.ID = gsub(" ", "", drug_to_target$DrugBank.ID) # remove whitespaces after the drug_to_target strsplit

#remove cell lines with replicates falling in more than one cluster and general data wrangling

drugs$CELL = gsub(" ", "", drugs$CELL) #remove whitespaces

drugs$CELL = gsub("-", "", drugs$CELL) #remove hyphens 


# Merge replicates on metabolome dataset 

#AZ_metab.intensities = cbind(AZ_metab.intensities1, AZ_metab.intensities2[4:248])

AZ_metab.intensities = AZ_metab.intensities1

AZ_metab.intensities$ionActive = NULL
AZ_metab.intensities$ionMz = NULL


colnames(AZ_metab.intensities)[2] = '184A1'
colnames(AZ_metab.intensities)[3] = '5637'
colnames(AZ_metab.intensities)[4] = '647V'
colnames(AZ_metab.intensities)[5] = '769P'
colnames(AZ_metab.intensities)[6] = '786O'
colnames(AZ_metab.intensities)[7] = '977'

rm(AZ_metab.intensities1, AZ_metab)



# transpose data frame 

metab.data$dataset = AZ_metab.intensities

metab.data$dataset = data.frame(t(metab.data$dataset), colname = T)

colnames(metab.data$dataset) = metab.data$dataset[1, ]

metab.data$dataset = metab.data$dataset[-1,]

metab.data$dataset = tibble::rownames_to_column(metab.data$dataset, "CELL")

metab.data$dataset$CELL = factor(metab.data$dataset$CELL)

rm(AZ_metab.intensities)

rm(db_to_pubchem_drug, drug_to_target)

metab.data$dataset$CELL = gsub("786O", "7860", metab.data$dataset$CELL) 

metab.data$dataset$"TRUE" <- NULL

#merge datapoints in which there is more than one concentration of compound per cell line

drugs$CELL = gsub("/ATCC", "", drugs$CELL) 

drugs$CELL = gsub(".Fine", "", drugs$CELL) 

drugs$CELL = gsub("/", "", drugs$CELL) 

drug_sensitivity_NCI = drugs %>% 
  group_by(NSC, CELL, CONCUNIT) %>%
  dplyr::mutate(aNLOGGI50 = mean(NLOGGI50) ) %>%
  ungroup()%>%
  dplyr::mutate(NSC = factor(NSC))%>%
  dplyr::mutate(CELL = factor(CELL))%>%
  group_by(NSC, CELL, CONCUNIT)%>%
  dplyr::slice(1)%>%
  ungroup()%>%
  dplyr::select(-c(INDN, TOTN, STDDEV))


rm(drugs, metab.filter, metab, db_to_pubchem_drug, drug_to_target)

#remove u (g/mL) from actinomycin NSC = 3053 and keep only [Molar]

drug_sensitivity_NCI <- drug_sensitivity_NCI[!(drug_sensitivity_NCI$NSC == '3053' &
                                                               drug_sensitivity_NCI$CONCUNIT == 'u ') ,]


#import kegg Cancer drugs

kEGG_Cancer_drug_set = read.csv('Drugs_MoA_Kegg_Antineoplastics.csv', stringsAsFactors = F)

tmp_store_DrugBankID_CAS <- list()

library(KEGGREST)

for (i in 1:length(kEGG_Cancer_drug_set$KEGG)){
  
  k = kEGG_Cancer_drug_set$KEGG[i]

  tmp_getDrug <- keggGet(k)

  tmp_getDrug_2 <- tmp_getDrug[[1]]$DBLINKS
  
  tmp_store_DrugBankID_CAS[i] <- list(c(ifelse(length(tmp_getDrug_2[grepl(pattern = 'CAS', x = tmp_getDrug_2)]) == 0 , 'NA', tmp_getDrug_2[grepl(pattern = 'CAS', x = tmp_getDrug_2)]), 
                                   ifelse(length(tmp_getDrug_2[grepl(pattern = 'PubChem', x = tmp_getDrug_2)]) == 0 , 'NA', tmp_getDrug_2[grepl(pattern = 'PubChem', x = tmp_getDrug_2)])))
  }

tmp_store_DrugBankID_CAS.df <- t(data.frame(tmp_store_DrugBankID_CAS))

kEGG_Cancer_drug_set <- cbind(kEGG_Cancer_drug_set, tmp_store_DrugBankID_CAS.df)

kEGG_Cancer_drug_set$`1` <- as.character(kEGG_Cancer_drug_set$`1`)

kEGG_Cancer_drug_set$`2` <- as.character(kEGG_Cancer_drug_set$`2`)

colnames(kEGG_Cancer_drug_set)[9:10] <- c('CAS','Pubchem')

kEGG_Cancer_drug_set$CAS <- gsub(pattern = "CAS: ", replacement = "", x = kEGG_Cancer_drug_set$CAS) 

kEGG_Cancer_drug_set$Pubchem <- gsub(pattern = "PubChem: ", replacement = "", x = kEGG_Cancer_drug_set$Pubchem) 

rm(k, i, tmp_getDrug, tmp_getDrug_2, tmp_store_DrugBankID_CAS, tmp_store_DrugBankID_CAS.df)

rownames(kEGG_Cancer_drug_set) <- NULL

#fix drugs IDs in kegg metadata  

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Nilotinib', 'CAS'] = '641571-10-0'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Thioguanine', 'CAS'] = '154-42-7'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Trifluridine, mixt', 'CAS'] = '70-00-8'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Panobinostat', 'CAS'] = '404950-80-7'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Afatinib', 'CAS'] = '439081-18-2'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Bosutinib', 'CAS'] = '380843-75-4'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Regorafenib', 'CAS'] = '755037-03-7'

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == 'Alectinib', 'CAS'] = '1256580-46-7'


# include experimental drugs into kegg_Cancer drug set

# [1] "Class"             "Target"            "Generic.name"      "ATC"               "KEGG"              "Products..USA."    "Products..Japan." 
# [8] "Indications..USA." "CAS"               "Pubchem"        


experimental_drugs <- data.frame(Class =             c("Fatty acid biosynthesis" ,"Fatty acid biosynthesis"      ,"Molecularly targeted agent"     ,"Fatty acid biosynthesis"    ,"OXPHOS"    , "Glutamine inhibitor" , "Glutamine inhibitor" , "OXPHOS"              , "OXPHOS"                      , "HIF-1alpha"          , "HSP90"       , "HIF-1alpha"),
                                 Target =            c("CPT1"                    ,"CPT1"                         ,"mTOR kinase inhibitor"          ,'ATP citrate lyase'          ,"MPC"       ,"Glutaminase GLS1"     , "Glutaminase GLS1"    , "ATP Synthase"        , "Mitochondrial complex I"     , "HIF-1alpha"          , "HSP90"       , "HIF-1alpha"),
                                 Generic.name =      c("Oxfenicine"              ,"Etomoxir"                     ,"Rapamycin"                      ,"MEDICA-16"                  ,"UK-5099"   ,"BPTES"                , "CB-839"              , "Oligomycin A"        , "Metformin"                   , "Panzem (2-ME)"       , "17-AAG"      , "YC-1"),  
                                 ATC =               c(NA),
                                 KEGG =              c(NA),
                                 Products..USA. =     c(NA),
                                 Products..Japan. =  c(NA),
                                 Indications..USA. =  c(NA),
                                 CAS =               c("32462-30-9"              ,"828934-41-4"                 ,"53123-88-9"                     ,"87272-20-6"                 ,"56396-35-1" ,"314045-39-1"          , "1439399-58-2"        , "579-13-5"             , "1115-70-4"                  , "362-07-2"            , "75747-14-7"  , "170632-47-0"),
                                 Pubchem =           c(NA))

kEGG_Cancer_drug_set <- rbind(kEGG_Cancer_drug_set, experimental_drugs)


# Include whether drug is in screen or not

kEGG_Cancer_drug_set$is_Screen <- 0

kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Erlotinib" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Lenvatinib" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Irinotecan" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Topotecan" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Clofarabine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Cladribine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Mercaptopurine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Fluorouracil" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Gemcitabine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Decitabine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Pemetrexed" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Methotrexate" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Docetaxel" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Paclitaxel" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Everolimus" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Temsirolimus" ,"is_Screen"] <- 0
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Chlormethine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Oxaliplatin" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Asparaginase" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Omacetaxine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "BPTES" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "CB-839" ,"is_Screen"] <- 0  #removed from 1 to 0, 20220613
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Oligomycin A" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Shikonin" ,"is_Screen"] <- 0  #this guy is not on screen, removed from 1 to 0, 20220613
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Panzem (2-ME)" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "17-AAG" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "YC-1" ,"is_Screen"] <- 1

# new additions 20220613
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "UK-5099" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "MEDICA-16" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Rapamycin" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Metformin" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Etomoxir" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Oxfenicine" ,"is_Screen"] <- 1
kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$Generic.name == "Trametinib" ,"is_Screen"] <- 1


# export Cancer drug indications as drug metadata 

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\drug-metadata")

tmp <- subset(kEGG_Cancer_drug_set, is_Screen == 1)

tmp = tmp[!(duplicated(tmp$Generic.name) #remove duplicated with multiple inducations 
                                              | duplicated(tmp$Generic.name, fromLast = F)), ]

write.csv(tmp,"drug_metadata_moa.csv")

#plot drug classes

setwd(path_fig)

tmp <- read.csv("drug_metadata_moa.csv")

ggplot(tmp, aes(Class))+
  geom_histogram(stat="count")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,vjust=0,hjust=1))->plt

ggsave('drugs_by_moa.pdf',plt,height = 5,width = 5)
ggsave('drugs_by_moa.png',plt,height = 5,width = 5)

# add CAS to drug_sensitivity_NCI



library(dplyr)

kEGG_Cancer_drug_set <- inner_join(NSC_to_CAS, kEGG_Cancer_drug_set, by = 'CAS')

#experimental_drugs[!(experimental_drugs$CAS %in% kEGG_Cancer_drug_set$CAS),]

kEGG_Cancer_drug_set = kEGG_Cancer_drug_set[!(duplicated(kEGG_Cancer_drug_set[1]) #remove duplicated with multiple inducations 
                                                            | duplicated(kEGG_Cancer_drug_set[1], fromLast = F)), ]


kEGG_Cancer_drug_set <- kEGG_Cancer_drug_set[kEGG_Cancer_drug_set$NSC %in% unique(drug_sensitivity_NCI$NSC),] #remove the ones that are not in drug_sens

kEGG_Cancer_drug_set = kEGG_Cancer_drug_set[!(duplicated(kEGG_Cancer_drug_set[5]) #remove duplicates, both drugs had corr 0.7
                                                             | duplicated(kEGG_Cancer_drug_set[5], fromLast = F)), ]


drug_sensitivity_NCI_cancer <- drug_sensitivity_NCI[drug_sensitivity_NCI$CELL %in% unique(metab.data$dataset$CELL),]

drug_sensitivity_NCI_cancer <- drug_sensitivity_NCI_cancer[drug_sensitivity_NCI_cancer$NSC %in% unique(kEGG_Cancer_drug_set$NSC),]

rm(NSC_to_CAS)

metab.data$data_ions <- AZ_metab.ions

rm(AZ_metab.ions)

drug_data$Kegg_metadata <- kEGG_Cancer_drug_set
drug_data$NCI60_GI50_cancer <- drug_sensitivity_NCI_cancer
drug_data$exp_drugs <- experimental_drugs

rm(drug_sensitivity_NCI, drug_sensitivity_NCI_cancer, kEGG_Cancer_drug_set, experimental_drugs)


# #remove "X001..." from "X001...NCIH226_MTX_4_t1" in the sample name.
# 
# metab.data$dataset$CELL <- 
#   
# lapply(strsplit(as.character(metab.data$dataset$CELL), split = "...", fixed = T), function(x) x[2])
# 
# 
# metab.data$dataset$CELL


