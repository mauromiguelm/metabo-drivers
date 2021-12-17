## This code analyse metabolomics data

# load packages and definitions -------------------------------------------

path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\metabolomics'
path_metabolomics_in <- '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed'

source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\plate_converter.R")
library(openxlsx)
library(dplyr)
library(rhdf5)


setwd(path_metabolomics_in)

dataContent<- h5ls("metabolomics_raw.h5")

data<- rhdf5::h5read(file = "metabolomics_raw.h5", '/data')

setwd(paste(path_data_file,'metabolomics', sep = "//"))

metadata <- read.csv("metadata_clean.csv")

drugs_in_screen <- c(unique(metadata$drug)[!(unique(metadata$drug) %in% c("PBS", "DMSO"))])

# calculating log2(FCs) for all data and for every drug.  --------

#define 96 well plates within each 384 wp
#P1-Q1 to P1-Q4 

quadrant_map <- lapply(1:4,get_quadrant_wells)

quadrant_map <- do.call(rbind,quadrant_map)

lapply(unique(metadata$source_plate), function(plate_idx){
  #plate_idx = "P1"
  tmp <- subset(metadata, source_plate == plate_idx)
  
  tmp <- merge(tmp, quadrant_map[,c("well384",'quadrant')], by = "well384")
  
  return(tmp)
  
})

#iterating over plates, cells, drug, and conc to calculate FCs

metadata$cell_plate <- paste(metadata$cell, metadata$source_plate)


lapply(unique(metadata$cell_plate), function(cell_plate_idx){
  #subset metadata to  drug & control groups
  #cell_plate_idx = 'MDAMB231 P1'
  #TODO isolate control by 96wp, now they are grouped in a large batch of 384
  print(cell_plate_idx)
  results_list <- list()
  
  tmp <- subset(metadata, cell_plate == cell_plate_idx, select = c("idx","drug", "cell", "conc", "source_plate")) #keep one time point with GR24.. all time points have the same value..
  
  tmp_control <- subset(tmp, drug == "DMSO")
  
  tmp_drug <- subset(tmp, drug %in% drugs_in_screen)
  
  tmp_drug$drug_conc <- paste(tmp_drug$drug, tmp_drug$conc)
  
  #for each drug and concentration, calculate fold changes to control
  
  drug_conc <- unique(paste(tmp_drug$drug, tmp_drug$conc))
  #drug_conc_idx <- "Irinotecan 0.11"
  lapply(drug_conc, function(drug_conc_idx){
    #iterate over each drug_conc
    
    tmp_drug_conc <- subset(tmp_drug,drug_conc == drug_conc_idx)

    lapply(1:nrow(data), function(metab_idx){
      #metab_idx = 1
      drug_metab <- median(data[metab_idx,tmp_drug_conc$idx])
      control_metab <- median(data[metab_idx,tmp_control$idx])
      return(list(paste(cell_plate_idx,drug_conc_idx, metab_idx), log2(drug_metab/control_metab)))
    }) -> results
    
    results_list <- append(results_list,results)
    
  })
  return(results_list)
}) -> metab_fcs

list <- unlist(metab_fcs, recursive = FALSE)
list <- unlist(list, recursive = FALSE)
data <- do.call(rbind, list)


# Compare metabolomics results for R/Sgroups, and see which one is better
# check which metabolites are driving the correlation, hopefully positive controls
# Correlate FC with GR50 and GR24
# For metabolic drugs: check if direct substrate/products involved in MoA are relating to any drug metric 
# Multiomics: Link protein/gene/mRNA levels to metabolomics/drug sensitivity
# Multiomics: distance_to_target analysis