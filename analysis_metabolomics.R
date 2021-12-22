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

ions <- rhdf5::h5read(file = "metabolomics_raw.h5", '/annotation')

ions <- data.frame(ions)

setwd(paste(path_data_file,'metabolomics', sep = "//"))

metadata <- read.csv("metadata_clean.csv")

drugs_in_screen <- c(unique(metadata$drug)[!(unique(metadata$drug) %in% c("PBS", "DMSO"))])

# calculating log2(FCs) for all data and for every drug.  --------

#iterating over plates, cells, drug, and conc to calculate FCs

metadata$cell_plate <- paste(metadata$cell, metadata$source_plate)

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))


lapply(unique(metadata$cell_plate), function(cell_plate_idx){
  #subset metadata to  drug & control groups
  #cell_plate_idx <- unique(metadata$cell_plate)[1]
  if (!any(grepl(cell_plate_idx,x=list.files()))){
    cell_plate_idx
    print(cell_plate_idx)
    results_list <- list()
    
    tmp <- subset(metadata, cell_plate == cell_plate_idx, select = c("idx","drug", "cell", "conc", "source_plate",'quadrant')) #keep one time point with GR24.. all time points have the same value..
    
    tmp_control <- subset(tmp, drug == "DMSO")
    
    print(paste("dimention of control vector =",dim(tmp_control)))
    
    tmp_drug <- subset(tmp, drug %in% drugs_in_screen)
    
    tmp_drug$drug_conc <- paste(tmp_drug$drug, tmp_drug$conc)
    
    #for each drug and concentration, calculate fold changes to control
    
    drug_conc <- unique(paste(tmp_drug$drug, tmp_drug$conc))
    #drug_conc_idx <- "Irinotecan 0.11"
    lapply(drug_conc, function(drug_conc_idx){
      #iterate over each drug_conc
      
      tmp_drug_conc <- subset(tmp_drug,drug_conc == drug_conc_idx)
      
      if(nrow(tmp_drug_conc)>0){
        
        lapply(1:nrow(data), function(metab_idx){
          #metab_idx = 1
          
          drug_metab <- data.frame("intensity"= data[metab_idx,tmp_drug_conc$idx])
          
          drug_metab$Group.1 <- tmp_drug_conc$quadrant
          
          control_metab <- aggregate(data[metab_idx,tmp_control$idx], by = list(tmp_control$quadrant),median)
          
          drug_fcs <- merge(drug_metab, control_metab, by = 'Group.1')
          
          drug_fcs <- log2(drug_fcs[,2]/ drug_fcs[,3])
          
          p_value <- t.test(drug_fcs, mu =0)$p.value
          
          return(list(paste(cell_plate_idx,drug_conc_idx, metab_idx), median(drug_fcs), p_value))
        })-> results
        
        return(results)
        
      }
      
      
    }) -> results
    
    metab_fcs <- unlist(results, recursive = FALSE)
    metab_fcs <- do.call(rbind, metab_fcs)
    
    metab_fcs <- data.frame(metab_fcs)
    
    metab_fcs<- cbind(data.frame(stringr::str_split_fixed(string = metab_fcs[,1]," ",n = 5)),metab_fcs[,2:3])
    
    metab_fcs[,6] <- as.numeric(metab_fcs[,6])
    metab_fcs[,7] <- as.numeric(metab_fcs[,7])
    
    write.csv(metab_fcs, paste0(cell_plate_idx,'.csv'))
    
  }else{
    NA
  }
})


# use linear regression to see if any metabolites have association with GR50
# lm(metab~GR50) across all cell lines and concentrations that we have filtered strong effects/unnefective concentrations
# do this wiht methotrexate as we have clear expectations


# Compare metabolomics results for R/Sgroups, and see which one is better

# plot the metabolites that are increased in R or S as a volcano plot



# check which metabolites are driving the correlation, hopefully positive controls
# Correlate FC with GR50 and GR24
# For metabolic drugs: check if direct substrate/products involved in MoA are relating to any drug metric 
# Multiomics: Link protein/gene/mRNA levels to metabolomics/drug sensitivity
# Multiomics: distance_to_target analysis