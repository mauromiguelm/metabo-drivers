## This code imports and clean metabolomics data

# load packages nad definitions -------------------------------------------

path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\metabolomics'
path_metabolomics_in <- '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed'

library(openxlsx)
source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\plate_converter.R")

#import drug source plate map

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\growthData")

source_plates <- data.frame(
  filenames = list.files(
    path = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\Growthdata\\",
    pattern = "[P][1-2][.txt]",
    recursive = T))

source_plates$sourceid <- ifelse(grepl(pattern = "P1", source_plates$filenames), "1MSP001", "2MSP001")

source_plates$batch <- ifelse(grepl(pattern = "20190513", source_plates$filenames), "batch_1", "batch_2")

setwd("..")

tmp <- read.xlsx("Description.xlsx")

source_plates$uniqueID <- NA

cell <- strsplit(x = as.character(source_plates$filenames), split = "/")

cell <- lapply(cell, "[[", 3)

cell <- strsplit(x = unlist(cell), split = "_")

plate <- unlist(lapply(cell, "[", 3))

cell <- unlist(lapply(cell, "[[", 1))

tmp <- tmp[match(x = cell, table = tmp$cell_name), ]

msp1 <- grepl(pattern = "1MSP", x = source_plates$sourceid)

source_plates$uniqueID[msp1] <- tmp$source1[msp1]

msp2 <- grepl(pattern = "2MSP", x = source_plates$sourceid)

source_plates$uniqueID[msp2] <- tmp$source2[msp2]

rm(tmp, msp1, msp2)

#import drug source plate exceptions

source("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\exceptions\\log_processing.r")

#import metabolomics data

setwd(path_metabolomics_in)

# # map 96wp metabolomics to 384 source -----------------------------------
#each 384 plate has four 96 source plates, to a total of 8x 96 well plates
#used for metabolomics


# #remove echo problems ---------------------------------------------------
#remove wells that we had problems with drug transfer from echo pipetting



# # correct pipetting mistakes that can be corrected ----------------------
#when transfering source plate to 96 wp, there were a few pipetting mistakes
#all pipetting mistakes were recorded and can now be corrected


#SF539	mistake intra data 384 H1 on P2Q3 row H of all CLs and 384
#P1 on P2Q3 row D of all CLs instead of original layout due to pipetting mistake

#TODO

#IGROV1	mistake intra data 384 H1 on P2Q3 row H of all CLs and 384
#P1 on P2Q3 row D of all CLs instead of original layout due to pipetting mistake; 

#TODO

#MDAMB231	mistake intra data 384 H1 on P2Q3 row H of all CLs and 384 P1
#on P2Q3 row D of all CLs instead of original layout due to pipetting mistake

#TODO

#HOP62	metabolomics: hop62 p1-q1:q4, A1 is H12, front == end

#TODO flip plates completely

#COLO205	metabolomics: colo 205 cl3_p2_q1 row f on e and e on f

#TODO switch rows 


# # remove pipetting mistakes that cannot be corrected  -------------------


# # filter ions based on CV -----------------------------------------------

# # filter samples based on CV --------------------------------------------


#export results