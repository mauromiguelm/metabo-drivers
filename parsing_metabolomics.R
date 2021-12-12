## This code imports and clean metabolomics data

# load packages and definitions -------------------------------------------

path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\metabolomics'
path_metabolomics_in <- '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed'

library(openxlsx)
library(dplyr)
library(rhdf5)
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

tmp <- tmp[4:nrow(tmp),] #remove experiments that had problems with automatic liquid transfer

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

#removed unused source plates

source_plates <- source_plates[7:nrow(source_plates),]

rm(tmp, msp1, msp2)

#import drug source plate exceptions

source("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\exceptions\\log_processing.r")

# #import metabolomics data -----------------------------------------------

setwd(path_metabolomics_in)

dataContent<- h5ls("metabolomics_raw.h5")

metadata <- rhdf5::h5read(file = "metabolomics_raw.h5", 'samples/name')

metadata <- strsplit(metadata, "_")

metadata <- data.frame(do.call(rbind, metadata))
metadata$X3 <- NULL
colnames(metadata) <- c('injseq', 'source_plate','well384','well96','cell','drug',"conc")

metadata$idx <- 1:nrow(metadata)

# # remove controls  ------------------------------------------------------

#remove samples used to equilibrate system
metadata <- metadata[!grepl('burnin',metadata$cell),]

#remove samples used as batch controls
metadata <- metadata[!grepl('hct15Mtx',metadata$drug)& !grepl('poolCtrl',metadata$drug) & !grepl('hct15Ctrl',metadata$drug) & !grepl('SolvCrtl',metadata$drug),]

# # correct pipetting mistakes that can be corrected ----------------------
#when transfering source plate to 96 wp, there were a few pipetting mistakes
#all pipetting mistakes were recorded and can now be corrected

#SF539, IGROV1, MDAMB231: mistake intra data 384 H1 on P2Q3 row H of all CLs and 384
#P1 on P2Q3 row D of all CLs instead of original layout due to pipetting mistake

#apply fix to each cell line

metadata = switch_row(data = metadata,'SF529',plate = "P2",row_flip1 = "H",row_flip2 = "P",cols = seq(1,24,2))
metadata = switch_row(data = metadata,'IGROV1',plate = "P2",row_flip1 = "H",row_flip2 = "P",cols = seq(1,24,2))
metadata = switch_row(data = metadata,'MDAMB231',plate = "P2",row_flip1 = "H",row_flip2 = "P",cols = seq(1,24,2))

#COLO205	metabolomics: colo 205 cl3_p2_q1 row f on e and e on f

metadata = switch_row(metadata, 'COLO205',plate = "P2",row_flip1 = "F",row_flip2 = "E",cols = seq(1,24,2))

#HOP62	metabolomics: hop62 p1-q1:q4, A1 is H12, front == end

metadata = invert_plate(data = metadata, cell = "HOP62",plate = "P1",wells = paste0(rep(LETTERS[1:16][c(T,F)], each = 12), seq(1,24,2)))

# #remove echo problems ---------------------------------------------------
#remove wells that we had problems with drug transfer from echo pipetting

metadata$cell_plate <- paste(metadata$cell,metadata$source_plate)
source_plates$cell_plate <- unlist(lapply(strsplit(source_plates$filenames,"/"),"[[",3))

source_plates$cell_plate <- unlist(lapply(strsplit(source_plates$cell_plate,".txt"),"[[",1))

source_plates$cell_plate <- unlist(lapply(strsplit(source_plates$cell_plate,"_"),function(x){paste(x[[1]],x[[3]])}))

exception_wells <- exceptions[[2]]

exception_wells$DrugNameTransferError <- gsub(",",";",exception_wells$DrugNameTransferError)

as.data.frame(do.call(rbind, apply(exception_wells, 1, function(x) {
  do.call(expand.grid, strsplit(x, ","))
}))) -> exception_wells

exception_wells$source_well <- paste(exception_wells$cellPlateBC, exception_wells$WellNameTransferError)

metadata <- metadata %>% dplyr::right_join(source_plates[,c('uniqueID','cell_plate')], by = 'cell_plate')

metadata$source_well <- paste(metadata$uniqueID,metadata$well384)

exclusions_report <- subset(exception_wells, cellPlateBC%in% unique(metadata$uniqueID))

#these exclusions were already removed when measuring metabolomics, but report anyways
metadata <- subset(metadata, !source_well %in% unique(exclusions_report$source_well))

setwd(path_data_file)

write.csv(exclusions_report, 'exclusion_report_automatic_pipetting_error.csv')

#export results

setwd(paste(path_data_file,'metabolomics', sep = "//"))

write.csv(metadata,"metadata_clean.csv")



