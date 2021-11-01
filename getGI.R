## This code calculates the GI values across all concentrations within a growth curve.


### THIS code only works with commit 

load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\growth_curves.RData")

library(dplyr); library(tidyr); library(ggplot2)



# GR metrics, their example code only use the endpoint and time at treatment, but I will try all time points

# we need to select time point zero as the time we treated the samples (as per their requirements)

devtools::load_all("C:\\Mauro_r_library\\GRmetrics")

devtools::load_all("C:\\Mauro_r_library\\tdsR")

time_treatment <- 0

drugs_in_screen <- c(unique(data_corrected$Drug[!(data_corrected$Drug %in% c("PBS", "DMSO", NA, "Water", "WATER", "H2O", "Dmso", "Control", "Ctrl", "exception"))]))

#matrix()   # columns will be drugs, rows will be cell lines

data_GRmetrics <- data_corrected

#cell_line <- "NCIH460"
#drug  <- "Docetaxel"
#time_treatment <- 24

data_Grmetrics_Ttm <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # preparing the data for one drug
data_Grmetrics_Ctr <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # getting the matching controls ready

data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Drug %in% drugs_in_screen) #I chose clofarabine since it has a pretty dose response curve
data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Drug == "DMSO" & Final_conc_uM %in% c(333, 367)) # I separated the drug from the control, as its easier to parse

# keep only matching time points

match_time_intervals <- intersect(data_Grmetrics_Ctr$Time, data_Grmetrics_Ttm$Time) #make sure both datasets cover the same time

time_treatment <- max(match_time_intervals[match_time_intervals < time_treatment])

data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time %in% match_time_intervals)
data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time %in% match_time_intervals)

# Keeping certain constraints on time, so we are sure conf is not limiting (eg. >80)
# also, GRmetrics require the time point at treatment, and not before. So zero will be time point slightly before treatment

data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time >= time_treatment & Time == 72 | Time == 0)
data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time >= time_treatment & Time == 72 | Time == 0)

data_Grmetrics_Ctr$Time <- data_Grmetrics_Ctr %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time
data_Grmetrics_Ttm$Time <- data_Grmetrics_Ttm %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time


#data_Grmetrics_Ctr$Time <- ifelse(data_Grmetrics_Ctr$Time == min(data_Grmetrics_Ctr$Time), 0, data_Grmetrics_Ctr$Time) #FIXME for each cell line get a specific time, instead of general
#data_Grmetrics_Ttm$Time <- ifelse(data_Grmetrics_Ttm$Time == min(data_Grmetrics_Ttm$Time), 0, data_Grmetrics_Ttm$Time)

#TODO select only one time point, eg. 60h, and then calculate the GRmetric across this time point.. could be even 48h

# Prepare data_Grmetrics_Ctr as example Case C

data_Grmetrics_Ctr$time = data_Grmetrics_Ctr$Time

tmp <- strsplit(x = rownames(data_Grmetrics_Ctr), "_")

tmp <- lapply(tmp, "[[",1)

tmp <- gsub(tmp, pattern = "/", replacement = "_")

tmp <- strsplit(x = tmp, "_")

tmp <- lapply(tmp, "[[",3)

data_Grmetrics_Ctr$cell_line <- tmp

data_Grmetrics_Ctr$agent = "-"

data_Grmetrics_Ctr$concentration = 0

data_Grmetrics_Ctr$cell_count <- data_Grmetrics_Ctr$Conf

data_Grmetrics_Ctr <- data_Grmetrics_Ctr[,c("cell_line", "agent","time", "concentration", "cell_count")]

rownames(data_Grmetrics_Ctr) <- NULL

# Prepare data_Grmetrics_Ttm as example Case C

head(data_Grmetrics_Ttm)

data_Grmetrics_Ttm$time = data_Grmetrics_Ttm$Time

tmp <- strsplit(x = rownames(data_Grmetrics_Ttm), "_")

tmp <- lapply(tmp, "[[",1)

tmp <- gsub(tmp, pattern = "/", replacement = "_")

tmp <- strsplit(x = tmp, "_")

tmp <- lapply(tmp, "[[",3)

data_Grmetrics_Ttm$cell_line <- tmp; rm(tmp)

data_Grmetrics_Ttm$agent <- (data_Grmetrics_Ttm$Drug)

data_Grmetrics_Ttm$concentration = data_Grmetrics_Ttm$Final_conc_uM

data_Grmetrics_Ttm$cell_count <- data_Grmetrics_Ttm$Conf

data_Grmetrics_Ttm <- data_Grmetrics_Ttm[,c("cell_line", "agent","time", "concentration", "cell_count")]

rownames(data_Grmetrics_Ttm) <- NULL

# merge drug df and control df

data_comb <- rbind(data_Grmetrics_Ctr, data_Grmetrics_Ttm)

#data_comb <-subset(data_comb, cell_line =="HT29" & agent %in% c("-", "Chlormethine", "BPTES"))

data_comb$cell_line <- as.character(data_comb$cell_line)

#library(GRmetrics)

output1 = GRmetrics::GRfit(inputData = data_comb, groupingVariables = 
                              c("cell_line", 'agent'), case = "C")

#calculate average GI, with sd

output1 <- output1 %>% group_by(cell_line, agent, concentration) %>% dplyr::summarize(GI_sd = sd(norm_rel_cell_count, na.rm = T)/sqrt(n()),
                                                                                      GI_mean = mean(norm_rel_cell_count, na.rm=T))



## import tdsR data ####

library(dplyr)

load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\tdsR_output.Rdata")

output_tdsR$experiment <- rownames(output_tdsR)

output1$experiment <- base::do.call(paste, output1[,c("cell_line", "agent")])

output_tdsR[,c("cell_line", "agent")] <- NULL

colnames(output_tdsR)[2] <- "tds"

output_tdsR$tds <- as.numeric(as.character(output_tdsR$tds))



output_comb <- output1 %>% dplyr::inner_join(output_tdsR, by = "experiment")

#### calculate corrected GI50 ########

tmp <- output_comb

tmp <- do.call(data.frame,lapply(tmp, function(x) replace(x, is.infinite(x),NA)))

data_GRmetrics <- data_corrected

#cell_line <- "NCIH460"
#drug  <- "Docetaxel"
#time_treatment <- 24

data_Grmetrics_Ttm <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # preparing the data for one drug
data_Grmetrics_Ctr <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # getting the matching controls ready

data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Drug %in% drugs_in_screen) #I chose clofarabine since it has a pretty dose response curve
data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Drug == "DMSO" & Final_conc_uM %in% c(333, 367)) # I separated the drug from the control, as its easier to parse

# keep only matching time points

match_time_intervals <- intersect(data_Grmetrics_Ctr$Time, data_Grmetrics_Ttm$Time) #make sure both datasets cover the same time

time_treatment <- max(match_time_intervals[match_time_intervals < time_treatment])

data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time %in% match_time_intervals)
data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time %in% match_time_intervals)

# Keeping certain constraints on time, so we are sure conf is not limiting (eg. >80)
# also, GRmetrics require the time point at treatment, and not before. So zero will be time point slightly before treatment

lapply(split(tmp, tmp$experiment), function(idx){
  
  #idx <- split(tmp, tmp$experiment)[[1]]
  
  if(idx$tds <= 60){
    
    
    data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, cell == idx$cell_line)
    data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, cell == idx$cell_line & Drug == idx$agent)
    
  
    data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time >= time_treatment & Time == 72 | Time == idx$tds)
    data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time >= time_treatment & Time == 72 | Time == idx$tds)
    
    data_Grmetrics_Ctr$Time <- ifelse(data_Grmetrics_Ctr$Time == min(data_Grmetrics_Ctr$Time), 0, data_Grmetrics_Ctr$Time) #FIXME for each cell line get a specific time, instead of general
    data_Grmetrics_Ttm$Time <- ifelse(data_Grmetrics_Ttm$Time == min(data_Grmetrics_Ttm$Time), 0, data_Grmetrics_Ttm$Time)
    
    
    data_Grmetrics_Ctr$time = data_Grmetrics_Ctr$Time
    
    tmp <- strsplit(x = rownames(data_Grmetrics_Ctr), "_")
    
    tmp <- lapply(tmp, "[[",1)
    
    tmp <- gsub(tmp, pattern = "/", replacement = "_")
    
    tmp <- strsplit(x = tmp, "_")
    
    tmp <- lapply(tmp, "[[",3)
    
    data_Grmetrics_Ctr$cell_line <- tmp
    
    data_Grmetrics_Ctr$agent = "-"
    
    data_Grmetrics_Ctr$concentration = 0
    
    data_Grmetrics_Ctr$cell_count <- data_Grmetrics_Ctr$Conf
    
    data_Grmetrics_Ctr <- data_Grmetrics_Ctr[,c("cell_line", "agent","time", "concentration", "cell_count")]
    
    rownames(data_Grmetrics_Ctr) <- NULL
    
    data_Grmetrics_Ttm$time = data_Grmetrics_Ttm$Time
    
    tmp <- strsplit(x = rownames(data_Grmetrics_Ttm), "_")
    
    tmp <- lapply(tmp, "[[",1)
    
    tmp <- gsub(tmp, pattern = "/", replacement = "_")
    
    tmp <- strsplit(x = tmp, "_")
    
    tmp <- lapply(tmp, "[[",3)
    
    data_Grmetrics_Ttm$cell_line <- tmp; rm(tmp)
    
    data_Grmetrics_Ttm$agent <- (data_Grmetrics_Ttm$Drug)
    
    data_Grmetrics_Ttm$concentration = data_Grmetrics_Ttm$Final_conc_uM
    
    data_Grmetrics_Ttm$cell_count <- data_Grmetrics_Ttm$Conf
    
    data_Grmetrics_Ttm <- data_Grmetrics_Ttm[,c("cell_line", "agent","time", "concentration", "cell_count")]
    
    rownames(data_Grmetrics_Ttm) <- NULL
    
    # merge drug df and control df
    
    data_comb <- rbind(data_Grmetrics_Ctr, data_Grmetrics_Ttm)
    
    #data_comb <-subset(data_comb, cell_line =="HT29" & agent %in% c("-", "Chlormethine", "BPTES"))
    
    data_comb$cell_line <- as.character(data_comb$cell_line)
    
    output1_corrected = GRmetrics::GRfit(inputData = data_comb, groupingVariables =
                                 c("cell_line", 'agent'), case = "C")
    
    
    return(output1_corrected)
    
  }else{
    
    output1_corrected <- NA
    
  }
  
}) -> tmp

output2 <- do.call(rbind, tmp)

output2$group <- do.call(paste, output2[,c("cell_line", "agent", "concentration")])

output2 <- output2 %>% group_by(group) %>% dplyr::summarize(GI_sd_new = sd(norm_rel_cell_count, na.rm = T)/sqrt(n()),
                                                                                      GI_mean_new = mean(norm_rel_cell_count, na.rm=T))



output2 <- output2[,c("group", "GI_sd_new","GI_mean_new")]

output1$group <- do.call(paste, output1[,c("cell_line", "agent", "concentration")])

output_gi <- output1 %>% left_join(output2, by = "group")

#FIXME sclae concentration over drugs

ggplot(output_gi, aes(x = GI_mean, y = GI_mean_new, col = (concentration)))+
  geom_point()

