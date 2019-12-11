#source data

library(GRmetrics); library(dplyr); library(tidyr)

source("\\\\d.ethz.ch/groups/biol/sysbc/sauer_1/users/Mauro/Cell_culture_data/190310_LargeScreen/analysis/parsing-growthCurves.R")

# GRmetrics ---------------------------------------------------------------


# GR metrics, their example code only use the endpoint and time at treatment, but I will try all time points

# we need to select time point zero as the time we treated the samples (as per their requirements)


time_treatment <- 24

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

data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time >= time_treatment & Time == 72 | Time == 24)
data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time >= time_treatment & Time == 72 | Time == 24)

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

data_comb <-subset(data_comb, cell_line =="HT29" & agent %in% c("-", "Chlormethine", "BPTES"))

data_comb$cell_line <- as.character(data_comb$cell_line)

output1 = GRfit(inputData = data_comb, groupingVariables = c("cell_line", 'agent'), case = "C")

GRdrawDRC(output1, points =F)

GRbox(output1, metric ='GR50', groupVariable = c('cell_line'), 
      pointColor = c("agent"))

GRbox(output1, metric ='GR50', groupVariable = c('agent'), 
      pointColor = c("cell_line"))

output2 <- GRmetrics::GRgetMetrics(output1)

cc <- scales::seq_gradient_pal("blue", "red")(seq(0,1,length.out= 18))

# run tdsR on dataset, and compare with 

### ### CONTINUE HERE ### ### 
### ### CONTINUE HERE ### ### 
### ### CONTINUE HERE ### ### 
### ### CONTINUE HERE ### ### 
### ### CONTINUE HERE ### ### 

#FIXME sometimes the fit return a p-value instead of NA, include p-value in the control flow

tdsR_fit(inputData = data_comb, groupingVariables = c("cell_line", "agent", "time"))

#cc <- grDevices::colorRampPalette(colors = c("blue", "red"))

data_figure <- subset(data_corrected, Drug %in% c("DMSO",drugs_in_screen))

data_figure <- data_figure[!(as.character(data_figure$Final_conc_uM) %in% c("33.3", "3.33")),]

data_figure$Final_conc_uM <- ifelse(data_figure$Final_conc_uM == 333, 0, data_figure$Final_conc_uM)


ggplot(subset(subset(data_figure, cell == "NCIH460"), Time > 24 & Time < 72 ) , aes(x = Time, y = Conf, color = factor(Final_conc_uM)))+
  geom_smooth()+
  facet_grid(~cell)+
  scale_color_manual(values = cc)+
  guides(color=guide_legend(title="Concentration (uM)"))


plot(output2$time, output2$GR50, xlab = "Time (h)", ylab = "GR50", col = output2$cell_line)



# # GR analysis w.r.t. time -----------------------------------------------


# heatmap GR over time

output2$drug_cell <- paste(output2$agent, output2$cell_line, sep = "_") 

output2$GR50 <- ifelse(output2$GR50 == Inf | output2$GR50 == -Inf, NA, output2$GR50)

output2_matrix <- output2

output2_matrix <- output2_matrix[,c(1,2,9)]

output2_matrix <- spread(output2_matrix, key = "agent", value = "GR50")

rownames(output2_matrix) <- output2_matrix[,1]

output2_matrix[,1] <- NULL

output2_matrix <- output2_matrix[,colSums(is.na(output2_matrix))<nrow(output2_matrix)]

library(gplots)

(as.matrix(output2_matrix))

ggplot2::ggplot(output2 ,aes(x = agent, y = cell_line, fill = log(GR50))) + geom_tile()+
  scale_fill_continuous(na.value = "grey")+
  theme(axis.text.x= element_text(angle = 90))

