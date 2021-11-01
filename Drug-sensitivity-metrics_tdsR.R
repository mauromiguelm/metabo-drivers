## This code takes a vector of growth and calculate GI50, the initial parsing of the data is the same as in the GRmetrics packege,
#  so it will be easier to incorporate later.


load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\growth_curves.RData")

library(dplyr); library(tidyr); library(ggplot2)

# GR metrics, their example code only use the endpoint and time at treatment, but I will try all time points

# we need to select time point zero as the time we treated the samples (as per their requirements)

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

#FIXME data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time > 0)
#FIXME data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time > 0)

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

drug = "Asparaginase"

data_comb <- rbind(data_Grmetrics_Ctr, data_Grmetrics_Ttm)

#data_comb <-subset(data_comb, cell_line =="T47D" & agent %in% c("-", "Pemetrexed"))

data_comb <-subset(data_comb, agent %in% c("-", drug))

data_comb$cell_line <- as.character(data_comb$cell_line)

data_comb <- subset(data_comb, time <= 100)

#data_comb <- subset(data_comb, cell_line == "T47D" & agent == "Methotrexate")

devtools::load_all("C:\\Mauro_r_library\\tdsR")

tmp <- tdsR_fit(inputData = data_comb,
                groupingVariables = c("cell_line", "agent"),
                upperLimitThreshold = 1,
                upperLimit = 0.8,
                timeTreatment = 0, 
                smoothData = T,
                orderConc = T)

output_params <- tdsR_getOutput(inputData = tmp, metric = "parameters")

output <- tdsR_getOutput(inputData = tmp, metric = "tdsR")

output <- data.frame(output)

tmp <- strsplit(rownames(output), split = " ")

output$cell_line <- unlist(lapply(tmp, "[[", 1))

output$agent <- unlist(lapply(tmp, "[[", 2))

# e <- ggplot(output, aes(x=reorder(agent,as.numeric(X1) , median), y = as.numeric(X1)))
# 
# e + geom_violin(aes(fill = agent))+
#   theme(legend.position = "n", axis.text.x = element_text(angle = 90, size = 13.1))+geom_jitter(
#     aes(),
#     position = position_jitter(0.2),
#     size = 1.5)





#### plot cell line growth curves #####

###FIXME to show Methotrexate
###FIXME to show Omacetaxine

#load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\tdsR_output.Rdata")

lapply(unique(data_comb$agent), function(idx){
  
  output_tdsR = output
  
  test <- subset(data_comb, agent == idx | agent == "-")
  
  test$key <- paste(test$cell_line, test$agent)
  
  output_tdsR$key <- paste(output_tdsR$cell_line, output_tdsR$agent)
  
  output_tdsR <- output_tdsR[,c("tds", "key")]
  
  test <- left_join(test, output_tdsR, by = "key")
  
  test$tds <- as.numeric(as.character(test$tds))
  
  #test <- subset(test, cell_line == "ACHN" | cell_line == "M14" )
  
  library(RColorBrewer)
  
  p <- ggplot(test, aes(x = time, y= cell_count, col = factor(concentration )))+
    geom_smooth()+
    geom_vline(aes(xintercept = tds))+
    facet_wrap(~cell_line, ncol = 3)+
    scale_colour_hue(h = c(270, 360))+
    theme_bw()+
    scale_x_continuous(n.breaks = 10)+
    geom_vline(xintercept = 0, linetype="dotted", size = 1)
  
  setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\tdsR_orderConc_zero-lassymp")
  
  ggsave(filename = paste0(idx, ".png"), p, device = "png")
  
})

output_tdsR = output

test <- subset(data_comb, agent == "Gemci" | agent == "-")

test$key <- paste(test$cell_line, test$agent)

output_tdsR$key <- paste(output_tdsR$cell_line, output_tdsR$agent)

output_tdsR <- output_tdsR[,c("tds", "key")]

test <- left_join(test, output_tdsR, by = "key")

test$tds <- as.numeric(as.character(test$tds))

#test <- subset(test, cell_line == "ACHN" | cell_line == "M14" )

library(RColorBrewer)

ggplot(test, aes(x = time, y= cell_count, col = factor(concentration )))+
  geom_smooth()+
  geom_vline(aes(xintercept = tds))+
  facet_wrap(~cell_line, ncol = 3)+
  scale_colour_hue(h = c(270, 360))+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)+
  geom_vline(xintercept = 0, linetype="dotted", size = 1)
  #+

  #scale_x_continuous(breaks = seq(-10,100,10), limits = c(-10,100))

##### plot output parameters ####

tmp <- subset(output_params, key == "ACHN Erlotinib")


### export for downstream processing

rm(list = ls()[!ls()%in% "output_tdsR"])

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data")

save(output_tdsR, file = "tdsR_output.RData")

