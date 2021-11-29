## This code takes a vector of growth and calculate GI50, the initial parsing of the data is the same as in the GRmetrics packege,
#  so it will be easier to incorporate later.



# load packages nad definitions -------------------------------------------



path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\growth_metrics'

load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\growth_curves_filtered.RData")

library(dplyr); library(tidyr); library(ggplot2);library(gplots)library(RColorBrewer)
library(plyr)

devtools::load_all("C:\\Users\\masierom\\polybox\\Programing\\GRmetrics_GI50")
devtools::load_all("C:\\Users\\masierom\\polybox\\Programing\\tdsR")

time_treatment <- 0

drugs_in_screen <- c(unique(data_corrected$Drug[!(data_corrected$Drug %in% c("PBS", "DMSO"))]))

# calculate growth inhibition 50 metric -----------------------------------

#split analysis between plate 1 and plate 2

lapply(c("P1",'P2'),function(plate){
  data_GRmetrics <- subset(data_corrected, source_plate == plate)
  
  
  #cell_line <- "NCIH460"
  #drug  <- "Docetaxel"
  
  data_Grmetrics_Ttm <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # preparing the data for one drug
  data_Grmetrics_Ctr <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # getting the matching controls ready
  
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Drug %in% drugs_in_screen) #I chose clofarabine since it has a pretty dose response curve
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Drug == "DMSO" & Final_conc_uM %in% c(367)) # I separated the drug from the control, as its easier to parse
  
  # keep only matching time points
  
  match_time_intervals <- intersect(data_Grmetrics_Ctr$Time, data_Grmetrics_Ttm$Time) #make sure both datasets cover the same time
  
  time_treatment <- max(match_time_intervals[match_time_intervals < time_treatment])
  
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time %in% match_time_intervals)
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time %in% match_time_intervals)
  
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time >= time_treatment & Time == 72 | Time == 0)
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time >= time_treatment & Time == 72 | Time == 0)
  
  data_Grmetrics_Ctr$Time <- data_Grmetrics_Ctr %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time
  data_Grmetrics_Ttm$Time <- data_Grmetrics_Ttm %>% group_by(cell, Drug) %>% mutate(Time = ifelse(Time == min(Time), 0, Time)) %>% ungroup() %>% .$Time
  
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
  
  data_comb$cell_line <- as.character(data_comb$cell_line)

  output1 = GRmetrics::GRfit(inputData = data_comb, groupingVariables =
                               c("cell_line", 'agent'), case = "C")
  
  #GRdrawDRC(output1, points =F)
  
  
  # GRbox(output1, metric ='GR50', groupVariable = c('agent'),
  #       pointColor = c("cell_line"))
  
  
  output<-GRmetrics::GRgetMetrics(output1)
  
  return(output)
  
  
}) -> output_GI50


output_GI50 <- do.call(rbind, output_GI50)

#save output

setwd(path_data_file)

write.csv(output_GI50, "outcomes_growth_inhibition.R")


#plot results

tmp <- subset(data_corrected,Time == 24, select= c('Drug', 'cell', "GR24", "Final_conc_uM"))

tmp <- tmp%>% group_by(cell, Drug, Final_conc_uM) %>% summarise(GR24 = mean(GR24))

tmp <- tmp%>% group_by(cell, Drug) %>%
  dplyr::slice(which.max(Final_conc_uM))

tmp$Final_conc_uM <- NULL

tmp <- pivot_wider(tmp, values_from = GR24, names_from = "cell")

tmp_names <- tmp$Drug

tmp$Drug <- NULL

rownames(tmp) <- tmp_names

setwd(paste(path_fig, sep = "\\"))
png(paste("GR24",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PuOr"), margins=c(16,16),
          cexRow = 1.5,cexCol = 1.5,na.color = 'grey', scale = 'col')
dev.off()

tmp <- output_GI50[,c("cell_line", "agent","GR50")]

tmp$GR50 <- log10(tmp$GR50)

tmp$GR50 <- ifelse(is.infinite(tmp$GR50),NA,tmp$GR50)


tmp <- pivot_wider(tmp, values_from = GR50, names_from = "cell_line")

tmp <- tmp[rowSums(is.na(tmp)) != ncol(tmp), ]

tmp_names <- tmp$agent

tmp$agent <- NULL

rownames(tmp) <- tmp_names


setwd(paste(path_fig, sep = "\\"))
png(paste("GR50",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PiYG"), margins=c(16,16),
          cexRow = 2.0,cexCol = 2.5,na.color = 'grey', scale = 'row',Rowv = 'none',Colv = 'none',keysize=0.75)
dev.off()



# #calculate percent inhibition for GR24 ----------------------------------



lapply(c("P1", "P2"), function(plate){ 
  # calculate results by plate
  #plate = "P1"
  tmp <- subset(data_corrected, Time == 0 & source_plate == plate, select = c("GR24","Drug", "cell", "Final_conc_uM", "source_plate")) #keep one time point with GR24.. all time points have the same value..
  
  tmp <- tmp %>% 
    dplyr::group_by(Drug, cell,Final_conc_uM) %>% 
    dplyr::summarise(median_GR24 = median(GR24))
  
  tmp <- tmp%>% dplyr::group_by(cell,Drug)%>% dplyr::arrange(Final_conc_uM) %>% dplyr::mutate(Final_conc_uM = sequence(n()))
  
  control <- subset(tmp, Drug == "DMSO")
  
  drug <- subset(tmp, Drug %in% drugs_in_screen)
  
  tmp <- dplyr::right_join(drug, control[,c("cell", "median_GR24")], by = "cell")
  
  tmp$percent_change_GR <- (tmp$median_GR24.x/tmp$median_GR24.y) *100
  
  return(tmp)
})-> output_GR24

output_GR24 <- do.call(rbind,output_GR24)

# save output 

setwd(path_data_file)

write.csv(output_GI50, "outcomes_GR24.R")

#plot results

tmp <- output_GR24

lapply(unique(tmp$cell), function(cell_idx){
  #plot results by cell
  #cell_idx = "MDAMB231"
  cell_data <- subset(tmp, cell == cell_idx, select=c(Drug, Final_conc_uM, percent_change_GR))
  
  cell_data <- pivot_wider(cell_data, values_from = percent_change_GR, names_from = Final_conc_uM)
  colnames(cell_data) <- paste("conc=",colnames(cell_data))
  
  tmp_name <- cell_data$`conc= Drug`
  
  cell_data$`conc= Drug` <- NULL
  
  rownames(cell_data) <- tmp_name
  
  setwd(paste(path_fig, sep = "\\"))
  png(paste(cell_idx,"percent_change_GR24.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(cell_data), trace="none", key=T,col = (RColorBrewer::brewer.pal(n=11, "RdYlBu")), margins=c(16,16),
            cexRow = 1.5,cexCol = 2.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none')
  dev.off()
})

col_breaks <- seq(0,120,length.out = 10)

col_idxs <- (RColorBrewer::brewer.pal(n=9, "RdYlBu"))

lapply(unique(tmp$Drug), function(drug_idx){
  #plot results by drug
  drug_data <- subset(tmp, Drug == drug_idx, select=c(cell, Final_conc_uM, percent_change_GR))
  
  threshold_low <- ifelse(drug_data$percent_change_GR<0,"<0","")
  
  threshold_high <- ifelse(drug_data$percent_change_GR>120,">120","")
  text_data <- cbind(drug_data[,c("cell", "Final_conc_uM")], data.frame(threshold_high))
  
  text_data$threshold_high <- paste0(text_data$threshold_high, threshold_low)
  
  drug_data$percent_change_GR <- ifelse(drug_data$percent_change_GR<0,0,drug_data$percent_change_GR)
  
  drug_data <- pivot_wider(drug_data, values_from = percent_change_GR, names_from = Final_conc_uM)
  colnames(drug_data) <- paste("conc=",colnames(drug_data))
  
  tmp_name <- drug_data$`conc= cell`
  
  drug_data$`conc= cell` <- NULL
  
  rownames(drug_data) <- tmp_name
  
  
  text_data <- pivot_wider(text_data, values_from = threshold_high, names_from = Final_conc_uM)
  colnames(text_data) <- paste("conc=",colnames(text_data))
  
  tmp_name <- text_data$`conc= cell`
  
  text_data$`conc= cell` <- NULL
  
  rownames(text_data) <- tmp_name
  
  
  
  setwd(paste(path_fig, sep = "\\"))
  png(paste("drug=",drug_idx,"percent_change_GR24.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3)
  dev.off()
  
  png(paste("drug=",drug_idx,"percent_change_GR24_cluster.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3,Colv='none')
  dev.off()
  
  
})


# #calculate tdsR ---------------------------------------------------------

lapply(c("P1","P2"), function(x){
  
  data_GRmetrics <- subset(data_corrected, source_plate == x)
  
  
  data_Grmetrics_Ttm <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # preparing the data for one drug
  data_Grmetrics_Ctr <- data_GRmetrics#[grepl(cell_line, rownames(data_GRmetrics)),] # getting the matching controls ready
  
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Drug %in% drugs_in_screen) #I chose clofarabine since it has a pretty dose response curve
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Drug == "DMSO") # I separated the drug from the control, as its easier to parse
  
  # keep only matching time points
  
  match_time_intervals <- unique(intersect(data_Grmetrics_Ctr$Time, data_Grmetrics_Ttm$Time)) #make sure both datasets cover the same time
  
  data_Grmetrics_Ttm <- subset(data_Grmetrics_Ttm, Time %in% match_time_intervals)
  data_Grmetrics_Ctr <- subset(data_Grmetrics_Ctr, Time %in% match_time_intervals)
  
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
  
  
  #data_comb <-subset(data_comb, cell_line =="T47D" & agent %in% c("-", "Pemetrexed"))
  
  data_comb <-subset(data_comb, agent %in% c("-", drugs_in_screen))
  
  data_comb$cell_line <- as.character(data_comb$cell_line)
  
  data_comb <- subset(data_comb, time <= 100 & time >=0)
  
  tmp <- tdsR_fit(inputData = data_comb,
                  groupingVariables = c("cell_line", "agent"),
                  upperLimit = 0.8,
                  upperLimitThreshold = 1.0,
                  timeTreatment = 0, 
                  smoothData = T,
                  orderConc = T)
  
  
  return(tmp)
})-> output_tdsR


output_tdsR <- lapply(output_tdsR,function(x){tdsR_getOutput(x, metric = "tdsR")} )

output_tdsR <- do.call(rbind, output_tdsR)

#save output 

setwd(path_data_file)

write.csv(output_tdsR, "outcomes_tdsR.R")


#plot results

output <- output_tdsR

tmp <- strsplit(rownames(output_tdsR), split = " ")

output$cell_line <- unlist(lapply(tmp, "[[", 1))

output$agent <- unlist(lapply(tmp, "[[", 2))

tmp <- subset(output, select= c(cell_line, agent, tds))

tmp <- pivot_wider(tmp, values_from = tds, names_from = "cell_line")

tmp_names <- tmp$agent

tmp$agent <- NULL

rownames(tmp) <- tmp_names

tmp_colour <- ifelse(tmp >= 24, 1,0)


setwd(paste(path_fig, sep = "\\"))
png(paste("tdsR_24h",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp_colour), trace="none", key=T,col = RColorBrewer::brewer.pal(n=3, "PuOr"), margins=c(16,16),
          cexRow = 1.5,cexCol = 1.5,na.color = 'grey')
dev.off()

setwd(paste(path_fig, sep = "\\"))
png(paste("tdsR",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PuOr"), margins=c(16,16),
          cexRow = 1.5,cexCol = 1.5,na.color = 'grey')
dev.off()


tmp <- subset(data_corrected, cell == "MALME3M" & Drug == "Omacetaxine")

ggplot(tmp, aes(Time, Conf, col = factor(Final_conc_uM)))+
  geom_point()


# #TODO calculate corrected GI50 ------------------------------------------
#TODO calculate corrected GI50




