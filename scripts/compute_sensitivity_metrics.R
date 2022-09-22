## This code takes growth curve from large screen and calculate drug metrics, 

# load packages nad definitions -------------------------------------------

path_metabolomics_in <- '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\metabolomicsData_processed'
#path_data_file = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data'
path_data_file = "C:\\Users\\mauro\\Documents\\phd_results"
#path_fig = '\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\growth_metrics'
path_fig = "C:\\Users\\mauro\\Documents\\phd_results"
load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\growth_curves_filtered.RData")

library(dplyr); library(tidyr); library(ggplot2);library(gplots);library(RColorBrewer)
library(plyr);library(viridis)

devtools::load_all("C:\\Users\\masierom\\polybox\\Programing\\GRmetrics_GI50")
library(tdsR)

time_treatment <- 0
cutoff_GR_max_growth_effect <- 50 #pairs of drug & conc with abs change < this number will be excluded
cutoff_GR_conc_min_CV <- 10 #pairs of drug & conc with CV < this number will be excluded
drugs_in_screen <- c(unique(data_corrected$Drug[!(data_corrected$Drug %in% c("PBS", "DMSO"))]))

#import log2fc metabolomics data

setwd(paste0(path_data_file,"\\metabolomics","\\log2fc"))

metab_fcs <- lapply(list.files(pattern = "_P"),read.csv)

metab_fcs <- do.call(rbind,metab_fcs)

metab_fcs$X <- NULL
names(metab_fcs) <- c("cell_line","source_plate",'drug','concentration','ionIndex','log2fc','pvalue')

#import cell metadata from depmap 

setwd("C:\\Users\\mauro\\Documents\\phd_results\\metadata_cells")

cell_metadata <- read.csv("metadata_cells_depmap.csv")

#import cleaned metadata

setwd(paste(path_data_file,'metabolomics', sep = "//"))
metadata <- read.csv("metadata_clean.csv")


# calculate growth inhibition 50 metrics -----------------------------------

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

write.csv(output_GI50, "outcomes_growth_inhibition50.csv")

#plot results for GR24

setwd(path_data_file)

output_GI50 <- read.csv("outcomes_growth_inhibition50.csv")

tmp <- subset(data_corrected,Time == 24, select= c('Drug', 'cell', "GR24", "Final_conc_uM"))

tmp <- tmp%>% dplyr::group_by(cell, Drug, Final_conc_uM) %>% dplyr::summarise(GR24 = mean(GR24))

tmp <- tmp%>% group_by(cell, Drug) %>%
  dplyr::slice(which.max(Final_conc_uM))

tmp$Final_conc_uM <- NULL

tmp <- pivot_wider(tmp, values_from = GR24, names_from = "cell")

tmp_names <- tmp$Drug

tmp$Drug <- NULL

rownames(tmp) <- tmp_names

setwd(paste(path_fig,"GR24", sep = "\\"))
png(paste("GR24",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PuOr"), margins=c(16,16),
          cexRow = 1.5,cexCol = 1.5,na.color = 'grey', scale = 'col')
dev.off()

#plot results for GI50

tmp <- output_GI50[,c("cell_line", "agent","GR50")]

tmp$GR50 <- log10(tmp$GR50)

tmp <- pivot_wider(tmp, values_from = GR50, names_from = "cell_line")

tmp <- data.frame(tmp)

tmp_names <- tmp$agent

tmp$agent <- NULL

rownames(tmp) <- tmp_names

tmp <- t(apply(tmp, 1,function(x) ifelse(is.infinite(x),max(x[is.finite(x)]*1.1,na.rm = T),x)))

tmp <- ifelse(is.infinite(tmp),NA,tmp)

tmp <- tmp[rowSums(is.na(tmp)) != ncol(tmp), ]

setwd(paste(path_fig,"GR50", sep = "\\"))
png(paste("GR50",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PiYG"), margins=c(16,16),
          cexRow = 1.0,cexCol = 1.5,scale = 'row')
dev.off()


# #calculate percent inhibition for GR24 ----------------------------------

lapply(c("P1", "P2"), function(plate){ 
  # calculate results by plate
  #plate = "P2"
  
  tmp <- subset(data_corrected, Time == 0 & source_plate == plate, select = c("GR24","Drug", "cell", "Final_conc_uM", "source_plate")) #keep one time point with GR24.. all time points have the same value..
  
  tmp <- tmp %>% 
    dplyr::group_by(Drug, cell,Final_conc_uM) %>% 
    dplyr::summarise(median_GR24 = median(GR24))
  
  lapply(unique(tmp$Drug), function(drug_idx){
    tmp <- subset(tmp, Drug == drug_idx)
    tmp$Final_conc_uM <- tmp %>%  dplyr::group_indices(Final_conc_uM)
    return(tmp)
  }) -> tmp
  
  tmp <- do.call(rbind,tmp)
  
  control <- subset(tmp, Drug == "DMSO")
  
  drug <- subset(tmp, Drug %in% drugs_in_screen)
  
  tmp <- dplyr::right_join(drug, control[,c("cell", "median_GR24")], by = "cell")
  
  tmp$percent_change_GR <- (tmp$median_GR24.x/tmp$median_GR24.y) *100
  
  return(tmp)
})-> output_GR24

output_GR24 <- do.call(rbind,output_GR24)

# save output 

setwd(path_data_file)

write.csv(output_GR24, "outcomes_GR24.csv")

#plot results

tmp <- output_GR24

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
  
  setwd(paste(path_fig,'RSgroups', sep = "\\"))
  png(paste("drug=",drug_idx,"percent_change_GR24.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3)
  dev.off()
  
  png(paste("drug=",drug_idx,"percent_change_GR24_cluster.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3,Colv='none')
  dev.off()
  
  
})


# #calculate percent inhibition for confluence ----------------------------------

lapply(c("P1", "P2"), function(plate){ 
  # calculate results by plate
  #plate = "P1"
  tmp <- subset(data_corrected, Time >=0 & source_plate == plate, select = c("Conf",'Time',"Drug", "cell", "Final_conc_uM", "source_plate")) #keep one time point with GR24.. all time points have the same value..
  
  tmp <- tmp %>% 
    dplyr::group_by(Drug, cell,Final_conc_uM,Time) %>% 
    dplyr::summarise(median_conf = median(Conf))
  
  tmp <- tmp%>% dplyr::group_by(cell,Drug,Time)%>% dplyr::arrange(Final_conc_uM) %>% dplyr::mutate(Final_conc_uM = sequence(n()))
  
  tmp$cell_time <- paste(tmp$cell, tmp$Time)
  
  control <- subset(tmp, Drug == "DMSO")
  
  drug <- subset(tmp, Drug %in% drugs_in_screen)
  
  tmp <- dplyr::right_join(drug, control[,c("cell_time", "median_conf")], by = "cell_time")
  
  tmp$percent_change_GR <- (tmp$median_conf.x/tmp$median_conf.y) *100
  
  return(tmp)
})-> output_conf

output_conf <- do.call(rbind,output_conf)

# save output 

setwd(path_data_file)

write.csv(output_conf, "outcomes_conf_change_to_control.csv")

# plot phenotypic results -------------------------------------------------

#plot results for GR24 as heatmap

tmp <- subset(data_corrected,Time == 24, select= c('Drug', 'cell', "GR24", "Final_conc_uM"))

tmp <- tmp%>% dplyr::group_by(cell, Drug, Final_conc_uM) %>% dplyr::summarise(GR24 = mean(GR24))

tmp <- tmp%>% group_by(cell, Drug) %>%
  dplyr::slice(which.max(Final_conc_uM))

tmp$Final_conc_uM <- NULL

tmp <- pivot_wider(tmp, values_from = GR24, names_from = "cell")

tmp_names <- tmp$Drug

tmp$Drug <- NULL

rownames(tmp) <- tmp_names

setwd(paste(path_fig,"GR24", sep = "\\"))
png(paste("GR24",".png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PuOr"), margins=c(16,16),
          cexRow = 1.5,cexCol = 1.5,na.color = 'grey', scale = 'col')
dev.off()

#plot results for GI50 as heatmap 

tmp <- output_GI50[,c("cell_line", "agent","IC50")]

tmp$IC50 <- log10(tmp$IC50)

tmp <- pivot_wider(tmp, values_from = IC50, names_from = "cell_line")

tmp <- data.frame(tmp)

tmp_names <- tmp$agent

tmp$agent <- NULL

rownames(tmp) <- tmp_names

tmp <- t(apply(tmp, 1,function(x) ifelse(is.infinite(x),max(x[is.finite(x)]*1,na.rm = T),x)))

tmp <- ifelse(is.infinite(tmp),NA,tmp)

tmp <- tmp[rowSums(is.na(tmp)) != ncol(tmp), ]

labels_heatmap <- cell_metadata[order(match(cell_metadata$cell_line_display_name, colnames(tmp))),]

library(RColorBrewer)

labels_heatmap$colors <- as.character(factor(labels_heatmap$lineage_1, labels = RColorBrewer::brewer.pal(length(unique(labels_heatmap$lineage_1)), "Spectral")))

setwd(paste0(path_fig,'\\phenotype_figures'))

drugs_of_heatmap <- rownames(tmp)

pdf("IC50_phenotypes_72_legend.pdf")
heatmap_dendogram <- heatmap.2(as.matrix(tmp), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "RdBu"), margins=c(16,16),
          cexRow = 0.85,cexCol = 0.85,scale = 'row',ColSideColors = labels_heatmap$colors)

legend(x="bottomright", legend=c(labels_heatmap[!duplicated(labels_heatmap$lineage_1),]$lineage_1), 
       fill=labels_heatmap[!duplicated(labels_heatmap$lineage_1),]$colors)
dev.off()

#import GR24 data

setwd(path_data_file)

output_GR24 <- read.csv("outcomes_GR24.csv")

tmp <- output_GR24

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
  
  setwd(paste(path_fig,'RSgroups', sep = "\\"))
  png(paste("drug=",drug_idx,"percent_change_GR24.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3)
  dev.off()
  
  png(paste("drug=",drug_idx,"percent_change_GR24_cluster.png",sep="_"),height = 1200,width = 1200)
  heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
            cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3,Colv='none')
  dev.off()
  
  
})

# plot variance for each pair of cell and drug

GR24_var <- output_GR24 %>% dplyr::group_by(Drug,cell) %>% dplyr::summarize(variation = (sd(percent_change_GR)/mean(percent_change_GR))*100)

GR24_var <- pivot_wider(GR24_var,values_from = variation, names_from = Drug)
GR24_var <- data.frame(GR24_var)
rownames(GR24_var) <- GR24_var$cell
GR24_var$cell <- NULL

labels_heatmap <- cell_metadata[order(match(cell_metadata$cell_line_display_name, rownames(GR24_var))),]


labels_heatmap$colors <- as.character(factor(labels_heatmap$lineage_1, labels = RColorBrewer::brewer.pal(length(unique(labels_heatmap$lineage_1)), "Spectral")))
library(RColorBrewer)

m=as.matrix((GR24_var))
plt <- heatmap(m,scale = "row",col = bluered(100), RowSideColors = labels_heatmap$colors)
dev.off()

# plot basal growth rates at time of sampling, drug == PBS

GR_basal_24h <- subset(data_corrected, Time == 24 & Drug %in% c("PBS"))

GR_basal_24h$cell <- factor(GR_basal_24h$cell, levels = labels_heatmap$cell_line_display_name[rev(plt$rowInd)])

ggplot(GR_basal_24h, aes(cell, GR24))+
  geom_boxplot()+
  scale_y_continuous(limits = c(0,0.07))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))

# plot avg basal growth rates at time of sampling, drug == PBS

GR24_basal_avg <- data.frame(GR_basal_24h %>% group_by(cell) %>% dplyr::summarize(avg_basal_GR24 = median(GR24)))

rownames(GR24_basal_avg) <- GR24_basal_avg$cell

GR24_basal_avg$cell <- NULL

GR24_basal_avg$tmp <- GR24_basal_avg$avg_basal_GR24


GR24_basal_avg <- GR24_basal_avg[order(rownames(GR24_basal_avg)),] 

GR24_basal_avg <- GR24_basal_avg[rownames(GR24_basal_avg)[heatmap_dendogram$colInd],]


setwd(paste0(path_fig,'\\phenotype_figures'))

pdf("GR_median_24h.pdf")
heatmap.2(as.matrix(GR24_basal_avg), trace="none", key=T,col = RColorBrewer::brewer.pal(n=11, "PiYG"), margins=c(16,16),
          cexRow = 0.85,cexCol = 0.85,scale = 'none',Rowv = "null",Colv = "null")

dev.off()

# plot basal growth rates at time of sampling, drug == all controls

GR_basal_24h <- subset(data_corrected, Time == 24 & Drug %in% c("PBS", "DMSO"))

GR_basal_24h$cell <- factor(GR_basal_24h$cell, levels = labels_heatmap$cell_line_display_name[rev(plt$rowInd)])

ggplot(GR_basal_24h, aes(cell, GR24, colour = Drug))+
  geom_boxplot()+
  scale_y_continuous(limits = c(0,0.07))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0, hjust = 1))


# plot variance across concentratin for each pair of cell and drug

GR24_var <- output_GR24 %>% dplyr::group_by(Drug,Final_conc_uM) %>% dplyr::summarize(variation = (sd(percent_change_GR)/mean(percent_change_GR))*100)

GR24_var <- pivot_wider(GR24_var,values_from = variation, names_from = Final_conc_uM)
GR24_var <- data.frame(GR24_var)
rownames(GR24_var) <- GR24_var$Drug
GR24_var$Drug <- NULL

GR24_var <-GR24_var[rev(na.omit(drugs_of_heatmap[heatmap_dendogram$rowInd],match(rownames(GR24_var) ))),] 

library(viridis)

setwd(paste0(path_fig,'\\phenotype_figures'))

pdf("variance_across_concentration.pdf")
heatmap.2(abs(as.matrix(GR24_var)), trace="none", key=T,col = viridis(20), margins=c(16,16),
          cexRow = 0.65,cexCol = 0.65,scale = 'row',Rowv = "null",Colv = "null")
dev.off()


# define groups R/S based on GR24 ----------------------------

# Remove GR24 outliers
# remove pairs of drug_conc that are inefective across all CCLs or that are too strong
#low = drug_conc without effect across all CCLs

#import GR24 data

setwd(path_data_file)

output_GR24 <- read.csv("outcomes_GR24.csv")

#calculate CV
GR24_outliers_low <- output_GR24 %>% dplyr::group_by(Drug, Final_conc_uM) %>% dplyr::summarize(variation = (sd(percent_change_GR)/mean(percent_change_GR))*100)

#label low
GR24_outliers_low$outliers <- ifelse(GR24_outliers_low$variation<=  cutoff_GR_conc_min_CV,"low", NA)

#save output low

setwd(path_data_file)

write.csv(GR24_outliers_low, "GR24_outliers_low.csv")

#define high outliers

GR24_outliers_high <- output_GR24

GR24_outliers_high$outliers <- ifelse(GR24_outliers_high$percent_change_GR<=  cutoff_GR_max_growth_effect,"high", NA)

#save output

setwd(path_data_file)

write.csv(GR24_outliers_high, "GR24_outliers_high.csv")


# plot exclusions  

exclusions_min_growth <-GR24_outliers_low[GR24_outliers_low$variation<=  cutoff_GR_conc_min_CV,]

tmp <- exclusions_min_growth

tmp$outliers <- NULL

tmp$X <- NULL

tmp$variation <- log10(abs(tmp$variation))

tmp <- pivot_wider(tmp, values_from = variation, names_from = Final_conc_uM)

colnames(tmp) <- paste("conc=",colnames(tmp))

tmp_name <- tmp$`conc= Drug`

tmp$`conc= Drug` <- NULL

rownames(tmp) <- tmp_name

exclusions <- ifelse(tmp <= log10(cutoff_GR_conc_min_CV), "off","")

setwd(paste(path_fig,'RSgroups', sep = "\\"))
png(paste("percent_change_CV.png",sep="_"),height = 1200,width = 1200)
heatmap.2(as.matrix(tmp), trace="none", key=T,col = rev(RColorBrewer::brewer.pal(n=11, "RdYlBu")), margins=c(16,16),
          cexRow = 1.5,cexCol = 2.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',cellnote = exclusions,notecol="black",notecex=2)
dev.off()

# create R/S groups based on GR24 results

#remove concentrations with no effect

filtered_data <- output_GR24
low_outliers <- subset(GR24_outliers_low, outliers == 'low')
rows_to_exclude <- !paste(filtered_data$Drug,filtered_data$Final_conc_uM)%in% paste(low_outliers$Drug, low_outliers$Final_conc_uM)
filtered_data <- filtered_data[rows_to_exclude,]

library(parallel)
numWorkers <- detectCores()-1

setwd(path_data_file)
cl <-makeCluster(numWorkers, type="PSOCK",outfile = "tmp_err.txt")

# iterate over thresholds and calculate pvals/fdr
get_random_threshold <- function(min,max, n_, dist = 40){
  #min = min value from vector to be sampled
  #max = max value from vector to be sampled
  #n = number of cutoffs
  #dist = min distance between cutoffs
  sample_distribution <- runif(n = n_*100, min = min, max = max) #generate uniform distribution
  cutoff_distribution <- NULL
  while(length(cutoff_distribution) <n_){
    #iterate until getting n cutoffs with distance > dist
    cutoffs <- sample(sample_distribution,size = 2)
    cutoffs <- sort(cutoffs)
    if(abs(cutoffs[1]-cutoffs[2])>dist){
      cutoff_distribution <- append(cutoff_distribution,list(cutoffs))
    }
  }
  return(cutoff_distribution)
}

lapply(unique(metab_fcs$drug), function(drug_idx){
  tmp <- subset(metab_fcs, drug == drug_idx)
  tmp$conc <- tmp %>%  dplyr::group_indices(concentration)
  return(tmp)
}) -> metab_fcs

metab_fcs <- do.call(rbind, metab_fcs)

iterate_over_thresholds <-  function(drug_idx,gr24_data,data_metab,get_random_threshold_){
  #drug_idx = '17-AAG'
  sub_data <- subset(gr24_data, Drug == drug_idx)
  min_value <- 20 #min(sub_data$percent_change_GR)
  max_value <- 90 #max(sub_data$percent_change_GR)
  thresholds <- get_random_threshold_(min=min_value,max=max_value,n_=1000)
  sub_metab_data <- data_metab[data_metab$drug==drug_idx,]
  sub_metab_data$cell_conc <- paste(sub_metab_data$cell_line,sub_metab_data$conc)
  
  results_comb <- data.frame()
  for(idx in 1:length(thresholds)){
    #idx = 6
    #idx = 7
    sub_data$group <-NA
    sub_data$group <- ifelse(sub_data$percent_change_GR <=thresholds[[idx]][1],'S',sub_data$group)
    sub_data$group <- ifelse(sub_data$percent_change_GR >=thresholds[[idx]][2],'R',sub_data$group)
    sub_data$group <- ifelse(is.na(sub_data$group),'I',sub_data$group)
    sub_data$group <- ifelse(sub_data$percent_change_GR <=20,NA,sub_data$group)
    sub_data$group <- ifelse(sub_data$percent_change_GR >=90,NA,sub_data$group)
    sub_data$cell_conc <- paste(sub_data$cell,sub_data$Final_conc_uM)
    #calculate stats between R/S groups
    data_s <- subset(sub_data, group=='I')
    data_s <- sub_metab_data[sub_metab_data$cell_conc %in% unique(data_s$cell_conc),]
    data_r <- subset(sub_data, group=='R')
    data_r <- sub_metab_data[sub_metab_data$cell_conc %in% unique(data_r$cell_conc),]
    
    lapply(unique(data_s$ionIndex),function(ion_idx){
      #ion_idx = 2
      if(min(length(unique(data_r$cell_conc)), length(unique(data_s$cell_conc)))>=2){
        pval <- t.test(data_s[data_s$ionIndex==ion_idx,'log2fc'],data_r[data_r$ionIndex==ion_idx,'log2fc'])[[3]]
        log2fc <- abs(median(data_s[data_s$ionIndex==ion_idx,'log2fc'])-median(data_r[data_r$ionIndex==ion_idx,'log2fc']))
        
        return(data.frame(drug = drug_idx,
                          t1=thresholds[[idx]][1],
                          t2=thresholds[[idx]][2],
                          log2fc=log2fc,
                          pval=pval))
      }
    }) ->results 
    
    results <- do.call(rbind, results)
    results_comb <- rbind(results_comb,results)
  }
  return(results_comb)

}

#iterate_over_thresholds(drug_idx,gr24_data,metadata_,data_metab,get_random_threshold_)  

parLapply(cl=cl,unique(filtered_data$Drug),iterate_over_thresholds,gr24_data = filtered_data,
                                                                data_metab = metab_fcs,
          get_random_threshold_ = get_random_threshold) -> out_thresholds


parallel::stopCluster(cl)
 
rm(cl)

# plotting threshold for each drug 

lapply(seq_along(out_thresholds),function(idx){
  
  #idx = 9
  sub_out <- out_thresholds[[idx]]
  drug <- unique(sub_out$drug)
  sub_out$log2fc <- abs(sub_out$log2fc)
  
  sub_out$t1t2 <- paste(sub_out$t1,sub_out$t2)
  
  #remove not significant ions
  sub_out <- subset(sub_out, pval<0.05)
  
  sub_out <- sub_out %>%dplyr::group_by(t1,t2) %>% dplyr::summarise(count_fc = sum(log2fc>log2(3/2)),
                                                             count_p = sum(pval<0.05),
                                                             mean_fc = mean(log2fc))
  
  setwd(paste0(path_fig,'\\iter_threshold'))
  
  ggplot(sub_out, aes(x=t1,y=t2,col = count_fc))+
    geom_point()+
    geom_point(data=sub_out[which(sub_out$count_fc == max(sub_out$count_fc)), ], colour="red", size=3)+
    scale_colour_viridis()->plt
  
  ggsave(paste0("ithreshold_drug=",drug,'_log2fc_pval_filter_RvsI_min20max80.png'),plt)
  
})

#determine the number of cell_conc_groups per threshold

lapply(seq_along(out_thresholds),function(idx){
  
  sub_out <- out_thresholds[[idx]]
  drug <- unique(sub_out$drug)
  sub_out$log2fc <- abs(sub_out$log2fc)
  
  sub_out$t1t2 <- paste(sub_out$t1,sub_out$t2)
  
  #remove not significant ions
  sub_out <- subset(sub_out, pval<0.05)
  
  sub_out <- sub_out %>%dplyr::group_by(t1,t2) %>% dplyr::summarise(count_fc = sum(log2fc>log2(3/2)),
                                                                    count_p = sum(pval<0.05),
                                                                    mean_fc = mean(log2fc))
  
    
  threshold <- sub_out[which(sub_out$count_fc == max(sub_out$count_fc)), ][1,]
  
  gr24_data <- subset(filtered_data,Drug == drug )
  
  sens <- sum(gr24_data$percent_change_GR <= threshold$t1)
  
  res <- sum(gr24_data$percent_change_GR >= threshold$t2)
  
  
  int <- sum(gr24_data$percent_change_GR < threshold$t2 &
        gr24_data$percent_change_GR > threshold$t1)
  
  out <- data.frame(drug,threshold$t1,threshold$t2,threshold$count_fc,sens,int,res)
  
  return(out)
  
}) -> counted_groups

  
counted_groups <- do.call(rbind,counted_groups)

# save results
setwd(path_data_file)

write.csv(file = 'counted_groups_ithreshold_GR24_RvsI_between20-90.csv',x = counted_groups)

save(list="out_thresholds",file='iter_threshold_GR24_RvsI_between20-90_diff40.Rdata')

#save filtered RS groups

setwd(path_data_file)

write.csv(filtered_data,'outcomes_GR24_RSgroups_filtered_RvsI.csv')

#plot R/S groups

tmp <- filtered_data

col_breaks <- seq(0,120,length.out = 10)

col_idxs <- (RColorBrewer::brewer.pal(n=9, "RdYlBu"))

lapply(unique(tmp$Drug), function(drug_idx){
  #plot results by drug
  
  #drug_idx = 'Methotrexate'
  print(drug_idx)
  drug_data <- subset(tmp, Drug == drug_idx, select=c(cell, Final_conc_uM, percent_change_GR))
  
  if(!all(is.na(drug_data$percent_change_GR))){
    
    threshold_high <- ifelse(drug_data$percent_change_GR>120,">120","")
    text_data <- cbind(drug_data[,c("cell", "Final_conc_uM")], data.frame(threshold_high))
    
    drug_data$percent_change_GR <- ifelse(drug_data$percent_change_GR<0,0,drug_data$percent_change_GR)
    
    drug_data <- pivot_wider(drug_data, values_from = percent_change_GR, names_from = Final_conc_uM)
    
    colnames(drug_data) <- paste("conc=",colnames(drug_data))
    
    drug_data <- drug_data[,colSums(is.na(drug_data))<nrow(drug_data)]
    
    tmp_name <- drug_data$`conc= cell`
    
    drug_data$`conc= cell` <- NULL
    
    rownames(drug_data) <- tmp_name
    
    if(ncol(drug_data)==1){
      # create a column of ones so heatmap can be generated
      
      drug_data$emptycol<- rep(1,length.out = nrow(drug_data))
      
      drug_data <- drug_data[rev(colnames(drug_data))]
      
      text_data$emoty_column <- ""
      
      text_data <- text_data[rev(colnames(text_data))]
    }
    
    rownames(drug_data) <- tmp_name
    
    text_data <- pivot_wider(text_data, values_from = threshold_high, names_from = Final_conc_uM)
    colnames(text_data) <- paste("conc=",colnames(text_data))
    text_data <- text_data[,colSums(is.na(text_data))<nrow(text_data)]
    
    tmp_name <- text_data$`conc= cell`
    
    text_data$`conc= cell` <- NULL
    
    rownames(text_data) <- tmp_name
    
    setwd(paste(path_fig,'RSgroups', sep = "\\"))
    png(paste("drug=",drug_idx,"percent_change_GR24_filtered.png",sep="_"),height = 1200,width = 1200)
    
    heatmap.2(as.matrix(drug_data), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
              cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = text_data,notecol="black",notecex=3)
    dev.off()
  
  }

})

#plot groups on filtered data

col_breaks <- seq(0,120,length.out = 10)

col_idxs <- (RColorBrewer::brewer.pal(n=9, "RdYlBu"))

lapply(unique(GR24_RSgroups$Drug), function(drug_idx){
  #plot results by drug on data filtered by too strong effects
  
  #drug_idx = "Clofarabine"
  print(drug_idx)
  drug_data <- subset(tmp, Drug == drug_idx, select=c(cell, Final_conc_uM, percent_change_GR))
  
  if(!all(is.na(drug_data$percent_change_GR))){
    
    text_data <- subset(GR24_RSgroups, Drug == drug_idx, select=c("cell", "Final_conc_uM", 'group'))
    
    drug_data$percent_change_GR <- ifelse(drug_data$percent_change_GR<0,0,drug_data$percent_change_GR)
    
    drug_data <- pivot_wider(drug_data, values_from = percent_change_GR, names_from = Final_conc_uM)
    
    colnames(drug_data) <- paste("conc=",colnames(drug_data))
    
    #drug_data <- drug_data[,colSums(is.na(drug_data))<nrow(drug_data)]   
    
    tmp_name <- drug_data$`conc= cell`
    
    drug_data$`conc= cell` <- NULL
    
    
    text_data <- pivot_wider(text_data, values_from = group, names_from = Final_conc_uM)
    colnames(text_data) <- paste("conc=",colnames(text_data))
    #text_data <- text_data[,colSums(is.na(text_data))<nrow(text_data)]
    
    tmp_name <- text_data$`conc= cell`
    
    text_data$`conc= cell` <- NULL
    
    rownames(text_data) <- tmp_name
    
    grouping_def <-apply(text_data,1, function(x){unique((x))})
    
    grouping_def <- unlist(lapply(grouping_def, function(x){if(length(x)>=2){na.omit(x)}else(x)}))
    
    R_idxs <- which(grouping_def == "R")
    S_idxs <- which(grouping_def =="S")
    
    R_group <- data.frame(drug_data[R_idxs,])
    rownames(R_group) <- tmp_name[R_idxs]
    R_text <- text_data[R_idxs,]
    S_group <- data.frame(drug_data[S_idxs,])
    rownames(S_group) <- tmp_name[S_idxs]
    S_text <- text_data[S_idxs,]
    
    
    if(nrow(S_group)==1){
      # create a column of ones so heatmap can be generated
      S_group <- rbind(S_group,S_group)
      S_group[1,] = rep(1, length.out = ncol(S_group)) 
    }
    
    if(nrow(R_group)==1){
      # create a column of ones so heatmap can be generated
      R_group <- rbind(R_group,R_group)
      R_group[1,] = rep(1, length.out = ncol(R_group)) 
    }
    
    
    rownames(drug_data) <- tmp_name
    
    if(all(c("R", "S") %in% unique(as.character(as.matrix(text_data))))){
      
      #plot Resistant
      
      setwd(paste(path_fig,"RSgroups", sep = "\\"))
      png(paste("drug=",drug_idx,"percent_change_GR24_groupCV_filtered_SEN.png",sep="_"),height = 1200,width = 1200)
      
      heatmap.2(as.matrix(R_group), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
                cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = R_text,notecol="black",notecex=3)
      dev.off()
      
      #plot Sensitive
      
  
      png(paste("drug=",drug_idx,"percent_change_GR24_groupCV_filtered_RES.png",sep="_"),height = 1200,width = 1200)
      
      heatmap.2(as.matrix(S_group), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
                cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = S_text,notecol="black",notecex=3)
      dev.off()
      
    }else{
      print(paste("drug=",drug_idx,"has not one of the groups"))
    }
    
  }
})


#plot results by drug on data not filtered by too strong effects

GR24_RSgroups_unfiltered <- subset(output_GR24, Drug %in% unique(GR24_RSgroups$Drug))

GR24_RSgroups_unfiltered$idx <- paste(GR24_RSgroups_unfiltered$Drug, GR24_RSgroups_unfiltered$cell)

tmp <- subset(GR24_RSgroups, select= c(cell, Drug, group))

tmp$idx <- paste(tmp$Drug, tmp$cell)

tmp <- tmp%>% dplyr::group_by(idx) %>% dplyr::slice(1)

GR24_RSgroups_unfiltered <- dplyr::left_join(GR24_RSgroups_unfiltered, tmp[,c("idx", "group")], by = "idx",keep=F)

col_breaks <- seq(0,120,length.out = 10)

col_idxs <- (RColorBrewer::brewer.pal(n=9, "RdYlBu"))

lapply(unique(GR24_RSgroups_unfiltered$Drug), function(drug_idx){
  #plot results by drug on data not filtered by too strong effects
  
  #drug_idx = 'Everolimus'
  print(drug_idx)
  drug_data <- subset(GR24_RSgroups_unfiltered, Drug == drug_idx, select=c(cell, Final_conc_uM, percent_change_GR))
  
  if(!all(is.na(drug_data$percent_change_GR))){
    
    text_data <- subset(GR24_RSgroups_unfiltered, Drug == drug_idx, select=c("cell", "Final_conc_uM", 'group'))
    
    drug_data$percent_change_GR <- ifelse(drug_data$percent_change_GR<0,0,drug_data$percent_change_GR)
    
    drug_data <- pivot_wider(drug_data, values_from = percent_change_GR, names_from = Final_conc_uM)
    
    colnames(drug_data) <- paste("conc=",colnames(drug_data))
    
    #drug_data <- drug_data[,colSums(is.na(drug_data))<nrow(drug_data)]   
    
    tmp_name <- drug_data$`conc= cell`
    
    drug_data$`conc= cell` <- NULL
    
    rownames(drug_data) <- tmp_name
    
    if(ncol(drug_data)==1){
      # create a column of ones so heatmap can be generated
      
      drug_data$emptycol<- rep(1,length.out = nrow(drug_data))
      
      drug_data <- drug_data[rev(colnames(drug_data))]
      
      text_data$empty_column <- ""
      
      text_data <- text_data[rev(colnames(text_data))]
    }
    
    
    rownames(drug_data) <- tmp_name
    
    text_data <- pivot_wider(text_data, values_from = group, names_from = Final_conc_uM)
    colnames(text_data) <- paste("conc=",colnames(text_data))
    #text_data <- text_data[,colSums(is.na(text_data))<nrow(text_data)]
    
    tmp_name <- text_data$`conc= cell`
    
    text_data$`conc= cell` <- NULL
    
    rownames(text_data) <- tmp_name
    
    grouping_def <-apply(text_data,1, function(x){unique((x))})
    
    grouping_def <- unlist(lapply(grouping_def, function(x){if(length(x)>=2){na.omit(x)}else(x)}))
    
    R_idxs <- which(grouping_def == "R")
    S_idxs <- which(grouping_def =="S")
    
    R_group <- data.frame(drug_data[R_idxs,])
    rownames(R_group) <- tmp_name[R_idxs]
    R_text <- text_data[R_idxs,]
    S_group <- data.frame(drug_data[S_idxs,])
    rownames(S_group) <- tmp_name[S_idxs]
    S_text <- text_data[S_idxs,]
    
    if(nrow(S_group)==1){
      # create a column of ones so heatmap can be generated
      S_group <- rbind(S_group,S_group)
      S_group[1,] = rep(1, length.out = ncol(S_group)) 
    }
    
    if(nrow(R_group)==1){
      # create a column of ones so heatmap can be generated
      R_group <- rbind(R_group,R_group)
      R_group[1,] = rep(1, length.out = ncol(R_group)) 
    }
    
    
    rownames(drug_data) <- tmp_name
    
    if(all(c("R", "S") %in% unique(as.character(as.matrix(text_data))))){
      #plot Resistant
      
      setwd(paste(path_fig,"RSgroups", sep = "\\"))
      
      png(paste("drug=",drug_idx,"percent_change_GR24_groupCV_RES.png",sep="_"),height = 1200,width = 1200)
      
      heatmap.2(as.matrix(R_group), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
                cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = R_text,notecol="black",notecex=3)
      dev.off()
      
      #plot Sensitive
      
    
      png(paste("drug=",drug_idx,"percent_change_GR24_groupCV_SEN.png",sep="_"),height = 1200,width = 1200)
      
      heatmap.2(as.matrix(S_group), trace="none", key=T,col = col_idxs, margins=c(16,16),breaks = col_breaks,
                cexRow = 2,cexCol = 3.5,na.color = 'grey', scale = 'none',Rowv = 'none',Colv = 'none',keysize=0.75,cellnote = S_text,notecol="black",notecex=3)
      dev.off()
      
    }else{
      print(paste("drug=",drug_idx,"has not one of the groups"))
    }
      
    }
})


# plot GR50 for each RS group


nonEdgy_trans <- function(){   scales::trans_new("nonEdgy",                     transform = function(x) {ifelse(is.finite(x), x, sign(x) * 6)},                     inverse = function(x) {ifelse(abs(x) == 6, sign(x) * Inf, x)}                     ) }

tmp <- subset(filtered_data, !is.na(percent_change_GR))
tmp <- subset(tmp, select= c(cell, Drug, group))

tmp$experiment <- paste(tmp$cell,tmp$Drug)

tmp <- tmp%>% dplyr::group_by(experiment) %>% dplyr::slice(1)

tmp <- output_GI50%>%dplyr::inner_join(tmp,by='experiment' )

lapply(unique(tmp$agent), function(drug_idx){
  
  drug_data <- subset(tmp, Drug == drug_idx)
  
  drug_data$group <- factor(drug_data$group, levels=c("R", "I", "S"))
  
  setwd(paste(path_fig,"RSgroups", sep = "\\"))
  
  plt <- ggplot(drug_data,aes(group, log10(GR50)))+
    geom_jitter(aes(colour = ifelse(!is.finite(GR50), 
                      "Infinite point", 'defined metrix')),
                width = 0.18)+
    scale_y_continuous(trans = "nonEdgy")+
    theme_light()+
    labs(colour="GR50 undefined")
  
  ggsave(paste("drug=",drug_idx,"grouping_by_GR50.png",sep="_"),plt,width = 6,height = 3)
  
  
})