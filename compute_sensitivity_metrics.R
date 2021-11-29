## This code takes a vector of growth and calculate GI50, the initial parsing of the data is the same as in the GRmetrics packege,
#  so it will be easier to incorporate later.


load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\growth_curves.RData")

library(dplyr); library(tidyr); library(ggplot2)

library(GRmetrics, lib.loc = "C:\\Mauro_r_library\\tdsR")

# GR metrics, their example code only use the endpoint and time at treatment, but I will try all time points

# we need to select time point zero as the time we treated the samples (as per their requirements)

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



#GRdrawDRC(output1, points =F)


# GRbox(output1, metric ='GR50', groupVariable = c('agent'),
#       pointColor = c("cell_line"))


output<-GRmetrics::GRgetMetrics(output1)

library(corrplot)

test <- as.matrix(output[,c(7:14, 17:23, 26:31)])

test <- cor(test, method = "spe")

corrplot(test, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45)

## import tdsR data ####


library(dplyr)

load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\tdsR_output.Rdata")

### visualize tds across drugs ###

ggplot(output_tdsR , aes(x = reorder(agent, tds, median), y = tds, fill = agent))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15), axis.text.y = element_text(size =15))+
  scale_y_continuous(n.breaks = 10)

### tdsR growth rate: fit across all agents ###

output_tdsR$experiment <- rownames(output_tdsR)

output_tdsR[,c("cell_line", "agent")] <- NULL

colnames(output_tdsR)[2] <- "tds"

output_tdsR$tds <- as.numeric(as.character(output_tdsR$tds))

output_comb <- output %>% dplyr::inner_join(output_tdsR, by = "experiment")

test <- as.matrix(output_comb[,c("ctrl_cell_doublings", "GI50", "tds")])

test <- cor(test, method = "spe")

corrplot(test, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, addCoef.col = T)


### tdsR growth rate: fit across each drug to see drug specific effects ###

output_comb <- output %>% dplyr::inner_join(output_tdsR, by = "experiment")

test <- as.matrix(output_comb[,c("ctrl_cell_doublings", "GI50", "tds")])

test <- cor(test, method = "spe")

corrplot(test, type = "upper", order = "hclust",
         tl.col = "black", tl.srt = 45, addCoef.col = T)

### tdsR: is there a correlation between GR and tdsR ###

data_cor <- do.call( rbind, lapply( split(output_comb, output_comb$agent),
                                   function(x) data.frame(group = x$agent[1], cor = cor(x$ctrl_cell_doublings, x$tds, method = "spe")) ) )

data_cor$cor <- round(data_cor$cor,digits = 2)

ggplot(data_cor, aes(x = 1,y = reorder(group,as.numeric(cor) , mean), fill = cor)) +
  geom_tile(aes(fill = cor)) +
  geom_text(aes(label = cor)) +
  ylab("")+
  scale_fill_gradient2(low = "green", high = "red", mid = "white")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))



### tdsR: is there a correlation between GI50 and tdsR ###

data_cor <- do.call( rbind, lapply( split(output_comb, output_comb$agent),
                                    function(x) data.frame(group = x$agent[1], cor = cor(x$GI50, x$tds, method = "spe")) ) )

data_cor$cor <- round(data_cor$cor,digits = 2)

ggplot(data_cor, aes(x = 1,y = reorder(group,as.numeric(cor) , mean), fill = cor)) +
  geom_tile(aes(fill = cor)) +
  geom_text(aes(label = cor)) +
  ylab("")+
  scale_fill_gradient2(low = "green", high = "red", mid = "white")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures\\tdsR")

lapply(split(output_comb, output_comb$agent),function(x){

  fig <- ggplot(x, aes(log(GI50), tds, col = cell_line))+
    geom_point()+
    theme_bw()+
    theme(legend.title=element_text(size=2),
          legend.text=element_text(size=2),
            aspect.ratio=1)



    ggsave(paste0(x$agent[1],".png"),plot = fig  , device = "png")

       })


library(plotly)

tmp <- subset(output_comb, agent == "Pemetrexed")

fig <- ggplot(tmp, aes((GI50), tds, col = cell_line))+
  geom_point()+
  theme_bw()+
  theme(legend.title=element_text(size=5),
        legend.text=element_text(size=5),
        aspect.ratio=1)

fig

ggplotly(fig)


### check how many drugs could not be fissed

lapply(split(output_comb, output_comb$agent), function(x){

  #x = split(output_comb, output_comb$agent)[[2]]

  out <- sum(x$GR50 == Inf)

  names(out) <- x$agent[1]

  return(out)


}) -> unfit

unfit_df <- t(data.frame(unfit))

unfit_df <- as.data.frame(unfit_df)

colnames(unfit_df)[1] <- "lack of log fit"

unfit_df$agent <- names(unfit)



######### heatmap gi50 vs. tds  ######

metric = "tds"

tmp <- do.call(data.frame,lapply(output_comb, function(x) replace(x, is.infinite(x),NA)))

tmp <- tmp[,c("cell_line", metric, "agent")]

tmp <- spread(data = tmp, key = "agent", value = get(metric))

rownames(tmp) <- tmp$cell_line

tmp$cell_line <- NULL

colfunc<-colorRampPalette(c("green","white","red"))

gplots::heatmap.2(as.matrix(tmp), col= colfunc(200), na.color = "grey", trace="none", 
                  scale = "n", margins=c(12,8), srtCol=75, cexCol = 1.2)

rm(tmp)


#### calculate corrected GI50 ########

tmp <- output_comb

selected_drugs <- unfit_df$agent[unfit_df$`lack of log fit` <=11]

tmp <- subset(tmp, agent%in% selected_drugs) 
 
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
    
    print(idx$cell_line)
    print(idx$agent)
  
    
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
  
    
    output_corrected<-GRmetrics::GRgetMetrics(output1_corrected)
    
    idx$GI50_corr <- output_corrected$GI50
    
    return(idx)
    
  }else{
    
    idx$GI50_corr <- NA
    
    return(idx)
  }
  
}) -> tmp

test <- do.call(rbind, tmp)


plot(log(test$GI50_corr), log(test$GI50),col = test$cell_line,
     xlab = "GI50 Corrected", ylab = "GI50")

test$dist <- sqrt(((test$GI50)^2 + (test$GI50_corr)^2))

ggplot(test , aes(x = reorder(agent, dist, median), y = dist, fill = agent))+
  geom_boxplot()+
  geom_point(position = position_jitterdodge())+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15), axis.text.y = element_text(size =15))+
  scale_y_continuous(n.breaks = 10)

  ggplot(subset(test, agent == "Oxaliplatin"), aes(x = log(GI50), y = log(GI50_corr)))+
  geom_point(aes(col = cell_line))

data_cor <- do.call( rbind, lapply( split(test, as.character(test$agent)),
                                    function(x) data.frame(group = as.character(x$agent)[1], cor = cor(x$GI50, x$GI50_corr, method = "spe", use = "complete.obs"))))

data_cor$cor <- round(data_cor$cor, 2)

ggplot(data_cor, aes(x = 1,y = reorder(group,as.numeric(cor) , mean), fill = cor)) +
  geom_tile(aes(fill = cor)) +
  geom_text(aes(label = cor)) +
  ylab("")+
  scale_fill_gradient2(low = "green", high = "red", mid = "green")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


data_cor <- do.call( rbind, lapply( split(test, as.character(test$agent)),
                                    function(x) data.frame(group = as.character(x$agent)[1],
                                                           wilc_test = wilcox.test(x$GI50, x$GI50_corr, alternative ="two.sided")[[3]])))

ggplot(data_cor, aes(x = 1, y = reorder(group,as.numeric(wilc_test) , mean))) +
  geom_tile(aes(fill = as.numeric(wilc_test <0.05))) +
  geom_text(aes(label = wilc_test <0.05)) +
  ylab("")+
  scale_fill_gradient2(low = "blue", high = "red")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))


######### heatmap gi50 vs. gi50 new vs. log2fc(giNew/giOld)  ######

test$log2fc <- log2(test$GI50_corr/test$GI50)

metric = "log2fc"

tmp <- test[,c("cell_line", metric, "agent")]

tmp <- do.call(data.frame,lapply(tmp, function(x) replace(x, is.infinite(x),NA)))

tmp <- spread(data = tmp, key = "agent", value = get(metric))

rownames(tmp) <- tmp$cell_line

tmp$cell_line <- NULL

colfunc<-colorRampPalette(c("green","white","red"))

gplots::heatmap.2(as.matrix(tmp), col= colfunc(200), na.color = "grey", trace="none", 
                  scale = "n", margins=c(12,8), srtCol=75, cexCol = 1.2)

rm(tmp)
