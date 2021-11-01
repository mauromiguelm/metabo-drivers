library(ggplot2)
library(plyr)
library(dplyr)


##### import data ##########

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\from_Andrei\\distances\\cropped")

files <- list.files(full.names = T, include.dirs = T, recursive = T, pattern = ".json")

files <- files[grepl("euclidean", files)]

lapply(files,function(fileName){
  
  cell_line <- strsplit(strsplit(fileName, split = "/")[[1]][3], split = "_")[[1]][1]
  
  drug <- strsplit(strsplit(fileName, split = "/")[[1]][4], split = "\\.")[[1]][1]
  
  data <- rjson::fromJSON(file = fileName)
  
  time <- data$time
  
  data$time <- NULL
  
  data <- reshape2::melt(data, value.name = "value")
  
  data$time <- time
  
  data$cell_line <- cell_line
  
  data$agent <- drug
  
  return(data)
  
  
}) -> data

data <- do.call(rbind, data)

#####  statistics #########

# calculate average control before time zero, and use that for tdsR::pheno-ml
# plot CV of distance before zero controls, to see if there are outliers for one cell or drug

data_stats <- data %>% group_by(cell_line, agent) %>% dplyr::filter(time < 0) %>%
  dplyr::summarise(mean = mean(value), median = median(value), n = n(), sd = sd(value))

data_stats$groups <-  do.call(paste, data_stats[,c(1,2)])

ggplot(data_stats, aes(y = mean,x = sd, group = agent, col = agent))+
  geom_point()


######## normalize all data points by baseline control mean #######

data_norm <- data 

data_norm$norm_groups <- do.call(paste, data[,!colnames(data) %in%c("value", "time")])

norm_groups <- unique(data_norm$norm_groups)

lapply(norm_groups, function(idx){
  
  #idx = norm_groups[1]
  
  data_sub <- subset(data_norm, norm_groups == idx)
  
  control_avg <- mean(base::subset(data_sub, time <0, value, drop = T))
  
  data_sub$value <- data_sub$value - control_avg
  
  return(data_sub)
  
}) -> data_norm

data_norm <- do.call(rbind, data_norm)

data_norm$norm_groups <- NULL 



#### calculate tdsR #######

source("C:\\Users\\masierom\\Downloads\\tdsr-pheno-ml\\core.R")

#c("value", "well", "concentration", "agent","time")

data_tds <- data_norm

colnames(data_tds)[2] <- "concentration"

groups <- do.call(paste, data_tds[,c(4,5)])

data_tds$groups <- groups

groups <- unique(groups)

#groups <- unique(groups)[grepl(pattern = "ACHN", x = unique(groups))]

lapply(unique(groups), function(idx){
  
  #idx = groups[1]
  
  data_sub  <- subset(data_tds, groups == idx)
  
    tmp <- tdsR_fit(inputData = data_sub, groupingVariables = c("cell_line", "agent"), 
                  case = "B",
                  timeTreatment = 0,
                  lowerLimit = 1,   
                  smoothData = T,orderConc = T, 
                  limitThreshold = 5)
  
  
}) -> data_out

tmp <- lapply(data_out, "[[",2)

tmp <- as.data.frame(do.call(rbind, tmp))

names(tmp)[1] <- "tds"

tmp$groups <- groups

data_comb <- left_join(tmp, data_tds, by = "groups")

#### plot tds of morphology ####

data_comb$tds <- ifelse(data_comb$tds >= 100, 100, data_comb$tds)

ggplot(data_comb , aes(x = reorder(agent, tds, median), y = tds, fill = agent))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15), axis.text.y = element_text(size =15))+
  scale_y_continuous(n.breaks = 10)

ggplot(data_comb, aes(x = agent, y= tds, group = agent))+
  geom_boxplot()+
  geom_vline(aes(xintercept = tds))+
  facet_wrap(~agent, ncol = 3)+
  scale_colour_hue(h = c(270, 360))+
  theme_bw()+
  scale_x_continuous(n.breaks = 10)+
  geom_vline(xintercept = 0, linetype="dotted", size = 1)

##### compare tdsR with growth curves tdsR, to see if its different to any drug  #####

load("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data\\tdsR_output.Rdata")

output_tdsR$source <- "growth rates"

tmp$source <- "pheno-ml"

data_tds_comb <- rbind(output_tdsR[,c("groups", "tds", "source")], tmp[,c("groups", "tds", "source")])

data_tds_comb$cell_line <- unlist(lapply(strsplit(as.character(data_tds_comb$groups), " ", fixed = T), "[[", 1))

data_tds_comb$agent <- unlist(lapply(strsplit(as.character(data_tds_comb$groups), " ", fixed = T), "[[", 2))

tmp <- strsplit(data_tds_comb$agent, split = "_")

tmp <- unlist(lapply(tmp,"[[", 1))

data_tds_comb$agent <- tmp

p1 <- ggplot(data_tds_comb , aes(x = reorder(agent, tds, median), y = tds, fill = source, colour = source))+
  geom_boxplot()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 15), axis.text.y = element_text(size =15))+
  scale_y_continuous(n.breaks = 10)

setwd()
ggsave()