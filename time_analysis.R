if(!require(chron)){
  install.packages("chron")
  library(chron)
}

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen")

raw_data <- readxl::read_xlsx("Description.xlsx")




# Preparing time vectors for 384 well plates ------------------------------

# Creating a time vector

time_data <- matrix(data = NA, 
                    nrow = nrow(raw_data), 
                    ncol = 13, 
                    dimnames = list(NULL, 
                                    c("cell_line", "date_start","treatment_day", "sampling_day",
                                      "time_treatment384_p1_end", "time_treatment384_p2_end",
                                      "time_seed_96_finish", "time_treatment_96p1_end", "time_treatment_96p2_end",
                                      "p1_sampling_start", "p1_sampling_end",
                                      "p2_sampling_start", "p2_sampling_end")))

time_data <- as.data.frame(time_data)

time_data[,"cell_line"] <- raw_data$cell_name

time_data[,"date_start"] <- raw_data$start_data

time_data[,"date_start"] <- as.Date(as.character(time_data$date_start), format = "%Y%m%d")

time_data[,"sampling_day"] <- raw_data$day_sampling

time_data[,"sampling_day"] <- time_data[,"date_start"] + as.numeric(time_data[,"sampling_day"])

time_data[,"treatment_day"] <- time_data[,"date_start"]

time_data[,"treatment_day"] <- time_data[,"date_start"] + (raw_data$day_sampling-1)

convert_time <- function(x){
  if(is.na(x)){
    return(NA)
    }else{
      lim <- (2-nchar(x)%%2)
      txt <- paste(substring(x,1,lim),substring(x,lim+1,nchar(x)),"00", sep = ":")
      return(txt)
    }
  }

#organizing SEEDING time

time_data[,"time_seed_96_finish"] <- unlist(lapply(raw_data$`96_finish_seeding`, convert_time))
time_data[,"time_seed_96_finish"] <- chron(times = time_data[,"time_seed_96_finish"]) 

tmp <- paste(time_data[,"date_start"], time_data[,"time_seed_96_finish"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"time_seed_96_finish"] <- tmp 

#organizing 384 TREATMENT time

time_data[,"time_treatment384_p1_end"] <- unlist(lapply(raw_data$time_treatment_384_p1_end, convert_time))
time_data[,"time_treatment384_p1_end"] <- chron(times = time_data[,"time_treatment384_p1_end"]) 

tmp <- paste(time_data[,"treatment_day"], time_data[,"time_treatment384_p1_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"time_treatment384_p1_end"] <- tmp 


time_data[,"time_treatment384_p2_end"] <- unlist(lapply(raw_data$time_treatment_384_p2_end, convert_time))
time_data[,"time_treatment384_p2_end"] <- chron(times = time_data[,"time_treatment384_p2_end"]) 


tmp <- paste(time_data[,"treatment_day"], time_data[,"time_treatment384_p2_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"time_treatment384_p2_end"] <- tmp 

#organizing 96 TREATMENT time

time_data[,"time_treatment_96p1_end"] <- unlist(lapply(raw_data$`time_treatment_96_p1_end***`, convert_time))
time_data[,"time_treatment_96p1_end"] <- chron(times = time_data[,"time_treatment_96p1_end"]) 

tmp <- paste(time_data[,"treatment_day"], time_data[,"time_treatment_96p1_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"time_treatment_96p1_end"] <- tmp 

time_data[,"time_treatment_96p2_end"] <- unlist(lapply(raw_data$`time_treatment_96_p2_end***`, convert_time))
time_data[,"time_treatment_96p2_end"] <- chron(times = time_data[,"time_treatment_96p2_end"]) 

tmp <- paste(time_data[,"treatment_day"], time_data[,"time_treatment_96p2_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"time_treatment_96p2_end"] <- tmp 

#organizing SAMPLING TIME 96wp 

time_data[,"p1_sampling_start"] <- unlist(lapply(raw_data$`p1_sampling start`, convert_time))
time_data[,"p1_sampling_start"] <- chron(times = time_data[,"p1_sampling_start"]) 

tmp <- paste(time_data[,"sampling_day"], time_data[,"p1_sampling_start"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"p1_sampling_start"] <- tmp 

time_data[,"p1_sampling_end"] <- unlist(lapply(raw_data$`p1_sampling end`, convert_time))
time_data[,"p1_sampling_end"] <- chron(times = time_data[,"p1_sampling_end"]) 

tmp <- paste(time_data[,"sampling_day"], time_data[,"p1_sampling_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"p1_sampling_end"] <- tmp


time_data[,"p2_sampling_start"] <- unlist(lapply(raw_data$`p2_sampling start`, convert_time))
time_data[,"p2_sampling_start"] <- chron(times = time_data[,"p2_sampling_start"]) 

tmp <- paste(time_data[,"sampling_day"], time_data[,"p2_sampling_start"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"p2_sampling_start"] <- tmp

time_data[,"p2_sampling_end"] <- unlist(lapply(raw_data$`p2_sampling end`, convert_time))
time_data[,"p2_sampling_end"] <- chron(times = time_data[,"p2_sampling_end"]) 

tmp <- paste(time_data[,"sampling_day"], time_data[,"p2_sampling_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"p2_sampling_end"] <- tmp

# creating list vector with approximate sampling time

# here the order is q1, q2, q3 and q4 for p1 and p2, respectivelly.



lapply(rownames(time_data), function(row){
  
  #row <- 8
  
  data <- time_data[row,]
  
  if(any(is.na(c(data[["p1_sampling_end"]], data[["p1_sampling_start"]])))){
    
    return(NA)
    
    }else{
      
      sampling_length <- difftime(data[["p1_sampling_end"]], data[["p1_sampling_start"]], units = "secs")
      
      sampling_length <- seq(0, as.numeric(sampling_length), length.out = 4)
      
      plate_sampling_start <- lapply(sampling_length, function(time)
      {data["p1_sampling_start"] + time})
      
      plate_sampling_start <- do.call(rbind, plate_sampling_start)[,1]
      
      return(plate_sampling_start)

      
    }
    
    
  }) -> time_sampling_96p1

lapply(as.numeric(rownames(time_data)), function(row){
  
  #row <- 4

  data <- time_data[row,]
  
  time_sampling <- time_sampling_96p1[[row]]
  
  if(any(is.na(c(data[["time_treatment384_p1_end"]], time_sampling)))){
    
    return(NA)
    
  }else{
  
  time_sampling <- time_sampling_96p1[[row]]
  
  treatment_time <- difftime(time_sampling, data$time_treatment_96p1_end, units = "h")
  
  return(treatment_time)
  
  }
  
}) -> time_diff_96

time_diff_96 <- data.frame(sampling_p1= unlist(time_diff_96))

## plate 2

lapply(rownames(time_data), function(row){
  
  #row <- 8
  
  data <- time_data[row,]
  
  if(any(is.na(c(data[["p2_sampling_end"]], data[["p2_sampling_start"]])))){
    
    return(NA)
    
  }else{
    
    sampling_length <- difftime(data[["p2_sampling_end"]], data[["p2_sampling_start"]], units = "secs")
    
    sampling_length <- seq(0, as.numeric(sampling_length), length.out = 4)
    
    plate_sampling_start <- lapply(sampling_length, function(time)
    {data["p2_sampling_start"] + time})
    
    plate_sampling_start <- do.call(rbind, plate_sampling_start)[,1]
    
    return(plate_sampling_start)
    
  }
  
  
}) -> time_sampling_96p2

unlist(lapply(as.numeric(rownames(time_data)), function(row){
  
  #row <- 4
  
  data <- time_data[row,]
  
  time_sampling <- time_sampling_96p2[[row]]
  
  if(any(is.na(c(data[["time_treatment384_p2_end"]], time_sampling)))){
    
    return(NA)
    
  }else{
    
    time_sampling <- time_sampling_96p2[[row]]
    
    treatment_time <-difftime(time_sampling, data$time_treatment_96p2_end, units = "h")
    
    
    return(treatment_time)
    
  }
  
})) -> time_diff_96$sampling_p2

# plotting times 384 incubation times ----------------------------------------------------------

plot = F # plot = T

if(plot == T){
  
  require(ggplot2)
  require(ggExtra)
  
  diff_times_384 <- data.frame(incubation_time_384p1 = as.numeric(time_data$time_treatment384_p1_end - time_data$time_seed_96_finish))
  diff_times_384$incubation_time_384p2 <- as.numeric(time_data$time_treatment384_p2_end - time_data$time_seed_96_finish)
  
  (p <- ggplot(diff_times_384, aes(x=incubation_time_384p1, y=incubation_time_384p2)) +
      geom_point() +
      theme(legend.position="none")+
      scale_x_continuous(breaks = seq(floor(min(diff_times_384, na.rm = T)), ceiling(max(diff_times_384, na.rm = T)),2))+
      scale_y_continuous(breaks = seq(floor(min(diff_times_384, na.rm = T)), ceiling(max(diff_times_384, na.rm = T)),2))+
      theme_bw())
  
  # Set relative size of marginal plots (main plot 10x bigger than marginals)
  (p <- ggMarginal(p, type="histogram", size=2))
  
  ggsave(filename = "paperReady_384_incubation_times.png", plot = p, device = "png",
         path = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures")
  
}


# plotting times 96 treatment times --------------------------------------

if(plot == T){
  
  require(ggplot2)
  require(ggExtra)
  
  (p <- ggplot(time_diff_96, aes(x= sampling_p1, y=sampling_p2)) +
      geom_point() +
      theme(legend.position="none")+
      scale_x_continuous(breaks = seq(floor(min(diff_times_384, na.rm = T)), ceiling(max(diff_times_384, na.rm = T)),2))+
      scale_y_continuous(breaks = seq(floor(min(diff_times_384, na.rm = T)), ceiling(max(diff_times_384, na.rm = T)),2))+
      theme_bw())
  
  # Set relative size of marginal plots (main plot 10x bigger than marginals)
  (p <- ggMarginal(p, type="histogram", size=2))
  
  ggsave(filename = "paperReady_384_incubation_times.png", plot = p, device = "png",
         path = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\figures")
  
}

# organizing data for output ----------------------------------

time_vectors <- list(plate384 = time_data,
                     samp_p1_96 = time_sampling_96p1,
                     samp_p2_96 = time_sampling_96p2)

rm(list = ls()[ls() != "time_vectors"])
