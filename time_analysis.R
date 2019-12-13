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
time_data[,"time_treatment_96p1_end"] <- chron(times = time_data[,"time_treatment96_p1_end"]) 

tmp <- paste(time_data[,"treatment_day"], time_data[,"time_treatment_96p1_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"time_treatment_96p1_end"] <- tmp 

time_data[,"time_treatment_96p2_end"] <- unlist(lapply(raw_data$`time_treatment_96_p2_end***`, convert_time))
time_data[,"time_treatment_96p2_end"] <- chron(times = time_data[,"time_treatment96_p2_end"]) 

tmp <- paste(time_data[,"treatment_day"], time_data[,"time_treatment_96p2_end"])

tmp <- ifelse(grepl(patter = "NA", tmp), NA, tmp)

tmp <- as.POSIXct(tmp)

time_data[,"time_treatment_96p2_end"] <- tmp 

#organizing SAMPLING TIME time

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

