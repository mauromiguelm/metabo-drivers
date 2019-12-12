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
                    ncol = 3, 
                    dimnames = list(NULL, 
                                    c("cell_line", "date_start", "sampling_day")))

time_data <- as.data.frame(time_data)

time_data[,"cell_line"] <- raw_data$cell_name

time_data[,"date_start"] <- raw_data$start_data

time_data[,"date_start"] <- as.Date(as.character(time_data$date_start), format = "%Y%m%d")

time_data[,"sampling_day"] <- raw_data$day_sampling

time_data[,"sampling_day"] <- time_data[,"date_start"] + as.numeric(time_data[,"sampling_day"])

time_data <- 

chron::chron(times = raw_data$`p1_sampling start`, format = "hm")





