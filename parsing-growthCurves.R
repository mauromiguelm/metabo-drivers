## import R functions ####

library(ggplot2); library(RColorBrewer)

# Importing plate times ####

source("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\analysis\\time_analysis.R")

# import functions cell culture #####

source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\CellCultureAnalysis.R")

source("C:\\Users\\masierom\\polybox\\Programing\\Project_exometabolites\\modelling_growth_curves.R")

#import plate 96-384 converter #####

source('C:/Users/masierom/polybox/Programing/96_to_384/Convert_96_to_384.R')

# Importing exceptions ####

source("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\exceptions\\log_processing.r")

## Importing source plates #####

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\growthData")

source_plates <- data.frame(
  filenames = list.files(
    path = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\Growthdata\\",
    pattern = "[P][1-2][.txt]",
    recursive = T))

# source_plates <- rbind(
#   data.frame(filenames = list.files ("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\growthData\\20190828\\results", 
#                                                    pattern = ".txt",
#                                                    recursive = T)))

source_plates$sourceid <- ifelse(grepl(pattern = "P1", source_plates$filenames), "1MSP001", "2MSP001")

source_plates$batch <- ifelse(grepl(pattern = "20190513", source_plates$filenames), "batch_1", "batch_2")

#import growth curves ####

fileNames <- as.character(source_plates$filenames)

fileNames<- fileNames[7:length(fileNames)] #remove CLs that do not have plate 2

data <- 
  
  lapply(fileNames, function(x) {
    
    #x = fileNames[10]
    
    df <-
      ReadIncuCyteData(read_platemap = F, FileName_IncuCyte = x,
                       Plate_size = 384, FileDirectory = getwd(),
                       time_output = "GMT")
    
    
    return(df)
    
  })

names(data) <- fileNames

#### load data ####

#load("C:/Users/mauro/Desktop/20200113_envir.2.RData")

### Create elapsed and sampling times based on time_vectors for the 384 plates ######

data <- lapply(names(data), function(filename){
  
  #filename <- names(data)[13]
  
  # filename = "20191029/results/OVCAR4_CL3_P1.txt"
  
  print(filename)
  
  start_str <- "results/"
  
  start_length <- nchar(start_str)
  
  start_pos <- regexpr(start_str, text = filename)[[1]][1] + start_length
  
  end_str <- ".txt" 
  
  end_pos <- regexpr(end_str, text = filename)[[1]][1] - 1
  
  cell_plate <- substr(filename, start_pos, end_pos)   
  
  cell_line <- strsplit(cell_plate, "_")[[1]][1]  
  
  plate <- strsplit(cell_plate, "_")[[1]][3]
  
  p1_treatment_end <-  time_vectors$plate384[time_vectors$plate384$cell_line == cell_line,
                                             "time_treatment_96p1_end"]
  
  p2_treatment_end <- time_vectors$plate384[time_vectors$plate384$cell_line == cell_line,
                                            "time_treatment_96p2_end"]
  
  
  if(plate == "P1"){
    
    data_elapsed <- data[[filename]]
    
    # for each plate, calculate the time in which it was sampled
    
    data_elapsed <- conv_384_96(data_elapsed) #get which 384 wells belong to each 96 quadrant
    
    #cell_line = "HCT15"
    
    quad_time_keys <-  time_vectors$samp_p1_96[[cell_line]]
    
    if(is.null(quad_time_keys)) stop("cannot find matching time for CL")
    
    samp_time <- quad_time_keys[data_elapsed$quart]
    
    samp_time_diff <- difftime(samp_time, 
                               as.POSIXct(as.character(p1_treatment_end), tz = "GMT"), units = "h")
    
    samp_time_diff <- round(as.numeric(samp_time_diff),2)
    
    data_elapsed$samp_time <- samp_time_diff
    
    # for each plate, create elapse time 
    
    data_elapsed$Time <- difftime(as.POSIXct(data_elapsed$Time, tz = "GMT"), 
                                  as.POSIXct(as.character(p1_treatment_end), tz = "GMT"), units = "h")
    
    data_elapsed$Time <- round(data_elapsed$Time,digits = 2)
    
    data_elapsed$Time <- as.numeric(data_elapsed$Time)
    
  }else if(plate == "P2"){
    
    data_elapsed <- data[[filename]]
    
    # for each plate, calculate the time in which it was sampled
    
    data_elapsed <- conv_384_96(data_elapsed) #get which 384 wells belong to each 96 quadrant
    
    #cell_line = "HCT15"
    
    quad_time_keys <-  time_vectors$samp_p2_96[[cell_line]]
    
    if(is.null(quad_time_keys)) stop("cannot find matching time for CL")
    
    samp_time <- quad_time_keys[data_elapsed$quart]
    
    samp_time_diff <- difftime(samp_time, 
                               as.POSIXct(as.character(p2_treatment_end), tz = "GMT"), units = "h")
    
    samp_time_diff <- round(as.numeric(samp_time_diff),2)
    
    data_elapsed$samp_time <- samp_time_diff
    
    data_elapsed$Time <- difftime(as.POSIXct(data_elapsed$Time, tz = "GMT"), 
                                  as.POSIXct(as.character(p2_treatment_end), tz = "GMT"), units = "h")
    
    data_elapsed$Time <- round(data_elapsed$Time,digits = 2)
    
    data_elapsed$Time <- as.numeric(data_elapsed$Time)
    
  }else(
    stop("problems with regex while parsing plate times")
  )
  
  return(data_elapsed)
  
})

names(data) <- fileNames

# removing bad wells based on growth before treatment #####

data_corrected <- data

skip_outlier <- c("20190513/results/HCT15_P1.txt", # For the HCT-15_P1, due to time points before drugs, we have a problem with outlier detection. Hence, I am using the unfilteed data
                  "20190828/results/ACHN_CL3_P1.txt", # only 24h of growth
                  "20190828/results/ACHN_CL3_P2.txt") # For the HCT-15_P1, due to time points before drugs, we have a problem with outlier detection. Hence, I am using the unfilteed data


data_skip <- data_corrected[skip_outlier]

data_corrected <- data_corrected[!(names(data_corrected) %in% skip_outlier)]

tmp_names <- names(data_corrected)

data_corrected <-
  
  lapply(names(data_corrected), FUN = function(x){
    
    print(x)
    
    # plate_idx = names(data_corrected)[14]   #FIXME the name contains / which is not allowed as a finame.. use regex to only get the real plate name or remove the backslashed
    
    #x = "20190925/results/SF539_CL1_P2.txt"
    
    plate_name = strsplit(x, split = "/")[[1]][3]
    
    x = data_corrected[[x]]
    
    good_wells <-
      filter_growth_outliers(plate_name = plate_name, data = x, time_control = 0, save_diag_plots = F,
                             save_plots_directory = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users/Mauro/cell_culture_data/190310_LargeScreen/figures")
    
    new_data <-
      x[!(x[,"Well"] %in% good_wells$wells_exception),]
    
    print(paste("Plate:", plate_name,
                "Exception-Number:",good_wells$wells_exception_number
    ))
    
    return(new_data)
    
  })


names(data_corrected) <- tmp_names

lapply(data_skip, function(plate){
  if(is.null(plate)){warning("the plate is null")}
})

data_corrected <- append(data_corrected, data_skip)

data_corrected <- data_corrected[fileNames]

rm(data_skip, tmp_names)

# Fit the confluence for each well, and return fitted confluence. Keep the same plate map structure. ####


#FIXME the functin below has to retunr all columns.. some are mising

# tmp <- data_corrected
# 
# data_corrected <- tmp

base::lapply(names(data_corrected), function(plate_name){
  
  #plate_name =  names(data_corrected)[5] #FIXME delete me
  
  print(plate_name)
  
  r.2.threshold = 0.8
  
  data_raw <- data_corrected[[plate_name]]
  
  plate_name <- strsplit(plate_name, split = "/")[[1]][3]
  
  grouping_vars <- "Well"
  
  metadata_cols <- colnames(data_raw)[!colnames(data_raw) %in% c("Time", "Conf")]
  
  metadata_df <- data_raw[,metadata_cols]
  
  metadata_df <- metadata_df %>% group_by(get(grouping_vars)) %>% slice(1) %>% ungroup()
  
  metadata_df[ncol(metadata_df)] <- NULL
  
  fitted_data <- 
    
    lapply(unique(data_raw$Well), function(idx_well){
      
      #idx_well <- unique(data_raw$Well)[5] #FIXME delete me 
      
      data_well_raw <- subset(data_raw[,c("Time", "Conf", "Well")], Well == idx_well)
      
      min_scan_Time <- ceiling(min(data_well_raw$Time))
      
      max_scan_Time <- floor(max(data_well_raw$Time))
      
      model_pred_poly <- get_growthMetrics(data_well_raw, degree = 7)[[3]]
      
      time_sequence <- seq(min_scan_Time ,max_scan_Time,1)
      
      model_well_pred <- predict(model_pred_poly, data.frame(Time = time_sequence))
      
      pred_conf <- data.frame(Well = idx_well, Time = time_sequence, Conf = model_well_pred)
      
      well_metrics <- summary(model_pred_poly)
      
      well_metrics <- data.frame(Well = idx_well, adj.r.2 = well_metrics$adj.r.squared)
      
      return(list(pred_conf, well_metrics))
      
    })
  
  # save diagnostic plots
  
  diagnostics_well <- base::do.call(rbind, lapply(fitted_data,  function(x) x[[2]]) ) #for each plate, plot diagnostic plots with the r-squared.. check if any of fits have ploblems
  
  save_plots_directory = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users/mauro/cell_culture_data/190310_largescreen/figures"
  
  filename = paste("fitting-qc", plate_name, gsub(" ", "_",as.character(Sys.time())), ".png", sep = "_")
  
  filename = gsub(pattern = ":", replacement = "_", x = filename)
  
  bad_wells = diagnostics_well$Well[diagnostics_well$adj.r.2 < r.2.threshold]
  
  #print(paste(as.character(bad_wells), plate_name), sep = "_")
  
  png(filename = paste(save_plots_directory, filename , sep = "/"))
  fig <- hist(diagnostics_well$adj.r.2, main = plate_name)
  text( paste0("adj.r.2 < 0.8 = ", length(bad_wells)) , x = fig$breaks[3], y = max(fig$counts))
  dev.off()
  
  # return a df with the right columns
  
  pred_conf <- base::do.call(rbind, lapply(fitted_data,  function(x) x[[1]]) ) #for each plate, plot diagnostic plots with the R-squared.. check if any of fits have ploblems
  
  pred_conf <- subset(pred_conf, !(Well %in% bad_wells))
  
  pred_conf <- inner_join(pred_conf, metadata_df, by = grouping_vars)
  
  return(pred_conf)
  
}) -> data_corrected

names(data_corrected) <- fileNames

# combining metadata of source.plates into growth_data ######

source_layout <- # import source plate layout data
  
  lapply(list.files("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\cpd_data\\large_screen_plate_layout",
                    pattern = "randomized",
                    full.names = T,
                    recursive = T),
         function(x) {
           
           #x =  "\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\cpd_data\\large_screen_plate_layout/190806_NewMSP_Layout/randomized_layout_1MSP_batch2.xls"
           
           tmp_data  = readxl::read_xls(x)
           
           tmp_data$Well = paste0(tmp_data$Row, tmp_data$Column)
           
           return(tmp_data)
           
         })

data_corrected <-
  
  lapply(names(data_corrected), function(x){
    
    stopifnot(require(dplyr))
    
    #x = names(data_corrected)[1] 
    
    tmp_data = data_corrected[[x]]
    
    match_source <- source_plates[source_plates$filenames == x,c("sourceid", "batch")] 
    
    if(match_source$batch == "batch_1"){
      
      if(grepl(pattern = "1MSP", x = match_source$sourceid)){
        
        tmp_data <- inner_join(tmp_data, source_layout[[1]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
        return(tmp_data)
        
      }else if(grepl(pattern = "2MSP", x =match_source$sourceid)){
        
        tmp_data <- inner_join(tmp_data, source_layout[[2]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
        return(tmp_data)
        
      }else{
        
        stop("some experimental plate could not be matched to a source plate layout.")
        
      }
      
    } else if(match_source$batch == "batch_2"){
      
      if(grepl(pattern = "1MSP", x = match_source$sourceid)){
        
        tmp_data <- inner_join( tmp_data, source_layout[[3]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
        return(tmp_data)
        
      }else if(grepl(pattern = "2MSP", x =match_source$sourceid)){
        
        tmp_data <- inner_join(  tmp_data, source_layout[[4]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
        
        return(tmp_data)
        
      }else{
        
        stop("some experimental plate could not be matched to a source plate layout.")
        
      }
      
    }else{
      
      stop("some information is wrong with the match_source$source_id")
      
    }
    
  } )

names(data_corrected) <- fileNames

# relabelling bad wells based on Echo transfer. #########
# In future I can relabel those with medium, to increase the number of controls, today as exception

data_corrected <-
  
  lapply(names(data_corrected), function(x){
    
    #x = "20191029/results/OVCAR4_CL3_P2.txt"
    
    source_id <- unlist(subset(source_plates, filenames == x, sourceid))
    
    exception_wells <- subset(exceptions, Destination.Plate.Name  == source_id, Destination.Well, drop = T)
    
    tmp_data <- data_corrected[[x]]
    
    tmp_data[(tmp_data$Well %in% exception_wells), c("Drug")] <- c("exception")
    
    tmp_data[(tmp_data$Well %in% exception_wells), c("Final_conc_uM")] <- c(333)
    
    tmp_data$cell<- strsplit(x, split = "/")[[1]][3]
    
    tmp_data$cell <- strsplit(tmp_data$cell, split = "_")[[1]][1]
    
    tmp <- strsplit(x, split = "_")[[1]][3]
    
    tmp_data$source_plate <- strsplit(tmp, split = ".", fixed = T)[[1]][1]
    
    return(tmp_data)
    
  })

names(data_corrected) <- fileNames

data_corrected <- do.call(rbind, data_corrected)

data_corrected <- data_corrected[!(is.na(data_corrected$Final_conc_uM)), ]

# removing DMSO that is not 33.3 and 3.33, as it causes many data analysis problems ######

data_corrected <- subset(data_corrected, !(Drug == "DMSO" & !(Final_conc_uM %in% c(333,367)))) #FIXME include DMSO that increased volume

rm(list = ls()[!ls()=="data_corrected"])

