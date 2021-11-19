## import R functions ####

library(ggplot2); library(RColorBrewer)

# Importing plate times ####  

source("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\analysis\\time_analysis.R")

# import functions cell culture #####

sapply(list.files("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\",full.names=T),source,.GlobalEnv)

# Importing exceptions ####

source("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\exceptions\\log_processing.r")

## paths & parameters ##

path_clean_data <- "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data"
path_fig <- "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users/mauro/cell_culture_data/190310_largescreen/figures"


correct_for_initial_seeding = T
skip_first_time = F
r.2.threshold = 0.8
slope_cutoff = 0.05
p.val_cutoff = 0.05
intercept_sd_cutoff = 2.5
poly_degree = 8


## Importing source plates #####

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\growthData")

source_plates <- data.frame(
  filenames = list.files(
    path = "\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\Growthdata\\",
    pattern = "[P][1-2][.txt]",
    recursive = T))

source_plates$sourceid <- ifelse(grepl(pattern = "P1", source_plates$filenames), "1MSP001", "2MSP001")

source_plates$batch <- ifelse(grepl(pattern = "20190513", source_plates$filenames), "batch_1", "batch_2")

setwd("..")

tmp <- read.xlsx("Description.xlsx")

source_plates$uniqueID <- NA

cell <- strsplit(x = as.character(source_plates$filenames), split = "/")

cell <- lapply(cell, "[[", 3)

cell <- strsplit(x = unlist(cell), split = "_")

plate <- unlist(lapply(cell, "[", 3))

cell <- unlist(lapply(cell, "[[", 1))

tmp <- tmp[match(x = cell, table = tmp$cell_name), ]

msp1 <- grepl(pattern = "1MSP", x = source_plates$sourceid)

source_plates$uniqueID[msp1] <- tmp$source1[msp1]

msp2 <- grepl(pattern = "2MSP", x = source_plates$sourceid)

source_plates$uniqueID[msp2] <- tmp$source2[msp2]

rm(tmp, msp1, msp2)

setwd("./growthData")

#import growth curves ####

fileNames <- as.character(source_plates$filenames)

fileNames<- fileNames[7:length(fileNames)] #remove CLs that do not have plate 2

data <- 
  
  lapply(fileNames, function(x) {
    
    #x = fileNames[1]
    print(fileNames)
    
    df <-
      ReadIncuCyteData(read_platemap = F, FileName_IncuCyte = x,
                       Plate_size = 384, FileDirectory = getwd(),
                       time_output = "GMT",skip_first_time=skip_first_time,
                       correct_init_seeding = correct_for_initial_seeding)
    
    
    return(df)
    
  })

names(data) <- fileNames

#### load data ####

### Create elapsed and sampling times based on time_vectors for the 384 plates ######

data <- lapply(names(data), function(filename){
  
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

tmp_names <- names(data_corrected)

data_corrected <-
  
  lapply(names(data_corrected), FUN = function(x){
    
    print(x)
    
    plate_name = strsplit(x, split = "/")[[1]][3]
    
    x = data_corrected[[x]]
    
    good_wells <-
      filter_growth_outliers(plate_name = paste(plate_name,"corr-init-seed=",correct_for_initial_seeding,'remove-first-time',skip_first_time,sep = "_"), data = x, time_control = 0, save_diag_plots = F,
                             slope_cutoff = slope_cutoff,
                             p.val_cutoff = p.val_cutoff,
                             intercept_sd_cutoff = intercept_sd_cutoff,
                             return_stats = T,
                             save_plots_directory = paste(path_fig,"exclusions-growth", sep = "\\"))
    
    
    new_data <-
      x[!(x[,"Well"] %in% good_wells[[1]]$wells_exception),]
    
    print(paste("Plate:", plate_name,
                "Exception-Number:",good_wells[[1]]$wells_exception_number
    ))
    
    return(list(new_data,good_wells[[2]]))
    
    })

#save summary of plate outliers to excel

setwd(paste0(path_clean_data,"/growth_data"))

write.csv(do.call(rbind,lapply(data_corrected,"[[",2)),
          paste("growth_data_exclusions_growth","corr-init-seed=",correct_for_initial_seeding,'remove-first-time',skip_first_time,".csv",sep = "_"))

data_corrected <- lapply(data_corrected,"[[",1)

names(data_corrected) <- tmp_names

data_unfitted <- data_corrected

for(poly_degree in poly_degree){
    
  data_corrected = data_unfitted
  
  # Fit the confluence for each well, and return fitted confluence. Keep the same plate map structure. ####
  
  exclusions_noise <- data.frame("plate"=character(),"exclusions"=integer())
  
  base::lapply(names(data_corrected), function(plate_name){
    
    print(plate_name)
    
    data_raw <- data_corrected[[plate_name]]
    
    plate_name <- strsplit(plate_name, split = "/")[[1]][3]
    
    grouping_vars <- "Well"
    
    metadata_cols <- colnames(data_raw)[!colnames(data_raw) %in% c("Time", "Conf")]
    
    metadata_df <- data_raw[,metadata_cols]
    
    metadata_df <- metadata_df %>% dplyr::group_by(get(grouping_vars)) %>% dplyr::slice(1) %>% ungroup()
    
    metadata_df[ncol(metadata_df)] <- NULL
    
    fitted_data <- 
      
      lapply(unique(data_raw$Well), function(idx_well){
        
        #idx_well <- unique(data_raw$Well)[5] #FIXME delete me 
        
        #idx_well = 'E11'
        
        data_well_raw <- subset(data_raw[,c("Time", "Conf", "Well")], Well == idx_well)
        
        min_scan_Time <- ceiling(min(data_well_raw$Time))
        
        max_scan_Time <- floor(max(data_well_raw$Time))
        
        model <- get_growthMetrics(data = data_well_raw, degree = poly_degree)
        
      
        time_12 <- model[[2]][4][which(abs(model[[2]][4] - 12) == min(abs(model[[2]][4] - 12))),]
        
        time_24 <- model[[2]][4][which(abs(model[[2]][4] - 24) == min(abs(model[[2]][4] - 24))),]
        
        time_34 <- model[[2]][4][which(abs(model[[2]][4] - 34) == min(abs(model[[2]][4] - 34))),]
        
        time_1h_window <- c(time_24-0.5, time_24+0.5)
        
        time_2h_window <-c(time_24-1, time_24+1)
        
        time_3h_window <- c(time_24-1.5, time_24+1.5)
        
        GR_12 <- subset(model[[2]],time==time_12,select = 'GR',drop = T)
        
        GR_24 <- subset(model[[2]],time==time_24,select = 'GR',drop = T)
        
        GR_34 <- subset(model[[2]],time==time_34,select = 'GR',drop = T)
        
        GR_1h_window <- mean(subset(model[[2]],time >=time_1h_window[1]&time <=time_1h_window[2],select = 'GR',drop = T))
        
        GR_2h_window <- mean(subset(model[[2]],time >=time_2h_window[1]&time <=time_2h_window[2],select = 'GR',drop = T))
        
        GR_3h_window <- mean(subset(model[[2]],time >=time_3h_window[1]&time <=time_3h_window[2],select = 'GR',drop = T))
        
        model_pred_poly <- model[[3]]
        
        time_sequence <- seq(min_scan_Time ,max_scan_Time,1)
        
        model_well_pred <- predict(model_pred_poly, data.frame(Time = time_sequence))
        
        pred_conf <- data.frame(Well = idx_well, Time = time_sequence, Conf = model_well_pred)
        
        pred_conf$GR12 <-         GR_12
        pred_conf$GR24 <-         GR_24
        pred_conf$GR34 <-         GR_34
        pred_conf$GR_1h_window <- GR_1h_window
        pred_conf$GR_2h_window <- GR_2h_window
        pred_conf$GR_3h_window <- GR_3h_window
        
      
        well_metrics <- summary(model_pred_poly)
        
        well_metrics <- data.frame(Well = idx_well, adj.r.2 = well_metrics$adj.r.squared)
        
        return(list(pred_conf, well_metrics))
        
      })
    
    # save diagnostic plots
    
    diagnostics_well <- base::do.call(rbind, lapply(fitted_data,  function(x) x[[2]]) ) #for each plate, plot diagnostic plots with the r-squared.. check if any of fits have ploblems
    
    filename = paste("fitting-qc", plate_name, gsub(" ", "_",as.character(Sys.time())), ".png", sep = "_")
    
    filename = gsub(pattern = ":", replacement = "_", x = filename)
    
    bad_wells = diagnostics_well$Well[diagnostics_well$adj.r.2 < r.2.threshold]
    
    exclusions_noise[plate_name,] = c(plate_name,length(bad_wells))
    
    #print(paste(as.character(bad_wells), plate_name), sep = "_")
    
    png(filename = paste(path_fig,"exclusions-noise", filename , sep = "/"))
    fig <- hist(diagnostics_well$adj.r.2, main = plate_name)
    text( paste0("adj.r.2 < 0.8 = ", length(bad_wells)) , x = fig$breaks[3], y = max(fig$counts))
    dev.off()
    
    # return a df with the conf + metadata columns
    
    pred_conf <- base::do.call(rbind, lapply(fitted_data,  function(x) x[[1]]) ) #for each plate, plot diagnostic plots with the R-squared.. check if any of fits have ploblems
    
    pred_conf <- subset(pred_conf, !(Well %in% bad_wells))
    
    pred_conf <- inner_join(pred_conf, metadata_df, by = grouping_vars)
    
    #return a df with 
    
    
    return(list(pred_conf, exclusions_noise))
    
  }) -> data_corrected
  
  #save exclusions based on r^2
  
  setwd(paste0(path_clean_data,"/growth_data"))
  write.csv(do.call(rbind,lapply(data_corrected,"[[",2)),paste0("growth_data_exclusions_noise","polydegree=",poly_degree,".csv"))
  data_corrected <- lapply(data_corrected,function(x){
    return(x[[1]])
    })
  names(data_corrected) <- fileNames
  
  # combining metadata of source.plates into growth_data ######
  
  #source_layout[[1:2]] one and two refer to the first batch (plate 1 and 2), which had problems in Echo transfer.
  #source_layout[[3:4]] three and four refer to the second GOOD batch (plate 3 and 4), which were used in the big screen.
  
  
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
  
  data_corrected <-
    
    lapply(names(data_corrected), function(x){
      
      #x = "20191106/results/LOXIMVI_CL2_P1.txt"
      
      source_id <- unlist(subset(source_plates, filenames == x, uniqueID))
      
      exception_wells <- subset(exceptions[[2]], cellPlateBC  == source_id, WellNameTransferError, drop = T)
      
      exception_wells <-  unlist(strsplit(exception_wells, split = ","))
      
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
  
  setwd(paste0(path_clean_data,"/growth_data"))
  
  #remove unsused wells (wells without drug in drug source plate).
  
  data_unused_wells <- data_corrected[(is.na(data_corrected$Drug)), ]
  
  data_unused_wells <- rbind(data_unused_wells, data_corrected[(is.na(data_corrected$Final_conc_uM)), ])
  
  data_corrected <- data_corrected[!(is.na(data_corrected$Drug)),]
  
  data_corrected <- data_corrected[!(is.na(data_corrected$Final_conc_uM)),]
  
  # removing DMSO that is not 33.3 and 3.33, as it causes many data analysis problems ######
  
  data_corrected <- subset(data_corrected, !(Drug == "DMSO" & !(Final_conc_uM %in% c(333,367))))
  
  data_corrected <- subset(data_corrected, !(Drug ==  "exception"))

}

#exclude groups with <3 reps

data_drugs <- subset(data_corrected,!Drug %in% c("PBS","DMSO") & Time == 24)

exclusions_reps <- data_drugs %>% group_by(cell, Drug, Final_conc_uM) %>% summarise(nreps = n())

exclusions_reps <- subset(exclusions_reps, nreps <3)

write.csv(exclusions_reps, "exclusion_nreps.csv")

groups_to_exclude <- paste(exclusions_reps$cell,exclusions_reps$Drug, exclusions_reps$Final_conc_uM,sep = "_")

groups <- paste(data_corrected$cell,data_corrected$Drug, data_corrected$Final_conc_uM,sep = "_")

data_corrected <- data_corrected[!groups %in% groups_to_exclude,]

#save final outcome

rm(list = ls()[!ls()=="data_corrected"])

setwd("\\\\d.ethz.ch\\groups\\biol\\sysbc\\sauer_1\\users\\Mauro\\Cell_culture_data\\190310_LargeScreen\\clean_data")

save(data_corrected, file = "growth_curves_filtered.RData")
