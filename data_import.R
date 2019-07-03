library(ggplot2)

# import data

source("C:\\Users\\masierom\\polybox\\Programing\\Tecan_\\CellCultureAnalysis.R")

setwd("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cell_culture_data\\190310_Pre_screen_test_CPDdilution\\20190513\\results")

fileNames <- list.files(full.names = F, pattern = ".txt")

data <- 

lapply(fileNames, function(x) {
  
  df <-
  ReadIncuCyteData(read_platemap = F, FileName_IncuCyte = x, Plate_size = 384, FileDirectory = getwd())
  
  df$Time <- as.numeric(df$Time)
  
  return(df)
  
  })

names(data) <- fileNames

neval <- function(month) {
  deparse(match.call()$month)
}


analyze_growth_zero <-
  function(data_source){
    
    
data_source <- filtered_data[filtered_data$CELL == cell & filtered_data$Well == well,c("Time", "Conf")]
    

# Modelling confluence over time ------------------------------------------
    
data <- list()

data$store_vars <- data.frame(time = data_source[1],
                              confluence = data_source[2])

names(data$store_vars) <- c("time", "confluence")

# variables 

poly_degree = 5 # degrees of polynomial fit

# Modelling confluence wtih a polynomial regression model ---------

data$growth_model <- lm(confluence ~ poly(time, degree = poly_degree, raw = T), data$store_vars)

# we can plot it to see if the fit was good

data$store_vars$fitted <- fitted(data$growth_model) #this function gives the predicted confluence based on the original time vector supplied

#plot(data$store_vars$time, data$store_vars$confluence, type = "p") # plotting confluence over time
#lines(data$store_vars$time, data$store_vars$fitted, type = "l", col = "red") #plotting fitted line

polynomial_deriv <- function(x) {  #this function takes the first derivative, or the growth rate equation
    
    terms <- coef(x)
    
    stopifnot(names(terms)[1]=="(Intercept)")
    
    filter_terms <- terms[-1] #remove intercept
    
    stopifnot(all(grepl("^poly", names(filter_terms))))
    
    degrees <- as.numeric(gsub(pattern = "poly\\(.*\\)", replacement = "", names(filter_terms)))
    
    terms_deriv <- setNames(c(filter_terms * degrees, 0), names(terms))
    
    terms_deriv[is.na(terms_deriv)] <- 0
    
    return(matrix(terms_deriv, ncol = 1))
    
  }
  
  data$store_vars$growth_rate <-
    model.matrix(data$growth_model) %*% polynomial_deriv(data$growth_model)
  
  return(diff(sign(data$store_vars$growth_rate)))
  
}


# fixing time problems in dataset ---------------------------------------

# HCT15 did not record some well, so I am correcting for this

data$HCT15_P1.txt$Time <- data$HCT15_P1.txt$Time + 
  (max(data$A549_P1.txt$Time) - max(data$HCT15_P1.txt$Time))

# filtering wells ---------------------------------------------------------

#removing bad wells


store_outlier <-

  sapply(names(data), USE.NAMES = T, function(x){
    
    plate_name = x
    
    x = data[[x]]
    
    good_wells <- 
    filter_growth_outliers(plate_name = plate_name,data = x, time_control = 24, save_diag_plots = T, 
                           save_plots_directory = "\\\\imsbnas.d.ethz.ch/sauer1/users/Mauro/Cell_culture_data/190310_Pre_screen_test_CPDdilution/figures")
    
    new_data <-
    x[!(x[,"Well"] %in% good_wells$wells_exception),]
    
    return(list(new_data, exception_number = good_wells$wells_exception_number))
    
    })

store_outlier[[1]]

#lapply(filtered_data, function(x) plot_multi_well(input_df =  x, max_timepoint = 24))


# removing wells that Echo did not transfer 

source("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cell_culture_data\\190310_Pre_screen_test_CPDdilution\\exceptions\\log_processing.R")

source_plates <- data.frame(fileNames = list.files("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cell_culture_data\\190310_Pre_screen_test_CPDdilution\\20190513\\results"))

source_plates$sourceID <- ifelse(grepl(pattern = "P1", source_plates$fileNames), "1MSP001", "2MSP001")

filtered_data <- 

  lapply(names(data), function(x){
    
    #x = "A549_P1.txt"
    
    source_ID <- unlist(subset(source_plates, fileNames == x, sourceID))
    
    exception_wells <- subset(exceptions, Destination.Plate.Name  == source_ID, Destination.Well, drop = T)
    
    tmp_data <- filtered_data[[x]]
    
    tmp_data <- tmp_data[!(tmp_data$Well %in% exception_wells),]
    
    tmp_data$CELL <- strsplit(x, split = "_")[[1]][1]
    
    tmp_data$Source_Plate <- strsplit(x, split = "_")[[1]][2]
    
    return(tmp_data)
    
  })


# combining metadata of source.plates into growth_data


source_layout <-# import source plate layout data
  
  lapply(list.files("\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cpd_data\\190508_MSP_Layout",
                    pattern = "randomized",
                    full.names = T), 
         function(x) {
           
           tmp_data  = readxl::read_xls(x)
           
           tmp_data$Well = paste0(tmp_data$Row, tmp_data$Column)
           
           return(tmp_data)
           
         })


filtered_data <-

lapply(filtered_data, function(x){
  
  stopifnot(require(dplyr))
  
  #x = filtered_data[[1]]
  
  if(grepl(pattern = "P1", x = x$Source_Plate[1])){
    
    tmp_data <- inner_join( data.frame(x), source_layout[[1]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
    
    return(tmp_data)
    
  }else if(grepl(pattern = "P2", x = x$Source_Plate[1])){
    
    tmp_data <- inner_join( data.frame(x), source_layout[[2]][,c("Well","Drug" , "Final_conc_uM")], by = "Well")
    
    return(tmp_data)
    
  }else{
    
    stop("Some experimental plate could not be matched to a source plate layout.")
    
  }
    
} )



# combine dataframes

filtered_data <- do.call(rbind, filtered_data)
filtered_data <- filtered_data[!(is.na(filtered_data$Drug)),]


(plot <- 

ggplot(subset(filtered_data, Source_Plate= "P1.txt"), aes(x = Time, y = Conf, color = factor(Final_conc_uM)))+
  geom_smooth(formula = y ~ poly(x, 6), method = "lm")+
  facet_wrap(~Drug + CELL, nrow = 5)+
  ggplot2::theme_minimal()+
  geom_vline(xintercept = 48, colour = "red"))
  
ggsave(plot = plot, filename = "AllDoseResponsePlots.png",
       path = "\\\\imsbnas.d.ethz.ch\\sauer1\\users\\Mauro\\Cell_culture_data\\190310_Pre_screen_test_CPDdilution\\figures", device = "png",
       width = 30, height = 15)



# calculating second derivative, where is first zero?


filtered_data <- subset(filtered_data, Source_Plate == "P1.txt" & Time <= 72)

growthSolution <- filtered_data

for(cell in unique(filtered_data$CELL)){
  
  for(well in unique(filtered_data$Well)){
    
    
    
    print(cell)
    print(well)
    
    if(nrow(filtered_data[filtered_data$CELL == cell & filtered_data$Well == well,c("Time", "Conf")])== 0){
      next()
    }
    
    store <- 
    
    analyze_growth_zero(
    
    filtered_data[filtered_data$CELL == cell & filtered_data$Well == well,c("Time", "Conf")]
    )
    
    
    store <- store[!(is.na(store))]
    
    store <- sum(!(store == 0))
    
    growthSolution[growthSolution$Well == well & growthSolution$CELL == cell,"Conf"] <- store
    
  }
  
}

growthSolution <- growthSolution %>% group_by(Well, CELL) %>% slice(1)

ggplot(growthSolution, aes(CELL, Drug, fill = Conf))+
  geom_tile()+
  scale_fill_gradient(low = "blue", high = "red")+
  
