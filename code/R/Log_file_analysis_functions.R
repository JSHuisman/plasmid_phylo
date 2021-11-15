#############################################
## Functions needed to analyse BEAST log files
##
## Author: Jana Huisman
## Last Edit: August 2018
#############################################

# For the ESS
library(coda)

# For the HPD intervals
library("HDInterval")

# different error measure
library('entropy')


#############################################
## Essential function to remove Burnin from
## the trace files
#############################################

removeBurnin <- function(df, burninFrac=0.1) {
  n <- dim(df)[1]
  return(df[-(1:ceiling(n*burninFrac)),])
}
######
prep_log_file <- function(log_file){
  df <- try(read.table(log_file, header=TRUE, sep='\t', skipNul = TRUE))
  #df <- try(read.table(log_file, header=TRUE, sep='\t'))
  if (inherits(df, "try-error")) print(log_file)
  df <- removeBurnin(df, burninFrac=0.1)
  return(df)
}

get_log_files <- function(cfg_folder){
  all_files <- list.files(path=cfg_folder, pattern = ".log$", 
                          recursive=FALSE, full.names = TRUE)
  #print(all_files[!file.info(all_files)$isdir])
  return(all_files[!file.info(all_files)$isdir])
}

get_single_log_file <- function(cfg_folder){
  return(get_log_files(cfg_folder)[1])
}

#############################################
## To calculate the effective sampling size (ESS)
## and test whether this is sufficient (i.e. more than 200)
#############################################

ESS_from_df <- function(df){
  compute_ESS_vector <- function(cname) effectiveSize(df[,cname])
  ESS_vector <- sapply(colnames(df), compute_ESS_vector)
  names(ESS_vector) <- colnames(df)
  return(ESS_vector)
}
######
ESS_from_log_file <- function(log_file){
  df <- prep_log_file(log_file)
  ESS_vector <- ESS_from_df(df)
  return(ESS_vector)
}

######
ESS_experiment_test <- function(cfg_folder) {
  
  # load the names of all runs (numbers from 1 to nrun)
  log_files <- get_log_files(cfg_folder)
  
  temp_df <- prep_log_file(log_files[1])
  
  # the matrix to save the ESS values
  ESS <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))
  
  for (file_index in 1:length(log_files)){
    ESS[file_index,] <- ESS_from_log_file(log_files[file_index])
  }
  
  return(data.frame(ESS))
}

######

cfg_folder_validity <- function(cfg_folders){
  
  log_files_present <- vector(mode="logical",length = length(cfg_folders))
  
  for (i in 1:length(cfg_folders)){
    single_log_file <- get_log_files(cfg_folders[i])[1]
    if (length(single_log_file)==0 ){
      log_files_present[i] <- FALSE
    } else{
      log_files_present[i] <- TRUE
      }
  }
  return(cfg_folders[log_files_present])
}



calc_mean_ESS_rates <- function(cfg_folders){
  cfg_folders <- cfg_folder_validity(cfg_folders)
  
  single_log_file <- get_log_files(cfg_folders[1])[1]

  single_log_df <- prep_log_file(single_log_file)
  mean_ESS_rates <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  
  for (i in 1:length(cfg_folders)){
    ESS_rates <- (ESS_experiment_test(cfg_folders[i]))  
    mean_ESS_rates[i,] <- sapply(ESS_rates,mean)
  }
  mean_ESS_rates <- data.frame(mean_ESS_rates)
  
  return(mean_ESS_rates)
}


#############################################
## Automatic calculation of the mean of means
##
#############################################

mean_of_log_file <- function(log_file){
  df <- prep_log_file(log_file)
  means <- colMeans(df)
  return(means)
}

mean_of_means <- function(cfg_folder){
  # load the names of all runs (numbers from 1 to nrun)
  log_files <- get_log_files(cfg_folder)
  temp_df <- prep_log_file(log_files[1])
  
  # the dataframe of lower and upper hpd for this variable
  means_matrix <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))

  for (file_index in 1:length(log_files)){
    means_matrix[file_index,] <- mean_of_log_file(log_files[file_index])
  }
  
  mean_of_means_vector = colMeans(means_matrix)
  
  return(mean_of_means_vector)
}

mean_of_means_all_cfgs <- function(cfg_folders){
  
  single_log_file <- get_log_files(cfg_folders[1])[1]
  single_log_df <- prep_log_file(single_log_file)
  
  mean_of_means_matrix <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  
  for (i in 1:length(cfg_folders)){
    mean_of_means_matrix[i,] <- mean_of_means(cfg_folders[i])
  }
  
  return(data.frame(mean_of_means_matrix))
}

mean_HPD_bounds_all_cfgs <- function(cfg_folders){
  single_log_file <- get_log_files(cfg_folders[1])[1]
  single_log_df <- prep_log_file(single_log_file)
  
  mean_HPD_bound_upper <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  mean_HPD_bound_lower <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  
  for (i in 1:length(cfg_folders)){
    HPD.df <- HPD_matrices_from_cfg_folder(cfg_folders[i])
    mean_HPD_bound_lower[i,] <- colMeans(HPD.df$lower)
    mean_HPD_bound_upper[i,] <- colMeans(HPD.df$upper)
  }
  
  return(list(lower=data.frame(mean_HPD_bound_lower), upper=data.frame(mean_HPD_bound_upper)))
}

#############################################
## Automatic calculation of the HPD interval (upper, lower) based on the trace
##
#############################################

HPD_from_df <- function(df){
  compute_HPD_matrix <- function(cname) hdi(df[,cname], credMass=0.95)
  HPD_matrix <- sapply(colnames(df), compute_HPD_matrix)
  names(HPD_matrix) <- colnames(df)
  return(HPD_matrix)
}

######
HPD_from_log_file <- function(log_file){
  df <- prep_log_file(log_file)
  HPD_matrix <- HPD_from_df(df)
  return(HPD_matrix)
}

######
HPD_matrices_from_cfg_folder <- function(cfg_folder){
  # load the names of all runs (numbers from 1 to nrun)
  log_files <- get_log_files(cfg_folder)
  temp_df <- prep_log_file(log_files[1])

  # the dataframe of lower and upper hpd for this variable
  HPD_lower_matrix <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))
  HPD_upper_matrix <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))
  
  for (file_index in 1:length(log_files)){
    result <- HPD_from_log_file(log_files[file_index])
    HPD_lower_matrix[file_index,] <- result[1,]
    HPD_upper_matrix[file_index,] <- result[2,]
  }
  
  return(list(lower=data.frame(HPD_lower_matrix), upper=data.frame(HPD_upper_matrix)))
}

######
HPD_per_variable <- function(HPD.df, variable){
  variable.df <- data.frame(Lower = HPD.df$lower[,variable], Upper = HPD.df$upper[,variable], stringsAsFactors = F)
  return(variable.df)
}

######
HPD_error <- function(HPD_lower_matrix, HPD_upper_matrix, true_values=true_val) {
  
  error_rate <- vector()
  for (variable in colnames(HPD_lower_matrix)){
    error_low<-apply(HPD_lower_matrix, MARGIN=1, function(x) {x[variable]>true_values[variable]})   
    error_up<-apply(HPD_upper_matrix, MARGIN=1, function(x) {x[variable]<true_values[variable]})   
    # hier wordt het samengevat - sum
    error_rate <- cbind(error_rate,sum(error_low|error_up)/length(HPD_lower_matrix[,1]))
  }
  colnames(error_rate)<- colnames(HPD_lower_matrix)
  return(error_rate)
}

######

HPD_error_all_cfgs <- function(cfg_folders){
  
  single_log_file <- get_log_files(cfg_folders[1])[1]
  single_log_df <- prep_log_file(single_log_file)
  
  HPD_error_rates <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  
  for (i in 1:length(cfg_folders)){
    HPD.df <- HPD_matrices_from_cfg_folder(cfg_folders[i])
    HPD_error_rates[i,] <- HPD_error(HPD.df$lower,HPD.df$upper, true_values=true_val[i,])
  }
  
  error_rate <- data.frame(HPD_error_rates)
  return(error_rate)
}



######
HPD_width_error <- function(HPD_lower_matrix, HPD_upper_matrix) {
  
  # the dataframe of lower and upper hpd for this variable
  HPD_width.df <- HPD_upper_matrix-HPD_lower_matrix
  error<-apply(HPD_width.df, MARGIN=2, mean)   
  
  return(error)
}

HPD_width_error_all_cfgs <- function(cfg_folders){
  
  single_log_file <- get_log_files(cfg_folders[1])[1]
  single_log_df <- prep_log_file(single_log_file)
  
  HPD_width_error_rates <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  
  for (i in 1:length(cfg_folders)){
    HPD.df <- HPD_matrices_from_cfg_folder(cfg_folders[i])
    HPD_width_error_rates[i,] <- HPD_width_error(HPD.df$lower,HPD.df$upper)
  }
  
  error_rate <- data.frame(HPD_width_error_rates)
  return(error_rate)
}


HPD_width_var_all_cfgs <- function(cfg_folders){

  HPD_width_error_rates <- vector('list', length=length(cfg_folders))
  
  for (i in 1:length(cfg_folders)){
    HPD.df <- HPD_matrices_from_cfg_folder(cfg_folders[i])
    HPD_width_error_rates[[i]] <-   data.frame(HPD.df$upper-HPD.df$lower)
  }
  
  return(HPD_width_error_rates)
}


#############################################

ESS_and_HPD_from_log_file <- function(log_file){
  df <- prep_log_file(log_file)
  ESS_vector <- ESS_from_df(df)
  HPD_matrix <- HPD_from_df(df)
  return(list(ESS=ESS_vector, HPD=HPD_matrix))
}
######
ESS_and_HPD_from_cfg_folder <- function(cfg_folder){
  # load the names of all runs (numbers from 1 to nrun)
  log_files <- get_log_files(cfg_folder)
  temp_df <- prep_log_file(log_files[1])
  
  # the dataframe of ESS, lower and upper hpd for this variable
  ESS <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))
  HPD_lower_matrix <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))
  HPD_upper_matrix <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))
  
  for (file_index in 1:length(log_files)){
    result <- ESS_and_HPD_from_log_file(log_files[file_index])
    ESS[file_index,] <- result$ESS
    HPD_lower_matrix[file_index,] <- result$HPD[1,]
    HPD_upper_matrix[file_index,] <- result$HPD[2,]
  }
  
  return(list(ESS=data.frame(ESS), lower=data.frame(HPD_lower_matrix), upper=data.frame(HPD_upper_matrix)))
}

#############################################

Entropy_error <- function(cfg_folder) {
  
  # load the names of all runs (numbers from 1 to nrun)
  log_files <- get_log_files(cfg_folder)
  temp_df <- prep_log_file(log_files[1])
  
  # the matrix of entropy
  Ent_matrix <- matrix(data=NA, nrow=length(log_files), ncol=length(colnames(temp_df)), dimnames = list(log_files = 1:length(log_files), variable = colnames(temp_df)))

  for (i in 1:length(log_files)){
    df <- prep_log_file(log_files[i])
      
    for (variable in colnames(df)){
        normalized = (df[,variable]-min(df[,variable]))/(max(df[,variable])-min(df[,variable]))
        ENT<-entropy.empirical(discretize(normalized, numBins = 50, r = c(0,1)), unit="log2")
        Ent_matrix[i,variable]<-ENT
      }
  }
  
  Ent.df = data.frame(Ent_matrix)
  error<-apply(Ent.df, MARGIN=2, mean)   

  return(error)
}



Entropy_error_all_cfgs <- function(cfg_folders){
  
  single_log_file <- get_log_files(cfg_folders[1])[1]
  single_log_df <- prep_log_file(single_log_file)
  
  Entropy_error_rates <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  
  for (i in 1:length(cfg_folders)){
    Entropy_error_rates[i,] <- Entropy_error(cfg_folders[i])
  }
  
  error_rate <- data.frame(Entropy_error_rates)
  return(error_rate)
}

#############################################
all_errors_all_cfgs <- function(cfg_folders){
  
  single_log_file <- get_log_files(cfg_folders[1])[1]
  single_log_df <- prep_log_file(single_log_file)

  HPD_error_rates <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  HPD_width_error_rates <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  Entropy_error_rates <- matrix(data=NA, nrow=length(cfg_folders), ncol=length(colnames(single_log_df)), dimnames = list(config = 1:length(cfg_folders), variable = colnames(single_log_df)))
  
  for (i in 1:length(cfg_folders)){
    HPD.df <- HPD_matrices_from_cfg_folder(cfg_folders[i])
    HPD_error_rates[i,] <- HPD_error(HPD.df$lower,HPD.df$upper, true_values=true_val[i,])
    HPD_width_error_rates[i,] <- HPD_width_error(HPD.df$lower,HPD.df$upper)
    Entropy_error_rates[i,] <- Entropy_error(cfg_folders[i])
  }
  
  return(list(HPD = data.frame(HPD_error_rates),HPD_width = data.frame(HPD_width_error_rates), Entropy = data.frame(Entropy_error_rates)))
}

