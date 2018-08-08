################################################################################
# function for UN comtrade concordence
################################################################################
getCtyName <- function(id.vec, concord.df){
  # Prepare output vector
  out.vec <- rep(NA, length(id.vec))
  matches <- match(id.vec, concord.df$`Country.Code`)
  out.vec <- concord.df[matches,]$`Country.Name.English`
  
  if(any(is.na(out.vec))){
    warning(paste("concordence missing: ", 
                  paste(id.vec[is.na(out.vec)], collapse = ", "), 
                  sep = ""))
  } else{
    #print("All variables successfully converted")
  }
  
  return(out.vec)
}


################################################################################
# function for cluster concordence across different models
################################################################################
concordClusters <- function(cluster.vec, folder, concord.table){
  
  folder.name <- str_replace(folder, "/", "")
  concord.table.sub <- concord.table %>%
    filter(folder == folder.name)
  
  if(nrow(concord.table.sub) == 0){
    stop("Folder not matched")
  }
  
  # Prepare output vector
  out.vec <- rep(NA, length(cluster.vec))
  matches <- match(cluster.vec, concord.table.sub$cluster)
  out.vec <- concord.table.sub[matches,]$cluster_set
  
  if(any(is.na(out.vec))){
    stop(paste("Cluster concordence missing: ", 
               paste(cluster.vec[is.na(out.vec)], collapse = ", "), 
               sep = ""))
  } else{
    #print("All variables successfully converted")
  }
  
  return(out.vec)
}


################################################################################
# function for cluster-name concordence across different models
################################################################################
concordClusterNames <- function(cluster.vec, folder, concord.table){
  
  folder.name <- str_replace(folder, "/", "")
  concord.table.sub <- concord.table %>%
    filter(folder == folder.name)
  
  # Prepare output vector
  out.vec <- rep(NA, length(cluster.vec))
  matches <- match(cluster.vec, concord.table.sub$cluster)
  out.vec <- as.character(concord.table.sub[matches,]$cluster_set_name)
  
  if(any(is.na(out.vec))){
    stop(paste("Cluster-Name concordence missing: ", 
               paste(cluster.vec[is.na(out.vec)], collapse = ", "), 
               sep = ""))
  } else{
    #print("All variables successfully converted")
  }
  
  return(out.vec)
}


################################################################################
# function for normalizing trade by products
################################################################################
genTF_row <- function(imp, exp){
  impexp <- which(dyad_data$importer_ISO == imp & dyad_data$exporter_ISO == exp)
  expimp <- which(dyad_data$importer_ISO == exp & dyad_data$exporter_ISO == imp)
  
  if(length(impexp) == 0){
    imp_exp <- rep(0, (ncol(dyad_data)-4))
  } else {
    imp_exp <- as.numeric(dyad_data[impexp,-c(1,2,3,ncol(dyad_data))])
  }
  
  if(length(expimp) == 0){
    exp_imp <- rep(0, (ncol(dyad_data)-4))
  } else {
    exp_imp <- as.numeric(dyad_data[expimp,-c(1,2,3,ncol(dyad_data))])
  }
  
  agree <- all(names(dyad_data[,-c(1,2,3,ncol(dyad_data))]) == names(total))
  
  if(agree){
    norm_impexp <- imp_exp/total
    norm_expimp <- exp_imp/total
    row <- c(norm_impexp, norm_expimp)
  } else {
    stop("*** product column names don't agree")
  }
  
  return(row)
  
}


################################################################################
# function for generating dyadic ID
################################################################################
genID <- function(imp_iso, exp_iso){
  imp.n <- as.numeric(imp_iso)
  exp.n <- as.numeric(exp_iso)
  if (imp.n <= exp.n){
    id <- paste(imp.n, exp.n, sep = "_")
  } else if (imp.n > exp.n) {
    id <- paste(exp.n, imp.n, sep = "_")
  } else {
    cat("error in dyads", imp.n, exp.n, "\n")
    stop("error in country id")
  }
  return(id)
}


################################################################################
# function for generating product description
################################################################################
getDesc <- function(code, version){
  if(version == "S3"){
    version <- "S4"
  }
  desc <- read.csv(DESC_FILE, colClasses = "character", header = TRUE)
  first2 <- substring(code, 1,2)
  second2 <- substring(code, 3,4)
  idx <- which(desc$code == first2 & desc$version == version)
  if(length(idx) == 0){
    idx <- which(desc$code == first2)
    if(length(idx) > 0){
      description <- as.character(paste(desc[idx, c("short_desc")][1], 
                                        second2, sep = "_"))
    } else {
      description <- code
    }
  } else {
    description <- as.character(paste(desc[idx, c("short_desc")], 
                                      second2, sep = "_"))
  }
  description <- as.character(description)
  if(description == "99"){
    description <- "Not classified in the SITC"
  }
  description <- gsub("\\_$", "", description)
  return(description)
}


################################################################################
# function to calculate sum + 1
################################################################################
sumPlusOne <- function(x){
  out <- sum(as.numeric(x), na.rm = TRUE) + 1
  return(out)
}
