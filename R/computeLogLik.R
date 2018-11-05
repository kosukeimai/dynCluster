#' A function to compute hold-out log-likelihood
#' 
#' @param n.cluster The number of clusters in the fitted model (integer).
#' @param n.period The number of periods in the fitted model (integer).
#' @param file.dir The file directory of dynCluster outputs.
#' @param test.df A dataframe that contains the hold-out data (D_ijtk and W_ijtk) 
#' created using the function computeDsWs.
#' @param test.run A logical statement indicating whether to conduct a quick test 
#' run based on only two dyads (100_104 and 156_840).
#' @return A list object that contains a scalar for the overall hold-out log-likelihood
#' and a dataframe that contains the hold-out log-likelihood for each dyad.
#' @examples
#' system.time(
#'   out <- computeLogLik(n.cluster = 3, n.period = 5, 
#'                        file.dir = "~/Dropbox/out/", test.df = df, 
#'                        test.run = FALSE)
#' )
computeLogLik <- function(n.cluster, n.period, file.dir, test.df, test.run = FALSE){
  
  ##############################################################################
  ## set up
  ##############################################################################
  # load packages
  pkg <- c("foreign",
           "tidyverse")
  lapply(pkg, require, character.only = TRUE)
  
  # create vector of clusters and periods (starts from zero in C++)
  cluster.vec <- (seq(1, n.cluster) - 1)
  cluster.vec
  period.vec <- (seq(1, n.decades, 1) - 1)
  period.vec
  
  
  ##############################################################################
  ## read in output: pie for last decade in fitted model (pie_zt-1)
  ##############################################################################
  # compile output into dataframe format
  pie.df <- map_df(cluster.vec, function(y) {
    file <- paste(file.dir, "/PI_cluster_T", (n.decades - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        mutate(cluster = y) %>%
        rename(pie_z = V1)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  

  ##############################################################################
  ## read in output: Transition Matrix (P_zt-1_zt)
  ##############################################################################
  PMAT.FILE <- paste(file.dir, "/CHECKPOINT_PMAT.txt", sep = "")
  PMAT.FILE
  p.matrix <- read.delim(PMAT.FILE, header = FALSE, check.names = FALSE, sep = "")
  p.matrix <- as.matrix(p.matrix)
  colnames(p.matrix) <- NULL
  rownames(p.matrix) <- NULL
  p.matrix
  rowSums(p.matrix)
  
  # sum rows of equivalent clusters (flipping problem)
  p.matrix.r1 <- p.matrix[1:n.cluster,]
  p.matrix.r1
  
  p.matrix.r2 <- p.matrix[(n.cluster + 1):(n.cluster*2),]
  p.matrix.r2
  
  p.matrix.r <- p.matrix.r1 + p.matrix.r2
  p.matrix.r
  rowSums(p.matrix.r)
  
  # normalize
  p.matrix.r.n <- p.matrix.r/rowSums(p.matrix.r)
  p.matrix.r.n
  rowSums(p.matrix.r.n)
  
  # split matrix and collapse
  p.matrix.c1 <- p.matrix.r.n[,1:n.cluster]
  p.matrix.c1
  
  p.matrix.c2 <- p.matrix.r.n[,(n.cluster + 1):(n.cluster*2)]
  p.matrix.c2
  
  # sum matrix
  p.matrix.collapse <- p.matrix.c1 + p.matrix.c2
  p.matrix.collapse
  rowSums(p.matrix.collapse)
  
  
  ##############################################################################
  ## read in output: q_kzt doesn't change over periods
  ##############################################################################
  # read into df dyads for each cluster 
  q.df <- map_df(cluster.vec, function(y) {
    
    file <- paste(file.dir, "/Q_cluster_T", (n.decades - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        mutate(cluster = y) %>%
        rename(sitc = V1,
               q_kz = V2)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  
  n.sitc.q <- n_distinct(q.df$sitc)
  n.sitc.q
  
  q.df <- q.df %>% 
    group_by(cluster) %>%
    mutate(direction = ifelse(row_number() <= n.sitc.q, "", ".1")) %>%
    ungroup() %>%
    mutate(sitc_dir = paste(sitc, direction, sep = "")) %>%
    select(-direction)
  
  
  ##############################################################################
  ## read in output: mu_kzt doesn't change over periods
  ##############################################################################
  # read into df dyads for each cluster 
  mu.df <- map_df(cluster.vec, function(y) {
    
    file <- paste(file.dir, "/MU_cluster_T", (n.decades - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        mutate(cluster = y) %>%
        rename(sitc = V1,
               mu_kz = V2)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  
  n.sitc.mu <- n_distinct(mu.df$sitc)
  n.sitc.mu
  
  mu.df <- mu.df %>% 
    group_by(cluster) %>%
    mutate(direction = ifelse(row_number() <= n.sitc.mu, "", ".1")) %>%
    ungroup() %>%
    mutate(sitc_dir = paste(sitc, direction, sep = "")) %>%
    select(-direction)
  
  
  ################################################################################
  ## read in output: sigma_kzt, doesn't change over periods
  ################################################################################
  # read into df dyads for each cluster 
  sigma.df <- map_df(cluster.vec, function(y) {
    
    file <- paste(file.dir, "/SIGMA_cluster_T", (n.decades - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        mutate(cluster = y) %>%
        rename(sitc = V1,
               sigma_kz = V2)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  
  n.sitc.sigma <- n_distinct(sigma.df$sitc)
  n.sitc.sigma
  
  sigma.df <- sigma.df %>% 
    group_by(cluster) %>%
    mutate(direction = ifelse(row_number() <= n.sitc.sigma, "", ".1")) %>%
    ungroup() %>%
    mutate(sitc_dir = paste(sitc, direction, sep = "")) %>%
    select(-direction)
  
  
  ##############################################################################
  ## merge estimation data
  ##############################################################################
  est.df <- left_join(q.df, mu.df, by = c("sitc_dir", "cluster"))
  est.df <- left_join(est.df, sigma.df, by = c("sitc_dir", "cluster"))
  
  est.df <- est.df %>%
    select(sitc_dir, sitc, cluster, q_kz, mu_kz, sigma_kz) %>%
    mutate(sitc = as.character(sitc))
  
  head(est.df)
  
  
  ##############################################################################
  ## compute hold-out sample log-likelihood
  ##############################################################################
  if(test.run){
    
    obs.df <- test.df %>% 
      filter(dyad_id == "100_104" | dyad_id == "156_840") %>%
      arrange(dyad_id, importer_ISO, exporter_ISO)
    
  } else{
    
    obs.df <- test.df
    
  }
  
  out.z1 <- map_df(cluster.vec, function(z1) {
    
    out.z2 <- map_df(cluster.vec, function(z2) {
      
      cat("Computing estimates for transitions from cluster ", z1, " to ", z2, "\n")
      
      # extract pie_zt-1 (scalar)
      pie.out <- pie.df %>% 
        filter(cluster == z1) %>%
        pull(pie_z)
      
      # extract transition probability (scalar)
      p.mat.out <- p.matrix.collapse[z1 + 1, z2 + 1]
      
      # merge observed data from hold-out period and estimates from fitted model (1250 rows) for specific cluster
      merge.df.z2 <- left_join(obs.df,
                               est.df %>% filter(cluster == z2) %>% select(-sitc), 
                               by = c("sitc_dir"))
      
      # compute and create vars
      merge.df.z2 <- merge.df.z2 %>%
        mutate(norm_density = dnorm(W_ijtk, mu_kz, sigma_kz)^(1 - D_ijtk),
               prod = (q_kz^D_ijtk)*((1 - q_kz)^(1 - D_ijtk))*norm_density) %>%
        group_by(dyad_id) %>%
        summarize(prod_2k = prod(prod, na.rm = TRUE)) %>%
        mutate(all_prod = pie.out*p.mat.out*prod_2k,
               cluster_z_T_ij = z2)
      
    })
    
    out.z2 <- out.z2 %>%
      mutate(cluster_z_T_ij_m1 = z1) %>%
      select(dyad_id, cluster_z_T_ij_m1, cluster_z_T_ij, all_prod) %>%
      ungroup() %>%
      arrange(dyad_id)
  })
  
  ll.out <- out.z1 %>%
    group_by(dyad_id) %>%
    summarize(sum_trans = sum(all_prod),
              log_lik_dyad = log(sum_trans))
  
  list(log.lik = sum(ll.out$log_lik_dyad[!is.infinite(ll.out$log_lik_dyad)], 
                     na.rm = TRUE),
       ll.df = ll.out)
  
}
