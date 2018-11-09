#' A function to compute hold-out log-likelihood
#' 
#' @param n.cluster The number of clusters in the fitted model (integer).
#' @param n.period The number of periods in the fitted model (integer).
#' @param file.dir The file directory of dynCluster outputs.
#' @param test.df A dataframe that contains the hold-out data (D_ijtk and W_ijtk) 
#' created using the function computeDsWs.
#' @param test.run A logical TRUE/FALSE statement indicating whether to conduct a quick test 
#' run based on only two dyads (100_104 and 156_840).
#' @return A list object that contains a scalar for the overall hold-out log-likelihood
#' and a dataframe that contains the hold-out log-likelihood for each dyad.
#' @import dplyr 
#' @import purrr
#' @import tidyr 
#' @import utils
#' @importFrom stats dnorm
#' @export
#' @examples
#' # set directories
#' test.dir <- "./example/toy/dynamic_1970.csv"
#' out.dir <- "./example/toy"
#' 
#' # extract D_ijtk and W_ijtk
#' df <- computeDsWs(test.dir)
#' 
#' # compute hold-out log-likelihood
#' system.time(
#'   out <- computeLogLik(n.cluster = 3, n.period = 10, 
#'                        file.dir = out.dir, test.df = df, 
#'                        test.run = TRUE)
#' )
computeLogLik <- function(n.cluster, n.period, file.dir, test.df, test.run){
  
  ##############################################################################
  ## set up
  ##############################################################################
  
  # create vector of clusters and periods (starts from zero in C++)
  cluster.vec <- (seq(1, n.cluster) - 1)
  cluster.vec
  period.vec <- (seq(1, n.period, 1) - 1)
  period.vec
  
  
  ##############################################################################
  ## read in output: pie for last decade in fitted model (pie_zt-1)
  ##############################################################################
  # compile output into dataframe format
  pie.df <- purrr::map_df(cluster.vec, function(y) {
    file <- paste(file.dir, "/PI_cluster_T", (n.period - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- utils::read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        dplyr::mutate(cluster = y) %>%
        dplyr::rename(pie_z = V1)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  

  ##############################################################################
  ## read in output: Transition Matrix (P_zt-1_zt)
  ##############################################################################
  PMAT.FILE <- paste(file.dir, "/CHECKPOINT_PMAT.txt", sep = "")
  PMAT.FILE
  p.matrix <- utils::read.delim(PMAT.FILE, header = FALSE, check.names = FALSE, sep = "")
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
  q.df <- purrr::map_df(cluster.vec, function(y) {
    
    file <- paste(file.dir, "/Q_cluster_T", (n.period - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- utils::read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        dplyr::mutate(cluster = y) %>%
        dplyr::rename(sitc = V1,
                      q_kz = V2)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  
  n.sitc.q <- dplyr::n_distinct(q.df$sitc)
  n.sitc.q
  
  q.df <- q.df %>% 
    dplyr::group_by(cluster) %>%
    dplyr::mutate(direction = ifelse(row_number() <= n.sitc.q, "", ".1")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sitc_dir = paste(sitc, direction, sep = "")) %>%
    dplyr::select(-direction)
  
  
  ##############################################################################
  ## read in output: mu_kzt doesn't change over periods
  ##############################################################################
  # read into df dyads for each cluster 
  mu.df <- purrr::map_df(cluster.vec, function(y) {
    
    file <- paste(file.dir, "/MU_cluster_T", (n.period - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- utils::read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        dplyr::mutate(cluster = y) %>%
        dplyr::rename(sitc = V1,
               mu_kz = V2)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  
  n.sitc.mu <- dplyr::n_distinct(mu.df$sitc)
  n.sitc.mu
  
  mu.df <- mu.df %>% 
    dplyr::group_by(cluster) %>%
    dplyr::mutate(direction = ifelse(row_number() <= n.sitc.mu, "", ".1")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sitc_dir = paste(sitc, direction, sep = "")) %>%
    dplyr::select(-direction)
  
  
  ################################################################################
  ## read in output: sigma_kzt, doesn't change over periods
  ################################################################################
  # read into df dyads for each cluster 
  sigma.df <- purrr::map_df(cluster.vec, function(y) {
    
    file <- paste(file.dir, "/SIGMA_cluster_T", (n.period - 1), "_Z", y, ".txt", 
                  sep = "")
    
    tryCatch({
      #if(file.exists(file)) {
      data <- utils::read.delim(file, header = FALSE, check.names = FALSE)
      data <- data %>% 
        dplyr::mutate(cluster = y) %>%
        dplyr::rename(sitc = V1,
               sigma_kz = V2)
    }, error = function(e) NULL) # if no dyads are assigned to a cluster-year, no corresponding txt file will be created, assign NULL and skip
    
  })
  
  n.sitc.sigma <- dplyr::n_distinct(sigma.df$sitc)
  n.sitc.sigma
  
  sigma.df <- sigma.df %>% 
    dplyr::group_by(cluster) %>%
    dplyr::mutate(direction = ifelse(row_number() <= n.sitc.sigma, "", ".1")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sitc_dir = paste(sitc, direction, sep = "")) %>%
    dplyr::select(-direction)
  
  
  ##############################################################################
  ## merge estimation data
  ##############################################################################
  est.df <- dplyr::left_join(q.df, mu.df, by = c("sitc_dir", "cluster"))
  est.df <- dplyr::left_join(est.df, sigma.df, by = c("sitc_dir", "cluster"))
  
  est.df <- est.df %>%
    dplyr::select(sitc_dir, sitc, cluster, q_kz, mu_kz, sigma_kz) %>%
    dplyr::mutate(sitc = as.character(sitc))
  
  #head(est.df)
  
  
  ##############################################################################
  ## compute hold-out sample log-likelihood
  ##############################################################################
  if(test.run){
    
    test.dyad <- test.df %>% 
      dplyr::pull(dyad_id) %>%
      unique()
      
    test.dyad <- test.dyad[1:2]
    
    obs.df <- test.df %>%
      dplyr::filter(dyad_id %in% test.dyad) %>%
      dplyr::arrange(dyad_id, importer_ISO, exporter_ISO)
    
    # obs.df <- test.df %>% 
    #   dplyr::filter(dyad_id == "100_104" | dyad_id == "156_840") %>%
    #   dplyr::arrange(dyad_id, importer_ISO, exporter_ISO)
    
  } else{
    
    obs.df <- test.df
    
  }
  
  out.z1 <- purrr::map_df(cluster.vec, function(z1) {
    
    out.z2 <- purrr::map_df(cluster.vec, function(z2) {
      
      cat("Computing estimates for transitions from cluster ", z1, " to ", z2, "\n")
      
      # extract pie_zt-1 (scalar)
      pie.out <- pie.df %>% 
        dplyr::filter(cluster == z1) %>%
        dplyr::pull(pie_z)
      
      # extract transition probability (scalar)
      p.mat.out <- p.matrix.collapse[z1 + 1, z2 + 1]
      
      # merge observed data from hold-out period and estimates from fitted model (1250 rows) for specific cluster
      merge.df.z2 <- dplyr::left_join(obs.df,
                                      est.df %>% dplyr::filter(cluster == z2) %>% select(-sitc), 
                                      by = c("sitc_dir"))
      
      # compute and create vars
      merge.df.z2 <- merge.df.z2 %>%
        dplyr::mutate(norm_density = stats::dnorm(W_ijtk, mu_kz, sigma_kz)^(1 - D_ijtk),
                      prod = (q_kz^D_ijtk)*((1 - q_kz)^(1 - D_ijtk))*norm_density) %>%
        dplyr::group_by(dyad_id) %>%
        dplyr::summarize(prod_2k = prod(prod, na.rm = TRUE)) %>%
        dplyr::mutate(all_prod = pie.out*p.mat.out*prod_2k,
                      cluster_z_T_ij = z2)
      
    })
    
    out.z2 <- out.z2 %>%
      dplyr::mutate(cluster_z_T_ij_m1 = z1) %>%
      dplyr::select(dyad_id, cluster_z_T_ij_m1, cluster_z_T_ij, all_prod) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(dyad_id)
  })
  
  ll.out <- out.z1 %>%
    dplyr::group_by(dyad_id) %>%
    dplyr::summarize(sum_trans = sum(all_prod),
                     log_lik_dyad = log(sum_trans))
  
  list(log.lik = sum(ll.out$log_lik_dyad[!is.infinite(ll.out$log_lik_dyad)], 
                     na.rm = TRUE),
       ll.df = ll.out)
  
}
