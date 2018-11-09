#' A function that computes D_ijtk and W_ijtk using the held-out trade data
#'
#' @param file.dir The file directory of the held-out data.
#' @return A dataframe that contains D_ijtk and W_ijtk in long form.
#' @import dplyr 
#' @import purrr
#' @import tidyr 
#' @import utils
#' @export
#' @examples
#' file.dir <- "./example/toy/dynamic_1970.csv"
#' df <- computeDsWs(file.dir)
computeDsWs <- function(file.dir){
  
  ##############################################################################
  ## read in raw trade data X_ijtk
  ##############################################################################
  cat("Reading in raw data\n")
  
  # read sample data
  sampleData <- utils::read.table(file.dir, 
                                  header = TRUE, 
                                  check.names = FALSE, nrow = 1)
  ncols <- ncol(sampleData)
  
  # read in all data
  df.r <- utils::read.table(file.dir,
                            header = TRUE,
                            check.names = FALSE,
                            colClasses = c(rep("character", 3),
                                           rep("numeric", (ncols - 3))))
  
  # # check
  # head(df.r)
  # glimpse(df.r)
  
  # change to long form
  df.l <- df.r %>%
    tidyr::gather(sitc, X_ijtk, 
                  -year, -importer_ISO, -exporter_ISO) %>%
    dplyr::arrange(year, importer_ISO, exporter_ISO)
  

  ##############################################################################
  ## create variables of interest
  # V_ijtk (directed-dyad-year-product trade share)
  # V_ijt_last_prod (directed-dyad-year-product trade share of last product)
  # D_ijtk (whether the trade share proportion V_ijtk is zero)
  # V_ijt_last_prod_c (V_ijt_last_prod + c)
  # W_ijtk (log of normalized trade share by last-column product)
  ##############################################################################
  
  # set small constant
  c <- 0.0001
  
  # compute variables
  cat("Computing variables of interest\n")
  
  df.out <- df.l%>%
    dplyr::group_by(importer_ISO, exporter_ISO) %>%
    dplyr::mutate(tot_direct_dyad_year_volume = sum(X_ijtk, na.rm = TRUE),
                  V_ijtk = X_ijtk/tot_direct_dyad_year_volume,
                  sum_direct_dyad_prod_share = sum(V_ijtk, na.rm = TRUE),
                  D_ijtk = if_else(V_ijtk == 0 | (X_ijtk == 0 & tot_direct_dyad_year_volume == 0), 1, 0),
                  V_ijt_last_prod = X_ijtk[sitc == names(df.r)[length(names(df.r))]],
                  V_ijt_last_prod_c = if_else(V_ijt_last_prod == 0, V_ijt_last_prod + c, V_ijt_last_prod), # consistent with C++ code
                  W_ijtk = if_else(V_ijtk == 0, NA_real_, log(V_ijtk/V_ijt_last_prod_c))) %>%
    dplyr::ungroup()
  
  # # check
  # head(df.out)
  # df.out %>% filter(importer_ISO == 840 & exporter_ISO == 156) %>% head(20)
  # df.out %>% filter(importer_ISO == 840 & exporter_ISO == 156) %>% tail(20)
  # df.out %>% filter(importer_ISO == 840 & exporter_ISO == 156 & 
  #                     sitc == names(df.r)[length(names(df.r))])
  # df.out %>% filter(sum_direct_dyad_prod_share != 1)
  # summary(df.out$D_ijtk)
  # summary(df.out$W_ijtk)
  
  # create vars
  df.out <- df.out %>% 
    dplyr::select(importer_ISO, exporter_ISO, sitc, D_ijtk, W_ijtk) %>%
    dplyr::mutate(dyad_id = ifelse(as.numeric(importer_ISO) < as.numeric(exporter_ISO), 
                                   paste(importer_ISO, exporter_ISO, sep = "_"), 
                                   paste(exporter_ISO, importer_ISO, sep = "_")),
                  sitc_dir = ifelse(as.numeric(importer_ISO) < as.numeric(exporter_ISO), 
                                    sitc, 
                                    paste(sitc, ".1", sep = "")),
                  sitc_dir_flip = ifelse(as.numeric(importer_ISO) > as.numeric(exporter_ISO),
                                         sitc,
                                         paste(sitc, ".1", sep = "")))
  
  return(df.out)
}
