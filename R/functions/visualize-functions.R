################################################################################
# functions to calculate IIT
################################################################################
calIIT <- function(imp, exp){
  if(imp + exp == 0){
    return(0)
  } else {
    return(1-(abs(imp - exp)/(imp + exp)))
  }
}

ClusterIIT <- function(TF_row, uniq_product_name, weights = NULL){
  ## if(length(TF_row)/2 != length(uniq_product_name)){
  ##     stop("length of TF should be 2 x length of unique products")
  ## }
  
  pnames <- gsub("\\.1$", "", names(TF_row))
  uniq.pnames <- unique(pnames)
  iit <- rep(NA, length(uniq.pnames))
  if(is.null(weights)){ # equal weights
    weights <- rep(1, length(uniq.pnames))
  }
  
  for(i in 1:length(uniq.pnames)){
    ## print(i)
    prod_i <- uniq.pnames[i]
    idx <- which(pnames == prod_i)
    iit[i] <- calIIT(as.numeric(TF_row[idx[1]]),
                     as.numeric(TF_row[idx[2]]))
    if(as.numeric(TF_row[idx[1]]) + as.numeric(TF_row[idx[2]]) == 0){
      ## remove zero trade products
      weights[i] <- 0 
    }
    
    ## weights[i] <- as.numeric(TF_row[idx[1]])+as.numeric(TF_row[idx[2]])
    
  }
  
  
  return(weighted.mean(iit, weights))
}


################################################################################
# cluster type scatterplot
################################################################################
cluster_type_scatterplot <- function(dyad.full.df,
                                     cluster.proportion.df, 
                                     product.description.df,
                                     dyad.trade.1d.df,
                                     palette, 
                                     sitc.1d.labels,
                                     base.size = 17,
                                     xlab.title,
                                     ylab.title,
                                     axis.limits = c(0, 1),
                                     circle.alpha = 0.9,
                                     circle.scale = c(1, 18)){
  
  
  # extract and normalize aggregated tf values (sitc 1st digit)
  p.1d.bar.all <- product.description.df %>%
    group_by(cluster, cty, sitc_num_1d) %>% # cty indicates direction of trade
    slice(1) %>%
    group_by(cluster, cty) %>%
    mutate(value_1d_tot = sum(value_1d),
           value_1d_norm = value_1d/value_1d_tot) %>%
    ungroup()
  
  # calculate IIT info for each cluster using cluster trade proportion file
  iit.df <- data.frame(cluster = seq(1, nrow(cluster.proportion.df), by = 1),
                       iit = apply(cluster.proportion.df, 1, ClusterIIT, 
                                   unique(names(cluster.proportion.df))))
  
  # merge with product composition data
  p.1d.bar.all <- left_join(p.1d.bar.all, iit.df, 
                            by = "cluster")
  
  # create facet variable for plot
  p.1d.bar.all <- p.1d.bar.all %>%
    mutate(facet_id = paste("Cluster ", cluster, 
                            ", IIT Index: ", round(iit, 2), 
                            sep = ""),
           # convert to factor
           sitc_num_1d = factor(sitc_num_1d, 
                                labels = sitc.1d.labels),
           # clean names
           cty = str_replace_all(cty, "first", "First"),
           cty = str_replace_all(cty, "second", "Second"))
  
  # get cluster membership for dyad in each year
  dyad.df.temp <- dyad.full.df %>%
    select(imp, exp, cluster, year) %>%
    mutate(dyad_id_n = ifelse(as.numeric(imp) < as.numeric(exp), 
                              paste(imp, exp, sep = "_"), 
                              paste(exp, imp, sep = "_"))) %>%
    select(-imp, -exp)
  
  # clean sitc 1-digit raw trade data
  dyad.trade.collapse <- dyad.trade.1d.df %>%
    select(-imp, -exp) %>% 
    mutate(year = as.numeric(year))
  
  # merge cluster membership with sitc 1d data
  dyad.trade.sitc1d <- left_join(dyad.df.temp, dyad.trade.collapse,
                                 by = c("dyad_id_n", "year"))
  
  # calculate average trade volume for dyads in same cluster
  dyad.trade.sitc1d.mean <- dyad.trade.sitc1d %>%
    group_by(cluster, year) %>%
    summarise_at(vars(contains("SITC0")), mean, na.rm = TRUE) %>%
    group_by(cluster) %>%
    summarise_at(vars(contains("SITC0")), mean, na.rm = TRUE) %>%
    ungroup()
  
  # reshape to long format
  size.wt <- dyad.trade.sitc1d.mean %>%
    gather(sitc1d, t_volume, -cluster) 
  
  # reshape proportion dataframe
  p.1d.scatter.df <- p.1d.bar.all %>%
    mutate(sitc1d = substr(sitc, 1, 7)) %>%
    select(cluster, cty, sitc1d, sitc_num_1d, value_1d_norm) %>%
    spread(cty, value_1d_norm) %>%
    mutate(cluster_name = paste("Cluster ", cluster, sep = ""))
  
  # merge data with trade proportion and volume
  p.1d.scatter.df <- left_join(p.1d.scatter.df, size.wt, 
                               by = c("cluster", "sitc1d"))
  
  # merge with IIT data
  p.1d.scatter.df <- left_join(p.1d.scatter.df, iit.df, 
                               by = "cluster")
  
  # reorder clusters so results are comparable across output versions
  p.1d.scatter.df <- p.1d.scatter.df %>%
    mutate(facet_id = paste("atop('", cluster_name,
                            "', atop('(IIT Index: ", round(iit, digits = 2), ")'))",
                            sep = ""),
           first_flip = First,
           second_flip = Second)
  
  # create weights by cluster
  p.1d.scatter.df <- p.1d.scatter.df %>%
    group_by(cluster) %>%
    mutate(max_t_volume_c = max(t_volume),
           t_volume_wt_c = t_volume/max_t_volume_c) %>%
    ungroup() %>%
    mutate(max_t_volume_all = max(t_volume),
           t_volume_wt_all = t_volume/max_t_volume_all) %>%
    filter(!is.na(cluster))
  
  # plot
  p.p.scatter.all.across <- ggplot(p.1d.scatter.df, 
                                   aes(x = first_flip, 
                                       y = second_flip, 
                                       group = as.factor(sitc_num_1d), 
                                       color = as.factor(sitc_num_1d),
                                       #size = t_volume_wt_all,
                                       #size = t_volume
                                       size = sqrt(t_volume))) +
    geom_abline(intercept = 0, slope = 1,
                linetype = "solid", size = 0.3) +
    geom_point(alpha = circle.alpha) + 
    scale_color_manual(values = palette,
                       #name = "SITC 1st Digit Products",
                       guide = guide_legend(override.aes = list(size = 4.5))) +
    scale_size(range = circle.scale, guide = FALSE) + 
    scale_x_continuous(limits = axis.limits) + 
    scale_y_continuous(limits = axis.limits) + 
    facet_grid(. ~ facet_id, labeller = label_parsed) +
    coord_fixed(ratio = 1) +
    xlab(xlab.title) + 
    ylab(ylab.title) + 
    theme_bw() + theme(plot.title = element_text(lineheight = .8, face = "bold", 
                                                 margin = margin(0, 0, 20, 0)),
                       axis.title.x = element_text(size = base.size,
                                                   margin = margin(20, 0, 0, 0)),
                       axis.title.y = element_text(size = base.size,
                                                   margin = margin(0, 20, 0, 0)),
                       axis.text.x = element_text(size = base.size - 3),
                       axis.text.y = element_text(size = base.size - 3),
                       strip.text = element_text(size = base.size + 5),
                       strip.background = element_blank(),
                       panel.border = element_rect(colour = "black"),
                       legend.key.size = unit(0.3, "inch"),
                       legend.key = element_rect(colour = "white"),
                       legend.text = element_text(size = base.size - 3),
                       legend.title = element_blank())
  
  return(p.p.scatter.all.across)
}


################################################################################
# function to create transparency for heatmaps
################################################################################
addTrans <- function(color,trans)
{
  library(RColorBrewer)
  library(plotrix)
  
  # This function adds transparancy to a color.
  # Define transparancy with an integer between 0 and 255
  # 0 being fully transparant and 255 being fully visable
  # Works with either color and trans a vector of equal length,
  # or one of the two of length 1.
  
  if (length(color)!=length(trans)&!any(c(length(color),length(trans))==1)) stop("Vector lengths not correct")
  if (length(color)==1 & length(trans)>1) color <- rep(color,length(trans))
  if (length(trans)==1 & length(color)>1) trans <- rep(trans,length(color))
  
  num2hex <- function(x)
  {
    hex <- unlist(strsplit("0123456789ABCDEF",split=""))
    return(paste(hex[(x-x%%16)/16+1],hex[x%%16+1],sep=""))
  }
  rgb <- rbind(col2rgb(color),trans)
  res <- paste("#",apply(apply(rgb,2,num2hex),2,paste,collapse=""),sep="")
  return(res)
}



################################################################################
## this heatmap function takes a matrix with color indices
################################################################################
myheatmap_matrix <- function(matrix,
                             x.names = NULL,
                             y.names = NULL,
                             col.twotone = TRUE, 
                             twotone.cutoff = 0,
                             col.desc = NULL,
                             main = NULL,
                             colors = NULL,
                             cex.main = 1.5, 
                             cex.x = 1.2, cex.y = .3, cex.legend = 1.5,
                             xlas = 2, ylas = 2,
                             plot.legend = TRUE, legend.labels = NULL, 
                             legend.inset = 0.06,
                             x.txt.by = 5){
  library(RColorBrewer)
  library(plotrix)
  
  if(is.null(x.names)){
    x.names = colnames(matrix)
    x.names = x.names[seq(1, length(x.names), by = x.txt.by)]
  }
  if(is.null(y.names)){
    y.names = rownames(matrix)
  }
  if(is.null(main)){
    main <- ""
  }
  
  n.y <- nrow(matrix)
  n.x <- ncol(matrix)
  
  if(is.null(colors)){
    stop("colors must be provided")
  }
  if(plot.legend == TRUE){
    par(mar = c(3, 10, 3, 0), cex.main = cex.main, xpd = TRUE) # set margins
  } else {
    par(mar = c(5, 12, 1, 1), cex.main = cex.main)
  }
  ### draw plot area
  plot(1, 1,
       if(plot.legend ==TRUE){
         #xlim = c(.5, n.x*1.25)
         xlim = c(.5, n.x + .5)
       } else {
         xlim = c(.5, n.x + .5)
       },
       #ylim = c(.5, n.y + .5),
       ylim = c(.5, n.y + 10),
       main = main,
       ## main = "HeatMap of Intra-industry Trade",
       xlab = "",
       ylab = "",
       type = "n", axes = FALSE)
  
  ## labels on y-axis
  mtext(y.names, side = 2, at = 1:n.y, las = ylas, line = -2,
        cex = cex.y)
  ## labels on x-axis
  mtext(x.names, side = 1, at = seq(1, length(colnames(matrix)), x.txt.by), 
        las = xlas, line = -4.5,
        cex = cex.x)
  
  ## draw rectangles: loop across products and across variables
  for (i in 1:n.y) {
    for (t in 1:n.x) {
      if(is.na(matrix[i,t])){
        polygon.col <- "lightgrey"
        #polygon.col <- "white"
      } else{
        color.type <- matrix[i,t]
        polygon.col <- colors[color.type]
      }
      #       polygon(c(.5 + t - 1, .5 + t, .5 + t, .5 + t - 1),
      #               c(.5 + i - 1, .5 + i - 1, .5 + i, .5 + i),
      #               # vertices of rectangle
      #               density = NA,
      #               # filled polygons
      #               border = "white",
      #               lwd = 0.1,
      #               #border = NA,
      #               # white border of rectangle
      #               col = polygon.col)
      rect(.5 + t - 1, 
           .5 + i - 1, 
           .5 + t,
           .5 + i,
           # vertices of rectangle
           density = NA,
           # filled polygons
           #border = "white",
           #lwd = 0.02,
           border = NA,
           # white border of rectangle
           col = polygon.col)
      # specify color depending on matrix
    }
  }
  
  if(plot.legend == TRUE){
    if(is.null(col.desc)){
      string <- seq(1, length(colors))
    } else {
      string <- col.desc
    }
    ## print(string)
    ## print(polygon.cols)
    # legend(n.x+.60, n.y,
    #        as.character(string),
    #        col=colors,
    #        lty=rep(1,length(string)), lwd=rep(5, length(string)),
    #        bty="n",
    #        cex=1)
    legend("top",
           as.character(string),
           col = colors,
           horiz = TRUE,
           lty = rep(1,length(string)), lwd = rep(5, length(string)),
           bty = "n",
           inset = legend.inset,
           cex = cex.legend)
  }
}



################################################################################
# function to create heatmap for products
################################################################################
HeatmapRedBlue <- function(matrix,
                           x.names = colnames(matrix), 
                           y.names = rownames(matrix),
                           main = "Title",
                           cex.main = 1, cex.x = 1.2, cex.y = .2, 
                           legend.cex = 1.4,
                           xlas = 2, ylas = 2,
                           xline = -1,
                           plot.legend = TRUE, legend.labels = NULL,
                           blue.cut = NULL, red.cut = NULL,
                           plot.mar = c(7, 8, 4, 3)){
  library(RColorBrewer)
  library(plotrix)
  
  ## This function takes a matrix with continuous values, and uses
  ## Red for positive and Blue for negative values: It makes values
  ## closer to zero lighter using addTrans function
  
  n.y <- nrow(matrix)
  n.x <- ncol(matrix)
  
  ## determining polygon.cols as a matrix with same dimension as matrix
  polygon.cols <- matrix(NA, nrow = n.y, ncol = n.x)
  if(is.null(blue.cut)){
    blue.cut <- as.numeric(quantile(matrix, probs = seq(0,1,.01), 
                                    na.rm = TRUE)[2])
    min.value <- as.numeric(min(c(matrix)))
  } else {
    min.value <- blue.cut
  }
  if(is.null(red.cut)){
    red.cut <- as.numeric(quantile(matrix, probs = seq(0,1,.01), 
                                   na.rm = TRUE)[100])
    max.value <- as.numeric(max(c(matrix)))
  } else {
    max.value <- red.cut
  }
  
  ## others to have values proportional to the cut values
  for(i in 1:n.y){
    for (j in 1:n.x){
      value <- matrix[i,j]
      if(!exists("value")){
        value <- 0
      }
      if(value>0){
        polygon.cols[i,j] <- addTrans("red", 255 * abs(value/red.cut))
      } else if (value < 0){
        polygon.cols[i,j] <- addTrans("blue", 255 * abs(value/blue.cut))
      } else {
        polygon.cols[i,j] <- addTrans("white", 255)
      }
    }
  }
  ## setting extreme values to have darkest 
  polygon.cols[matrix<blue.cut] <- addTrans("blue", 255)
  polygon.cols[matrix>red.cut] <- addTrans("red", 255)
  
  if(plot.legend == TRUE){
    par(mar = plot.mar - c(0, 0, 0, 2), cex.main = cex.main) # set margins
    ## par(mar = c(2, 2, 2, 2), cex.main=cex.main) # set margins
  } else {
    par(mar = plot.mar, cex.main = cex.main)
  }
  ### draw plot area
  if(n.x < 3){
    plot(1, 1,
         if(plot.legend == TRUE){
           xlim = c(0.5, n.x + .7)
         } else {
           xlim = c(0.5, n.x + 0.5)
         },
         ylim = c(0.5, n.y + 0.5),
         main = main,
         ## main = "HeatMap of Intra-industry Trade",
         xlab = "",
         ylab = "",
         type = "n", axes = FALSE)
  } else {
    plot(1, 1,
         if(plot.legend ==TRUE){
           #xlim = c(0.5, n.x*1.25)
           xlim = c(0.5, n.x*1.35)
         } else {
           xlim = c(0.5, n.x + 0.5)
         },
         ylim = c(0.5, n.y + 0.5),
         main = main,
         ## main = "HeatMap of Intra-industry Trade",
         xlab = "",
         ylab = "",
         type = "n", axes = FALSE)
    
  }
  ## labels on y-axis
  mtext(y.names, side = 2, at = 1:n.y, las = ylas, line = -1,
        cex = cex.y)
  ## labels on x-axis
  mtext(x.names, side = 1, at = 1:n.x, las = xlas, line = xline,
        cex = cex.x)
  
  ## draw rectangles: loop across products and across variables
  for (i in 1:n.y) {
    for (t in 1:n.x) {
      if(is.na(matrix[i,t])){
        matrix[i,t] <- 0
      }
      
      polygon(c(0.5 + t - 1, 0.5 + t, 
                0.5 + t, 0.5 + t - 1),
              c(0.5 + i - 1, 0.5 + i - 1, 
                0.5 + i, 0.5 + i),
              # vertices of rectangle
              density = NA,
              # filled polygons
              ## border = "white",
              border = NA,
              # white border of rectangle
              col = polygon.cols[i,t])
      # specify color depending on matrix
    }
  }
  
  if(plot.legend == TRUE){
    if(is.null(legend.labels)){
      legend.labels <- c(round(min.value,2), round(max.value,2))
    }
    if(n.x < 3){
      color.legend(n.x+.65, n.y*.4, n.x+.7, n.y*.9,  legend = legend.labels,
                   c(polygon.cols), gradient = "y", align = "rb",
                   cex = legend.cex)
    } else {
      
      reds <- colorRampPalette(c("white", "red"))
      reds.legend <- reds(length(which(matrix>0)))
      print(length(which(matrix>0)))
      reds.legend <- reds(round(1000*(max.value/abs(max.value-min.value))))
      
      blues <- colorRampPalette(c("blue", "white"))
      ## blues.legend <- blues(round(100*(min.value/abs(max.value-min.value))))
      blues.legend <- blues(round(1000*(abs(min.value)/abs(max.value-min.value))))
      ## blues.legend <- blues(length(which(matrix<0)))
      legend.cols <- c(blues.legend, reds.legend)
      
      if(n.x == 3){
        color.legend(#n.x*1.18, n.y*.38, n.x*1.2, n.y*.86,  
          n.x*1.24, n.y*.38, n.x*1.26, n.y*.86,  
          legend = legend.labels,
          rect.col = legend.cols,
          gradient = "y", align = "rb",
          cex = legend.cex)
      } else {
        color.legend(#n.x*1.18, n.y*.38, n.x*1.2, n.y*.86,  
          n.x*1.18, n.y*.38, n.x*1.2, n.y*.86,  
          legend = legend.labels,
          rect.col = legend.cols,
          gradient = "y", align = "rb",
          cex = legend.cex) 
      }
    }
  }
}


################################################################################
# product heatmap
################################################################################
product_heatmap <- function(data,
                            cluster.order = seq(1, n.cluster*2, 1),
                            filename,
                            product.description,
                            product.txt.cex = 1.2,
                            industry.label.loc = NULL,
                            sitc.1d.dash.shift = 0.1,
                            cluster.label.hjust = 1.1,
                            cluster.label.cex = 2,
                            pdf.width = 10, 
                            pdf.height = 10,
                            legend.cex = 1.75,
                            legend.label.high.txt = "Higher Trade\nProportion",
                            legend.label.low.txt = "Lower Trade\nProportion",
                            legend.label.high.y = 400,
                            legend.label.low.y = 155){
  # add cluster number to df
  data.cl <- data %>%
    mutate(cluster = seq(1, nrow(data), by = 1))
  
  prod.order <- names(data)[1:(ncol(data)/2)]
    
  # reshape
  tf <- data.cl %>%
    gather(sitc, value, -cluster) %>%
    mutate(direction = ifelse(grepl("\\.1", sitc), "second", "first"),
           sitc = str_replace_all(sitc, "\\.1", "")) %>%
    group_by(sitc) %>%
    mutate(mean = mean(value),
           demean = value - mean) %>%
    ungroup() %>%
    select(cluster, sitc, direction, demean) %>%
    mutate(cluster_dir = paste("Cluster", cluster, direction, sep = "_")) %>%
    select(-cluster, -direction) %>%
    mutate(sitc = factor(sitc, levels = prod.order)) %>%
    spread(cluster_dir, demean) 
  
  # flip directions to be consistent across versions
  tf <- tf[, c(1, 1 + cluster.order)]
  
  # convert to matrix for plotting
  tf.dev <- tf %>%
    select(-sitc)
  
  tf.dev <- as.matrix(tf.dev)
  
  # set axis labels
  y.names <- tf$sitc
  y.names <- str_replace_all(y.names, "SITC0_", "")
  x.names <- c()
  
  # plot
  pdf(file = filename,
      width = pdf.width, 
      height = pdf.height)  

  # plot main heatmap
  HeatmapRedBlue(tf.dev , x.names = x.names, 
                 #y.names=y.names
                 y.names = NULL,
                 ## main="Deviation from Mean Level of Trade",
                 main = "",
                 ## legend.labels=c("", ""),
                 xlas = 1, cex.main = 2, xline = 2, cex.x = 1.8#,
                 #blue.cut = -0.01192352, red.cut = 0.07726225
  ) 
  
  if(is.null(industry.label.loc)){
    
    # set product description location
    products <- tf$sitc
    industry.cut <- grep("^SITC0_1", products)[1]
    industry.cut <- c(industry.cut, grep("^SITC0_2", products)[1])
    industry.cut <- c(industry.cut, grep("^SITC0_3", products)[1])
    industry.cut <- c(industry.cut, grep("^SITC0_4", products)[1])
    industry.cut <- c(industry.cut, grep("^SITC0_5", products)[1])
    industry.cut <- c(industry.cut, grep("^SITC0_6", products)[1])
    industry.cut <- c(industry.cut, grep("^SITC0_7", products)[1])
    industry.cut <- c(industry.cut, grep("^SITC0_8", products)[1])
    industry.cut <- c(industry.cut, grep("^SITC0_9", products)[1])
    
    # plot horizontal dash separater and product label
    #sitc.1d.dash.shift <- 0.1
    
    for(k in 1:length(industry.cut)){
      # draw dash lines
      lines(c(-2, (ncol(tf.dev) + .5)), 
            c(industry.cut[k] + sitc.1d.dash.shift, 
              industry.cut[k] + sitc.1d.dash.shift),
            lty = 2, lwd = 1.5, col = "black")
      if(k==1){
        loc <- (industry.cut[k] + sitc.1d.dash.shift)/2
      } else {
        loc <- industry.cut[k-1] + ((industry.cut[k]-industry.cut[k-1])/2)
      }
      # draw product description
      mtext(descr[k],
            side = 2, at = loc, las = 2, line = 7, adj = 0, cex = product.txt.cex)
    }
  } else{
    loc <- industry.label.loc
    
    for(k in 1:length(industry.label.loc)){
      lines(c(-2, (ncol(tf.dev) + .5)), 
            c(industry.label.loc[k] + sitc.1d.dash.shift, 
              industry.label.loc[k] + sitc.1d.dash.shift),
            lty = 2, lwd = 1.5, col = "black")
      # if(k==1){
      #   loc <- (industry.label.loc[k] + sitc.1d.dash.shift)/2
      # } else {
      #   loc <- industry.label.loc[k-1] + ((industry.cut[k]-industry.cut[k-1])/2)
      # }
      # draw product description
      mtext(descr[k],
            side = 2, at = loc[k], las = 2, line = 7, adj = 0, cex = product.txt.cex)
    }
  }
  
  # plot x labels
  x.names <- rep(1:nrow(data))
  x.names <- paste("Cluster ", x.names, sep = "")

  for(i in 1:(ncol(tf.dev)/2)){
     loc <- 2*(i-1) + cluster.label.hjust
    
     # add cluster labels
     mtext(x.names[i],
           side = 1, at = loc, line = .1, adj = 0, cex = cluster.label.cex)
      
     # add vertical white line to separate clusters
     abline(v = i*2 + .5, lwd = 5, col = "white")

    }

  # add legend labels
  text((ncol(tf.dev)*1.26), legend.label.high.y, legend.label.high.txt, 
       cex = legend.cex, offset = 0)
  text((ncol(tf.dev)*1.26), legend.label.low.y, legend.label.low.txt, 
       cex = legend.cex, offset = 0)
  
  dev.off()

}


################################################################################
# trend plot by proportion of dyads
################################################################################
cluster_trend_ndyads <- function(data, 
                                 palette = NULL,
                                 linetype = NULL,
                                 base.txt.size = 16,
                                 xlab.title,
                                 ylab.title,
                                 label.v.shift.vec = NULL,
                                 label.h.shift.vec = NULL){
  
  # calculate number of dyads per cluster-year and proportions
  df.out <- data %>%
    group_by(year) %>%
    mutate(n_dyads_year = length(dyad_id_n)) %>%
    group_by(cluster, year) %>%
    summarise(n_dyads = length(cluster),
              n_dyads_year = head(n_dyads_year, 1)) %>%
    ungroup() %>%
    mutate(per_dyads = n_dyads/n_dyads_year, 
           cluster_name = paste("Cluster ", cluster, sep = ""))
  
  # extract data for annotations
  txt.df <- df.out %>%
    filter(year == max(df.out$year)) %>% 
    select(cluster, year, per_dyads, cluster_name) %>%
    mutate(per_dyads = per_dyads + label.v.shift.vec,
           year = year + label.h.shift.vec) %>%
    arrange(desc(per_dyads))
  
  # plot
  p.ndyad.year <- ggplot(df.out, 
                         aes(x = year, 
                             y = per_dyads, 
                             group = as.factor(cluster),
                             color = as.factor(cluster))) +
    geom_text(data = txt.df,
              label = txt.df$cluster_name,
              size = base.txt.size - 11.5) +
    scale_x_continuous(xlab.title,
                       breaks = seq(min(df.out$year), 
                                    max(df.out$year), 
                                    by = 10), 
                       limits = c(min(df.out$year) - 1, 
                                  max(df.out$year) + 1)) + 
    scale_y_continuous(ylab.title,
                       breaks = seq(0, 1, 
                                    by = 0.25),
                       limits = c(0, 1)) + 
    theme_bw() + theme(plot.title = element_text(lineheight = .8, 
                                                 face = "bold", 
                                                 vjust = 1.5),
                       axis.title.x = element_text(size = base.txt.size,
                                                   margin = margin(20, 0, 0, 0)),
                       axis.title.y = element_text(size = base.txt.size, 
                                                   margin = margin(0, 10, 0, 20)),
                       axis.text.x = element_text(size = base.txt.size - 2),
                       axis.text.y = element_text(size = base.txt.size - 2),
                       #panel.grid.minor = element_blank(),
                       legend.position = "none") 
  
  if(is.null(palette) & is.null(linetype)){
    p.ndyad.year <- p.ndyad.year + 
      geom_line() 
    
  } else if(!is.null(palette) & is.null(linetype)){
    p.ndyad.year <- p.ndyad.year + 
      geom_line() + 
      scale_colour_manual(values = palette) 
    
  } else if(is.null(palette) & !is.null(linetype)){
    
    p.ndyad.year <- p.ndyad.year + 
      geom_line(aes(linetype = as.factor(cluster))) + 
      scale_linetype_manual(values = linetype) 
    
  } else if(!is.null(palette) & !is.null(linetype)) {
    p.ndyad.year <- p.ndyad.year + 
      geom_line(aes(linetype = as.factor(cluster))) + 
      scale_colour_manual(values = palette) +
      scale_linetype_manual(values = linetype) 
  }
  
  return(p.ndyad.year)
}


################################################################################
# trend plot by proportion of world trade
################################################################################
cluster_trend_tvolume <- function(data, 
                                  palette = NULL,
                                  linetype = NULL,
                                  base.txt.size = 16,
                                  xlab.title,
                                  ylab.title,
                                  label.v.shift.vec = NULL,
                                  label.h.shift.vec = NULL){
  
  # calculate cluster-year proportion of trade to World
  df.out <- data %>%
    group_by(year) %>%
    mutate(world_trade_year = sum(value)) %>%
    group_by(cluster, year) %>%
    mutate(tvolume_cy = sum(value),
           total_volume = tvolume_cy/world_trade_year) %>%
    slice(1) %>%
    ungroup() %>%
    select(cluster, year, tvolume_cy, world_trade_year, total_volume) %>%
    mutate(year = as.numeric(year),
           cluster_name = paste("Cluster ", cluster, sep = ""))
  
  # extract data for annotations
  txt.df <- df.out %>%
    filter(year == min(df.out$year)) %>% 
    select(cluster, year, total_volume, cluster_name) %>%
    mutate(total_volume = total_volume + label.v.shift.vec,
           year = year + label.h.shift.vec) %>%
    arrange(desc(total_volume))
  
  # plot
  p.tvolume.year <- ggplot(df.out, 
                           aes(x = year, 
                               y = total_volume, 
                               group = as.factor(cluster),
                               color = as.factor(cluster))) +
    geom_text(data = txt.df,
              label = txt.df$cluster_name,
              size = base.txt.size - 11.5) +
    scale_x_continuous(xlab.title, 
                       breaks = seq(min(df.out$year), 
                                    max(df.out$year), 
                                    by = 10), 
                       limits = c(min(df.out$year) - 1, 
                                  max(df.out$year) + 1)) + 
    scale_y_continuous(ylab.title,
                       breaks = seq(0, 1, 
                                    by = 0.25),
                       limits = c(0, 1)) + 
    theme_bw() + theme(plot.title = element_text(lineheight = .8, 
                                                 face = "bold", 
                                                 vjust = 1.5),
                       axis.title.x = element_text(size = base.txt.size,
                                                   margin = margin(20, 0, 0, 0)),
                       axis.title.y = element_text(size = base.txt.size, 
                                                   margin = margin(0, 10, 0, 20)),
                       axis.text.x = element_text(size = base.txt.size - 2),
                       axis.text.y = element_text(size = base.txt.size - 2),
                       #panel.grid.minor = element_blank(),
                       legend.position = "none") 
  
  if(is.null(palette) & is.null(linetype)){
    p.tvolume.year <- p.tvolume.year + 
      geom_line() 
    
  } else if(!is.null(palette) & is.null(linetype)){
    p.tvolume.year <- p.tvolume.year + 
      geom_line() + 
      scale_colour_manual(values = palette) 
    
  } else if(is.null(palette) & !is.null(linetype)){
    
    p.tvolume.year <- p.tvolume.year + 
      geom_line(aes(linetype = as.factor(cluster))) + 
      scale_linetype_manual(values = linetype) 
    
  } else if(!is.null(palette) & !is.null(linetype)) {
    p.tvolume.year <- p.tvolume.year + 
      geom_line(aes(linetype = as.factor(cluster))) + 
      scale_colour_manual(values = palette) +
      scale_linetype_manual(values = linetype) 
  }
  
  return(p.tvolume.year)
}


################################################################################
# partner heatmap
################################################################################
partner_heatmap <- function(data,
                            cty.vec,
                            palette,
                            plot.legend = TRUE,
                            title.cex = 3.5,
                            title.hjust = 1.03,
                            legend.cex = 3,
                            legend.inset = 0.03,
                            x.txt.by = 1,
                            file.path = NULL,
                            pdf.name = "heatmap_",
                            pdf.width = 15,
                            pdf.height = 40){
  
  # subset vars, and convert clusters
  member.df <- data %>%
    select(dyad_id_n, iso3n_imp, iso3n_exp, cty_imp, cty_exp, year, cluster) 

  # plot
  walk(cty, function(x){
    # set palette
    heatmap.palette <- palette
    heatmap.cty <- x
    heatmap.desc <- paste("Cluster ", 
                          seq(1, max(data$cluster), by = 1), 
                          sep = "")
    
    cat("Plotting heatmap for", heatmap.cty, "\n")
    
    iso3n <- countrycode(heatmap.cty, "country.name", "iso3n")
    
    # reshape dyad-year long to wide format
    member.df <- member.df %>%
      filter(iso3n_imp == iso3n | iso3n_exp == iso3n) %>%
      mutate(partner = ifelse(iso3n_imp != iso3n, iso3n_imp, iso3n_exp),
             partner_name = ifelse(iso3n_imp != iso3n, cty_imp, cty_exp)) %>%
      distinct(dyad_id_n, year, .keep_all = TRUE) 
    
    # calculate order of partners based on number of partners in cluster 1, 2, etc.
    if(n.cluster == 3) {
      count <- member.df %>%
        count(cluster, partner_name) %>%
        arrange(cluster, desc(n)) %>%
        spread(cluster, n) %>% 
        rename(n_c1 = `1`,
               n_c2 = `2`,
               n_c3 = `3`) %>%
        ungroup() %>%
        arrange(desc(n_c3), desc(n_c2), desc(n_c1)) 
      
    } else if(n.cluster == 7){
      count <- member.df %>% 
        count(cluster, partner_name) %>%
        arrange(cluster, desc(n)) %>%
        spread(cluster, n)
      
      if(!(7 %in% names(count))){
        count$`7` <- NA
      }
      
      count <- count %>%
        rename(n_c1 = `1`,
               n_c2 = `2`,
               n_c3 = `3`,
               n_c4 = `4`,
               n_c5 = `5`,
               n_c6 = `6`,
               n_c7 = `7`) %>%
        arrange(desc(n_c7), desc(n_c6), desc(n_c5),
                desc(n_c4), desc(n_c3), desc(n_c2), desc(n_c1))
      
    } else if(n.cluster == 15){
      count <- member.df %>% 
        count(cluster, partner_name) %>%
        arrange(cluster, desc(n)) %>%
        spread(cluster, n)
      
      if(!(1 %in% names(count))){
        count$`1` <- NA
      } 
      if (!(2 %in% names(count))){
        count$`2` <- NA
      }
      if (!(3 %in% names(count))){
        count$`3` <- NA
      }
      if (!(4 %in% names(count))){
        count$`4` <- NA
      }
      if (!(5 %in% names(count))){
        count$`5` <- NA
      }
      if (!(6 %in% names(count))){
        count$`6` <- NA
      }
      if (!(7 %in% names(count))){
        count$`7` <- NA
      }
      if (!(8 %in% names(count))){
        count$`8` <- NA
      }
      if (!(9 %in% names(count))){
        count$`9` <- NA
      }
      if (!(10 %in% names(count))){
        count$`10` <- NA
      }
      if (!(11 %in% names(count))){
        count$`11` <- NA
      }
      if (!(12 %in% names(count))){
        count$`12` <- NA
      }
      if (!(13 %in% names(count))){
        count$`13` <- NA
      }
      if (!(14 %in% names(count))){
        count$`14` <- NA
      }
      if (!(15 %in% names(count))){
        count$`15` <- NA
      }
      
      count <- count %>%
        rename(n_c1 = `1`,
               n_c2 = `2`,
               n_c3 = `3`,
               n_c4 = `4`,
               n_c5 = `5`,
               n_c6 = `6`,
               n_c7 = `7`,
               n_c8 = `8`,
               n_c9 = `9`,
               n_c10 = `10`,
               n_c11 = `11`,
               n_c12 = `12`,
               n_c13 = `13`,
               n_c14 = `14`,
               n_c15 = `15`) %>%
        arrange(desc(n_c15), desc(n_c14), desc(n_c13), desc(n_c12), desc(n_c11),
                desc(n_c10), desc(n_c9), desc(n_c8), desc(n_c7), desc(n_c6), 
                desc(n_c5), desc(n_c4), desc(n_c3), desc(n_c2), desc(n_c1))
      
    } else{
      stop("Model not supported")
    }
    
    # reshape to wide format
    member.df.wide <- member.df %>%
      select(-iso3n_imp, -iso3n_exp, -cty_imp, -cty_exp, -partner, -dyad_id_n) %>%
      spread(year, cluster) 
    
    # reorder rows based on count
    member.df.wide <- member.df.wide[match(rev(count$partner_name), 
                                           member.df.wide$partner_name),]
    
    
    # shorten full names for plotting
    member.df.wide$partner_name <- str_replace_all(member.df.wide$partner_name, 
                                                   "South Georgia and the South Sandwich Islands", 
                                                   "South Georgia and Sandwich Islands")
    member.df.wide$partner_name <- str_replace_all(member.df.wide$partner_name,
                                                   "United States Minor Outlying Islands", 
                                                   "US Minor Outlying Islands")
    member.df.wide$partner_name <- str_replace_all(member.df.wide$partner_name,
                                                   "Saint Vincent and the Grenadines", 
                                                   "St. Vincent and the Grenadines")
    row.names(member.df.wide) <- member.df.wide$partner_name
    
    # drop vars
    member.df.wide <- member.df.wide %>% 
      select(-partner_name)
    
    # convert to matrix
    member.matrix <- as.matrix(member.df.wide)
    
    # drop partners with any missing
    member.matrix <- na.omit(member.matrix)
    
    # plot heatmap
    pdf(file = paste(file.path, pdf.name,
                     heatmap.cty, ".pdf", sep = ""), 
        width = pdf.width, height = pdf.height)
    
    myheatmap_matrix(matrix = member.matrix,
                     colors = heatmap.palette,
                     col.desc = heatmap.desc,
                     main = "",
                     cex.main = 1.5, cex.x = 2, cex.y = 0.9, 
                     plot.legend = plot.legend,
                     xlas = 2, ylas = 2,
                     cex.legend = legend.cex,
                     legend.inset = legend.inset,
                     x.txt.by = x.txt.by)
    
    ## text(ncol(heatmap.df())/2, nrow(heatmap.df())*1.12,
    text((ncol(member.matrix) + 1)/2, 
         #nrow(member.matrix)*1.08,
         #paste("Cluster Membership: \n Partners of ", toupper(heatmap.cty), sep=""),
         nrow(member.matrix)*title.hjust,
         toupper(heatmap.cty),
         cex = title.cex)
    dev.off()
  })
}


