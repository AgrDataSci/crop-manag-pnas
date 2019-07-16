####################################################
# Script used to model the performance of crop varieties with climatic variables
# a forward selection with Plackett-Luce trees is performed to select the best
# generalizable model for crop performance forecast
####################################################

library("tidyverse")
library("Matrix")
library("reshape2")
library("magrittr")
library("gosset")
library("caret")
library("psychotools")
library("PlackettLuce")
library("BradleyTerryScalable")
library("partykit")
library("igraph")
library("foreach")
library("doParallel")
library("abind")
library("dismo")
library("rgeos")
library("gstat")
library("raster")
library("alphahull")

# .............................................................
# .............................................................
# Read data ####

# define number of cores for parallelisation
n_cpu <- 3

# read data
df <- read_csv("data/tricot_data.csv", na = c("NA",""))

# crop species that will be analysed 
crop <-  sort(unique(df$crop))


# .............................................................
# .............................................................
# This analysis will evaluate the performance of Plackett-Luce models under following models:
# 1. PLT-null without variables
# 2. PLT-design with variables related to experiental design 
# 3. PLT-climate with climatic variables
# 4. PLT-clim+loc with variables from climate and location
approach <- c("PLT-null","PLT-design","PLT-climate","PLT-clim+loc")

# .............................................................
# .............................................................
# Run models ####
set.seed(123)

for(m in seq_along(crop)){
  
  # .............................................................
  # .............................................................
  # Subset data and check rankings and explanatory variables ####
  
  cat("\n Starting analysis:", toupper(crop[m]), " \n Time:", date(), "\n")
  
  # Subset data for the m crop
  mydata <- df[df$crop == crop[m] , ]
  
  # Create folder for outputs
  output <- paste0("output/", crop[m], "/")
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  
  # reclassify factors levels to fit in the subseted data
  mydata[c("soil","season")] <- lapply(mydata[c("soil","season")], 
                                       function(x) as.factor(as.character(x)))
  
  # check consistency of data and remove inconsistent rows
  if(crop[m] == "commonbean"){
    #remove NA's in rankings for overall_performance vs local
    keep <- !is.na(mydata$var_a) & !is.na(mydata$var_b) & !is.na(mydata$var_c)
    mydata <- mydata[keep,]
  }
  
  if(crop[m] == "wheat"){
    #for wheat
    #remove items (varieties) tested in only 1 season
    #here we remove items, not the entire row, rows with more than 1 NA will be removed next
    rm_item <- as.data.frame(table(mydata$variety_a, mydata$season))
    rm_item[rm_item > 0] <- 1
    rm_item <- aggregate(rm_item[,"Freq"], by=list(rm_item[,"Var1"]), sum)
    rm_item <- as.vector(c(as.character(rm_item[rm_item$x<2, 1 ]), "HUW468"))
    
    for(i in 1:length(rm_item)){
      mydata$variety_a <- ifelse(mydata$variety_a == rm_item[i], NA, mydata$variety_a)
      mydata$variety_b <- ifelse(mydata$variety_b == rm_item[i], NA, mydata$variety_b)
      mydata$variety_c <- ifelse(mydata$variety_c == rm_item[i], NA, mydata$variety_c)
    }
    
    #remove rankings with more than 1 NA per observer
    keep <- NULL
    for(i in 1:nrow(mydata)){
      if(sum(is.na(mydata[i,c("variety_a","variety_b","variety_c")])) > 1) keep[i] <- FALSE else keep[i] <- TRUE
    }
    
    #remove wrong evaluations best == worst in overall_performance
    keep <- ifelse(mydata$best==mydata$worst, FALSE, keep)
    
    mydata <- mydata[keep, ]
  }
  
  # .............................................................
  # .............................................................
  # Check correlation between overall performance and yield
  # this information is only available for wheat (IND) and common beans (NIC)
  if(crop[m] != "durumwheat"){
    # rankings from overall performance
    overall <- mydata[,c("variety_a","variety_b","variety_c",
                         "characteristic","best","worst")]
    # rankings from yield
    yield <- mydata[,c("variety_a","variety_b","variety_c",
                       "yield","best_yield","worst_yield")]
    # make sure that both datasets has the same names
    names(yield) <- names(overall)
    
    # remove NA's in yield
    ykeep <- !is.na(yield$best) & !is.na(yield$worst) &  yield$best != yield$worst
    
    # keep those rows with no NA's
    yield <- yield[ykeep,]
    overall <- overall[ykeep,]
    
    # get rankings for yield
    yield <- to_rankings(yield, 
                         items = c("variety_a","variety_b","variety_c"),
                         input = c("best","worst"),
                         type = "tricot",
                         grouped.rankings = TRUE)
    
    # parsed matrix of rankings 
    yield <- yield[1:length(yield), , as.grouped_rankings = FALSE]
    
    # get rankings for overall performance
    overall <- to_rankings(overall, 
                           items = c("variety_a","variety_b","variety_c"),
                           input = c("best","worst"),
                           type = "tricot",
                           grouped.rankings = TRUE)
    
    # parsed matrix of rankings 
    overall <- overall[1:length(overall), , as.grouped_rankings = FALSE]
    
    # Calculate Kendall's correlation
    # export this output
    ktau <- kendallTau(yield, overall)
    
    capture.output(ktau, file = paste0(output, "kendall_tau_yield_vs_overall.txt"))
  }
  
  # .............................................................
  # .............................................................
  # Take explanatory variables in a separate dataframe
  covar <- cbind(mydata[, c("season","soil","lat","lon",
                             "xy","yx","planting_day","planting_date")],
                 mydata[, c(30:(ncol(mydata)-4))])
  
  covar <- covar[,!names(covar) %in% c("maxNT_VEG","planting_date")]
  covar <- covar[,!grepl("ETo_", names(covar))]
  covar <- covar[,!grepl("sd", names(covar))]
  covar <- covar[,!grepl("WB_", names(covar))]
  
  # Check variance in explanatory variables
  varout <- caret::nearZeroVar(covar)
  cat("\n Removing these variables with near zero variance: \n", 
      sort(names(covar)[varout]),"\n")
  
  # Remove variables with near zero variance
  covar <- covar[,-varout]
  
  # Take rankings from mydata 
  if (crop[m] != "durumwheat") {
    mydata <- cbind(mydata[, c("variety_a","variety_b","variety_c",
                               "characteristic","best","worst",
                               "overall_vs_local","var_a","var_b","var_c")])
  }
  
  if (crop[m] == "durumwheat") {
    mydata <- cbind(mydata[, c("package","variety_a","variety_b","variety_c","variety_d",
                               "rank_variety_a","rank_variety_b","rank_variety_c","rank_variety_d")])
  }
  
  #Take the number of rows in this dataset
  n <- nrow(mydata)
  
  # Generate a Plackett-Luce grouped rankings
  ## for common beans we add the comparison with the local var as partial rankings
  ## in this approach the three variaties from tricot are compared together in the same rank
  ## then each variety is compared with local as an individual rank (additional.rank)
  if(crop[m] == "commonbean"){
    G <- to_rankings(data = mydata,
                     items = c("variety_a","variety_b","variety_c"),
                     input = c("best","worst"),
                     type = "tricot",
                     add.rank = mydata[c("var_a","var_b","var_c")],
                     grouped.rankings = TRUE)
  }
  
  ## information about local variety is not available for wheat
  ## then a simple grouped rankings is created
  if(crop[m] == "wheat"){
    G <- to_rankings(data = mydata,
                     items = c("variety_a","variety_b","variety_c"),
                     input = c("best","worst"),
                     type = "tricot",
                     grouped.rankings = TRUE)
  }
  
  ## durum data is another format
  ## reorganise the dataset to convert into grouped rankings
  if(crop[m] == "durumwheat"){
    
    #convert data in long format
    R <- bind_cols(melt(mydata[,c(1,2:5)], id.vars = c("package")),
                   melt(mydata[,c(1,6:9)], id.vars = c("package")))
    
    #keep item names and rankings
    R %<>% 
      dplyr::select(c(1,3,6)) %>%
      set_colnames(c("package","var","rank")) %>%
      arrange(., package)
    
    #check the number of observations per variety
    #remove those tested in less than 20 farms 
    rmitem <- R %>% 
      group_by(var) %>%  
      count(var) %>%
      filter(n < 20) %>%
      dplyr::select(var) %>%
      as.matrix() %>%
      as.vector()
    
    #take variety names
    items <- sort(unique(R[,"var"]))
    
    #match those to remove
    rmitem <- match(rmitem, items)
    
    #convert dataset from long format into a sparsed matrix
    R <- sparseMatrix(i = rep(1:n, each = 4),
                      j = match(R[,"var"], items),
                      x = R[,"rank"],
                      dimnames = list(as.integer(unique(R[,"package"])), items ))
    #sparsed matrix into a ordinary matrix from R package base
    R <- as.matrix(R)
    #remove varieties selected previously
    R <- R[,-rmitem]
    #remove rows with less tham 2 items
    R <- R[-which(rowSums(R > 0) < 2), ]
    
    #row names here are the packages codes, get this to merge rankings and explanatory variables
    keep <- row.names(R)
    
    #add package codes to covar
    covar <- cbind(package = mydata[,"package"],covar)
    
    #subset covar
    covar <- covar[covar$package %in% keep , ]
    #reorder by packages
    covar %<>% arrange(., package)
    
    #define new n
    n <- nrow(R)
    
    #rename rows the rankings
    row.names(R) <- 1:n
    row.names(covar) <- row.names(R)
    
    pack <- mydata[,"package"]
    
    #convert to a grouped rankings
    G <- grouped_rankings(as.rankings(R), seq_len(n))
    
    #remove packages codes
    covar <- covar[,-1]
    
  }
  
  cat("This analysis will use",n,"observations.\n")
  
  # Check network
  source("script/helper_01_check_network.R")
  
  
  # merge grouped rankings with explanatory variables
  mydata <- cbind(G, covar)
  
  # .............................................................
  # .............................................................
  # Set parameters for forward selection ####
  
  # Define folds based on which season this crop was evaluated
  folds <- as.integer(as.factor(as.character(mydata$season)))
  # number of folds
  k <- max(folds)
  # minsize
  minsize <- round((n * 0.3), -1)
  # bonferroni = TRUE
  bonferroni <- TRUE
  # mean method
  mean.method <- "foldsize"
  # akaike.weights
  akaike.weights <- TRUE
  # alpha
  alpha <- 0.01
  # npseudo
  npseudo <- 1.5

  # .............................................................
  # .............................................................
  # Run forward cross-validation ####
  # mod <- forward(G ~ .,
  #                data = mydata, 
  #                k = k,
  #                folds = folds,
  #                select.by = "deviance",
  #                akaike.weights = akaike.weights,
  #                ncores = n_cpu,
  #                mean.method = mean.method,
  #                minsize = minsize, 
  #                alpha = alpha,
  #                bonferroni = bonferroni)
  
  best_model <- as.formula(mod$raw$call)
  best_model <- all.vars(best_model)[-1]
  
  # .............................................................
  # .............................................................
  # Run blocked cross-validation ####
  
  # Define list of explanatory variables to use in each approach
  mydata$P1 <- rep(1, n)
  
  attrib <- list(zero = c("P1"),
                 desi = c("lon","lat","xy","yx","soil","planting_day"),
                 clim = best_model,
                 clsp = c(best_model, c("lon","lat","yx") ))
  
  attrib <- lapply(attrib, function(x){
    
    x <- paste(x, collapse = " + ")
    
    x <- paste0("G ~ ", x)
    
    x <- as.formula(x)
    
  })
  
  blocked_cv <- lapply(attrib, function(x){
    crossvalidation(x, 
                    data = mydata, 
                    k = k, 
                    folds = folds,
                    minsize = minsize,
                    alpha = alpha, 
                    bonferroni = bonferroni, 
                    npseudo = npseudo,
                    mean.method = mean.method)
  })
  
  # combine model coefficients
  blocked_cv <- lapply(blocked_cv, function(x){
    bind_cols(call = x$raw$call, 
              x$coeffs)
  })
  
  blocked_cv <- bind_rows(blocked_cv)
  
  # .............................................................
  # .............................................................
  # Run average season using historical data ####
  cat("Starting average season with blocked cross-validation", date() ,"\n")
  
  #load climatology data
  load(paste0("data/",crop[m],"/",
              crop[m],"_climatology.rda"))
  #define number of predictions to be generated n.years * n.estimated planting dates
  npreds <- length(climatology) * length(climatology[[1]])
  
  #list into array
  nr <- nrow(climatology[[1]][[1]])
  nc <- ncol(climatology[[1]][[1]])
  
  climatology %<>%
    unlist(.) %>%
    array(., dim = c(nr, nc, npreds),
          dimnames = list(1:nr, names(climatology[[1]][[1]]), 1:npreds ))
  
  # durumwheat have a particular order, filter data to fit this order
  if(crop[m]=="durumwheat"){
    
    tf <- as.data.frame(cbind(pack, climatology[[1]][[1]]))
    
    keep <- tf %<>%
      filter(., pack %in% keep ) %>%
      rownames_to_column() %>%
      arrange(., pack) %>%
      dplyr::select(1) %>%
      as.matrix() %>%
      as.vector()
    
    rm(tf)
  }
  
  
  pb <- txtProgressBar(min = 1, max = npreds, style = 3)
  
  avg_season <- NULL
  
  for(a in seq_len(npreds)){
    
    #temporary dataframe with G and explanatory variables
    a_df <- cbind(G, as.data.frame(climatology[keep, , a]))
    
    model <-  crossvalidation(attrib$clim, 
                              data = mydata, 
                              k = k, 
                              folds = folds,
                              minsize = minsize,
                              alpha = alpha, 
                              bonferroni = bonferroni, 
                              npseudo = npseudo,
                              mean.method = mean.method)
    
   
    
    estimators <- model[["coeffs"]]
    
    avg_season <- rbind(avg_season, estimators)
    
    setTxtProgressBar(pb, a)
    
  }
  
  close(pb)
  
  cat("End of average season \n")
  
  avg_season <- as.data.frame(colMeans(avg_season))
  
  avg_season <- as_tibble(t(avg_season))
  
  avg_season <- mutate(avg_season, 
                       call = "G ~ average_season")
  
  blocked_cv <- bind_rows(blocked_cv, avg_season)
  
  # write matrix as a .csv file
  write_csv(blocked_cv, paste0(output, "performance_PLT_blocked_crossvalidation.csv"))
  
  
  # .............................................................
  # .............................................................
  # Run k-fold cross-validation for the different model approaches ####
  cat("Starting k-fold cross-validation",date(),"\n")
  #Define number of folds
  k <- 10
  #Define samples
  folds <- sample(rep(1:k, times = ceiling(n/k), 
                      length.out=n), replace=FALSE)
  
  
  kfolds_cv <- lapply(attrib, function(x){
    crossvalidation(x, 
                    data = mydata, 
                    k = k, 
                    folds = folds,
                    minsize = minsize,
                    alpha = alpha, 
                    bonferroni = bonferroni, 
                    npseudo = npseudo,
                    mean.method = "equal")
  })
  
  # combine model coefficients
  kfolds_cv <- lapply(kfolds_cv, function(x){
    bind_cols(call = x$raw$call, 
              x$coeffs)
  })
  
  kfolds_cv <- bind_rows(kfolds_cv)
  
  #write matrix as a csv file
  write_csv(kfolds_cv, paste0(output, "performance_PLT_10fold_crossvalidation.csv"))
  
  
  # .............................................................
  # .............................................................
  # Fit the best model from k-fold and generate Plackett-Luce plots ####
  
  #Fit pltree with best model
  cat("Fit pltree with best model \n")
  tree <- pltree(attrib$clim,
                 data = mydata, 
                 alpha = alpha, 
                 minsize = minsize,
                 npseudo = npseudo)
  
  #export summary of fitted pltree
  capture.output(tree, summary(tree), file = here(output, "pltree_clim.txt") )
  
  #Write this tree as a .svg file
  svg(filename = here(output, "PLTree.svg"), 
      width=15, 
      height=15, 
      pointsize=12)
  partykit::plot.modelparty(tree) 
  dev.off()
  
  #Generate plots with error bars using qvcal
  labels <- NULL
  
  plots <- plot_nodes(tree, labels = labels, font.size = c(20, 23))
  #define file names
  nodepaths <- paste0(here(output, paste0("PLTree_", names(plots) ,".svg")))
  #put names in a list
  nodepaths <- as.list(nodepaths)
  #define names of each element in this list
  names(nodepaths) <- names(plots)
  
  #write plots as .svg file 
  h <- length(itemnames) + 1.5
  mapply(function(X, Y){
    
    ggsave(filename = Y, plot = X,
           dpi = 600, width = 15, height = h, units = "cm")
    
  }, X = plots, Y = nodepaths )
  
  # .............................................................
  # .............................................................
  #Make a histogram to show the distribution of ####
  ## explanatory variables among blocks (seasons)
  
  #get names of vars within nodes
  vars <- partykit:::.list.rules.party(tree)
  vars <- strsplit(vars, "&|<=|[>]")
  mvalue <- which.max(as.vector(unlist(lapply(vars, length))))
  vars <- vars[[mvalue]]
  vars <- vars[seq(1, length(vars)-1, 2)]
  vars <- gsub(" ", "", vars)
  
  #get split values
  ni <- nodeids(tree)
  ni_terminal <- nodeids(tree, terminal = TRUE)
  ni_inner <- ni[!ni %in% ni_terminal]
  breaks <- sapply(ni_inner, function(X){
    split_node(node_party(tree[[X]]))$breaks} )
  
  for(i in unique(vars)){
    
    b <- breaks[ vars%in% vars[i] ]
    
    h <- ggplot(mydata, aes(mydata[, vars[i] ], fill = season)) +
      geom_histogram(bins = 50) +
      geom_vline(xintercept = b, col = "red") + 
      labs(x = "", y = "Count") +
      scale_fill_brewer(palette = "GnBu", direction = -1, name = "") +
      scale_x_continuous(expand = expand_scale(mult = c(0, .05))) +
      scale_y_continuous(expand = expand_scale(mult = c(0, .01))) +
      theme_classic() +
      theme(axis.text = element_text(size=20, colour="black"),
            axis.text.x = element_text(size=20, angle = 0, 
                                       hjust=0.5,vjust=1, face="plain"),
            axis.text.y = element_text(size=20, angle = 0, 
                                       hjust=1, vjust=0.5, face="plain"),
            axis.title=element_text(size=20,face="bold"),
            axis.line = element_line(),
            legend.position=c(.15, .8),
            legend.text=element_text(size=20),
            legend.background = element_rect(colour = NA),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.background = element_blank())
    
    ggsave(here(output, paste0("histogram_",vars[i],".svg")), 
           plot = h, dpi = 600,
           width = 25, height = 20, units = "cm")
    
  }
  
  
  # .............................................................
  # .............................................................
  # Calculate worst regret ####
  #Regret is the difference with the best variety in each node
  WR <- worstRegret(tree)
  
  if(!is.null(labels)) {WR$items <- labels}
  
  write.csv(WR, here(output, "worst_regret.csv"))
  
  
  # .............................................................
  # .............................................................
  # Plot network connections among items ####
  R <- mydata$G
  
  R <- R[1:length(R), , as.grouped_rankings = FALSE]
  
  R <- as.rankings(R)
  
  adj <- PlackettLuce::adjacency(R)
  
  adj <- as.vector(adj)
  
  adj <- t(matrix(adj, nrow = ncol(R), ncol = ncol(R)))
  
  dimnames(adj) <- list(dimnames(R)[[2]], dimnames(R)[[2]])
  
  adj <- btdata(adj, return_graph = TRUE)
  
  svg(filename = here(output, "connections.svg"), 
      width=30, 
      height=30, 
      pointsize=55)
  plot.igraph(adj$graph, vertex.size = 30, edge.arrow.size = 0.1) 
  dev.off()
  
  
  # .............................................................
  # .............................................................
  # Check if elevation play a role in the performance of varieties in Ethiopia ####
  if(crop[m]=="durumwheat"){
    nodes <- partykit::nodeids(tree, terminal = TRUE)
    
    #get the estimates from nodes
    node <- list()
    
    for(i in seq_along(nodes)){
      node[[i]] <- tree[[ nodes[i] ]]$node$info$object %>%
        itempar(., vcov = TRUE, alias = TRUE) %>%
        qvcalc::qvcalc.itempar(.) %>%
        extract2(2) %>%
        as.matrix() %>%
        as_tibble() %>%
        mutate(items = items)
    }
    
    node <- bind_cols(items = node[[1]]$items, n1 = node[[1]]$estimate, n2 = node[[2]]$estimate )
    
    #read the passport data
    pass <- read_csv("data/durumwheat/passportdata_eth.csv")
    
    pass %<>%
      mutate(items = accession) %>%
      inner_join(. , node, by = "items" ) %>%
      filter(. , !is.na(x) & !is.na(y)) %>%
      dplyr::select(. , items, x , y , n1, n2, adm1, year_charac)
    
    #get the elevation from origin of items
    #read the elevation data
    #source http://gisweb.ciat.cgiar.org/TRMM/SRTM_Resampled_250m/
    r <- stack("data/durumwheat/eth_elev_srtm.tif")
    
    e <- extract(r, pass[c("x","y")], buffer = 1000)
    
    names(e) <- paste0(pass$items,".")
    
    e <- unlist(e, use.names = TRUE)
    
    e <- as.data.frame(cbind(items = str_split(names(e), 
                                               pattern = "[.]", simplify = TRUE)[,1],
                             elev = e))
    
    rownames(e) <- 1:nrow(e)
    
    #add to the main data
    pass %<>%
      inner_join(., e , by = "items", all.y= TRUE ) %>%
      mutate(. , dif = n1 - n2,
             group = ifelse(n1 >= mean(n1) & n2 <= mean(n2), 1,
                            ifelse(n1 <= mean(n1) & n2 >= mean(n2), 2, NA)),
             elev = as.numeric(as.character(elev)))
    
    ttest <- with(pass,
                  t.test(elev ~ group))
    
  }
  
  
  
}