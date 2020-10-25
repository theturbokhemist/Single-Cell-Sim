Generate_results_object3 <- function(target_expression, standard_object, nth_neighbor = 20, cutoff_percent = 0.85) {
  
  reference_expression <- standard_object$expression_ref
  
  #Intialize distance_df
  ncol_df <- 3
  nrow_df <- nrow(target_expression)
  
  distance_df <- data.frame(matrix(ncol = ncol_df , nrow = nrow_df))
  colnames(distance_df) <- c("Model_ID","Nearest_Neighbor_Euc", "Neighbor_Model_ID")
  
  distance_df[, "Model_ID"] <- 1:nrow_df
  
  
  #Calc distances
  dists <- as.matrix(pdist(target_expression, target_expression))
  
  #Generate list of distances
  for (i in 1:nrow(distance_df)) {
    
    distance_df[i, "Neighbor_Model_ID"] <- order(dists[,i])[(nth_neighbor + 1)]
    
    distance_df[i, "Nearest_Neighbor_Euc"] <- dists[,i][distance_df[i, "Neighbor_Model_ID"]]
    
  }
  
  cutoff <- distance_df[, "Nearest_Neighbor_Euc"][order(distance_df[, "Nearest_Neighbor_Euc"])[floor(cutoff_percent*nrow(distance_df))]]
  
  
  #Distribution of the distance to the nth nearest neighbor for each model
  distribution_target_plot <- ggplot(data = distance_df, aes(Nearest_Neighbor_Euc)) + geom_histogram(binwidth = 0.01) + ylim(0, 500) + 
    geom_vline(xintercept = cutoff, color = "red", size = .3)
  
  
  #Remove all target models that have nth neighbor farther than cutoff
  distance_df_filtered <- distance_df[which(distance_df[,"Nearest_Neighbor_Euc"] <= cutoff),]
  
  target_expression_filtered <- target_expression[distance_df_filtered[,"Model_ID"],]
  
 
  #Overlay Plots
  
  overlay_unfiltered <- ggplot() + 
    geom_point(data = as.data.frame(reference_expression), aes(A, B, color = "DF1")) +
    geom_point( data = as.data.frame(target_expression), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
    geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
    labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
    scale_colour_manual(values = c(DF1 = "black", DF2 = "red", DF3 = "yellow"), 
                        labels = c(DF1 = "Reference", DF2 = "Target", DF3 = "Centers"),
                        name = "Data") + 
    guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))
  
  
  overlay_filtered<- ggplot() + 
    geom_point(data = as.data.frame(reference_expression), aes(A, B, color = "DF1")) +
    geom_point( data = as.data.frame(target_expression_filtered), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
    geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
    labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
    scale_colour_manual(values = c(DF1 = "black", DF2 = "red", DF3 = "yellow"), 
                        labels = c(DF1 = "Reference", DF2 = "Target_Filtered", DF3 = "Centers"),
                        name = "Data") + 
    guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))
  
  overlay_list_target <- list(unfiltered = overlay_unfiltered, filtered = overlay_filtered)
  
  
  #########################################Comparing to reference expression
  #Intialize ref_to_target_dist.df

  ncol <- 3
  nrow <- nrow(reference_expression)
  
  ref_to_target_dist.df <- data.frame(matrix(ncol = ncol , nrow = nrow))
  colnames(ref_to_target_dist.df) <- c("Ref_Model_ID","Nearest_Target_Euc", "Target_ID")
  
  ref_to_target_dist.df[, "Ref_Model_ID"] <- 1:nrow
  
  #Distance of each reference model to target model
  dists2 <- as.matrix(pdist(reference_expression, target_expression_filtered))
  
  
  i <- 1
  for (i in 1:nrow(ref_to_target_dist.df)) {
    
    ref_to_target_dist.df[i, "Target_ID"] <- order(dists2[i, ])[nth_neighbor]
    
    ref_to_target_dist.df[i, "Nearest_Target_Euc"] <- dists2[i,][ref_to_target_dist.df[i, "Target_ID"]]
    
  }
  
  #Distribution of the distance to the nth nearest neighbor for each model
  distribution_reference_plot <- ggplot(data = ref_to_target_dist.df, aes(Nearest_Target_Euc)) + geom_histogram(binwidth = 0.01) + ylim(0, 5000) + 
    geom_vline(xintercept = cutoff, color = "red", size = .3)
  
  
  #Remove all target models that have nth neighbor farther than cutoff
  
  kept_models <- which(ref_to_target_dist.df[,"Nearest_Target_Euc"] <= cutoff)
  
  
  ref_to_target_dist.df_filtered <- ref_to_target_dist.df[kept_models,]
  
  reference_expression_filtered <- reference_expression[ref_to_target_dist.df_filtered[,"Ref_Model_ID"],]
  
  #Overlay Plots
  
  overlay_ref <- ggplot() + 
    geom_point(data = as.data.frame(reference_expression), aes(A, B, color = "DF1")) +
    geom_point( data = as.data.frame(reference_expression_filtered), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
    geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
    labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
    scale_colour_manual(values = c(DF1 = "black", DF2 = "red", DF3 = "green"), 
                        labels = c(DF1 = "Reference", DF2 = "Filtered Reference", DF3 = "Centers"),
                        name = "Data") + 
    guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))
  
  
  overlay_ref.target<- ggplot() + xlim(-3.5,2) + ylim(-3.5, 2) +
    geom_point(data = as.data.frame(reference_expression_filtered), aes(A, B, color = "DF1")) +
    geom_point( data = as.data.frame(target_expression_filtered), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
    geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
    labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
    scale_colour_manual(values = c(DF1 = "black", DF2 = "red", DF3 = "green"), 
                        labels = c(DF1 = "Filtered Reference", DF2 = "Filtered Target", DF3 = "Centers"),
                        name = "Data") + 
    guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))
  
  overlay_list_ref <- list(unfiltered = overlay_ref, filtered = overlay_ref.target)
  
  
  ###For 2 genes
  
  #Distances Distribution
  # Distance_density_plot <- ggplot() +geom_density(data = distance_df, aes(Nearest_Ref_Neighbor_Euc)) + 
  #   labs(title="Distance to Nearest Reference Neighbor") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Euclidean Distance")
  
  
  #geom_vline(xintercept = cutoff, color = "red") +
  
  filtered_params_df <- as.data.frame(sracipeParams(standard_object$reference_rset)[kept_models,])
  
  
  output <- list(cutoff = cutoff, distance_df = distance_df, distance_df_filtered = distance_df_filtered, target_expression_filtered = target_expression_filtered,
                 distribution_target_plot = distribution_target_plot, overlay_list_target = overlay_list_target, ref_to_target_dist.df = ref_to_target_dist.df, 
                 ref_to_target_dist.df_filtered = ref_to_target_dist.df_filtered, reference_expression_filtered = reference_expression_filtered, 
                 distribution_reference_plot = distribution_reference_plot, overlay_list_ref = overlay_list_ref, filtered_params_df = filtered_params_df)
  output
}

test3 <- Generate_results_object3(target_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z, standard_object = standard_object_TS, nth_neighbor = 20, 
                                  cutoff_percent = 0.85)


test3$distribution_target_plot

test3$distribution_reference_plot

dim(test3$dists2)

test3$distribution_target_plot

do.call("grid.arrange", test3$overlay_list_target)

do.call("grid.arrange", test3$overlay_list_ref)
#13928499
#########################################################################################Resample reference models from previous reference models
resample_reference <- function(reference_object, identified_reference_parameters) {
  
  topology <- reference_object$topology
  numModels <- nrow(reference_object$expression_ref)
  
  resampled_reference_rset <- sracipeSimulate(topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = 0.02)
  
  resampled_reference_params_df <- sracipeParams(resampled_reference_rset) 
  
  
  ####
  nparams <- ncol(identified_reference_parameters)
  param_names <- colnames(resampled_reference_params_df)
  hill_params <- grep(pattern = "N", param_names)
  
  nmodels_resampled <- nrow(resampled_reference_params_df)
  nmodels_filtered <- nrow(identified_reference_parameters)
  
  ###
  
  
  for (i in 1:nparams) {
    
    resampled_reference_params_df[, i] <- sample(identified_reference_parameters[,i], nmodels_resampled, replace=TRUE)
    
  }
  
  #resimulate
  sracipeParams(resampled_reference_rset) <- resampled_reference_params_df
  
  resampled_reference_rset <- sracipeSimulate(resampled_reference_rset, integrate = TRUE, genParams = FALSE)
  
  #renormalize
  
  resampled_expr_reference_log <- log2(t(assay(resampled_reference_rset)))
  
  resampled_rset_expr_log_z <- sweep(resampled_expr_reference_log,
                                     2, reference_object$means_ref, FUN = "-")
  
  resampled_rset_expr_log_z <- sweep(resampled_rset_expr_log_z,
                                     2, reference_object$sds_ref, FUN = "/")
  
  #rotate
  resampled_simulated_rset_expr_PCs <- resampled_rset_expr_log_z %*% reference_object$prcomp_ref$rotation
  
  reference_object$resampled_expression <- resampled_rset_expr_log_z
  
  reference_object$resampled_expression_PCs <- resampled_simulated_rset_expr_PCs
  
  reference_object$resampled_reference_rset <- resampled_reference_rset
  
  reference_object
  
}


resampled_standard_object <- reference_resampled(standard_object = standard_object_TS, results_object = test3)  
  
ggplot() + 
  geom_point(data = as.data.frame(standard_object_TS$expression_ref), aes(A, B, color = "DF1")) +
  geom_point( data = as.data.frame(resampled_standard_object$expression_ref), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
  labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "black", DF2 = "red", DF3 = "yellow"), 
                      labels = c(DF1 = "Reference", DF2 = "Target", DF3 = "Centers"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))
  
test5 <- Generate_results_object3(target_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z, standard_object = resampled_standard_object, nth_neighbor = 20, 
                                  cutoff_percent = 0.85)

do.call("grid.arrange", test5$overlay_list_target)

do.call("grid.arrange", test5$overlay_list_ref)

test6 <- extract_param_ranges(results_object = test5, standard_object = standard_object_TS, nbins = 10, cutoff = 0.80)

do.call("grid.arrange", test6$plots_list)
#########################################################################################
extract_param_ranges <- function(results_object, standard_object, nbins = 10, cutoff = 0.9) {
  
  nparams <- ncol(results_object$filtered_params_df)
  param_names <- colnames(sracipeParams(standard_object$reference_rset))
  hill_params <- grep(pattern = "N", param_names)
  
  ###
  uniform_counts <- nrow(results_object$filtered_params_df)/nbins
  
  
  counts_tables <- vector(mode = "list", length = nparams)
  names(counts_tables) <- param_names
  
  cutoff_tables <- vector(mode = "list", length = nparams)
  names(cutoff_tables) <- param_names
  
  plots_list <- vector(mode = "list", length = nparams)
  
  extracted_param_ranges <- vector(mode = "list", length = nparams)
  names(extracted_param_ranges) <- param_names
  ###
  
  for (i in 1:nparams) {
    
    if (i %in% hill_params) {
      
      uniform_counts_hill <- nrow(results_object$filtered_params_df)/length(standard_object$standard_state[[i]][1]:standard_object$standard_state[[i]][2])
      
      table <- table(results_object$filtered_params_df[,i])
      
      counts_tables[[i]] <- table
      
      ###
      
      cutoff_tables[[i]] <- table/uniform_counts_hill
      
    } else {
      
      bin_width <- (standard_object$standard_state[[i]][2] - standard_object$standard_state[[i]][1])/nbins
      
      breaks_vec <- seq(standard_object$standard_state[[i]][1], standard_object$standard_state[[i]][2], by = bin_width)
      
      table <- table(cut(results_object$filtered_params_df[,i], breaks = breaks_vec, include.lowest = TRUE))
      
      counts_tables[[i]] <- table
      
      ###
      
      cutoff_tables[[i]] <- table/uniform_counts
      
    }
    
  }
  
  i <- 1
  
  for (i in 1:length(plots_list)) {
    
    plots_list[[i]] <- ggplot(as.data.frame(cutoff_tables[[i]]), aes(x = Var1, weight = Freq)) + geom_bar() + geom_hline(yintercept = cutoff, color = "red") +
      ggtitle(param_names[i]) + xlab("Param Value") +ylab("Ratio")
    
  }
  
  #extract params
  i <- 1
  for (i in 1:length(cutoff_tables)) {
    
    df <- as.data.frame(cutoff_tables[[i]])
    cutoff_rows <- which(df[,"Freq"] >= cutoff)
    
    df <- df[cutoff_rows,]
    
    ranges_as.char <- as.character(df[, "Var1"])
    
    extract_numerics <- as.numeric(unlist(regmatches(ranges_as.char,
                                                     gregexpr("[[:digit:]]+\\.*[[:digit:]]*",ranges_as.char))))
    
    param_range <- range(extract_numerics)
    
    extracted_param_ranges[[i]] <- param_range
    
  }
  
  
  
  
  output <- list(counts_tables = counts_tables, cutoff_tables = cutoff_tables, plots_list = plots_list, extracted_param_ranges = extracted_param_ranges)
  
  
}


test4 <- extract_param_ranges(results_object = test3, standard_object = standard_object_TS, nbins = 10, cutoff = 1)
test4$extracted_param_ranges
do.call("grid.arrange", test4$plots_list)
########################################################################################

test <- as.matrix(pdist(experimental_object_TS.G_A.l$simulated_rset_expr_log_z, experimental_object_TS.G_A.l$simulated_rset_expr_log_z))
test[1:10, 1:10]
order(test[,1])[(21)]

test3 <- Generate_results_object3(target_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z, standard_object = standard_object_TS, nth_neighbor = 20, 
                                  cutoff = 0.3)
test3$overlay_list
do.call("grid.arrange", test3$overlay_list)
cutoff3 <- 0.3

ggplot(data = test3, aes(Nearest_Neighbor_Euc)) + geom_histogram(binwidth = 0.001) + ylim(0, 50) + geom_vline(xintercept = cutoff3, color = "red", size = .3)

test3_filter <- test3[which(test3[,1] <= cutoff3),]

test3_expr <- experimental_object_TS.G_A.l$simulated_rset_expr_log_z[as.numeric(rownames(test3_filter)),]

experimental_object_TS.G_A.l$simulated_rset_expr_log_z[1:6,]
test3_expr[1:5,]

ggplot() + 
  geom_point(data = as.data.frame(test3_expr), aes(A, B)) + xlim(-3.5, 1) + ylim(-3,2)


ggplot() + 
  geom_point(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(A, B, color = "DF1")) +
  geom_point( data = as.data.frame(test3_expr), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
  labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "black", DF2 = "red", DF3 = "yellow"), 
                      labels = c(DF1 = "Reference", DF2 = "Target", DF3 = "Centers"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))


##############################################
Generate_results_object4 <- function(experimental_expression, num_neighbors = 5) {
  
  
  #Intialize distance_df
  ncol_df <- 2
  nrow_df <- nrow(experimental_expression)
  
  distance_df <- data.frame(matrix(ncol = ncol_df , nrow = nrow_df))
  colnames(distance_df) <- c("Mean_Nearest_Neighbors", "Model_ID")
  distance_df[,"Model_ID"] <- 1:nrow_df
  
  neighbors_list <- vector(mode = "list", length = nrow_df)
  
  #Calc distances
  dists <- as.matrix(pdist(experimental_expression, experimental_expression))
  
  #Generate list of distances
  for (i in 1:nrow(distance_df)) {
    
    nearest_neighbors <- order(dists[,i])[2:(num_neighbors + 1)]
    neighbors_list[[i]] <- nearest_neighbors
    
    distance_df[i, 1] <- mean(dists[,i][nearest_neighbors])
    
  }
  
  
  ###For 2 genes
  
  #Distances Distribution
  # Distance_density_plot <- ggplot() +geom_density(data = distance_df, aes(Nearest_Ref_Neighbor_Euc)) + 
  #   labs(title="Distance to Nearest Reference Neighbor") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Euclidean Distance")
  
  
  #geom_vline(xintercept = cutoff, color = "red") +
  output <- list(distance_df = distance_df, neighbors_list = neighbors_list)
}

test4 <- Generate_results_object4(experimental_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z, num_neighbors = 20)
test4 <- test4$distance_df

cutoff4 <- 0.1

ggplot(data = test4, aes(Mean_Nearest_Neighbors)) + geom_histogram(binwidth = 0.001) + ylim(0, 50) + geom_vline(xintercept = cutoff4, color = "red", size = .3)

test4_filter <- test4[which(test4[,1] <= cutoff4),]

test4_expr <- experimental_object_TS.G_A.l$simulated_rset_expr_log_z[as.numeric(rownames(test4_filter)),]

experimental_object_TS.G_A.l$simulated_rset_expr_log_z[1:6,]
test4_expr[1:5,]

ggplot() + 
  geom_point(data = as.data.frame(test4_expr), aes(A, B))

ggplot() + 
  geom_point(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(A, B, color = "DF1")) +
  geom_point( data = as.data.frame(test4_expr), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
  labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "black", DF2 = "blue", DF3 = "red"), 
                      labels = c(DF1 = "Reference", DF2 = "Target", DF3 = "Centers"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))



ggplot() + 
  geom_point(data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B))

ggplot() + 
  geom_point(data = as.data.frame(standard_object_TS$expression_ref), aes(A, B))


ggplot() + 
  geom_point(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(A, B, color = "DF1")) +
  geom_point( data = as.data.frame(test3_expr), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8)
