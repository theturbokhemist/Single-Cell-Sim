##############################Functions#################################################################################################################################


##############################generate_reference_object#########################################################################################################################
generate_reference_object <- function(topology, reference_parameter_ranges, numModels = 10000, integrateStepSize = 0.02) {
  
  reference_rset <- sracipeSimulate(topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  state <- reference_parameter_ranges
  
  #hill gene indexes
  hill_genes <- grep(pattern = "N", colnames(sracipeParams(reference_rset)))
  
  #assign vaues from uniform samples distributions to each model of reference rset.  
  for (i in 1:length(state)) {
    
    sracipeParams(reference_rset)[,i] <- runif(numModels, min= as.numeric(state[[i]][1]), max= as.numeric(state[[i]][2]))  
    
  }
  
  for (k in 1:length(hill_genes)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(reference_rset)[,hill_genes[k]] <- round(runif(numModels, min= as.numeric(state[[hill_genes[k]]][1]) - 0.5, max= as.numeric(state[[hill_genes[k]]][2]) + 0.5))  
    
  }
  
  reference_params.df <- sracipeParams(reference_rset)
  
  ref_expression_rset <- sracipeSimulate(reference_rset, integrate = TRUE, genParams = FALSE)
  
  expr_reference_log <- log2(t(assay(ref_expression_rset)))
  
  means_reference <- colMeans(expr_reference_log)
  
  sds_reference <- apply(expr_reference_log, 2, sd)
  
  #z normalize
  expr_reference_log_z <- sweep(expr_reference_log, 
                                2, means_reference, FUN = "-")
  
  expr_reference_log_z <- sweep(expr_reference_log_z, 
                                2, sds_reference, FUN = "/")
  
  #find eigenvector,loading scores
  prcomp_reference <- prcomp(expr_reference_log_z, center = FALSE, scale. = FALSE)
  
  #list of all the ref information
  reference_rset_object <- list(reference_rset = ref_expression_rset, reference_params.df = reference_params.df, 
                                reference_parameter_ranges = reference_parameter_ranges, expression_ref = expr_reference_log_z,
                                prcomp_ref = prcomp_reference,
                                means_ref = means_reference, sds_ref = sds_reference, topology = topology)
  
  reference_rset_object
  
}

##############################generate_target_object#########################################################################################################################
generate_target_object <- function(target_parameter_ranges = target_parameter_ranges, reference_object, numModels = 10000, integrateStepSize = 0.02, save_rset = FALSE) {
  
  state <- target_parameter_ranges
  
  parameter_names <- colnames(sracipeParams(reference_object$reference_rset))
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  simulated_rset <- sracipeSimulate(reference_object$topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  for (t in 1:length(state)) {
    
    sracipeParams(simulated_rset)[,t] <- runif(numModels, min = state[[t]][[1]], max = state[[t]][[2]])  
    
  }
  
  #sample and round all the hill coefficients
  for (u in 1:length(hill_params)) {
    #you need to subtract 0.5 from lower bound and add 0.5 to upper bound in order to get a uniform distribution
    sracipeParams(simulated_rset)[,hill_params[u]] <- round(runif(numModels, min= round(state[[hill_params[u]]][[1]]) - 0.5, max= round(state[[hill_params[u]]][[2]]) + 0.5))  
    
  }
  target_params.df <- sracipeParams(simulated_rset)
  
  simulated_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE)
  
  simulated_rset_expr_log <- log2(t(assay(simulated_rset)))
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log,
                                     2, reference_object$means_ref, FUN = "-")
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log_z,
                                     2, reference_object$sds_ref, FUN = "/")
  
  #rotate
  simulated_rset_expr_PCs <- simulated_rset_expr_log_z %*% reference_object$prcomp_ref$rotation
  
  
  expression_object <- list(simulated_rset_expr_log_z = simulated_rset_expr_log_z, simulated_rset_expr_PCs = simulated_rset_expr_PCs)
  
  if (save_rset == TRUE) {
    
    expression_object <- list(simulated_rset_expr_log_z = simulated_rset_expr_log_z, simulated_rset_expr_PCs = simulated_rset_expr_PCs, simulated_rset = simulated_rset,
                              arguments = list(numModels = numModels, integrateStepSize = integrateStepSize, state = state), target_params.df = target_params.df)
    
  }
  expression_object
}


##############################identify_reference_models#########################################################################################################################
identify_reference_models <- function(target_expression, reference_object, nth_neighbor = 20, cutoff_percent = 0.85, n_closest_neighbors = NULL, centers, 
                                      weighting_factors = NULL) {
  
  #defining colors, determing number of models to sample in each cluster
  n_ref_models <- nrow(reference_object$expression_ref)
  
  if (!is.null(weighting_factors)) {
    
    num_clust <- length(weighting_factors)
    
    colors <- c("yellow", "springgreen1", "turqiose1", "purple")[1:num_clust]
    colors <- factor(colors, levels = colors)
    
    sample_vec <- c()
    
    for (i in 1:num_clust) {
      
      
      sample_vec <- c(sample_vec, floor(weighting_factors[i]*n_ref_models))
      
    }
    
    sample_vec[num_clust] <- n_ref_models - sum(sample_vec[1:(num_clust-1)])
    names(sample_vec) <- colors
    
  }

  

  
  
  if (!is.null(reference_object$resampled_expression)) {
    
    reference_expression <- reference_object$resampled_expression
    
    resampled_overlay <- ggplot() + 
      geom_point(data = as.data.frame(reference_object$expression_ref), aes(A, B, color = "DF1")) +
      geom_point( data = as.data.frame(reference_expression), aes(A, B, color = "DF2"), shape=1) +
      geom_point( data = as.data.frame(target_expression), aes(A, B, color = "DF3"), shape = 16, alpha = 2) +
      geom_point(data = as.data.frame(centers), aes(A, B, color="DF4"), shape=13, size=7, alpha = 3) + 
      labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
      scale_colour_manual(values = c(DF1 = "cadetblue3", DF2 = "red", DF3 = "black", DF4 = "green4"), 
                          labels = c(DF1 = "Original Reference", DF2 = "Current Reference", DF3 = "Target", DF4 = "Centers"),
                          name = "Data") + 
      guides(color = guide_legend(override.aes = list(shape = c(16, 1, 16, 13), size = c(1.5, 1.5, 1.5, 7)), title.hjust = .5)) 
    
    
  } else {
    
    reference_expression <- reference_object$expression_ref
    
  }
  
  ###expression plot
  density.ref.target.plot <- ggplot() + geom_density(data = as.data.frame(target_expression), aes(A, color = "A"), linetype = "dashed") +
    geom_density(data = as.data.frame(target_expression), aes(B, color = "B"), linetype = "dashed") +
    geom_density(data = as.data.frame(reference_expression), aes(A, color = "C")) +
    geom_density(data = as.data.frame(reference_expression), aes(B, color = "D")) +
    labs(title="Expression Density")  +  theme(plot.title = element_text(hjust = 0.5)) +   
    scale_colour_manual("Gene", values = c("blue", "red", "blue", "red"), labels = c(A = "Target(A)", B = "Target(B)", C = "Reference(A)", D = "Reference(B)"))
  
  
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
  # distance_df_filtered <- distance_df[which(distance_df[,"Nearest_Neighbor_Euc"] <= cutoff),]
  # 
  # target_expression_filtered <- target_expression[distance_df_filtered[,"Model_ID"],]
  
  
  #ref.target_overlay
  
  ref.target_overlay <- ggplot() + 
    geom_point( data = as.data.frame(reference_expression), aes(A, B, color = "DF2"), shape=1) +
    geom_point( data = as.data.frame(target_expression), aes(A, B, color = "DF3"), shape = 16, alpha = 2) +
    geom_point(data = as.data.frame(centers), aes(A, B, color="DF4"), shape=13, size=7, alpha = 3) + 
    labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
    scale_colour_manual(values = c(DF2 = "red", DF3 = "black", DF4 = "green4"), 
                        labels = c( DF2 = "Current Reference", DF3 = "Target", DF4 = "Centers"),
                        name = "Data") + 
    guides(color = guide_legend(override.aes = list(shape = c(1, 16, 13), size = c(1.5, 1.5, 7)), title.hjust = .5)) 
  
  
  #########################################Comparing to reference expression
  #Intialize ref_to_target_dist.df
  
  ncol <- 3
  nrow <- nrow(reference_expression)
  
  ref_to_target_dist.df <- data.frame(matrix(ncol = ncol , nrow = nrow))
  colnames(ref_to_target_dist.df) <- c("Ref_Model_ID","Nearest_Target_Euc", "Target_ID")
  
  ref_to_target_dist.df[, "Ref_Model_ID"] <- 1:nrow
  
  #Distance of each reference model to target model
  dists2 <- as.matrix(pdist(reference_expression, target_expression))
  
  
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
  
  #####keep reference that is close to targets
  
  if (!is.null(n_closest_neighbors)) {
    
    neighbors_list <- vector(mode = "list", length = ncol(dists2))
    
    i <- 1
    for (i in 1:length(neighbors_list)) {
      
      neighbors_list[[i]] <- order(dists2[,i])[1:n_closest_neighbors]
      
    }
    
    closest_neighbors <- unlist(neighbors_list)
    
    kept_models <- c(kept_models, closest_neighbors)
    
    kept_models <- kept_models[-which(duplicated(kept_models))]
    
  }
  
  
  ref_to_target_dist.df_filtered <- ref_to_target_dist.df[kept_models,]
  
  reference_expression_filtered <- reference_expression[ref_to_target_dist.df_filtered[,"Ref_Model_ID"],]
  

  
  ####density
  
  density.kept_ref.target.plot <- ggplot() + geom_density(data = as.data.frame(target_expression), aes(A, color = "A"), linetype = "dashed") +
    geom_density(data = as.data.frame(target_expression), aes(B, color = "B"), linetype = "dashed") +
    geom_density(data = as.data.frame(reference_expression_filtered), aes(A, color = "C")) +
    geom_density(data = as.data.frame(reference_expression_filtered), aes(B, color = "D")) +
    labs(title="Expression Density")  +  theme(plot.title = element_text(hjust = 0.5)) +
    scale_colour_manual("Gene", values = c("blue", "red", "blue", "red"), labels = c(A = "Target(A)", B = "Target(B)", C = "Kept Reference(A)", D = "Kept Reference(B)"))
  
  
    
  ###determinig which kept models belong to which cluster if weighting factors present
  
  if (!is.null(weighting_factors)) {
  
  weighted_IDs <- knn(train = centers, test = reference_expression_filtered, cl = colors, k = 1)
  kept_expression_with_clusters <- data.frame(reference_expression_filtered, colors = weighted_IDs)
  
  kept_ref.target_overlay <- ggplot() + 
    geom_point(data = as.data.frame(reference_expression), aes(A, B, color = "DF1"), shape=1) +
    geom_point( data = as.data.frame(kept_expression_with_clusters), aes(A, B, color = colors), shape=16) +
    geom_point( data = as.data.frame(target_expression), aes(A, B, color = "DF2"), shape = 16, alpha = 2) +
    geom_point(data = as.data.frame(centers), aes(A, B, color="DF3"), shape=13, size=7, alpha = 3) + 
    labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
    scale_colour_manual(values = c(DF1 = "red", DF2 = "black", DF3 = "blue", levels(colors)), 
                        labels = c(DF1 = "Current Reference", DF2 = "Target", DF3 = "Centers", yellow =  "Kept Reference Cluster 1", springgreen1 = "Kept Reference Cluster 2"),
                        name = "Data") + 
    guides(color = guide_legend(override.aes = list(shape = c(1, 16, 13, 16, 16), size = c(1.5, 1.5, 7, 1.5, 1.5)), title.hjust = .5))
  
  if (!is.null(reference_object$resampled_expression)) {
    
    total_overlay <- ggplot() + 
      geom_point(data = as.data.frame(reference_object$expression_ref), aes(A, B, color = "DF1"), shape = 16, alpha = 1) +
      geom_point( data = as.data.frame(reference_expression), aes(A, B, color = "DF2"), shape=1, alpha = 2) +
      geom_point( data = as.data.frame(kept_expression_with_clusters), aes(A, B, color = colors), shape = 16, alpha = 3) +
      geom_point( data = as.data.frame(target_expression), aes(A, B, color = "DF3"), shape = 16, alpha = 4) +
      geom_point(data = as.data.frame(centers), aes(A, B, color="DF4"), shape=13, size=7, alpha = 5) + 
      labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
      scale_colour_manual(values = c(DF1 = "cadetblue3", DF2 = " red", DF3 = "black", DF4 = "blue", levels(colors)), 
                          labels = c(DF1 = "Original Reference", DF2 = "Current Reference", DF3 = "Kept Target", DF4 = "Centers",
                                     yellow =  "Kept Reference Cluster 1", springgreen1 = "Kept Reference Cluster 2"), name = "Data") + 
      guides(color = guide_legend(override.aes = list(shape = c(16, 1, 16, 13, 16,16), size = c(1.5, 1.5, 1.5, 7, 1.5, 1.5)), title.hjust = .5)) 
    
  }
  
  
  } else {
    
    #resampled_ref.target_overlay
    kept_ref.target_overlay <- ggplot() + 
      geom_point(data = as.data.frame(reference_expression), aes(A, B, color = "DF1"), shape=1) +
      geom_point( data = as.data.frame(reference_expression_filtered), aes(A, B, color = "DF2"), shape=16) +
      geom_point( data = as.data.frame(target_expression), aes(A, B, color = "DF3"), shape = 16, alpha = 2) +
      geom_point(data = as.data.frame(centers), aes(A, B, color="DF4"), shape=13, size=7, alpha = 3) + 
      labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
      scale_colour_manual(values = c(DF1 = "red", DF2 = "yellow", DF3 = "black", DF4 = "green4"), 
                          labels = c(DF1 = "Current Reference", DF2 = "Kept Reference", DF3 = "Target", DF4 = "Centers"),
                          name = "Data") + 
      guides(color = guide_legend(override.aes = list(shape = c(1, 16, 16, 13), size = c(1.5, 1.5, 1.5, 7)), title.hjust = .5))
    
    
    if (!is.null(reference_object$resampled_expression)) {
      
      total_overlay <- ggplot() + 
        geom_point(data = as.data.frame(reference_object$expression_ref), aes(A, B, color = "DF1"), shape = 16, alpha = 1) +
        geom_point( data = as.data.frame(reference_expression), aes(A, B, color = "DF2"), shape=1, alpha = 2) +
        geom_point( data = as.data.frame(reference_expression_filtered), aes(A, B, color = "DF3"), shape = 16, alpha = 3) +
        geom_point( data = as.data.frame(target_expression), aes(A, B, color = "DF4"), shape = 16, alpha = 4) +
        geom_point(data = as.data.frame(centers), aes(A, B, color="DF5"), shape=13, size=7, alpha = 5) + 
        labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
        scale_colour_manual(values = c(DF1 = "cadetblue3", DF2 = " red", DF3 = "yellow", DF4 = "black", DF5 = "green4"), 
                            labels = c(DF1 = "Original Reference", DF2 = "Current Reference", DF3 = "Kept Reference", DF4 = "Target", DF5 = "Centers"),
                            name = "Data") + 
        guides(color = guide_legend(override.aes = list(shape = c(16, 1, 16, 16, 13), size = c(1.5, 1.5, 1.5, 1.5, 7)), title.hjust = .5)) 
      
    }
  }

  
  

  
  ###For 2 genes
  
  #Distances Distribution
  # Distance_density_plot <- ggplot() +geom_density(data = distance_df, aes(Nearest_Ref_Neighbor_Euc)) + 
  #   labs(title="Distance to Nearest Reference Neighbor") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Euclidean Distance")
  
  
  #geom_vline(xintercept = cutoff, color = "red") +
  
  if (!is.null(reference_object$resampled_reference_rset)) {
    
    
    filtered_params_df <- as.data.frame(sracipeParams(reference_object$resampled_reference_rset)[kept_models,])
    
  } else {
    
    filtered_params_df <- as.data.frame(sracipeParams(reference_object$reference_rset)[kept_models,])
    
  }
  
  
  
  output <- list(cutoff = cutoff, reference_expression = reference_expression, distance_df = distance_df, 
                 distribution_target_plot = distribution_target_plot, ref.target_overlay = ref.target_overlay, ref_to_target_dist.df = ref_to_target_dist.df, 
                 ref_to_target_dist.df_filtered = ref_to_target_dist.df_filtered, reference_expression_filtered = reference_expression_filtered, 
                 distribution_reference_plot = distribution_reference_plot, kept_ref.target_overlay = kept_ref.target_overlay, filtered_params_df = filtered_params_df,
                 density.ref.target.plot = density.ref.target.plot, density.kept_ref.target.plot = density.kept_ref.target.plot)
  
  if (!is.null(reference_object$resampled_expression)) {
    
    output$resampled_overlay <- resampled_overlay
    output$total_overlay <- total_overlay
  }
  
  if (!is.null(weighting_factors)) {
    
    output$weighted_IDs <- weighted_IDs 
    output$weightings <- sample_vec
    
  }
  
  output
}

##############################resample_reference#########################################################################################################################
resample_reference <- function(reference_object, identified_reference_parameters, weighted_IDs = NULL, weightings) {
  
  topology <- reference_object$topology
  numModels <- nrow(reference_object$expression_ref)
  
  resampled_reference_rset <- sracipeSimulate(topology, integrate = FALSE,
                                              numModels = numModels, genParams = TRUE,
                                              integrateStepSize = 0.02)
  
  resampled_reference_params_df <- sracipeParams(resampled_reference_rset) 
  
  if (!is.null(weighted_IDs)) {
    
    levels <- levels(weighted_IDs)
    
    weighted_IDs_list <- vector(mode = "list", length = length(levels))
    
    for (i in 1:length(weighted_IDs_list)) {
      
      weighted_IDs_list[[i]] <- which(weighted_IDs == levels[[i]])
      
    }
  

  ###
  
  
  for (i in 1:nparams) {
    
    samples <- c()
    
    for (j in 1:length(weighted_IDs_list)) {
      
      samples <- c(samples, sample(identified_reference_parameters[,i][weighted_IDs_list[[j]]], weightings[j], replace=TRUE))
      
      
      
    }
    
    resampled_reference_params_df[, i] <- samples
    
  }
  } else {
    
    nparams <- ncol(identified_reference_parameters)
    param_names <- colnames(resampled_reference_params_df)
    hill_params <- grep(pattern = "N", param_names)
    
    nmodels_resampled <- nrow(resampled_reference_params_df)
    nmodels_filtered <- nrow(identified_reference_parameters)
    
    ###
    
    for (i in 1:nparams) {
      
      resampled_reference_params_df[, i] <- sample(identified_reference_parameters[,i], nmodels_resampled, replace=TRUE)
      
    } 
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
  
  reference_object$resampled_reference_params.df <- sracipeParams(resampled_reference_rset)
  
  reference_object
}

##############################extract_parameters#########################################################################################################################
extract_parameters2 <- function(identified_reference_parameters, reference_object, nbins = 10, cutoff = 0.8) {
  
  nparams <- ncol(identified_reference_parameters)
  param_names <- colnames(sracipeParams(reference_object$reference_rset))
  hill_params <- grep(pattern = "N", param_names)
  
  ###
  uniform_counts <- nrow(identified_reference_parameters)/nbins
  
  
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
      
      uniform_counts_hill <- nrow(identified_reference_parameters)/length(reference_object$reference_parameter_ranges[[i]][1]:reference_object$reference_parameter_ranges[[i]][2])
      
      table <- table(identified_reference_parameters[,i])
      
      counts_tables[[i]] <- table
      
      ###
      
      cutoff_tables[[i]] <- table/uniform_counts_hill
      
    } else {
      
      bin_width <- (reference_object$reference_parameter_ranges[[i]][2] - reference_object$reference_parameter_ranges[[i]][1])/nbins
      
      breaks_vec <- seq(reference_object$reference_parameter_ranges[[i]][1], reference_object$reference_parameter_ranges[[i]][2], by = bin_width)
      
      table <- table(cut(identified_reference_parameters[,i], breaks = breaks_vec, include.lowest = TRUE))
      
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

extract_parameters <- function(parameters.df, reference_object, nbins = 10) {
  
  nparams <- ncol(parameters.df)
  param_names <- colnames(sracipeParams(reference_object$reference_rset))
  hill_params <- grep(pattern = "N", param_names)
  
  nmodels <- nrow(parameters.df)
  ###
  uniform.models_per_bin <- nmodels/nbins
  
  
  counts_tables <- vector(mode = "list", length = nparams)
  names(counts_tables) <- param_names
  
  percent_tables <- vector(mode = "list", length = nparams)
  names(percent_tables) <- param_names
  
  # cutoff_tables <- vector(mode = "list", length = nparams)
  # names(cutoff_tables) <- param_names
  
  plots_list <- vector(mode = "list", length = nparams)
  
  extracted_param_ranges <- vector(mode = "list", length = nparams)
  names(extracted_param_ranges) <- param_names
  ###
  
  for (i in 1:nparams) {
    
    if (i %in% hill_params) {
      
      table <- table(parameters.df[,i])
      
      counts_tables[[i]] <- table
      
      percent_tables[[i]] <- 100*table/nmodels
      
      # uniform.models_per_bin.hill <- nmodels/length(reference_object$reference_parameter_ranges[[i]][1]:reference_object$reference_parameter_ranges[[i]][2])
      # cutoff_tables[[i]] <- table/uniform_counts_hill
      
    } else {
      
      bin_width <- (reference_object$reference_parameter_ranges[[i]][2] - reference_object$reference_parameter_ranges[[i]][1])/nbins
      
      breaks_vec <- seq(reference_object$reference_parameter_ranges[[i]][1], reference_object$reference_parameter_ranges[[i]][2], by = bin_width)
      
      table <- table(cut(parameters.df[,i], breaks = breaks_vec, include.lowest = TRUE))
      
      counts_tables[[i]] <- table
      percent_tables[[i]] <- 100*table/nmodels
      
      ###
      
      # cutoff_tables[[i]] <- table/uniform.models_per_bin
      
    }
    
  }
  
  
  i <- 1
  
  for (i in 1:length(plots_list)) {
    
    if (i %in% hill_params) {
      
      nbins_hill <- length(reference_object$reference_parameter_ranges[[i]][1]:reference_object$reference_parameter_ranges[[i]][2])
      
      plots_list[[i]] <- ggplot(as.data.frame(percent_tables[[i]]), aes(x = Var1, weight = Freq)) + geom_bar() + geom_hline(yintercept = 100/nbins_hill, color = "red") +
        ggtitle(param_names[i]) + xlab("Param Value") +ylab("Ratio")
      
    } else {
      
      plots_list[[i]] <- ggplot(as.data.frame(percent_tables[[i]]), aes(x = Var1, weight = Freq)) + geom_bar() + geom_hline(yintercept = 100/nbins, color = "red") +
        ggtitle(param_names[i]) + xlab("Param Value") +ylab("Ratio")
      
    }
    
    # plots_list[[i]] <- ggplot(as.data.frame(cutoff_tables[[i]]), aes(x = Var1, weight = Freq)) + geom_bar() + geom_hline(yintercept = cutoff, color = "red") +
    #   ggtitle(param_names[i]) + xlab("Param Value") +ylab("Ratio")
    
  }
  
  
  
  output <- list(counts_tables = counts_tables, percent_tables = percent_tables, plots_list = plots_list)
  
  
}



##############################calc_ks###############################################################################################################
calc_ks <- function(ref_expression, target_expression) {
  
  #score_vec
  score_vec <- vector(length = ncol(ref_expression))
  
  #calc ks
  
  for (w in 1:length(score_vec)) {
    
    ks_test <- dgof::ks.test(x = target_expression[,w], y = ref_expression[,w])
    score_vec[w] <- as.numeric(ks_test$statistic)
    
  }
  
  
  scores_object <- list(score_vec = score_vec, score_vec_sum = sum(score_vec))
}
##############################optimize_parameters###############################################################################################################

optimize_parameters <- function(reference_object, target_object, num_iterations = 5, nth_neighbor = 20, cutoff_percent = 0.85, nbins = 10, n_closest_neighbors = NULL, 
                                weighting_factors = FALSE, center_type = "reference", num_clust = 2) {
  
  #Param info
  param_names <- colnames(sracipeParams(reference_object$reference_rset))
  nparams <- length(param_names)
  
  n_ref_models <- nrow(reference_object$expression_ref)
  
  #finding centers, defining colors, determing weighting factors
  
  if (center_type == "reference") {
    
    #centers
    kmeans <- kmeans(reference_object$expression_ref, centers = num_clust, nstart = 25)
    centers <- kmeans$centers
    
    #colors
    colors <- c("yellow", "springgreen1", "turqiose1", "purple")[1:num_clust]
    colors <- factor(colors, levels = colors)
        
    #weighting factors
    wf_vec <- c()
    sample_vec <- c()
    
    knn <- knn(train = centers, test = target_object$simulated_rset_expr_log_z, cl = colors, k = 1)

    knn_freq <- table(knn)

    for (i in 1:num_clust) {
      
      wf_vec <- c(wf_vec, knn_freq[i]/sum(knn_freq))
      
      sample_vec <- c(sample_vec, floor(wf_vec[i]*n_ref_models))
  
    }
    names(wf_vec) <- colors
    
    sample_vec[num_clust] <- n_ref_models - sum(sample_vec[1:(num_clust-1)])
    names(sample_vec) <- colors
    
  } else if (center_type == "target") {
    
    kmeans <- kmeans(target_object$simulated_rset_expr_log_z, centers = num_clust, nstart = 25)
    centers <- kmeans$centers
    
    #colors
    colors <- c("yellow", "springgreen1", "turqiose1", "purple")[1:num_clust]
    colors <- factor(colors, levels = colors)
    
    #weighting factors
    wf_vec <- c()
    sample_vec <- c()
    
    knn <- knn(train = centers, test = target_object$simulated_rset_expr_log_z, cl = colors, k = 1)
    knn_freq <- table(knn)
    for (i in 1:num_clust) {
      
      wf_vec <- c(wf_vec, knn_freq[i]/sum(knn_freq))
      
      sample_vec <- c(sample_vec, floor(wf_vec[i]*n_ref_models))

    }
    names(wf_vec) <- colors
    
    sample_vec[num_clust] <- n_ref_models - sum(sample_vec[1:(num_clust-1)])
    names(sample_vec) <- colors
  }
  

  
  ###
  identified_ref_models_list <- vector(mode = "list", length = num_iterations)
  
  extracted_params_list <- vector(mode = "list", length = num_iterations)
  extracted_ref_params_list <- vector(mode = "list", length = (num_iterations + 1))
  extracted_ref_params_list[[1]] <- extract_parameters(parameters.df = reference_object$reference_params.df, reference_object = reference_object, nbins = nbins)
  
  reference_objects_list <- vector(mode = "list", length = num_iterations)
  
  ##ks dynamics
  
  dynamics.ref.ks <- vector(length = (num_iterations + 1))
  dynamics.ref.ks[1] <- calc_ks(ref_expression = reference_object$expression_ref, target_expression = target_object$simulated_rset_expr_log_z)$score_vec_sum
  
  dynamics.kept_ref.ks <- vector(length = num_iterations)
  
  ##model number dynamics
  dynamics.frac_models_kept <- vector(length = num_iterations)
  
  #cluster ratio dynamics
  dynamics.cluster_ratio <- vector(length = num_iterations)
  
  dynamics.cluster_ratio_ref <- vector(length = (num_iterations + 1))
  
  knn <- knn(train = centers, test = reference_object$expression_ref, cl = colors, k = 1)
  knn_freq <- table(knn)
  ratio <- 100*knn_freq[1]/knn_freq[2]
  
  dynamics.cluster_ratio_ref[1] <- ratio
  
  
  #reference parameters list
  ref_params_list <- vector(mode = "list", length = (num_iterations + 1))
  ref_params_list[[1]] <- reference_object$reference_params.df
  ##
  
  for (i in 1:num_iterations) {
    
    reference_objects_list[[i]] <- reference_object
    
    ##################identify_reference_models#########
    
    if (weighting_factors == TRUE) {
      
      identified_ref_models <- identify_reference_models(target_expression = target_object$simulated_rset_expr_log_z,
                                                         reference_object = reference_object, nth_neighbor = nth_neighbor, cutoff_percent = cutoff_percent, 
                                                         n_closest_neighbors = n_closest_neighbors, centers = centers, weighting_factors = wf_vec )
      
    } else {
    
    identified_ref_models <- identify_reference_models(target_expression = target_object$simulated_rset_expr_log_z,
                                                       reference_object = reference_object, nth_neighbor = nth_neighbor, cutoff_percent = cutoff_percent, 
                                                       n_closest_neighbors = n_closest_neighbors, centers = centers)
    
    }
    
    identified_ref_models_list[[i]] <- identified_ref_models
    
    
    dynamics.kept_ref.ks[i] <- calc_ks(ref_expression = identified_ref_models$reference_expression_filtered,
                                       target_expression = target_object$simulated_rset_expr_log_z)$score_vec_sum 
    
    
    dynamics.frac_models_kept[i] <- nrow(identified_ref_models$reference_expression_filtered)/nrow(identified_ref_models$reference_expression)
    
    
    #cluster ratio dynamics kept models
    knn <- knn(train = centers, test = identified_ref_models$reference_expression_filtered, cl = factor(c(1,2)), k = 1)
    knn_freq <- table(knn)
    ratio <- 100*knn_freq[1]/knn_freq[2]
    
    dynamics.cluster_ratio[i] <- ratio
    
    
    ###################extract_parameters#########
    
    extracted_parameters <- extract_parameters(parameters.df = identified_ref_models$filtered_params_df, reference_object = reference_object, nbins = nbins)
    
    
    # extracted_parameters <- extract_parameters2(identified_reference_parameters = identified_ref_models$filtered_params_df, reference_object = reference_object,
    #                                            nbins = nbins, cutoff = 1)
    
    extracted_params_list[[i]] <- extracted_parameters
    
    ###################resampled_reference_object#########
    
    if (weighting_factors == TRUE) { 
      
    reference_object <- resample_reference(reference_object = reference_object, identified_reference_parameters = identified_ref_models$filtered_params_df, 
                                           weighted_IDs = reference_object$weighted_IDs, weightings = reference_object$weightings)
    
    } else {
      
      reference_object <- resample_reference(reference_object = reference_object, identified_reference_parameters = identified_ref_models$filtered_params_df)
      
    }
    
    dynamics.ref.ks[i + 1] <- calc_ks(ref_expression = reference_object$resampled_expression, target_expression = target_object$simulated_rset_expr_log_z)$score_vec_sum
    
    extracted_ref_params_list[[i + 1]] <- extract_parameters(parameters.df = reference_object$resampled_reference_params.df, reference_object = reference_object, nbins = nbins)
    
    ref_params_list[[i + 1]] <- reference_object$resampled_reference_params.df
    
    
    #cluster ratio dynamics reference models
    knn <- knn(train = centers, test = reference_object$resampled_expression, cl = factor(c(1,2)), k = 1)
    knn_freq <- table(knn)
    ratio <- 100*knn_freq[1]/knn_freq[2]
    
    dynamics.cluster_ratio_ref[i + 1] <- ratio
    
  }
  
  
  #dynamics
  
  df <- data.frame(Iteration = paste0("I", 0:num_iterations),
                   value=dynamics.ref.ks)
  
  df$Iteration <- factor(df$Iteration, levels = df$Iteration)
  
  dynamics.ref.ks_plot <- ggplot(data=df, aes(x=Iteration, y=value, group=1)) + geom_line() + geom_point() + 
    labs(title="KS Score of Target vs Reference ", x ="Iteration", y = "KS Score")
  
  
  ###################
  
  df <- data.frame(Iteration = paste0("I", 1:num_iterations),
                   value=dynamics.kept_ref.ks)
  
  df$Iteration <- factor(df$Iteration, levels = df$Iteration)
  
  dynamics.kept_ref.ks_plot <- ggplot(data=df, aes(x=Iteration, y=value, group=1)) + geom_line() + geom_point() + labs(title="KS Score of Target vs Kept Reference ",
                                                                                                                       x ="Iteration", y = "KS Score")
  
  ###################
  
  df <- data.frame(Iteration = paste0("I", 1:num_iterations),
                   value=dynamics.frac_models_kept)
  
  df$Iteration <- factor(df$Iteration, levels = df$Iteration)
  
  dynamics.frac_models_kept_plot <- ggplot(data=df, aes(x=Iteration, y=value, group=1)) + geom_line() + geom_point() + labs(title="Fraction of Reference Models Kept",
                                                                                                                            x ="Iteration", y = "Fraction of Reference Models")
  ####
  dynamics_plots_list <- list(dynamics.ref.ks_plot = dynamics.ref.ks_plot, dynamics.kept_ref.ks_plot = dynamics.kept_ref.ks_plot,
                              dynamics.frac_models_kept_plot = dynamics.frac_models_kept_plot)
  
  
  #####kept_reference_params dynamics
  param_dynamics_dfs <- vector(mode = "list", length = nparams)
  names(param_dynamics_dfs) <- param_names
  
  param_dynamics_plots <- vector(mode = "list", length = nparams)
  names(param_dynamics_plots) <- param_names
  
  i <- 1
  
  for (i in 1:nparams) {
    
    dfs_list <- list()
    
    for (j in 1:num_iterations) {
      
      df <- as.data.frame(extracted_params_list[[j]]$percent_tables[[i]])
      df$Iteration <- paste0("I",j)
      colnames(df) <- c("Bin", "Percent_of_Models", "Iteration")
      df$Bin <- as.character(df$Bin)
      
      df <- df[, c("Bin", "Iteration", "Percent_of_Models")]
      
      dfs_list <- append(dfs_list, list(df))
      
    }
    
    param_dynamics_dfs[[i]] <- do.call("rbind", dfs_list)
    param_dynamics_dfs[[i]]$Iteration <- factor(param_dynamics_dfs[[i]]$Iteration, levels = unique(param_dynamics_dfs[[i]]$Iteration))
    
    param_dynamics_plots[[i]] <- ggplot(param_dynamics_dfs[[i]], aes(x=Iteration, y= Percent_of_Models, group=Bin)) +
      geom_line(aes(color=Bin)) + geom_point(aes(color=Bin)) + 
      ggtitle(paste0("Parameter: ", param_names[i])) + xlab("Iteration") + ylab("Percent of Kept Models")
    
  }
  
  
  #####reference_params dynamics
  ref_param_dynamics_dfs <- vector(mode = "list", length = nparams)
  names(ref_param_dynamics_dfs) <- param_names
  
  ref_param_dynamics_plots <- vector(mode = "list", length = nparams)
  names(ref_param_dynamics_plots) <- param_names
  
  i <- 1
  
  for (i in 1:nparams) {
    
    dfs_list <- list()
    
    for (j in 1:length(extracted_ref_params_list)) {
      
      df <- as.data.frame(extracted_ref_params_list[[j]]$percent_tables[[i]])
      df$Iteration <- paste0("I",(j-1))
      colnames(df) <- c("Bin", "Percent_of_Models", "Iteration")
      df$Bin <- as.character(df$Bin)
      
      df <- df[, c("Bin", "Iteration", "Percent_of_Models")]
      
      dfs_list <- append(dfs_list, list(df))
      
    }
    
    ref_param_dynamics_dfs[[i]] <- do.call("rbind", dfs_list)
    ref_param_dynamics_dfs[[i]]$Iteration <- factor(ref_param_dynamics_dfs[[i]]$Iteration, levels = unique(ref_param_dynamics_dfs[[i]]$Iteration))
    
    ref_param_dynamics_plots[[i]] <- ggplot(ref_param_dynamics_dfs[[i]], aes(x=Iteration, y= Percent_of_Models, group=Bin)) +
      geom_line(aes(color=Bin)) + geom_point(aes(color=Bin)) + 
      ggtitle(paste0("Parameter: ", param_names[i])) + xlab("Iteration") + ylab("Percent of Reference Models")
    
  }
  
  #####parameter correlation
  corr_dynamics_dfs <- vector(mode = "list", length = nparams)
  names(corr_dynamics_dfs) <- param_names
  
  corr_dynamics_plots <- vector(mode = "list", length = nparams)
  names(corr_dynamics_plots) <- param_names
  
  i <- 1
  
  for (i in 1:nparams) {
    
    dfs_list <- list()
    
    
    for (j in 1:num_iterations) {
      
      corr <- cor(identified_ref_models_list[[j]]$filtered_params_df)
      
      df <- data.frame(Iteration = rep(paste0("I",j), nparams), CoCof = corr[,i], Parameter = param_names)
      
      df <- df[-i,]
      
      dfs_list <- append(dfs_list, list(df))
      
    }
    
    corr_dynamics_dfs[[i]] <- do.call("rbind", dfs_list)
    corr_dynamics_dfs[[i]]$Iteration <- factor(corr_dynamics_dfs[[i]]$Iteration, levels = unique(corr_dynamics_dfs[[i]]$Iteration))
    
    corr_dynamics_plots[[i]] <- ggplot(corr_dynamics_dfs[[i]], aes(x=Iteration, y= CoCof, group=Parameter)) +
      geom_line(aes(color=Parameter)) + geom_point(aes(color=Parameter)) + 
      ggtitle(paste0("Parameter: ", param_names[i])) + xlab("Iteration") + ylab("Correlation Coefficient") + ylim(-0.25, 0.25)
    
  }
  
  
  #####parameter correlation resampled reference
  corr_ref_dynamics_dfs <- vector(mode = "list", length = nparams)
  names(corr_ref_dynamics_dfs) <- param_names
  
  corr_ref_dynamics_plots <- vector(mode = "list", length = nparams)
  names(corr_ref_dynamics_plots) <- param_names
  
  i <- 1
  
  for (i in 1:nparams) {
    
    dfs_list <- list()
    
    
    for (j in 1:length(ref_params_list)) {
      
      corr <- cor(ref_params_list[[j]])
      
      df <- data.frame(Iteration = rep(paste0("I",(j-1)), nparams), CoCof = corr[,i], Parameter = param_names)
      
      df <- df[-i,]
      
      dfs_list <- append(dfs_list, list(df))
      
    }
    
    corr_ref_dynamics_dfs[[i]] <- do.call("rbind", dfs_list)
    corr_ref_dynamics_dfs[[i]]$Iteration <- factor(corr_ref_dynamics_dfs[[i]]$Iteration, levels = unique(corr_ref_dynamics_dfs[[i]]$Iteration))
    
    corr_ref_dynamics_plots[[i]] <- ggplot(corr_ref_dynamics_dfs[[i]], aes(x=Iteration, y= CoCof, group=Parameter)) +
      geom_line(aes(color=Parameter)) + geom_point(aes(color=Parameter)) + 
      ggtitle(paste0("Parameter: ", param_names[i])) + xlab("Iteration") + ylab("Correlation Coefficient") + ylim(-0.04, 0.04)
    
  }
  
  
  
  
  ###cluster ratio dynamics plots
  knn_target <- knn(train = centers, test = target_object$simulated_rset_expr_log_z, cl = factor(c(1,2)), k = 1)
  knn_target_freq <- table(knn_target)
  target_ratio <- 100*knn_target_freq[1]/knn_target_freq[2]
  
  #kept
  df <- data.frame(Iteration = paste0("I", 1:num_iterations),
                   value=dynamics.cluster_ratio)
  
  df$Iteration <- factor(df$Iteration, levels = df$Iteration)
  
  dynamics.cluster_ratio.plot <- ggplot(data=df, aes(x=Iteration, y=value,group=1)) + geom_line() + geom_point() + 
    labs(title="Cluster Ratio of Kept Models", x ="Iteration", y = "Cluster Ratio") + geom_hline(yintercept=target_ratio, color = "red", size=1)
  
  #ref
  df <- data.frame(Iteration = paste0("I", 0:num_iterations),
                   value=dynamics.cluster_ratio_ref)
  
  df$Iteration <- factor(df$Iteration, levels = df$Iteration)
  
  dynamics.cluster_ratio.ref.plot <- ggplot(data=df, aes(x=Iteration, y=value,group=1)) + geom_line() + geom_point() + 
    labs(title="Cluster Ratio of Reference", x ="Iteration", y = "Cluster Ratio") + geom_hline(yintercept=target_ratio, color = "red", size=1)
  
  
  
  
  dynamics_params <- list(param_dynamics_dfs = param_dynamics_dfs, param_dynamics_plots = param_dynamics_plots, ref_param_dynamics_dfs = ref_param_dynamics_dfs, 
                          ref_param_dynamics_plots = ref_param_dynamics_plots, corr_dynamics_dfs = corr_dynamics_dfs, corr_dynamics_plots = corr_dynamics_plots, 
                          corr_ref_dynamics_dfs = corr_ref_dynamics_dfs, corr_ref_dynamics_plots = corr_ref_dynamics_plots)
  
  
  dynamics <- list(dynamics.ref.ks= dynamics.ref.ks, dynamics.kept_ref.ks = dynamics.kept_ref.ks, dynamics.frac_models_kept = dynamics.frac_models_kept, 
                   dynamics_plots_list = dynamics_plots_list, dynamics_params = dynamics_params, 
                   dynamics.cluster_ratio.plot = dynamics.cluster_ratio.plot, dynamics.cluster_ratio.ref.plot = dynamics.cluster_ratio.ref.plot)
  
  
  
  
  
  
  results <- list(identified_ref_models_list = identified_ref_models_list, extracted_params_list = extracted_params_list, reference_objects_list = reference_objects_list, 
                  dynamics = dynamics)
  
  
}



