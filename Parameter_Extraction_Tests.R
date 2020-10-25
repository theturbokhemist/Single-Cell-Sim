topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))

reference_parameter_ranges <- list(c(1,100), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100))

target_parameter_ranges_TS.G_A.l <- list(c(1,10), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100))


target_parameter_ranges_TS.G_A.m <- list(c(45,55), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100))

target_parameter_ranges_TS_G_A.h_K_A.l <- list(c(80,100), c(1,100), c(0.1,0.3), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100))

#generate ref
reference_object_TS.30k <- generate_reference_object(topology = topology_TS, reference_parameter_ranges = reference_parameter_ranges, numModels = 30000)
reference_object_TS.10k <- generate_reference_object(topology = topology_TS, reference_parameter_ranges = reference_parameter_ranges, numModels = 10000)

standard_expression_TS.kmeans <- kmeans(reference_object_TS.30k$expression_ref, centers = 2, nstart = 25)
standard_expression_TS.kmeans$cluster


target_object_TS_G_A.h_K_A.l_kmeans <- kmeans(target_object_TS_G_A.h_K_A.l$simulated_rset_expr_log_z, centers = standard_expression_TS.kmeans$centers, nstart = 1)
target_object_TS_G_A.h_K_A.l_kmeans$centers

df <- target_object_TS_G_A.h_K_A.l$simulated_rset_expr_log_z
a <- knn(train = standard_expression_TS.kmeans$centers, test = target_object_TS_G_A.h_K_A.l$simulated_rset_expr_log_z, cl = factor(c("yellow","springgreen1")), k = 1)

df <- data.frame(df, Cluster = a)

df <- cbind(df, as.data.frame(a))
colnames(df) <- c("A", "B", "Cluster")  
df <- as.data.frame(df)

ggplot() + geom_point(data = df, aes(A, B, color = factor(Cluster)), shape=1)  + 
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color = "red"), shape=13, size=7, alpha = 5) +
  scale_color_manual(values=c('blue','black', "red"))+labs(color='Cluster')

a <- factor(c("yellow","springgreen1"), levels = c("yellow","springgreen1"))

levels(a)
which(a == factor(c("yellow","springgreen1"))[1])

df_main <- data.frame(reference_object_TS.30k$expression_ref, color = "cadetblue3", shape = "16", size = "1.5")

a <- knn(train = standard_expression_TS.kmeans$centers, test = target_object_TS_G_A.h_K_A.l$simulated_rset_expr_log_z, cl = factor(c("yellow","springgreen1")), k = 1)
df <- data.frame(target_object_TS_G_A.h_K_A.l$simulated_rset_expr_log_z, color = a, shape = "16", size = "1.5")

df_main <- rbind(df_main, df)

df <- data.frame(target_object_TS_G_A.h_K_A.l$simulated_rset_expr_log_z, color = "black", shape = "16", size = "1.5")
df_main <- rbind(df_main, df)

df <- data.frame(standard_expression_TS.kmeans$centers, color = "blue", shape = "13", size = "7")
df_main <- rbind(df_main, df)

ggplot(df_main, aes(A, B, color = color, shape = shape, size = size)) + geom_point() +
  scale_shape_manual(values=c(16, 16, 16, 16, 7)) + 
  scale_color_manual(values=c('cadetblue3','yellow', 'springgreen1', "black", "blue"))+
  theme(legend.position="top")

ggplot() + 
  geom_point(data = as.data.frame(reference_object_TS.30k$expression_ref), aes(A, B, color = "DF1"), shape = 16, alpha = 1) +
  geom_point( data = as.data.frame(test$reference_objects_list[[2]]$resampled_expression), aes(A, B, color = "DF2"), shape=1, alpha = 2) +
  geom_point(data = df, aes(A, B, color = factor(Cluster)), shape=16, alpha = 3) + 
  geom_point( data = as.data.frame(target_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B, color = "DF3"), shape = 16, alpha = 4) +
  # geom_point( data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color = "DF4"), shape = 16, alpha = 5) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF4"), shape=13, size=7, alpha = 5) +
  scale_color_manual(values=c("cadetblue3","red", "black","blue", "yellow", "springgreen1"),
                     labels = c(DF1 = "Reference", DF2 = "Resampled Reference", DF3 = "Target", DF4 = "Centers", a = "Kept Reference Cluster 1", Z= "Z"),
                     name = "Data") +
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 16, 13, 16, 16), size = c(1.5, 1.5, 1.5, 7, 1.5, 1.5)), title.hjust = .5))
  
#turqiose1 and and purple

colors <- c("yellow", "springgreen1", "turqiose1", "purple")[1:2]
colors <- factor(colors, levels = colors)

  
weighted_IDs <- knn(train = standard_expression_TS.kmeans$centers, test = test3$identified_ref_models_list[[1]]$reference_expression_filtered, cl = colors, k = 1)
kept_expression_with_clusters <- data.frame(test3$identified_ref_models_list[[1]]$reference_expression_filtered, colors = weighted_IDs)

ggplot() + 
  geom_point(data = as.data.frame(reference_object_TS.30k$expression_ref), aes(A, B, color = "DF1"), shape=1) +
  geom_point( data = as.data.frame(kept_expression_with_clusters), aes(A, B, color = colors), shape=16) +
  geom_point( data = as.data.frame(target_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B, color = "DF2"), shape = 16, alpha = 2) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=13, size=7, alpha = 3) + 
  labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "red", DF2 = "black", DF3 = "blue", levels(colors)), 
                      labels = c(DF1 = "Current Reference", DF2 = "Target", DF3 = "Centers", yellow =  "Kept Reference Cluster 1", springgreen1 = "Kept Reference Cluster 2"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(1, 16, 13, 16, 16), size = c(1.5, 1.5, 7, 1.5, 1.5)), title.hjust = .5))



ggplot(data = df, aes(A, B), shape=1) + geom_point(fill = Cluster)


ggplot() + geom_point(data = df, aes(A, B, fill = Cluster), shape=1, color = "NA")  + 
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="Red"), shape=13, size=7, alpha = 5)


geom_point( data = as.data.frame(target_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B, color = "DF3"), shape = 16, alpha = 2)

table(target_object_TS_G_A.h_K_A.l_kmeans$cluster)
table(df$Cluster)



#generate target

target_object_TS.G_A.l <- generate_target_object(target_parameter_ranges = target_parameter_ranges_TS.G_A.l, reference_object = reference_object_TS.30k,
                                                 numModels = 2000, save_rset = T)


target_object_TS.G_A.m <- generate_target_object(target_parameter_ranges = target_parameter_ranges_TS.G_A.m, reference_object = reference_object_TS.30k,
                                                 numModels = 1000, save_rset = T)


target_object_TS_G_A.h_K_A.l <- generate_target_object(target_parameter_ranges = target_parameter_ranges_TS_G_A.h_K_A.l, reference_object = reference_object_TS.30k,
                                                 numModels = 1000, save_rset = T)
#####correlation
cor(target_object_TS_G_A.h_K_A.l$target_params.df)

a <- data.frame(a = target_object_TS_G_A.h_K_A.l$target_params.df$G_A, b = target_object_TS_G_A.h_K_A.l$target_params.df$K_A)
  
  ggscatter(a, x = "a", y = "b", 
            add = "reg.line", conf.int = TRUE, 
            cor.coef = TRUE, cor.method = "pearson",
            xlab = "G_A", ylab = "K_A")

a <- data.frame(a = test$identified_ref_models_list[[1]]$filtered_params_df$G_A, b = test$identified_ref_models_list[[1]]$filtered_params_df$K_A)
b <- data.frame(a = test$identified_ref_models_list[[20]]$filtered_params_df$G_A, b = test$identified_ref_models_list[[20]]$filtered_params_df$K_A)

# library("ggpubr")
a <- ggscatter(a, x = "a", y = "b", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "G_A", ylab = "K_A")

b <- ggscatter(b, x = "a", y = "b", 
                add = "reg.line", conf.int = TRUE, 
                cor.coef = TRUE, cor.method = "pearson",
                xlab = "G_A", ylab = "K_A")

a <- list(a, b)

do.call("grid.arrange", a)


##################identify_reference_models#########

identified_ref_models <- identify_reference_models(target_expression = target_object_TS.G_A.l$simulated_rset_expr_log_z,
                                                   reference_object = reference_object_TS.30k, nth_neighbor = 20, cutoff_percent = .9, n_closest_neighbors = 10)

identified_ref_models$distribution_target_plot
identified_ref_models$kept_ref.target_overlay

identified_ref_models$distribution_reference_plot
do.call("grid.arrange", identified_ref_models$overlay_list_ref)

###################extract_parameters#########

extracted_parameters_test <- extract_parameters2(parameters.df = target_object_TS.G_A.l$target_params.df, reference_object = reference_object_TS.30k, 
                                                 nbins = 10)

# extracted_parameters <- extract_parameters(identified_reference_parameters = identified_ref_models$filtered_params_df, reference_object = reference_object_TS.30k,
#                                            nbins = 10, cutoff = 0.8)

df <- as.data.frame(extracted_parameters_test$percent_tables$N_B_A)

df$Iteration <- paste0("I",0)
do.call("grid.arrange", extracted_parameters_test$plots_list)


rbind(extracted_parameters_test$percent_tables$G_A, extracted_parameters_test$percent_tables$G_B)  
###########
extract_parameters2 <- function(parameters.df, reference_object, nbins = 10) {
  
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
  
  #extract param
  i <- 1
  for (i in 1:length(percent_tables)) {
    
    df <- as.data.frame(percent_tables[[i]])
    cutoff_rows <- which(df[,"Freq"] >= cutoff)
    
    df <- df[cutoff_rows,]
    
    ranges_as.char <- as.character(df[, "Var1"])
    
    extract_numerics <- as.numeric(unlist(regmatches(ranges_as.char,
                                                     gregexpr("[[:digit:]]+\\.*[[:digit:]]*",ranges_as.char))))
    
    param_range <- range(extract_numerics)
    
    extracted_param_ranges[[i]] <- param_range
    
  }
  
  
  # i <- 1
  # for (i in 1:length(cutoff_tables)) {
  #   
  #   df <- as.data.frame(cutoff_tables[[i]])
  #   cutoff_rows <- which(df[,"Freq"] >= cutoff)
  #   
  #   df <- df[cutoff_rows,]
  #   
  #   ranges_as.char <- as.character(df[, "Var1"])
  #   
  #   extract_numerics <- as.numeric(unlist(regmatches(ranges_as.char,
  #                                                    gregexpr("[[:digit:]]+\\.*[[:digit:]]*",ranges_as.char))))
  #   
  #   param_range <- range(extract_numerics)
  #   
  #   extracted_param_ranges[[i]] <- param_range
  #   
  # }
  # output <- list(counts_tables = counts_tables, cutoff_tables = cutoff_tables, plots_list = plots_list, extracted_param_ranges = extracted_param_ranges)
  # 
  
  
  output <- list(counts_tables = counts_tables, percent_tables = percent_tables, plots_list = plots_list)
  
  
}



##################resample_reference#########

resampled_reference_object <- resample_reference(reference_object = reference_object_TS.30k, identified_reference_parameters = identified_ref_models$filtered_params_df)

resampled_reference_object$resampled_reference_rset



##################identify_reference_models again#########

identified_ref_models2 <- identify_reference_models(target_expression = target_object_TS.G_A.l$simulated_rset_expr_log_z,
                                                    reference_object = resampled_reference_object, nth_neighbor = 20, cutoff_percent = .85)


identified_ref_models2$resampled_overlay

identified_ref_models2$

do.call("grid.arrange", identified_ref_models2$overlay_list_target)

do.call("grid.arrange", identified_ref_models2$overlay_list_ref)

dim(identified_ref_models2$filtered_params_df)
dim(identified_ref_models$filtered_params_df)

identified_ref_models2$

###################extract_parameters#########

extracted_parameters2 <- extract_parameters(identified_reference_parameters = identified_ref_models2$filtered_params_df, reference_object = resampled_reference_object,
                                            nbins = 10, cutoff = 0.8)

do.call("grid.arrange", extracted_parameters2$plots_list)



##################resample_reference again#########

resampled_reference_object2 <- resample_reference(reference_object = reference_object_TS.20k, identified_reference_parameters = identified_ref_models2$filtered_params_df)

resampled_reference_object2$resampled_expression



############
test <- optimize_parameters(reference_object = reference_object_TS.30k, target_object = target_object_TS_G_A.h_K_A.l, num_iterations = 15, cutoff_percent = .9,
                            nth_neighbor = 20, nbins = 10, n_closest_neighbors = 20, weighting_factors = FALSE, center_type = "target", num_clust = 2)

test3 <- optimize_parameters(reference_object = reference_object_TS.30k, target_object = target_object_TS_G_A.h_K_A.l, num_iterations = 15, cutoff_percent = .9,
                            nth_neighbor = 20, nbins = 10, n_closest_neighbors = 20, weighting_factors = TRUE, center_type = "target", num_clust = 2)

test3$identified_ref_models_list[[15]]$total_overlay

do.call("grid.arrange", test3$dynamics$dynamics_params$corr_dynamics_plots)
do.call("grid.arrange", test3$dynamics$dynamics_params$corr_ref_dynamics_plots)


test3$dynamics$dynamics.cluster_ratio.plot
test3$dynamics$dynamics.cluster_ratio.ref.plot
cor(test$reference_objects_list[[2]]$resampled_reference_params.df)

test$identified_ref_models_list[[1]]$ref.target_overlay
test$identified_ref_models_list[[1]]$distribution_target_plot

test$identified_ref_models_list[[1]]$distribution_reference_plot
test$identified_ref_models_list[[1]]$kept_ref.target_overlay

test$identified_ref_models_list[[2]]$resampled_overlay
test$identified_ref_models_list[[2]]$distribution_reference_plot
test$identified_ref_models_list[[2]]$total_overlay

test$identified_ref_models_list[[1]]$density.ref.target.plot
test$identified_ref_models_list[[1]]$density.kept_ref.target.plot

do.call("grid.arrange", test3$extracted_params_list[[5]]$plots_list)

test$identified_ref_models_list[[40]]$resampled_overlay
test$identified_ref_models_list[[40]]$distribution_reference_plot
test$identified_ref_models_list[[40]]$total_overlay

do.call("grid.arrange", test$extracted_params_list[[40]]$plots_list)

test$identified_ref_models_list[[1]]$density.ref.target.plot
test$identified_ref_models_list[[2]]$density.kept_ref.target.plot



test3$dynamics$dynamics_plots_list$dynamics.ref.ks_plot
test3$dynamics$dynamics_plots_list$dynamics.kept_ref.ks_plot
test$dynamics$dynamics_plots_list$dynamics.frac_models_kept_plot

do.call("grid.arrange", test$extracted_params_list[[3]]$plots_list)


do.call("grid.arrange", test$dynamics$dynamics_params$param_dynamics_plots)
do.call("grid.arrange", test$dynamics$dynamics_params$ref_param_dynamics_plots)

gtools::mixedorder(c("I0", "I11", "I1"))

ggplot(test$dynamics$param_dynamics_dfs[[1]], aes(x=Iteration, y= Percent_of_Models, group=Bin)) +
  geom_line(aes(color=Bin)) +
  geom_point(aes(color=Bin)) + ggtitle("Plot of length \n by dose") +
  xlab("Dose (mg)") + ylab("Teeth length")


df1 <- data.frame(iteration = paste0("I", 0:5),
                 value=test$dynamics$dynamics.ref.ks)

ggplot(data=df1, aes(x=iteration, y=value, group=1)) +
  geom_line() + geom_point()

ggplot() + geom_density(data = as.data.frame(test$identified_ref_models_list[[1]]$reference_expression_filtered), aes(A, color = "A"), linetype = "dashed") +
           geom_density(data = as.data.frame(test$identified_ref_models_list[[1]]$reference_expression_filtered), aes(B, color = "B")) +
           labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
           scale_colour_manual(values = c(A = "black", B = "red"), 
                      labels = c(A = "Reference", B = "Resampled Reference"),
                      name = "Data")

ggplot() + 
  geom_point(data = as.data.frame(reference_object_TS.30k$expression_ref), aes(A, B, color = "DF1"), shape = 1) +
  geom_point( data = as.data.frame(test$reference_objects_list[[2]]$resampled_expression), aes(A, B, color = "DF2"), shape=1) +
  geom_point( data = as.data.frame(target_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B, color = "DF3"), shape = 16, alpha = 2) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF4"), shape=13, size=7, alpha = 3) + 
  labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "black", DF2 = "light blue", DF3 = "red", DF4 = "purple"), 
                      labels = c(DF1 = "Reference", DF2 = "Resampled Reference", DF3 = "Target", DF4 = "Centers"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(1, 1, 16, 13), size = c(1.5, 1.5, 1.5, 7)), title.hjust = .5)) 
# + xlim(-3, 3) + ylim(-3, 3)
ggplot() + 
  geom_point(data = as.data.frame(reference_object_TS.30k$expression_ref), aes(A, B, color = "DF1"), shape = 16, alpha = 1) +
  geom_point( data = as.data.frame(test$reference_objects_list[[2]]$resampled_expression), aes(A, B, color = "DF2"), shape=1, alpha = 2) +
  geom_point( data = as.data.frame(test$identified_ref_models_list[[2]]$reference_expression_filtered), aes(A, B, color = "DF5"), shape = 16, alpha = 1) +
  geom_point( data = as.data.frame(target_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B, color = "DF3"), shape = 16, alpha = 4) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF4"), shape=13, size=7, alpha = 5) + 
  labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "cadetblue3", DF2 = " red", DF3 = "black", DF4 = "green4", DF5 = "yellow"), 
                      labels = c(DF1 = "Reference", DF2 = "Resampled Reference", DF3 = "Target", DF4 = "Centers", DF5 = "Kept Reference"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 16, 13, 16), size = c(1.5, 1.5, 1.5, 7, 1.5)), title.hjust = .5)) 



results.9_20 <- optimize_parameters(reference_object = reference_object_TS.30k, target_object = target_object_TS.G_A.l, num_iterations = 100, cutoff_percent = .9,
                                    nth_neighbor = 20)

results.9_20.G_A.m <- optimize_parameters(reference_object = reference_object_TS.30k, target_object = target_object_TS.G_A.m, num_iterations = 50, cutoff_percent = .9,
                                    nth_neighbor = 20, nbins = 20)


results$identified_ref_models_list[[1]]$overlay_target

results.9_20$identified_ref_models_list[[100]]$resampled_overlay

do.call("grid.arrange", results.9_20$identified_ref_models_list[[2]]$overlay_list_ref)

do.call("grid.arrange", results.9_20$extracted_params_list[[100]]$plots_list)

############################################
results$identified_ref_models_list[[1]]$distribution_target_plot

results.9_20.G_A.m$identified_ref_models_list[[50]]$resampled_overlay
results.9_20.G_A.m$identified_ref_models_list[[1]]$overlay_target
results.9_20.G_A.m$identified_ref_models_list[[1]]$overlay_target

do.call("grid.arrange", results.9_20.G_A.m$identified_ref_models_list[[50]]$overlay_list_ref)

do.call("grid.arrange", results.9_20.G_A.m$extracted_params_list[[50]]$plots_list)


ggplot() +  geom_point(data = as.data.frame(reference_object_TS.30k$expression_ref), aes(A, B))
  
############################################
results.95_20_G_A.h_K_A.l <- optimize_parameters(reference_object = reference_object_TS.30k, target_object = target_object_TS_G_A.h_K_A.l, num_iterations = 40, cutoff_percent = .95,
                                          nth_neighbor = 20, nbins = 20)

results.95_20_G_A.h_K_A.l$identified_ref_models_list[[40]]$resampled_overlay

results.95_20_G_A.h_K_A.l$identified_ref_models_list[[1]]$overlay_target
do.call("grid.arrange", results.95_20_G_A.h_K_A.l$identified_ref_models_list[[1]]$overlay_list_ref)

do.call("grid.arrange", results.95_20_G_A.h_K_A.l$extracted_params_list[[40]]$plots_list)


#ssh -Y d.gordin@login.discovery.neu.edu
