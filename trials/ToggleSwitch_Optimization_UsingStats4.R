Generate_results_object2 <- function(standard_object, experimental_expression) {
  
  nparams <- ncol(sracipeParams(standard_object$reference_rset))
  param_names <- colnames(sracipeParams(standard_object$reference_rset))
  
  standard_expression.kmeans <- kmeans(standard_object$expression_ref, centers = 2, nstart = 25)
  
  #Intialize distance_df
  ncol_df <- 2
  nrow_df <- nrow(experimental_expression)
  
  distance_df <- data.frame(matrix(ncol = ncol_df , nrow = nrow_df))
  colnames(distance_df) <- c("Nearest_Ref_Neighbor_Euc", "Model_ID")

  
  #Calc distances
  dists <- as.matrix(pdist(standard_object$expression_ref, experimental_expression))
  
  #Generate list of distances
  for (i in 1:nrow(distance_df)) {
    
    distance_df[i, 2] <- order(dists[,i])[1]
    
    distance_df[i, 1] <- dists[,i][distance_df[i, 2]]
    
  }
  

  ###For 2 genes
  
  #Distances Distribution
  # Distance_density_plot <- ggplot() +geom_density(data = distance_df, aes(Nearest_Ref_Neighbor_Euc)) + 
  #   labs(title="Distance to Nearest Reference Neighbor") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Euclidean Distance")
  
  
  #geom_vline(xintercept = cutoff, color = "red") +
  output <- distance_df
}


Distance_density_plot <- ggplot() +geom_histogram(data = test, aes(Nearest_Ref_Neighbor_Euc)) + 
  labs(title="Distance to Nearest Reference Neighbor") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Euclidean Distance")

test <- Generate_results_object2(standard_object = standard_object_TS_100k, experimental_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z)

test <- test[which(test[,1] <= 0.1),]
dim(test)


ggplot(data = test, aes(Nearest_Ref_Neighbor_Euc)) + geom_histogram(binwidth = 0.001) + ylim(0, 200)


nbins <- 100
breaks <- seq(0, 0.1, by = 0.1/nbins)

df1 <- as.data.frame(table(cut(test[,1], breaks = breaks, include.lowest = TRUE))) 
df1[, 1] <- breaks[2:(nbins + 1)]

ggplot(df1, aes(x = Var1, weight = Freq)) + geom_bar() + ylim(0, 200)

### 


test2 <- Generate_results_object2(standard_object = standard_object_TS_250k, experimental_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z)

test2 <- test2[which(test2[,1] <= 0.1),]
dim(test2)


ggplot(data = test2, aes(Nearest_Ref_Neighbor_Euc)) + geom_histogram(binwidth = 0.001) + ylim(0, 200)

nbins <- 100
breaks <- seq(0, 0.1, by = 0.1/nbins)

df2 <- as.data.frame(table(cut(test2[,1], breaks = breaks, include.lowest = TRUE)))
df2[, 1] <- breaks[2:(nbins + 1)]

ggplot(df2, aes(x = Var1, weight = Freq)) + geom_bar() + ylim(0, 200)

# ggplot(data = test2, aes(Nearest_Ref_Neighbor_Euc)) + geom_bar() +
#   scale_x_binned(n.breaks = 50) +
#   labs(title="Distance to Nearest Reference Neighbor") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Euclidean Distance")
# 
# Distance_density_plot
  
  filtered_models <- which(distance_df[,"Distance_nearest_exp_neighbor"] <= cutoff)
  filtered_params_df <- as.data.frame(sracipeParams(standard_object$reference_rset)[filtered_models,])
  
  #param plots n
  param_plots_list <- vector(mode = "list", length = nparams)
  names(param_plots_list) <- param_names
  
  for (j in 1:nparams) {
    
    df <- filtered_params_df
    colnames(df)[j] <- "placeholder"
    
    plot <- ggplot(df, aes(x = placeholder)) + geom_line(stat = "density")
    plot <- plot + labs(title=param_names[j], x="Param Values", y = "Density") + theme(plot.title = element_text(hjust = 0.5))
    
    
    param_plots_list[[j]] <- plot 
    
  }
  
  
  
  filtered_standard_expr <- as.data.frame(standard_object$expression_ref[filtered_models,])
  
  Expression_plot <- ggplot() + geom_line(data = filtered_standard_expr, aes(A, color = "A"), stat = "density") + 
    geom_line(data = filtered_standard_expr, aes(B, color = "B"), stat = "density") +
    scale_colour_manual(values=c("A"="blue", "B"="red"), name="Genes") + labs(title="Filtered Standard Expression") + 
    theme(plot.title = element_text(hjust = 0.5)) + ylab("Expression")
  
  Expression_overlay_plot <- ggplot() + 
    geom_point(data = as.data.frame(standard_object$expression_ref), aes(A, B, color = "DF1")) +
    geom_point( data = filtered_standard_expr, aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
    geom_point(data = as.data.frame(standard_expression.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
    labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
    scale_colour_manual(values = c(DF1 = "black", DF2 = "blue", DF3 = "red"), 
                        labels = c(DF1 = "Standard", DF2 = "Filtered Standard", DF3 = "Centers"),
                        name = "Data") + 
    guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))
  
  
  output <- list(std.vs.exp_distance_matrix = dists, shortest_distance_df = distance_df, filtered_standard_expr = filtered_standard_expr, filtered_models = filtered_models,
                 filtered_params_df = filtered_params_df, Distance_density_plot = Distance_density_plot, Expression_plot  = Expression_plot,
                 Expression_overlay_plot = Expression_overlay_plot, param_plots_list = param_plots_list)
  output
  
  
  
}


##########################################
extract_param_ranges <- function(results_object, standard_object, nbins = 10, cutoff = 0.9) {
  
  nparams <- ncol(results_object$filtered_params_df)
  param_names <- colnames(results_object$filtered_params_df)
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
      
      table <- table(Results_object_TS_G_A.l_co005$filtered_params_df[,i])
      
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
##########################################

a <- extract_param_ranges(results_object = Results_object_TS_G_A.l_co005, standard_object = standard_object_TS_100k, nbins = 10, cutoff = 0.77)
do.call("grid.arrange", a$plots_list)
a$extracted_param_ranges

ggplot(as.data.frame(a$cutoff_tables[[1]]), aes(x = Var1, weight = Freq)) + geom_bar() + geom_hline(yintercept = 0.9, color = "red")

table(Results_object_TS_G_A.l_co005$filtered_params_df[,"N_A_B"])
length(Results_object_TS_G_A.l_co005$filtered_params_df[,"N_A_B"])
a <- cut(Results_object_TS_G_A.l_co005$filtered_params_df[,"G_A"], breaks = 10)
a <- cut(Results_object_TS_G_A.l_co005$filtered_params_df[,"G_A"], breaks = c, include.lowest = TRUE)
b <- cut(c(1,100,Results_object_TS_G_A.l_co005$filtered_params_df[,"G_A"]), breaks = 10)

table(a)
table(b)

table(cut(c(1.1, 100), breaks = c))
cut.default(c(0,10), breaks = 10, labels = NULL, include.lowest = FALSE, right = FALSE)


as.data.frame(a$cutoff_tables[[1]])


which(as.data.frame(a$cutoff_tables[[1]])[,"Freq"] >= 0.8)

b <- as.character(as.data.frame(a$cutoff_tables[[1]])[1,1])
as.numeric(gsub("([0-9]+).*$", "\\1", b))

as.numeric(sapply(strsplit(as.character(as.data.frame(a$cutoff_tables[[1]])[2,1]), " "), "[[", 1))


regmatches(c(b), gregexpr("[[:digit:]]+", c(b)))

as.numeric(unlist(regmatches(b,
                             gregexpr("[[:digit:]]+\\.*[[:digit:]]*",b))))


b <- as.data.frame(a$cutoff_tables[[1]])
cutoff_rows <- which(b[,"Freq"] >= 0.8)

b <- b[cutoff_rows,]
c <- as.character(b[,"Var1"])
d <- as.numeric(unlist(regmatches(c,
                             gregexpr("[[:digit:]]+\\.*[[:digit:]]*",c))))
range(d)



#################################
dists_target <- as.matrix(pdist(experimental_object_TS.G_A.l$simulated_rset_expr_log_z, experimental_object_TS.G_A.l$simulated_rset_expr_log_z))
dists_target <- as.matrix(dist(experimental_object_TS.G_A.l$simulated_rset_expr_log_z))
