##############################Generate_results_object###########################################################################################################
# install.packages("pdist")
# library(pdist)

Generate_results_object <- function(standard_object, experimental_expression, cutoff = 0.075) {
  
  nparams <- ncol(sracipeParams(standard_object$reference_rset))
  param_names <- colnames(sracipeParams(standard_object$reference_rset))
  
  standard_expression.kmeans <- kmeans(standard_object$expression_ref, centers = 2, nstart = 25)
  
  #Intialize distance_df
  ncol_df <- 2
  nrow_df <- nrow(standard_object$expression_ref)
  
  distance_df <- data.frame(matrix(ncol = ncol_df , nrow = nrow_df))
  distance_df[,ncol_df] <- 1:nrow_df
  
  #Calc distances
  dists <- as.matrix(pdist(standard_object$expression_ref, experimental_expression))
  
  #Generate list of distances
  for (i in 1:nrow(distance_df)) {
    
    distance_df[i, 1] <- dists[i,][order(dists[i,])[1]]
    
  }
  
  colnames(distance_df) <- c("Distance_nearest_exp_neighbor", "ModelID")
  
  
  ###For 2 genes
  
  #Distances Distribution
  Distance_density_plot <- ggplot() +geom_density(data = distance_df, aes(Distance_nearest_exp_neighbor)) + geom_vline(xintercept = cutoff, color = "red") +
    labs(title="Distance to Nearest Experimental Neighbor") + theme(plot.title = element_text(hjust = 0.5)) + xlab("Euclidean Distance")
  
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


###############################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))

standard_state_TS <- list(c(1,100), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100)) #standard

standard_object_TS_100k <- generate_standard_object(topology = topology_TS, standard_state = standard_state_TS, numModels = 100000)

standard_object_TS_250k <- generate_standard_object(topology = topology_TS, standard_state = standard_state_TS, numModels = 250000)

standard_expr_A <- standard_object_TS_100k$expression_ref[, "A"]

standard_expression_TS.kmeans <- kmeans(standard_object_TS_100k$expression_ref, centers = 2, nstart = 25)


experimental_state_TS.G_A.l <- list(c(1,10), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100)) #experimental

experimental_object_TS.G_A.l <- generate_experimental_object(experimental_state = experimental_state_TS.G_A.l, standard_rset_object = standard_object_TS_100k, numModels = 1000,
                                                             integrateStepSize = 0.02, save_rset = TRUE)
experimental_expr_A <- experimental_object_TS.G_A.l$simulated_rset_expr_log_z[, "A"]


ks_test <- dgof::ks.test(x = standard_expr_A , y = experimental_expr_A)
ks_test$data.name

group <- c(rep("Standard", length(standard_expr_A)), rep("Experimental", length(experimental_expr_A)))
dat <- data.frame(KSD = c(standard_expr_A,experimental_expr_A), group = group)

cdf1 <- ecdf(standard_expr_A) 
cdf2 <- ecdf(experimental_expr_A) 
# find min and max statistics to draw line between points of greatest distance
minMax <- seq(min(standard_expr_A, experimental_expr_A),
              max(standard_expr_A, experimental_expr_A), length.out=length(standard_expr_A)) 
x0 <- minMax[which( abs(cdf1(minMax) - cdf2(minMax)) == max(abs(cdf1(minMax) - cdf2(minMax))) )] 
y0 <- cdf1(x0) 
y1 <- cdf2(x0) 

ggplot(dat, aes(x = KSD, group = group, color = group))+
  stat_ecdf(size=1) +
  theme_bw(base_size = 28) +
  theme(legend.position ="top") +
  xlab("Expression") +
  ylab("ECDF") +
  #geom_line(size=1) +
  geom_segment(aes(x = x0[1], y = y0[1], xend = x0[1], yend = y1[1]),
               linetype = "dashed", color = "red") +
  geom_point(aes(x = x0[1] , y= y0[1]), color="red", size=8) +
  geom_point(aes(x = x0[1] , y= y1[1]), color="red", size=8) +
  ggtitle("K-S Test: Standard vs Experimental") +
  theme(legend.title=element_blank()) + theme(plot.title = element_text(hjust = 0.5))


ggplot(dat, aes(x = KSD, group = group, color = group)) + geom_density() + theme_bw(base_size = 14) +
  theme(legend.position ="top") +
  xlab("Expression") +
  ylab("Density") +  ggtitle("Standard vs Experimental Expression of A") + theme(plot.title = element_text(hjust = 0.5))

#Standard expression 
ggplot(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(A, B)) + geom_point() + 
  labs(title="Standard Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

ggplot() + geom_line(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(A, color = "A"), stat = "density") + 
  geom_line(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(B, color = "B"), stat = "density") +
  scale_colour_manual(values=c("A"="blue", "B"="red"), name="Genes") + labs(title="Standard Expression") + 
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Expression")

#param_filtered_standard_expression 
param_filtered_standard_expression <- as.data.frame(standard_object_TS_100k$expression_ref[which(sracipeParams(standard_object_TS_100k$reference_rset)[,"G_A"] <= 10),])

ggplot(data = param_filtered_standard_expression, aes(A, B)) + geom_point() + 
  labs(title="Param-filtered Standard Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

ggplot() + geom_line(data = param_filtered_standard_expression, aes(A, color = "A"), stat = "density") + 
  geom_line(data = param_filtered_standard_expression, aes(B, color = "B"), stat = "density") +
  scale_colour_manual(values=c("A"="blue", "B"="red"), name="Genes") + labs(title="Param-filtered Standard Expression") + 
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Expression")

ggplot() + 
  geom_point(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(A, B, color = "DF1")) +
  geom_point( data = param_filtered_standard_expression, aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
  labs(title="Expression Overlay") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "black", DF2 = "blue", DF3 = "red"), 
                      labels = c(DF1 = "Standard", DF2 = "Param-Filtered Standard", DF3 = "Centers"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))

#Experimental expression 
ggplot(data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B)) + geom_point() + 
  labs(title="Experimental Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

ggplot() + geom_line(data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, color = "A"), stat = "density") + 
  geom_line(data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z), aes(B, color = "B"), stat = "density") +
  scale_colour_manual(values=c("A"="blue", "B"="red"), name="Genes") + labs(title="Experimental Expression") + 
  theme(plot.title = element_text(hjust = 0.5)) + xlab("Expression")

ggplot() + 
  geom_point(data = as.data.frame(standard_object_TS_100k$expression_ref), aes(A, B, color = "DF1")) +
  geom_point( data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, B, color = "DF2"), shape = 1, alpha = 0.8) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="DF3"), shape=4, size=8) + 
  labs(title="Expression Overlay") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values = c(DF1 = "black", DF2 = "blue", DF3 = "red"), 
                      labels = c(DF1 = "Standard", DF2 = "Experimental Expression", DF3 = "Centers"),
                      name = "Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))

#Distance-filtered expression
Results_object_TS_G_A.l_co065 <- Generate_results_object(standard_object = standard_object_TS_100k, 
                                                         experimental_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z, cutoff = 0.065)

ggplot(data = Results_object_TS_G_A.l_co065$filtered_standard_expr, aes(A, B)) + geom_point() + 
  labs(title="Distance-Filtered Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

Results_object_TS_G_A.l_co065$Expression_overlay_plot

Results_object_TS_G_A.l_co065

do.call("grid.arrange", Results_object_TS_G_A.l_co065$param_plots_list)

ggplot() + geom_density(data = as.data.frame(Results_object_TS_G_A.l_co065$shortest_distance_df[which(Results_object_TS_G_A.l_co065$shortest_distance_df[,1] <= 0.2),]),
                        aes(Distance_nearest_exp_neighbor))


Results_object_TS_G_A.l_co005 <- Generate_results_object(standard_object = standard_object_TS_100k, 
                                                         experimental_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z, cutoff = 0.005)

# ggplot(data = Results_object_TS_G_A.l_co055$filtered_standard_expr, aes(A, B)) + geom_point() + 
#   labs(title="Distance-Filtered Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

Results_object_TS_G_A.l_co055$filtered_params_df

Results_object_TS_G_A.l_co01$Distance_density_plot
do.call("grid.arrange", Results_object_TS_G_A.l_co005$param_plots_list)

ggplot(data = Results_object_TS_G_A.l_co005$filtered_params_df, aes(G_A)) + geom_histogram() + stat_bin(bins = 4) + 
  scale_x_continuous(breaks = c(0, 100))

ggplot(data = Results_object_TS_G_A.l_co005$filtered_params_df, aes(N_B_A)) +
geom_bar() +
  scale_x_binned(n.breaks = 10)


Results_object_TS_G_A.l_co005$Distance_density_plot


ggplot(Results_object_TS_G_A.l_co01$filtered_params_df) +
  geom_histogram(aes(G_A), breaks = c(0, 100), color = "black") + 
  scale_x_continuous(breaks = c(0, 100))

which(Results_object_TS_G_A.l_co01$filtered_params_df[,"G_A"] < 1)

Results_object_TS_G_A.l_co01$filtered_params_df
###
test2 <- binarize.kMeans(vect = Results_object_TS_G_A.l_co065$filtered_params_df[,"G_A"])
plot(cbind(test2@originalMeasurements, test2@binarizedMeasurements))

test <- binarize.BASC(vect = Results_object_TS_G_A.l_co065$filtered_params_df[,"G_A"])
plotStepFunctions(test)
