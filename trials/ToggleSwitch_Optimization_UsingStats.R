##############################Functions########################################################################################################################################################
##############################generate_standard_object########################################################################################################################################################
generate_standard_object <- function(topology, standard_state, numModels = 10000, integrateStepSize = 0.02) {
  
  reference_rset <- sracipeSimulate(topology, integrate = FALSE,
                                    numModels = numModels, genParams = TRUE,
                                    integrateStepSize = integrateStepSize)
  
  state <- standard_state
  
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
  reference_rset_object <- list(reference_rset = ref_expression_rset, standard_state = standard_state, expression_ref = expr_reference_log_z, prcomp_ref = prcomp_reference,
                                means_ref = means_reference, sds_ref = sds_reference, topology = topology)
  
  reference_rset_object
  
}


##############################generate_experimental_object########################################################################################################################################################
generate_experimental_object <- function(experimental_state = experimental_state, standard_rset_object, numModels = 10000, integrateStepSize = 0.02, save_rset = FALSE) {
  
  state <- experimental_state
  
  parameter_names <- colnames(sracipeParams(standard_rset_object$reference_rset))
  
  hill_params <- grep(pattern = "N", parameter_names)
  
  simulated_rset <- sracipeSimulate(standard_rset_object$topology, integrate = FALSE,
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
  
  simulated_rset <- sracipeSimulate(simulated_rset, integrate = TRUE, genParams = FALSE)
  
  simulated_rset_expr_log <- log2(t(assay(simulated_rset)))
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log,
                                     2, standard_rset_object$means_ref, FUN = "-")
  
  simulated_rset_expr_log_z <- sweep(simulated_rset_expr_log_z,
                                     2, standard_rset_object$sds_ref, FUN = "/")
  
  #rotate
  simulated_rset_expr_PCs <- simulated_rset_expr_log_z %*% standard_rset_object$prcomp_ref$rotation
  
  
  expression_object <- list(simulated_rset_expr_log_z = simulated_rset_expr_log_z, simulated_rset_expr_PCs = simulated_rset_expr_PCs)
  
  if (save_rset == TRUE) {
    
    expression_object <- list(simulated_rset_expr_log_z = simulated_rset_expr_log_z, simulated_rset_expr_PCs = simulated_rset_expr_PCs, simulated_rset = simulated_rset,
                              arguments = list(numModels = numModels, integrateStepSize = integrateStepSize, state = state))
    
  }
  expression_object
}


##############################Generate_distances_list########################################################################################################################################################
Generate_distances_list <- function(standard_object, experimental_expression, dist.method = "Euclidean") {
  
  #Intialize distance_df
  ncol_df <- 2
  nrow_df <- nrow(standard_object$expression_ref)
  
  distance_df <- data.frame(matrix(ncol = ncol_df , nrow = nrow(experimental_expression)))
  distance_df[,ncol_df] <- 1:nrow_df
  
  #Initialize list of distances
  distance_df_list <- vector(mode = "list", length = nrow_df)
  
  #Experimental ncol and nrow
  nrow_experimental_expr <- nrow(experimental_expression)
  ncol_experimental_expr <- ncol(experimental_expression)
  
  for (i in 1:nrow_df) {
    
    if (i%%100 == 0) {
      print(i)
    }
    
    for (j in 1:nrow_experimental_expr) {
      
      if (dist.method == "Euclidean") {
        
        distance_df[j, 1] <- dist(rbind(standard_object$expression_ref[i,], experimental_expression[j,]))[1]
        
      } else if (dist.method == "Manhattan") {
        
        vec <- vector(length = ncol_experimental_expr)
        
        for (k in 1:ncol_experimental_expr) {
          
          vec[k] <- abs(standard_object$expression_ref[i, k] - experimental_expression[j, k])
          
        }
        
        distance_df[j, 1] <- sum(vec)
      }
    }
    distance_df_list[[i]] <- distance_df
    
  }
  output <- distance_df_list
  output
}
  

##############################Rank_Models########################################################################################################################################################


##############################Testing########################################################################################################################################################

topology_TS <- data.frame(Source = c("A", "B"), Target = c("B", "A"), Type = c(2,2))

standard_state_TS <- list(c(1,100), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100)) #standard

standard_object_TS <- generate_standard_object(topology = topology_TS, standard_state = standard_state_TS, numModels = 20000)

#nearest neighbor distance
standard_expr_TS_knn.G_A.l <- spatstat::nndist(X = standard_object_TS$expression_ref, k = 1)
plot(density(standard_expr_TS_knn.G_A.l))
as.data.frame(standard_object_TS$expression_ref[which(standard_expr_TS_knn.G_A.l <= 0.008),])
ggplot() + geom_point(data = as.data.frame(standard_object_TS$expression_ref[which(standard_expr_TS_knn.G_A.l <= 0.008),]), aes(A,B))


# library(ggplot2)
####Plots
ggplot(data = as.data.frame(standard_object_TS$expression_ref),aes(A, B)) + geom_point() + labs(title="Standard Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

ggplot() + geom_density(data = as.data.frame(standard_object_TS$expression_ref), aes(A, color = "A")) + geom_density(data = as.data.frame(standard_object_TS$expression_ref), aes(B, color = "B")) +
  scale_colour_manual(values=c("A"="blue", "B"="red"), name="Genes") + labs(title="Standard Expression Density") + theme(plot.title = element_text(hjust = 0.5))

###kmeans
set.seed(123)

standard_expression_TS.kmeans <- kmeans(standard_object_TS$expression_ref, centers = 2, nstart = 25)
standard_expression_TS.kmeans$centers
standard_expression_TS.kmeans.df <- as.data.frame(cbind(standard_object_TS$expression_ref, standard_expression_TS.kmeans$cluster ))
colnames(standard_expression_TS.kmeans.df)[3] <- "Cluster"
standard_expression_TS.kmeans.df[, "Cluster"] <- as.factor(standard_expression_TS.kmeans.df[, "Cluster"])

###Plots_fancy
standard_expression_TS.scatterPlot <- ggplot(standard_expression_TS.kmeans.df,aes(A, B, color = Cluster)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.0006,0.095), legend.justification=c(0,1)) + labs(title="Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + 
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B), color='black',size=8, shape=4 )

# Marginal density plot of x (top panel)
Adensity <- ggplot(standard_expression_TS.kmeans.df, aes(A, fill=Cluster)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none") + labs(title="Expression Density of A") + theme(plot.title = element_text(hjust = 0.5))

# Marginal density plot of y (right panel)
Bdensity <- ggplot(standard_expression_TS.kmeans.df, aes(B, fill=Cluster)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none") + labs(title="Expression Density of B") + theme(plot.title = element_text(hjust = 0.5))


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

# library("gridExtra")
grid.arrange(Adensity, blankPlot, standard_expression_TS.scatterPlot, Bdensity, 
             ncol=2, nrow=2)

# plot(density(standard_object_TS$expression_ref[, "A"]), col = "red", main = "Expression Density of Standard")
# # lines(density(experimental_object_TS.G_A.l$simulated_rset_expr_log_z[, "B"]), col = "green")
# par(new = TRUE)
# plot(density(standard_object_TS$expression_ref[, "B"]), col = "green", axes = FALSE, xlab = "", ylab = "", main = "")
# 
# legend(-1.9,0.4, legend=c("A", "B"),
#        col=c("red", "green"), lty=1:1, cex=0.8)


# library(spatstat)
##############################G_A########################################################################################################################################################

experimental_state_TS.G_A.l <- list(c(1,10), c(1,100), c(0.1,1), c(0.1, 1), c(1,55), c(1,55), c(1,5), c(1,5), c(1,100), c(1,100)) #experimental

experimental_object_TS.G_A.l <- generate_experimental_object(experimental_state = experimental_state_TS.G_A.l, standard_rset_object = standard_object_TS, numModels = 20000, integrateStepSize = 0.02, save_rset = TRUE)



#Plots
ggplot(data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z),aes(A, B)) + geom_point() + 
  labs(title="Experimental Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

ggplot() + geom_density(data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z), aes(A, color = "A")) + 
  geom_density(data = as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z), aes(B, color = "B")) +
  scale_colour_manual(values=c("A"="blue", "B"="red"), name="Genes") + labs(title="Experimental Expression Density") + theme(plot.title = element_text(hjust = 0.5))

# library(fpc)
db <- fpc::dbscan(experimental_object_TS.G_A.l$simulated_rset_expr_log_z, eps = 0.15, MinPts = 5)
# Plot DBSCAN results
plot(db, experimental_object_TS.G_A.l, main = "DBSCAN", frame = FALSE)
###kmeans
set.seed(123)

experimental_expression_TS.kmeans <- kmeans(experimental_object_TS.G_A.l$simulated_rset_expr_log_z, centers = 2, nstart = 25)
experimental_expression_TS.kmeans.df <- as.data.frame(cbind(experimental_object_TS.G_A.l$simulated_rset_expr_log_z, experimental_expression_TS.kmeans$cluster ))
colnames(experimental_expression_TS.kmeans.df)[3] <- "Cluster"
experimental_expression_TS.kmeans.df[, "Cluster"] <- as.factor(experimental_expression_TS.kmeans.df[, "Cluster"])


###Plots_fancy
ggplot(experimental_expression_TS.kmeans.df,aes(A, B)) +geom_point() +  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B), color='red',size=8, shape=4 ) +
  labs(title="Experimental Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5))


#overlap
  ggplot() + 
  geom_point(data = standard_expression_TS.kmeans.df, aes(A, B, color = "Standard")) +
  geom_point( data = experimental_expression_TS.kmeans.df, aes(A, B, color = "Experimental"), shape = 1, alpha = 0.8) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B, color="Centers"), shape=4, size=8) + 
  labs(title="Expression Scatterplot") +  theme(plot.title = element_text(hjust = 0.5), legend.position = c(0, 0), legend.justification = c("left", "bottom")) +
  scale_colour_manual(values=c("Standard"="black", "Experimental"="blue", "Centers"="red"), name="Data") + 
  guides(color = guide_legend(override.aes = list(shape = c(16, 1, 4), size = c(1.5, 1.5, 8)), title.hjust = .5))


experimental_expression_TS.scatterPlot <- ggplot(experimental_expression_TS.kmeans.df,aes(A, B, color = Cluster)) + 
  geom_point() + 
  scale_color_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position=c(0.0006,0.095), legend.justification=c(0,1)) + labs(title="Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5))

# Marginal density plot of x (top panel)
Adensity <- ggplot(experimental_expression_TS.kmeans.df, aes(A)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none") + labs(title="Expression Density of A") + theme(plot.title = element_text(hjust = 0.5))

# Marginal density plot of y (right panel)
Bdensity <- ggplot(experimental_expression_TS.kmeans.df, aes(B, fill=Cluster)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = c('#999999','#E69F00')) + 
  theme(legend.position = "none") + labs(title="Expression Density of B") + theme(plot.title = element_text(hjust = 0.5))

ggplot(experimental_expression_TS.kmeans.df, aes(A)) + 
  geom_density(alpha=.5,  color = "green") + geom_density(data = experimental_expression_TS.kmeans.df, aes(B), alpha=.5, color = "red") + theme(legend.position="top")


blankPlot <- ggplot()+geom_blank(aes(1,1))+
  theme(plot.background = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_blank(),
        axis.ticks = element_blank()
  )

# library("gridExtra")
grid.arrange(Adensity, blankPlot, standard_expression_TS.scatterPlot, Bdensity, 
             ncol=2, nrow=2)


###Using Cutoff
# experimental_expr_TS_knn.G_A.l <- spatstat::nndist(X = experimental_object_TS.G_A.l$simulated_rset_expr_log_z, k = 1)

plot(density(experimental_expr_TS_knn.G_A.l), main = "Distance to Nearest Neighbor"); abline(v=0.01, col = "red")

experimental_expr_TS.G_A.l_filtered <- as.data.frame(experimental_object_TS.G_A.l$simulated_rset_expr_log_z[which(experimental_expr_TS_knn.G_A.l <= 0.01),])

distances_list_TS.G_A.l <- Generate_distances_list(standard_object = standard_object_TS, experimental_expression = experimental_object_TS.G_A.l$simulated_rset_expr_log_z)

#Plots
ggplot(data = experimental_expr_TS.G_A.l_filtered, aes(A, B)) + geom_point() + 
  labs(title="Experimental Expression Scatterplot with nn Distance Cutoff of 0.01") + theme(plot.title = element_text(hjust = 0.5)) + xlim(-3.5, 2) + ylim(-3.5, 2)

ggplot() + geom_density(data = experimental_expr_TS.G_A.l_filtered, aes(A, color = "A")) + 
  geom_density(data = experimental_expr_TS.G_A.l_filtered, aes(B, color = "B")) +
  scale_colour_manual(values=c("A"="blue", "B"="red"), name="Genes") + labs(title="Experimental Expression Density with nn Distance Cutoff of 0.01") + theme(plot.title = element_text(hjust = 0.5))


plot(density(experimental_expr_TS.G_A.l_filtered[, "A"]), col = "red", main = "Expression Density of Experimental with knn cutoff of 0.01")
par(new = TRUE)
plot(density(experimental_expr_TS.G_A.l_filtered[, "B"]), col = "green", axes = FALSE, xlab = "", ylab = "", main = "")
legend(-2.1,0.5, legend=c("A", "B"),
       col=c("red", "green"), lty=1:1, cex=0.8)

ggplot(experimental_expr_TS.G_A.l_filtered,aes(A, B)) +geom_point() +  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B), color='red',size=8, shape=4 ) +
  labs(title="Experimental Expression Scatterplot with nn cutoff of 0.01") + theme(plot.title = element_text(hjust = 0.5)) 


ggplot(experimental_expr_TS.G_A.l_filtered, aes(A)) + 
  geom_density(alpha=.5,  color = "green") + geom_density(data = experimental_expr_TS.G_A.l_filtered, aes(B), alpha=.5, color = "red")


ggplot() + 
  geom_point(data = standard_expression_TS.kmeans.df, aes(A, B), color = "black") +
  geom_point( data = experimental_expr_TS.G_A.l_filtered, aes(A, B), color = "blue", shape = 1, alpha = 0.8) +
  geom_point(data = as.data.frame(standard_expression_TS.kmeans$centers), aes(A, B), color='red',size=8, shape=4 ) + 
  labs(title="Expression Scatterplot") + theme(plot.title = element_text(hjust = 0.5))