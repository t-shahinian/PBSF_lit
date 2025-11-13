library("Rtsne")
library("umap")
library("FNN")
library("igraph")
library("bluster")
library("parallel")
library("pdist")


plot_RGB_tSNE=function(location, latent_dat,pointsize=2,textsize=15){
  
  # suppressMessages(require(Rtsne))
  
  info = as.data.frame(location)
  colnames(info) = c("sdimx","sdimy")
  
  PCvalues = latent_dat
  
  tsne <- Rtsne(PCvalues,dims=3,check_duplicates = FALSE)
  r = (tsne$Y[,1]-min(tsne$Y[,1]))/(max(tsne$Y[,1])-min(tsne$Y[,1]))
  g = (tsne$Y[,2]-min(tsne$Y[,2]))/(max(tsne$Y[,2])-min(tsne$Y[,2]))
  b = (tsne$Y[,3]-min(tsne$Y[,3]))/(max(tsne$Y[,3])-min(tsne$Y[,3]))
  x =  info$sdimx
  y =  info$sdimy
  dat = data.frame(x,y,r,g,b)
  p1=ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
    geom_point(size=pointsize) +
    scale_color_identity()+
    #ggtitle(paste0("RGB tSNE"))+
    theme_void()+
    theme(#plot.title = element_text(size = textsize),
          text = element_text(size = textsize),
          #axis.title = element_text(face="bold"),
          #axis.text.x=element_text(size = 22) ,
          legend.position = "bottom")
  
  return(list("RGB"=dat,"figure"=p1))
}

plot_RGB_UMAP=function(location, latent_dat,pointsize=2,textsize=15){
  
  # suppressMessages(require(umap))
  
  info = as.data.frame(location)
  colnames(info) = c("sdimx","sdimy")
  
  PCvalues = latent_dat
  
  umap <- umap(PCvalues,n_components = 3)
  r = (umap$layout[,1]-min(umap$layout[,1]))/(max(umap$layout[,1])-min(umap$layout[,1]))
  g = (umap$layout[,2]-min(umap$layout[,2]))/(max(umap$layout[,2])-min(umap$layout[,2]))
  b = (umap$layout[,3]-min(umap$layout[,3]))/(max(umap$layout[,3])-min(umap$layout[,3]))
  x =  info$sdimx
  y = info$sdimy
  dat = data.frame(x,y,r,g,b)
  p1=ggplot(data=dat, aes(x=x, y=y, col=rgb(r,g,b))) +
    geom_point(size=pointsize) +
    scale_color_identity()+
    #ggtitle(paste0("RGB UMAP"))+
    theme_void()+
    theme(#plot.title = element_text(size = textsize),
          text = element_text(size = textsize),
          #axis.title = element_text(face="bold"),
          #axis.text.x=element_text(size = 22) ,
          legend.position = "bottom")
  
  return(list("RGB"=dat,"figure"=p1))
}


#' @title Obtain clustering cluster labels through louvain method.
#' @description This function performs louvain clustering on input low dimensional components.
#' @param clusternum The desired number of clusters the user wants to obtain.
#' @param latent_dat A d by n matrix of low dimensional components, d is number of PCs, n is number of spots.
#' @param knearest An integers, number of nearest neighbors for KNN graph construction in louvain clustering.
#' @return The cluster labels.
#'
#' @export
louvain_clustering = function(clusternum, latent_dat, knearest=100){
  set.seed(1234)
  # suppressMessages(require(FNN))
  # suppressMessages(require(igraph))
  # suppressMessages(require(bluster))
  
  PCvalues = latent_dat
  info.spatial = as.data.frame(t(PCvalues))
  colnames(info.spatial) =  paste0("factor", 1:nrow(PCvalues))
  knn.norm = FNN::get.knn(as.matrix(t(PCvalues)), k = knearest)
  knn.norm = data.frame(from = rep(1:nrow(knn.norm$nn.index),
                                   k=knearest), 
                        to = as.vector(knn.norm$nn.index), 
                        weight = 1/(1 + as.vector(knn.norm$nn.dist)))
  nw.norm = igraph::graph_from_data_frame(knn.norm, directed = FALSE)
  nw.norm = igraph::simplify(nw.norm)
  lc.norm = igraph::cluster_louvain(nw.norm)
  t <- proc.time()
  merged <- bluster::mergeCommunities(nw.norm, lc.norm$membership, 
                                      number=clusternum)
  proc.time() - t
  clusterlabel = as.character(as.integer(as.factor(paste0("cluster",merged))))
  return("cluster_label"=clusterlabel)
}


#' @title Visualize cluster labels on locations.
#' @description This function visualizes cluster labels on locations.
#' @param location A n by k matrix of spot locations.
#' @param clusterlabel A vector of cluster labels for spots.
#' @param pointsize An integer, the point size of each spot.
#' @param textsize An integer, the text size in the legend.
#' @param title_in A character string, the title you want to display at the top of the figure.
#' @param color_in A vector of colors for each cluster.
#' @param legend A character string, the position of the figure legend. Select from "top", "bottom","left" or "right".
#' @return A ggplot object.
#' @export
plot_cluster = function(location, clusterlabel, pointsize=3, text_size=15, 
                        title_in, color_in, legend="none"){
  cluster = clusterlabel
  loc_x=location[,1]
  loc_y=location[,2]
  datt = data.frame(cluster, loc_x, loc_y)
  p = ggplot(datt, aes(x = location[,1], y = location[,2], color = cluster)) +
    geom_point( alpha = 1,size=pointsize) +
    scale_color_manual(values = color_in)+
    ggtitle(paste0(title_in))+
    theme_void()+
    theme(plot.title = element_text(size = text_size,  face = "bold"),
          text = element_text(size = text_size),
          #axis.title = element_text(face="bold"),
          #axis.text.x=element_text(size = 15) ,
          legend.position =legend)
  p
}

get_PCA = function(expr,PCnum){
  
  n=dim(expr)[2]
  k = dim(expr)[1]
  output_sub_mean=matrix(0,k,n)
  for(i_k in 1:k){
    output_sub_mean[i_k,]=expr[i_k,]-mean(expr[i_k,])
  }
  svd_output_sub_mean=svd(output_sub_mean)
  A_ini=svd_output_sub_mean$u[,1:PCnum]
  Z_pca = t(A_ini) %*% output_sub_mean
  return(Z_pca)
}


#' @import pdist
fx_1NN = function(i,location_in){
  # library(pdist)
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  return(min(line_i))
}

#' @title Calculate CHAOS score to measure clustering performance.
#' @description CHAOS score measures the spatial continuity of the detected spatial domains.
#' Lower CHAOS score indicates better spatial domian clustering performance.
#' @param clusterlabel Cluster labels.
#' @param location A n by k matrix of spatial locations.
#' @return A numeric value for CHAOS score.
#'
#' @import parallel
#'
#' @export
fx_CHAOS = function(clusterlabel, location){
  # require(parallel)
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  matched_location = scale(matched_location)
  dist_val = rep(0,length(unique(clusterlabel)))
  count = 0
  for(k in unique(clusterlabel)){
    count = count + 1
    location_cluster = matched_location[which(clusterlabel == k),]
    if(length(location_cluster)==2){next}
    #require(parallel)
    results = mclapply(1:dim(location_cluster)[1], fx_1NN, 
                       location_in=location_cluster,mc.cores = 5)
    dist_val[count] = sum(unlist(results))
  }
  dist_val = na.omit(dist_val)
  return(sum(dist_val)/length(clusterlabel))
  
}


fx_PAS = function(clusterlabel, location){
  # require(parallel)
  
  matched_location=location
  NAs = which(is.na(clusterlabel))
  if(length(NAs>0)){
    clusterlabel=clusterlabel[-NAs]
    matched_location = matched_location[-NAs,]
  }
  
  results = mclapply(1:dim(matched_location)[1], fx_kNN, 
                     location_in=matched_location,k=10,cluster_in=clusterlabel, 
                     mc.cores = 5)
  return(sum(unlist(results))/length(clusterlabel))
}


fx_kNN = function(i,location_in,k,cluster_in){
  #library(pdist)
  line_i = rep(0,dim(location_in)[1])
  line_i = pdist(location_in[i,],location_in[-i,])@dist
  ind = order(line_i)[1:k]
  cluster_use = cluster_in[-i]
  if(sum(cluster_use[ind] != cluster_in[i])>(k/2)){
    return(1)
  }else{
    return(0)
  }
}
