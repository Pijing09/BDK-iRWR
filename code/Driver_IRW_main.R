Driver_IRWR_main <- function(network,mutation,tum_exp,seeds,root_path,damping){
  source(paste(root_path,"/code/net_mat2vec.R",sep = ""))
  source(paste(root_path,"/code/Katz_cent.R",sep = ""))
  source(paste(root_path,"/code/normalization.R",sep = ""))
  source(paste(root_path,"/code/BDK_iRWR_GeneScore.R",sep = ""))
  
  ###  calculate betweenness of the network
  library(igraph)
  network_vec <- net_mat2vec(network)    
  new_graphDiffCo_PPI <- make_graph(network_vec,directed = T)
  diff_co_PPI_bet <- betweenness(new_graphDiffCo_PPI,directed = T,normalized = T)
  diff_co_PPI_bet <- as.matrix(diff_co_PPI_bet)
  
  ###  caculate katz centrality of the network
  diff_co_PPI_katz <- Katz_cent(t(network))
  
  ###  obtain the interacting genes and normalization
  BKD_net_i_genes <- Reduce(intersect,list(v1 = rownames(diff_co_PPI_bet),v2 = rownames(diff_co_PPI_katz),
                                           v3 = rownames(network)))
  bet_norm_i <- as.matrix(normalization(diff_co_PPI_bet[BKD_net_i_genes,]))
  katz_norm_i <- as.matrix(normalization(diff_co_PPI_katz[BKD_net_i_genes,]))
  bet_katz <- as.matrix((bet_norm_i + katz_norm_i)/2)
  
  ###  obtain the bet_katz of the seeds as random jumping probility
  bet_katz_seeds <- bet_katz
  bet_katz_seeds[1:length(bet_katz_seeds),] <- 0
  seeds_i <- intersect(seeds,BKD_net_i_genes)
  bet_katz_seeds[seeds_i,] <- bet_katz[seeds_i,]
  
  ###  obtain the mean of tumor expression data as initialization
  tum_exp_i <- tum_exp[BKD_net_i_genes,]
  tum_exp_mean_i <- as.matrix(apply(tum_exp_i,1,mean))
  tum_exp_mean_norm_i <- normalization(tum_exp_mean_i)
  network_i <- network[BKD_net_i_genes,BKD_net_i_genes]
  rand_prob <- as.matrix(normalization(bet_katz_seeds))
  diff_co_PPI_seeds_S <- Driver_IRW_GeneScore(network_i,tum_exp_mean_norm_i,rand_prob,damping,epsilon = 1e-8,maxit=1000)
  
  ## selecting the mutated genes from 
  mut_sum <- rowSums(mutation)
  mut_genes <- rownames(mutation)[which(mut_sum > 0)]
  mut_IRW_genes_i <- intersect(mut_genes,rownames(diff_co_PPI_seeds_S))
  IRW_mut_score_n <- normalization(diff_co_PPI_seeds_S[mut_IRW_genes_i,])
  IRW_mut_score_rank <- IRW_mut_score_n[order(IRW_mut_score_n,decreasing = T)]
  return(IRW_mut_score_rank)
}

