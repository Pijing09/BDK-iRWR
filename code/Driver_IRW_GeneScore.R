Driver_IRW_GeneScore <- function(NewNetwork,prior_score,rand_p,damping,epsilon = 0.00000001,maxit=500){
  AdjMat = t(NewNetwork)
  Adj_colsum=apply(AdjMat,2,sum)
  trans_Mat = matrix(0,nrow(AdjMat),ncol(AdjMat))
  
  degreeMat = matrix(rep(Adj_colsum,length(Adj_colsum)),nrow = length(Adj_colsum))   
  trans_degree = degreeMat * AdjMat
  trans_degree_colsum = apply(trans_degree,2,sum)
  
  for (i in 1:ncol(AdjMat)){ 
    if (trans_degree_colsum[i]!=0){
      trans_Mat[,i] = 0.85*(trans_degree[,i]/trans_degree_colsum[i]) + 0.15*(1/nrow(AdjMat))
    }
    else {
      trans_Mat[,i]=matrix(1/nrow(AdjMat),nrow(AdjMat),1)
    }  
  }
  
  node_num_col = ncol(trans_Mat)
  node_score_norm = prior_score/(sum(prior_score))
  rand_p_norm = rand_p/sum(rand_p)
  node_score_t = matrix(0,node_num_col,1)
  
  for(k in 1:maxit){
    
    node_score_t = damping*trans_Mat%*%node_score_norm + (1-damping)* rand_p_norm
    node_score_t = node_score_t/sum(node_score_t)
    wucha = sqrt(sum((node_score_t - node_score_norm)^2))
    node_score_norm = node_score_t
    message(paste("the",k,"th iteration"))
    if(wucha <= epsilon) break
	
  }
  node_score = (node_score_t - min(node_score_t))/(max(node_score_t) - min(node_score_t))
  return(node_score)
}