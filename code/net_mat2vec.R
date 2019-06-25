#### transform the 0-1 matrix to the specific format in igraph package
net_mat2vec  <-  function(net_matrix){
  edge = list()
  for (i in 1:nrow(net_matrix)){
    edge[[i]] = names(which(net_matrix[i,] == 1))
  }
  names(edge) = rownames(net_matrix)
  
  edge_mat = matrix(0,1,2)
  for(i in 1:nrow(net_matrix)){
	if (length(edge[[i]]) != 0){    
	   edge_mat_i = cbind(names(edge)[i],edge[[i]])
	   edge_mat = rbind(edge_mat,edge_mat_i)
	}       
  }
  edge_mat = edge_mat[-1,]
  edge_v = list(t(edge_mat))
  edge_v = unlist(edge_v)
}


