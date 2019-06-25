MutFreScore <- function(mutationMat){
  mut_fre = apply(mutationMat,1,sum)
  mut_fre_score = (mut_fre - min(mut_fre))/(max(mut_fre) - min(mut_fre))
  return(mut_fre_score)
}