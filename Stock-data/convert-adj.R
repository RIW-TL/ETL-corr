matrix_to_edgelist_vectorized <- function(adj_matrix, node_names) {
  # 获取上三角矩阵的索引
  upper_tri_indices <- which(upper.tri(adj_matrix) & adj_matrix == 1, arr.ind = TRUE)
  
  # 创建边列表
  edges <- data.frame(
    Source = upper_tri_indices[, 1],
    Target = upper_tri_indices[, 2],
    Type = "Undirected"
  )
  
  # 如果有节点名称，使用名称而不是索引
  if (!is.null(node_names)) {
    edges$Source <- node_names[edges$Source]
    edges$Target <- node_names[edges$Target]
  }
  
  return(edges)
}

node_names = read.table("clipboard",header = FALSE)[,1]

Ome_ECM_hat[which(Ome_ECM_hat != 0)] = 1
to_ECM = matrix_to_edgelist_vectorized(Ome_ECM_hat, node_names)
write.csv(to_ECM,file = "ECM.csv")

Ome_SECM_hat[which(Ome_SECM_hat != 0)] = 1
to_SECM = matrix_to_edgelist_vectorized(Ome_SECM_hat, node_names)
write.csv(to_SECM,file = "SECM.csv")

Ome_single_hat[which(Ome_single_hat != 0)] = 1
to_single = matrix_to_edgelist_vectorized(Ome_single_hat, node_names)
write.csv(to_single,file = "single.csv")

Ome_pool_hat[which(Ome_pool_hat != 0)] = 1
to_pool = matrix_to_edgelist_vectorized(Ome_pool_hat, node_names)
write.csv(to_pool,file = "Pool.csv")


submat = Ome_ECM_hat[1:35,36:67]
submat[abs(submat) < 0.01] = 0
Ome_ECM_hat[1:35,36:67] = submat

submat = Ome_ECM_hat[1:35,68:104]
submat[abs(submat) < 0.005] = 0
Ome_ECM_hat[1:35,68:104] = submat

submat = Ome_ECM_hat[36:67,68:104]
submat[abs(submat) < 0.005] = 0
Ome_ECM_hat[36:67,68:104] = submat

Ome_ECM_hat[36:67,1:35] = t(Ome_ECM_hat[1:35,36:67])
Ome_ECM_hat[68:104,1:35] = t(Ome_ECM_hat[1:35,68:104])
Ome_ECM_hat[68:104,36:67] = t(Ome_ECM_hat[36:67,68:104])

#-----------------------------------------------------------------


edge = function(Ome_hat){
  
  p = ncol(Ome_hat)
  
  density.in = c()
  density.in[1] = length(which(Ome_hat[1:35,1:35] != 0))/(35^2)
  density.in[2] = length(which(Ome_hat[36:67,36:67] != 0))/(32^2)
  density.in[3] = length(which(Ome_hat[68:104,68:104] != 0))/(37^2)

  
  density.cross = c()
  density.cross[1] = length(which(Ome_hat[1:35,36:67] != 0))/(35*32)
  density.cross[2] = length(which(Ome_hat[1:35,68:104] != 0))/(35*37)
  density.cross[3] = length(which(Ome_hat[36:67,68:104] != 0))/(32*37)
  
  return(list(density.in = density.in, density.cross = density.cross))
}



edge.detect.result = matrix(0,6,6)
rownames(edge.detect.result) = 
  c("Target-Glasso","Trans-ECM-Glasso","Trans-SECM-Glasso",
                                 "Pooled-Glasso","CLIME","Trans-CLIME")
colnames(edge.detect.result) = 
  c("class-1","class-2","class-3","class1-2","class1-3","class2-3")

a = edge(Ome_single_hat)
edge.detect.result[1,] = c(a$density.in, a$density.cross)

a = edge(Ome_ECM_hat)
edge.detect.result[2,] = c(a$density.in, a$density.cross)

a = edge(Ome_SECM_hat)
edge.detect.result[3,] = c(a$density.in, a$density.cross)

a = edge(Ome_pool_hat)
edge.detect.result[4,] = c(a$density.in, a$density.cross)

a = edge(Theta.initial)
edge.detect.result[5,] = c(a$density.in, a$density.cross)

a = edge(Omega.trans.clime.hat)
edge.detect.result[6,] = c(a$density.in, a$density.cross)
