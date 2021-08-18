

# ADD CAUTCHY SCHWARZ INEQUALITY



A <- t(matrix(c(1,2,3,4,
                1,2,3,4,
                1,2,3,4,
                1,2,3,4,
                -1,0,1,-3,
                7,1,-2,-3,
                1,2,3,4,
                1,2,3,4,
                1,2,3,4,
                1,2,3,4), nrow = 4))


A <- t(matrix(c(1,-3,5,
                0,1,7,
                2,-4,28,
                3,-13,0,
                4,-14,-3,
                5,-9,30,
                6,2,173), nrow = 3))

A <- t(matrix(c(1,2,3, 4,-2,
                -1,0,1,4,-7,
                7,1,-2,4,1,
                1,0,10,2,0,
                213,1,8,6,4), nrow = 5))

A <- t(matrix(c(1,2,3,
                -1,0,1,
                7,1,-2), nrow = 3))

singular_Value_Decomposition <- function(A){
  # A = P [D 0] Q'
  P
}

rank_Matrix <- function(A){
  
  # if A is mxn with m != n, then the larger (be it rows or cols) will be lin dep!
  if (ncol(A) > nrow(A)){
    A = t(A)
  }
  
  # lin.dep indices
  CS_indices <- unique(linDep_Cautchy_Schwartz(A)$j)
  print(CS_indices)
  # if CS_indices is empty, A is full-rank, hence rank = nrow(A) = ncol(A)
  if (is.null(CS_indices)){
    message("Matrix is full-rank!")
    return(nrow(A))
  } else {
    # not full rank situation, due to transposition, nrow(A) > ncol(A)
    rank = nrow(A)-length(CS_indices)
    message("Matrix has order ", nrow(A),"x",ncol(A), " and rank ", rank)
    return(rank)
  }
}


rank_Matrix(A)

linDep_Cautchy_Schwartz <- function(A){
  # Consists of checking whether <u,v> >= ||u|| ||v||
  # Strict equality indicate linear dependence, i.e. <u,v> = ||u|| ||v||
  
  
  # By transposing the matrix here, the algo needs to know whether the idx
  # output is meant for the row or col space!!!
  
  # # Check whether there are more rows or cols, then choose shorter space
  # # This reduces the # of iterations in the double loop
  # if(ncol(A) > nrow(A)){
  #   A = t(A)
  # } 
  
  # output container
  linDepIdx = c()
  
  # Cauchy-Schwarzt Loop
  for (i in 1:nrow(A)){
    for (j in i:nrow(A)){
      if (i != j){
        # compute norms and dot prod of the span of the row/col space
        u_norm = sqrt(sum(A[i,]^2))
        v_norm = sqrt(sum(A[j,]^2))
        u_v = A[i,] %*% A[j,]
        
        
        print("it")
        print(round(u_norm) * round(v_norm))
        print(u_v)
        
        
        # check for strict equality to find lin. dependent rows/cols
        # if (isTRUE(all.equal(u_v == u_norm * v_norm, round(u_v) == round(u_norm) * round(v_norm)))){
        if (u_v == u_norm * v_norm){
          # mirrored duplicates should be excluded
          linDepIdx = rbind.data.frame(linDepIdx,c(i,j)) 
        }
      }
    }
  }
  
  # output construction
  if (is.null(nrow(linDepIdx))){
    message("No linear dependent cols/rows were found!")
    return(NULL)
  } else {
    colnames(linDepIdx) = c("i","j")
    message("Linear combinations of cols/rows were found!")
    return(linDepIdx)
  }
}

linDep_Cautchy_Schwartz(A)

adjugate <- function(A){
  
  n = nrow(A)
  m = ncol(A)
  
  if (n != m){
    message("Matrix must be square!")
    return(NULL)
  }
  
  # create emtpy cofactor matrix
  C = matrix(NA,nrow = n, ncol = m)
  
  # populate the cofactor matrix
  for (i in 1:n){
    for (j in 1:m){
      C[i,j] = (-1)^(i+j)*det(A[-i,-j])
    }
  }
  
  return(t(C))
}

generalized_Inverse <- function(A){
  
  ### step 1 : Find a LIN submatrix of order rxr 
  message("Step 1 : Find a LIN submatrix W of order rxr in A!")
  
  # All idxs in the second col should be lin. dep. with some vecs in the first col
  # Hence, removing those entries should guarantee NON-singularity of A
  rmIdx <- unique(linDep_Cautchy_Schwartz(A)$j)
  
  # Check whether A is full rank anyways
  if (is.null(rmIdx)){
    W = A
  } else {
    # make matrix square
    diff_row_col = abs((nrow(A)-length(rmIdx))-ncol(A))
    W <- A[-rmIdx, -diff_row_col]
  }
  
  # check assumptions for non-singularity
  if (det(W) == 0){
    message("Hope that never happens ;)")
    return(NULL)
  }
  
  ### step 2 : (W^-1)^T
  message("Step 2 : (W^-1)^T!")
  
  # create adjugate matrix - adjugate() doesnt work for 2x2
  if (dim(W)[1] == 2){
    adj_W = -1*W
    diag(adj_W) = rev(diag(W))
  } else {
    adj_W = adjugate(W)
  }
  
  transp_inv_W <- t(adj_W/det(W))
  
  ### step 3 : replace elements of A with (W^-1)^T
  message("Step 3 : Project rows/cols of (W^-1)^T into a a zero matrix A0 of order dim(A)!")
  
  # Again, if A is non-singular from scratch skip Step 3
  if (is.null(rmIdx)){
    A0 = transp_inv_W
  } else{
    A0 <- 0*A
    A0[-rmIdx, -diff_row_col] = transp_inv_W
    
  }
  
  
  ### step 4 : t(A)
  message("Step 3 : G = A0^T!")
  
  G <- t(A0)
  
  return(G)
}

check_Penrose_Cond <- function(A,
                               G,
                               all_Penrose_check = F,
                               digits = 2){
  
  message("Disclaimer: All matrices are transformed to pure integer matrices first, so consider that...\n")
  
  if (all_Penrose_check == F){
    Penrose_1 = all(round(A%*%G%*%A,2) == A)
    message("AGA=A is ", Penrose_1)
    return(Penrose_1)
  } else{
    Penrose_1 = all(round(A%*%G%*%A) == A)
    message("AGA = A is ", Penrose_1)
    
    Penrose_2 = all(round(G%*%A%*%G,digits = digits) == round(G,digits = digits))
    message("GAG = G is ", Penrose_2)
    
    Penrose_3 = all(round(A%*%G) == t(round(A%*%G)))
    message("AG = (AG)' is ", Penrose_3)
    
    Penrose_4 = all(round(G%*%A) == t(round(G%*%A)))
    message("GA = (GA)' is ", Penrose_4)
    
    return(all(Penrose_1,Penrose_2,Penrose_3,Penrose_4))
  }
}


G <- generalized_Inverse(A)

check_Penrose_Cond(A,G, all_Penrose_check = T)

round(G%*%A%*%G,2) == round(G,2)

