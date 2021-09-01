#------------------------------------------------------------------------------#
# Linear Algebra Functions Nobody Asked For -------------------------------
#------------------------------------------------------------------------------#

# None of the functions below have been thourougly tested!


# create linear independent vectors aka a basis ---------------------------

create_basis <- function(dim, lower = 0, upper = 9){
  basis = matrix(sample(lower:upper,dim^2), nrow = dim)
  
  while (det(basis) == 0){
    basis = matrix(sample(lower:upper,dim^2), nrow = dim)
  }
  
  return(basis)
  
  # output list of vecs instead of matrix
  out = list()
  for (rowcol in 1:dim){
    out[[rowcol]] = basis[,rowcol]
  }
  return(out)
}

# Is Positive Definite ----------------------------------------------------

isPosDef <- function(A){# test 1 - symmetry
  
  # test 1 - is symmetric?
  test1 = all(round(A, 5) == round(t(A),5))
  
  message("Test 1 - Matrix is symmetric? ", test1)
  
  # test 2 - positive eigenvalues
  test2 = prod(Re(eigen(A)$value)) > 0
  
  message("Test 2 - All (real) eigenvalues are positive? ", test2)
  
  # test 3 - positive upper left submatrices
  dimA = nrow(A)
  upLeftDets = c(A[1,1])
  if (A[1,1] < 0){
    test3 = FALSE
  } else { 
    for (i in 2:dimA){
      upLeftDet = det(A[1:i,1:i])
      if (upLeftDet < 0){
        test3 = FALSE
      } else {
        upLeftDets = c(upLeftDets,det(A[1:i,1:i]))
      }
    }
  }
  test3 = all(upLeftDets > 0)
  
  message("Test 3 - Determinant of upper left submatrices are > 0? ", test3)
  
  # output  
  return(all(test1,test2,test3))
}

# Matrix Power ------------------------------------------------------------

mpow <- function(A, n){
  # uglllyyyy
  n = n-1
  
  if(n == 0){
    return(A)
  } else if (n == -1){
    # return identity matrix
    I <- matrix(0,nrow = nrow(A), ncol = ncol(A))
    diag(I) = rep(1,nrow(A))
    return(I)
  }
  
  A_to_the_i <- A
  for (i in 1:n){
    A_to_the_i <- A%*%A_to_the_i
  }; return(A_to_the_i)
}


# Canonical Form (A = UDU^-1 ----------------------------------------------

canonical_form <- function(A){
  # create UDU
  
  # create U
  eign = eigen(A)
  U = eign$vectors
  
  # create D
  D <- matrix(0,nrow = length(eign$values), ncol = length(eign$values))
  diag(D) = eign$values
  
  # create U-1
  U_inv = adjugate(U)*1/det(U)
  
  # output
  res <- list(U = U, D = D, U_inv = U_inv)
  return(res) 
  
  # message(all(round(U%*%D%*%U_inv) == A))
}



# Fast Exponentiation Using Canonical Form --------------------------------

fastExp <- function(A,p){
  UDU = canonical_form(A)
  
  # fast exponentiation of eigenvalues
  D <- UDU$D
  D_to_the_p = diag(D)^p
  
  # calculate UD^pU^-1
  # U doesnt have to be included in the power calc, because A = UDU^-1, thus A^2 = UDU^-1UDU^-1, which is UDDU^-1
  diag(D) = D_to_the_p
  
  return(UDU$U%*%D%*%UDU$U_inv)
}


# Single Cautchy-Schwartz Ineq Test ---------------------------------------

linDep_Cautchy_Schwartz <- function(u,v, tol = 1e-05){
  # Consists of checking whether <u,v> >= ||u|| ||v||
  # Strict equality indicate linear dependence, i.e. <u,v> = ||u|| ||v||
  
  
  # By transposing the matrix here, the algo needs to know whether the idx
  # output is meant for the row or col space!!!
  
  # # Check whether there are more rows or cols, then choose shorter space
  # # This reduces the # of iterations in the double loop
  # if(ncol(A) > nrow(A)){
  #   A = t(A)
  # } 
  
  # compute norms of the vectors
  u_norm = sqrt(sum(u^2))
  v_norm = sqrt(sum(v^2))
  u_v_norm = u_norm * v_norm
  
  # dot product of vectors
  u_v = drop(v %*% u)

  

  # check for strict equality to find lin. dependent rows/cols
  if (abs(u_v-u_v_norm) <= tol){
    return(FALSE)
  } else {
    return(TRUE)
  }
  
}



# Cautchy-Schwartz-Inequality to Identify Linear Dependent Cols/Rows ------


linDep_Cautchy_Schwartz_Matrix <- function(A){
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
        
        # check for strict equality to find lin. dependent rows/cols
        if (u_v == u_norm * v_norm ){
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
 


# Create an Adjugate Matrix (transpose of a cofactor matrix) --------------

adjugate <- function(A){
  
  n = nrow(A)
  m = ncol(A)
  
  if (n != m){
    message("Matrix must be square!")
    return(NULL)
  }
  
  # create emtpy cofactor matrix
  C = matrix(NA,nrow = n, ncol = m)
  
  # special case for A of order 2x2
  if (n == 2){
    C = (-1)*A
    diag(C) = rev(diag(A))
    return(t(C))
  }
  
  # populate the cofactor matrix
  for (i in 1:n){
    for (j in 1:m){
      C[i,j] = (-1)^(i+j)*det(A[-i,-j])
    }
  }
  
  return(t(C))
}

# Create generalized inverse from a mxn matrix ----------------------------

# does only work for obvious cases of lin dep that are actually identified
# by the Cautchy-Schwartz Inequality (e.g. [2,4] = 2*[1,2]), but not for
# compositve linear dependencies (e.g. [4,0] = [2,-2] + [2,2]).


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
  
  # create adjoint matrix - adjoint() doesnt work for 2x2
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


# Check Penrose Condition of Inverses -------------------------------------

check_Penrose_Cond <- function(A,
                               G,
                               all_Penrose_check = F,
                               digits = 2){
  
  message("Disclaimer: All matrices are transformed to pure integer matrices first, so consider that...\n")
  
  if (all_Penrose_check == F){
    Penrose_1 = all(round(A%*%G%*%A,2) == A)
    message("AGA=A is ", Penrose_1)
    return(Penrose_1)
    }else{
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


# Inverse of Matrix (Square) ----------------------------------------------

inverse <- function(A){
  # Check for squareness
  if (nrow(A) != ncol(A)){
    message("Matrix is not invertible")
    return(NULL)
  }
  
  # check for singularity
  if (det(A) == 0){
    message("Matrix is not invertible")
    return(NULL)
  }
  
  # invert!
  adjA = adjugate(A)
  return(adjA*det(A)^-1)
}


#  Orthogonalize Matrix (AAT or ATA) --------------------------------------

orthogonalize <- function(A){
  preQ = A%*%t(A)
  Q = eigen(preQ, symmetric = T)$vectors
  for (row in 1:nrow(Q)){
    # Normalize Rows
    Q[row,] = 1/sqrt(c(Q[row,]%*%Q[row,]))*c(Q[row,])
  }
  return(Q)
}

# Rank --------------------------------------------------------------------

rank_Matrix <- function(A, tol = 1e-12){
  return(sum(abs(Re(eigen(A)$values)) > tol))
}


# SVD ---------------------------------------------------------------------


singular_Value_Decomposition <- function(A){
  # create P (left signular vectors)
  P <- orthogonalize(A)
  
  # create Q (right-singular vectors)
  Q <- orthogonalize(t(A))
  
  # get eigenvalues
  eigenVal <- eigen(A%*%t(A))$values
  
  # create diagonal matrix of signular values
  D <- matrix(0,nrow = length(eigenVal), ncol = length(eigenVal))
  diag(D) <- sqrt(eigenVal)
  
  # output generation
  res <- list(P = P,D = D, Q = Q)
  
  message("Warning, there may be sign ambiguity in singular vectors!!!")
  
  return(res)
}



# Plotting Eigenvectors from 2x2 Matrix -----------------------------------

plotEigenVec <- function(A,
                         offset = 1,
                         plotBasisVecs = T,
                         plotSpan = T,
                         plotTransBasis = T){
  
  # assumptions for function:
  if (all(dim(A) != 2)){
    message("A is not 2x2. Exiting...")
    return(NULL)
  }
  
  par(pty="s")
  
  maxMat = max(A) + offset
  
  plot(1, type = "n",                        
       xlab = "", ylab = "",
       xlim = c(-maxMat, maxMat), ylim = c(-maxMat, maxMat),
       main = paste("Eigenvalues = ",eigen(A)$values), cex.main = 0.75)
  grid()
  
  # basis vectors
  if (plotBasisVecs){
    arrows(0,0,0,1, length = 0.05)
    arrows(0,0,1,0, length = 0.05)
  }
  
  # span basis vectors
  arrows(0,-1*2*maxMat,0,1*2*maxMat,
         length = 0.05, col = rgb(0,0,0,0.3), code = 0, lty = 2)
  arrows(-1*2*maxMat,0,1*2*maxMat,0,
         length = 0.05, col = rgb(0,0,0,0.3), code = 0, lty = 2)
  
  # transformed basis vectors
  if (plotTransBasis){
    arrows(0,0,A[1,1],A[2,1], length = 0.05, col = rgb(0,0,0.8,0.7))
    arrows(0,0,A[1,2],A[2,2], length = 0.05, col = rgb(0,0,0.8,0.7))
  }
  
  # span transformed basis vectors
  if (plotSpan){
    arrows(-A[1,1]*maxMat,-A[2,1]*maxMat,
           A[1,1]*maxMat,A[2,1]*maxMat,
           length = 0.05, col = rgb(0,0,0.8,0.3), code = 0, lty = 2)
    arrows(-A[1,2]*maxMat,-A[2,2]*maxMat,
           A[1,2]*maxMat,A[2,2]*maxMat,
           length = 0.05, col = rgb(0,0,0.8,0.3), code = 0, lty = 2)
  }
  
  
  # eigenvectors
  eign = eigen(A)
  
  if (any(Im(eign$values) != 0)){
    message("Some eigenvalues are complex, but only real values are considered!\n")
  }
  
  # imagenary part of eigenvectors is stripped
  eigenVec_1 = Re(eign$vectors[,1])*Re(eign$values[1]) 
  eigenVec_2 = Re(eign$vectors[,2])*Re(eign$values[2])
  
  # eigenvectors plotting
  arrows(0,0,eigenVec_1[1],eigenVec_1[2], length = 0.05, col = rgb(0.8,0,0.8,0.7))
  arrows(0,0,eigenVec_2[1],eigenVec_2[2], length = 0.05, col = rgb(0.8,0,0,0.7))
  
  # span eigenvectors plotting
  if (plotSpan){
    arrows(-eigenVec_1[1]*maxMat,-eigenVec_1[2]*maxMat,
           eigenVec_1[1]*maxMat,eigenVec_1[2]*maxMat,
           length = 0.05, col = rgb(0.8,0,0.8,0.3), lty = 2) # aes
    arrows(-eigenVec_2[1]*maxMat,-eigenVec_2[2]*maxMat,
           eigenVec_2[1]*maxMat,eigenVec_2[2]*maxMat,
           length = 0.05, col = rgb(0.8,0,0,0.3), lty = 2) # aes
  }
  
  # plot origin
  points(0,0, pch = 16, cex = 0.7,)
  
  # legend
  legend("topleft", c("Basis", "Transformed Basis","Scaled Eigenvectors"),
         col = c(rgb(0,0,0,1),rgb(0,0,0.8,0.7),rgb(0.8,0,0.8,0.7)), 
         pch = c(16,16,16),
         inset=c(1,0), xpd=TRUE, horiz=F, bty="n",
         cex = 0.75)
  
  # restore par settings to default
  par(mfrow = c(1,1))
}


# Plot Matrix Transformation Ax = y ---------------------------------------

plotMatrixTransformation <- function(A,v,
                                     offset = 1,
                                     plotBasisVecs = T,
                                     splitPlot = T){
  # assumptions for function:
  if (all(dim(A) != 2)){
    message("A is not 2x2. Exiting...")
    return(NULL)
  }
  
  par(mfrow = c(1,2),pty="s")
  
  maxMat = max(A) + offset
  
  plot(1, type = "n",                        
       xlab = "", ylab = "",
       xlim = c(-maxMat, maxMat), ylim = c(-maxMat, maxMat),
       main = "Before Transformation", cex.main = 0.75)
  grid()
  
  # basis vectors
  if (plotBasisVecs){
    arrows(0,0,0,1, length = 0.05)
    arrows(0,0,1,0, length = 0.05)
  }
  
  # span basis vectors
  arrows(0,-1*2*maxMat,0,1*2*maxMat,
         length = 0.05, col = rgb(0,0,0,0.3), code = 0, lty = 2)
  arrows(-1*2*maxMat,0,1*2*maxMat,0,
         length = 0.05, col = rgb(0,0,0,0.3), code = 0, lty = 2)
  
  # input vector
  arrows(0,0,v[1],v[2], length = 0.05, col = rgb(0.8,0,0,0.8))
  text(v[1],v[2]+round(offset/2),paste("[",v[1],v[2],"]"),
       col = rgb(0.8,0,0,0.8), cex = 0.75)
  
  # plot origin and tip of vec
  points(0,0, pch = 16, cex = 0.7,)
  points(0,0, pch = 16, cex = 0.7,)
  
  
  # legend
  legend("bottomright", c("Basis","Vector x (Ax=y)"),
         col = c(rgb(0,0,0.8,0.7),rgb(0.8,0,0.8,0.7)), 
         pch = c(16,16), inset=c(0,1), xpd=TRUE, horiz=T, bty="n",
         cex = 0.75)
  
  # second plot with transformation
  
  # find new boundaries 
  v_trans = A%*%v
  if (max(v_trans) > maxMat){
    maxMat = max(v_trans) + offset
  }
  
  if (splitPlot){
    # emtpy plot
    plot(1, type = "n",                        
         xlab = "", ylab = "",
         xlim = c(-maxMat, maxMat), ylim = c(-maxMat, maxMat),
         main = "After Transformation", cex.main = 0.75)
    
    # formula for getting angle between vectors
    angle <- function(v1,v2){
      acos( sum(v1*v2) / ( sqrt(sum(v1*v1)) * sqrt(sum(v2*v2)) ) )
    }
    
    # find angle
    angle1 = angle(c(1,0),A[,1])
    angle2 = -angle(c(0,1),A[,2])
    
    
    # apply grid (swipe through grid in intervals (seq) and draw abline with slope = angle)
    sapply(seq(-maxMat*10,maxMat*10,by=1), function(inter) abline(a=inter, 
                                                                  b=tan(angle1),
                                                                  lty=3,
                                                                  col="lightgray"))
    sapply(seq(-maxMat*10,maxMat*10,by=4), function(inter) abline(a=inter, 
                                                                  b=tan(angle2+pi/2),
                                                                  lty=3,
                                                                  col="lightgray"))
  }
  
  
  # transformed basis vectors
  if (plotBasisVecs){
    arrows(0,0,A[1,1],A[2,1], length = 0.05, col = rgb(0,0,0.8,0.7))
    arrows(0,0,A[1,2],A[2,2], length = 0.05, col = rgb(0,0,0.8,0.7))
  }
  
  # span transformed basis vectors
  
  arrows(-A[1,1]*maxMat,-A[2,1]*maxMat,
         A[1,1]*maxMat,A[2,1]*maxMat,
         length = 0.05, col = rgb(0,0,0.8,0.3), code = 0, lty = 2)
  arrows(-A[1,2]*maxMat,-A[2,2]*maxMat,
         A[1,2]*maxMat,A[2,2]*maxMat,
         length = 0.05, col = rgb(0,0,0.8,0.3), code = 0, lty = 2)
  
  
  
  # input vector transformed
  arrows(0,0,v_trans[1],v_trans[2], length = 0.05, col = rgb(0.5,0,0.8,0.8))
  text(v_trans[1],v_trans[2]+round(offset/2),paste("[",v_trans[1],v_trans[2],"]"),
       col = rgb(0.5,0,0.8,0.8), cex = 0.75)
  
  # legend
  legend("bottomright", c("Basis","Vector y (Ax=y)"),
         col = c(rgb(0,0,0.8,0.7),rgb(0.8,0,0.8,0.7)), 
         pch = c(16,16), inset=c(0,1), xpd=TRUE, horiz=T, bty="n",
         cex = 0.75)
  
  # restore par settings to default
  par(mfrow = c(1,1))
}

plotMatrixTransformation(A,c(1,3),offset = 2)

