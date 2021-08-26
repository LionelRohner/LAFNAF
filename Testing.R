

# Some matrices

# lin dep
A <- t(matrix(c(1,2,3,4,
                1,2,3,4,
                1,2,3,4,
                1,2,3,4),nrow = 4))

# lin dep rect
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

# lin dep (row space)
A <- t(matrix(c(1,-3,5,
                0,1,7,
                2,-4,28,
                3,-13,0,
                4,-14,-3,
                5,-9,30,
                6,2,173), nrow = 3))

# full rank 5x5
A <- t(matrix(c(1,2,3, 4,-2,
                -1,0,1,4,-7,
                7,1,-2,4,1,
                1,0,10,2,0,
                213,1,8,6,4), nrow = 5))

# full rank 3x3
A <- t(matrix(c(1,2,3,
                -1,0,1,
                7,1,-2), nrow = 3))

# lin dep 3x3 (Row 1 + Row 2 == Row 3)
A <- t(matrix(c(2,2,2,
                -2,2,-2,
                0,4,0), nrow = 3))

# full rank - all positive entries
A <- t(matrix(c(2,0,2,
                3,4,5,
                17,13,0), nrow = 3))

# from p.126 (matrix book) --> works for UDU-1
A <- t(matrix(c(-1,-2,-2,
                1,2,1,
                -1,-1,-1), nrow = 3))

# full rank - positive 
A <- t(matrix(c(3,1,
                0,2), nrow = 2))

# zero
A <- t(matrix(c(0,0,
                0,0), nrow = 2))

# rotation - eigenvalues +- i
A <- t(matrix(c(0,-1,
                1,0), nrow = 2))



  
fastExp <- function(A){
  
}

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
  
  # message(all(round(U%*%D%*%U_inv) == A))
  
  return(res) 
}

canonical_form(A)

# SVD --> Add sign checker ------------------------------------------------


# create P
P <- orthogonalize(A)

# create Q
Q <- orthogonalize(t(A))

eigenVal <- eigen(P)$values

D <- matrix(0,nrow = length(eigenVal), ncol = length(eigenVal))
diag(D) <- sqrt(eigenVal)

P%*%D%*%t(Q)

singular_Value_Decomposition <- function(A){
  
}




