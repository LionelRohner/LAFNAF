

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


# symmetric + pos def
A <- t(matrix(c(2,2,1,
                2,5,1,
                1,1,2), nrow = 3))


isPosDef(A)

# full rank - positive 
A <- t(matrix(c(3,1,
                0,2), nrow = 2))

# zero
A <- t(matrix(c(0,0,
                0,0), nrow = 2))

# rotation - eigenvalues +- i
A <- t(matrix(c(0,-1,
                1,0), nrow = 2))












# Is Positive Definite ----------------------------------------------------

isPosDef <- function(A){# test 1 - symmetry
  
  # test 1 - is symmetric?
  test1 = all(round(A, 5) == round(t(A),5))
  
  message("Test 1 - Matrix is symmetric? ", test1)
  
  # test 2 - positive eigenvalues
  test2 = prod(eigen(A)$value) > 0
  
  message("Test 2 - All eigenvalues are positive? ", test2)
  
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




