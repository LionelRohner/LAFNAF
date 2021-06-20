

library(plm)

library(RConics)


A <- t(matrix(c(1,2,3,4,
                1,2,3,4,
                -1,0,1,-3,
                7,1,-2,-3), nrow = 4))

A <- t(matrix(c(4,1,2,0,
              1,1,5,15,
              3,1,3,5), ncol = 3))


A <- t(matrix(c(3,1,3,5,
                1,1,1,1,
                1,1,1,1), ncol = 3))




# step 1 : Find a LIN submatrix of order rxr 


# cols = setdiff(1:ncol(A),detect.lindep(A,suppressPrint = T))
# rows = setdiff(1:nrow(A),detect.lindep(t(A), suppressPrint = T))

q <- qr(A)
A_LIN <- A[,q$pivot[seq(q$rank)]]



# THATS REALLY BAD

if (nrow(A_LIN) > ncol(A_LIN)){
  rows <- 1:ncol(A_LIN)
  cols <- 1:ncol(A_LIN)
} else {
  cols <- 1:nrow(A_LIN)
  rows <- 1:nrow(A_LIN)
}

# should be rxr
W <- A_LIN[rows,cols]

# step 2 : (W^-1)^T

# create adjoint matrix - adjoint() doesnt work for 2x2
if (dim(W)[1] == 2){
  adj_W = -1*W
  diag(adj_W) = rev(diag(W))
} else {
  adj_W = adjoint(W)
}

transp_inv_W <- t(adj_W/det(W))

# step 3 : replace elements of A with (W^-1)^T
A_mod <- 0*A
A_mod[rows,cols] = transp_inv_W


# step 4 : t(A)
G <- t(A_mod)

# Test
A%*%G%*%A
