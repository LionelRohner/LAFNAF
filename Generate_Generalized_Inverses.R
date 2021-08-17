

# ADD CAUTCHY SCHWARZ INEQUALITY


library(plm)

library(RConics)


A <- t(matrix(c(1,2,3,4,
                1,2,3,4,
                -1,0,1,-3,
                7,1,-2,-3), nrow = 4))




linDep_Cautchy_Schwartz <- function(A){
  # Consists of checking whether <u,v> >= ||u|| ||v||
  # Strict equality indicate linear dependence, i.e. <u,v> = ||u|| ||v||
  
  # output
  linDepIdx = data.frame(i = numeric(), j = numeric())
  
  # Cauchy-Schwarzt Loop
  for (i in 1:nrow(A)){
    for (j in i:nrow(A)){
      if (i != j){
        
        # get norms and dot prod
        u_norm = sqrt(sum(A[i,]^2))
        v_norm = sqrt(sum(A[j,]^2))
        u_v = A[i,] %*% A[j,]
        
        # check for strict equality
        if (u_v == u_norm * v_norm ){
          # mirrored duplicates should be excluded
          linDepIdx = rbind(linDepIdx,c(i,j)) 
        }
      }
    }
  }
  colnames(linDepIdx) = c("i","j")
  
  if (nrow(linDepIdx) == 0){
    message("No linear dependent cols/rows were found!")
    return(linDepIdx)
  } else {
    message("Linear combinations of cols/rows were found!")
    return(linDepIdx)
  }
}
  

for (i in 1:nrow(A)){
  for (j in i:nrow(A)){
  message(i,", ",j)}}

linDep_Cautchy_Schwartz(A)
  

A[-1,-1]

# It should only find one

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
