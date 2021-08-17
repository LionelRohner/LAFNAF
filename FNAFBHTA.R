#------------------------------------------------------------------------------#
# Functions That Nobody Asked For, But Here They Are ----------------------
#------------------------------------------------------------------------------#

# Matrix Power ------------------------------------------------------------

mpow <- function(A, n){
  A_to_the_i <- A
  for (i in 1:n-1){
    A_to_the_i <- A%*%A_to_the_i
  }; return(A_to_the_i)
}


# Cautchy-Schwartz-Inequality to Identify Linear Dependent Cols/Ro --------

linDep_Cautchy_Schwartz <- function(A){
  # Consists of checking whether <u,v> >= ||u|| ||v||
  # Strict equality indicate linear dependence, i.e. <u,v> = ||u|| ||v||
  
  # output
  linDepIdx = data.frame(i = numeric(), j = numeric())
  
  # Cauchy-Schwarzt Loop
  for (i in 1:nrow(A)){
    # for (j in i:nrow(A))
    for (j in 1:nrow(A)){
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

