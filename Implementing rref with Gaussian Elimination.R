# find sub matrix that is a basis

# lin dep (row space)
A <- t(matrix(c(2,2,2,
                4,4,4,
                1,-3,5,
                0,1,7,
                2,-4,28,
                1,-2,12,
                0,0,1,
                0,0,1e-10,
                -1,1,3), nrow = 3))

nc = ncol(A)
nr = nrow(A)

# remove parallel vectors
rmIdx = unique(linDep_Cautchy_Schwartz_Matrix(A)$j)

# set small values to zero
ifelse(A <= 1e-10, 0, A)

# put zero-vector and parallel vectors at the bottom
if (!is.null(rmIdx)){
  A[rmIdx,] <- rep(0,nc)
  A <- append_to_bottom(A,rmIdx)
} 



firstPivot = which(A[,1]==1)[1]

if (!is.null(firstPivot)){
  A <- swap(A,firstPivot,1)
}

nextPivot = 2
idx = c()
for (i in 1:nr){
  for (j in 1:nc){
    if(i != 1 && j == nextPivot && A[i,j] == 1 && A[i,j-1] == 0){
      nextPivot = nextPivot + 1
      idx = c(idx,i)
    }
  }
}

A[idx,]


find_linDep_submatrix <- function(A){}








A <- t(matrix(c(2,2,2,
                4,4,4,
                6,6,6,
                0,1,7), nrow = 3))

A <- t(matrix(c(2,2,2,2,
                4,4,4,4,
                6,6,6,6,
                8,8,8,8,
                0,1,7,-1), nrow = 4))


nc = ncol(A)
nr = nrow(A)

# lin dep vectors must span the smaller of the row or column space of A
maxRank = min(nc,nr)

if (nc > nr){
  A = t(A)
}


all_combos = combn(1:nr,maxRank)

cnt = 0
for (combo in 1:ncol(all_combos)){
  subMatrix = A[all_combos[,combo],]
  if(det(subMatrix) != 0){
    break
  } 
  cnt = cnt + 1
}



if (cnt == ncol(all_combos)){
  while(det(subMatrix) == 0){
    
    if (exists("subsubmatrix")){
      rmd = sort(sample(1:nrow(subsubmatrix), maxRank - 1))
      subMatrix = subsubmatrix[rmd,rmd]
    } else {
      rmd = sort(sample(1:nrow(subMatrix), maxRank - 1))
      subMatrix = subMatrix[rmd,rmd]  
    }
  
    nc = ncol(subMatrix)
    nr = nrow(subMatrix)
    
    maxRank = min(nc,nr)
    
    all_combos = combn(1:nr,maxRank)
    
    cnt = 0
    for (combo in 1:ncol(all_combos)){
      subsubMatrix = subMatrix[all_combos[,combo],]
      if(det(subsubMatrix) != 0){
        break
      } 
      cnt = cnt + 1
    }
  }
}

# how to retrace the original indices
