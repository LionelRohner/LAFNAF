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
#
A <- t(matrix(c(1,-3,5,
                0,1,7,
                0,0,1), nrow = 3))


#
A <- t(matrix(c(6,3,2,
                12,6,-2,
                24,12,-1), nrow = 3))

# 
A <- t(matrix(c(6,12,24,
                3,6,12,
                0,0,0,
                2,-2,1,
                0,0,0,
                1,1,1,
                0,0,0), nrow = 3))

# Auxilliaries ------------------------------------------------------------

# check if column vector below pivot is all zero (used in rref)
col_Is_All_Zero = function(A,currentCol){
  return(all(A[-c(1:currentCol),currentCol] == 0))
}


# perform row operations (used in rref)
gaussian_Elimination <- function(A, idxNonzero, currentCol){
  for (ele in idxNonzero){
    A[ele,] = A[ele,]-A[ele,currentCol]*A[currentCol,]
  }
  return(A)
}

# RM if only used in swap_Zero_Vectors
find_Zero_Vectors = function(A){
  return(which(rowSums(sqrt(A^2)) == 0))
}

# used to put zero vectors at the bottom if found in matrix
swap_Zero_Vectors <- function(A){
  
  zeroVecIdx = find_Zero_Vectors(A)
  
  if (length(zeroVecIdx) == 0){
    return(A)
  } else {
    return(append_to_bottom(A,zeroVecIdx))
  }
}

remove_Parallel_Vectors <- function(A){
  # find parallel vectors
  rmIdx = unique(linDep_Cautchy_Schwartz_Matrix(A)$j)
  
  # put zero-vector and parallel vectors at the bottom
  if (!is.null(rmIdx)){
    A[rmIdx,] <- rep(0,nc)
    return(append_to_bottom(A,rmIdx))
  } else {
    return(A)
  }
}



# -------------------------------------------------------------------------

# 1.) Setup variables
nc = ncol(A)
nr = nrow(A)

# first pivot must be on first column
currentCol = 1
idxPivot = 1

# 1.1.) remove parallel vectors
rmIdx = unique(linDep_Cautchy_Schwartz_Matrix(A)$j)

# 1.2.) set small values to zero
ifelse(abs(A) <= 1e-10, 0, A)

# put zero-vector and parallel vectors at the bottom
if (!is.null(rmIdx)){
  A[rmIdx,] <- rep(0,nc)
  A <- append_to_bottom(A,rmIdx)
} 

# 2.) find first pivot (if exists, else normalize first row to first element)
firstPivotIdx = which(A[,1]==1)[1]

# check if a pivot exists?
if (!is.na(firstPivotIdx)){
  A <- swap(A,firstPivotIdx,1)
  
} else {
  # Create a pivot (normalize first row)
  A[1,] = (1/A[1,1]) * A[1,]
}

A

# # check whether all entries of first column == 0
# col_Is_All_Zero = function(A,idxPivot,currentCol){
#   return(!any(A[-c(idxPivot),currentCol] == 0))
# }


# 2.1.) Check if column are all 0 except for the row with the pivot
if (!col_Is_All_Zero(A,currentCol = currentCol)){
  
  # find nonzero elements in col. vector
  idxNonzero = which(A[-idxPivot,currentCol]!=0) + currentCol
  
  # Perform Elementary Row OPs
  A = gaussian_Elimination(A,idxNonzero,currentCol)
} 


A

# 2.2.) Reorganize matrix such that zero vectors are at the bottom, rm parallel vectors
A = swap_Zero_Vectors(A)

A = remove_Parallel_Vectors(A)



A2 = A



A = A2

# 3.) find potential next pivots and swap rows if exist
nextPivot = 2
currentCol = 2
currentRow = 1


# make loop easier
for (col in currentCol:nc){
  
  # skip if col is zero
  if (all(A[-c(1:currentRow),col] == 0)){
    next
  }

  idxPivot = which(A[-c(1:currentRow),col] == 1)[1]
  
  # if no leading ones, create one
  if (is.na(idxPivot)){

    # if leading one does not exist, create one
    A[currentRow+1,] = (1/A[currentRow+1,col]) * A[currentRow+1,]
    currentRow = currentRow + 1
    
    # 2.1.) Check if column are all 0 except for the row with the pivot
    if (!col_Is_All_Zero(A,currentCol = col)){

      # find nonzero elements in col. vector
      idxNonzero = which(A[-c(1:currentRow),col]!=0) + col

      # Perform Elementary Row OPs
      A = gaussian_Elimination(A,idxNonzero,col)
    }
    
    next  
  }
  
  # if a leading one exists, swap and continue
  if (all(A[idxPivot,1:col-1] == 0)){
    # swap
    A <- swap(A,idxPivot,currentRow)
    
    # 2.1.) Check if column are all 0 except for the row with the pivot
    if (!col_Is_All_Zero(A,currentCol = col)){

      # find nonzero elements in col. vector
      idxNonzero = which(A[-c(1:currentRow),col]!=0) + col

      # Perform Elementary Row OPs
      A = gaussian_Elimination(A,idxNonzero,col)
    }
    
    # prep next it
    currentRow = currentRow + 1
  }
}


# detect pivots
A_test = A[-find_Zero_Vectors(A),]
nrTest = nrow(A_test)

if (all(diag(A_test) == 1)){
  # return
  diag(nrTest)
}

for (row in 1:nrTest){
  for (col in 1:nc){
    if (A_test[row, col] != 1){
      if (col == nc & row != nc){
        message("Failed")
      }
    }
  }
}

# WHAT WAS THAT ABOUT? ----------------------------------------------------

# FIND RANK?

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
