

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

# full rank - positive
A <- t(matrix(c(3,1,
                0,2), nrow = 2))

# zero
A <- t(matrix(c(0,0,
                0,0), nrow = 2))

# rotation - eigenvalues +- i
A <- t(matrix(c(0,-1,
                1,0), nrow = 2))


plotMatrixTransformation <- function(A,v,
                                     offset = 1,
                                     plotBasisVecs = T){
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
    maxMat = max(v_trans)
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
  
  
  # legend
  legend("bottomright", c("Basis","Vector y (Ax=y)"),
         col = c(rgb(0,0,0.8,0.7),rgb(0.8,0,0.8,0.7)), 
         pch = c(16,16), inset=c(0,1), xpd=TRUE, horiz=T, bty="n",
         cex = 0.75)
  
  # restore par settings to default
  par(mfrow = c(1,1))
}
  

v = c(1,1)
plotMatrixTransformation(A,v, offset = 1)



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




