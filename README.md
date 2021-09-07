# LAFNAF
**L**inear **A**lgebra **F**unctions That **N**obody **A**sked **F**or, But Here They Are.

Originally, I wrote these functions to verify some linear algebra exercises. As the number of functions grew, I created this repository. The functions are not efficient (it's plain R-code) nor have they been thoroughly tested. Besides, all these functions have been implemented more efficiently in other packages or even in base R. It is just a fun exercise for me to get a better understanding of linear algebra.

### Main Functions:

* **create_Basis:** Create linear independent vectors, i.e. a basis for any real vector-space.
*  **is_Pos_Def:** Checks whether a matrix **A** is positive definite.
*  **mPow:** Calculate powers of some matrix **A**.
*  **canonical_Form:** Compute the decomposition of a matrix **A** into **UDU^-1**, where **U** is the matrix of eigenvectors and the diagonal matrix **D**, the similar canonical form, containing the eigenvalues of **A**.
*  **fastExp:** Fast exponentiation using the similar canonical form.
*  **linDep_Cautchy_Schwartz:** Checks for the Cautchy-Schwartz-Inequality between two vectors or between rows/cols of matrix to check for linear dependence.
*  **adjugate:** Creates an adjugate matrix of some matrix **A**.
*  **generalized_Inverse:** Should create generalized inverses of any matrix **A**, however, the part that should find a nonsingular submatrix in **A** does not always work.
*  **check_Penrose_Cond:** Check for the four Penrose condition for matrix inverses, i.e. A = AMA, M = MAM, AM = (AM)^T, MA = (MA)^T
*  **inverse:** Creates an inverse of a nonsingular matrix **A**.
*  **orthogonalize:** Create an orthogonal matrix from some matrix **A**. Used for my SVD implementation (computations of AA^T and A^TA).
*  **rank_Matrix:** Compute rank of matrix.
*  **singular_Value_Decomposition:** Does not work yet, I need to solve the problem of the sign ambiguity of eigenvector calculations.
*  **find_Zero_Vectors:** Find zero-vectors in a matrix using vector norm.
*  **rref:** Reduce some matrix **A** to reduced row echelon form. Not thoroughly tested!!!
  
### Plot Functions:
*  **plot_EigenVec:** Plots eigenvalues of a 2x2 matrix.
*  **plot_Matrix_Transformation:** Plots vectors **x** and **y** as in **Ax=y**, where **A** is 2x2.

### Auxiliary Functions:
* **compare_Floats:** Compare equality of floats. Principle: Compare absolute difference of two floats to a tolerance value (default = 1e-06). 
*  **swap:** Swap rows in a some matrix **A**. Used in **rref**.
*  **add_To_Bottom:** Add rows to the bottom of a matrix **A**.
*  **col_Is_All_Zero:** Check (TRUE/FALSE) whether a column is all zeros. Used in **rref**.
*  **gaussian_Elimination:** Performs Gaussian elimination on non-zero rows of a matrix. Used in **rref**.
*  **swap_Zero_Vectors:** Combination of **find_Zero_Vectors** and **add_to_bottom** used to rearrange rows in **rref**.
*  **remove_Parallel_Vectors:** Remove parallel vectors (rows or cols) in a matrix based on the Cautchy Schwartz Equality. Used in **rref**.
