# LAFNAF
<ins>L</ins>inear <ins>A</ins>lgebra <ins>F</ins>unctions That <ins>N</ins>obody <ins>A</ins>sked <ins>F</ins>or, But Here They Are.

Originally, I wrote these functions to verify some linear algebra exercises. As the number of functions grew, I created this repository. The functions are not efficient (it's plain R-code) nor have they been thoroughly tested. Besides, all these functions have been implemented more efficiently in other packages or even in base R. It is just a fun exercise for me to get a better understanding of linear algebra.

### To Do List:
* The **singular_Value_Decomposition** function does not work correctly as the eigenvectors produced by R have a random sign, thus the matrix product **PDQ^-1** (corresponds to **UΣV^-1**) does not reconstitute **A**. Solution: Implement sign-flip function from https://digital.library.unt.edu/ark:/67531/metadc900575/m2/1/high_res_d/920802.pdf.
* ~~The **rref** function does not return the reduced row echelon form, since it does not set free variables, which are located above a pivot to zero. Solution: Repeat Gaussian Elimination for free variables.~~
* ~~The **generalized_Inverse** function, I am still trying to find a way to find a nonsingular matrix in **A** without using the builtin function **qr**. Once rref works fine, this can be used to find rank r and consequetnly to find a submatrix by checking, which rxr submatrix has a determinant > 0.~~
* ~~Once **rref** works, rewrite **rank_Matrix**.~~

### Main Functions:

* **create_Basis:** Create a set of linear independent vectors that span a real vector-space.
*  **is_Pos_Def:** Checks whether a matrix **A** is positive definite.
*  **mPow:** Calculate powers of some matrix **A**.
*  **canonical_Form:** Compute the decomposition of a matrix **A** into **UDU^-1**, where **U** is the matrix of eigenvectors and the diagonal matrix **D**, the similar canonical form, containing the eigenvalues of **A**.
*  **fastExp:** Fast exponentiation using the similar canonical form.
*  **linDep_Cautchy_Schwartz:** Checks for the Cautchy-Schwartz-Inequality between two vectors or between rows/cols of matrix to check for linear dependence.
*  **adjugate:** Creates an adjugate matrix of some matrix **A**.
*  **generalized_Inverse:** Create generalized inverses of any matrix **A**. Needs more testing.
*  **check_Penrose_Cond:** Check for the four Penrose condition for matrix inverses, i.e. A = AMA, M = MAM, AM = (AM)^T, MA = (MA)^T
*  **inverse:** Creates an inverse of a nonsingular matrix **A**.
*  **orthogonalize:** Create an orthogonal matrix from some matrix **A**. Used for my SVD implementation (computations of AA^T and A^TA).
*  **rank_Matrix:** Compute rank of matrix.
*  **singular_Value_Decomposition:** Does not work yet, I need to solve the problem of the sign ambiguity of eigenvector calculations.
*  **find_Zero_Vectors:** Find zero-vectors in a matrix using vector norm.
*  **ref:** Reduce some matrix **A** to row echelon form. Used for **rank_Matrix**, because its quicker than **rref**.
*  **rref:** Reduce some matrix **A** to reduced row echelon form. Validated to some extent... Needs more testing.

  
### Plot Functions:
*  **plot_EigenVec:** Plots eigenvalues of a 2x2 matrix.
*  **plot_Matrix_Transformation:** Plots vectors **x** and **y** as in **Ax=y**, where **A** is 2x2.

### Auxiliary Functions:
* **compare_Floats:** Compare equality of floats. Principle: Compare absolute difference of two floats to a tolerance value (default = 1e-06). 
*  **swap:** Swap rows in a some matrix **A**. Used in **rref**.
*  **add_To_Bottom:** Add rows to the bottom of a matrix **A**. Used in **rref**.
*  **col_Is_All_Zero:** Check (TRUE/FALSE) whether a column is all zeros. Used in **rref**.
*  **gaussian_Elimination:** Performs Gaussian elimination on non-zero rows of a matrix. Used in **rref**.
*  **swap_Zero_Vectors:** Combination of **find_Zero_Vectors** and **add_to_bottom** used to rearrange rows in **rref**.
*  **remove_Parallel_Vectors:** Remove parallel vectors (rows or cols) in a matrix based on the Cautchy Schwartz Equality. Used in **rref**.
