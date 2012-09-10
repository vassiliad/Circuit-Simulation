////////////////////////////////////////////////////////////////////////////////
// File: unit_lower_triangular.c                                              //
// Routines:                                                                  //
//    Unit_Lower_Triangular_Solve                                             //
//    Unit_Lower_Triangular_Inverse                                           //
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
//  void Unit_Lower_Triangular_Solve(double *L, double *B, double x[], int n) //
//                                                                            //
//  Description:                                                              //
//     This routine solves the linear equation Lx = B, where L is an n x n    //
//     unit lower triangular matrix.  (Only the subdiagonal part of the matrix//
//     is addressed.)  The diagonal is assumed to consist of 1's and is not   //
//     addressed.                                                             //
//     The algorithm follows:                                                 //
//                          x[0] = B[0], and                                  //
//            x[i] = B[i] - (L[i][0] * x[0]  + ... + L[i][i-1] * x[i-1]),     //
//     for i = 1, ..., n-1.                                                   //
//                                                                            //
//  Arguments:                                                                //
//     double *L   Pointer to the first element of the unit lower triangular  //
//                 matrix.                                                    //
//     double *B   Pointer to the column vector, (n x 1) matrix, B.           //
//     double *x   Pointer to the column vector, (n x 1) matrix, x.           //
//     int     n   The number of rows or columns of the matrix L.             //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double A[N][N], B[N], x[N];                                            //
//                                                                            //
//     (your code to create matrix A and column vector B)                     //
//     Unit_Lower_Triangular_Solve(&A[0][0], B, x, n);                        //
//     printf(" The solution is \n");                                         //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Unit_Lower_Triangular_Solve(double *L, double B[], double x[], int n)
{
   int i, k;

//         Solve the linear equation Lx = B for x, where L is a unit lower
//         triangular matrix.                                      
  
   x[0] = B[0];
   for (k = 1, L += n; k < n; L += n, k++) 
      for (i = 0, x[k] = B[k]; i < k; i++) x[k] -= x[i] * *(L + i);
}


////////////////////////////////////////////////////////////////////////////////
//  void Unit_Lower_Triangular_Inverse(double *L,  int n)                     //
//                                                                            //
//  Description:                                                              //
//     This routine calculates the inverse of the unit lower triangular       //
//     matrix L.  Only the subdiagonal part of the matrix is addressed.       //
//     The diagonal is assumed to consist of 1's and is not addressed.        //
//     The algorithm follows:                                                 //
//        Let M be the inverse of L, then L M = I,                            //
//          M[i][j] = -( L[i][j] M[j][j] + ... + L[i][i-1] M[i-1][j] ),       //
//     for i = 1, ..., n-1, j = 0, ..., i - 1.                                //
//                                                                            //
//  Arguments:                                                                //
//     double *L   On input, the pointer to the first element of the matrix   //
//                 whose unit lower triangular elements form the matrix which //
//                 is to be inverted. On output, the lower triangular part is //
//                 replaced by the inverse.  The diagonal and superdiagonal   //
//                 elements are not modified.                                 //
//     int     n   The number of rows and/or columns of the matrix L.         //
//                                                                            //
//  Return Values:                                                            //
//     void                                                                   //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double L[N][N];                                                        //
//                                                                            //
//     (your code to create the matrix L)                                     //
//     Unit_Lower_Triangular_Inverse(&L[0][0], N);                            //
//     printf(" The inverse is \n");                                          //
//           ...                                                              //
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
void Unit_Lower_Triangular_Inverse(double *L, int n)
{
   int i, j, k;
   double *p_i, *p_j, *p_k;

//         Invert the subdiagonal part of the matrix L row by row where
//         the diagonal elements are assumed to be 1.0.

   for (i = 1, p_i = L + n; i < n; i++, p_i += n) {
      for (j = 0, p_j = L; j < i; p_j += n, j++) {
         *(p_i + j) = -*(p_i + j);
         for (k = j + 1, p_k = p_j + n; k < i; k++, p_k += n)
            *(p_i + j) -= *(p_i + k) * *(p_k + j);
      }
   }
}
