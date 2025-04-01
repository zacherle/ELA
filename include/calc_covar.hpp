
// Typedef for a matrix stored in column‐major order.
// A matrix is represented by a 1D vector of doubles. For a matrix with nr rows and nc columns,
// the element (i, j) (with 0-based indexing) is at index: i + j * nr.
typedef std::vector<double> Matrix;


namespace calc_covar {

    // Subroutine eigv_sym4
    // Computes the eigenvalues and eigenvectors of a 4x4 symmetric matrix using LAPACK's dsyev.
    // The matrix 'a' (dimension 4x4) is stored in column‐major order.
    // On exit, if successful, a contains the orthonormal eigenvectors and w holds the eigenvalues.
    void eigv_sym4(Matrix &a, double w[4]);

    // Function cov_matrix1
    // Computes the covariance matrix.
    // Parameters:
    //   m     : INTEGER, number of rows for resw, w, and Gw.
    //   n     : INTEGER, number of columns for Gw and dimension of the covariance matrix (n x n).
    //   resw  : REAL(real64) array of dimension m.
    //   Gw    : REAL(real64) matrix of dimension (m, n) stored in column‐major order.
    //   w     : REAL(real64) array of dimension m.
    // Returns:
    //   cov_matrix1 : REAL(real64) matrix of dimension (n, n) stored in column‐major order.
    Matrix cov_matrix1(int m, int n, const std::vector<double> &resw, const std::vector<double> &w, const Matrix &Gw);

    // Subroutine get_errellipse
    // Computes the error ellipse parameters from a 4x4 covariance matrix 'co' stored in column‐major order.
    // Output parameters:
    //   dxer, dyer : The square roots of the absolute values of co(1,1) and co(2,2).
    //   dzer, dter : The square roots of the absolute values of co(3,3) and co(4,4).
    //   l1, l2     : The lengths of the major and minor axes of the error ellipse.
    //   theta      : The orientation angle of the error ellipse (in degrees).
    void get_errellipse(const Matrix &co, double &dxer, double &dyer, double &dzer, double &dter,
                        double &l1, double &l2, double &theta);

} // end namespace calc_covar


