#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <numeric>
#include <lapacke.h>

// Define constants from layers and const_raddeg modules
const double model_error = 0.2;    // value assumed from layers module
const double reading_error = 0.1;  // value assumed from layers module
const double RAD2DEG = 180.0 / M_PI; // constant from const_raddeg module

// Typedef for a matrix stored in column‐major order.
// A matrix is represented by a 1D vector of doubles. For a matrix with nr rows and nc columns,
// the element (i, j) (with 0-based indexing) is at index: i + j * nr.
typedef std::vector<double> Matrix;

// Inline helper functions for accessing and setting matrix elements in column‐major order.
inline double get(const Matrix &A, int nr, int i, int j) {
    return A[i + j * nr];
}
inline void set(Matrix &A, int nr, int i, int j, double value) {
    A[i + j * nr] = value;
}

// Helper function: Multiply two square matrices A and B (dimension n x n) stored in column‐major order.
Matrix matmul_square(const Matrix &A, const Matrix &B, int n) {
    Matrix C(n * n, 0.0);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++) {
                sum += get(A, n, i, k) * get(B, n, k, j);
            }
            set(C, n, i, j, sum);
        }
    }
    return C;
}

// Helper function: Compute matrix product GG = MATMUL(TRANSPOSE(Gw), Gw)
// Gw is an m x n matrix stored in column‐major order. The result GG is n x n.
Matrix matmul_transpose(const Matrix &Gw, int m, int n) {
    Matrix GG(n * n, 0.0);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            double sum = 0.0;
            for (int k = 0; k < m; k++) {
                // In Gw, element (k, i) is at index k + i * m, (k, j) is at k + j * m.
                sum += Gw[k + i * m] * Gw[k + j * m];
            }
            set(GG, n, i, j, sum);
        }
    }
    return GG;
}

namespace calc_covar {

    // Subroutine eigv_sym4
    // Computes the eigenvalues and eigenvectors of a 4x4 symmetric matrix using LAPACK's dsyev.
    // The matrix 'a' (dimension 4x4) is stored in column‐major order.
    // On exit, if successful, a contains the orthonormal eigenvectors and w holds the eigenvalues.
    void eigv_sym4(Matrix &a, double w[4])
    {
        const int n = 4;
        const int lda = n;
        int info;
        int lwork = -1;
        double work_query;
        // Workspace query for optimal lwork
        info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'L', n, a.data(), lda, w, &work_query, lwork);
        lwork = static_cast<int>(work_query);
        std::vector<double> work(lwork);
        // Evaluate eigenvalues and eigenvectors
        info = LAPACKE_dsyev_work(LAPACK_COL_MAJOR, 'V', 'L', n, a.data(), lda, w, work.data(), lwork);
        if (info > 0) {
            std::cout << "The algorithm failed to compute eigenvalues." << std::endl;
            std::exit(EXIT_FAILURE);
        }
    }

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
    Matrix cov_matrix1(int m, int n, const std::vector<double> &resw, const std::vector<double> &w, const Matrix &Gw)
    {
        //--- Begin function cov_matrix1 ---
        // For this implementation, n must equal 4 as eigv_sym4 operates on a 4x4 matrix.
        if (n != 4) {
            std::cerr << "cov_matrix1 error: n must be 4." << std::endl;
            std::exit(EXIT_FAILURE);
        }
        const double kappa = 1.0e-8;
        // Compute GG = MATMUL(TRANSPOSE(Gw), Gw)
        Matrix GG = matmul_transpose(Gw, m, n); // GG is n x n.
        Matrix eigV = GG;  // Copy of GG to be used for eigen-decomposition.
        double eig_l[4];
        // Compute eigenvalues and eigenvectors of eigV (which is GG)
        eigv_sym4(eigV, eig_l);
        // Adjust eigenvalues if their ratio is below kappa.
        double min_eig = eig_l[0], max_eig = eig_l[0];
        for (int i = 1; i < 4; i++) {
            if (eig_l[i] < min_eig)
                min_eig = eig_l[i];
            if (eig_l[i] > max_eig)
                max_eig = eig_l[i];
        }
        if (min_eig / max_eig < kappa) {
            for (int i = 0; i < 4; i++) {
                eig_l[i] = eig_l[i] + kappa;
            }
        }
        // Build pmat such that pmat(j,:) = eigV(:,j)/eig_l(j)
        // Here, we set the j-th row of pmat from the j-th column of eigV.
        Matrix pmat(4 * 4, 0.0);
        for (int j = 0; j < 4; j++) {
            for (int i = 0; i < 4; i++) {
                // In eigV, element (i,j) is accessed as get(eigV, 4, i, j)
                // We assign this to pmat at position (j,i).
                double val = get(eigV, 4, i, j) / eig_l[j];
                set(pmat, 4, j, i, val);
            }
        }
        // Compute covariance matrix in eigenvalue form: covM = MATMUL(eigV, pmat)
        Matrix covM = matmul_square(eigV, pmat, 4);
        // Compute sumw = sum of elements of w
        double sumw = std::accumulate(w.begin(), w.end(), 0.0);
        // Count number of positive elements in w
        int marrv = std::count_if(w.begin(), w.end(), [](double val) { return val > 0.0; });
        // Compute variance factor varD from resw and w.
        double dot = 0.0;
        for (int i = 0; i < m; i++) {
            dot += resw[i] * resw[i];
        }
        double varD = dot / sumw;
        if (marrv > 4) {
            varD = varD * marrv / (marrv - 4);
        } else {
            varD = varD * marrv;
        }
        // reading error: ensure varD is not below reading_error².
        if (varD < reading_error * reading_error) {
            varD = reading_error * reading_error;
        }
        // model error is added to the variance.
        varD = varD + model_error * model_error;
        // Multiply the covariance matrix by the variance factor.
        for (auto &val : covM) {
            val *= varD;
        }
        return covM;
        //--- End function cov_matrix1 ---
    }

    // Subroutine get_errellipse
    // Computes the error ellipse parameters from a 4x4 covariance matrix 'co' stored in column‐major order.
    // Output parameters:
    //   dxer, dyer : The square roots of the absolute values of co(1,1) and co(2,2).
    //   dzer, dter : The square roots of the absolute values of co(3,3) and co(4,4).
    //   l1, l2     : The lengths of the major and minor axes of the error ellipse.
    //   theta      : The orientation angle of the error ellipse (in degrees).
    void get_errellipse(const Matrix &co, double &dxer, double &dyer, double &dzer, double &dter,
                        double &l1, double &l2, double &theta)
    {
        //--- Begin subroutine get_errellipse ---
        // Compute the determinant from the 2x2 upper-left block of co.
        double deter = get(co, 4, 0, 0) * get(co, 4, 1, 1) - get(co, 4, 0, 1) * get(co, 4, 1, 0);
        double d11 = get(co, 4, 1, 1) / deter;
        double d22 = get(co, 4, 0, 0) / deter;
        double d21 = -get(co, 4, 1, 0) / deter;
        theta = std::atan(2 * d21 / (d11 - d22)) / 2.0;
        double cos_theta = std::cos(theta);
        double sin_theta = std::sin(theta);
        double al = d11 * cos_theta * cos_theta + 2 * d21 * cos_theta * sin_theta + d22 * sin_theta * sin_theta;
        double bl = d11 * sin_theta * sin_theta - 2 * d21 * cos_theta * sin_theta + d22 * cos_theta * cos_theta;
        l1 = std::sqrt(1.0 / al);
        l2 = std::sqrt(1.0 / bl);
        // Convert theta to degrees.
        theta = theta * RAD2DEG;
        // Adjust if l2 is greater than l1.
        if (l2 > l1) {
            double tl = l1;
            l1 = l2;
            l2 = tl;
            theta = theta + 90.0;
        }
        dxer = std::sqrt(std::fabs(get(co, 4, 0, 0)));
        dyer = std::sqrt(std::fabs(get(co, 4, 1, 1)));
        dzer = std::sqrt(std::fabs(get(co, 4, 2, 2)));
        dter = std::sqrt(std::fabs(get(co, 4, 3, 3)));
        //--- End subroutine get_errellipse ---
    }

} // end namespace calc_covar

// Example main function demonstrating usage (can be removed or modified as needed)
int main() {
    // Example usage of eigv_sym4:
    // Define a 4x4 symmetric matrix (column-major order)
    Matrix A = {
        4.0, 1.0, 2.0, 0.0,  // first column
        1.0, 3.0, 0.0, 1.0,  // second column
        2.0, 0.0, 2.0, 1.0,  // third column
        0.0, 1.0, 1.0, 1.0   // fourth column
    };
    double w[4];
    calc_covar::eigv_sym4(A, w);
    std::cout << "Eigenvalues:" << std::endl;
    for (int i = 0; i < 4; i++) {
        std::cout << w[i] << std::endl;
    }

    // Example usage of cov_matrix1:
    // Assume m = 5 and n = 4.
    int m = 5, n = 4;
    // Example resw and w vectors.
    std::vector<double> resw = {0.5, 0.6, 0.7, 0.8, 0.9};
    std::vector<double> warr = {1.0, 2.0, 3.0, 4.0, 5.0};
    // Example Gw: an m x n matrix stored in column-major order.
    // For simplicity we initialize Gw with 5 rows and 4 columns.
    Matrix Gw = {
        // Column 1:
        1.0, 2.0, 3.0, 4.0, 5.0,
        // Column 2:
        2.0, 3.0, 4.0, 5.0, 6.0,
        // Column 3:
        3.0, 4.0, 5.0, 6.0, 7.0,
        // Column 4:
        4.0, 5.0, 6.0, 7.0, 8.0
    };
    Matrix covMat = calc_covar::cov_matrix1(m, n, resw, warr, Gw);
    std::cout << "Covariance Matrix (4x4):" << std::endl;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            std::cout << get(covMat, n, i, j) << " ";
        }
        std::cout << std::endl;
    }

    // Example usage of get_errellipse:
    // Use covMat as the 4x4 covariance matrix "co".
    double dxer, dyer, dzer, dter, l1, l2, theta;
    calc_covar::get_errellipse(covMat, dxer, dyer, dzer, dter, l1, l2, theta);
    std::cout << "Error Ellipse Parameters:" << std::endl;
    std::cout << "dxer = " << dxer << ", dyer = " << dyer << std::endl;
    std::cout << "dzer = " << dzer << ", dter = " << dter << std::endl;
    std::cout << "l1 (major axis) = " << l1 << ", l2 (minor axis) = " << l2 << std::endl;
    std::cout << "theta = " << theta << " degrees" << std::endl;

    return 0;
}

