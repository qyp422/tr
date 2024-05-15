#include <iostream>
#include <cmath>

const int MAX_SIZE = 100;

// 矩阵乘法
void matrixMultiply(const double A[MAX_SIZE][MAX_SIZE], const double B[MAX_SIZE][MAX_SIZE], double C[MAX_SIZE][MAX_SIZE], int m, int n, int p) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < p; ++j) {
            C[i][j] = 0;
            for (int k = 0; k < n; ++k) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

// 转置矩阵
void transpose(const double A[MAX_SIZE][MAX_SIZE], double At[MAX_SIZE][MAX_SIZE], int m, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            At[i][j] = A[j][i];
        }
    }
}

// 矩阵的特征值分解
void eigenvalueDecomposition(const double A[MAX_SIZE][MAX_SIZE], double eigenvalues[MAX_SIZE], int n) {
    // 这里使用简化的方法，假设 A 是实对称矩阵，直接对角化
    // 实际上，特征值分解是一个复杂的数值计算问题
    // 对于大型矩阵，可以使用更复杂的算法，如雅可比方法
    for (int i = 0; i < n; ++i) {
        eigenvalues[i] = A[i][i];
    }
}

// 主函数，奇异值分解
void svd(const double A[MAX_SIZE][MAX_SIZE], double U[MAX_SIZE][MAX_SIZE], double S[MAX_SIZE], double Vt[MAX_SIZE][MAX_SIZE], int m, int n) {
    // A 的转置矩阵
    double At[MAX_SIZE][MAX_SIZE];
    transpose(A, At, m, n);

    // A * A^T 的乘积
    double AAt[MAX_SIZE][MAX_SIZE];
    matrixMultiply(A, At, AAt, m, n, m);

    // 计算 A * A^T 的特征值
    double eigenvalues[MAX_SIZE];
    eigenvalueDecomposition(AAt, eigenvalues, m);

    // 计算奇异值
    for (int i = 0; i < m; ++i) {
        S[i] = std::sqrt(std::abs(eigenvalues[i])); // 奇异值是特征值的平方根
    }

    // 计算 Vt
    // 这里假设 V 是从特征向量构成的，实际上需要更复杂的计算
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < m; ++j) {
            Vt[i][j] = At[i][j] / S[j];
        }
    }

    // 计算 U
    // U = A * V * S^(-1)
    double S_inv[MAX_SIZE][MAX_SIZE] = {0};
    for (int i = 0; i < m; ++i) {
        if (S[i] != 0) {
            S_inv[i][i] = 1 / S[i];
        }
    }
    double AV[MAX_SIZE][MAX_SIZE];
    matrixMultiply(A, Vt, AV, m, n, m);
    matrixMultiply(AV, S_inv, U, m, m, m);
}

// 打印矩阵
void printMatrix(const double A[MAX_SIZE][MAX_SIZE], int m, int n) {
    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
}

int main() {
    // 例子矩阵 A
    double A[MAX_SIZE][MAX_SIZE] = {
        {1, 2, 3},
        {4, 5, 6},
        {7, 8, 9}
    };

    // 存储奇异值分解结果的矩阵
    double U[MAX_SIZE][MAX_SIZE], Vt[MAX_SIZE][MAX_SIZE];
    double S[MAX_SIZE] = {0};

    // 进行奇异值分解
    svd(A, U, S, Vt, 3, 3);

    // 打印结果
    std::cout << "U:\n";
    printMatrix(U, 3, 3);

    std::cout << "S:\n";
    for (int i = 0; i < 3; ++i) {
        std::cout << S[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "Vt:\n";
    printMatrix(Vt, 3, 3);

    return 0;
}
