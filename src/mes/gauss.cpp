#include <vector>
#include "matrix.h"

std::vector<double> Gauss(const Fem::Matrix& A, std::vector<double>& b) {
    int n = A.getRows();
    Fem::Matrix augmentedMatrix = A;
    for (int i = 0; i < n; ++i) {
        augmentedMatrix[i].push_back(b[i]);
    }
    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (std::abs(augmentedMatrix[k][i]) > std::abs(augmentedMatrix[maxRow][i])) {
                maxRow = k;
            }
        }
        std::swap(augmentedMatrix[i], augmentedMatrix[maxRow]);
        double pivot = augmentedMatrix[i][i];
        for (int j = i; j <= n; ++j) {
            augmentedMatrix[i][j] /= pivot;
        }
        for (int k = i + 1; k < n; ++k) {
            double factor = augmentedMatrix[k][i];
            for (int j = i; j <= n; ++j) {
                augmentedMatrix[k][j] -= factor * augmentedMatrix[i][j];
            }
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        b[i] = augmentedMatrix[i][n];
        for (int j = i + 1; j < n; ++j) {
            b[i] -= augmentedMatrix[i][j] * b[j];
        }
    }

    return b;
}
