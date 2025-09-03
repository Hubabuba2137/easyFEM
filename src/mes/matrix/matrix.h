#ifndef FEM_MATRIX_H
#define FEM_MATRIX_H

#include <vector>
#include <iostream>
#include <stdexcept>
#include <math.h>

namespace Fem {

    struct Matrix {
        std::vector<std::vector<double>> data;
        int rows;
        int cols;

        Matrix(int r, int c);
        Matrix(int r, int c, const std::vector<std::vector<double>>& init_data);

        int getRows() const;
        int getCols();

        std::vector<double>& operator[](int index);

        Matrix operator+(const Matrix& other) const;
        Matrix operator-(const Matrix& other) const;
        Matrix operator*(const Matrix& other) const;
        Matrix operator*(const float val) const;

        double determinant() const;
        Matrix transpose() const;
        bool isSquare() const;
        double trace() const;
        Matrix inverse();
        Matrix getMinor(int row, int col) const;
        Matrix multiply(float val);
    };

    std::ostream& operator<<(std::ostream& os, Matrix& mat);
    std::istream& operator>>(std::istream& is, Matrix& mat);

} // namespace Fem

#endif // FEM_MATRIX_H
