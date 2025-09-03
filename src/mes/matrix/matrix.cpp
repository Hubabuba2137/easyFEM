#include <vector>
#include <iostream>
#include <stdexcept>
#include <math.h>

#include "matrix.h"

namespace Fem
{
        //constructor
        Matrix::Matrix(int r, int c) : rows(r), cols(c) {
            if (r <= 0 || c <= 0) {
                throw std::invalid_argument("Dimensions must be positive integers.");
            }
            data.resize(r, std::vector<double>(c, 0.0));
        }

        //constructor with initializing values
        Matrix::Matrix(int r, int c, const std::vector<std::vector<double>>& init_data) : rows(r), cols(c) {
            if (r <= 0 || c <= 0) {
                throw std::invalid_argument("Dimensions must be positive integers.");
            }
            if (init_data.size() != r || init_data[0].size() != c) {
                throw std::invalid_argument("Initial data dimensions do not match specified dimensions.");
            }
            data = init_data;
        }

        int Matrix::getRows() const{
            return rows; 
        }

        int Matrix::getCols() { 
            return cols; 
        }

        //indexing element at given place []
        std::vector<double>& Matrix::operator[](int index) {
            if (index < 0 || index >= rows) {
                throw std::out_of_range("Row index out of bounds.");
            }
            return data[index];
        }

        //matrix addition (A+B)
        Matrix Matrix::operator+(const Matrix& other) const {
            if (rows != other.rows || cols != other.cols) {
                throw std::invalid_argument("Matrices must have the same dimensions for addition.");
            }

            Matrix result(rows, cols);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    result[i][j] = data[i][j] + other.data[i][j];
                }
            }
            return result;
        }

        //catrix subtraction (A-B)
        Matrix Matrix::operator-(const Matrix& other) const {
            if (rows != other.rows || cols != other.cols) {
                throw std::invalid_argument("Matrices must have the same dimensions for subtraction.");
            }

            Matrix result(rows, cols);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    result[i][j] = data[i][j] - other.data[i][j];
                }
            }
            return result;
        }

        //matrix multiplication (A*B)
        Matrix Matrix::operator*(const Matrix& other) const {
            if (cols != other.rows) {
                throw std::invalid_argument("Matrices cannot be multiplied. Columns of first matrix must match rows of second.");
            }

            Matrix result(rows, other.cols);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < other.cols; ++j) {
                    double sum = 0;
                    for (int k = 0; k < cols; ++k) {
                        sum += data[i][k] * other.data[k][j];
                    }
                    result[i][j] = sum;
                }
            }
            return result;
        }

        Matrix Matrix::operator*(const float val) const{
            Matrix result(rows, cols);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    result[i][j] = result[i][j]*val;
                }
            }
            return result;
        }

        //calculatinjg determinant
        double Matrix::determinant() const {
            if (!isSquare()) {
                throw std::invalid_argument("Determinant can only be calculated for square matrices.");
            }

            if (rows == 1) {
                return data[0][0];
            } else if (rows == 2) {
                return data[0][0] * data[1][1] - data[0][1] * data[1][0];
            } else {
                double det = 0;
                for (int col = 0; col < cols; ++col) {
                    Matrix minor = getMinor(0, col);
                    det += pow(-1, col) * data[0][col] * minor.determinant();
                }
                return det;
            }
        }

        //transpose the matrix
        Matrix Matrix::transpose() const {
            Matrix result(cols, rows);
            for (int i = 0; i < cols; ++i) {
                for (int j = 0; j < rows; ++j) {
                    result[i][j] = data[j][i];
                }
            }
            return result;
        }

        bool Matrix::isSquare() const { return rows == cols; }

        //calculate sum of diagonal elements
        double Matrix::trace() const {
            if (!isSquare()) {
                throw std::invalid_argument("Trace can only be calculated for square matrices.");
            }

            double sum = 0;
            for (int i = 0; i < rows; ++i) {
                sum += data[i][i];
            }
            return sum;
        }

        //inversing matrix
        Matrix Matrix::inverse() {
            //check if the matrix is square and has non-zero determinant
            if (!this->isSquare()) {
                throw std::invalid_argument("Inverse only exists for square matrices.");
            }

            double det = this->determinant();
            if (det == 0.0) {
                throw std::invalid_argument("Matrix is singular; inverse does not exist.");
            }

            int size = this->getRows();
            std::vector<std::vector<double>> cofactors(size, std::vector<double>(size));

            //calculate the cofactor matrix
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    Matrix minor = this->getMinor(i, j);
                    double minor_det = minor.determinant();
                    cofactors[i][j] = pow(-1, i + j) * minor_det;
                }
            }

            //transpose the cofactor matrix to get adjugate (adjoint)
            Matrix adjugate(size, size);
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    adjugate[j][i] = cofactors[i][j];
                }
            }

            //multiply each element of the adjugate by 1/determinant
            Matrix inverse(size, size);
            double inv_det = 1.0 / det;
            for (int i = 0; i < size; ++i) {
                for (int j = 0; j < size; ++j) {
                    inverse[i][j] = adjugate[i][j] * inv_det;
                }
            }

            return inverse;
        }

        //helper function for greating aaaa minor
        Matrix Matrix::getMinor(int row, int col) const {
            std::vector<std::vector<double>> minor_data;

            for (int i = 0; i < rows; ++i) {
                if (i == row) continue;
                std::vector<double> minor_row;
                for (int j = 0; j < cols; ++j) {
                    if (j == col) continue;
                    minor_row.push_back(data[i][j]);
                }
                minor_data.push_back(minor_row);
            }

            return Matrix(rows - 1, cols - 1, minor_data);
        }

        Matrix Matrix::multiply(float val){
            Matrix result(rows, cols);
            for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                    result[i][j] = result[i][j]*val;
                }
            }
            return result;
        }
}

    //overload << operator for printing
std::ostream& Fem::operator<<(std::ostream& os, Matrix& mat) {
    for (int i = 0; i < mat.getRows(); ++i) {
        for (int j = 0; j < mat.getCols(); ++j) {
            os << mat.data[i][j] << "   ";
        }
        os << std::endl;
    }
    return os;
}

    //overload >> operator for input
std::istream& Fem::operator>>(std::istream& is, Matrix& mat) {
    for (int i = 0; i < mat.getRows(); ++i) {
        for (int j = 0; j < mat.getCols(); ++j) {
            is >> mat.data[i][j];
        }
    }
    return is;
}