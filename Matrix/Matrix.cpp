// Matrix.cpp : Console Application where the user can create matricies and do matrix mamth
// Resource used to make this: https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
#include "stdafx.h"

#pragma once
#ifndef MATRIX_CPP
#define MATRIX_CPP
#include "Matrix.h"
#include <iostream>
#include <stdexcept>  
#include <limits>  

#define SWAP(x,y) (x^=y),(y^=x),(x^=y)
using namespace std;

template<typename T>
Matrix<T>:: Matrix(const std::vector<std::vector<T>>& _matrix) {
	matrix = _matrix;
	rows = matrix.size();
	cols = matrix[0].size();
}

template<typename T>
Matrix<T>::Matrix(unsigned _rows, unsigned _cols, const T& _initial) {
	matrix.resize(_rows);
	for (unsigned i = 0; i < matrix.size(); i++) {
		matrix[i].resize(_cols, _initial);
	}
	rows = _rows;
	cols = _cols;
}

template<typename T>
Matrix<T>::Matrix(const Matrix& matrix2) {
	matrix = matrix2.matrix;
	rows = matrix2.get_rows();
	cols = matrix2.get_cols();
}

template<typename T>
Matrix<T>::~Matrix() {}


//MATRIX o MATRIX OVERLOADS
//Addition
template<typename T>
Matrix<T> Matrix<T>:: operator+(const Matrix<T>& matrix2){
	if (matrix2.get_rows() != get_rows() && matrix2.get_cols() != get_cols) {
		throw invalid_argument("Cannot add matricies");
	}
	Matrix results(rows, cols, 0.0);

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			results(i, j) = this->matrix[i][j] + matrix2(i, j);
		}
	}
	return results;
}

//Addition and equal
template<typename T>
Matrix<T>& Matrix<T>:: operator+=(const Matrix<T>& matrix2) {
	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			this->matrix[i][j] += matrix2(i, j);
		}
	}
	return *this;
}

//Subtraction
template<typename T>
Matrix<T> Matrix<T>:: operator-(const Matrix<T>& matrix2) {
	if (matrix2.get_rows() != get_rows() && matrix2.get_cols() != get_cols) {
		throw invalid_argument("Cannot subtract matricies");
	}

	Matrix results(rows, cols, 0.0);

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			results(i, j) = this->matrix[i][j] - matrix2(i, j);
		}
	}
	return c;
}

//Subtraction and equal
template<typename T>
Matrix<T>& Matrix<T>:: operator-=(const Matrix<T>& matrix2) {

	Matrix result = (*this) - matrix2;
	(*this) = result;
	return result;
}

//Multiplication
template<typename T>
Matrix<T> Matrix<T>:: operator*(const Matrix<T>& matrix2) {
	if (matrix2.get_rows() != get_cols()) {
		throw invalid_argument("Cannot multiply matricies");
	}

	unsigned cols = matrix2.get_cols();
	Matrix result(rows, matrix2.get_cols(), 0.0);

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			for (unsigned k = 0; k <= rows; k++) {
				result(i, j) += this->matrix[i][k] * matrix2(k, j);
			}
		}
	}
	return result;
}

//Multiplication and equal
template<typename T>
Matrix<T>& Matrix<T>:: operator*=(const Matrix<T>& matrix2) {
	Matrix result = (*this) * matrix2;
	(*this) = result;
	return *this;
}


//SCALAR o MATRIX OVERLOADS
//Scalar Addition
template<typename T>
Matrix<T> Matrix<T>:: operator+(const T& scalar) {
	Matrix result(rows, cols, 0.0); 

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			result(i, j) = this->matrix[i][j] + scalar;
		}
	}
	return result;
}

//Scalar Multiplication
template<typename T>
Matrix<T> Matrix<T>:: operator*(const T& scalar) {
	Matrix result(rows, cols, 0.0);

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			result(i, j) = this->matrix[i][j] * scalar;
		}
	}
	return result;
}

//Scalar Subtraction
template<typename T>
Matrix<T> Matrix<T>:: operator-(const T& scalar) {
	Matrix result(rows, cols, 0.0);

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			result(i, j) = this->matrix[i][j] - scalar;
		}
	}
	return result;
}

//Scalar Division
template<typename T>
Matrix<T> Matrix<T>:: operator/(const T& scalar) {
	Matrix result(rows, cols, 0.0);

	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			result(i, j) = this->matrix[i][j] / scalar;
		}
	}
	return result;
}

//FUNCTIONS
template<typename T>
std::vector<T> Matrix<T>::diagonal() {
	std::vector<T> result(rows, 0.0);

	for (unsigned i = 0; i < rows; i++) {
		result[i] = this->matrix[i][i];
	}

	return result;
}

template<typename T>
T& Matrix<T>:: operator()(const unsigned &row, const unsigned &col) {
	return this->matrix[row][col];
}

template<typename T>
const T& Matrix<T>::operator()(const unsigned &row, const unsigned &col) const {
	return this->matrix[row][col];
}

template<typename T>
unsigned Matrix<T>::get_rows() const {
	return this->rows;
}

template<typename T>
unsigned Matrix<T>::get_cols() const {
	return this->cols;
}

template<typename T>
void Matrix<T>::printMatrix() const {
	for (unsigned i = 0; i < rows; i++) {
		for (unsigned j = 0; j < cols; j++) {
			std::cout << (this->matrix[i][j]) << ", ";
		}
		std::cout << std::endl;
	}
}

//Returns Reduced Row Echelon Form of matrix.
template<typename T>
Matrix<T> Matrix<T>::rref() {
	Matrix result((*this)); //Copies matrix to result
	
	
	unsigned lead = 0;
	unsigned rowCount = result.get_rows();
	unsigned colCount = result.get_cols();

	//RREF
	while (lead < rowCount) {
		float d, m;

		for (int r = 0; r < rowCount; r++) {
			d = result(lead,lead);
			m = result(r, lead) / d;

			for (int c = 0; c < colCount; c++) {
				if (r == lead)
					result(r, c) /= d;
				else
					result(r, c) -= result(lead, c)*m;
			}
		}
		lead++;
	}

	//Changes '-0' to '0'
	for (unsigned i = 0; i < rowCount; i++) {
		for (unsigned j = 0; j < colCount; j++) {
			if (result(i, j) == -0) {
				result(i, j) = 0;
			}
		}
	}

	return result;
}


template<typename T>
std::vector<T> Matrix<T>::getRow(unsigned i)const {
	std::vector<T> aRow;
	aRow.resize(cols);
	for (unsigned j = 0; j < cols; j++) {
		aRow[j] = this->matrix[i][j];
	}
	return aRow;
}

template<typename T>
void Matrix<T>::replaceRow(unsigned rowNum, const std::vector<T>& newRow) {
	for (unsigned j = 0; j < cols; j++) {
		this->matrix[rowNum][j] = newRow[j];
	}
}

template<typename T>
T trace() {
	if (rows != cols)
		throw exception("Cannot trace non-square matrix");

	T value = 0;
	for (unsigned i = 0; i < rows; i++) {
		value += matrix[i][i];
	}
	return trace;
}

//MAIN FUNCTION
int main(int argc, char **argv) {
	//std::vector<std::vector<double>> v4 = ;
	Matrix<double> mat1({ { 1.0,1.0,1.0 },
						{ 2.0,2.0,2.0 } });
	Matrix<double> mat2({ { 1.0,2.0 },
						  { 1.0,2.0 },
						  { 1.0,2.0 }});

	Matrix<double> mat3= mat1*mat2;

	//mat3.printMatrix();
	mat3 += mat3;
	//mat3.printMatrix();
	mat3 = mat3 - 6;
	mat3 = mat3 * 4;
	mat3 = mat3 / 2;
	//mat3.printMatrix();


	Matrix<double> mat4 ({ {1,2,-1,-4},
					{2,3,-1,-11},
					{-2,0,-3,22} });

	mat4 = mat4.rref();
	mat4.printMatrix();

	return 0;
}

#endif // !__MATRIX_CPP
