#ifndef MATRIX_H
#define MATRIX_H
#include <vector>


template <typename T> class Matrix{
private:
	std::vector<std::vector<T>> matrix;
	unsigned rows;
	unsigned cols;

public:
	Matrix(unsigned _rows, unsigned _cols, const T& _initial);
	Matrix(const std::vector<std::vector<T>>& _matrix);
	Matrix(const Matrix<T>& matrix2);
	~Matrix();

	//Matrix o Matrix
	Matrix<T> operator+(const Matrix<T>& matrix2);
	Matrix<T>& operator+=(const Matrix<T>& matrix2);
	Matrix<T> operator-(const Matrix<T>& matrix2);
	Matrix<T>& operator-=(const Matrix<T>& matrix2);
	Matrix<T> operator*(const Matrix<T>& matrix2);
	Matrix<T>& operator*=(const Matrix<T>& matrix2);
	Matrix<T> transpose();
	//Scalar o Matrix 
	Matrix<T> operator+(const T& scalar);
	Matrix<T> operator-(const T& scalar);
	Matrix<T> operator*(const T& scalar);
	Matrix<T> operator/(const T& scalar);
	T& operator()(const unsigned &row, const unsigned &col);
	const T& operator()(const unsigned& row, const unsigned& col) const;

	//Methods
	unsigned get_rows() const;
	unsigned get_cols() const;
	std::vector<T> diagonal();
	void printMatrix() const;
};

#include "Matrix.cpp"
#endif