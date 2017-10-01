#ifndef MATRIX_H
#define MATRIX_H
#include <vector>


template <typename T> class Matrix{
private:
	std::vector<std::vector<T>> matrix;
	unsigned rows; //num of rows
	unsigned cols; //numm of cols

public:
	//CONSTUCTORS
	Matrix(unsigned _rows, unsigned _cols, const T& _initial);
	Matrix(const std::vector<std::vector<T>>& _matrix);
	Matrix(const Matrix<T>& matrix2);
	~Matrix();

	//OPERATOR OVERLOADS
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

	//Gets element of matrix
	T& operator()(const unsigned &row, const unsigned &col);
	const T& operator()(const unsigned& row, const unsigned& col) const;

	//FUNCTIONS
	unsigned get_rows() const; //returns num of rows
	unsigned get_cols() const; //returns num of cols
	std::vector<T> diagonal(); //
	void printMatrix() const; //prints matrix in console
	Matrix<T> rref(); //Returns reduced row echelon form of matrix
	std::vector<T> getRow(unsigned i) const; //gets a row of matrix
	void replaceRow(unsigned rowNum, const std::vector<T>& newRow); //replaces row of matrix
	T trace(); //Get trace of matrix
	T det(); //Get matrix determint 
	Matrix<T> inverse(); //Get matrix inverse

};

#include "Matrix.cpp"
#endif