#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

class Matrix {
public:
    int n;                                // количество строк
    int m;                                // количество столбцов

    vector<vector<float>> mat;

    Matrix operator+(Matrix b);
    Matrix operator-(Matrix b);
    Matrix operator*(Matrix b);
    Matrix operator*(int b);

    friend Matrix operator*(int a, Matrix b);

    Matrix();

    Matrix(int i, int j);         // матрица, заполненная нулями

    Matrix(int i, int j, vector<vector<float>> s);  
                                  // матрица из вектора векторов

    Matrix GaussianElimination(); // метод Гаусса,
                                  // либо получение нулевой строки/столбца

    float scalar(Matrix b);       // скалярное произведение векторов

    float maxnorm();

    void Set_Elems();             // ввод значений элементов матрицы

    float norm();                 // норма матрицы

    float Frobenius_norm();       // норма Фробениуса(евклидова норма) 

    float angle(Matrix b);        // угол в радианах

    int rank();                   // ранг матрицы

    Matrix minr(int, int);        // минор матрицы

    Matrix Transpose();           // транспонированная матрица

    Matrix Hadamard_product(Matrix b);

    friend ostream &operator<<(ostream &out, Matrix a);
    // вывод матрицы в поток

    friend ifstream& operator>>(ifstream& fin, Matrix& matr); 

    Matrix Inverse();             // обратная матрица

    void write(ofstream &);       // запись в бинарный файл

    void read(ifstream &);        // чтение бинарного файла
};

class IdentityMatrix : Matrix {
public:
    IdentityMatrix(int n) : Matrix(n, n) {
        vector<vector<double>> b;
        for (int i = 0; i < n; i++) mat[i][i] = 1;
        }
};

class DiagonalMatrix : Matrix{
public:
    DiagonalMatrix(int i, int j): Matrix(i, j) {};
};

class UpperMatrix : Matrix{
    UpperMatrix(int i, int j) : Matrix(i, j) {};
};

class LowerMatrix : Matrix{
    LowerMatrix(int i, int j) : Matrix(i,j) {};
};

class SymmetricMatrix : Matrix{
    SymmetricMatrix(int i) : Matrix(i, i) {};
};


float Matrix_trace(const Matrix A); // след матрицы

float det(const Matrix A);          // определитель

class DifferentSizesException {
public:
    string what = "Matrixes have different sizes";
};

class MultiplySizeException {
public:
    string what = "Column size in the 1st matrix != line size in the 2nd matrix";
};

class MatrixLineExceprion {
public:
    string what = "There is no such line";
};

class NotSquareMatrix {
public:
    string what = "Pushing not squared matrix";
};

class NotMatrix {
public:
    string what = "Pushing not matrix";
};

class SingularMatrixException {
public:
    string what = "Invertible matrix.";
};

class NotSymmetricMatrixException {
public:
    string what = "Entered not symmetrix matrix";
};

class NotUpperTriangleMatrixException {
public:
    string what = "Entered not upper triangle matrix";
};

class NotLowerTriangleMatrixException {
public:
    string what = "Entered not lower triangle matrix";
};

class WrongSizeHadamardException {
public:
    string what = "Wrong size for Hadamard product";
};

class NotVector {
public:
    string what = "This is not a vector";
};

class AngleError {
public:
    string what = "Angle not defined";
};

class NoOppositeMatrixException {
public:
    string what = "No opposite matrix! det(matrix) = 0!";
};

class FileNotFoundException {
public:
    string what = "Can't open/find file";
};
