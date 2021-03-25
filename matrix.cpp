#include "matrix.h"

using namespace std;

Matrix::Matrix(int i, int j) {
    n = i;
    m = j;
    vector<float> str;
    for (int l = 0; l < i; l++) {
        for (int k = 0; k < j; k++)
            str.push_back(0);
        mat.push_back(str);
    }
}

Matrix::Matrix() : Matrix(1, 1){
  
}

Matrix Matrix::operator+(Matrix b) {
        if (n != b.n or m != b.m) {
            throw DifferentSizesException();
        }

        Matrix c(n, m);
        for (int l = 0; l < n; l++) {
            for (int k = 0; k < m; k++) {
                c.mat[l][k] = mat[l][k] + b.mat[l][k];
            }
        }
        return c;
}

Matrix Matrix::operator-(Matrix b) {
        if (n != b.n or m != b.m) {
            throw DifferentSizesException();
        }

        Matrix c(n, m);
        
        for (int l = 0; l < n; l++) {
            for (int k = 0; k < m; k++) {
                c.mat[l][k] = mat[l][k] - b.mat[l][k];
            }
        }
        return c;
}


Matrix Matrix::operator*(Matrix b) {
        if (n != b.m) {
            throw DifferentSizesException();
        }

        Matrix c(n, m);
        for (int l = 0; l < n; l++) {
            for (int k = 0; k < m; k++) {
                float tmp = 0;
                for (int s = 0; s < n; s++) {
                    tmp += mat[l][s] * b.mat[s][k];
                }
                c.mat[l][k] = tmp;

            }
        }
        return c;
}


Matrix Matrix::operator*(int b) {


    Matrix c(n, m);
    for (int l = 0; l < n; l++) {
        for (int k = 0; k < m; k++) {
            c.mat[l][k] = b * mat[l][k];
        }
    }
    return c;
}

Matrix operator*(int a, Matrix b) {
    Matrix c(b.n, b.m);
    for (int l = 0; l < b.n; l++) {
        for (int k = 0; k < b.m; k++) {
            c.mat[l][k] = a * b.mat[l][k];
        }
    }
    return c;
}

Matrix Matrix::Hadamard_product(Matrix b) {
        if (n != b.n or m != b.m)
            throw WrongSizeHadamardException();
        Matrix c(b.n, b.m);
        for (int l = 0; l < b.n; l++) {
            for (int k = 0; k < b.m; k++) {
                c.mat[l][k] = mat[l][k] * b.mat[l][k];
            }
        }
        return c;
}

float Matrix_trace(const Matrix A) {
    int min = 0;
    if (A.n < A.m) {
        min = A.n;
    } else {
        min = A.m;
    }
    float trace = 0;

    for (int i = 0; i < min; i++) {
        trace += A.mat[i][i];
    }
    return trace;
}

Matrix Matrix::GaussianElimination() {
    Matrix matr(n, m, mat);
    vector<float> tmp;
    float max;
    int maxi;
    int min;
    min = matr.m < matr.n ? matr.m : matr.n;
    for (int j = 0; j < min; j++) {
        max = 0;
        maxi = j;
        for (int s = j; s < n; s++) {
            if (abs(matr.mat[s][j]) > max) {
                max = abs(matr.mat[s][j]);
                maxi = s;
            }
        }
        if (max == 0) {
            return matr;
        }
        vector<float> tmp = matr.mat[j];
        matr.mat[j] = matr.mat[maxi];
        matr.mat[maxi] = tmp;

        for (int k = j + 1; k < n; k++) {
            float r = matr.mat[k][j];
            for (int l = j; l < m; l++) {
                matr.mat[k][l] = matr.mat[k][l] - matr.mat[j][l] / matr.mat[j][j] * r;

                mat[k][l];
            }
        }

    }
    return matr;
}


ostream &operator<<(ostream &out, Matrix a) {
        if (!out)
            throw FileNotFoundException();
        for (int i = 0; i < a.n; i++) {
            for (int j = 0; j < a.m; j++) {
                out << a.mat[i][j] << " ";
            }
            out << '\n';
        }
        return out;
}

Matrix::Matrix(int i, int j, vector<vector<float>> s) {
    n = i;
    m = j;
    for (auto v: s) {
        mat.push_back(v);
    }
}

float det(const Matrix A) {
        if (A.n != A.m) {
            throw NotSquareMatrix();;
        }
        Matrix M = A;
        M = M.GaussianElimination();
        float D = 1;
        for (int i = 0; i < M.n; i++) {
            D *= M.mat[i][i];
        }
        return D;
}

float Matrix::scalar(Matrix b) {
        if (not((m == 1 or n == 1) and (m == b.m and n == b.n))) {
            throw MultiplySizeException();
        }
        if (m == 1) {
            float S = 0;
            for (int i = 0; i < n; i++) {
                S += mat[i][0] * b.mat[i][0];
            }
            return S;
        } else {
            float S = 0;
            for (int j = 0; j < m; j++) {
                S += mat[0][j] * b.mat[0][j];
            }
            return S;
        }
}

float Matrix::norm() {
    return sqrt((this->scalar(*this)));
}

float Matrix::maxnorm() {
        if (not((m == 1 or n == 1))) {
            throw NotVector();
        }
        if (m == 1) {
            float S = abs(mat[0][0]);
            for (int i = 0; i < n; i++) {
                if (S < abs(mat[i][0])) {
                    S = abs(mat[i][0]);
                }
            }
            return S;
        } else {
            float S = abs(mat[0][0]);
            for (int j = 0; j < m; j++) {
                if (S < abs(mat[0][j])) {
                    S = abs(mat[0][j]);
                }

            }
            return S;
        }
}

float Matrix::Frobenius_norm() {
    float s = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            s += mat[i][j] * mat[i][j];
        }
    return sqrt(s);
}

float Matrix::angle(Matrix b) {
        if (this->norm() == 0 or b.norm() == 0)throw AngleError();
        float anglecos;
        anglecos = (this->scalar(b)) / (this->norm() * b.norm());
        float angle;
        angle = acos(anglecos);
        return angle;
}

int Matrix::rank() {
    Matrix M = this->GaussianElimination();
    int r = 0;
    for (int i = 0; i < n; i++) {
        bool iszero = false;
        for (int j = 0; j < m; j++) {
            if (abs(M.mat[i][j]) > 0.01) iszero = true;
        }
        r = iszero == true ? r + 1 : r;
    }
    return r;
}

Matrix Matrix::Transpose() {
    Matrix M(m, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            M.mat[j][i] = mat[i][j];
        }
    return M;
}


Matrix Matrix::minr(int k, int l) {
        if (k < 1 or k > n or l < 1 or l > m)
            throw MatrixLineExceprion();

        k -= 1;
        l -= 1;
        Matrix M(n - 1, m - 1);
        int del = 0;
        int i = 0, j = 0, q = 0, r = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                if (i == k or j == l) { continue; }
                if (i > k and j < l)
                    M.mat[i - 1][j] = mat[i][j];
                else if (i < k and j < l)
                    M.mat[i][j] = mat[i][j];
                if (i < k and j > l)
                    M.mat[i][j - 1] = mat[i][j];
                if (i > k and j > l)
                    M.mat[i - 1][j - 1] = mat[i][j];
            }
        return M;
}

Matrix Matrix::Inverse() {
        Matrix T = this->Transpose();
        float d = det(T);
        double epsilon = 0.001;
        if (abs(d) < epsilon) throw NoOppositeMatrixException();
        if (T.rank() < n) throw SingularMatrixException();
        Matrix I(n, m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                I.mat[i][j] = pow(-1, i + j) * det(T.minr(i + 1, j + 1)) / d;
            }
        return I;
}

ifstream& operator>>(ifstream& fin,  Matrix& mat) {
    string line;
     int n = 0;
     int m = 0;
    double val;
    vector<vector<float>> vvtmp;
    if (fin.is_open()) {
        while (getline(fin, line)) {
            m = 0;
            n++;
            vector<float> vtmp;
            istringstream iss(line);
            while (!iss.eof()) {
                m++;
                iss >> val;
                vtmp.push_back(val);
            }
            vvtmp.push_back(vtmp);
            line.clear();
        }
        mat.mat = vvtmp;
        mat.n = n;
        mat.m = m;
    }
    else    throw FileNotFoundException();
    
    return fin;
}

void Matrix::write(ofstream &a) {
        if (!a)
            throw FileNotFoundException();
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                void *k = new float[1];
                float r;
                r = mat[i][j];
                k = (void *) &r;
                a.write((char *) k, sizeof(float));
            }
}

void Matrix::read(ifstream &a) {
        if (!a.is_open())throw FileNotFoundException();
        void *k = new float[1];
        float r;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++) {
                a.read((char *) k, sizeof(float));
                mat[i][j] = ((float *) k)[0];
            }
}

void Matrix::Set_Elems() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            cout << "Enter value arr["<< i << "][" << j << "]:" << endl;
            cin >> mat[i][j];
        }
    }
  cout << '\n';
}