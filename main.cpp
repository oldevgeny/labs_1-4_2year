#include "matrix.h"

using namespace std;

int main() {
  try{
    cout << "ЛАБОРАТОРНАЯ РАБОТА 1.\n" << endl;


    Matrix A(2, 2);
    //A.Set_Elems();
    A.mat[0][0] = 4;
    A.mat[0][1] = 2;
    A.mat[1][0] = 3;
    A.mat[1][1] = 1;

    Matrix B(2, 2);
    //B.Set_Elems();
    B.mat[0][0] = 4;
    B.mat[0][1] = 1;
    B.mat[1][0] = 1;
    B.mat[1][1] = 3;

    Matrix C(2, 2);

    cout << "A:\n" << A << endl;

    cout << "B:\n" << B << endl;

    C = A * B;
    cout << "A * B: \n" << C << endl;

    C = A + B;
    cout << "A + B: \n" << C << endl;

    C = A - B;
    cout << "A - B: \n" << C << endl;

    C = A * 2;
    cout << "A * q: (q = 2)\n" << C << endl;

    C = A.Hadamard_product(B);
    cout << "Произведение Адамара: \n" << C << endl;



    cout << "ЛАБОРАТОРНАЯ РАБОТА 2.\n" << endl;


    cout << "След матрицы A: \n" << Matrix_trace(A) << '\n' << endl;


    vector<vector<float>> s{{1, -2, 3},
                            {4, 0, 6},
                            {-7, 8, 9}};
    Matrix S(3, 3, s);
    cout << "Матрица S:\n" << S << endl;

    cout << "Определитель матрицы S методом Гаусса:" << endl;
    
    Matrix s1 = S.GaussianElimination();
    cout << s1 << "Определитель матрицы S равен " << det(s1) << ".\n" << endl;


    vector<vector<float>> b{{1, 2, 4}};
    vector<vector<float>> c{{1, 3, 2}};
    Matrix k1(1, 3, b);
    Matrix k2(1, 3, c);
    cout << "Вектор b = " << k1 << "Вектор c = " << k2 << endl;
    float sk = k1.scalar(k2);
    cout << "Скалярное произведение векторов b и c = " << sk << "." << endl;

    cout << "Евклидова норма вектора b = " << k1.norm() << "." << endl;
    cout << "Максимальная норма вектора b = " << k1.maxnorm() << "." << endl;
    cout << "Евклидова норма вектора с = " << k2.norm() << "." << endl;
    cout << "Максимальная норма вектора с = " << k2.maxnorm() << "." << endl;

    cout << "Норма Фробениуса матрицы A = " << A.Frobenius_norm() << "." << endl;



    cout << "ЛАБОРАТОРНАЯ РАБОТА 3.\n" << endl;


    cout << "Угол между векторами b и c = " << k1.angle(k2) << ".\n" << endl;

    cout << "Ранг матрицы A = " << A.rank() << "." << endl;
    cout << "Ранг матрицы S = " << S.rank() << ".\n" << endl;

    cout << "Матрица, обратная A:\n" << A.Inverse() << endl;
    cout << "Проверка: (A*A^(-1))\n" << A * A.Inverse() << endl;

    cout << "Транспонированная матрица A:\n" << A.Transpose() << endl;



    cout << "ЛАБОРАТОРНАЯ РАБОТА 4.\n" << endl;


    cout << "Читаем файл input.txt ..." << endl;
    ifstream f;
    f.open("input.txt");
    Matrix R;
    f >> R;
    cout << "Выводим прочитанное: " << endl;
    cout << R << endl;

    ofstream out;
    out.open("output.txt");
    cout << "Записываем вектор b в файл output.txt ..." << endl;
    out << k2;
    cout << "Записываем матрицу S в файл output.txt ..." << endl;
    out << S;
    out.close();

    ofstream outstrm("file.bin", ios::binary);
    cout << "Записываем вектор b в бинарный файл ..." << endl;
    k1.write(outstrm);
    cout << "Записываем матрицу S в бинарный файл ..." << endl;
    S.write(outstrm);
    outstrm.close();

    cout << "Читаем бинарный файл ..." << endl;
    ifstream instr("file.bin", ios::binary);
    k1.read(instr);
    cout << "Выводим считанный из бинарного файла вектор:" << endl;
    cout << k1;
    S.read(instr);
    cout << "Выводим считанную из бинарного файла матрицу:" << endl;
    cout << S << endl;
    instr.close();
  }
  catch(DifferentSizesException error){
    cout << error.what << endl;
  }
  catch(MultiplySizeException error){
    cout << error.what << endl;
  }
  catch(MatrixLineExceprion error){
    cout << error.what << endl;
  }
  catch(NotSquareMatrix error){
    cout << error.what << endl;
  }
  catch(NotMatrix error){
    cout << error.what << endl;
  }
  catch(SingularMatrixException error){
    cout << error.what << endl;
  }
  catch(NotSymmetricMatrixException error){
    cout << error.what << endl;
  }
  catch(NotUpperTriangleMatrixException error){
    cout << error.what << endl;
  }
  catch(NotLowerTriangleMatrixException error){
    cout << error.what << endl;
  }
  catch(WrongSizeHadamardException error){
    cout << error.what << endl;
  }
  catch(NotVector error){
    cout << error.what << endl;
  }
  catch(AngleError error){
    cout << error.what << endl;
  }
  catch(NoOppositeMatrixException error){
    cout << error.what << endl;
  }
  catch(FileNotFoundException error){
    cout << error.what << endl;
  }
}