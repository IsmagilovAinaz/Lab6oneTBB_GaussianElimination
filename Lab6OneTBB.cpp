#include <iostream>
#include <tbb/tbb.h>
#include <chrono>
using namespace std;

int Rrand(int min, int max) {
    return 1 + rand() % (max - min + 1) + min;
}

double** CreateMatrix(int rows, int colums) {
    double** arr = new double* [rows];
    for (int i = 0; i < rows; i++) {
        arr[i] = new double[colums];
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < colums; j++) {
            arr[i][j] = 0;
        }
    }
    return arr;
}

double** MatrixFill(double** arr, int rows, int colums) {
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < colums; j++) {
            arr[i][j] = Rrand(-100, 100);
        }
    }
    return arr;
}

double* ArrFill(double* arr, int len) {
    for (int i = 0; i < len; i++) {
        arr[i] = Rrand(-100, 100);
    }
    return arr;
}

double* ArrFillZeros(double* arr, int len) {
    for (int i = 0; i < len; i++) {
        arr[i] = 0;
    }
    return arr;
}

void PrintMatrix(double** matrix, int rows, int colums) {
    if (matrix != NULL) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < colums; j++) {
                cout << matrix[i][j] << " ";
            }
            cout << endl;
        }
    }
    else {
        cout << "Matrix is empty" << endl;
    }
}

void PrintArr(double* arr, int n) {
    for (int i = 0; i < n; i++) {
        cout << arr[i] << " ";
    }
    cout << endl;
}

double* GaussTBB(double** A, double* Y, int n) {
    double* X = new double[n];

    for (int k = 0; k < n; k++)
    {
        tbb::parallel_for(tbb::blocked_range<int>(k + 1, n), [&](const tbb::blocked_range<int>& range) {
            for (int j = range.begin(); j < range.end(); j++) {
                double m = A[j][k] / A[k][k];
                for (int i = k; i < n; i++) {
                    A[j][i] = A[j][i] - m * A[k][i];
                }
                Y[j] = Y[j] - m * Y[k];
            }
        });
    }
    for (int k = n - 1; k >= 0; k--) {
        X[k] = Y[k];
        for (int i = k + 1; i < n; i++) {
            X[k] = X[k] - A[k][i] * X[i];

        }
        X[k] = X[k] / A[k][k];
    }
    return X;
}

double* FindY(double** A, double* X, int n) {
    double* Y = new double[n];
    for (int i = 0; i < n; i++) {
        Y[i] = 0;
        for (int j = 0; j < n; j++) {
            Y[i] += A[i][j] * X[j];
        }
    }
    return Y;
}

double* Gauss(double** A, double* Y, int n) {
    double* X = new double[n];
    for (int k = 0; k < n; k++)
    {
        for (int j = k + 1; j < n; j++) {
            double m = A[j][k] / A[k][k];
            for (int i = k; i < n; i++) {
                A[j][i] = A[j][i] - m * A[k][i];
            }
            Y[j] = Y[j] - m * Y[k];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        X[k] = Y[k];
        for (int i = k + 1; i < n; i++) {
            X[k] = X[k] - A[k][i] * X[i];

        }
        X[k] = X[k] / A[k][k];
    }

    return X;
}

int ArrEquial(double* X1, double* X2, int n, double eps) {
    int c = 0;
    for (int i = 0; i < n; i++) {
        if (std::abs(X1[i] - X2[i]) > eps) {
            c++;
        }
    }
    return c;
}

int EnterNumOfUnknown() {
    int nX = 0;
    do {
        cout << "Enter the Number of unknown quantities in the system: ";
        cin >> nX;
        if (cin.fail() || nX <= 0) {
            cin.clear();
            cin.ignore((numeric_limits<streamsize>::max)(), '\n');
            cout << "Incorrect input!" << endl;
        }
    } while (nX <= 0);
    return nX;
}

int main()
{
    //tbb::global_control control(tbb::global_control::max_allowed_parallelism, 6);
    int nX;
    nX = EnterNumOfUnknown();

    double** A = CreateMatrix(nX, nX);
    for (int i = 0; i < nX; i++) {
        A[i] = new double[nX];
    }
    for (int i = 0; i < nX; i++) {
        for (int j = 0; j < nX; j++) {
            A[i][j] = 0;
        }
    }
    double* X = new double[nX];
    double* X2 = new double[nX];
    double* Y = new double[nX];
    A = MatrixFill(A, nX, nX);
    X2 = ArrFill(X2, nX);
    //cout << "Matrix A:" << endl;
    //PrintMatrix(A, nX, nX);
    Y = FindY(A, X2, nX);
    //PrintArr(Y, nX);
    //cout << "Matrix Y:" << endl;
    //PrintArr(Y, nX);
    std::cout << "OneTBB solution" << std::endl;
    auto begin = std::chrono::steady_clock::now();
    X = GaussTBB(A, Y, nX);
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "fuction time: " << elapsed_ms.count() << " ms\n";
    int f = ArrEquial(X, X2, nX, 0.001);
    std::cout << "Number of mismatches:" << f <<std::endl;
    std::cout << "Direct solution" << std::endl;
    begin = std::chrono::steady_clock::now();
    X = Gauss(A, Y, nX);
    end = std::chrono::steady_clock::now();
    elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
    std::cout << "fuction time: " << elapsed_ms.count() << " ms\n";
    f = ArrEquial(X, X2, nX, 0.001);
    std::cout << "Number of mismatches:" << f << std::endl;
    for (int i = 0; i < nX; i++) {
        delete[] A[i];
    }
    delete[] A;
    delete[] X;
    delete[] Y;
}
