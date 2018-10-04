#include<bits/stdc++.h>
#include<eigen3/Eigen/Dense>


using namespace std;
using namespace Eigen;

MatrixXf mult(const MatrixXf &A, const MatrixXf &B) {
    int N = A.rows();
    MatrixXf C = MatrixXf::Zero(N, N);

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                C(i, j) += A(i, k) * B(k, j);
            }
        }
    }


//    cout << C << endl;
    return C;
}

MatrixXf mult_rec_helper(MatrixXf A, MatrixXf B) {
    if (A.rows() == 1) {
        return A * B;

    } else {
        int half = A.rows() / 2;


        MatrixXf A11 = A.block(0, 0, half, half);
        MatrixXf A12 = A.block(0, half, half, half);
        MatrixXf A21 = A.block(half, 0, half, half);
        MatrixXf A22 = A.block(half, half, half, half);

        MatrixXf B11 = B.block(0, 0, half, half);
        MatrixXf B12 = B.block(0, half, half, half);
        MatrixXf B21 = B.block(half, 0, half, half);
        MatrixXf B22 = B.block(half, half, half, half);


        MatrixXf C11 = mult_rec_helper(A11, B11) + mult_rec_helper(A12, B21);
        MatrixXf C12 = mult_rec_helper(A11, B12) + mult_rec_helper(A12, B22);
        MatrixXf C21 = mult_rec_helper(A21, B11) + mult_rec_helper(A22, B21);
        MatrixXf C22 = mult_rec_helper(A21, B12) + mult_rec_helper(A22, B22);


        MatrixXf C = MatrixXf::Zero(A.rows(), A.cols());
        C.block(0, 0, half, half) = C11;
        C.block(0, half, half, half) = C12;
        C.block(half, 0, half, half) = C21;
        C.block(half, half, half, half) = C22;

        return C;


    }

}

MatrixXf mult_rec(MatrixXf &A, MatrixXf &B) {
    int N = A.rows();
    MatrixXf C = MatrixXf::Zero(N, N);
    C = mult_rec_helper(A, B);
//    cout << C << endl;
    return C;
}

MatrixXf strassen_helper(MatrixXf A, MatrixXf B) {
    if (A.rows() == 1) {
        return A * B;

    } else {
        int half = A.rows() / 2;


        MatrixXf A11 = A.block(0, 0, half, half);
        MatrixXf A12 = A.block(0, half, half, half);
        MatrixXf A21 = A.block(half, 0, half, half);
        MatrixXf A22 = A.block(half, half, half, half);

        MatrixXf B11 = B.block(0, 0, half, half);
        MatrixXf B12 = B.block(0, half, half, half);
        MatrixXf B21 = B.block(half, 0, half, half);
        MatrixXf B22 = B.block(half, half, half, half);


        MatrixXf M1 = strassen_helper((A11 + A22), (B11 + B22));
        MatrixXf M2 = strassen_helper((A21 + A22), B11);
        MatrixXf M3 = strassen_helper((A11), (B12 - B22));
        MatrixXf M4 = strassen_helper((A22), (B21 - B11));
        MatrixXf M5 = strassen_helper((A11 + A12), (B22));
        MatrixXf M6 = strassen_helper((A21 - A11), (B11 + B12));
        MatrixXf M7 = strassen_helper((A12 - A22), (B21 + B22));


        MatrixXf C11 = M1 + M4 - M5 + M7;
        MatrixXf C12 = M3 + M5;
        MatrixXf C21 = M2 + M4;
        MatrixXf C22 = M1 - M2 + M3 + M6;


        MatrixXf C = MatrixXf::Zero(A.rows(), A.cols());
        C.block(0, 0, half, half) = C11;
        C.block(0, half, half, half) = C12;
        C.block(half, 0, half, half) = C21;
        C.block(half, half, half, half) = C22;

        return C;


    }
}
MatrixXf strassen(const MatrixXf &A, const MatrixXf &B) {
    int N = A.rows();
    MatrixXf C = MatrixXf::Zero(N, N);
    C = strassen_helper(A, B);




    return C;
}

int main() {


    srand(time(0));
    cout << setprecision(6) << setfill(' ');

    for (int i = 1; i < 9; i++) {
        int N = 1 << i;
        cout << "Matrix size = " << N << endl;
        MatrixXd AA = MatrixXd::Random(N, N);
        MatrixXd BB = MatrixXd::Random(N, N);
        MatrixXd ans = AA * BB;
        MatrixXf A = AA.cast<float>();
        MatrixXf B = BB.cast<float>();

        auto start = std::chrono::steady_clock::now();
        MatrixXf W = mult(A, B);
        auto finish = std::chrono::steady_clock::now();
        cout << setw(24) << " " << setw(15)
             << "Time (s)" << setw(20) << "Error (l2-norm)" << endl;
        cout << setw(24) << "Naive iterative " << setw(15)
             << std::chrono::duration_cast<std::chrono::duration<double> >
                     (finish - start).count() << setw(20) << (W.cast<double>() - ans).norm() << endl;

        start = std::chrono::steady_clock::now();
        MatrixXf X = mult_rec(A, B);
        finish = std::chrono::steady_clock::now();
        cout << setw(24) << "Naive recursive " << setw(15)
             << std::chrono::duration_cast<std::chrono::duration<double> >
                     (finish - start).count() << setw(20) << (X.cast<double>() - ans).norm() << endl;
        start = std::chrono::steady_clock::now();
        MatrixXf Y = strassen(A, B);
        finish = std::chrono::steady_clock::now();
        cout << setw(24) << "Strassen recursive " << setw(15)
             << std::chrono::duration_cast<std::chrono::duration<double> >
                     (finish - start).count() << setw(20) << (Y.cast<double>() - ans).norm() << endl;
        start = std::chrono::steady_clock::now();
        MatrixXf Z = A * B;
        finish = std::chrono::steady_clock::now();
        cout << setw(24) << "Eigen built-in " << setw(15)
             << std::chrono::duration_cast<std::chrono::duration<double> >
                     (finish - start).count() << setw(20) << (Z.cast<double>() - ans).norm() << "\n\n\n";
    }
}

