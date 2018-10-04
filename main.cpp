#include<bits/stdc++.h>
#include<eigen3/Eigen/Dense>

void mult_rec_helper(MatrixXf &matrix, MatrixXf &matrix1, MatrixXf c, int i, int i1, int i2);

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

MatrixXf mult_rec(MatrixXf &A,MatrixXf &B) {
    int N = A.rows();
    MatrixXf C = MatrixXf::Zero(N, N);


     mult_rec_helper(A,B,C,0,0,0);
}

void mult_rec_helper(MatrixXf &A, MatrixXf &B, MatrixXf C, int i, int j, int k) {
    if(A.rows()==1){
        C(i,j)=A(i,k)*B(k,j);

    }else{
        int half = A.rows()/2;
        mult_rec_helper(A.block(0,0,half,half), B.block(), MatrixXf C, int i, int j, int k);
        mult_rec_helper(MatrixXf &A, MatrixXf &B, MatrixXf C, int i, int j, int k);
        mult_rec_helper(MatrixXf &A, MatrixXf &B, MatrixXf C, int i, int j, int k);
        mult_rec_helper(MatrixXf &A, MatrixXf &B, MatrixXf C, int i, int j, int k);


    }

}




MatrixXf strassen(const MatrixXf &A, const MatrixXf &B) {
    int N = A.rows();
    MatrixXf C(N, N);

    //TODO: Point (e)

    return C;
}

int main() {
//    const int rows = 2;
//    const int cols = 2;
//    MatrixXf A(rows, cols);
//    A << 1., 2., 3.,
//            4., 5.;
//
//    MatrixXf B(rows, cols);
//    B << 1., 2., 3.,
//            4., 5.;
//
//    mult(A, B);


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
