#include <iostream>
#include <cmath>
#include <fstream>
#include "armadillo"
#include <iomanip>
#include <Windows.h>
#include <ctime>
#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace arma;


void solve_c(int n,double h,double* x,double* f, double *u, double *v);
void make_file(int n, double *u, double *v);
vec solve_LU(int n, double *x,double *f);
vec solve_LU_d(int n, double *x,double *f);
void make_file_LU(int n, double *x, vec v_LU);
int time(int n, double h,double *x,double *u, double *v,double *f, vec v_LU);
double* Reerror(double *a, double *b, int n);
void make_file_err(double* err_pointer, int n);


//main
int main(int argc, char* argv[])
{

    int n=0;

    if (argc <= 2) {
        cout << "Bad Usage: " << argv[0] <<
            " read also output file and n (int) on same line" << endl;
        exit(1);
    }
    else {
        n = atoi(argv[argc - 1]);  // last element
    }

    int i;
    // handling flag
    bool solved = 0;

    // array u & v& x & f
    double *u = new double[n + 2];
    double *v = new double[n + 2];
    double *x = new double[n + 2];
    double *f = new double[n + 1];
    // length of step
    double h;
    h = 1.0 / (n + 1);
    // initialize x elements
    for (i = 1; i <= n; i++) {
        x[i] = h*i;
    }
   x[0]=0;
   x[n+1]=1;

    // boundary condition
    v[0] = 0;
    v[n + 1] = 0;
    u[0] = 0;
    u[n+1] = 0;

    // fill f,the right part of equation
    for (i = 1; i <= n; i++) {
        f[i] = h*h * 100 * exp(-10 * x[i]);
    }
    f[0]=0;

    for( i = 1; i <= n; i++){
        u[i] = 1.0-(1-exp(-10))*x[i]-exp(-10*x[i]);
    }

    // a vector for lu
    vec v_LU = zeros<vec>(n + 2, 1);

    //pointer to store relative error
    double* relative_error_pointer = NULL;

    //maximum error
    double Max_relative_error;

    //flags
    //-s : general solution(general matrix LU decomposition)
    //-sc : optimised algorithm for tridiagonal matrix (task c)
    //-sLU : LU decomposition(matrix set as tridiagonal)
    //-err : get error

    // check flag
    for (i = 1; i < argc; i++) {
        if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("t") != string::npos)) {
             //Time flag (-t):
            time( n, h,x,u, v,f,v_LU);
            solved = 1;
        }
        if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("s") != string::npos) || argc == 2) {
            // -s
            solved = 1;
            if (string(argv[i]).find("sLU") != string::npos) {
                // LU decomposition (-sLU)
                v_LU.subvec(1, n) = solve_LU_d(n, x,f);
                make_file_LU(n, u, v_LU);

            }
            else if(string(argv[i]).find("sc") != string::npos){
                // c (-sc)
                solve_c(n, h, x, f, u, v);
                make_file(n, u, v);
                solved = 1;

            }
            else {
                // LU decomposition (-s)
                v_LU.subvec(1, n) = solve_LU(n, x,f);
                make_file_LU(n, u, v_LU);

            }
        }


        if (solved == 1 && (string(argv[i]).find("-") == 0 && string(argv[i]).find("err") != string::npos)) {
            //error calculate (-err)
            // calculate relative error and find max
            relative_error_pointer = new double[n+1];
            relative_error_pointer = Reerror(v, u, n);
            mat relative_error(1, n + 1);

            for (int i = 1; i <= n; i++)
            {
                relative_error(0, i) = relative_error_pointer[i];
            }

            Max_relative_error = relative_error_pointer[0];

            cout << "maximum error : " << Max_relative_error << endl;
            //make file
            make_file_err(relative_error_pointer, n);
        }
    }


    cout << "number of floating point numbers for optimized algorithm = " << 3 * (n - 1) << endl;
    cout << "number of floating point numbers for general algorithm = " << pow(n ,3)<< endl;
    cout << "check result in text file" << endl;
    system("pause");

    return 0;
}

//define function for optimized algorithm
void solve_c(int n,double h,double* x,double* f, double *u, double *v)
{
    int i;

    //double* x;
   // x = new double[n+2];
    double* d;
    d = new double[n+1];
    double* diag_temp = new double[n+1];
    double* a = new double[n+1];
    double* c = new double[n+1];
   // double* f;
   // f = new double[n+1];
   // double h;

    // length of step
   // h = 1.0 / (n + 1);

    // initialize d elements
    for (i = 1; i <= n; i++) {
        d[i] = 2;
        a[i] = -1;
        c[i] = -1;
    }
    c[n] = 0;
    a[1] = 0;

    double d_temp = d[1];
    v[1] = f[1]/d_temp;
    for(i = 2;i<=n;i++){
        diag_temp[i] = c[i-1]/d_temp;
        d_temp = d[i] - a[i]*diag_temp[i];
        v[i] = (f[i] - v[i-1]*a[i])/d_temp;
    }

    // backward substitution
    for(i = n-1; i >=1;i--){
        v[i] -= diag_temp[i+1]*v[i+1];
    }
}

//make file for u[i] and v[i] with optimized algorithm
void make_file(int n, double *u, double *v)
{
    int i;

    char filename[30];

    fstream myfile;
    sprintf_s(filename, "solution_c_n%d.txt", n);
    myfile.open(filename, ios::out);
    myfile << "u[i]        v[i]" << endl;
    for (i = 1; i <= n; i++) {
        myfile << u[i] << "         " << v[i] << endl;
    }
    myfile.close();

}

//define function for relative error calculation and finding max relative error(Yisha)
//a is for v[i] and b for u[i](closed-form)
double* Reerror(double *a, double *b, int n)
{
    double* reerr;
    reerr = new double[n + 1];
    for (int i = 1; i <= n; i++) {
        reerr[i] = log10(fabs((a[i] - b[i]) / b[i]));

        if (i == 1) {
            reerr[0] = reerr[i];
        }
        else {
            if (reerr[0] < reerr[i]) {
                reerr[0] = reerr[i];
            }
        }
    }
    return reerr;
}


void make_file_err(double* err_pointer, int n) {
    int i;

    char filename[20];

    fstream myfile;
    sprintf_s(filename, "error_list_n%d.txt", n);
    myfile.open(filename, ios::out);
    myfile << "maximum error : " << err_pointer[0] << endl;
    myfile << "error list : " << endl;
    for (i = 1; i <= n; i++) {
        myfile << setw(10) << err_pointer[i] << endl;
    }
    myfile.close();
}




vec solve_LU(int n, double *x,double* f)
{
    int i, j;
  //  double h;
    mat A = zeros<mat>(n, n);
    vec F(n);
    mat L;
    mat U;
    vec z;
    // fill matrix A
    for (i = 0; i<n; i++) {
        for (j = 0; j<n; j++) {
            A(i, j) = i + j + 3;
        }
    }
    
   // fill vec F to use function LU
   for (i = 1; i <= n; i++) {
       F[i-1] = f[i];
   }


    // Uv = z, Lz = f
    //nuber of floating operations = 4/3*n^3
    //round off
    if (n > 1000) {
        cout << "number of floating operations = " << pow(n, 3) << endl;
    }
    else {
        cout << "number of floating operations = " << ceilf(((float)4 / (float)3) * (float)(pow(n, 3))) << endl;
    }

    lu(L, U, A); // 2/3*n^2
    z = solve(L, F);
    return solve(U, z);
}




vec solve_LU_d(int n, double *x,double *f)
{
    int i, j;
    //double h;
    mat A = zeros<mat>(n, n);
    vec F(n);
    mat L;
    mat U;
    vec z;
    // fill matrix A
    for (i = 0; i<n; i++) {
        for (j = 0; j<n; j++) {
            if (i == j) {
                A(i, j) = 2;
            }
            else if (fabs(i - j) == 1) {
                A(i, j) = -1;
            }
        }
    }

   // fill f,the right part of equation
   for (i = 1; i <= n; i++) {
       F[i-1] = f[i];
   }
  // f[0]=0;

    // Uv = z, Lz = f
    lu(L, U, A);
    z = solve(L, F);
    return solve(U, z);
}


void make_file_LU(int n, double *u, vec v_LU)
{
    int i;
    char filename[30];

    fstream myfile;
    sprintf_s(filename, "solution_LUd_n%d.txt", n);
    myfile.open(filename, ios::out);
    myfile << "u[i]        v[i]" << endl;
    for (i = 1; i <= n ; i++) {
        myfile << u[i] << "      " << v_LU[i] << endl;
    }
    myfile.close();
}
int time(int n, double h,double *x,double *u, double *v,double *f, vec v_LU)
{
    high_resolution_clock::time_point start_c;
    high_resolution_clock::time_point start_LU;
    high_resolution_clock::time_point end_c;
    high_resolution_clock::time_point end_LU;
    start_c = high_resolution_clock::now();
    solve_c(n, h, x, f, u, v);
    end_c = high_resolution_clock::now();
    duration<double> time_span_c = end_c - start_c;
    cout << "n = " << n << endl;
    cout << "optimised solution took me " << setw(10) << time_span_c.count() << setw(8) << "seconds" << endl;
    start_LU = high_resolution_clock::now();
    v_LU = solve_LU_d(n, x,f);
    end_LU = high_resolution_clock::now();
    duration<double> time_span_LU = end_LU - start_LU;
    cout << "LU decomposition took me " << setw(12) << time_span_LU.count() << setw(8) << "seconds" << endl;
    return 0;
}
