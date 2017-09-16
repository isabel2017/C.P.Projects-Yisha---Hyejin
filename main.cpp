/*
**     Project1: a) and c)
**     The algorithm for solving the tridiagonal matrix
**     Relative error calculation and max value is found.
**     equation is implemented (requiering O(8n) FLOPS).
**     Results are written to textfile and read by a python
**     script (project1_b_plot.py) to make plots.
*/
# include <iostream>
# include <fstream>
# include <iomanip>
# include <cmath>
#include<armadillo>

using namespace std;
using namespace arma;
ofstream ofile;

// Declaring two functions that will be used:
double Solution(double x) {return 1.0-(1-exp(-10))*x-exp(-10*x);}// closed form solution
double f(double x) {return 100*exp(-10*x);}

//decalring relative error calculation Function
double* Reerror(double *a, double *b, int n)
{
    double* reerr;
    reerr= new double[n+1];
    for (int i=1;i<=n;i++){
        reerr[i]=log10(fabs((a[i]-b[i])/b[i]));
    }

//    reerr = new double[n+1];
//    for (int i=1;i<=n;i++){
//        reerr[i]=log10(fabs((a[i]-b[i])/b[i]));
//        if(i != 1){
//            if(reerr[0] < reerr[i]){
//                reerr[0] = reerr[i];
//            }
//        }
//        else{
//            reerr[0] = reerr[i];
//        }
//    }
    return reerr;
}


// Main program reads filename and n from command line:
int main(int argc, char* argv[]) {

    // Declaration of initial variables:
    char *outfilename;
    int n;
   // double m[4];
   // m[0] = (double)atoi(argv[3]);

    // Read in output file and n,
    // abort if there are too few command-line arguments:
    if( argc <= 2 ){
      cout << "Bad Usage: " << argv[0] <<
          " read also output file and n (int) on same line" << endl;
      exit(1);
    }
    else{
      outfilename = argv[1]; // first command line argument.
      n = atoi(argv[2]);  // second command line argument.
    }

    // Constants of the problem:
    double h = 1.0/(n+1.0);
    double *x = new double[n+2];
    double *b_twidd = new double[n+1]; // construction with n+1 points to make
                                       // indexing close to mathematics.
    b_twidd[0] = 0;

    // The constituents of the tridiagonal matrix A:
    // Zeroth element not needed, but included to make indexing easy:
    int *a = new int[n+1];
    int *b = new int[n+1];
    int *c = new int[n+1];

    // Temporal variabel in Gaussian elimination:
    double *diag_temp = new double[n+1];

    // Real solution and approximated one:
    double *u = new double[n+2]; // Analytical solution
    double *v = new double[n+2]; // Numerical solution
    // Including extra points to make the indexing easy:
    u[0] = 0;
    v[0] = 0;

    // Filling up x-array, making x[0] = 0 and x[n+1] = 1:
    for (int i=0; i<=n+1; i++) {
        x[i] = i*h;
        // Could print results to check:
        //cout << "x = " << x[i] << " and " << "h^2*f(x) = " << h*h*f(x[i]) << endl;
    }

    // Filling up b_twiddle array, i.e. right hand side of equation:
    for (int i=1; i<=n; i++) {
        b_twidd[i] = h*h*f(x[i]);
        // Could print here to check:
        //cout << "b_twidd = " << b_twidd[i] << "for x = " << x[i] << endl;
        u[i] = Solution(x[i]);
        //cout << "u = " << u[i] << " for x = " << x[i] <<  endl;
        b[i] = 2;
        a[i] = -1;
        c[i] = -1;
    }
    c[n] = 0;
    a[1] = 0;

    // Algorithm for finding v:
    // a(i)*v(i-1) + b(i)*v(i) + c(i)*v(i+1) = b_twidd(i)
    // Row reduction; forward substitution:
    double b_temp = b[1];
    v[1] = b_twidd[1]/b_temp;//???
    for (int i=2;i<=n;i++) {
       // Temporary value needed also in next loop:
       diag_temp[i] = c[i-1]/b_temp;
       // Temporary diagonal element:
       b_temp = b[i] - a[i]*diag_temp[i];
       // Updating right hand side of matrix equation:
       v[i] = (b_twidd[i]-v[i-1]*a[i])/b_temp;
    }

    // Row reduction; backward substition:
    for (int i=n-1;i>=1;i--) {
        v[i] -= diag_temp[i+1]*v[i+1];
    }

    // calculate relative error and find max
    double* relative_error_pointer = new double[n];
    relative_error_pointer= Reerror(u,v,n);
    mat relative_error(1,n+1);

    for(int i=1;i<=n;i++)
    {
        relative_error(0,i)=relative_error_pointer[i];
    }

     double Max_relative_error=max(max(relative_error));


    // Open file and write results to file:
    ofile.open(outfilename);
    ofile << setiosflags(ios::showpoint | ios::uppercase);
    ofile << "       x:             u(x):          v(x):      relative_error(x):         Max_relative_error(x):   " << endl;
    for (int i=1;i<=n;i++) {
       ofile << setw(15) << setprecision(8) << x[i];
       ofile << setw(15) << setprecision(8) << u[i];
       ofile << setw(15) << setprecision(8) << v[i];
       ofile << setw(15) << setprecision(8) << relative_error[i];
       if(i == 1){
           ofile << setw(30) << setprecision(8) << Max_relative_error<<endl;
       }
       else{
           ofile << endl;
       }

    }

    ofile.close();



    delete [] x;
    delete [] b_twidd;
    delete [] a;
    delete [] b;
    delete [] c;
    delete [] u;
    delete [] v;
    delete [] relative_error_pointer;
    return 0;
}
