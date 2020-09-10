//#include "armadillo"
#include "lib.h"
#include "time.h"
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
using namespace std;

const int DRAW_N = 10;
const int RUNS = 10; // amount of times to solve all

/*
solve the symmetric tridiagonal system with diagonal d, off-diagonal e and rhs
g.

The arrays d2 and g2 are used for intermediate results and are assumed to be
allocated outside this function(this will also make proper timing of this
function easier) u, which is also allocated outside this function, will contain
the solution

e: double[n-1] off-diagonal
d: double[n]   diagonal
d2:double[n]   intermediate diagonal
g: double[n]   rhs
g2:double[n]   intermediate rhs
u: double[n]   solution
n: int         system size
*/
void solve_toeplitz(double *e, double *d, double *d2, double *g, double *g2,
                    double *u, int n) {
  // initial conditions
  d2[0] = d[0];
  g2[0] = g[0];
  // forward substitution
  for (int i = 1; i < n; i++) {
    d2[i] = d[i] - e[i - 1] * e[i - 1] / d2[i - 1];
    g2[i] = g[i] - e[i - 1] * g2[i - 1] / d2[i - 1];
  }

  // backward subs
  for (int i = n - 2; i >= 1; i--) {
    u[i] = (g2[i] - e[i] * u[i + 1]) / d2[i];
  }
}

void solve_specific(double *d2, double *g, double *g2, double *u, int n) {
  // initial conditions
  d2[0] = -2;
  g2[0] = g[0];
  // forward substitution
  for (int i = 1; i < n; i++) {
    d2[i] = -(i + 1) / ((double)i);
    g2[i] = g[i] - g2[i - 1] / d2[i - 1];
  }
  // backwards substitution
  for (int i = n - 2; i >= 1; i--) {
    u[i] = (g2[i] - u[i + 1]) / d2[i];
  }
}

// write a (double) array to file
void write_array_to_file(ostream &file, double *array, int n) {
  int draw_distance = n / min(n, 50);
  for (int i = 0; i < n; i++) {
    if (i % draw_distance == 0)
      file << to_string(array[i]) << " ";
  }
  file << to_string(array[n - 1]) << "\n";
}

double f(double x) { return 100 * exp(-10 * x); }

/*
Lifted verbatim from
https://xoax.net/cpp/ref/cpp_examples/incl/mean_med_mod_array/
in the first example
*/
double median(double daArray[], int iSize) {
  // Allocate an array of the same size and sort it.
  double *dpSorted = new double[iSize];
  for (int i = 0; i < iSize; ++i) {
    dpSorted[i] = daArray[i];
  }
  for (int i = iSize - 1; i > 0; --i) {
    for (int j = 0; j < i; ++j) {
      if (dpSorted[j] > dpSorted[j + 1]) {
        double dTemp = dpSorted[j];
        dpSorted[j] = dpSorted[j + 1];
        dpSorted[j + 1] = dTemp;
      }
    }
  }

  // Middle or average of middle values in the sorted array.
  double dMedian = 0.0;
  if ((iSize % 2) == 0) {
    dMedian = (dpSorted[iSize / 2] + dpSorted[(iSize / 2) - 1]) / 2.0;
  } else {
    dMedian = dpSorted[iSize / 2];
  }
  delete[] dpSorted;
  return dMedian;
}

double compute_error(double *u_approx, double *u_exact, int n) {
  double epsilon;
  epsilon = 0;
  for (int i = 1; i < (n - 1); i++) {
    epsilon = max(epsilon, abs((u_approx[i] - u_exact[i]) / u_exact[i]));
  }
  return log10(epsilon);
}

/*
solve the system using all three methods,
if writing = true, the function writes to the file data_n=<n>.txt.
keeps a record of time usage for all three methods in times, which is assumed
to be able to hold three doubles.
*/
void solve_all(int n, double *times, bool writing) {
  double *e, *d, *d2, *g, *g2, *u, *u_exact;
  double x, relative_error;
  double h = 1 / ((double)n);
  int draw_distance = n / min(n, 50);

  // allocate e, d, g
  e = new double[n - 1];
  d = new double[n];
  d2 = new double[n];
  g = new double[n];
  g2 = new double[n];
  u = new double[n];
  u_exact = new double[n];

  // open file for writing, write n and log10h, descriptors in file removed for
  // convenience of reading later
  ofstream myfile;
  if (writing) {
    myfile.open("data_n=" + to_string(n) + ".txt");
    myfile << "" << to_string(n) << "\n";
    myfile << "" << to_string(log10(h)) << "\n";
  }

  // fill e, d, g and compute u_exact
  for (int i = 0; i <= n - 2; i++) {
    x = 0 + i * h;
    if (i % draw_distance == 0 && writing)
      myfile << to_string(x) << " ";

    e[i] = 1;
    d[i] = -2;
    g[i] = -h * h * f(x);
    u_exact[i] = 1 - (1 - exp(-10)) * x - exp(-10 * x);
  }
  // fix endpoints
  u_exact[n - 1] = 0;
  g[n - 1] = -h * h * f(1);
  d[n - 1] = -2;
  if (writing) { // write the last entries of the x row in the file
    myfile << to_string(1 - draw_distance * h);
    myfile << to_string(1) << "\n";
    // write u_exact to file after x
    write_array_to_file(myfile, u_exact, n);
  }

  //*****************************SOLVING****************************************
  // solve toeplitz
  auto start = chrono::high_resolution_clock::now();
  solve_toeplitz(e, d, d2, g, g2, u, n);
  auto finish = chrono::high_resolution_clock::now();
  times[0] =
      chrono::duration_cast<chrono::microseconds>(finish - start).count();
  relative_error = compute_error(u, u_exact, n);
  if (writing) {
    write_array_to_file(myfile, u, n);
    myfile << "" << to_string(relative_error) << "\n";
  }

  // cout << "elapsed time(allocation and initialization): "
  //       << (1000000 * (start - t0) / CLOCKS_PER_SEC) << " ms\n";
  /*cout << "standard method \n";
  cout << "Maximum relative error: " << relative_error << "\n";
  cout << "elapsed time: "
       << chrono::duration_cast<chrono::microseconds>(finish - start).count()
       << " ms\n";*/

  /// solve specific
  start = chrono::high_resolution_clock::now();
  solve_specific(d2, g, g2, u, n);
  finish = chrono::high_resolution_clock::now();
  relative_error = compute_error(u, u_exact, n);
  times[1] =
      chrono::duration_cast<chrono::microseconds>(finish - start).count();
  if (writing) {
    write_array_to_file(myfile, u, n);
    myfile << "" << to_string(relative_error) << "\n";
  }
  /*cout << "specific method \n";
  cout << "Maximum relative error: " << compute_error(u, u_exact, n) << "\n";
  cout << "elapsed time: "
       << chrono::duration_cast<chrono::microseconds>(finish - start).count()
       << " ms\n";*/

  // LU method, only viable for "small" n
  if (n <= 1000) {
    // allocate  to matrix and index array
    int *index = new int[n];
    double dd;
    double **matrix = new double *[n];
    for (int i = 0; i < n; i++) {
      matrix[i] = new double[n];
      for (int j = 0; j < n; j++) {
        matrix[i][j] = 0;
      }
      matrix[i][i] = -2;
      if (i > 0)
        matrix[i][i - 1] = 1;
      if (i < n)
        matrix[i][i + 1] = 1;
    }
    start = chrono::high_resolution_clock::now();
    ludcmp(matrix, n, index, &dd);
    lubksb(matrix, n, index, g);
    finish = chrono::high_resolution_clock::now();
    // cout << "Maximum relative error: " << compute_error(u, u_exact, n) <<
    // "\n";
    relative_error = compute_error(g, u_exact, n);
    times[2] =
        chrono::duration_cast<chrono::microseconds>(finish - start).count();
    if (writing) {
      write_array_to_file(myfile, g, n);
      myfile << "" << to_string(relative_error) << "\n";
    }
    // deallocate matrix
    for (int i = 0; i < n; i++) {
      delete[] matrix[i];
    }
    delete[] matrix;
  }

  // deallocate arrays
  delete[] e, d, d2, g, g2, u, u_exact;
  myfile.close();
}

int main(int argc, char *argv[]) {
  // declare variables
  int n;
  double times[3]; // an array conatining timing for all three methods
  double times_standard[RUNS];
  double times_specific[RUNS];
  double times_LU[RUNS];

  // cout << "This is c++ version " << __cplusplus << "\n";

  // ask for n
  cout << "input steps n:\n";
  cin >> n;
  // for each n, solve with all methods RUNS times,
  for (int i = 1; i <= n; i++) {
    cout << "------------------n=" << to_string(int(pow(10, i)))
         << "------------------\n";
    for (int run = 0; run < RUNS; run++) {
      if (run == RUNS - 1)                  // if last run
        solve_all(pow(10, i), times, true); // last run writes to file
      else
        solve_all(pow(10, i), times, false); // all other do not write to file
      // extract timing data from times
      times_standard[run] = times[0];
      times_specific[run] = times[1];
      times_LU[run] = times[2];
    }
    // take median of times and print
    cout << "median time for standard method: " << median(times_standard, RUNS)
         << " ms\n";
    cout << "median time for specific method: " << median(times_specific, RUNS)
         << " ms\n";
    cout << "median time for LU method: " << median(times_LU, RUNS) << " ms\n";
  }
}

/*
COMPILATION

g++ FYS3150_project_1.cpp lib.cpp -O0
*/

/*


>g++ -o FYS3150_project_1 FYS3150_project_1.cpp lib.cpp
>FYS3150_project_1.exe
input steps n:
7
------------------n=10------------------
median time for standard method: 0 ms
median time for specific method: 0 ms
median time for LU method: 0 ms
------------------n=100------------------
median time for standard method: 0 ms
median time for specific method: 0 ms
median time for LU method: 1985 ms
------------------n=1000------------------
median time for standard method: 0 ms
median time for specific method: 0 ms
median time for LU method: 2.91087e+006 ms
------------------n=10000------------------
median time for standard method: 496.5 ms
median time for specific method: 0 ms
median time for LU method: 3.26535e+006 ms
------------------n=100000------------------
median time for standard method: 6483 ms
median time for specific method: 2989.5 ms
median time for LU method: 3.26535e+006 ms
------------------n=1000000------------------
median time for standard method: 54438.5 ms
median time for specific method: 28424.5 ms
median time for LU method: 3.26535e+006 ms
------------------n=10000000------------------
median time for standard method: 1.28924e+006 ms
median time for specific method: 289226 ms
median time for LU method: 3.26535e+006 ms




O3 optimizations
>g++ -o FYS3150_project_1_O3 FYS3150_project_1.cpp lib.cpp -O3
>FYS3150_project_1_O3.exe
input steps n:
7
------------------n=10------------------
median time for standard method: 0 ms
median time for specific method: 0 ms
median time for LU method: 0 ms
------------------n=100------------------
median time for standard method: 0 ms
median time for specific method: 0 ms
median time for LU method: 996.5 ms
------------------n=1000------------------
median time for standard method: 0 ms
median time for specific method: 0 ms
median time for LU method: 1.05431e+006 ms
------------------n=10000------------------
median time for standard method: 0 ms
median time for specific method: 497 ms
median time for LU method: 1.1627e+006 ms
------------------n=100000------------------
median time for standard method: 6521 ms
median time for specific method: 3048.5 ms
median time for LU method: 1.1627e+006 ms
------------------n=1000000------------------
median time for standard method: 40891 ms
median time for specific method: 15458.5 ms
median time for LU method: 1.1627e+006 ms
------------------n=10000000------------------
median time for standard method: 2.17688e+006 ms
median time for specific method: 839194 ms
median time for LU method: 1.1627e+006 ms



*/
