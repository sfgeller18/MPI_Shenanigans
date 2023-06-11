#include <iostream>
#include <cmath>
#include <mpi.h>
#include <gsl/gsl_sf_sqrt.h>

void JacobiCalc(float* A, float* x, float* b, int n) {
  double x_temp = new double[n];
  for (int i = 0; i < n; ++i) {
        double local_sum = 0.0;
        for (int j = 0; j < n; ++j) {
            if (j != i) {
                local_sum += A[i * n + j] * x[j];
            }
        }
  double general_sum;
  MPI_Allreduce(&local_sum, &general_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  x_temp[i] = (b[i] - general_Sum) / A[i * n + i];
    }
    for (int i = 0; i < n; ++i) {
        x[i] = x_temp[i];
    }
    delete[] x_temp;  
}

bool tolerance_check(float* A, float* x, float* b, int n, double e) {
  bool check = true;
  double local_res = 0.0;
  for (int i = 0; i < n; ++i) {
      double b_i = 0.0;
      for (int j = 0; j < n; ++j) {
          b_i += A[i * n + j] * x[j];
      }
      local_res += (b_i - b[i])*(b_i - b[i]);
  }
  double general_res;
  MPI_Allreduce(&res, &general_res, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  if (gsl_sf_sqrt(general_res) < e) {check = false;}
  return check;
}

float* JacobiMethod(float* A, float* x, float* b, int n, int iter, double e) {
  bool run = true;
  while (i<iter && run = true) {
    JacobiCalc(A, x, b, n);
    run = tolerance_check(A, x, b, n, e);
    ++i;
  }
  return x;
}
int main {
    MPI_Init(null, null);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    //Whatever you are looking to compute here
}
