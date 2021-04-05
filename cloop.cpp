#include <omp.h>
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>

// Store all Livermore loops in a vector for easier construcion of
// automated tests for each loop.
std::vector<std::function<void(void)>> loop_vec;

// Set global constants.
const int NUM_OF_THREADS = 8; // Number of threads used.
const int NUM_OF_ITER = 100; // Number of iterations.

// Get an arithmetic sequence vector using a set of
// given size and differential.
template <class T>
std::vector<T> GetVec(const int size, const float diff)
{
  std::vector<T> temp(size);
  #pragma omp parallel for
  for (int i = 0; i < temp.size(); ++i) {
    temp[i] = i*diff;
  }
  return temp;
}

// Get a 2D arithmetic sequence vector.
template <class T>
std::vector<std::vector<T>>
GetVec2D(const int size, const float diff)
{
  std::vector<std::vector<T>> temp(size, std::vector<T>(size, 0));
  #pragma omp parallel for collapse(2)
  for (int i = 0; i < temp.size(); ++i) {
    for (int j = 0; j < temp[0].size(); ++j) {
      temp[i][j] = i*j*diff;
    }
  }
  return temp;
}

// Get a 3D uniform vector.
template <class T>
std::vector<std::vector<std::vector<float>>>
GetUniformVec3D(const int x_size, const int y_size,
  const int z_size, int num)
{
  return std::vector<std::vector<std::vector<float>>>
    (x_size, std::vector<std::vector<float>>
      (y_size, std::vector<float>(z_size, num)));
}

// Get the sum of a vector using an OpenMP approach.
template <class T> 
double GetVecSum(const std::vector<T> input)
{
  double sum = 0;
  #pragma omp parallel for reduction(+:sum)
  for (int i = 0; i < input.size(); ++i) {
    sum += input[i];
  }
  return sum;
}

// Initialize the vector containing all Livermore loops.
void LoopInitialize()
{
  // Kernel 1 -- Hydro Fragment
  loop_vec.push_back([]() -> void {
    int i = 0, n = 1'000'000;
    float q = 0.05, r = 0.02, t = 0.01;
    auto x{GetVec<float>(n, 0.001)};
    auto y{GetVec<float>(n, 0.0003)};
    auto z{GetVec<float>(n, 0.0005)};

    do {
      #pragma omp parallel for
      for (int k = 0; k < n; k++) {
        x[k] = q + y[k]*(r*z[k + 10] + t*z[k + 11]);
      }
    } while (++i < NUM_OF_ITER);
    auto sum = GetVecSum(x);
  });
  
  // Kernel 2 -- ICCG Excerpt (Incomplete Cholesky Conjugate Gradient)
  loop_vec.push_back([]() -> void {
    int ipntp, ipnt, i, ii;
    int n = 10'000, j = 0;
    auto x{GetVec<float>(2*n, -0.002)};
    auto v{GetVec<float>(2*n, 0.0007)};
    do {
      ii = n;
      ipntp = 0;
      do {
        ipnt = ipntp;
        ipntp += ii;
        ii /= 2;
        i = ipntp ;
        // This loop is not suitable for omp parallel for
        // since x[i] and x[k] have dependency on the
        // previous loop.
        for (int k = ipnt + 1; k < ipntp; k = k + 2) {
          ++i;
          x[i] = x[k] - v[k]*x[k - 1] - v[k + 1]*x[k + 1];
        }
      } while (ii > 0);
    } while (++j < NUM_OF_ITER);
    // Loop No.2 can still benefit from the paralleled sum.
    auto sum = GetVecSum(x);
  });

  // Kernel 3 -- Inner Product
  loop_vec.push_back([]() -> void{
    double res = 0;
    int i = 0, n = 1'000'000;
    auto x{GetVec<float>(n, 0.0001)};
    auto z{GetVec<float>(n, 0.0015)};

    do {
      #pragma omp parallel for reduction(+:res)
      for (int k = 0; k < n; ++k) {
        res += z[k]*x[k];
      }
    } while (++i < NUM_OF_ITER);
  });

  // Kernel 4 -- Banded Linear Equations
  loop_vec.push_back([]() -> void {
    int j = 0, n = 1'000'000;
    int m = (1001 - 7)/2;
    auto x{GetVec<float>(n, 0.01)};
    auto y{GetVec<float>(n, 0.0035)};

    do {
      for (int k = 6; k < 10'000; k += m) {
        int lw = k - 6;
        float temp = x[k - 1];
        #pragma omp parallel for reduction(-:temp)
        for (int j = 4; j < n; j += 5) {
          temp -= x[lw]*y[j];
          ++lw;
        }
        x[k - 1] = y[4]*temp;
      }
    } while (++j < NUM_OF_ITER);
    auto sum = GetVecSum(x);
  });

  // Kernel 5 -- Tri-Diagonal Elimination (Below Diagonal)
  loop_vec.push_back([]() -> void {
    int j = 0, n = 10'000;
    std::vector<float> x(n, 0);
    auto y{GetVec<float>(n, 0.0000305)};
    auto z{GetVec<float>(n, 0.0000023)};
    do {
      // This loop can't be parallelized for x[i] has
      // dependency on x[i - 1].
      for (int i = 1; i < n; ++i) {
        x[i] = z[i]*(y[i] - x[i - 1]);
      }
    } while (++j < NUM_OF_ITER);
    auto sum = GetVecSum(x);
  });

  // Kernel 6 -- General Linear Recurrence Equations
  loop_vec.push_back([]() -> void {
    int j = 0, n = 100;
    auto w{GetVec<float>(n, 0.00000012)};
    auto b{GetVec2D<float>(n, 0.000007)};
    do {
      // Not suitable for parallelization since w[i]
      // have dependency on w[(i-k)-1].
      for (int i = 1; i < n; ++i) {
        w[i] =  0.0100;
        for (int k = 0; k < i; ++k) {
          w[i] += b[k][i]*w[(i - k) - 1];
        }
      }
    } while (++j < NUM_OF_ITER);
    auto sum = GetVecSum(w);
  });

  // Kernel 7 -- Equation of State Fragment
  loop_vec.push_back([]() -> void {
    int i = 0, n = 1'000'000;
    float q = 0.5, r = 0.2, t = 0.1;
    auto x{GetVec<float>(n, 0.0001)};
    auto y{GetVec<float>(n, 0.00023)};
    auto z{GetVec<float>(n, 0.0016)};
    auto u{GetVec<float>(n, 0.0021)};

    do {
      #pragma omp parallel for
      for (int k = 0; k < n; ++k) {
        x[k] = u[k] + r*(z[k] + r*y[k]) + t*(u[k + 3] 
          + r*(u[k + 2] + r*u[k + 1]) + t*(u[k + 6]
          + q*(u[k + 5] + q*u[k + 4])));
      }
    } while (++i < NUM_OF_ITER);
    auto sum = GetVecSum(x);
  });

  // Kernel 8 -- ADI Integration
  loop_vec.push_back([]() -> void {
    int j = 0, n = 100'000;
    auto u1{GetUniformVec3D<float>(2, n + 1, 4, 1.25)};
    auto u2{GetUniformVec3D<float>(2, n + 1, 4, 0.75)};
    auto u3{GetUniformVec3D<float>(2, n + 1, 4, 1.05)};
    std::vector<float> du1(n);
    std::vector<float> du2(n);
    std::vector<float> du3(n);
    int nl1 = 0, nl2 = 1;
    float a11 = 1.1, a12 = 1.2, a13 = 1.3;
    float a21 = 2.1, a22 = 2.2, a23 = 2.3;
    float a31 = 3.1, a32 = 3.2, a33 = 3.3, sig = 1.0;

    do {
      #pragma omp parallel for collapse(2)
      for (int kx = 1; kx < 3; ++kx) {
        for (int ky = 1; ky < n; ++ky) {
          du1[ky] = u1[nl1][ky + 1][kx] - u1[nl1][ky - 1][kx];
          du2[ky] = u2[nl1][ky + 1][kx] - u2[nl1][ky - 1][kx];
          du3[ky] = u3[nl1][ky + 1][kx] - u3[nl1][ky - 1][kx];
          u1[nl2][ky][kx] = u1[nl1][ky][kx] + a11*du1[ky] + a12*du2[ky]
            + a13*du3[ky] + sig*(u1[nl1][ky][kx + 1] - 2.0*u1[nl1][ky][kx]
            + u1[nl1][ky][kx - 1]);
          u2[nl2][ky][kx] = u2[nl1][ky][kx] + a21*du1[ky] + a22*du2[ky]
            + a23*du3[ky] + sig*(u2[nl1][ky][kx + 1] - 2.0*u2[nl1][ky][kx]
            + u2[nl1][ky][kx - 1]);
          u3[nl2][ky][kx] = u3[nl1][ky][kx] + a31*du1[ky] + a32*du2[ky]
            + a33*du3[ky] + sig*(u3[nl1][ky][kx + 1] - 2.0*u3[nl1][ky][kx]
            + u3[nl1][ky][kx - 1]);
        }
      }
    } while (++j < NUM_OF_ITER);
  });
}

int main(int argc, char const *argv[])
{
  // Initialize all Livermore loops.
  LoopInitialize();

  // Automated test for each loop in the loop vector.
  for (int i = 0; i < loop_vec.size(); ++i) {
    std::cout << "Loop No." << i + 1 << " is executing.\n";

    // Single thread execution.
    omp_set_num_threads(1);
    auto single_elapsed_time = -omp_get_wtime();
    loop_vec[i]();
    single_elapsed_time += omp_get_wtime();
    // std::cout << "The single thread version takes "
    //   << single_elapsed_time << " seconds.\n";

    // Multi-thread execution.
    omp_set_num_threads(NUM_OF_THREADS);
    auto multi_elapsed_time = -omp_get_wtime();
    loop_vec[i]();
    multi_elapsed_time += omp_get_wtime();

    // Calculate speedup of parallelization.
    // std::cout << "The paralleled version takes "
    //   << multi_elapsed_time << " seconds.\n";
    std::cout << "Speedup of loop No." << i + 1 << " is: "
      << single_elapsed_time/multi_elapsed_time << "\n\n";
  }

  return 0;
}
