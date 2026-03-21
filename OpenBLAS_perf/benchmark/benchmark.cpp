#include <iostream>
#include <vector>
#include <chrono>
#include <cmath>
#include <cstring>
#include <omp.h>

extern "C" {
#include <cblas.h>
}

#include "../src/trsm.h"


double time_my(int n, int m, double* A, double* B) {
    auto t1 = std::chrono::high_resolution_clock::now();

    trsm_double_parallel(Left, Lower, NoTrans, NonUnit, n, m, 1.0, A, B);

    auto t2 = std::chrono::high_resolution_clock::now();

    double time =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

    if (time < 0) time = 0;
    return time;
}

double time_blas(int n, int m, double* A, double* B) {
    auto t1 = std::chrono::high_resolution_clock::now();

    cblas_dtrsm(
        CblasRowMajor,
        CblasLeft,
        CblasLower,
        CblasNoTrans,
        CblasNonUnit,
        n, m, 1.0,
        A, n,
        B, m
    );

    auto t2 = std::chrono::high_resolution_clock::now();

    double time =
        std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1).count();

    if (time < 0) time = 0;
    return time;
}


double geo_mean(const std::vector<double>& v) {
    double log_sum = 0;
    int count = 0;

    for (double x : v) {
        if (x > 0) {
            log_sum += std::log(x);
            count++;
        }
    }

    return std::exp(log_sum / count);
}

int main() {

    int n = 4000; 
    int m = 4000;

    std::vector<int> threads = {1,4,8,16};

    std::vector<double> A(n*n);
    std::vector<double> B(n*m);
    std::vector<double> Bcopy(n*m);

    
    for (int i = 0; i < n; i++)
        for (int j = 0; j <= i; j++)
            A[i*n + j] = (i == j) ? 2.0 : 1.0;

    for (int i = 0; i < n*m; i++)
        B[i] = 1.0;

    for (int t : threads) {

    std::cout << "\n==== Потоки: " << t << " ====\n";

    
    std::string env = "OPENBLAS_NUM_THREADS=" + std::to_string(t);
    putenv(const_cast<char*>(env.c_str()));

    std::vector<double> perf;

    for (int i = 0; i < 10; i++) {

        memcpy(Bcopy.data(), B.data(), sizeof(double)*n*m);
        double my = time_my(n, m, A.data(), Bcopy.data());

        memcpy(Bcopy.data(), B.data(), sizeof(double)*n*m);
        double blas = time_blas(n, m, A.data(), Bcopy.data());

        double p = (blas / my) * 100.0;
        perf.push_back(p);

        std::cout << "Запуск " << i+1
                  << " | Моё время: " << my << " сек"
                  << " | OpenBLAS: " << blas << " сек"
                  << " | Производительность: " << p << " %\n";
    }

    double gm = geo_mean(perf);

    std::cout << ">>> Среднее геометрическое: " << gm << " %\n";
}

    return 0;
}