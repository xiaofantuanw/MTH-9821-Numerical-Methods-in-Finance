#include "LinearSolvers.hpp"
#include <iostream>
#include <iomanip>


int main() {
    LinearSolvers ls;
    int n = 14;
    double omega = 1.15;
    mat A = mat::Zero(n, n);
    vec b(n);
    vec x0 = vec::Ones(n);
    double tolerance = 1e-6;
    std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6) ;
    // Fill A and b as per the problem's specification
    for (int i = 0; i < n; ++i) {
        A(i, i) = 2; // Diagonal elements
        if (i > 0) {
            A(i, i - 1) = -1; // Sub-diagonal elements
        }
        if (i < n - 1) {
            A(i, i + 1) = -1; // Super-diagonal elements
        }
        b(i) = i * i; // b vector
    }

    // Solve using Jacobi iteration with residual-based stopping criterion
    auto [x_res_jacobi, iterations_res_jacobi] = ls.jacobi_banded(A, b, x0, tolerance, 1, 1, residual);
    std::cout << "Jacobi solution with residual-based criterion:\n" << x_res_jacobi << "\nIterations: " << iterations_res_jacobi << std::endl;

    // Solve using Jacobi iteration with consecutive approximation stopping criterion
    auto [x_con_jacobi, iterations_con_jacobi] = ls.jacobi_banded(A, b, x0, tolerance, 1, 1, consecutive);
    std::cout << "Jacobi solution with consecutive approximation criterion:\n" << x_con_jacobi << "\nIterations: " << iterations_con_jacobi << std::endl;

    // Solve using Gauss-Seidel iteration with residual-based stopping criterion
    auto [x_res_gs, iterations_res_gs] = ls.gauss_seidel(A, b, x0, tolerance, residual);
    std::cout << "Gauss-Seidel solution with residual-based criterion:\n" << x_res_gs << "\nIterations: " << iterations_res_gs << std::endl;

    // Solve using Gauss-Seidel iteration with consecutive approximation stopping criterion
    auto [x_con_gs, iterations_con_gs] = ls.gauss_seidel(A, b, x0, tolerance, consecutive);
    std::cout << "Gauss-Seidel solution with consecutive approximation criterion:\n" << x_con_gs << "\nIterations: " << iterations_con_gs << std::endl;

    // Solve using SOR iteration with residual-based stopping criterion
    auto [x_res_sor, iterations_res_sor] = ls.sor(A, b, x0, tolerance, omega, residual);
    std::cout << "SOR solution with residual-based criterion:\n" << x_res_sor << "\nIterations: " << iterations_res_sor << std::endl;

    // Solve using SOR iteration with consecutive approximation stopping criterion
    auto [x_con_sor, iterations_con_sor] = ls.sor(A, b, x0, tolerance, omega, consecutive);
    std::cout << "SOR solution with consecutive approximation criterion:\n" << x_con_sor << "\nIterations: " << iterations_con_sor << std::endl;

    // Loop over the desired w values (from 1.02 to 1.98, step 0.02)
    for (omega = 1.02; omega <= 2.00; omega += 0.02) {
        // Solve using SOR iteration with residual-based stopping criterion
        auto [_, iterations_sor] = ls.sor(A, b, x0, tolerance, omega, residual);
        std::cout << "w = " << omega
            << ", Iterations to convergence: " << iterations_sor << std::endl;
        //std::cout << iterations_sor << std::endl;
    }

    return 0;
}