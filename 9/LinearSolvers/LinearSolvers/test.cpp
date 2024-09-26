#include "LinearSolvers.hpp"
#include <iostream>


int main() {
    LinearSolvers ls;
    mat A(3, 3);
    vec b(3);
    double tolerance = 1e-6;
    StoppingCriterion criterion = residual;

    A << 4, -1, 0,
        -1, 4, -1,
        0, -1, 3;
    b << 1, 2, 2;

    vec x0(3);
    x0 << 0, 0, 0;

    auto [x, iterations] = ls.jacobi_banded(A, b, x0, tolerance,1,1, criterion);

    std::cout << "The approximate solution is:\n" << x << std::endl;
    std::cout << "Number of iterations: " << iterations << std::endl;

    return 0;
}