#include "LinearSolvers.hpp"
#include <cmath>

vec LinearSolvers::forward_subst(const mat& L, const vec& b)
{
    int n = L.rows();
    vec x(n);

    for (int i = 0; i < n; ++i) {
        x(i) = b(i);
        for (int j = 0; j < i; ++j) {
            x(i) -= L(i, j) * x(j);
        }
        x(i) /= L(i, i);
    }

    return x;
}




vec LinearSolvers::backward_subst(const mat& U, const vec& b)
{
    int n = U.rows();
    vec x(n);

    for (int i = n - 1; i >= 0; --i) {
        x(i) = b(i);
        for (int j = i + 1; j < n; ++j) {
            x(i) -= U(i, j) * x(j);
        }
        x(i) /= U(i, i);
    }

    return x;
}


std::tuple<mat, mat > LinearSolvers::lu_no_pivoting(mat A)
{
    int n = A.rows();
    mat L = mat::Identity(n, n);
    mat U = A;

    for (int k = 0; k < n; ++k) {
        for (int i = k + 1; i < n; ++i) {
            L(i, k) = U(i, k) / U(k, k);
            for (int j = k; j < n; ++j) {
                U(i, j) -= L(i, k) * U(k, j);
            }
        }
    }

    return { L, U };
}


std::tuple<permutation, mat, mat > LinearSolvers::lu_row_pivoting(mat A)
{
    int n = A.rows();
    permutation P(n);  // Initialize permutation matrix of size n
    P.setIdentity();
    mat L = mat::Identity(n, n);
    mat U = A;

    for (int k = 0; k < n; ++k) {
        // Find the row with the maximum pivot
        int maxRow = k;
        double maxPivot = std::abs(U(k, k));
        for (int i = k + 1; i < n; ++i) {
            double pivot = std::abs(U(i, k));
            if (pivot > maxPivot) {
                maxPivot = pivot;
                maxRow = i;
            }
        }

        // Check if we need to apply a row interchange
        if (k != maxRow) {
            // Swap the rows in U matrix
            U.row(k).swap(U.row(maxRow));

            // Swap the rows in permutation matrix
            P.applyTranspositionOnTheRight(k, maxRow);

            // Swap the rows in L matrix
            if (k > 0) {
                L.block(k, 0, 1, k).swap(L.block(maxRow, 0, 1, k));
            }
        }

        // Perform the elimination as usual
        for (int i = k + 1; i < n; ++i) {
            L(i, k) = U(i, k) / U(k, k);
            for (int j = k + 1; j < n; ++j) {
                U(i, j) -= L(i, k) * U(k, j);
            }
        }
    }

    // Reset the lower part of U to zero
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j) {
            U(i, j) = 0.0;
        }
    }

    return std::make_tuple(P, L, U);
}

mat LinearSolvers::cholesky(mat A)
{
    int n = A.rows();
    mat U = mat::Zero(n, n);

    for (int i = 0; i < n; ++i) {
        for (int j = i; j < n; ++j) {
            double sum = 0.0;
            if (i == j) { // Diagonal elements
                for (int k = 0; k < i; ++k) {
                    sum += std::pow(U(k, i), 2);
                }
                double val = A(i, i) - sum;
                if (val <= 0) {
                    throw std::runtime_error("Matrix is not positive definite.");
                }
                U(i, i) = std::sqrt(val);
            }
            else { // Off-diagonal elements
                for (int k = 0; k < i; ++k) {
                    sum += U(k, i) * U(k, j);
                }
                U(i, j) = (A(i, j) - sum) / U(i, i);
            }
        }
    }

    return U;
}


std::tuple<vec, uint> LinearSolvers::gauss_seidel(const mat& A, const vec& b, vec x,
                                   const double tolerance,
                                   const StoppingCriterion criterion) {
    uint iterations = 0;
    int n = A.rows();
    vec x_old = x;

    while (true) {
        x_old = x;
        
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A(i, j) * x(j);
                }
            }
            x(i) = (b(i) - sum) / A(i, i);
        }

        ++iterations;

        // Stopping criterion check
        if (criterion == consecutive) {
            if ((x - x_old).norm() < tolerance) {
                break;
            }
        } else {
            if ((b - A * x).norm() < tolerance) {
                break;
            }
        }
    }

    return std::make_tuple(x, iterations);
}

std::tuple<vec, uint> LinearSolvers::jacobi(const mat& A, const vec& b, vec x,
    const double tolerance,
    const StoppingCriterion criterion) {
    uint iterations = 0;
    int n = A.rows();
    vec x_old = x; // Previous iteration's solution

    while (true) {
        x_old = x;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) {
                    sum += A(i, j) * x_old(j);
                }
            }
            x(i) = (b(i) - sum) / A(i, i);
        }

        ++iterations;

        // Stopping criterion check
        if (criterion == consecutive) {
            if ((x - x_old).norm() < tolerance) {
                break;
            }
        }
        else {
            if ((b - A * x).norm() < tolerance) {
                break;
            }
        }

        // Optional: Add a maximum iteration count to prevent infinite loops
        // if the method does not converge.
    }

    return std::make_tuple(x, iterations);
}


std::tuple<vec, uint> LinearSolvers::sor(const mat& A, const vec& b, vec x,
    const double tolerance, const double omega,
    const StoppingCriterion criterion) {
    if (omega <= 0.0 || omega >= 2.0) {
        throw std::invalid_argument("Omega must be between 0 and 2 for SOR to converge.");
    }

    uint iterations = 0;
    int n = A.rows();
    vec x_old = x;

    while (true) {
        x_old = x;

        for (int i = 0; i < n; ++i) {
            double sigma = 0.0;
            for (int j = 0; j < i; ++j) {
                sigma += A(i, j) * x(j);
            }
            for (int j = i + 1; j < n; ++j) {
                sigma += A(i, j) * x_old(j);
            }
            x(i) = (1 - omega) * x_old(i) + (omega / A(i, i)) * (b(i) - sigma);
        }

        ++iterations;

        // Stopping criterion check
        if (criterion == consecutive) {
            if ((x - x_old).norm() < tolerance) {
                break;
            }
        }
        else {
            if ((b - A * x).norm() < tolerance) {
                break;
            }
        }

        // Optional: Add a maximum iteration count to prevent infinite loops
        // if the method does not converge.
    }

    return std::make_tuple(x, iterations);
}

std::tuple<vec, uint> LinearSolvers::gauss_seidel_banded(const mat& A, const vec& b, vec x,
    const double tolerance, const int lower_bandwidth, const int upper_bandwidth,
    const StoppingCriterion criterion) {
    uint iterations = 0;
    int n = A.rows();
    vec x_old = x;

    while (true) {
        x_old = x;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            int j_min = std::max(i - lower_bandwidth, 0);
            int j_max = std::min(i + upper_bandwidth, n - 1);
            for (int j = j_min; j <= j_max; ++j) {
                if (j != i) {
                    sum += A(i, j) * x(j);
                }
            }
            x(i) = (b(i) - sum) / A(i, i);
        }

        ++iterations;

        // Stopping criterion check
        if (criterion == consecutive) {
            if ((x - x_old).norm() < tolerance) {
                break;
            }
        }
        else {
            if ((b - A * x).norm() < tolerance) {
                break;
            }
        }
    }

    return std::make_tuple(x, iterations);
}

std::tuple<vec, uint> LinearSolvers::jacobi_banded(const mat& A, const vec& b, vec x,
    const double tolerance, const int lower_bandwidth, const int upper_bandwidth,
    const StoppingCriterion criterion) {
    uint iterations = 0;
    int n = A.rows();
    vec x_old = x;

    while (true) {
        x_old = x;

        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            int j_min = std::max(i - lower_bandwidth, 0);
            int j_max = std::min(i + upper_bandwidth, n - 1);
            for (int j = j_min; j <= j_max; ++j) {
                if (j != i) {
                    sum += A(i, j) * x_old(j);
                }
            }
            x(i) = (b(i) - sum) / A(i, i);
        }

        ++iterations;

        // Stopping criterion check
        if (criterion == consecutive) {
            if ((x - x_old).norm() < tolerance) {
                break;
            }
        }
        else {
            if ((b - A * x).norm() < tolerance) {
                break;
            }
        }
    }

    return std::make_tuple(x, iterations);
}

std::tuple<vec, uint> LinearSolvers::sor_banded(const vec& early_exercise, const mat& A, const vec& b, vec x,
    const double tolerance, const double omega,
    const int lower_bandwidth, const int upper_bandwidth,
    const StoppingCriterion criterion) {
    if (omega <= 0.0 || omega >= 2.0) {
        throw std::invalid_argument("Omega must be between 0 and 2 for SOR to converge.");
    }

    uint iterations = 0;
    int n = A.rows();
    vec x_old = x; // Previous iteration's solution

    while (true) {
        x_old = x;

        for (int i = 0; i < n; ++i) {
            double sigma = 0.0;
            int j_min = std::max(i - lower_bandwidth, 0);
            int j_max = std::min(i + upper_bandwidth, n - 1);
            for (int j = j_min; j < i; ++j) {
                sigma += A(i, j) * x(j);
            }
            for (int j = i + 1; j <= j_max; ++j) {
                sigma += A(i, j) * x_old(j);
            }
            x(i) = (1 - omega) * x_old(i) + (omega / A(i, i)) * (b(i) - sigma);
            if (x(i) < early_exercise(i))
                x(i) = early_exercise(i);
        }
        

        ++iterations;

        // Stopping criterion check
        if (criterion == consecutive) {
            if ((x - x_old).norm() < tolerance) {
                break;
            }
        }
        else {
            if ((b - A * x).norm() < tolerance) {
                break;
            }
        }

        // Optional: Add a maximum iteration count to prevent infinite loops
        // if the method does not converge.
    }

    return std::make_tuple(x, iterations);
}
