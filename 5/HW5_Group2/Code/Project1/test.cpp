#include "american_pricers.h"
#include "mc_regression.h"
#include "regression.h"
#include <random>
#include <iomanip>
#include <iostream>
#include <tuple>

using namespace montecarlo;
using namespace montecarlo::regression;
using namespace linalg::regression;
// Initialize the LCG with std::minstd_rand
std::minstd_rand lcg,lcg2;
std::normal_distribution<double> normal_dist(0.0, 1.0); // mean = 0, standard deviation = 1



int main() {
    // Parameters
    std::cout << std::fixed << std::setprecision(6);
    double spot = 40.5;
    double strike = 44;
    double r = 0.04;
    double q = 0.0;
    double sigma = 0.2;
    double maturity = 0.5;
    int M = 6;
    int N = 10000;
    int K = 100;
    double dt = maturity / M;
    // Generate paths
    //Eigen::ArrayXXd paths = generate_paths(spot, r, q, sigma, maturity, M, N);
    //generate paths

    for (int degree = 2; degree <= 10; degree++)
    {
        lcg.seed(100);  // Seed the LCG
        lcg2.seed(200);
        double avg_result = 0, se = 0;
        double avg_result_fw = 0, se_fw = 0;
        double avg_result_fw_out = 0, se_fw_out = 0;
        arr results(K), results_fw(K), results_fw_out(K);
        for (int k = 0; k < K; k++)
        {
            arr2 paths(N, M);
            paths.col(0).fill(spot);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M - 1; j++) {
                    std::normal_distribution<double> z(0.0, 1.0);
                    paths(i, j + 1) = paths(i, j) * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * z(lcg));
                }
            }
            arr2 paths_out(N, M);
            paths_out.col(0).fill(spot);
            for (int i = 0; i < N; i++) {
                for (int j = 0; j < M - 1; j++) {
                    std::normal_distribution<double> z(0.0, 1.0);
                    paths_out(i, j + 1) = paths_out(i, j) * exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * z(lcg2));
                }
            }
            // ... rest of the experiment and option pricing ...
            HermiteMCRegression p1(M, degree, true, true);
            //PolynomialMCRegression p1(M, degree, true);
            double result = regression_pricer_backward(spot, paths, -1, strike, maturity, r, p1, Longstaff_Schwartz);
            //std::cout << result << std::endl;
            avg_result += result / K;
            results(k) = result;
            double result_fw = regression_pricer_forward(spot, paths, -1, strike, maturity, r, p1);
            avg_result_fw += result_fw / K;
            results_fw(k) = result_fw;
            double result_fw_out = regression_pricer_forward(spot, paths_out, -1, strike, maturity, r, p1);
            //std::cout << result_fw_out << std::endl;
            avg_result_fw_out += result_fw_out / K;
            results_fw_out(k) = result_fw_out;
        }
        for (int k = 0; k < K; k++)
        {
            se += (results(k) - avg_result)*(results(k) - avg_result)/K;
            se_fw_out += (results_fw_out(k) - avg_result_fw_out) * (results_fw_out(k) - avg_result_fw_out) / K;
            se_fw += (results_fw(k) - avg_result_fw) * (results_fw(k) - avg_result_fw) / K;
        }
        se = sqrt(se) / sqrt(K);
        se_fw = sqrt(se_fw) / sqrt(K);
        se_fw_out = sqrt(se_fw_out) / sqrt(K);

        std::cout << "Longstaff, Poly, degree= " <<degree<<": "<< avg_result<<"; Standard Error: "<<se<<std::endl;
        std::cout << "Longstaff, Poly, forward, degree= " << degree << ": " << avg_result_fw << "; Standard Error: " << se_fw << std::endl;
        std::cout << "Longstaff, Poly, forward, out, degree= " << degree << ": " << avg_result_fw_out << "; Standard Error: " << se_fw_out << std::endl<<std::endl;
    }
}