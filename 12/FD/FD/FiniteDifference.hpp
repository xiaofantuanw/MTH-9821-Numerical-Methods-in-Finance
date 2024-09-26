#ifndef FiniteDifference_HPP
#define FiniteDifference_HPP
#include <vector>
#include "EuropeanOption.hpp"

class FiniteDifference
{
private:
    //Characteristics of the option (It is a put option)
    // t, S, K, T, sigma, r, q
    //double m_t;// Current time
    double m_S0;// Spot price
    double m_K;// Strike price
    double m_T;// Maturity
    double m_sigma;// Volatility
    double m_r;// Const interest rate
    double m_q;// Dividend rate
    //bool Call;//Call option or not
    int M;// Number of time steps
    double alpha_tmp;//Alpha
    double m_B;//Barrier

    //Internal Variables
    double tau_final;
    double x_left;
    double x_right;
    double delta_tau;
    int N;
    double delta_x;
    double alpha_real;
    double a;
    double b;
    double x_compute;
    std::vector<std::vector<double>> Results;
    std::vector<double> x_grid;
    std::vector<double> tau_grid;
    int N_left;
    int N_right;
    double delta_x_tmp;

public:
    FiniteDifference(double S, double K, double B, double T, double sigma, double r, double q, int time_step, double alpha);

    //Boundary conditions
    double f(double x);
    double g_left(double tau);
    double g_right(double tau);
    void print_intermediate_result() const;
    void Forward_Euler();
    void Backward_Euler();
    void Crank_Nicolson();
    //Two ways to find out the value of the option
    double uvalue() const;
    double Option_Value() const;
    double Delta() const;
    double Gamma() const;
    double Theta() const;
    double pointwise_error() const;
};



#endif // !FiniteDifference_HPP
