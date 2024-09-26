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
    double m_tdiv;//Time of paying dividend
    double m_sigma;// Volatility
    double m_r;// Const interest rate
    double m_q;// Dividend rate
    //bool Call;//Call option or not
    int M1;// Number of time steps before dividend payment
    int M2;//Number of time steps after dividend payment
    double alpha_tmp;//Alpha

    //Internal Variables
    double tau_final;
    double tau_div;
    double x_left_new;
    double x_right_new;
    double x_left;
    double x_right;
    double delta_tau_1;
    double delta_tau_2_tmp;
    double delta_tau_2;
    int N;
    int Nleft;
    int Nright;
    double delta_x;
    double alpha_real;
    double a;
    double b;
    double x_compute;
    double x_compute_new;
    int index;
    double xi;
    double xi_plus;
    double Si;
    double Si_plus;
    std::vector<std::vector<double>> Results;
    std::vector<double> x_grid;
    std::vector<double> tau_grid;
    std::vector<double> tau_grid_new;

public:
    FiniteDifference(double S, double K, double T, double tdiv, double sigma, double r, double q, int time_step, double alpha);

    //Boundary conditions
    double f(double x);
    double g_left(double tau);
    double g_right(double tau);
    double g_right_new(double tau);
    void print_intermediate_result() const;
    void Forward_Euler();
    //void Backward_Euler();
    void Crank_Nicolson();
    //Two ways to find out the value of the option
    double u_value() const;
    double Option_value() const;
    double Delta() const;
    double Gamma() const;
    double Theta() const;
    double RMS_error() const;
};



#endif // !FiniteDifference_HPP
