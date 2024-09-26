#include "FiniteDifference.hpp"
#include "EuropeanOption.hpp"
#include "LinearSolvers.hpp"
#include <cmath>
#include <iostream>
#include <vector>

void tridiagonalLUDecomposition(const std::vector<double>& a,
	const std::vector<double>& b,
	const std::vector<double>& c,
	std::vector<double>& l,
	std::vector<double>& u) {
	int n = a.size();
	l.resize(n, 0);
	u.resize(n, 0);

	u[0] = a[0];

	for (int i = 0; i < n - 1; i++) {
		l[i] = b[i] / u[i];
		u[i + 1] = a[i + 1] - l[i] * c[i];
	}
}

// Forward substitution
std::vector<double> forwardSubstitution(const std::vector<double>& l,
	const std::vector<double>& b) {
	int n = b.size();
	std::vector<double> y(n, 0.0);
	y[0] = b[0];

	for (int i = 1; i < n; i++) {
		y[i] = b[i] - l[i - 1] * y[i - 1];
	}

	return y;
}

// Backward substitution
std::vector<double> backwardSubstitution(const std::vector<double>& u,
	const std::vector<double>& c,
	const std::vector<double>& y) {
	int n = y.size();
	std::vector<double> x(n, 0.0);
	x[n - 1] = y[n - 1] / u[n - 1];

	for (int i = n - 2; i >= 0; i--) {
		x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
	}

	return x;
}

std::vector<double> solveTridiagonalSystem(const std::vector<double>& a,
	const std::vector<double>& b,
	const std::vector<double>& c,
	const std::vector<double>& l,
	const std::vector<double>& u,
	const std::vector<double>& rhs) {
	// Solve Ly = rhs using forward substitution
	std::vector<double> y = forwardSubstitution(l, rhs);

	// Solve Ux = y using backward substitution
	std::vector<double> x = backwardSubstitution(u, c, y);

	return x;
}


FiniteDifference::FiniteDifference(double S, double K, double B, double T, double sigma, double r, double q, int time_step, double alpha) : m_B(B), m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), M(time_step), alpha_tmp(alpha)
{
//Initialize the values
	x_compute = log(m_S0 / m_K);
	tau_final = T * sigma * sigma / 2.0;
	delta_tau = tau_final / M;
	delta_x_tmp = sqrt(delta_tau / alpha_tmp);
	x_left = log(B / K);
	N_left = (int)((x_compute - x_left) / delta_x_tmp);
	delta_x = (x_compute - x_left) / N_left;
	alpha_real = delta_tau / delta_x / delta_x;
	x_right = log(m_S0 / m_K) + (r - q - sigma * sigma / 2.) * T + 3. * sigma * sqrt(T);
	N_right = (int)((x_right - x_compute) / delta_x)+1;
	N = N_left + N_right;
	x_right = x_compute + N_right * delta_x;

	delta_x = (x_right - x_left) / N;
	
	a = (r - q) / sigma / sigma - 0.5;
	b = ((r - q) / sigma / sigma + 0.5) * ((r - q) / sigma / sigma + 0.5) + 2 * q / sigma / sigma;

	//x and tau grids
	for (int i = 0; i <= N; i++)
	{
		x_grid.push_back( x_left + i * delta_x);
	}

	for (int i = 0; i <= M; i++)
	{
		tau_grid.push_back( i * delta_tau);
	}

}

void FiniteDifference::print_intermediate_result() const
{
	std::cout << M << ", " << alpha_real << ", " << x_left << ", " << x_right << ", " << N << ", " << delta_x << ", " << delta_tau << std::endl;
	if (M == 4)
	{
		for (int i = 0; i <= M; i++)
		{
			for (int j = 0; j <= N; j++)
			{
				std::cout << Results[i][j] << ",";
			}
			std::cout << std::endl;
		}
	}
	
}

double FiniteDifference::f(double x)
{
	if (1 - exp(x) < 0)
	{
		return m_K * exp(a * x) * (-1 + exp(x));
	}
	return 0;
}

double FiniteDifference::g_left(double tau)
{
	return 0;
}

double FiniteDifference::g_right(double tau)
{
	return m_K * exp(a * x_right + b * tau) * (exp(x_right-2*m_q*tau/m_sigma/m_sigma)-exp(-2*m_r*tau/m_sigma/m_sigma) );
}

void FiniteDifference::Forward_Euler()
{
	//The first time level
	std::vector<double> this_level(N + 1);
	for (int i = 0; i <= N; i++)
	{
		this_level[i] = this->f(x_grid[i]);
	}
	Results.push_back(this_level);
	//The next time levels
	for (int t = 1; t <= M; t++)
	{
		std::vector<double> last_level(this_level);
		this_level[0] = g_left(tau_grid[t]);
		this_level[N] = g_right(tau_grid[t]);
		for (int i = 1; i < N; i++)
		{
			this_level[i] = alpha_real * last_level[i - 1] + alpha_real * last_level[i + 1] + (1 - 2 * alpha_real) * last_level[i];
		}
		Results.push_back(this_level);
	}
	
}

double FiniteDifference::uvalue() const
{
	return Results[M][N_left];
}

double FiniteDifference::Option_Value() const
{
	return exp(-a*x_compute-b*tau_final)*uvalue();
}

double FiniteDifference::Gamma() const
{
	double S_minus = m_K * exp(x_compute - delta_x);
	double S_plus = m_K * exp(x_compute + delta_x);
	double S_0 = m_K * exp(x_compute);
	double x_minus = x_compute - delta_x;
	double x_plus = x_compute + delta_x;
	double V_minus = exp(-a * x_minus - b * tau_final) * Results[M][N_left - 1];
	double V_plus = exp(-a * x_plus - b * tau_final) * Results[M][N_left + 1];
	double V_0=  exp(-a * x_compute - b * tau_final) * Results[M][N_left];
	return ((S_0 - S_minus) * V_plus - (S_plus - S_minus) * V_0 + (S_plus - S_0) * V_minus) / ((S_0 - S_minus) * (S_plus - S_0) * (S_plus - S_minus) / 2);
}

double FiniteDifference::Delta() const
{
	double S_minus = m_K * exp(x_compute - delta_x);
	double S_plus = m_K * exp(x_compute + delta_x);
	double S_0 = m_K * exp(x_compute);
	double x_minus = x_compute - delta_x;
	double x_plus = x_compute + delta_x;
	double V_minus = exp(-a * x_minus - b * tau_final) * Results[M][N_left - 1];
	double V_plus = exp(-a * x_plus - b * tau_final) * Results[M][N_left + 1];
	double V_0 = exp(-a * x_compute - b * tau_final) * Results[M][N_left];
	return (V_plus - V_minus) / (S_plus - S_minus);
}

double FiniteDifference::Theta() const
{
	double Option_minus = exp(-a * x_compute - b * (tau_final - delta_tau)) * Results[M - 1][N_left];
	return (Option_Value() - Option_minus) / (2 * delta_tau / m_sigma / m_sigma);
}

double FiniteDifference::pointwise_error() const
{
	EuropeanOption eu1(0, m_S0, m_K, m_T, m_sigma, m_r, m_q);
	EuropeanOption eu2(0, m_B*m_B/m_S0, m_K, m_T, m_sigma, m_r, m_q);
	double accu = eu1.Call() - eu2.Call() * pow(m_B / m_S0, 2 * a);
	return abs(accu - Option_Value());
}


void FiniteDifference::Backward_Euler()
{


	//The first time level
	std::vector<double> this_level(N + 1);
	for (int i = 0; i <= N; i++)
	{
		this_level[i] = this->f(x_grid[i]);
	}
	Results.push_back(this_level);
	//Do the LU decomposition
	std::vector<double> l, u;
	std::vector<double> a(N - 1, 1 + 2 * alpha_real);
	std::vector<double> b(N - 2, -alpha_real);
	std::vector<double> c(N - 2, -alpha_real);
	tridiagonalLUDecomposition(a, b, c, l, u);
	//The next time levels
	for (int t = 1; t <= M; t++)
	{
		std::vector<double> last_level(this_level);
		this_level[0] = g_left(tau_grid[t]);
		this_level[N] = g_right(tau_grid[t]);
		std::vector<double> rhs(N - 1);
		for (int i = 0; i < N - 1; i++)
		{
			rhs[i] = last_level[i + 1];
			if (i == 0)
				rhs[i] += alpha_real * this_level[0];
			if (i == N - 2)
				rhs[i] += alpha_real * this_level[N];
		}
		auto x = solveTridiagonalSystem(a, b, c, l, u, rhs);
		for (int i = 1; i < N; i++)
		{
			this_level[i] = x[i-1];
		}
		Results.push_back(this_level);
	}

}



void FiniteDifference::Crank_Nicolson()
{

	LinearSolvers ls;
	//The first time level
	std::vector<double> this_level(N + 1);
	for (int i = 0; i <= N; i++)
	{
		this_level[i] = this->f(x_grid[i]);
	}
	Results.push_back(this_level);
	//Do the LU decomposition
	//std::vector<double> l, u;
	//std::vector<double> a(N - 1, 1 +  alpha_real);
	//std::vector<double> b(N - 2, -0.5*alpha_real);
	//std::vector<double> c(N - 2, -0.5*alpha_real);
	//tridiagonalLUDecomposition(a, b, c, l, u);
	//The next time levels
	double tolerance = 1e-6;
	double omega = 1.2;
	mat A = mat::Zero(N-1, N-1);
	// Fill A and b as per the problem's specification
	for (int i = 0; i < N-1; ++i) {
		A(i, i) = 1+alpha_real; // Diagonal elements
		if (i > 0) {
			A(i, i - 1) = -0.5*alpha_real; // Sub-diagonal elements
		}
		if (i < N - 2) {
			A(i, i + 1) = -0.5*alpha_real; // Super-diagonal elements
		}
	}

	for (int t = 1; t <= M; t++)
	{
		std::vector<double> last_level(this_level);
		vec ll(N - 1);
		this_level[0] = g_left(tau_grid[t]);
		this_level[N] = g_right(tau_grid[t]);

		vec rhs(N - 1);
		for (int i = 0; i < N - 1; i++)
		{
			ll[i] = last_level[i];
			rhs[i] = last_level[i + 1] * (1 - alpha_real) + last_level[i] * alpha_real / 2 + last_level[i + 2] * alpha_real / 2;

			if (i == 0)
				rhs[i] += alpha_real/2 * this_level[0];
			if (i == N - 2)
				rhs[i] += alpha_real/2 * this_level[N];
		}
		//vec zz = vec::Zero(N - 1);
		auto res = ls.sor(A, rhs, ll, tolerance, omega, consecutive);
		auto x = std::get<0>(res);
		for (int i = 1; i < N; i++)
		{
			this_level[i] = x[i - 1];
		}
		Results.push_back(this_level);
	}

}