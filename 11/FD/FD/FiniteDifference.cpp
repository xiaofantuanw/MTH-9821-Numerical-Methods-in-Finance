#include "FiniteDifference.hpp"
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


FiniteDifference::FiniteDifference(double S, double K, double T,double tdiv, double sigma, double r, double q, int time_step, double alpha) : m_tdiv(tdiv), m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), M1(time_step), alpha_tmp(alpha)
{
//Initialize the values
	x_compute = log(S / K) + log(1 - q);
	x_compute_new = log(S / K);
	tau_final = T * sigma * sigma / 2.0;
	tau_div = (T - tdiv) * sigma * sigma / 2.0;
	delta_tau_1 = tau_div / M1;
	double x_left_tmp = log(m_S0 / m_K) + (r - sigma *sigma / 2.) * T - 3. * sigma * sqrt(T);
	double x_right_tmp = log(m_S0 / m_K) + (r - sigma *sigma / 2.) * T + 3. * sigma * sqrt(T);
	delta_x = sqrt(delta_tau_1 / alpha_tmp);
	Nleft = (int)((x_compute - x_left_tmp) / delta_x) + 1;
	Nright = (int)((x_right_tmp-x_compute) / delta_x) + 1;
	//std::cout << Nleft << ',' << Nright << std::endl;
	N = Nleft+Nright;
	x_left = x_compute - Nleft * delta_x;
	x_right = x_compute + Nright * delta_x;
	a = (r) / sigma / sigma - 0.5;
	b = ((r) / sigma / sigma + 0.5) * ((r) / sigma / sigma + 0.5);

	//x and tau grids
	for (int i = 0; i <= N; i++)
	{
		x_grid.push_back( x_left + i * delta_x);
	}

	for (int i = 0; i <= M1; i++)
	{
		tau_grid.push_back( i * delta_tau_1);
	}
	
	x_left_new = x_left - log(1 - q);
	x_right_new = x_right - log(1 - q);
	delta_tau_2_tmp = alpha_tmp * delta_x * delta_x;
	M2 = (int)((tau_final - tau_div) / delta_tau_2_tmp) + 1;
	delta_tau_2 = (tau_final - tau_div) / M2;
	alpha_real = delta_tau_2 / delta_x / delta_x;
	//std::cout << M2 << std::endl;
	for (int i = 0; i <= M2; i++)
	{
		tau_grid_new.push_back(i * delta_tau_2 + tau_div);
	}
	//Compute values
	//x_compute = log(m_S0 / m_K);
	//for (int i = 0; i < N; i++)
	//{
		//if (x_grid[i + 1] > x_compute)
		//{
			//index = i;
			//xi = x_grid[i];
			//xi_plus = x_grid[i + 1];
			//Si = m_K * exp(x_grid[i]);
			//Si_plus = m_K * exp(x_grid[i + 1]);
			//break;
		//}
	//}
}

void FiniteDifference::print_intermediate_result() const
{
	std::cout << "M_1: " << M1 << std::endl;
	std::cout << "M_2: " << M2<<std::endl;
	std::cout << "alpha_2: " << alpha_real<<std::endl;
	std::cout << "N: " << N<<std::endl;
	std::cout << "x_left: " << x_left << std::endl;
	std::cout << "x_right: " << x_right << std::endl;
	std::cout << "x_left,new: " << x_left_new << std::endl;
	std::cout << "x_right,new: " << x_right_new << std::endl;
	std::cout << "tau_div: " << tau_div << std::endl;
	std::cout << "delta_tau_1: " << delta_tau_1 << std::endl;
	std::cout << "delta_tau_2: " << delta_tau_2 << std::endl;
	std::cout << "delta_x: " << delta_x << std::endl;
	//std::cout << alpha_tmp << ',' << alpha_real <<','<<N<< std::endl;
	//if (M1 == 4)
	//{
	//	for (int i = 0; i < size(Results); i++)
	//	{
	//		for (int j = 0; j <= N; j++)
	//		{
	//			std::cout << Results[i][j] << ",";
	//		}
	//		std::cout << std::endl;
	//	}
	//}
}

double FiniteDifference::f(double x)
{
	if (exp(x) - 1 > 0)
	{
		return m_K * exp(a * x) * (exp(x)-1);
	}
	return 0;
}

double FiniteDifference::g_left(double tau)
{
	return 0;
}

double FiniteDifference::g_right(double tau)
{
	return m_K*exp(a*x_right+b*tau)*(exp(x_right)-exp(-2*m_r*tau/m_sigma/m_sigma));
}

double FiniteDifference::g_right_new(double tau)
{
	return m_K * exp(a * (x_right-log(1-m_q)) + b * tau) * (exp(x_right) - exp(-2 * m_r * tau / m_sigma / m_sigma));
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
	for (int t = 1; t <= M1; t++)
	{
		std::vector<double> last_level(this_level);
		this_level[0] = g_left(tau_grid[t]);
		this_level[N] = g_right(tau_grid[t]);
		for (int i = 1; i < N; i++)
		{
			this_level[i] = alpha_tmp * last_level[i - 1] + alpha_tmp * last_level[i + 1] + (1 - 2 * alpha_tmp) * last_level[i];
		}
		Results.push_back(this_level);
	}

	//The dividend payment level
	std::vector<double> last_level(this_level);
	for (int i = 0; i < N + 1; i++)
	{
		this_level[i] = pow(1 - m_q, -a) * last_level[i];
	}
	Results.push_back(this_level);

	//The time levels before payment
	for (int t = 1; t <= M2; t++)
	{
		std::vector<double> last_level(this_level);
		this_level[0] = g_left(tau_grid_new[t]);
		this_level[N] = g_right_new(tau_grid_new[t]);
		for (int i = 1; i < N; i++)
		{
			this_level[i]=alpha_real * last_level[i - 1] + alpha_real * last_level[i + 1] + (1 - 2 * alpha_real) * last_level[i];
		}
		Results.push_back(this_level);
	}
	
}

double FiniteDifference::u_value() const
{
	return Results[M1 + M2 + 1][Nleft];
}

double FiniteDifference::Option_value() const
{
	double u = this->u_value();
	return exp(-a * x_compute_new - b * tau_final) * u;
	
}

double FiniteDifference::Gamma() const
{
	double S_minus = m_K * exp(x_compute_new - delta_x);
	double S_0 = m_K * exp(x_compute_new);
	double S_plus = m_K * exp(x_compute_new + delta_x);
	double V_minus = exp(-a * (x_compute_new - delta_x) - b * tau_final) * Results[M1 + M2 + 1][Nleft - 1];
	double V_0 = this->Option_value();
	double V_plus = exp(-a * (x_compute_new + delta_x) - b * tau_final) * Results[M1 + M2 + 1][Nleft + 1];
	double gamma = ((S_0 - S_minus) * V_plus - (S_plus - S_minus) * V_0 + (S_plus - S_0) * V_minus) / ((S_0 - S_minus) * (S_plus - S_0) * (S_plus - S_minus) / 2);
	return gamma;
}

double FiniteDifference::Delta() const
{
	double S_minus = m_K * exp(x_compute_new - delta_x);
	double S_0 = m_K * exp(x_compute_new);
	double S_plus = m_K * exp(x_compute_new + delta_x);
	double V_minus = exp(-a * (x_compute_new - delta_x) - b * tau_final) * Results[M1 + M2 + 1][Nleft - 1];
	double V_0 = this->Option_value();
	double V_plus= exp(-a * (x_compute_new + delta_x) - b * tau_final) * Results[M1 + M2 + 1][Nleft + 1];
	return (V_plus - V_minus) / (S_plus - S_minus);
}

double FiniteDifference::Theta() const
{
	double delta_t = 2 * delta_tau_2 / m_sigma / m_sigma;
	double V_plus = exp(-a * x_compute_new - b * (tau_final - delta_tau_2)) * Results[M1 + M2][Nleft];
	return (V_plus - this->Option_value()) / delta_t;
}

double FiniteDifference::RMS_error() const
{
	double Vk, Vaccu,err=0,count=0;
	for (int i = 0; i <= N; i++)
	{
		Vk=(exp(-a * x_grid[i] - b * tau_final) * Results[M1][i]);
		Vaccu=(EuropeanOption(0, m_K*exp(x_left+i*delta_x), m_K, m_T, m_sigma, m_r, m_q).Put());
		if (Vaccu / m_K / exp(x_left) > 0.00001)
		{
			err += (Vk - Vaccu) * (Vk - Vaccu) / Vaccu / Vaccu;
			count++;
		}
	}
	return sqrt(err / count);

	
}


void FiniteDifference::Crank_Nicolson()
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
	std::vector<double> d(N - 1, 1 +  alpha_tmp);
	std::vector<double> b(N - 2, -0.5*alpha_tmp);
	std::vector<double> c(N - 2, -0.5*alpha_tmp);
	tridiagonalLUDecomposition(d, b, c, l, u);
	//The next time levels
	for (int t = 1; t <= M1; t++)
	{
		std::vector<double> last_level(this_level);
		this_level[0] = g_left(tau_grid[t]);
		this_level[N] = g_right(tau_grid[t]);
		std::vector<double> rhs(N - 1);
		for (int i = 0; i < N - 1; i++)
		{
			rhs[i] = last_level[i + 1] * (1 - alpha_tmp) + last_level[i] * alpha_tmp / 2 + last_level[i + 2] * alpha_tmp / 2;

			if (i == 0)
				rhs[i] += alpha_tmp/2 * this_level[0];
			if (i == N - 2)
				rhs[i] += alpha_tmp/2 * this_level[N];
		}
		auto x = solveTridiagonalSystem(d, b, c, l, u, rhs);
		for (int i = 1; i < N; i++)
		{
			this_level[i] = x[i - 1];
		}
		Results.push_back(this_level);
	}
	
	//The dividend payment level
	std::vector<double> last_level(this_level);
	for (int i = 0; i < N + 1; i++)
	{
		this_level[i] = pow(1 - m_q, -a) * last_level[i];
	}
	Results.push_back(this_level);

	//The next levels
	//Do the LU decomposition for a different alpha
	std::vector<double> d0(N - 1, 1 + alpha_real);
	std::vector<double> b0(N - 2, -0.5 * alpha_real);
	std::vector<double> c0(N - 2, -0.5 * alpha_real);
	tridiagonalLUDecomposition(d0, b0, c0, l, u);
	//The next time levels
	for (int t = 1; t <= M2; t++)
	{
		std::vector<double> last_level(this_level);
		this_level[0] = g_left(tau_grid_new[t]);
		this_level[N] = g_right_new(tau_grid_new[t]);
		std::vector<double> rhs(N - 1);
		for (int i = 0; i < N - 1; i++)
		{
			rhs[i] = last_level[i + 1] * (1 - alpha_real) + last_level[i] * alpha_real / 2 + last_level[i + 2] * alpha_real / 2;

			if (i == 0)
				rhs[i] += alpha_real / 2 * this_level[0];
			if (i == N - 2)
				rhs[i] += alpha_real / 2 * this_level[N];
		}
		auto x = solveTridiagonalSystem(d0, b0, c0, l, u, rhs);
		for (int i = 1; i < N; i++)
		{
			this_level[i] = x[i - 1];
		}
		Results.push_back(this_level);
	}



}