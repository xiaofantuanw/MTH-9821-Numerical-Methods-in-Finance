#include "FiniteDifference.hpp"
#include <cmath>
#include <iostream>

FiniteDifference::FiniteDifference(double S, double K, double T, double sigma, double r, double q, int time_step, double alpha) : m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), M(time_step), alpha_tmp(alpha)
{
//Initialize the values
	tau_final = T * sigma * sigma / 2.0;
	delta_tau = tau_final / M;
	x_left = log(m_S0 / m_K) + (r - q - sigma *sigma / 2.) * T - 3. * sigma * sqrt(T);
	x_right = log(m_S0 / m_K) + (r - q - sigma *sigma / 2.) * T + 3. * sigma * sqrt(T);
	N = (int)((x_right - x_left) / sqrt(delta_tau / alpha_tmp));
	delta_x = (x_right - x_left) / N;
	alpha_real = delta_tau / delta_x / delta_x;
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

	//Compute values
	x_compute = log(m_S0 / m_K);
	for (int i = 0; i < N; i++)
	{
		if (x_grid[i + 1] > x_compute)
		{
			index = i;
			xi = x_grid[i];
			xi_plus = x_grid[i + 1];
			Si = m_K * exp(x_grid[i]);
			Si_plus = m_K * exp(x_grid[i + 1]);
			break;
		}
	}
}

void FiniteDifference::print_intermediate_result() const
{
	//std::cout << alpha_tmp << ',' << alpha_real <<','<<N<< std::endl;
	for (int i = 0; i <= M; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			std::cout << Results[i][j] << ",";
		}
		std::cout << std::endl;
	}
}

double FiniteDifference::f(double x)
{
	if (1 - exp(x) > 0)
	{
		return m_K * exp(a * x) * (1 - exp(x));
	}
	return 0;
}

double FiniteDifference::g_left(double tau)
{
	return m_K * exp(a * x_left + b * tau) * (exp(-2 * m_r * tau / m_sigma / m_sigma) - exp(x_left - 2 * m_q * tau / m_sigma / m_sigma));
}

double FiniteDifference::g_right(double tau)
{
	return 0;
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

double FiniteDifference::Value1() const
{
	double Vi = exp(-a * xi - b * tau_final) * Results[M][index];
	double Vi_plus = exp(-a * xi_plus - b * tau_final) * Results[M][index + 1];
	double V = ((Si_plus - m_S0) * Vi + (m_S0 - Si) * Vi_plus) / (Si_plus - Si);
	return V;
}

double FiniteDifference::Value2() const
{
	double u_compute = ((xi_plus - x_compute) * Results[M][index] + (x_compute - xi) * Results[M][index + 1]) / (xi_plus - xi);
	double V = exp(-a * x_compute - b * tau_final) * u_compute;
	return V;
}

double FiniteDifference::Gamma() const
{
	int i = index;
	double Si_minus = m_K * exp(x_grid[i - 1]);
	double Si_plus2 = m_K * exp(x_grid[i + 2]);
	double Vi = exp(-a * xi - b * tau_final) * Results[M][i];
	double Vi_plus = exp(-a * xi_plus - b * tau_final) * Results[M][i + 1];
	double Vi_minus = exp(-a * x_grid[i - 1] - b * tau_final) * Results[M][i - 1];
	double Vi_plus2 = exp(-a * x_grid[i + 2] - b * tau_final) * Results[M][i + 2];
	double gamma = (((Vi_plus2 - Vi_plus) / (Si_plus2 - Si_plus)) - ((Vi - Vi_minus) / (Si - Si_minus))) / (Si_plus2 + Si_plus - Si - Si_minus) * 2;
	return gamma;
}

double FiniteDifference::Delta() const
{
	int i = index;
	double Vi = exp(-a * xi - b * tau_final) * Results[M][i];
	double Vi_plus = exp(-a * xi_plus - b * tau_final) * Results[M][i + 1];
	return (Vi_plus - Vi) / (Si_plus - Si);
}

double FiniteDifference::Theta() const
{
	double delta_t = 2 * delta_tau / m_sigma / m_sigma;
	double Vi_last = exp(-a * xi - b * (tau_final - delta_tau)) * Results[M - 1][index];
	double Vi_plus_last = exp(-a * xi_plus - b * (tau_final - delta_tau)) * Results[M - 1][index+1];
	double V_last = ((Si_plus - m_S0) * Vi_last + (m_S0 - Si) * Vi_plus_last) / (Si_plus - Si);
	double theta = (V_last - this->Value1()) / delta_t;
	return theta;
}

double FiniteDifference::RMS_error() const
{
	double Vk, Vaccu,err=0,count=0;
	for (int i = 0; i <= N; i++)
	{
		Vk=(exp(-a * x_grid[i] - b * tau_final) * Results[M][i]);
		Vaccu=(EuropeanOption(0, m_K*exp(x_left+i*delta_x), m_K, m_T, m_sigma, m_r, m_q).Put());
		if (Vaccu / m_K / exp(x_left) > 0.00001)
		{
			err += (Vk - Vaccu) * (Vk - Vaccu) / Vaccu / Vaccu;
			count++;
		}
	}
	return sqrt(err / count);

	
}