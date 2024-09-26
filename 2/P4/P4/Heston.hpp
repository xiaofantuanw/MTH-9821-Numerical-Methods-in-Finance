#ifndef Heston_HPP
#define Heston_HPP
#include <vector>
#include "RandomGenerators.hpp"

class Heston
{
private:
	int m_m;
	int m_n;
	double m_lam;
	double m_Vtilde;
	double m_eta;
	double m_rho;
	double m_r;
	double m_V0;
	double m_S0;
	double m_K;
	double m_T;
	double opt_value;
	RandomGenerators rg;
public:
	Heston(int m,int n,double lam,double Vtilde,double eta,double rho,double r,double V0,double S0, double T,double K)
	{
		m_m = m;
		m_K = K;
		m_n = n;
		m_lam = lam;
		m_Vtilde = Vtilde;
		m_eta = eta;
		m_rho = rho;
		m_r = r;
		m_V0 = V0;
		m_S0 = S0;
		m_T = T;
		rg = RandomGenerators();
		opt_value = 0;
	}
	std::vector<double> alg()
	{
		std::vector <double> nm = rg.Gen_Normal(m_m * m_n);
		double dt = m_T / m_m;
		int k = 0;
		double Vhat = 0.;
		double Phat = 0.;
		double V, S;
		for (int i = 0; i < m_n; i++)
		{
			V = m_V0;
			S = m_S0;
			for (int j = 0; j < m_m; j++)
			{
				double z1 = nm[k];
				k++;
				double z2 = nm[k];
				k++;
				double V1 = 0;
				if (V > 0)
				{
					V1 = V;
				}
				S = S * exp((m_r - V1 / 2.) * dt+sqrt(V1 * dt) * z1);
				V = V1 - m_lam * (V1 - m_Vtilde) * dt + m_eta * sqrt(V1 * dt) * (m_rho * z1 + sqrt(1 - m_rho * m_rho) * z2);

			}
			Vhat += V / m_n;
			double S1 = 0;
			if (m_S0 - S > 0)
			{
				S1 = m_S0-S;
			}
			Phat += exp(-m_r * m_T) * (S1) / m_n;
		}
		std::vector<double> result ;
		result.push_back(Vhat);
		result.push_back(Phat);
		std::cout << "This is the wrong V (Variance): "<<Vhat << std::endl<<"This is the right V (Value of option): "<< Phat << std::endl;
		opt_value = Phat;
		return result;
	}

	//Implied Volatility given put or call
	double find_vol(bool Call_Put)
	{
		double value = opt_value;
		int maxit = 100;
		double precision = 1.0e-7;
		double sigma = 0.5;
		for (int i = 0; i < maxit; i++)
		{
			EuropeanOption tmp(0, m_S0, m_K, m_T, sigma, m_r, 0);
			double price = 0,vega=0;
			if (Call_Put == 1)
			{
				price = tmp.Call();
				vega = tmp.VegaCall();
			}
			else
			{
				price = tmp.Put();
				vega = tmp.VegaPut();
			}
			double diff = value - price;
			if (abs(diff) < precision)
			{
				std::cout << "Implied Volatility: " << sigma << std::endl;
				return sigma;
			}
				
			sigma += diff / vega;//Newton's Method
		}
		return 0;
	}
};

#endif // !Heston_HPP
