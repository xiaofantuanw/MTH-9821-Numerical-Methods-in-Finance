#include "BinomialTree.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

BinomialTree::BinomialTree(double S, double K, double T, double sigma, double r, double q, int n, bool call, bool eu) :  m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), N(n), Call(call), Eu(eu)
{
	dt = m_T / N;
	u = exp(m_sigma * sqrt(dt));
	//std::cout << u << std::endl;
	d = 1 / u;
	disc = exp(-m_r * dt);
	p = (exp((m_r - m_q) * dt) - d) / (u - d);
	this->Simulate_Tree();
}

double BinomialTree::CallPricing(double S) const
{
	return std::max(S-m_K, 0.);
}

double BinomialTree::PutPricing(double S) const
{
	return std::max(m_K-S, 0.);
}

void BinomialTree::Simulate_Tree()
{
	//double dt = m_T / N;
	//double u = exp(m_sigma * sqrt(dt));
	//std::cout << u << std::endl;
	//double d = 1 / u;
	//double disc = exp(-m_r * dt);
	//double p = (exp((m_r - m_q) * dt) - d) / (u - d);
	//std::cout << p << std::endl;
	std::vector<double> init(N + 1, 0);
	//std::cout << init.size() << std::endl;
	if (Call == 1)
	{
		for (int i = 0; i <= N; i++)
		{
			init[i] = CallPricing(m_S0 * pow(u, N - i) * pow(d, i));
		}
	}
	else
	{
		for (int i = 0; i <= N; i++)
		{
			init[i] = PutPricing(m_S0 * pow(u, N - i) * pow(d, i));
		}
	}
	Nodes.push_back(init);
	for (int j = (N - 1); j >= 0; j--)
	{
		std::vector<double> tmp(j + 1, 0);
		double curr_price = m_S0 * pow(u, j);  // Initialize curr_price for i=0
		for (int i = 0; i <= j; i++)
		{
			tmp[i] = disc * (p * Nodes[N - j - 1][i] + (1 - p) * Nodes[N - j - 1][i + 1]);

			if (!Eu)
			{
				double exercise_value = (Call ? (curr_price - m_K) : (m_K - curr_price));
				exercise_value = std::max(exercise_value, 0.0);
				tmp[i] = std::max(tmp[i], exercise_value);

				// Update curr_price for next iteration
				curr_price *= (d / u);
			}
		}
		Nodes.push_back(tmp);
	}
}

double BinomialTree::Price() const
{
	return Nodes[N][0];
}

double BinomialTree::Delta() const
{
	return (Nodes[N - 1][0] - Nodes[N - 1][1]) / (u * m_S0 - d * m_S0);
}

double BinomialTree::Gamma() const
{
	double delta1 = (Nodes[N - 2][0] - Nodes[N - 2][1]) / (u * u * m_S0 - u * d * m_S0);
	double delta2= (Nodes[N - 2][1] - Nodes[N - 2][2]) / (u * d * m_S0 -d * d * m_S0);
	return (delta1 - delta2) * 2 / (u * u * m_S0 - d * d * m_S0);
}

double BinomialTree::Theta() const
{
	return (Nodes[N - 2][1] - Nodes[N][0]) /2 / dt;
}