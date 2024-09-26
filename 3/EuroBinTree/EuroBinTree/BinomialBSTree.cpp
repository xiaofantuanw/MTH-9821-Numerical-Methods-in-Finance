#include "BinomialBSTree.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>
#include "EuropeanOption.hpp"

double BinomialBSTree::Pricing(double S) const
{
	EuropeanOption eu(0, S, m_K, dt, m_sigma, m_r, m_q);
	if (Call == 1)
	{
		return eu.Call();
	}
	else
	{
		return eu.Put();
	}
}

BinomialBSTree::BinomialBSTree(double S, double K, double T, double sigma, double r, double q, int n, bool call) : m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), N(n-1), Call(call)
{
	dt = m_T / n;
	u = exp(m_sigma * sqrt(dt));
	//std::cout << u << std::endl;
	d = 1 / u;
	disc = exp(-m_r * dt);
	p = (exp((m_r - m_q) * dt) - d) / (u - d);
	this->Simulate_Tree();
}

void BinomialBSTree::Simulate_Tree()
{
	//double dt = m_T / N;
	//double u = exp(m_sigma * sqrt(dt));
	//std::cout << u << std::endl;
	//double d = 1 / u;
	//double disc = exp(-m_r * dt);
	//double p = (exp((m_r - m_q) * dt) - d) / (u - d);
	//std::cout << p << std::endl;
	std::vector<double> init(N+1, 0);
	//std::cout << init.size() << std::endl;
	
	for (int i = 0; i <= N; i++)
	{
		init[i] = Pricing(m_S0 * pow(u, N - i) * pow(d, i));
	}
	Nodes.push_back(init);
	for (int j = (N - 1); j >= 0; j--)
	{
		std::vector<double> tmp(j + 1, 0);
		for (int i = 0; i <= j; i++)
		{
			tmp[i] = disc * (p * Nodes[N - j - 1][i] + (1 - p) * Nodes[N - j - 1][i + 1]);
		}
		Nodes.push_back(tmp);
	}
}

double BinomialBSTree::Price() const
{
	return Nodes[N][0];
}

double BinomialBSTree::Delta() const
{
	return (Nodes[N - 1][0] - Nodes[N - 1][1]) / (u * m_S0 - d * m_S0);
}

double BinomialBSTree::Gamma() const
{
	double delta1 = (Nodes[N - 2][0] - Nodes[N - 2][1]) / (u * u * m_S0 - u * d * m_S0);
	double delta2 = (Nodes[N - 2][1] - Nodes[N - 2][2]) / (u * d * m_S0 - d * d * m_S0);
	return (delta1 - delta2) * 2 / (u * u * m_S0 - d * d * m_S0);
}

double BinomialBSTree::Theta() const
{
	return (Nodes[N - 2][1] - Nodes[N][0]) / 2 / dt;
}