#include "TrinomialBSTreeRichardson.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

TrinomialBSTreeRichardson::TrinomialBSTreeRichardson(double S, double K, double T, double sigma, double r, double q, int n, bool call, bool eu) : m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), N(n), Call(call), Eu(eu)
{
	bt1 = new TrinomialBSTree(S, K, T, sigma, r, q, n, call, eu);
	bt2 = new TrinomialBSTree(S, K, T, sigma, r, q, n/2, call, eu);
}

double TrinomialBSTreeRichardson::Price() const
{
	return (2*bt1->Price() - bt2->Price());
}

double TrinomialBSTreeRichardson::Delta() const
{
	return (2*bt1->Delta() - bt2->Delta()) ;
}

double TrinomialBSTreeRichardson::Gamma() const
{
	return (2*bt1->Gamma() - bt2->Gamma()) ;
}

double TrinomialBSTreeRichardson::Theta() const
{
	return (2*bt1->Theta() - bt2->Theta());
}