#include "BinomialBSTreeRichardson.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

BinomialBSTreeRichardson::BinomialBSTreeRichardson(double S, double K, double T, double sigma, double r, double q, int n, bool call) : m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), N(n), Call(call)
{
	bt1 = new BinomialBSTree(S, K, T, sigma, r, q, n, call);
	bt2 = new BinomialBSTree(S, K, T, sigma, r, q, n/2, call);
}

double BinomialBSTreeRichardson::Price() const
{
	return (2*bt1->Price() - bt2->Price());
}

double BinomialBSTreeRichardson::Delta() const
{
	return (2*bt1->Delta() - bt2->Delta()) ;
}

double BinomialBSTreeRichardson::Gamma() const
{
	return (2*bt1->Gamma() - bt2->Gamma()) ;
}

double BinomialBSTreeRichardson::Theta() const
{
	return (2*bt1->Theta() - bt2->Theta());
}