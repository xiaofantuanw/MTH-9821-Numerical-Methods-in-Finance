#include "BinomialAvgTree.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

BinomialAvgTree::BinomialAvgTree(double S, double K, double T, double sigma, double r, double q, int n, bool call) : m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), N(n), Call(call)
{
	bt1 = new BinomialTree(S, K, T, sigma, r, q, n, call);
	bt2 = new BinomialTree(S, K, T, sigma, r, q, n+1, call);
}

double BinomialAvgTree::Price() const
{
	return (bt1->Price() + bt2->Price()) / 2;
}

double BinomialAvgTree::Delta() const
{
	return (bt1->Delta() + bt2->Delta()) / 2;
}

double BinomialAvgTree::Gamma() const
{
	return (bt1->Gamma() + bt2->Gamma()) / 2;
}

double BinomialAvgTree::Theta() const
{
	return (bt1->Theta() + bt2->Theta()) / 2;
}