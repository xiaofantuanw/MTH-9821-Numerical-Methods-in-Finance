#include "BinomialAvgTree.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>
#include <vector>

BinomialAvgTree::BinomialAvgTree(double S, double K, double T, double sigma, double r, double q, int n, bool call, bool eu) : m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q), N(n), Call(call), Eu(eu)
{
	bt1 = new BinomialTree(S, K, T, sigma, r, q, n, call, eu);
	bt2 = new BinomialTree(S, K, T, sigma, r, q, n+1, call, eu);
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