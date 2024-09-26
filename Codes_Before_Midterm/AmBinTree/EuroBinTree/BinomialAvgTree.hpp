#ifndef BinomialAvgTree_HPP
#define BinomialAvgTree_HPP
#include "BinomialTree.hpp"

class BinomialAvgTree
{
private:
    // Characteristics of the option and the tree
    // t, S, K, T, sigma, r, q
    //double m_t;// Current time
    double m_S0;// Spot price
    double m_K;// Strike price
    double m_T;// Maturity
    double m_sigma;// Volatility
    double m_r;// Const interest rate
    double m_q;// Dividend rate
    bool Call;
    int N;// Number of time steps
    bool Eu; // European/American
    BinomialTree* bt1;
    BinomialTree* bt2;

public:
    BinomialAvgTree(double S, double K, double T, double sigma, double r, double q, int n, bool call, bool eu);
    virtual ~BinomialAvgTree() = default;
    // Black-Scholes Price
    double Price() const;
    //double Put() const;
    // Greeks
    double Delta() const;
    double Gamma() const;
    double Theta() const;
};


#endif // !BinomialAvgTree_HPP

