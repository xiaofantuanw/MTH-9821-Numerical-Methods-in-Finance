#ifndef BinomialBSTreeRichardson_HPP
#define BinomialBSTreeRichardson_HPP
#include "BinomialBSTree.hpp"

class BinomialBSTreeRichardson
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
    BinomialBSTree* bt1;
    BinomialBSTree* bt2;

public:
    BinomialBSTreeRichardson(double S, double K, double T, double sigma, double r, double q, int n, bool call);
    virtual ~BinomialBSTreeRichardson() = default;
    // Black-Scholes Price
    double Price() const;
    //double Put() const;
    // Greeks
    double Delta() const;
    double Gamma() const;
    double Theta() const;
};


#endif // !BinomialAvgTree_HPP


