#ifndef BinomialBSTree_HPP
#define BinomialBSTree_HPP
#include <vector>

class BinomialBSTree
{
private:
    //Characteristics of the option and the tree
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
    std::vector<std::vector<double>> Nodes;//Prices at the nodes of the tree


    //internal variables
    double u;
    double d;
    double disc;
    double pu;
    double pm;
    double pd;
    double dt;
	double Pricing(double S) const;//Price of a call option given the stock price
	//double PutPricing(double S) const;//Price of a put option given the stock price
public:
    BinomialBSTree(double S, double K, double T, double sigma, double r, double q, int N, bool call);
    virtual ~BinomialBSTree() = default;

    //Simulate the Call/Put tree
    void Simulate_Tree();
    //void Simulate_Put_Tree();

    // Black-Scholes Price
    double Price() const;
    //double Put() const;
    // Greeks
    double Delta() const;
    double Gamma() const;
    double Theta() const;
};

#endif // !BinomialBSTree_HPP
