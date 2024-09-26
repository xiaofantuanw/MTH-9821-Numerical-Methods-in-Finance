#ifndef TrinomialBSTree_HPP
#define TrinomialBSTree_HPP
#include <vector>

class TrinomialBSTree
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
    bool Eu; // European/American
    std::vector<std::vector<double>> Nodes;//Prices at the nodes of the tree


    //internal variables
    double u;
    double d;
    double disc;
    double pu;
    double pd;
    double pm;
    double dt;
	double Pricing(double S) const;//Price of a call option given the stock price
	//double PutPricing(double S) const;//Price of a put option given the stock price
public:
    TrinomialBSTree(double S, double K, double T, double sigma, double r, double q, int N, bool call, bool eu);
    virtual ~TrinomialBSTree() = default;

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

#endif // !TrinomialBSTree_HPP
