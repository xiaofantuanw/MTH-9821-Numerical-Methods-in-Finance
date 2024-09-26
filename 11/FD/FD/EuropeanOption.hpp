#ifndef EuropeanOption_hpp
#define EuropeanOption_hpp

#include <functional>

class EuropeanOption {
private:
    // t, S, K, T, sigma, r, q
    double m_t;// Current time
    double m_S0;// Spot price
    double m_K;// Strike price
    double m_T;// Maturity
    double m_sigma;// Volatility
    double m_r;// Const interest rate
    double m_q;// Dividend rate

    //Values needed for pricing
    double d1;
    double d2;
    double Zd1;
    double Zd2;
    double Nd1;
    double Nd2;
    double q_disc; // Discount of dividend
    double r_disc; // Discount of interest rate

    double PDF_Normal(double t) const;   // PDF of a standard normal variable
    double CDF_Normal(double t) const; // CDF of a standard normal variable


public:
    EuropeanOption(double t, double S, double K, double T, double sigma, double r, double q);
    ~EuropeanOption() = default;

    std::function<double(double, double)> CallPayoff() const;
    std::function<double(double, double)> PutPayoff() const;

    // Black-Scholes Price
    double Call() const;
    double Put() const;

    // Greeks
    double DeltaCall() const;
    double DeltaPut() const;
    double GammaCall() const;
    double GammaPut() const;
    double ThetaCall() const;
    double ThetaPut() const;
    double VegaCall() const;
    double VegaPut() const;
    double RhoCall() const;
    double RhoPut() const;
};

#endif /* EuropeanOption_hpp */