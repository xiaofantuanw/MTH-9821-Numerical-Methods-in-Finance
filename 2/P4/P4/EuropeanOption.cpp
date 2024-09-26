#define _USE_MATH_DEFINES
#include "EuropeanOption.hpp"
#include <cmath>
#include <iostream>
#include <algorithm>

double EuropeanOption::CDF_Normal(double t) const 
{
    return std::erfc(-t / std::sqrt(2.)) / 2.;
}

double EuropeanOption::PDF_Normal(double t) const 
{
    return std::exp(-t * t / 2.) / std::sqrt(2. * M_PI);
}

EuropeanOption::EuropeanOption(double t, double S, double K, double T, double sigma, double r, double q) : m_t(t), m_S0(S), m_K(K), m_T(T), m_sigma(sigma), m_r(r), m_q(q) 
{
    d1 = (log(S / K) + (r - q + sigma * sigma / 2.) * (T - t)) / (sigma * std::sqrt(T - t));
    d2 = d1 - sigma * std::sqrt(T - t);

    Zd1 = this->PDF_Normal(d1);
    Zd2 = this->PDF_Normal(d2);
    Nd1 = this->CDF_Normal(d1);
    Nd2 = this->CDF_Normal(d2);

    q_disc = std::exp(-q * (T - t));
    r_disc = std::exp(-r * (T - t));
}

double EuropeanOption::Call() const 
{
    return m_S0 * q_disc * Nd1 - m_K * r_disc * Nd2;
}

double EuropeanOption::Put() const 
{
    return -m_S0 * q_disc * (1 - Nd1) + m_K * r_disc * (1 - Nd2);
}


double EuropeanOption::DeltaCall() const 
{
    return q_disc * Nd1;
}
double EuropeanOption::DeltaPut() const 
{
    return -q_disc * (1- Nd1);
}

double EuropeanOption::GammaCall() const 
{
    return q_disc / (m_S0 * m_sigma * std::sqrt(m_T - m_t)) * Zd1;
}
double EuropeanOption::GammaPut() const 
{
    return this->GammaCall();
}

double EuropeanOption::ThetaCall() const 
{
    double theta = -(m_S0 * m_sigma * q_disc) / (2 * std::sqrt(m_T - m_t)) * Zd1+ m_q * m_S0 * q_disc * Nd1 - m_r * m_K * r_disc * Nd2;

    return theta;
}

double EuropeanOption::ThetaPut() const 
{
    double theta = -(m_S0 * m_sigma * q_disc) / (2 * std::sqrt(m_T - m_t)) * Zd1-m_q * m_S0 * q_disc * (1 - Nd1)+m_r * m_K * r_disc * (1 - Nd2);

    return theta;
}

double EuropeanOption::VegaCall() const 
{
    return m_S0 * q_disc * std::sqrt(m_T - m_t) * Zd1;
}

double EuropeanOption::VegaPut() const 
{
    return this->VegaCall();
}

double EuropeanOption::RhoCall() const 
{
    return m_K * (m_T - m_t) * r_disc * Nd2;
}

double EuropeanOption::RhoPut() const 
{
    return -m_K * (m_T - m_t) * r_disc * (1 - Nd2);
}
