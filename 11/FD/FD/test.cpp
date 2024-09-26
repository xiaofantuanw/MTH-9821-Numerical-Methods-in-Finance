#include "EuropeanOption.hpp"
#include "FiniteDifference.hpp"
#include <iostream>
#include <iomanip>

int main()
{
	std::cout << std::fixed << std::setprecision(6);
	double sigma = 0.25;
	double S = 53;
	double q = 0.01;
	double T = 1.0;
	double tdiv = 5.0 / 12.0;
	double K = 50;
	double r = 0.04;
	int M = 4;
	double alpha = 4;
	
	EuropeanOption bs(0, S, K, T, sigma, r, q);
	//double accu_value = bs.Put();
	std::cout << bs.Call() << std::endl;
	std::cout << bs.DeltaCall() << std::endl << bs.GammaCall() << std::endl << bs.ThetaCall()<<std::endl;
	
	
	for (int i = 1; i <= 4; i++)
	{
		FiniteDifference fd(S, K, T, tdiv, sigma, r, q, M, alpha);
		fd.Crank_Nicolson();
		fd.print_intermediate_result();
		std::cout << "u_value: " << fd.u_value()<<std::endl;
		std::cout << "Option value: " << fd.Option_value()<<std::endl;
		std::cout << "Delta: " << fd.Delta() << std::endl;
		std::cout << "Gamma: " << fd.Gamma() << std::endl;
		std::cout << "Theta: " << fd.Theta() << std::endl;
		M *= 4;
	}
		
	

}