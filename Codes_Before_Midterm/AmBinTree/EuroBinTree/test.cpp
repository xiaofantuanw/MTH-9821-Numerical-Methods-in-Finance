#include <iostream>
#include <iomanip>
#include "EuropeanOption.hpp"
#include "BinomialTree.hpp"
#include "BinomialAvgTree.hpp"
#include "BinomialBSTreeRichardson.hpp"
#include "BinomialBSTree.hpp"

int main()
{
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6);
	double S = 40.5, K=44, T=0.5, sigma=0.2, r=0.04, q=0.0;

	//BinomialAvgTree BatExact(S, K, T, sigma, r, q, 10000, 0, 0);

	//std::cout << "Data of Black-Scholes Put" << std::endl;
	//double V_exact = BatExact.Price();
	//double Delta_exact = BatExact.Delta();
	//double Gamma_exact = BatExact.Gamma();
	//double Theta_exact = BatExact.Theta();
	//std::cout <<"V_exact: "<< V_exact << "; ";
	//std::cout << "Delta_exact: " << Delta_exact << "; ";
	//std::cout << "Gamma_exact: " << Gamma_exact << "; ";
	//std::cout << "Theta_exact: " << Theta_exact << "\n";
	//BinomialTree b1(S, K, T, sigma, r, q, 1280,0);
	int N = 10;
	std::cout << "Binomial Tree" << '\n';
	for (int k = 0; k <= 7; k++)
	{
		
		BinomialTree b1(S, K, T, sigma, r, q, N , 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		//double err = abs(V_N - V_exact);
		//std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N <<  ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N <<  ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << '\n';
		N *= 2;
	}

	N = 10;
	std::cout << '\n'<< "Binomial Average Tree" << '\n';
	for (int k = 0; k <= 7; k++)
	{
		
		BinomialAvgTree b1(S, K, T, sigma, r, q, N, 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		//double err = abs(V_N - V_exact);
		//std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' ;
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' ;
		double Theta_N = b1.Theta();
		std::cout << Theta_N ;
		std::cout << "\n";
		N *= 2;
	}

	N = 10;
	std::cout << '\n' << "Binomial Black-Scholes Tree" << '\n';
	for (int k = 0; k <= 7; k++)
	{
		
		BinomialBSTree b1(S, K, T, sigma, r, q, N, 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' ;
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' ;
		double Theta_N = b1.Theta();
		std::cout << Theta_N;
		std::cout << "\n";
		N *= 2;
	}

	N = 10;
	std::cout << '\n' << "Binomial Black-Scholes Tree with Richardson" << '\n';
	for (int k = 0; k <= 7; k++)
	{
		
		BinomialBSTreeRichardson b1(S, K, T, sigma, r, q, N, 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' ;
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' ;
		double Theta_N = b1.Theta();
		std::cout << Theta_N;
		std::cout << "\n";
		N *= 2;
	}
}