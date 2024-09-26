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
	double S = 41, K=39, T=1.0, sigma=0.25, r=0.03, q=0.005;
	EuropeanOption eu1(0, S, K, T, sigma, r, q);
	//std::cout << "Hello world." << std::endl;
	std::cout << "Data of Black-Scholes Put" << std::endl;
	double V_BS = eu1.Put();
	double Delta_BS = eu1.DeltaPut();
	double Gamma_BS = eu1.GammaPut();
	double Theta_BS = eu1.ThetaPut();
	std::cout <<"V_BS: "<< V_BS << "; ";
	std::cout << "Delta_BS: " << Delta_BS << "; ";
	std::cout << "Gamma_BS: " << Gamma_BS << "; ";
	std::cout << "Theta_BS: " << Theta_BS << "\n";
	//BinomialTree b1(S, K, T, sigma, r, q, 1280,0);
	int N = 10;
	for (int k = 0; k <= 7; k++)
	{
		BinomialTree b1(S, K, T, sigma, r, q, N , 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_BS);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_BS) << ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_BS) << ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_BS);
		std::cout << "\n" ;
		N *= 2;
	}

	//N = 10;
	//for (int k = 0; k <= 7; k++)
	//{
	//	BinomialAvgTree b1(S, K, T, sigma, r, q, N, 0);
	//	std::cout << N << ',';
	//	double V_N = b1.Price();
	//	std::cout << V_N << ',';
	//	double err = abs(V_N - V_BS);
	//	std::cout << err << ',' << err * N << ',' << err * N * N << ',';
	//	double Delta_N = b1.Delta();
	//	std::cout << Delta_N << ',' << abs(Delta_N - Delta_BS) << ',';
	//	double Gamma_N = b1.Gamma();
	//	std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_BS) << ',';
	//	double Theta_N = b1.Theta();
	//	std::cout << Theta_N << ',' << abs(Theta_N - Theta_BS);
	//	std::cout << "\n";
	//	N *= 2;
	//}

	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		BinomialBSTree b1(S, K, T, sigma, r, q, N, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_BS);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_BS) << ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_BS) << ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_BS);
		std::cout << "\n";
		N *= 2;
	}

	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		BinomialBSTreeRichardson b1(S, K, T, sigma, r, q, N, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_BS);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_BS) << ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_BS) << ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_BS);
		std::cout << "\n";
		N *= 2;
	}
}