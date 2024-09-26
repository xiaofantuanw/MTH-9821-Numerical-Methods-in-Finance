#include <iostream>
#include <iomanip>
#include "EuropeanOption.hpp"
#include "TrinomialTree.hpp"
#include "TrinomialAvgTree.hpp"
#include "TrinomialBSTreeRichardson.hpp"
#include "TrinomialBSTree.hpp"

int main()
{
	std::cout << std::setiosflags(std::ios::fixed) << std::setprecision(6);
	//double S = 54, K=50, T=1.0, sigma=0.29, r=0.0375, q=0.01;
	double S = 33, K = 32, T = 10./12., sigma = 0.24, r = 0.045, q = 0.02;

	TrinomialAvgTree BatExact(S, K, T, sigma, r, q, 10000, 0, 0);

	std::cout << "Data of Black-Scholes Put" << std::endl;
	double V_exact = BatExact.Price();
	double Delta_exact = BatExact.Delta();
	double Gamma_exact = BatExact.Gamma();
	double Theta_exact = BatExact.Theta();
	std::cout <<"V_exact: "<< V_exact << "; ";
	std::cout << "Delta_exact: " << Delta_exact << "; ";
	std::cout << "Gamma_exact: " << Gamma_exact << "; ";
	std::cout << "Theta_exact: " << Theta_exact << "\n";
	//TrinomialTree b1(S, K, T, sigma, r, q, 1280,0);
	int N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialTree b1(S, K, T, sigma, r, q, N , 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n" ;
		N *= 2;
	}

	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialAvgTree b1(S, K, T, sigma, r, q, N, 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n";
		N *= 2;
	}

	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialBSTree b1(S, K, T, sigma, r, q, N, 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n";
		N *= 2;
	}
	
	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialBSTreeRichardson b1(S, K, T, sigma, r, q, N, 0, 0);
		std::cout << N << ',';
		double V_N = b1.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n";
		N *= 2;
	}


	// Variance Reduction
	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialTree b1(S, K, T, sigma, r, q, N, 0, 0);
		TrinomialTree b2(S, K, T, sigma, r, q, N, 0, 1);
		EuropeanOption eu1(0, S, K, T, sigma, r, q);

		std::cout << N << ',';
		double V_N = b1.Price() + eu1.Put() - b2.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta() + eu1.DeltaPut() - b2.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma() + eu1.GammaPut() - b2.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta() + eu1.ThetaPut() - b2.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n";
		N *= 2;
	}

	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialAvgTree b1(S, K, T, sigma, r, q, N, 0, 0);
		TrinomialAvgTree b2(S, K, T, sigma, r, q, N, 0, 1);
		EuropeanOption eu1(0, S, K, T, sigma, r, q);

		std::cout << N << ',';
		double V_N = b1.Price() + eu1.Put() - b2.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta() + eu1.DeltaPut() - b2.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma() + eu1.GammaPut() - b2.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta() + eu1.ThetaPut() - b2.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n";
		N *= 2;
	}

	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialBSTree b1(S, K, T, sigma, r, q, N, 0, 0);
		TrinomialBSTree b2(S, K, T, sigma, r, q, N, 0, 1);
		EuropeanOption eu1(0, S, K, T, sigma, r, q);

		std::cout << N << ',';
		double V_N = b1.Price() + eu1.Put() - b2.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta() + eu1.DeltaPut() - b2.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma() + eu1.GammaPut() - b2.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta() + eu1.ThetaPut() - b2.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n";
		N *= 2;
	}

	N = 10;
	for (int k = 0; k <= 7; k++)
	{
		TrinomialBSTreeRichardson b1(S, K, T, sigma, r, q, N, 0, 0);
		TrinomialBSTreeRichardson b2(S, K, T, sigma, r, q, N, 0, 1);
		EuropeanOption eu1(0, S, K, T, sigma, r, q);

		std::cout << N << ',';
		double V_N = b1.Price() + eu1.Put() - b2.Price();
		std::cout << V_N << ',';
		double err = abs(V_N - V_exact);
		std::cout << err << ',' << err * N << ',' << err * N * N << ',';
		double Delta_N = b1.Delta() + eu1.DeltaPut() - b2.Delta();
		std::cout << Delta_N << ',' << abs(Delta_N - Delta_exact) << ',';
		double Gamma_N = b1.Gamma() + eu1.GammaPut() - b2.Gamma();
		std::cout << Gamma_N << ',' << abs(Gamma_N - Gamma_exact) << ',';
		double Theta_N = b1.Theta() + eu1.ThetaPut() - b2.Theta();
		std::cout << Theta_N << ',' << abs(Theta_N - Theta_exact);
		std::cout << "\n";
		N *= 2;
	}
}