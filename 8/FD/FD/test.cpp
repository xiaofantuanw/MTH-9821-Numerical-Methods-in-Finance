#include "EuropeanOption.hpp"
#include "FiniteDifference.hpp"
#include <iostream>
#include <iomanip>

int main()
{
	std::cout << std::fixed << std::setprecision(6);
	double sigma = 0.32;
	double S = 36;
	double q = 0.02;
	double T = 0.75;
	double K = 40;
	double r = 0.045;
	int M = 4;
	double alpha = 5;
	
	EuropeanOption bs(0, S, K, T, sigma, r, q);
	double accu_value = bs.Put();
	std::cout << bs.DeltaPut()<<','<< bs.GammaPut() <<','<< bs.ThetaPut() << std::endl;
	
	for (int i = 1; i <= 6; i++)
	{
		FiniteDifference fd(S, K, T, sigma, r, q, M, alpha);
		fd.Crank_Nicolson();
		double fd_put_value = fd.Value1();
		double fd_put_value_2 = fd.Value2();
		double error_pointwise = abs(fd_put_value - accu_value);
		double error_pointwise_2 = abs(fd_put_value_2 - accu_value);
		double err_RMS = fd.RMS_error();
		double delta = fd.Delta();
		double gamma = fd.Gamma();
		double theta = fd.Theta();
		if (M == 4)
		{
			fd.print_intermediate_result();
		}
		std::cout << M << ',' << error_pointwise << ',' << error_pointwise_2 <<','<<err_RMS<<','<<delta << ',' << gamma << ',' << theta << std::endl;
		M *= 4;
	}

}