#include "EuropeanOption.hpp"
#include "FiniteDifference.hpp"
#include <iostream>
#include <iomanip>

int main()
{
	std::cout << std::fixed << std::setprecision(6);

	//double sigma = 0.28;
	//double S = 37;
	//double q = 0.02;
	//double T = 0.75;
	//double K = 40;
	//double r = 0.04;
	double S = 40.5, K = 44, T = 0.5, sigma = 0.2, r = 0.04, q = 0.0;
	int M = 4;
	double alpha = 5;
	
	double accu_value = 0;
	//std::cout << bs.DeltaPut()<<','<< bs.GammaPut() <<','<< bs.ThetaPut() << std::endl;
	
	for (int i = 1; i <= 4; i++)
	{
		FiniteDifference fd(S, K, T, sigma, r, q, M, alpha,0);
		fd.Crank_Nicolson();
		double fd_put_value = fd.Value1();
		//double fd_put_value_2 = fd.Value2();
		//double error_pointwise = abs(fd_put_value - accu_value);
		//double error_pointwise_2 = abs(fd_put_value_2 - accu_value);
		//double err_RMS = fd.RMS_error();
		double delta = fd.Delta();
		double gamma = fd.Gamma();
		double theta = fd.Theta();
		//if (M == 4)
		//{
		//	fd.print_intermediate_result();
		//}
		std::cout << M << ',' << fd_put_value <<','<<delta << ',' << gamma << ',' << theta << std::endl;
		M *= 4;
	}

}