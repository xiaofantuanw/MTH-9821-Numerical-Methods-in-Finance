#include "EuropeanOption.hpp"
#include "FiniteDifference.hpp"
#include <iostream>
#include <iomanip>

int main()
{
	std::cout << std::fixed << std::setprecision(6);
	double sigma = 0.28;
	double S = 42;
	double q = 0.015;
	double T = 7. / 12.;
	double K = 40;
	double B = 36;
	double r = 0.04;
	int M = 4;
	double alpha = 4;
	
	EuropeanOption bs(0, S, K, T, sigma, r, q);
	double accu_value = bs.Put();
	std::cout << bs.DeltaPut()<<','<< bs.GammaPut() <<','<< bs.ThetaPut() << std::endl;
	
	for (int i = 1; i <= 4; i++)
	{
		FiniteDifference fd(S, K, B, T, sigma, r, q, M, alpha);
		fd.Crank_Nicolson();
		//fd.print_intermediate_result();
		std::cout << fd.uvalue() <<", " <<fd.Option_Value()<<", "<<fd.pointwise_error()<<", "<<fd.Delta()<<", "<<fd.Gamma()<<", "<<fd.Theta()<< std::endl;
		//double fd_put_value = fd.Value1();
		//double fd_put_value_2 = fd.Value2();
		//double error_pointwise = abs(fd_put_value - accu_value);
		//double error_pointwise_2 = abs(fd_put_value_2 - accu_value);
		//double err_RMS = fd.RMS_error();
		//double delta = fd.Delta();
		//double gamma = fd.Gamma();
		//double theta = fd.Theta();
		//if (M == 4)
		//{
		//	fd.print_intermediate_result();
		//}
		//std::cout << M << ',' << error_pointwise << ',' << error_pointwise_2 <<','<<err_RMS<<','<<delta << ',' << gamma << ',' << theta << std::endl;
		M *= 4;
	}

}