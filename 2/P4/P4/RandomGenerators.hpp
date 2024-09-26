#ifndef RandomGenerators_HPP
#define RandomGenerators_HPP
#include <vector>
#include <iostream>

class RandomGenerators
{
private:
	long long seed;
	//int n;

public:
	RandomGenerators()
	{
		seed = 1;
		//n = 100;
	}
	RandomGenerators(int s, int num)
	{
		seed = s;
		//n = num;
	}
	std::vector<double> Gen_Uniform(int n)
	{
		std::vector<double> result;
		long long tmp = seed;
		for (int i = 0; i < n; i++)
		{
			tmp = (tmp * 39373) % ((long long)(pow(2, 31)) - 1);
			double f = (double)tmp / (pow(2, 31) - 1);
			result.push_back(f);
			//std::cout << f << std::endl;
		}
		return result;
	}

	std::vector<double> Gen_Normal(int n)
	{
		std::vector<double> result;
		std::vector<double> u = this->Gen_Uniform(10*n);
		int k = 0;
		for (int i = 0; i < n; i++)
		{
			double x = 10.,u1,u2;
			while (x > 1.)
			{
				u1 = u[k];
				u1 = 2 * u1 - 1;
				k++;
				u2 = u[k];
				u2 = 2 * u2 - 1;
				k++;
				x = u1 * u1 + u2 * u2;
			}
			double y = sqrt(-2. * log(x) / x);
			double z1 = u1 * y;
			double z2 = u2 * y;
			result.push_back(z1);
			result.push_back(z2);
			//std::cout << z1 << std::endl << z2 << std::endl;
		}
		return result;
	}
};

#endif // !RandomGenerators_HPP
