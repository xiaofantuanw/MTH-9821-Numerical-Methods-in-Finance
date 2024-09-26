#include "EuropeanOption.hpp"
#include "RandomGenerators.hpp"
#include "Heston.hpp"




int main()
{
	for (int k = 0; k < 6; k++)
	{
		Heston h(175, 500 * ((int)pow(2, k)),4,0.1225,0.25,-0.15,0.05,0.09,50.0,0.5,50);
		h.alg();
		h.find_vol(0);
	}
}