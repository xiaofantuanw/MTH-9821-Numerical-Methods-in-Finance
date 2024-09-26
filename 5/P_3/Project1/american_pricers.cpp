//Tengxiao Fan
//This is backward and forward regression methods
#include "american_pricers.h"
#include <iostream>

namespace montecarlo
{
	double regression_pricer_backward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
		const double interest_rate, regression::MCRegression& mc_regression, const MonteCarloRegressionMethod& method)
	{
        int N = paths.rows();
        int M = paths.cols();
        double dt = maturity / M;

        arr2 E = (w * (paths - strike)).max(arr2::Zero(N, M));
        arr2 C = arr2::Zero(N, M);
        arr2 P = arr2::Zero(N, M);
        P.col(M - 1) = E.col(M - 1);

        for (int t = M - 2; t > 0; --t) {
            arr St = paths.col(t);
            arr Pt_delta = P.col(t + 1) * exp(-interest_rate * dt);

            arr Ct = mc_regression.fit_predict(t, St, Pt_delta);
            C.col(t) = Ct;

            if (method == Tsitsiklis_VanRoy) 
            {
                P.col(t) = E.col(t).max(C.col(t));
            }
            else
            {
                P.col(t) = (E.col(t) >= C.col(t)).select(E.col(t), Pt_delta);
            }
            
        }

        // For t = 0
        arr Pt_delta_0 = P.col(1) * exp(-interest_rate * dt);
        C.col(0) = mc_regression.fit_predict_at_0(Pt_delta_0);
        P.col(0) = E.col(0).max(C.col(0));
        //for (int i = 0; i < N; ++i) {
        //    std::cout << E(i, M-1) << std::endl;
        //}
        return P.col(0).mean();
    }
	


    double regression_pricer_forward(const double spot, const arr2& paths, const int w, const double strike, const double maturity,
        const double interest_rate, const regression::MCRegression& mc_regression)
    {
        int M = int(paths.cols());
        int N = int(paths.rows());
        double dt = maturity / M;

        arr2 C = arr2::Zero(N, M);
        arr2 P = arr2::Zero(N, M);
        arr2 E = (w * (paths - strike)).max(arr2::Zero(N, M));
       
        C.col(0).fill(mc_regression.predict_at_0());
        for (int t = 1; t < M - 1; t++) 
        {
            auto Ct = mc_regression.predict(t, paths.col(t));
            C.col(t) = Ct;
        }
        C.col(M - 1) = E.col(M - 1);

        for (int n = 0; n < N; n++) 
        {
            for (int j = 0; j < M; j++) 
            {
                if (E(n, j) >= C(n, j)) 
                {
                    P(n, 0) = E(n, j) * exp(-interest_rate * j * dt);
                    break;
                }
            }
        }
        return P.col(0).mean();
    }
}