#ifndef LinearSolvers_HPP
#define LinearSolvers_HPP
#include <Eigen/Dense>
#include <tuple>


typedef Eigen::VectorXd vec;
typedef Eigen::MatrixXd mat;
typedef Eigen::PermutationMatrix<-1, -1, unsigned int> permutation;
typedef unsigned int uint;
enum StoppingCriterion { consecutive, residual };

class LinearSolvers
{
public:
	//Constructors and Distructor
	LinearSolvers() = default;
	LinearSolvers(LinearSolvers& source) = default;
	~LinearSolvers() = default;

	//Forward and Backward substitution
	vec forward_subst(const mat& L, const vec& b);
	vec backward_subst(const mat& U, const vec& b);

	//LU with or without pivoting
	std::tuple<mat, mat > lu_no_pivoting(mat A);
	std::tuple<permutation, mat, mat > lu_row_pivoting(mat A);

	//Cholesky Decomposition
	mat cholesky(mat A);

	//Iteration method
	std::tuple<vec, uint> gauss_seidel(const mat& A, const vec& b, vec x, const double tolerance, const StoppingCriterion criterion);
	std::tuple<vec, uint> jacobi(const mat& A, const vec& b, vec x, const double tolerance, const StoppingCriterion criterion);
	std::tuple<vec, uint> sor(const mat& A, const vec& b, vec x, const double tolerance, const double omega, const StoppingCriterion criterion);

	//Iteration method for banded matrix
	std::tuple<vec, uint> gauss_seidel_banded(const mat& A, const vec& b, vec x, const double tolerance, const int lower_bandwidth, const int upper_bandwidth, const StoppingCriterion criterion);
	std::tuple<vec, uint> jacobi_banded(const mat& A, const vec& b, vec x, const double tolerance, const int lower_bandwidth, const int upper_bandwidth, const StoppingCriterion criterion);
	std::tuple<vec, uint> sor_banded(const mat& A, const vec& b, vec x, const double tolerance, const double omega, const int lower_bandwidth, const int upper_bandwidth, const StoppingCriterion criterion);

};


#endif // !LinearSolvers_HPP
