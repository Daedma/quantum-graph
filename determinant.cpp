#include "QuantumGraph.hpp"
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <matplot/matplot.h>
#include <iostream>

void plot_det()
{
	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	auto numeric_det = [&graph2](double lambda) {return std::sqrt(lambda) * graph2.characteristicDeterminant(lambda);};

	auto eigenvalue = [](uint32_t n) {
		using boost::math::double_constants::pi;
		using boost::math::double_constants::root_pi;
		using boost::math::double_constants::quarter_pi;
		using boost::math::double_constants::one_div_pi;
		constexpr double B1 = quarter_pi;
		constexpr double B2 = -quarter_pi;
		constexpr double B3 = -root_pi / 3;
		double s_n = 0.5 * B1 + 0.25 * B2 + 0.25 * B3;
		return one_div_pi * std::asin(std::sqrt(2 / 3)) + n + s_n / n;
		};

	auto asymptotic_det = [](double lambda) {
		using boost::math::double_constants::pi;
		using boost::math::double_constants::root_pi;
		using boost::math::double_constants::quarter_pi;
		constexpr double B1 = quarter_pi;
		constexpr double B2 = -quarter_pi;
		constexpr double B3 = -root_pi / 3;
		double rho = std::sqrt(lambda);
		double sin_rho = std::sin(rho * pi);
		double sin_rho2 = sin_rho * sin_rho;
		double cos_rho = std::cos(rho * pi);

		return sin_rho * (2 - 3 * sin_rho2) + (3 * (B1 + B2 + B3) * pi * sin_rho2 * cos_rho - (B2 + B3) * pi * cos_rho) / rho;
		};

	auto lambda = matplot::linspace(100, 200, 100);
	matplot::plot(lambda, matplot::transform(lambda, asymptotic_det), lambda, matplot::transform(lambda, numeric_det));
	matplot::legend({ "num", "asympt" });
	matplot::xlabel("lambda");
	matplot::ylabel("Charachteristic function");
	matplot::show();
}

void test_eigenvalues_asympt(uint32_t lower, uint32_t upper)
{
	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};

	auto eigenvalue = [](uint32_t n) {
		using boost::math::double_constants::pi;
		using boost::math::double_constants::root_pi;
		using boost::math::double_constants::quarter_pi;
		using boost::math::double_constants::one_div_pi;
		constexpr double B1 = quarter_pi;
		constexpr double B2 = -quarter_pi;
		constexpr double B3 = -root_pi / 3;
		double s_n = 0.5 * B1 + 0.25 * B2 + 0.25 * B3;
		double rho = one_div_pi * std::asin(std::sqrt(2. / 3)) + n + s_n / n;
		double lambda = rho * rho;
		return lambda;
		};

	for (uint32_t n = lower; n < upper; ++n)
	{
		std::cout << "n : " << n << " | det : " << graph2.characteristicDeterminant(eigenvalue(n)) << std::endl;
	}
}

int main(int argc, char const* argv[])
{
	test_eigenvalues_asympt(500, 600);

	return 0;
}
