#include "QuantumGraph.hpp"
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include <matplot/matplot.h>

int main(int argc, char const *argv[])
{
	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	auto numeric_det = [&graph2](double lambda) {return graph2.characteristicDeterminant(lambda);};

	auto asymptotic_det = [](double lambda) {
		using boost::math::double_constants::pi;
		double rho = std::sqrt(lambda);
		return std::sin(rho*pi)/rho * (2 - 3*sin(rho*pi)*sin(rho*pi));
	};

	auto lambda = matplot::linspace(1000, 1100, 100);
	matplot::plot(lambda, matplot::transform(lambda, asymptotic_det), lambda, matplot::transform(lambda, numeric_det));
	matplot::legend({"num", "asympt"});
	matplot::xlabel("lambda");
	matplot::ylabel("Charachteristic function");
	matplot::show();
	return 0;
}
