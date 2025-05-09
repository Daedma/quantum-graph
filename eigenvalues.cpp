#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "QuantumGraph.hpp"


double get_null_ev(int n, int k) noexcept
{
	using boost::math::double_constants::one_div_pi;
	static const double main_term = one_div_pi * std::asin(std::sqrt(2. / 3.));
	double ev;
	switch (k)
	{
	case 1:
		ev = main_term + n;
		break;
	case 2:
		ev = -main_term + n + 1;
		break;
	case 3:
		ev = n + 1;
		break;
	default:
		ev = NAN;
		break;
	}
	return ev * ev;
}


int main(int argc, char const* argv[])
{
	constexpr double lowerBound = 99.0;
	constexpr double upperBound = 200.0;
	constexpr double step = 5.;
	constexpr double error = 1.e-8;
	constexpr size_t maxIter = 1000;

	std::vector<double> analytical;
	for (int i = 0;;++i)
	{
		double ev = get_null_ev(i / 3, i % 3 + 1);
		assert(!std::isnan(ev));
		if (ev >= lowerBound && ev <= upperBound)
		{
			analytical.emplace_back(ev);
		}
		if (ev > upperBound)
		{
			break;
		}
	}

	QuantumGraph graph{
        [](double) {return 0.;},
        [](double) {return 0.;},
        [](double) {return 0.;}
    };

	std::vector<double> numerical = graph.calcEigenvalues(lowerBound, upperBound, step, error, maxIter);

	double maxDiff = 0.0;
	size_t eigenvaluesCount = std::min(analytical.size(), numerical.size());

	std::cout << std::setprecision(17);

	for (size_t i = 0; i != eigenvaluesCount; ++i)
	{
		double diff = std::abs(analytical[i] - numerical[i]);
		if (diff > maxDiff)
		{
			maxDiff = diff;
		}
		std::cout << analytical[i] << ' ' << numerical[i] << ' ' << diff << '\n';
	}

	std::cout << "Max difference: " << std::scientific << maxDiff << '\n';
}
