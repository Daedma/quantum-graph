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
	constexpr double lowerBound = 999999.0;
	constexpr double upperBound = 1100000.0;

	constexpr double errorStart = 1.e-6;
	constexpr double errorEnd = 1.e+0;
	constexpr double errorStep = 1.e+1;

	constexpr double step = 100.;
	constexpr size_t maxIter = 10000;

	std::vector<double> analytical;
	for (int i = 0;; ++i)
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

	std::string filename = "experiments/error-control-test-" +
		std::to_string(static_cast<int>(lowerBound)) + "-" +
		std::to_string(static_cast<int>(upperBound)) + ".csv";
	std::ofstream ofs(filename);
	ofs << "Error,MaxDifference\n";

	for (double error = errorStart; error <= errorEnd; error *= errorStep)
	{
		QuantumGraph graph{
			[](double) { return 0.; },
			[](double) { return 0.; },
			[](double) { return 0.; }
		};

		std::vector<double> numerical = graph.calcEigenvalues(lowerBound, upperBound, 2, step, error, maxIter);

		double maxDiff = 0.0;
		size_t eigenvaluesCount = std::min(analytical.size(), numerical.size());

		for (size_t i = 0; i != eigenvaluesCount; ++i)
		{
			double diff = std::abs(analytical[i] - numerical[i]);
			if (diff > maxDiff)
			{
				maxDiff = diff;
			}
		}

		ofs << error << ',' << maxDiff << '\n';

		std::clog << "Error: " << error << ", Max difference: " << maxDiff << '\n';
	}

	ofs.close();
	std::clog << "Results written to " << filename << '\n';
	return 0;
}
