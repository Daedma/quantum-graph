#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <boost/math/constants/constants.hpp>
#include "QuantumGraph.hpp"


QuantumGraph genGraph(double q) noexcept
{
	return QuantumGraph{
		[q](double) { return q; },
		[q](double) { return q; },
		[q](double) { return q; }
	};
}

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

std::vector<double> getEigenvalues(double lowerBound, double upperBound, double q)
{
	std::vector<double> eigenvalues;
	for (int i = 0;; ++i)
	{
		double ev = get_null_ev(i / 3, i % 3 + 1) + q;
		if (ev >= lowerBound && ev <= upperBound)
		{
			eigenvalues.emplace_back(ev);
		}
		if (ev > upperBound)
		{
			break;
		}
	}
	return eigenvalues;
}


int main(int argc, char const* argv[])
{
	constexpr double errorStart = 1.e-10;
	constexpr double errorEnd = 1.e-2;
	constexpr double errorStep = 1.e+1;
	
	constexpr double step = 0.1;
	constexpr size_t maxIter = 10000;
	
	if (argc != 4)
	{
		std::cerr << "Usage: " << argv[0] << " <number of graphs> <lower bound> <upper bound>\n";
		return 1;
	}
	int num_of_graphs = std::atoi(argv[1]);
	double lowerBound = std::atoi(argv[2]);
	double upperBound = std::atoi(argv[3]);


	std::string filename = "experiments/error-control-gentest-" +
		std::to_string(static_cast<int>(lowerBound)) + "-" +
		std::to_string(static_cast<int>(upperBound)) + ".csv";
	std::ofstream ofs(filename);
	ofs << "Graph,Error,MaxDifference\n";

	std::uniform_real_distribution<double> dist{ -1000., 0. };
	std::mt19937 gen{ std::random_device{}() };

	for (size_t i = 0; i != num_of_graphs; ++i)
	{
		double q = dist(gen);

		QuantumGraph graph = genGraph(q);
		std::vector<double> analytical = getEigenvalues(lowerBound, upperBound, q);

		for (double error = errorStart; error <= errorEnd; error *= errorStep)
		{
			std::vector<double> numerical = graph.calcEigenvalues(lowerBound, upperBound, 2, step, error, maxIter);

			double maxDiff = 0.0;
			size_t eigenvaluesCount = std::min(analytical.size(), numerical.size());

			for (size_t j = 0; j != eigenvaluesCount; ++j)
			{
				double diff = std::abs(analytical[j] - numerical[j]);
				if (diff > maxDiff)
				{
					maxDiff = diff;
				}
			}

			ofs << i << ',' << error << ',' << maxDiff << '\n';

			std::clog << "Graph: " << i << ", Error: " << error << ", Max difference: " << maxDiff << '\n';
		}
	}

	ofs.close();
	std::clog << "Results written to " << filename << '\n';
	return 0;
}
