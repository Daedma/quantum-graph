#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "QuantumGraph.hpp"


double get_null_ev(int n, int k) noexcept
{
	double ev;
	switch (k)
	{
	case 1:
		ev = n;
		break;
	case 2:
		ev = n + 0.5;
		break;
	default:
		ev = NAN;
		break;
	}
	return ev * ev;
}


int main(int argc, char const* argv[])
{
	constexpr double lowerBound = -1.0;
	constexpr double upperBound = 200.0;
	constexpr double step = 0.1;
	constexpr double error = 1.e-7;
	constexpr size_t maxIter = 10000;

	std::vector<double> analytical;
	for (int i = 0;;++i)
	{
		double ev = i * i;
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
		[](double) {return 0.;},
		QuantumGraph::eBoundaryCondition::Neumann,
		QuantumGraph::eBoundaryCondition::Neumann,
		QuantumGraph::eBoundaryCondition::Neumann
	};

	// std::vector<double> numerical = graph.calcEigenvalues(lowerBound, upperBound, step, error, maxIter);
	std::vector<double> numerical = graph.calcEigenvalues(lowerBound, upperBound, 2, step, error, maxIter);

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
