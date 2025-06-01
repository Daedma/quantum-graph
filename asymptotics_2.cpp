#include <fstream>
#include <iostream>
#include <iomanip>
#include "QuantumGraph.hpp"

int main(int argc, char const* argv[])
{
	constexpr double lowerBound = -10.;
	constexpr double upperBound = 10000.;
	constexpr double step = 0.1;
	constexpr double error = 1.e-6;
	constexpr size_t windowSize = 3;

	QuantumGraph graph1(
		[](double x) { return x * x; },
		[](double x) { return x + 1; },
		[](double x) { return std::sin(x); }
	);

	std::vector<double> eigenvalues1 = graph1.calcEigenvalues(lowerBound, upperBound, windowSize, step, error);

	QuantumGraph graph2(
		[](double x) { return std::cos(x); },
		[](double x) { return std::sin(-x) * std::cos(x); },
		[](double x) { return x; }
	);

	std::vector<double> eigenvalues2 = graph2.calcEigenvalues(lowerBound, upperBound, windowSize, step, error);

	QuantumGraph graph3(
		[](double x) { return 1.0; },
		[](double x) { return 2.0; },
		[](double x) { return 3.0; }
	);

	std::vector<double> eigenvalues3 = graph3.calcEigenvalues(lowerBound, upperBound, windowSize, step, error);

	std::string filename = "experiments/eigenvalues-" + std::to_string(int(lowerBound)) + "-" + std::to_string(int(upperBound)) + ".csv";
	std::ofstream ofs{ filename };
	ofs << std::setprecision(16);
	ofs << "graph#1";
	for (auto i : eigenvalues1)
	{
		ofs << ", " << i;
	}
	ofs << "\ngraph#2";
	for (auto i : eigenvalues2)
	{
		ofs << ", " << i;
	}
	ofs << "\ngraph#3";
	for (auto i : eigenvalues3)
	{
		ofs << ", " << i;
	}

	std::clog << "Program finished";
	return 0;
}
