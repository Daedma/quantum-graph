#include <fstream>
#include "QuantumGraph.hpp"

int main(int argc, char const* argv[])
{
	constexpr double upperBound = 1000.;
	constexpr double lowerBound = -2.;
	constexpr double step = 0.1;
	constexpr size_t depth = 30;
	constexpr size_t numberOfNobes = 3000;

	QuantumGraph graph1{
		[](double) {return 0.;},
		[](double) {return 0.;},
		[](double) {return 0.;}
	};
	graph1.setNumberOfNodes(numberOfNobes);
	std::vector<double> eigenvalues1 = graph1.calcEigenvalues(lowerBound, upperBound, step, depth);

	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	graph2.setNumberOfNodes(numberOfNobes);
	std::vector<double> eigenvalues2 = graph2.calcEigenvalues(lowerBound, upperBound, step, depth);

	QuantumGraph graph3{
		[](double x) {return sin(x);},
		[](double x) {return exp(-x);},
		[](double x) {return exp(-x * x);}
	};
	graph3.setNumberOfNodes(numberOfNobes);
	std::vector<double> eigenvalues3 = graph3.calcEigenvalues(lowerBound, upperBound, step, depth);

	std::ofstream ofs{ "asymptotics-2-out.csv" };
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
	return 0;
}
