#include <iostream>
#include <cmath>
#include <boost/multi_array.hpp>
#include "QuantumGraph.hpp"

void checkSequences(const std::vector<double> seq1, const std::vector<double> seq2)
{
	boost::multi_array<double, 3> permutations(boost::extents[3][3][std::max(seq1.size(), seq2.size()) / 3 + 2]);
	for (size_t i = 0; i != 3; ++i)
	{
		for (size_t j = 0; j != 3; ++j)
		{
			for (size_t n = 0; 3 * n < std::min(seq1.size() - i, seq2.size() - j); ++n)
			{
				permutations[i][j][n] = (abs(sqrt(seq1[3 * n + i]) - sqrt(seq2[3 * n + j])) * n);
			}
		}
	}
	boost::multi_array<bool, 3> decreasing(boost::extents[3][3][std::max(seq1.size(), seq2.size()) / 3 + 2]);
	for (size_t i = 0; i != 3; ++i)
	{
		for (size_t j = 0; j != 3; ++j)
		{
			decreasing[i][j][0] = false;
			for (size_t n = 0; n + 2 < decreasing[i][j].size(); ++n)
			{
				decreasing[i][j][n] = (abs(permutations[i][j][n + 1] - permutations[i][j][n + 2]) < abs(permutations[i][j][n] - permutations[i][j][n + 1]));
			}
		}
	}
	for (size_t i = 0; i != 3; ++i)
	{
		for (size_t j = 0; j != 3; ++j)
		{
			std::cout << "i=" << i << ", j=" << j << " : ";
			for (size_t n = 0; n != permutations[i][j].size(); ++n)
			{
				std::cout << permutations[i][j][n] << ' ';
			}
			std::cout << "\ndecreasing : ";
			for (size_t n = 0; n != decreasing[i][j].size(); ++n)
			{
				std::cout << decreasing[i][j][n] << ' ';
			}
			std::cout << "\n";
		}
	}
}

int main()
{
	constexpr double upperBound = 3000.;

	QuantumGraph graph1{
		[](double) {return 0.;},
		[](double) {return 0.;},
		[](double) {return 0.;}
	};
	std::vector<double> eigenvalues0 = graph1.calcEigenvalues(0, upperBound, 0.2, 30);

	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	std::vector<double> eigenvaluesGraph2 = graph2.calcEigenvalues(0, upperBound, 0.2, 30);

	QuantumGraph graph3{
		[](double x) {return sin(x);},
		[](double x) {return exp(-x);},
		[](double x) {return exp(-pow(x, 2));}
	};
	std::vector<double> eigenvaluesGraph3 = graph3.calcEigenvalues(0, upperBound, 0.2, 30);

	std::cout << "Graph 0 & Graph 2 :\n";
	checkSequences(eigenvalues0, eigenvaluesGraph2);

	std::cout << "Graph 0 & Graph 3 :\n";
	checkSequences(eigenvalues0, eigenvaluesGraph3);

	return 0;
}