#include <iostream>
#include <cmath>
#include "QuantumGraph.hpp"

int main()
{
	QuantumGraph graph1{
		[](double) {return 0.;},
		[](double) {return 0.;},
		[](double) {return 0.;}
	};
	std::cout << "0, 0, 0 : ";
	for (auto i : graph1.getEigenvalues(0, 100))
	{
		std::cout << i << ' ';
	}

	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	std::cout << "\nx, -x, -sqrt(x) : ";
	for (auto i : graph2.getEigenvalues(0, 100))
	{
		std::cout << i << ' ';
	}

	QuantumGraph graph3{
		[](double x) {return sin(x);},
		[](double x) {return exp(-x);},
		[](double x) {return exp(pow(-x, 2));}
	};
	std::cout << "\nsin(x), exp(-x), exp(pow(-x, 2)) : ";
	for (auto i : graph3.getEigenvalues(0, 100))
	{
		std::cout << i << ' ';
	}
	return 0;
}