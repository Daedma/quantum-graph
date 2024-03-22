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
	std::vector<double> eigenvalues0 = graph1.getEigenvalues(0, 200);

	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	std::vector<double> eigenvaluesGraph2 = graph2.getEigenvalues(0, 200);

	QuantumGraph graph3{
		[](double x) {return sin(x);},
		[](double x) {return exp(-x);},
		[](double x) {return exp(pow(-x, 2));}
	};
	std::vector<double> eigenvaluesGraph3 = graph3.getEigenvalues(0, 200);

	for (size_t n = 0; n < std::min(eigenvalues0.size() + 2, eigenvaluesGraph2.size()); n += 3)
	{
		std::cout << abs(sqrt(eigenvalues0[n + 2]) - sqrt(eigenvaluesGraph2[n])) * ((n + 1) / 3 + 1) << " ";
	}
	std::cout << std::endl;

	for (size_t n = 1; n < std::min(eigenvalues0.size() + 2, eigenvaluesGraph2.size()); n += 3)
	{
		std::cout << abs(sqrt(eigenvalues0[n + 2]) - sqrt(eigenvaluesGraph2[n])) * (n / 3 + 1) << " ";
	}
	std::cout << std::endl;

	for (size_t n = 2; n < std::min(eigenvalues0.size() + 2, eigenvaluesGraph2.size()); n += 3)
	{
		std::cout << abs(sqrt(eigenvalues0[n + 2]) - sqrt(eigenvaluesGraph2[n])) * ((n - 1) / 3 + 1) << " ";
	}
	std::cout << std::endl;

	for (size_t n = 0; n < std::min(eigenvalues0.size() + 2, eigenvaluesGraph3.size()); n += 3)
	{
		std::cout << abs(sqrt(eigenvalues0[n + 2]) - sqrt(eigenvaluesGraph3[n])) * ((n + 1) / 3 + 1) << " ";
	}
	std::cout << std::endl;

	for (size_t n = 1; n < std::min(eigenvalues0.size() + 2, eigenvaluesGraph3.size()); n += 3)
	{
		std::cout << abs(sqrt(eigenvalues0[n + 2]) - sqrt(eigenvaluesGraph3[n])) * (n / 3 + 1) << " ";
	}
	std::cout << std::endl;

	for (size_t n = 2; n < std::min(eigenvalues0.size() + 2, eigenvaluesGraph3.size()); n += 3)
	{
		std::cout << abs(sqrt(eigenvalues0[n + 2]) - sqrt(eigenvaluesGraph3[n])) * ((n - 1) / 3 + 1) << " ";
	}
	std::cout << std::endl;

	return 0;
}