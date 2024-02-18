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
	for (auto i : graph1.getEigenvalues(0, 10))
	{
		std::cout << i << ' ';
	}
	return 0;
}