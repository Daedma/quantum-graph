#include "QuantumGraph.hpp"
#include <cmath>
#include <iostream>
#include <boost/math/constants/constants.hpp>
#include <matplot/matplot.h>

void plotEigenfunction(const QuantumGraph& graph, const char* style)
{
	std::vector<double> eigenvalues = graph.calcEigenvalues(0, 5);
	if (eigenvalues.empty())
	{
		return;
	}
	QuantumGraph::GraphFunction eigenfunction = graph.calcEigenfunction(eigenvalues[0]);
	if (eigenfunction)
	{
		auto func1 = [eigenfunction](double x) {return eigenfunction({ x, 0, 0 })[0];};
		auto func2 = [eigenfunction](double x) {return eigenfunction({ 0, x, 0 })[1];};
		auto func3 = [eigenfunction](double x) {return eigenfunction({ 0, 0, x })[2];};
		matplot::fplot(func1, std::array{ 0., boost::math::double_constants::pi }, style);
		matplot::hold(matplot::on);
		matplot::fplot(func2, std::array{ 0., boost::math::double_constants::pi }, style);
		matplot::hold(matplot::on);
		matplot::fplot(func3, std::array{ 0., boost::math::double_constants::pi }, style);
		matplot::hold(matplot::on);
	}
}

int main(int argc, char const* argv[])
{
	QuantumGraph graph1{
		[](double) {return 0.;},
		[](double) {return 0.;},
		[](double) {return 0.;}
	};
	plotEigenfunction(graph1, "r");

	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	plotEigenfunction(graph2, "b");

	QuantumGraph graph3{
		[](double x) {return sin(x);},
		[](double x) {return exp(-x);},
		[](double x) {return exp(-pow(x, 2));}
	};
	plotEigenfunction(graph3, "g");
	matplot::show();
}
