#include "QuantumGraph.hpp"
#include <cmath>
#include <array>
#include <boost/math/constants/constants.hpp>
#include <matplot/matplot.h>

void plotEigenfunction(const QuantumGraph& graph,
	const std::array<const char*, 3>& styles,
	const std::array<const char*, 3>& legends)
{
	std::vector<double> eigenvalues = graph.calcEigenvalues(0, 5);
	if (eigenvalues.empty())
	{
		return;
	}
	QuantumGraph::GraphFunction eigenfunction = graph.calcEigenfunction(eigenvalues[3]);
	if (eigenfunction)
	{
		auto func1 = [eigenfunction](double x) {return eigenfunction({ x, 0, 0 })[0];};
		auto func2 = [eigenfunction](double x) {return eigenfunction({ 0, x, 0 })[1];};
		auto func3 = [eigenfunction](double x) {return eigenfunction({ 0, 0, x })[2];};
		auto p1 = matplot::fplot(func1, std::array{ 0., boost::math::double_constants::pi }, styles[0]);
		p1->display_name(legends[0]);
		matplot::hold(matplot::on);
		auto p2 = matplot::fplot(func2, std::array{ 0., boost::math::double_constants::pi }, styles[1]);
		p2->display_name(legends[1]);
		matplot::hold(matplot::on);
		auto p3 = matplot::fplot(func3, std::array{ 0., boost::math::double_constants::pi }, styles[2]);
		p3->display_name(legends[2]);
		matplot::hold(matplot::on);
		auto lgd = matplot::legend(matplot::on);
		lgd->location(matplot::legend::general_alignment::bottom);
		lgd->num_columns(8);
	}
}

int main(int argc, char const* argv[])
{
	QuantumGraph graph1{
		[](double) {return 0.;},
		[](double) {return 0.;},
		[](double) {return 0.;}
	};
	plotEigenfunction(graph1, { "r", "--r", ":r" }, { "q_0_1", "q_0_2", "q_0_3" });

	QuantumGraph graph2{
		[](double x) {return x;},
		[](double x) {return -x;},
		[](double x) {return -sqrt(x);}
	};
	plotEigenfunction(graph2, { "b", "--b", ":b" }, { "q_1_1", "q_1_2", "q_1_3" });

	QuantumGraph graph3{
		[](double x) {return sin(x);},
		[](double x) {return exp(-x);},
		[](double x) {return exp(-pow(x, 2));}
	};
	plotEigenfunction(graph3, { "g", "--g", ":g" }, { "q_2_1", "q_2_2", "q_2_3" });
	matplot::ylabel("q_i_j(x_j)");
	matplot::xlabel("x_j\nX\nX\nX\nX\nX");
	matplot::grid(true);
	matplot::show();
}
