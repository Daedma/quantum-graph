#include <iostream>
#include <cmath>
#include <array>
#include <boost/numeric/odeint.hpp>
#include <boost/math/constants/constants.hpp>
#include <matplot/matplot.h>
#include "QuantumGraph.hpp"

using state_type = std::array<double, 2>;

void plotEigenfunction(const QuantumGraph& graph, const char* style)
{
	std::vector<double> eigenvalues = graph.calcEigenvalues(0, 2);
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

auto getAnaliticalS(double lambda)
{
	return [lambda](double x) {
		return sin(sqrt(lambda) * x) / sqrt(lambda);
		};
}

auto getAnaliticalC(double lambda)
{
	return [lambda](double x) {
		return cos(sqrt(lambda) * x);
		};
}

auto getAnaliticalEigenfunction(double lambda, size_t k)
{
	namespace constants = boost::math::double_constants;
	std::array<double, 3> coeffs;
	switch (k)
	{
	case 1:
		coeffs = { 0., 1., -1. };
		break;
	case 2:
		coeffs = { constants::root_two, 1., 1. };
		break;
	case 3:
		coeffs = { -constants::root_two, 1., 1. };
		break;
	default:
		break;
	}
	return[rho = sqrt(lambda), coeffs](const std::array<double, 3>& x)->std::array<double, 3>
	{
		return {
			coeffs[0] * cos(rho * x[0]),
			coeffs[1] * sin(rho * x[1]),
			coeffs[2] * sin(rho * x[2])
		};
	};
}

double getEigenvalue(size_t n, size_t k)
{
	namespace constants = boost::math::double_constants;
	constexpr double term = constants::one_div_pi * 0.955316618125;
	switch (k)
	{
	case 1:
		return n * n;
	case 2:
		return pow(term + n, 2);
	case 3:
		return pow(-term + n + 1, 2);
	default:
		break;
	}
	return NAN;
}

// 1. Сравнить численное решение при q1 = q2 = q3 = 0 с аналитическим, которое Вы находили.
// Вывести несколько первых собственных значений, найденных численно и аналитически, посчитать 
// разницу между ними, а также изобразить на графике компоненты собственных функций.

// 2. Представить в таблице собственные значения Graph 0, Graph 1, Graph 3. Проследить за тем, 
// что для малых значений n они могут сильно отличаться, а при больших номерах n все три 
// последовательности становятся близкими.

int main(int argc, char const* argv[])
{
	QuantumGraph graph1{
		[](double) {return 0.;},
		[](double) {return 0.;},
		[](double) {return 0.;}
	};
	std::vector<double> eigenvaluesNum = graph1.calcEigenvalues(0, 100, 0.2, 30);
	std::vector<double> eigenvaluesAnal(eigenvaluesNum.size() + 1); // first eigenvalue is missed
	for (size_t i = 0; i != eigenvaluesAnal.size(); ++i)
	{
		eigenvaluesAnal[i] = getEigenvalue(i / 3, i % 3 + 1);
	}
	std::cout << "Numerical: ";
	for (auto i : eigenvaluesNum)
	{
		std::cout << i << ", ";
	}
	std::cout << "\nAnalitical: ";
	for (auto i : eigenvaluesAnal)
	{
		std::cout << i << ", ";
	}
	std::cout << "\nDifference: ";
	for (size_t i = 0; i != eigenvaluesNum.size(); ++i)
	{
		std::cout << abs(eigenvaluesAnal[i + 1] - eigenvaluesNum[i]) << ", ";
	}

	matplot::fplot([func = getAnaliticalEigenfunction(eigenvaluesAnal[1], 2)](double x) {return func({ x, 0., 0. })[0];}, std::array{ 0., boost::math::double_constants::pi }, "r");
	matplot::hold(matplot::on);
	matplot::fplot([func = getAnaliticalEigenfunction(eigenvaluesAnal[1], 2)](double x) {return func({ 0., x, 0. })[1];}, std::array{ 0., boost::math::double_constants::pi }, "r");
	matplot::hold(matplot::on);
	matplot::fplot([func = getAnaliticalEigenfunction(eigenvaluesAnal[1], 2)](double x) {return func({ 0., 0., x })[2];}, std::array{ 0., boost::math::double_constants::pi }, "r");
	matplot::hold(matplot::on);

	plotEigenfunction(graph1, "b");
	matplot::show();
	return 0;
}
