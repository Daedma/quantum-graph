#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <matplot/matplot.h>
#include <iostream>
#include <functional>
#include <array>
#include <utility>
#include <vector>
#include <set>
#include <algorithm>
#include <numeric>

using state_type = std::array<double, 2>;

static const double PI = 3.141592653589793238462643383;

auto create_sturm_liouville_ode(std::function<double(double)> q, double lambda) noexcept
{
	return [q, lambda](const state_type& x, state_type& dxdt, double t) {
		dxdt[0] = x[1];
		dxdt[1] = q(t) * x[0] - lambda * x[0];
		};
}

std::pair<std::vector<double>, std::vector<double>>
sturm_liouville_solution(std::function<double(double)> q, double lambda, const state_type& init_conditions, size_t num_point, std::pair<double, double> bounds = { 0., PI })
{
	assert(num_point > 1 && bounds.first < bounds.second);
	double step = (bounds.second - bounds.first) / (num_point - 1);
	std::vector<double> x(num_point);
	x[0] = bounds.first;
	std::generate(x.begin() + 1, x.end(), [t = bounds.first, step]() mutable {return t += step;});
	std::vector<double> y(num_point);
	y[0] = init_conditions[0];
	state_type curx = init_conditions;
	auto sl_sys = create_sturm_liouville_ode(q, lambda);
	boost::numeric::odeint::runge_kutta4<state_type> rk;
	for (size_t i = 1; i != num_point; ++i)
	{
		rk.do_step(sl_sys, curx, x[i], step);
		y[i] = curx[0];
	}
	return { std::move(x), std::move(y) };
}

std::pair<double, double> get_CS(double q, double lambda, double x) noexcept
{
	if (q > lambda)
	{
		return {
			0.5 * exp(sqrt(q - lambda) * x) + 0.5 * exp(-sqrt(q - lambda) * x),
			1. / (2 * sqrt(q - lambda)) * exp(sqrt(q - lambda) * x) - 1. / (2 * sqrt(q - lambda)) * exp(-sqrt(q - lambda) * x)
		};
	}
	if (q < lambda)
	{
		return {
			cos(sqrt(lambda - q) * x),
			1. / sqrt(lambda - q) * sin(sqrt(lambda - q) * x)
		};
	}
	return {
		1.,
		x
	};
}

int main()
{
	size_t n = 1000;

	double q = 2, lambda = 1;
	auto CS = get_CS(q, lambda, PI);

	auto slS1 = sturm_liouville_solution([q](double x) {return q;}, lambda, { 1, 0 }, n);
	auto slC1 = sturm_liouville_solution([q](double x) {return q;}, lambda, { 0, 1 }, n);

	std::cout << "              |  numerical  |  analytical | q | lambda\n";
	std::cout << "C(PI, lambda) | " << std::setw(12) << slS1.second[n - 1]
		<< '|' << std::setw(12) << CS.first << " | " << q << " | " << lambda << '\n';
	std::cout << "S(PI, lambda) | " << std::setw(12) << slC1.second[n - 1]
		<< '|' << std::setw(12) << CS.second << " | " << q << " | " << lambda << '\n';

	q = 1, lambda = 2;
	CS = get_CS(q, lambda, PI);;
	auto slS2 = sturm_liouville_solution([q](double x) {return q;}, lambda, { 1, 0 }, n);
	auto slC2 = sturm_liouville_solution([q](double x) {return q;}, lambda, { 0, 1 }, n);

	std::cout << "C(PI, lambda) | " << std::setw(12) << slS2.second[n - 1]
		<< '|' << std::setw(12) << CS.first << " | " << q << " | " << lambda << '\n';
	std::cout << "S(PI, lambda) | " << std::setw(12) << slC2.second[n - 1]
		<< '|' << std::setw(12) << CS.second << " | " << q << " | " << lambda << '\n';

	q = 1, lambda = 1;
	CS = get_CS(q, lambda, PI);;
	auto slS3 = sturm_liouville_solution([q](double x) {return q;}, lambda, { 1, 0 }, n);
	auto slC3 = sturm_liouville_solution([q](double x) {return q;}, lambda, { 0, 1 }, n);

	std::cout << "C(PI, lambda) | " << std::setw(12) << slS3.second[n - 1]
		<< '|' << std::setw(12) << CS.first << " | " << q << " | " << lambda << '\n';
	std::cout << "S(PI, lambda) | " << std::setw(12) << slC3.second[n - 1]
		<< '|' << std::setw(12) << CS.second << " | " << q << " | " << lambda << '\n';

	matplot::plot(slS1.first, slS1.second, "r",
		slC1.first, slC1.second, "g",
		slS2.first, slS2.second, "b",
		slC2.first, slC2.second, "y",
		slS3.first, slS3.second, "c",
		slC3.first, slC3.second, "k");
	matplot::show();
	return 0;
}