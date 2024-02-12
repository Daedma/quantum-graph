#include <boost/numeric/odeint.hpp>
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

int main()
{
	auto sl1 = sturm_liouville_solution([](double x) {return 2;}, 2, { 1, 0 }, 1000);
	auto sl2 = sturm_liouville_solution([](double x) {return 2;}, 2, { 0, 1 }, 1000);
	matplot::plot(sl1.first, sl1.second, "r", sl2.first, sl2.second, "g");
	matplot::show();
	return 0;
}