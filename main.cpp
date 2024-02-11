#include <boost/numeric/odeint.hpp>
#include <iostream>
#include <functional>
#include <array>
#include <utility>

auto create_sturm_liouville_ode(std::function<double(double)> q, double lambda) noexcept
{
	using state_type = std::array<double, 2>;
	return [q, lambda](const state_type& x, state_type& dxdt, double t) {
		dxdt[0] = x[1];
		dxdt[1] = q(t) * x[1] - lambda * x[1];
		};
}

int main()
{
	std::cout << "Hello world\n" << std::endl;
}