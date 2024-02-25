#include <iostream>
#include <cmath>
#include <array>
#include <boost/numeric/odeint.hpp>

using state_type = std::array<double, 2>;

static constexpr double PI = 3.141592653589793;

auto getSLSystem(double lambda)
{
	return [lambda](const state_type& x, state_type& dxdt, double t) {
		dxdt[0] = x[1];
		dxdt[1] = -lambda * x[0];
		};
}

state_type getAnaliticalS(double lambda)
{
	return
	{
		sin(sqrt(lambda) * PI) / sqrt(lambda),
		cos(sqrt(lambda) * PI)
	};
}

state_type getAnaliticalC(double lambda)
{
	return
	{
		cos(sqrt(lambda) * PI),
		-sqrt(lambda) * sin(sqrt(lambda) * PI)
	};
}

double getAnaliticalDeterminant(double lambda)
{
	state_type S = getAnaliticalS(lambda);
	state_type C = getAnaliticalC(lambda);
	return
		C[1] * S[0] * S[0] +
		C[0] * S[1] * S[0] +
		C[0] * S[0] * S[1];
}

state_type getSolutionValues(double lambda, const state_type& initConditions, size_t numPoints)
{
	double step = PI / (numPoints - 1);
	state_type curx = initConditions;
	auto slSystem = getSLSystem(lambda);
	boost::numeric::odeint::runge_kutta4<state_type> rk;
	for (size_t i = 0; i != numPoints; ++i)
	{
		rk.do_step(slSystem, curx, step * i, step);
	}
	return curx;
}

state_type getNumericalS(double lambda)
{
	return getSolutionValues(lambda, { 0, 1 }, 1000);
}

state_type getNumericalC(double lambda)
{
	return getSolutionValues(lambda, { 1, 0 }, 1000);
}

double getNumericalDeterminant(double lambda)
{
	state_type S = getNumericalS(lambda);
	state_type C = getNumericalC(lambda);
	return
		C[1] * S[0] * S[0] +
		C[0] * S[1] * S[0] +
		C[0] * S[0] * S[1];
}

void print(const char* label, double v1, double v2, double v3, double v4)
{
	std::cout <<
		std::setw(15) << label << "|" <<
		std::setw(15) << v1 << "|" <<
		std::setw(15) << v2 << "|" <<
		std::setw(15) << v3 << "|" <<
		std::setw(15) << v4 << "\n";

}

int main()
{
	std::cout << "               | analitical    | numerical     | difference    | lambda\n";
	for (size_t n = 1; n != 4; ++n)
	{
		double lambda = n * n;
		state_type AS = getAnaliticalS(lambda);
		state_type AC = getAnaliticalC(lambda);
		state_type NS = getNumericalS(lambda);
		state_type NC = getNumericalC(lambda);
		double ADet = getAnaliticalDeterminant(lambda);
		double NDet = getNumericalDeterminant(lambda);
		print("S (PI, lambda)", AS[0], NS[0], abs(AS[0] - NS[0]), lambda);
		print("S'(PI, lambda)", AS[1], NS[1], abs(AS[1] - NS[1]), lambda);
		print("C (PI, lambda)", AC[0], NC[0], abs(AC[0] - NC[0]), lambda);
		print("C'(PI, lambda)", AC[1], NC[1], abs(AC[1] - NC[1]), lambda);
		print("Det", ADet, NDet, abs(ADet - NDet), lambda);
	}
}