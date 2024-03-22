#include "QuantumGraph.hpp"
#include <boost/numeric/odeint.hpp>
#include <Mathter/Decompositions/DecomposeSVD.hpp>

QuantumGraph::state_type QuantumGraph::getSolutionValues(size_t edge, double lambda, const state_type& initConditions, size_t numPoints) const
{
	double step = PI / (numPoints - 1);
	state_type curx = initConditions;
	auto slSystem = create_sturm_liouville_ode(m_potentials[edge - 1], lambda);
	boost::numeric::odeint::runge_kutta4<state_type> rk;
	for (size_t i = 0; i != numPoints - 1; ++i)
	{
		rk.do_step(slSystem, curx, step * i, step);
	}
	return curx;
}

std::array<double, 3> QuantumGraph::getNullSpaceBasis(const matrix& m) const
{
	constexpr double EPS = 1.e-10;
	auto [S, V, D] = mathter::DecomposeSVD(m);
	D = mathter::Transpose(D);
	for (size_t i = 0; i != 3; ++i)
	{
		if (abs(V(i, i)) <= EPS)
		{
			return { D(0, i), D(1, i), D(2, i) };
		}
	}
	return std::array<double, 3>{NAN, NAN, NAN};
}

