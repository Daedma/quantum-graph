#include "QuantumGraph.hpp"
#include <boost/numeric/odeint.hpp>

QuantumGraph::state_type QuantumGraph::getSolutionValues(size_t edge, double lambda, const state_type& initConditions, size_t numPoints)
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

