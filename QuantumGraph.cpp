#include "QuantumGraph.hpp"
#include <boost/numeric/odeint.hpp>
#include <Mathter/Decompositions/DecomposeSVD.hpp>
#include <algorithm>
#include <iostream>

std::vector<double> QuantumGraph::calcEigenvalues(double lowerBound, double higherBound, double step, double binarySearchDepth) const
{
	std::vector<double> eigenvalues;
	for (; lowerBound < higherBound; lowerBound += step)
	{
		double detLeft = characteristicDeterminant(lowerBound);
		double detRight = characteristicDeterminant(lowerBound + step);
		if (detRight * detLeft <= 0)
		{
			double left = lowerBound, right = lowerBound + step;
			for (size_t it = 0; it != binarySearchDepth; ++it)
			{
				double midlle = (left + right) / 2;
				double detMiddle = characteristicDeterminant(midlle);
				if (detMiddle * detLeft <= 0)
				{
					right = midlle;
					detRight = detMiddle;
				}
				else
				{
					left = midlle;
					detLeft = detMiddle;
				}
			}
			eigenvalues.emplace_back((left + right) / 2);
		}
	}
	return eigenvalues;
}

QuantumGraph::GraphFunction QuantumGraph::calcEigenfunction(double lambda, double tolerance) const
{
	std::array<double, 3> fsos = getNullSpaceBasis(getSystemMatrix(lambda), tolerance);
	if (std::any_of(fsos.cbegin(), fsos.cend(), std::isnan<double>))
	{
		return nullptr;
	}
	std::vector<double> edge1Values = getSolutionValues(1, lambda, { 1, 0 }, numberOfNodes);
	std::vector<double> edge2Values = getSolutionValues(1, lambda, { 0, 1 }, numberOfNodes);
	std::vector<double> edge3Values = getSolutionValues(1, lambda, { 0, 1 }, numberOfNodes);

	return[y = std::array{ std::move(edge1Values), std::move(edge2Values), std::move(edge3Values) }, size = numberOfNodes, fsos]
	(std::array<double, 3> x)->std::array<double, 3> {
		constexpr double PI_INV = 1.0 / 3.14159265358979323846;

		double x1 = x[0] * PI_INV;
		double x2 = x[1] * PI_INV;
		double x3 = x[2] * PI_INV;

		size_t index1 = static_cast<size_t>(x1 * (size - 1));
		size_t index2 = static_cast<size_t>(x2 * (size - 1));
		size_t index3 = static_cast<size_t>(x3 * (size - 1));

		double y1 = y[0][index1];
		double y2 = y[1][index2];
		double y3 = y[2][index3];

		double fraction1 = x1 * (size - 1) - index1;
		double fraction2 = x2 * (size - 1) - index2;
		double fraction3 = x3 * (size - 1) - index3;

		double interpolatedValue1 = (1 - fraction1) * y1 + fraction1 * (index1 + 1 < size ? y[0][index1 + 1] : y1);
		double interpolatedValue2 = (1 - fraction2) * y2 + fraction2 * (index2 + 1 < size ? y[1][index2 + 1] : y2);
		double interpolatedValue3 = (1 - fraction3) * y3 + fraction3 * (index3 + 1 < size ? y[2][index3 + 1] : y3);

		return { fsos[0] * interpolatedValue1, fsos[1] * interpolatedValue2, fsos[2] * interpolatedValue3 };
	};
}

QuantumGraph::StateType QuantumGraph::getSolutionValuesAtPI(size_t edge, double lambda, const StateType& initConditions, size_t numPoints) const
{
	double step = PI / (numPoints - 1);
	StateType curx = initConditions;
	auto slSystem = createSturmLiouvilleODE(m_potentials[edge - 1], lambda);
	boost::numeric::odeint::runge_kutta4<StateType> rk;
	for (size_t i = 0; i != numPoints - 1; ++i)
	{
		rk.do_step(slSystem, curx, step * i, step);
	}
	return curx;
}

std::vector<double> QuantumGraph::getSolutionValues(size_t edge, double lambda, const StateType& initConditions, size_t numPoints) const
{
	double step = PI / (numPoints - 1);
	std::vector<double> values(numPoints);
	values[0] = initConditions[0];
	StateType curx = initConditions;
	auto slSystem = createSturmLiouvilleODE(m_potentials[edge - 1], lambda);
	boost::numeric::odeint::runge_kutta4<StateType> rk;
	for (size_t i = 0; i != numPoints - 1; ++i)
	{
		rk.do_step(slSystem, curx, step * i, step);
		values[i + 1] = curx[0];
	}
	return values;
}

std::array<double, 3> QuantumGraph::getNullSpaceBasis(const Matrix33& matrix, double tolerance) noexcept
{
	auto [S, V, D] = mathter::DecomposeSVD(matrix);
	for (size_t i = 0; i != 3; ++i)
	{
		if (abs(V(i, i)) <= tolerance)
		{
			return { D(i, 0), D(i, 1), D(i, 2) };
		}
	}
	return std::array<double, 3>{NAN, NAN, NAN};
}

