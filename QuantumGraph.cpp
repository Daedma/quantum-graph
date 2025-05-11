#include "QuantumGraph.hpp"

#include <boost/numeric/odeint.hpp>
#include <Mathter/Decompositions/DecomposeSVD.hpp>

#include <algorithm>

std::vector<double> QuantumGraph::calcEigenvalues(double lowerBound, double higherBound, double step, double error, size_t maxIter) const
{
	std::vector<double> eigenvalues;
	for (; lowerBound < higherBound; lowerBound += step)
	{
		double detLeft = characteristicDeterminantSignSafe(lowerBound, error, true);
		double detRight = characteristicDeterminantSignSafe(lowerBound + step, error, true);
		if (detRight * detLeft <= 0)
		{
			double detError = error;
			double left = lowerBound, right = lowerBound + step;
			for (size_t it = 0; it != maxIter && (right - left) > 2 * error; ++it)
			{
				double midlle = (left + right) * 0.5;
				double detMiddle = characteristicDeterminantSignSafe(midlle, detError, false);
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
			eigenvalues.emplace_back((left + right) * 0.5);

		}
	}
	return eigenvalues;
}

std::vector<double> QuantumGraph::calcEigenvalues(double lowerBound, double higherBound, size_t windowSize, double initialStep, double error, size_t maxIter) const
{
	std::vector<double> eigenvalues;

	std::vector<double> window(windowSize, initialStep);
	size_t windowIndex = 0;
	double lastEigenvalue = lowerBound;
	double currentStep = initialStep;

	for (; lowerBound < higherBound; lowerBound += currentStep)
	{
		double detLeft = characteristicDeterminantSignSafe(lowerBound, error, true);
		double detRight = characteristicDeterminantSignSafe(lowerBound + currentStep, error, true);
		if (detRight * detLeft <= 0)
		{
			double detError = error;
			double left = lowerBound, right = lowerBound + currentStep;
			for (size_t it = 0; it != maxIter && (right - left) > 2 * error; ++it)
			{
				double midlle = (left + right) * 0.5;
				double detMiddle = characteristicDeterminantSignSafe(midlle, detError, false);
				if (detMiddle == 0.)
				{
					left = right = midlle;
					break;
				}
				else if (detMiddle * detLeft <= 0)
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
			double eigenvalue = (left + right) * 0.5;
			double diff = eigenvalue - lastEigenvalue;
			if (diff > 2 * error)
			{
				window[windowIndex] = diff;
				windowIndex = (windowIndex + 1) % windowSize;
				currentStep = *std::min_element(window.cbegin(), window.cend());
				eigenvalues.emplace_back(eigenvalue);
				lastEigenvalue = eigenvalue;
			}
		}
	}
	return eigenvalues;
}

std::vector<QuantumGraph::GraphFunction> QuantumGraph::calcEigenfunction(double lambda, size_t numberOfNodes, double error) const
{
	// fsos - fundametal system of solutions
	std::vector<std::array<double, 3>> fsos = getNullSpaceBasis(getSystemMatrix(lambda, error * 1.e-2), error);
	if (fsos.empty())
	{
		return {};
	}
	std::vector<double> edge1Values = getSolutionValues(1, lambda, { 1, 0 }, numberOfNodes);
	std::vector<double> edge2Values = getSolutionValues(2, lambda, { 0, 1 }, numberOfNodes);
	std::vector<double> edge3Values = getSolutionValues(3, lambda, { 0, 1 }, numberOfNodes);

	std::vector<GraphFunction> eigenfunctions;
	eigenfunctions.reserve(fsos.size());

	for (const auto& vec : fsos)
	{
		eigenfunctions.emplace_back(createEigenfunction(edge1Values, edge2Values, edge3Values, vec));
	}

	return eigenfunctions;
}

QuantumGraph::StateType QuantumGraph::getSolutionValuesAtPI(size_t edge, double lambda, const StateType& initConditions, double error, size_t baseNumPoints) const
{
	double initialStep = PI / (baseNumPoints - 1);
	StateType x = initConditions;
	auto system = createSturmLiouvilleODE(m_potentials[edge - 1], lambda);
	auto stepper = boost::numeric::odeint::make_dense_output(error, error, boost::numeric::odeint::runge_kutta_dopri5<StateType>());
	size_t steps = boost::numeric::odeint::integrate_const(stepper, system, x, 0., PI, initialStep);
	return x;
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

QuantumGraph::GraphFunction QuantumGraph::createEigenfunction(const std::vector<double>& edge1Values, const std::vector<double>& edge2Values, const std::vector<double>& edge3Values, const std::array<double, 3>& fsos) const
{
	return[y = std::array{ edge1Values, edge2Values, edge3Values }, fsos]
	(std::array<double, 3> x)->std::array<double, 3> {
		constexpr double PI_INV = 1.0 / 3.14159265358979323846;

		size_t size = std::min({ y[0].size(), y[1].size(), y[2].size() });

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

std::vector<std::array<double, 3>> QuantumGraph::getNullSpaceBasis(const Matrix33& matrix, double tolerance) noexcept
{
	std::vector<std::array<double, 3>> basis;
	auto [S, V, D] = mathter::DecomposeSVD(matrix);
	for (size_t i = 0; i != 3; ++i)
	{
		if (abs(V(i, i)) <= tolerance)
		{
			basis.push_back({ D(i, 0), D(i, 1), D(i, 2) });
		}
	}
	return basis;
}

double QuantumGraph::characteristicDeterminant(double lambda, double error) const
{
	// Начальная погрешность для вычисления C1, S2, S3
	double initialError = error / 9.0;

	// Первое вычисление значений с начальной погрешностью
	StateType sol1 = getSolutionValueAtPI(1, lambda, initialError);
	StateType sol2 = getSolutionValueAtPI(2, lambda, initialError);
	StateType sol3 = getSolutionValueAtPI(3, lambda, initialError);

	// Вычисление определителя
	double result = sol1[1] * sol2[0] * sol3[0]
		+ sol1[0] * sol2[1] * sol3[0]
		+ sol1[0] * sol2[0] * sol3[1];

	// Оценка погрешности на основе частных производных
	double derivativeSol11 = std::abs(sol2[0] * sol3[0]);
	double derivativeSol21 = std::abs(sol1[0] * sol3[0]);
	double derivativeSol31 = std::abs(sol1[0] * sol2[0]);

	double derivativeSol10 = std::abs(sol2[1] * sol3[0] + sol2[0] * sol3[1]);
	double derivativeSol20 = std::abs(sol1[1] * sol3[0] + sol1[0] * sol3[1]);
	double derivativeSol30 = std::abs(sol1[1] * sol2[0] + sol1[0] * sol2[1]);

	double derivativeSum = derivativeSol11 + derivativeSol21 + derivativeSol31
		+ derivativeSol20 + derivativeSol30 + derivativeSol10;

	// Корректировка погрешности, если оценка превышает заданную погрешность
	if (3.0 * derivativeSum * initialError > error)
	{
		double correctedError = error / (3.0 * derivativeSum);

		// Повторное вычисление значений с скорректированной погрешностью
		sol1 = getSolutionValueAtPI(1, lambda, correctedError);
		sol2 = getSolutionValueAtPI(2, lambda, correctedError);
		sol3 = getSolutionValueAtPI(3, lambda, correctedError);

		// Повторное вычисление определителя
		result = sol1[1] * sol2[0] * sol3[0]
			+ sol1[0] * sol2[1] * sol3[0]
			+ sol1[0] * sol2[0] * sol3[1];
	}

	return result;
}