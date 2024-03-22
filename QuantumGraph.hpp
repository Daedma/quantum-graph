#pragma once
#include <functional>
#include <array>
#include <vector>
#include <Mathter/Matrix.hpp>

class QuantumGraph
{
	static constexpr double PI = 3.141592653589793;

	static constexpr size_t NUMBER_OF_NODES = 1000;

	using state_type = std::array<double, 2>;

	using matrix = mathter::Matrix<double, 3, 3, mathter::eMatrixOrder::PRECEDE_VECTOR>;

private:
	std::array<std::function<double(double)>, 3> m_potentials;

public:
	using graph_function = std::function<std::array<double, 3>(double, double, double)>;

	QuantumGraph(std::function<double(double)> q1, std::function<double(double)> q2, std::function<double(double)> q3) noexcept :
		m_potentials({ q1, q2, q3 })
	{}

	std::vector<double> getEigenvalues(double lowerBound, double higherBound)
	{
		constexpr double STEP = 0.1;
		constexpr size_t BINARY_SEARCH_DEPTH = 30;

		std::vector<double> eigenvalues;

		for (; lowerBound < higherBound; lowerBound += STEP)
		{
			double detLeft = characteristicDeterminant(lowerBound);
			double detRight = characteristicDeterminant(lowerBound + STEP);
			if (detRight * detLeft <= 0)
			{
				double left = lowerBound, right = lowerBound + STEP;
				for (size_t it = 0; it != BINARY_SEARCH_DEPTH; ++it)
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




// private:
	state_type getCosValueAtPi(size_t edge, double lambda) const
	{
		return getSolutionValues(edge, lambda, { 1, 0 }, NUMBER_OF_NODES);
	}

	state_type getSinValueAtPi(size_t edge, double lambda) const
	{
		return getSolutionValues(edge, lambda, { 0, 1 }, NUMBER_OF_NODES);
	}

	state_type getSolutionValues(size_t edge, double lambda, const state_type& initConditions, size_t numPoints) const;

	auto create_sturm_liouville_ode(std::function<double(double)> q, double lambda) const noexcept
	{
		return [q, lambda](const state_type& x, state_type& dxdt, double t) {
			dxdt[0] = x[1];
			dxdt[1] = q(t) * x[0] - lambda * x[0];
			};
	}

	double characteristicDeterminant(double lambda) const
	{
		state_type C1 = getCosValueAtPi(1, lambda);
		state_type S2 = getSinValueAtPi(2, lambda);
		state_type S3 = getSinValueAtPi(3, lambda);

		return C1[1] * S2[0] * S3[0]
			+ C1[0] * S2[1] * S3[0]
			+ C1[0] * S2[0] * S3[1];
	}

	matrix getSystemMatrix(double lambda) const
	{
		state_type C1 = getCosValueAtPi(1, lambda);
		state_type S2 = getSinValueAtPi(2, lambda);
		state_type S3 = getSinValueAtPi(3, lambda);
		return {
			C1[0], -S2[0], 0.,
			0., S2[0], -S3[0],
			C1[1], S2[1], S3[1]
		};
	}

	std::array<double, 3> getNullSpaceBasis(const matrix& m) const;

};
