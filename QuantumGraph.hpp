#pragma once
#include <functional>
#include <array>
#include <vector>
#include <Mathter/Matrix.hpp>

class QuantumGraph
{
	static constexpr double PI = 3.141592653589793;

	using StateType = std::array<double, 2>;

	using Matrix33 = mathter::Matrix<double, 3, 3, mathter::eMatrixOrder::PRECEDE_VECTOR>;

private:
	size_t numberOfNodes = 1000;

	std::array<std::function<double(double)>, 3> m_potentials;

public:
	using GraphFunction = std::function<std::array<double, 3>(const std::array<double, 3>&)>;

	QuantumGraph(std::function<double(double)> q1, std::function<double(double)> q2, std::function<double(double)> q3) noexcept :
		m_potentials({ q1, q2, q3 })
	{}

	void setNumberOfNodes(size_t num) noexcept { numberOfNodes = num; }

	double getNumberOfNodes() const noexcept { return numberOfNodes; }

	std::vector<double> calcEigenvalues(double lowerBound, double higherBound, double step = 0.1, size_t binarySearchDepth = 30) const;

	GraphFunction calcEigenfunction(double lambda, double tolerance) const;

	GraphFunction calcEigenfunction(double lambda) const
	{
		return calcEigenfunction(lambda, abs(characteristicDeterminant(lambda) * 100.));
	}

private:
	StateType getCosValueAtPI(size_t edge, double lambda) const
	{
		return getSolutionValuesAtPI(edge, lambda, { 1, 0 }, numberOfNodes);
	}

	StateType getSinValueAtPI(size_t edge, double lambda) const
	{
		return getSolutionValuesAtPI(edge, lambda, { 0, 1 }, numberOfNodes);
	}

	StateType getSolutionValuesAtPI(size_t edge, double lambda, const StateType& initConditions, size_t numPoints) const;

	static auto createSturmLiouvilleODE(std::function<double(double)> q, double lambda) noexcept
	{
		return [q, lambda](const StateType& x, StateType& dxdt, double t) {
			dxdt[0] = x[1];
			dxdt[1] = q(t) * x[0] - lambda * x[0];
			};
	}

	double characteristicDeterminant(double lambda) const
	{
		StateType C1 = getCosValueAtPI(1, lambda);
		StateType S2 = getSinValueAtPI(2, lambda);
		StateType S3 = getSinValueAtPI(3, lambda);

		return C1[1] * S2[0] * S3[0]
			+ C1[0] * S2[1] * S3[0]
			+ C1[0] * S2[0] * S3[1];
	}

	Matrix33 getSystemMatrix(double lambda) const
	{
		StateType C1 = getCosValueAtPI(1, lambda);
		StateType S2 = getSinValueAtPI(2, lambda);
		StateType S3 = getSinValueAtPI(3, lambda);

		return {
			C1[0], -S2[0], 0.,
			0., S2[0], -S3[0],
			C1[1], S2[1], S3[1]
		};
	}

	std::vector<double> getSolutionValues(size_t edge, double lambda, const StateType& initConditions, size_t numPoints) const;

	static std::array<double, 3> getNullSpaceBasis(const Matrix33& matrix, double tolerance = 1.e-10) noexcept;
};
