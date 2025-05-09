#pragma once
#include <functional>
#include <array>
#include <vector>
#include <optional>
#include <Mathter/Matrix.hpp>

// TODO : Распараллелить вычисление C1, S2, S3 при расчете определителя
// TODO : Распараллелить поиск собственных значений
// TODO : Добавить обработку случая кратных собственных значений
// TODO : Добавить параметризацию граничных условий
// TODO : Подумать над адаптивным шагом при поиске собственных значений


class QuantumGraph
{
	static constexpr double PI = 3.141592653589793;

	using StateType = std::array<double, 2>;

	using Matrix33 = mathter::Matrix<double, 3, 3, mathter::eMatrixOrder::PRECEDE_VECTOR>;

private:
	std::array<std::function<double(double)>, 3> m_potentials;

public:
	using GraphFunction = std::function<std::array<double, 3>(const std::array<double, 3>&)>;

	QuantumGraph(std::function<double(double)> q1, std::function<double(double)> q2, std::function<double(double)> q3) noexcept :
		m_potentials({ q1, q2, q3 })
	{}

	std::vector<double> calcEigenvalues(double lowerBound, double higherBound, double step, double error, size_t maxIter) const;

	GraphFunction calcEigenfunction(double lambda, size_t numOfNodes, double tolerance) const;

	double characteristicDeterminant(double lambda, double error) const;

private:
	StateType getCosValueAtPI(size_t edge, double lambda, double error) const
	{
		return getSolutionValuesAtPI(edge, lambda, { 1, 0 }, error, 100);
	}

	StateType getSinValueAtPI(size_t edge, double lambda, double error) const
	{
		return getSolutionValuesAtPI(edge, lambda, { 0, 1 }, error, 100);
	}

	StateType getSolutionValuesAtPI(size_t edge, double lambda, const StateType& initConditions, double error, size_t baseNumPoints) const;

	static auto createSturmLiouvilleODE(std::function<double(double)> q, double lambda) noexcept
	{
		return [q, lambda](const StateType& x, StateType& dxdt, double t) {
			dxdt[0] = x[1];
			dxdt[1] = q(t) * x[0] - lambda * x[0];
			};
	}

	Matrix33 getSystemMatrix(double lambda) const
	{
		constexpr double error = 1.e-3;
		StateType C1 = getCosValueAtPI(1, lambda, error);
		StateType S2 = getSinValueAtPI(2, lambda, error);
		StateType S3 = getSinValueAtPI(3, lambda, error);

		return {
			C1[0], -S2[0], 0.,
			0., S2[0], -S3[0],
			C1[1], S2[1], S3[1]
		};
	}

	std::vector<double> getSolutionValues(size_t edge, double lambda, const StateType& initConditions, size_t numPoints) const;

	static std::array<double, 3> getNullSpaceBasis(const Matrix33& matrix, double tolerance = 1.e-10) noexcept;
};
