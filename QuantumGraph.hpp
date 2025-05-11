#pragma once
#include <functional>
#include <array>
#include <vector>
#include <iostream>
#include <Mathter/Matrix.hpp>

// TODO : Добавить параметризацию граничных условий

class QuantumGraph
{
	static constexpr double PI = 3.141592653589793;

	using StateType = std::array<double, 2>;

	using Matrix33 = mathter::Matrix<double, 3, 3, mathter::eMatrixOrder::PRECEDE_VECTOR>;

public:
	enum class eBoundaryCondition
	{
		Dirichlet,
		Neumann
	};

	static constexpr std::array<eBoundaryCondition, 3> dirichletOnly = { eBoundaryCondition::Dirichlet, eBoundaryCondition::Dirichlet, eBoundaryCondition::Dirichlet };

	static constexpr std::array<eBoundaryCondition, 3> neumannOnly = { eBoundaryCondition::Neumann, eBoundaryCondition::Neumann, eBoundaryCondition::Neumann };

private:
	std::array<std::function<double(double)>, 3> m_potentials;

	std::array<eBoundaryCondition, 3> m_boundaryConditions;

public:
	using GraphFunction = std::function<std::array<double, 3>(const std::array<double, 3>&)>;

	QuantumGraph(std::function<double(double)> q1, std::function<double(double)> q2, std::function<double(double)> q3,
		eBoundaryCondition bc1, eBoundaryCondition bc2, eBoundaryCondition bc3) noexcept :
		m_potentials({ q1, q2, q3 }), m_boundaryConditions({ bc1, bc2, bc3 })
	{}

	QuantumGraph(std::function<double(double)> q1, std::function<double(double)> q2, std::function<double(double)> q3) noexcept :
		QuantumGraph(q1, q2, q3, eBoundaryCondition::Neumann, eBoundaryCondition::Dirichlet, eBoundaryCondition::Dirichlet)
	{}

	std::vector<double> calcEigenvalues(double lowerBound, double higherBound, double step, double error, size_t maxIter = 1000) const;

	std::vector<double> calcEigenvalues(double lowerBound, double higherBound, size_t windowSize, double initialStep, double error, size_t maxIter = 1000) const;

	std::vector<double> calcEigenvalues(double lowerBound, double higherBound, double error = 1.e-6) const
	{
		return calcEigenvalues(lowerBound, higherBound, 2, 0.1, error, 1000);
	}

	std::vector<GraphFunction> calcEigenfunction(double lambda, size_t numOfNodes, double error) const;

	double characteristicDeterminant(double lambda, double error) const;

private:
	StateType getCosValueAtPI(size_t edge, double lambda, double error) const
	{
		return getSolutionValuesAtPI(edge, lambda, { 1, 0 }, error, 3);
	}

	StateType getSinValueAtPI(size_t edge, double lambda, double error) const
	{
		return getSolutionValuesAtPI(edge, lambda, { 0, 1 }, error, 3);
	}

	StateType getSolutionValuesAtPI(size_t edge, double lambda, const StateType& initConditions, double error, size_t baseNumPoints) const;

	StateType getSolutionValueAtPI(size_t edge, double lambda, double error) const
	{
		StateType initConditions = m_boundaryConditions[edge - 1] == eBoundaryCondition::Neumann ? StateType{ 1, 0 } : StateType{ 0, 1 };
		return getSolutionValuesAtPI(edge, lambda, initConditions, error, 3);
	}

	static auto createSturmLiouvilleODE(std::function<double(double)> q, double lambda) noexcept
	{
		return [q, lambda](const StateType& x, StateType& dxdt, double t) {
			dxdt[0] = x[1];
			dxdt[1] = q(t) * x[0] - lambda * x[0];
			};
	}

	Matrix33 getSystemMatrix(double lambda, double error) const
	{
		StateType sval1 = getSolutionValueAtPI(1, lambda, error);
		StateType sval2 = getSolutionValueAtPI(2, lambda, error);
		StateType sval3 = getSolutionValueAtPI(3, lambda, error);

		return {
			sval1[0], -sval2[0], 0.,
			0., sval2[0], -sval3[0],
			sval1[1], sval2[1], sval3[1]
		};
	}

	std::vector<double> getSolutionValues(size_t edge, double lambda, const StateType& initConditions, size_t numPoints) const;

	double characteristicDeterminantSignSafe(double lambda, double& error, bool saveError = true) const
	{
		constexpr size_t maxIter = 100;
		size_t it = 0;

		double localError = error;
		double* errorToUse = saveError ? &localError : &error;

		double result = characteristicDeterminant(lambda, *errorToUse);
		while ((result - *errorToUse) * (result + *errorToUse) < 0 && it != maxIter)
		{
			*errorToUse *= 0.5;
			result = characteristicDeterminant(lambda, *errorToUse);
			++it;
		}

		if (it == maxIter)
		{
			result = 0.;
		}
		
		return result;
	}

	double characteristicDeterminantSignSafe(double lambda, const double& error) const
	{
		double localError = error;
		return characteristicDeterminantSignSafe(lambda, localError, true);
	}

	GraphFunction createEigenfunction(const std::vector<double>& edge1Values, const std::vector<double>& edge2Values, const std::vector<double>& edge3Values, const std::array<double, 3>& fsos) const;

	static std::vector<std::array<double, 3>> getNullSpaceBasis(const Matrix33& matrix, double tolerance = 1.e-10) noexcept;
};
