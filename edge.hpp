#include <boost/numeric/odeint.hpp>
#include <functional>
#include <array>
#include <utility>
#include <vector>
#include <algorithm>
#include <numeric>

class Edge
{
	std::function<double(double)> _potential = nullptr;

	double _leftBound = 0.;

	double _rightBound = 3.141592653589793;

public:
	Edge(std::function<double(double)> potential, double leftBound = 0., double rightBound = 3.141592653589793) noexcept :
		_potential(potential), _leftBound(leftBound), _rightBound(rightBound)
	{}

	double getLeftBound() const noexcept
	{
		return _leftBound;
	}

	double getRightBound() const noexcept
	{
		return _rightBound;
	}

	std::function<double(double)> getPotential() const noexcept
	{

	}
};