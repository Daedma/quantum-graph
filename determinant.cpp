#include "QuantumGraph.hpp"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/math/constants/constants.hpp>

static constexpr double PI = boost::math::double_constants::pi;
using state_type = std::array<double, 2>;


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

double getAnalyticalDeterminant(double lambda)
{
	state_type S = getAnaliticalS(lambda);
	state_type C = getAnaliticalC(lambda);
	return
		C[1] * S[0] * S[0] +
		C[0] * S[1] * S[0] +
		C[0] * S[0] * S[1];
}

void test_determinant(double lower, double upper, double step)
{
	QuantumGraph graph{
		[](double) { return 0.; },
		[](double) { return 0.; },
		[](double) { return 0.; }
	};

	double maxDiscrepancy = 0;
	for (double error = 1.e-10; error <= 1.e-2; error *= 10)
	{
		std::cout << "Error: " << error << '\n';
		for (double lambda = lower; lambda <= upper; lambda += step)
		{
			double analytical = getAnalyticalDeterminant(lambda);
			double numerical = graph.characteristicDeterminant(lambda, error);
			double difference = std::abs(analytical - numerical);

			double disperancy = difference / error;
			if (disperancy > maxDiscrepancy)
			{
				maxDiscrepancy = disperancy;
				std::cout << "\tMax disperancy: " << maxDiscrepancy << '\n';
			}
		}
	}
	std::cout << "Final max disperancy: " << maxDiscrepancy << '\n';
}

int main()
{
	test_determinant(0, 10000, 0.1);
	return 0;
}
