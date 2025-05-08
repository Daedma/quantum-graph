#include <fstream>
#include <iostream>
#include <iomanip>
#include <random>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "QuantumGraph.hpp"

auto gen_potencial_sin(double a) noexcept
{
	using boost::math::double_constants::pi;
	return [a](double x) { return a * pi * sin(x); };
}

QuantumGraph gen_graph(double a1, double a2, double a3) noexcept
{
	return QuantumGraph{
		gen_potencial_sin(a1),
		gen_potencial_sin(a2),
		gen_potencial_sin(a3)
	};
}

double get_null_ev(int n, int k) noexcept
{
	using boost::math::double_constants::one_div_pi;
	static const double main_term = one_div_pi * std::asin(std::sqrt(2. / 3.));
	double ev;
	switch (k)
	{
	case 1:
		ev = main_term + n;
		break;
	case 2:
		ev = -main_term + n + 1;
		break;
	case 3:
		ev = n + 1;
		break;
	default:
		ev = NAN;
		break;
	}
	return ev * ev;
}


int main(int argc, char const* argv[])
{
	constexpr double lowerBound = 24600.;
	constexpr double upperBound = 25100.;
	constexpr double step = 50.;
	constexpr size_t depth = 25;
	constexpr size_t numberOfNobes = 6000;

	if (argc != 2)
	{
		std::cerr << "Usage: " << argv[0] << " <number of graphs>\n";
		return 1;
	}
	int num_of_graphs = std::atoi(argv[1]);

	const std::string filename = "ev3-" + std::to_string(static_cast<int>(lowerBound)) + "-" +
		std::to_string(static_cast<int>(upperBound)) + "-" + std::to_string(static_cast<int>(step)) + "-" +
		std::to_string(depth) + ".csv";

	std::ofstream ofs{ filename };
	ofs << std::setprecision(17);
	ofs << "a1, a2, a3, eigenvalues,\n";
	ofs << "0., 0., 0.,";
	for (int i = 0;;++i)
	{
		double ev = get_null_ev(i / 3, i % 3 + 1);
		assert(!std::isnan(ev));
		if (ev >= lowerBound && ev <= upperBound)
		{
			ofs << ev << ',';
		}
		if (ev > upperBound)
		{
			break;
		}
	}
	ofs << '\n';
	ofs.close();

	std::uniform_real_distribution<double> dist{ -5., 5. };
	std::mt19937 gen{ std::random_device{}() };

	for (int i = 0; i < num_of_graphs; ++i)
	{
		auto a1 = dist(gen);
		auto a2 = dist(gen);
		auto a3 = dist(gen);
		auto graph = gen_graph(a1, a2, a3);
		graph.setNumberOfNodes(numberOfNobes);
		auto ev = graph.calcEigenvalues(lowerBound, upperBound, step, depth);
		std::ofstream ofs{ filename, std::ios::app };
		ofs << std::setprecision(17);
		ofs << a1 << ',' << a2 << ',' << a3 << ',';
		for (auto e : ev)
			ofs << e << ',';
		ofs << '\n';
		ofs.close();
	}

	std::clog << "Program finished\n";
	return 0;
}
