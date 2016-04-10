#include <iostream>
#include <iomanip>
#include <utility_functions.hpp>

int main()
{
	std::cout << "        " << std::left << std::setw(15) << "Difference"
			  << "Result\n";

	bool result = true;
	log_exp_base = 10;

	// double -> double
	double temp_d = type_caster<double, double>(100);
	std::cout << "Result: " << std::left << std::setw(15) << temp_d - 100 << (result &= (AlmostEqualRelative(temp_d, 100.0))) << '\n';

	// int -> double
	temp_d = type_caster<int, double>(2);
	std::cout << "Result: " << std::left << std::setw(15) << temp_d - 100 << (result &= (AlmostEqualRelative(temp_d, 100.0))) << '\n';

	// double -> int
	int temp_i = type_caster<double, int>(100);
	std::cout << "Result: " << std::left << std::setw(15) << temp_i - 2 << (result &= (temp_i == 2)) << '\n';

	// int -> int
	temp_i = type_caster<int, int>(100);
	std::cout << "Result: " << std::left << std::setw(15) << temp_i - 100 << (result &= (temp_i == 100)) << '\n';

	std::cout << "Final result: " << result << '\n';

	return !result;
}