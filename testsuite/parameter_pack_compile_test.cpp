#include <parameter_pack.hpp>

const int col_width = 9;

int main()
{
	typedef double fp_type;

	std::vector<dna_array<fp_type, 5>> E_p{
		// A     C     G     T
		{ 1.00, 0.00, 0.00, 0.00 }, // 0
		{ 1.00, 0.00, 0.00, 0.00 }, // 1
		{ 0.40, 0.00, 0.25, 0.35 }, // 2
		{ 0.00, 0.38, 0.62, 0.00 }, // 3
		{ 0.00, 0.00, 1.00, 0.00 }, // 4
		{ 0.00, 0.00, 1.00, 0.00 }, // 5
		{ 0.53, 0.00, 0.00, 0.47 }, // 6
		{ 0.65, 0.35, 0.00, 0.00 }, // 7
		{ 0.65, 0.35, 0.00, 0.00 }, // 8
	};
	std::vector<fp_type> M_D_p{
		0.00, // 0 -> 1
		0.00, // 1 -> 2
		0.35, // 2 -> 3
		0.38, // 3 -> 4
		0.00, // 4 -> 5
		0.00, // 5 -> 6
		0.00, // 6 -> 7
		0.00 // 7 -> 8
	};
	std::vector<fp_type> D_D_p{
		0.00, // 0 -> 1
		0.00, // 1 -> 2
		0.00, // 2 -> 3
		0.00, // 3 -> 4
		0.98, // 4 -> 5
		0.42, // 5 -> 6
		0.00, // 6 -> 7
		0.00 // 7 -> 8
	};

	fp_type substitution_rate = 0.005;

	fp_type gap_open = 0.005;
	fp_type gap_extend = 0.19;

	fp_type insert_open = 0.17;
	fp_type insert_extend = 0.18;

	fp_type left_clip_open = 0.01;
	fp_type left_clip_extend = 0.90;
	fp_type right_clip_open = 0.01;
	fp_type right_clip_extend = 0.90;

	parameter_pack<fp_type> fp_pp;
	fp_pp.set_parameters(
		E_p,
		M_D_p,
		D_D_p,
		background_rates<double>{
			substitution_rate,
			gap_open,
			gap_extend,
			insert_open,
			insert_extend,
			1.0 / 11,
			left_clip_open,
			left_clip_extend,
			right_clip_open,
			right_clip_extend });
	fp_pp.display_parameters();

	parameter_pack<int32_t> int_pp;
	int_pp.set_parameters(
		E_p,
		M_D_p,
		D_D_p,
		background_rates<double>{
			substitution_rate,
			gap_open,
			gap_extend,
			insert_open,
			insert_extend,
			1.0 / 11,
			left_clip_open,
			left_clip_extend,
			right_clip_open,
			right_clip_extend });
	int_pp.display_parameters();

	parameter_pack<fp_type> fp_pp_convert1(int_pp);
	parameter_pack<int32_t> int_pp_convert2(fp_pp_convert1);

	std::cout << "Displaying parameters for doubly converted parameter pack:\n";
	fp_pp_convert1.display_parameters(false);
	int_pp_convert2.display_parameters(false);
	return 0;
}