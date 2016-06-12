#include <vector>

#include "dna_array.hpp"
#include "hmmalign.hpp"

int main()
{
	typedef double fp_type;

	std::vector<dna_array<fp_type, 5>> E_p{
		// A     C     G     T
		{ 1.00, 0.00, 0.00, 0.00 }, // 0
		{ 1.00, 0.00, 0.00, 0.00 }, // 1
		{ 0.40, 0.00, 0.25, 0.35 }, // 2
		{ 0.00, 0.38, 0.62, 0.00 }, // 3
	};
	std::vector<fp_type> M_D_p{
		0.00,
		0.00,
		0.00
	};
	std::vector<fp_type> D_D_p{
		0.00,
		0.00,
		0.00
	};

	fp_type substitution_rate = 0.005;

	fp_type gap_open = 0.005;

	fp_type insert_open = 0.005;
	fp_type insert_extend = 0.02;

	fp_type left_clip_open = 0.01;
	fp_type left_clip_extend = 0.90;
	fp_type right_clip_open = 0.01;
	fp_type right_clip_extend = 0.90;

	fp_type gap_extend = gap_open / (1 - insert_open - 1.0 / 12 - right_clip_open);

	reference_genome<double> parameters;
	parameters.set_parameters(
		E_p,
		M_D_p,
		D_D_p,
		background_rates{
			substitution_rate,
			gap_open,
			gap_extend,
			insert_open,
			insert_extend,
			1.0 / 12,
			left_clip_open,
			left_clip_extend,
			right_clip_open,
			right_clip_extend });

	std::vector<std::string> string_vector_1base{
		"A", "C", "G", "T"
	};

	std::vector<std::string> string_vector_2base{
		"AA", "AC", "AG", "AT",
		"CA", "CC", "CG", "CT",
		"GA", "GC", "GG", "GT",
		"TA", "TC", "TG", "TT"
	};

	std::vector<std::string> string_vector_3base{
		"AAA", "AAC", "AAG", "AAT",
		"ACA", "ACC", "ACG", "ACT",
		"AGA", "AGC", "AGG", "AGT",
		"ATA", "ATC", "ATG", "ATT",

		"CAA", "CAC", "CAG", "CAT",
		"CCA", "CCC", "CCG", "CCT",
		"CGA", "CGC", "CGG", "CGT",
		"CTA", "CTC", "CTG", "CTT",

		"GAA", "GAC", "GAG", "GAT",
		"GCA", "GCC", "GCG", "GCT",
		"GGA", "GGC", "GGG", "GGT",
		"GTA", "GTC", "GTG", "GTT",

		"TAA", "TAC", "TAG", "TAT",
		"TCA", "TCC", "TCG", "TCT",
		"TGA", "TGC", "TGG", "TGT",
		"TTA", "TTC", "TTG", "TTT"
	};

	std::vector<std::string> string_vector_4base{
		"AAAA", "AAAC", "AAAG", "AAAT",
		"AACA", "AACC", "AACG", "AACT",
		"AAGA", "AAGC", "AAGG", "AAGT",
		"AATA", "AATC", "AATG", "AATT",
		"ACAA", "ACAC", "ACAG", "ACAT",
		"ACCA", "ACCC", "ACCG", "ACCT",
		"ACGA", "ACGC", "ACGG", "ACGT",
		"ACTA", "ACTC", "ACTG", "ACTT",
		"AGAA", "AGAC", "AGAG", "AGAT",
		"AGCA", "AGCC", "AGCG", "AGCT",
		"AGGA", "AGGC", "AGGG", "AGGT",
		"AGTA", "AGTC", "AGTG", "AGTT",
		"ATAA", "ATAC", "ATAG", "ATAT",
		"ATCA", "ATCC", "ATCG", "ATCT",
		"ATGA", "ATGC", "ATGG", "ATGT",
		"ATTA", "ATTC", "ATTG", "ATTT",

		"CAAA", "CAAC", "CAAG", "CAAT",
		"CACA", "CACC", "CACG", "CACT",
		"CAGA", "CAGC", "CAGG", "CAGT",
		"CATA", "CATC", "CATG", "CATT",
		"CCAA", "CCAC", "CCAG", "CCAT",
		"CCCA", "CCCC", "CCCG", "CCCT",
		"CCGA", "CCGC", "CCGG", "CCGT",
		"CCTA", "CCTC", "CCTG", "CCTT",
		"CGAA", "CGAC", "CGAG", "CGAT",
		"CGCA", "CGCC", "CGCG", "CGCT",
		"CGGA", "CGGC", "CGGG", "CGGT",
		"CGTA", "CGTC", "CGTG", "CGTT",
		"CTAA", "CTAC", "CTAG", "CTAT",
		"CTCA", "CTCC", "CTCG", "CTCT",
		"CTGA", "CTGC", "CTGG", "CTGT",
		"CTTA", "CTTC", "CTTG", "CTTT",

		"GAAA", "GAAC", "GAAG", "GAAT",
		"GACA", "GACC", "GACG", "GACT",
		"GAGA", "GAGC", "GAGG", "GAGT",
		"GATA", "GATC", "GATG", "GATT",
		"GCAA", "GCAC", "GCAG", "GCAT",
		"GCCA", "GCCC", "GCCG", "GCCT",
		"GCGA", "GCGC", "GCGG", "GCGT",
		"GCTA", "GCTC", "GCTG", "GCTT",
		"GGAA", "GGAC", "GGAG", "GGAT",
		"GGCA", "GGCC", "GGCG", "GGCT",
		"GGGA", "GGGC", "GGGG", "GGGT",
		"GGTA", "GGTC", "GGTG", "GGTT",
		"GTAA", "GTAC", "GTAG", "GTAT",
		"GTCA", "GTCC", "GTCG", "GTCT",
		"GTGA", "GTGC", "GTGG", "GTGT",
		"GTTA", "GTTC", "GTTG", "GTTT",

		"TAAA", "TAAC", "TAAG", "TAAT",
		"TACA", "TACC", "TACG", "TACT",
		"TAGA", "TAGC", "TAGG", "TAGT",
		"TATA", "TATC", "TATG", "TATT",
		"TCAA", "TCAC", "TCAG", "TCAT",
		"TCCA", "TCCC", "TCCG", "TCCT",
		"TCGA", "TCGC", "TCGG", "TCGT",
		"TCTA", "TCTC", "TCTG", "TCTT",
		"TGAA", "TGAC", "TGAG", "TGAT",
		"TGCA", "TGCC", "TGCG", "TGCT",
		"TGGA", "TGGC", "TGGG", "TGGT",
		"TGTA", "TGTC", "TGTG", "TGTT",
		"TTAA", "TTAC", "TTAG", "TTAT",
		"TTCA", "TTCC", "TTCG", "TTCT",
		"TTGA", "TTGC", "TTGG", "TTGT",
		"TTTA", "TTTC", "TTTG", "TTTT"
	};

	std::cout << std::fixed << std::setprecision(7);
	fp_type sum = 0, result;

	hmmalign<double> my_hmm;

	for (const auto& i : { string_vector_1base, string_vector_2base, string_vector_3base, string_vector_4base })
	{
		for (const auto& j : i)
		{
			result = my_hmm.Lik(parameters, j);
			std::cout << j << ":\t" << result << '\n';
			sum += result;
		}
		std::cout << '\n';
	}
	std::cout << "Sum:\t" << sum << '\n';

	return 0;
}