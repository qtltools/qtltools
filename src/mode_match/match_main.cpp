/*Copyright (C) 2015 Olivier Delaneau, Halit Ongen, Emmanouil T. Dermitzakis
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.*/

#include "match_data.h"

void match_main(vector < string > & argv) {
	match_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF format.")
		("bam", boost::program_options::value< string >(), "Sequence data in BAM/SAM format.")
		("reg", boost::program_options::value< string >()->default_value(""), "Genomic region(s) to be processed.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mFilters\33[0m");
	opt_parameters.add_options()
		("filter-mapping-quality", boost::program_options::value< unsigned int >()->default_value(10), "Minimal phred mapping quality for a read to be considered.")
		("filter-base-quality", boost::program_options::value< unsigned int >()->default_value(5), "Minimal phred quality for a base to be considered.")
		("filter-binomial-pvalue", boost::program_options::value< double >()->default_value(0.05, "0.05"), "Binomial p-value threshold below which a het genotype is considered as exhibiting allelic imbalance.")
		("filter-minimal-coverage", boost::program_options::value< unsigned int >()->default_value(10), "Minimal coverage for a genotype to be considered.")
		("filter-imputation-qual", boost::program_options::value< double >()->default_value(0.90, "0.90"), "Minimal imputation information score for a variant to be considered.")
		("filter-imputation-prob", boost::program_options::value< double >()->default_value(0.99, "0.99"), "Minimal posterior probability for a genotype to be considered.")
		("filter-keep-duplicates", "Keep duplicate sequencing reads in the process.");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	boost::program_options::variables_map options;
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [match] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("MATCHING SAMPLES BETWEEN GENOTYPES AND SEQUENCES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("vcf")) vrb.error("Genotype data needs to be specified with --vcf [file.vcf]");
	if (!D.options.count("bam")) vrb.error("Sequence data needs to be specified with --bam [file.bam]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");

	//TO DO CHECK PARAMETER VALUES
	D.param_min_mapQ = D.options["filter-mapping-quality"].as < unsigned int > ();
	D.param_min_baseQ = D.options["filter-base-quality"].as < unsigned int > ();
	D.param_min_cov = D.options["filter-minimal-coverage"].as < unsigned int > ();
	D.param_min_pval = D.options["filter-binomial-pvalue"].as < double > ();
	D.param_min_gp = D.options["filter-imputation-prob"].as < double > ();
	D.param_min_iq = D.options["filter-imputation-qual"].as < double > ();
	D.param_dup_rd = (D.options.count("filter-keep-duplicates") != 0);
	vrb.bullet("Mapping quality >= " + stb.str(D.param_min_mapQ));
	vrb.bullet("Base quality >= " + stb.str(D.param_min_baseQ));
	vrb.bullet("Coverage >= " + stb.str(D.param_min_cov));
	vrb.bullet("Binomial p-value threshold = " + stb.str(D.param_min_pval));
	vrb.bullet("Genotype probability >= " + stb.str(D.param_min_gp));
	vrb.bullet("Imputation quality >= " + stb.str(D.param_min_iq));
	vrb.bullet("Keep duplicate reads = " + stb.str(D.param_dup_rd));

	//------------------------------------------
	// 5. READ FILES / INITIALIZE / RUN ANALYSIS
	//------------------------------------------

	D.processBasicOptions();
	D.readSampleFromVCF(D.options["vcf"].as < string > ());
	D.mergeSampleLists();
	D.readGenotypes(D.options["vcf"].as < string > (), D.options["reg"].as < string > ());
	D.readSequences(D.options["bam"].as < string > ());
	D.writeOutput(D.options["out"].as < string > ());
}
