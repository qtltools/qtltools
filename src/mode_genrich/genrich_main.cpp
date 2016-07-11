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

#include "genrich_data.h"

void genrich_main(vector < string > & argv) {
	genrich_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("tss", boost::program_options::value< string >(), "Phenotype file used for the QTL mapping.")
		("qtl", boost::program_options::value< string >(), "List of QTLs in BED format.")
		("ref", boost::program_options::value< string >(), "1000 Genomes genotypes in VCF/BCF format.")
		("gwas", boost::program_options::value< string >(), "List of GWAS hits in BED format.")
		("out", boost::program_options::value< string >(), "Output filename.");

	boost::program_options::options_description opt_param ("\x1B[32mParameters\33[0m");
	opt_param.add_options()
		("threshold-maf", boost::program_options::value< double >()->default_value(0.01), "MAF filter for sites in 1000 Genomes.")
		("threshold-ld", boost::program_options::value< double >()->default_value(0.5), "Consider that a GWAS hit and a QTL belong to the same signal when r2 >= arg.")
		("bin-distance", boost::program_options::value < unsigned int > ()->default_value(5000), "Maximal distance to assume 2 QTL are within the same bin.")
		("bin-maf", boost::program_options::value < double > ()->default_value(0.02), "Maximal frequency difference to assume 2 QTL are within the same bin.")
		("permute", boost::program_options::value < unsigned int > ()->default_value(1000), "Number of null sets of variants to be samples.");

	D.option_descriptions.add(opt_files).add(opt_param);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [genrich] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("COMPUTING GWAS ENRICHMENT FOR QTLs");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("ref")) vrb.error("Please specify a variant collection with --ref [file.vcf]");
	if (!D.options.count("qtl")) vrb.error("Please specify a QTL list with --qtl [qtl.bed]");
	if (!D.options.count("gwas")) vrb.error("Please specify a GWAS hits list with --gwas [hits.bed]");
	if (!D.options.count("tss")) vrb.error("Please specify a list of reference positions with --tss [phenotype.txt]");
	if (!D.options.count("out")) vrb.error("Please specify an output file with --out [out.txt]");

	//------------------
	// 5. SET PARAMETERS
	//------------------
	D.threshold_ld = D.options["threshold-ld"].as < double > ();
	D.threshold_maf = D.options["threshold-maf"].as < double > ();
	D.bin_distance = D.options["bin-distance"].as < unsigned int > ();
	D.bin_maf = D.options["bin-maf"].as < double > ();
	D.n_permutations = D.options["permute"].as < unsigned int > ();
	//TO BE DONE: check parameter values within reasonable range.
	vrb.bullet("LD threshold = " + stb.str(D.threshold_ld));
	vrb.bullet("MAF threshold = " + stb.str(D.threshold_maf));
	vrb.bullet("Distance binning = " + stb.str(D.bin_distance));
	vrb.bullet("MAF binning = " + stb.str(D.bin_maf));
	vrb.bullet("#permutations = " + stb.str(D.n_permutations));

	//--------------
	// 6. READ FILES
	//--------------
	D.processBasicOptions();
	D.readSampleFromVCF(D.options["ref"].as < string > ());
	D.mergeSampleLists();
	D.readPhenotypes(D.options["tss"].as < string > ());

	D.readReferenceGenotypes(D.options["ref"].as < string > ());

	D.readQTL(D.options["qtl"].as < string > ());
	D.readGWAS(D.options["gwas"].as < string > ());

	//----------------
	// 7. RUN ANALYSIS
	//----------------
	D.binningAllVariants();
	D.overlapGWASandQTL(D.options["out"].as < string > ());
}
