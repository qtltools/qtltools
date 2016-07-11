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

#include "trans_data.h"

void trans_main(vector < string > & argv) {
	trans_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();  //Mandatory
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF format.")
		("bed", boost::program_options::value< string >(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< string >(), "Covariates in TXT format.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("normal", "Quantile normalize phenotype data.")
		("window", boost::program_options::value< double >()->default_value(5e6, "5e6"), "Cis-window of variants to be excluded.")
		("threshold", boost::program_options::value< double >()->default_value(1e-5, "1e-5"), "P-value threshold below which hits are reported.")
		("bins", boost::program_options::value< unsigned int >()->default_value(1000), "Number of bins to use to categorize all p-values above --threshold.");

	boost::program_options::options_description opt_modes ("\x1B[32mAnalysis type\33[0m");
	opt_modes.add_options()
		("nominal", "MODE1: NOMINAL PASS [Pvalues are not adjusted].")
		("adjust", boost::program_options::value< string >(), "MODE2: ADJUSTED PASS [Pvalues are adjusted].")
		("permute", "MODE3: PERMUTATION PASS [Permute all phenotypes once].")
		("sample", boost::program_options::value< unsigned int >(), "MODE4: PERMUTATION PASS [Permute randomly chosen phenotypes multiple times].");
    
    boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
    opt_parallel.add_options()
        ("chunk", boost::program_options::value< vector < int > >()->multitoken(), "Specify which chunk needs to be processed");

	D.option_descriptions.add(opt_files).add(opt_parameters).add(opt_modes).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [trans] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("MAPPING QTL IN TRANS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("vcf")) vrb.error("Genotype data needs to be specified with --vcf [file.vcf]");
	if (!D.options.count("bed")) vrb.error("Phenotype data needs to be specified with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");
	int nMode = D.options.count("nominal") + D.options.count("adjust") + D.options.count("permute") + D.options.count("sample");
	if (nMode != 1) vrb.error("Please, specify only one of these options [--nominal, --adjust, --permute, --sample]");
    if (D.options.count("chunk") && D.options.count("sample")) vrb.error("--chunk cannot be combined with --sample");

	//---------
	// 5. MODES
	//---------

	//MODE1: NOMINAL PASS NON ADJUSTED
	if (D.options.count("nominal")) {
		D.mode = TRANS_MODE1;
		vrb.bullet("TASK: Perform a full nominal pass, do not adjust p-values");
	}

	//MODE2: NOMINAL PASS ADJUSTED
	if (D.options.count("adjust")) {
		D.mode = TRANS_MODE2;
		vrb.bullet("TASK: Test and adjust p-values using [" + D.options["adjust"].as < string > () +"]");
	}

	//MODE3: PERMUTATION PASS
	if (D.options.count("permute")) {
		D.mode = TRANS_MODE3;
		vrb.bullet("TASK: Permute all phenotype once and test");
	}

	//MODE4: PERMUTATION PASS
	if (D.options.count("sample")) {
		D.mode = TRANS_MODE4;
		vrb.bullet("TASK: Permute randomly chosen phenotypes " + stb.str(D.options["sample"].as < unsigned int > ()) + " times and test");
	}

	//--------------
	// 6. SET PARAMS
	//--------------
	if (D.options["window"].as < double > () <= 0 || D.options["window"].as < double > () > 1e9) vrb.error ("Incorrect cis-window size");
	vrb.bullet("Cis-window size is " + stb.str((int)D.options["window"].as < double > ()) + " bp");
	D.cis_window = D.options["window"].as < double > ();
    if (D.options.count("chunk")) {
        vector < int > nChunk = D.options["chunk"].as < vector < int > > ();
        if (nChunk.size() != 2 || nChunk[0] > nChunk[1]) vrb.error("Incorrect --chunk arguments!");
        vrb.bullet("Chunk = [" + stb.str(nChunk[0]) + "/" + stb.str(nChunk[1]) + "]");
    }
	D.n_bins = D.options["bins"].as < unsigned int > ();
	vrb.bullet("#bins = " + stb.str(D.n_bins));
	vrb.bullet("threshold = " + stb.str(D.options["threshold"].as < double > ()));
    
    
    D.processBasicOptions();
    if (D.options.count("chunk")) {
        D.scanPhenotypes(D.options["bed"].as < string > ());
        D.setPhenotypeLines(D.options["chunk"].as < vector < int > > ()[0], D.options["chunk"].as < vector < int > > ()[1]);
        D.clear();
    }

	//---------------------------
	// 7. READ FILES & INITIALIZE
	//---------------------------
    //D.processBasicOptions();
	D.readSampleFromVCF(D.options["vcf"].as < string > ());
	D.checkSampleInBED(D.options["bed"].as < string > ());
	if (D.options.count("cov")) D.checkSampleInCOV(D.options["cov"].as < string > ());

	D.readPhenotypes(D.options["bed"].as < string > ());
	D.imputePhenotypes();

	if (D.options.count("cov")) {
		D.readCovariates(D.options["cov"].as < string > ());
		D.residualizePhenotypes();
	}
	if (D.options.count("normal")) D.normalTranformPhenotypes();

	if (D.options.count("permute")) D.shufflePhenotypes();
	if (D.options.count("sample")) D.samplePhenotypes(D.options["sample"].as < unsigned int > ());
	if (D.options.count("adjust")) D.buildNullDistribution(D.options["adjust"].as < string > ());
	D.getCorrelationThreshold(D.options["threshold"].as < double > ());

	D.normalizePhenotypes();

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.runTransPass(D.options["vcf"].as < string > (), D.options["out"].as < string > ());
}
