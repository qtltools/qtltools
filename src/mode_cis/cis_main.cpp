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

#define _DECLARE_TOOLBOX_HERE
#include "cis_data.h"

void cis_main(vector < string > & argv) {
	cis_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF/BED format.")
		("bed", boost::program_options::value< string >(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< string >(), "Covariates in TXT format.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("normal", "Normal transform the phenotypes.")
		("window", boost::program_options::value< unsigned int >()->default_value(1000000), "Size of the cis-window.");

	boost::program_options::options_description opt_modes ("\x1B[32mAnalysis type\33[0m");
	opt_modes.add_options()
		("permute", boost::program_options::value< int >(), "MODE1: PERMUTATION PASS.")
		("nominal", boost::program_options::value< string >(), "MODE2: NOMINAL PASS.")
		("mapping", boost::program_options::value< string >(), "MODE3: MAPPING PASS.");

	boost::program_options::options_description opt_aggr ("\x1B[32mPhenotype aggregation methods\33[0m");
	opt_aggr.add_options()
		("grp-best", "Correct for multiple phenotypes within a group.")
		("grp-pca1", "Run PCA on phenotypes within a group and use PC1 for association testing.")
		("grp-mean", "Average phenotypes within a group and use the results for association testing.");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("chunk", boost::program_options::value< vector < int > >()->multitoken(), "Specify which chunk needs to be processed")
		("region", boost::program_options::value< string >(), "Region of interest.")
		("region-pair", boost::program_options::value< vector < string > >()->multitoken(), "Pairs of interest.");

	D.option_descriptions.add(opt_files).add(opt_parameters).add(opt_modes).add(opt_aggr).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [cis] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("MAPPING QTL IN CIS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------BEST
	if (!D.options.count("vcf")) vrb.error("Genotype data needs to be specified with --vcf [file.vcf]");
	if (!D.options.count("bed")) vrb.error("Phenotype data needs to be specified with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");
	int nParallel = D.options.count("chunk") + D.options.count("region") + D.options.count("region-pair");
	if (nParallel > 1) vrb.error("Please, specify one of these options [--region, --chunk, --region-pair]");
	int nMode = D.options.count("permute") + D.options.count("mapping") + D.options.count("nominal");
	if (nMode != 1) vrb.error("Please, specify only one of these options [--permute, --nominal, --mapping]");
	string outFile = D.options["out"].as < string > ();

	//---------
	// 5. MODES
	//---------

	//ONLY MODE1: PERMUTATION PASS
	if (D.options.count("permute")) {
		D.mode = CIS_PERM;
		if (D.options["permute"].as < int >() < 10) vrb.error("Number of permutations is incorrect.");
		vrb.bullet("TASK: Perform " + stb.str(D.options["permute"].as < int >()) + " permutations");
		D.n_permutations = D.options["permute"].as < int >();
	}
	//ONLY MODE2: NOMINAL PASS
	if (D.options.count("nominal")) {
		D.mode = CIS_NOMI;
		string argument = D.options["nominal"].as < string >();
		if (stb.numeric(argument)) {
			D.threshold = atof(argument.c_str());
			if (D.threshold  <= 0 ||  D.threshold> 1.0) vrb.error("Significance threshold is outside of the range ]0,1]");
			else vrb.bullet("TASK: Report all nominal associations with p <= " + stb.str(D.threshold));
		} else vrb.bullet("TASK: Report all nominal associations with p below thresholds in [" + argument +"]");
	}
	//ONLY MODE3: MAPPING PASS
	if (D.options.count("mapping")) {
		D.mode = CIS_COND;
		vrb.bullet("TASK: Perform a conditional pass");
	}

	//--------------
	// 6. SET PARAMS
	//--------------
	if (D.options.count("cov")) vrb.bullet("Linear model: Phe ~ Var + Cov");
	else vrb.bullet("Linear model: Phe ~ Var");
	if (D.options["window"].as < unsigned int > () > 1000000000) vrb.error("Incorrect cis-window size!");
	vrb.bullet("Cis-window size is " + stb.str(D.options["window"].as < unsigned int > ()) + " bp");
	D.cis_window = D.options["window"].as < unsigned int > ();
	if (D.options.count("chunk")) {
		vector < int > nChunk = D.options["chunk"].as < vector < int > > ();
		if (nChunk.size() != 2 || nChunk[0] > nChunk[1] || nChunk[0] < 0) vrb.error("Incorrect --chunk arguments!");
		vrb.bullet("Chunk = [" + stb.str(nChunk[0]) + "/" + stb.str(nChunk[1]) + "]");
	} else if(D.options.count("region")) vrb.bullet("Region = [" + D.options["region"].as < string > () +"]");
	else if(D.options.count("region-pair")) {
		vector < string > regions = D.options["region-pair"].as < vector < string > > ();
		vrb.bullet("Region for genotypes = [" + regions[0] +"]");
		vrb.bullet("Region for phenotypes = [" + regions[1] +"]");
		D.full_test = true;
	}

	//---------------------------
	// 7. SET AGGREGATION METHODS
	//---------------------------
	int n_aggregation_methods = D.options.count("grp-best") + D.options.count("grp-pca1") + D.options.count("grp-mean");
	if (n_aggregation_methods > 1) vrb.error("Only one of the --grp-XXX options is allowed");
	if (D.options.count("grp-best")) {
		vrb.bullet("Phenotypes are regrouped within groups [method: best]");
		D.grp_mode = GRP_BEST;
	} else if (D.options.count("grp-pca1")) {
		vrb.bullet("Phenotypes are regrouped within groups [method: pca1]");
		D.grp_mode = GRP_PCA1;
	} else if (D.options.count("grp-mean")) {
		vrb.bullet("Phenotypes are regrouped within groups [method: mean]");
		D.grp_mode = GRP_MEAN;
	} else {
		D.grp_mode = GRP_NONE;
	}

	//--------------
	// 8. SET REGION
	//--------------
	if (D.options.count("chunk")) {
		vector < int > nChunk = D.options["chunk"].as < vector < int > > ();
		if (nChunk[0] == 0) {
			D.writeHeader(outFile);
			return;
		} else {
			D.scanPhenotypes(D.options["bed"].as < string > ());
			D.setPhenotypeRegion(nChunk[0] - 1, nChunk[1]);
			D.clear();
			D.regionGenotype = genomic_region(D.regionPhenotype, D.options["window"].as < unsigned int > ());
		}
    } else if (D.options.count("region")){
        if (!D.setPhenotypeRegion(D.options["region"].as < string > ())) vrb.error("Impossible to interpret region [" + D.options["region"].as < string > () + "]");
        D.regionGenotype = genomic_region(D.regionPhenotype, D.options["window"].as < unsigned int > ());
    } else if (D.options.count("region-pair")) {
    	vector < string > regions = D.options["region-pair"].as < vector < string > > ();
    	if (!D.setPhenotypeRegion(regions[0])) vrb.error("Impossible to interpret region [" + regions[0] + "]");
    	if (!D.setGenotypeRegion(regions[1])) vrb.error("Impossible to interpret region [" + regions[1] + "]");
    	D.regionGenotype = genomic_region(D.regionGenotype, D.options["window"].as < unsigned int > ());
    	D.regionPhenotype = genomic_region(D.regionPhenotype, D.options["window"].as < unsigned int > ());
    	D.full_test = true;
    	D.grp_mode = GRP_BEST;
    }
	

	//---------------
	// 9. READ FILES
	//---------------
	D.processBasicOptions();
	D.readSampleFromBED(D.options["bed"].as < string > ());										//Read samples in BED
   	htsFile * fp = hts_open(D.options["vcf"].as < string > ().c_str(),"r");						
   	if (fp->format.format == sam) D.readSampleFromBED(D.options["vcf"].as < string > ());
	else D.readSampleFromVCF(D.options["vcf"].as < string > ());
	hts_close(fp);

	if (D.options.count("cov")) D.readSampleFromCOV(D.options["cov"].as < string > ());			//Read samples in COV
	D.mergeSampleLists();																		//Merge all sample lists
	D.readPhenotypes(D.options["bed"].as < string > ());										//Read data in BED
	D.readGenotypes(D.options["vcf"].as < string > ());											//Read data in VCF
	if (D.options.count("cov")) D.readCovariates(D.options["cov"].as < string > ());
	if (D.options.count("mapping")) D.readThresholds(D.options["mapping"].as < string > ());
	if (D.options.count("nominal")) {
		string argument = D.options["nominal"].as < string >();
		if (!stb.numeric(argument)) D.readThresholds(argument);
	}

	//------------------------
	// 10. INITIALIZE ANALYSIS
	//------------------------
	D.imputeGenotypes();
	D.imputePhenotypes();
	if (D.options.count("cov")) D.residualizePhenotypes();
	D.collapsePhenotypes();
	if (D.options.count("normal")) D.normalTransformPhenotypes();

	//-----------------
	// 11. RUN ANALYSIS
	//-----------------
	switch (D.mode) {
	case CIS_PERM: D.runPermutationPass(outFile); break;
	case CIS_NOMI: D.runNominalPass(outFile); break;
	case CIS_COND: D.runConditionalPass(outFile); break;
	}
}
