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

#include "correct_data.h"

void correct_main(vector < string > & argv) {
	correct_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF/BED format.")
		("bed", boost::program_options::value< string >(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< string >(), "Covariates in TXT format.")
		("qtl", boost::program_options::value< vector < string > >()->multitoken(), "QTL covariates.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_param ("\x1B[32mParameters\33[0m");
	opt_param.add_options()
		("normal", "Normal transformation fo the data.");

	D.option_descriptions.add(opt_files).add(opt_param);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [correct] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CORRECTING GENOTYPES OR PHENOTYPES FOR COVARIATES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if ((D.options.count("vcf") + D.options.count("bed")) != 1) vrb.error("One input file has to be specified using either --vcf [file.vcf] or --bed [file.bed]");
	if ((D.options.count("cov") + D.options.count("normal")) < 1) vrb.error("At least one data transformation has to be specified with --cov [file.cov] and/or --normal");

	//---------
	// 5. MODES
	//---------

	//MODE1: correcting phenotypes
	if (D.options.count("bed")) {
		D.mode = CORRECT_BED;
		if (D.options.count("cov")) { vrb.bullet("Correct phenotypes using linear model: Phe ~ Cov"); D.residualize = true; }
		if (D.options.count("normal")) { vrb.bullet("Quantile normal transform phenotypes"); D.normalize = true; }
	}
	//MODE2: correcting genotypes
	if (D.options.count("vcf")) {
		D.mode = CORRECT_VCF;
		if (D.options.count("cov")) { vrb.bullet("Correct genotypes using linear model: Phe ~ Cov"); D.residualize = true; }
		if (D.options.count("normal")) { vrb.bullet("Quantile normal transform genotypes"); D.normalize = true; }
	}

	//---------------------------
	// 6. READ FILES & INITIALIZE
	//---------------------------
	D.processBasicOptions();
	if (D.options.count("bed")) D.readSampleFromBED(D.options["bed"].as < string > ());
	if (D.options.count("vcf")) D.readSampleFromVCF(D.options["vcf"].as < string > ());
	if (D.options.count("cov")) D.readSampleFromCOV(D.options["cov"].as < string > ());
	D.mergeSampleLists();
	if (D.options.count("cov")) D.readCovariates(D.options["cov"].as < string > ());
	if (D.options.count("qtl")) {
		vector < string > files = D.options["qtl"].as < vector < string > > ();
		assert(files .size() == 2);
		D.readQTLCovariates(files[0], files[1]);
	}

	//----------------
	// 7. RUN ANALYSIS
	//----------------
	if (D.options.count("bed")) D.processBED(D.options["bed"].as < string > (), D.options["out"].as < string > ());
	if (D.options.count("vcf")) D.processVCF(D.options["vcf"].as < string > (), D.options["out"].as < string > ());
}
