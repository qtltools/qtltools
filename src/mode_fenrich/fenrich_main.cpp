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

#include "fenrich_data.h"

void fenrich_main(vector < string > & argv) {
	fenrich_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("qtl", boost::program_options::value< string >(), "QTL list in TXT format.")
		("bed", boost::program_options::value< string >(), "Functional annotations in BED format.")
		("tss", boost::program_options::value< string >(), "TSS positions.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("permute", boost::program_options::value< unsigned int >()->default_value(1000), "Permutation number to empirically assess significance.");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [fenrich] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("MEASURE ENRICHMENT OF QTL WITHIN ANNOTATIONS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("qtl")) vrb.error("QTL data needs to be specified with --qtl [file.txt]");
	if (!D.options.count("tss")) vrb.error("TSS data needs to be specified with --tss [file.bed]");
	if (!D.options.count("bed")) vrb.error("Annotation data needs to be specified with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");

	//--------------
	// 6. SET PARAMS
	//--------------
	D.n_permutation = D.options["permute"].as < unsigned int > ();
	vrb.bullet("#permutations = " + stb.str(D.options["permute"].as < unsigned int > ()));

	//---------------------------
	// 7. READ FILES & INITIALIZE
	//---------------------------
	D.processBasicOptions();
	D.readAnnotation(D.options["bed"].as < string > ());
	D.readTSS(D.options["tss"].as < string > ());
	D.readQTL(D.options["qtl"].as < string > ());

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.mapAnnotation2QTL();
	D.runEnrichmentPass(D.options["out"].as < string > ());
}

