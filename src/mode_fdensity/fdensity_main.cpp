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

#include "fdensity_data.h"

void fdensity_main(vector < string > & argv) {
	fdensity_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bed", boost::program_options::value< string >(), "Functional annotations in BED format.")
		("qtl", boost::program_options::value< string >(), "QTL positions.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("window", boost::program_options::value< int >()->default_value(1000000), "Window size around TSS in bp.")
		("bin", boost::program_options::value< int >()->default_value(1000), "Bin size in bp.");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [fdensity] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("ANNOTATION DENSITY ARROUND QTLs");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("qtl")) vrb.error("QTL data needs to be specified with --qtl [file.bed]");
	if (!D.options.count("bed")) vrb.error("Annotation data needs to be specified with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");

	//--------------
	// 6. SET PARAMS
	//--------------
	D.window = D.options["window"].as < int > ();
	D.bin = D.options["bin"].as < int > ();
	vrb.bullet("window = " + stb.str(D.options["window"].as < int > ()));
	vrb.bullet("bin = " + stb.str(D.options["bin"].as < int > ()));

	//---------------------------
	// 7. READ FILES & INITIALIZE
	//---------------------------
	D.processBasicOptions();
	D.readAnnotation(D.options["bed"].as < string > ());
	D.readQTL(D.options["qtl"].as < string > ());
	D.buildIntervalTrees();

	//----------------
	// 8. RUN ANALYSIS
	//----------------
	D.runDensityCalculation(D.options["out"].as < string > ());
}

