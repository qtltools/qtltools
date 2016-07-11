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

#include "bamstat_data.h"

void bamstat_main(vector < string > & argv) {
	bamstat_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("bam", boost::program_options::value< string >(), "Sequence data in BAM/SAM format.")
		("bed", boost::program_options::value< string >(), "Annotation data in BED format.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mFilters\33[0m");
	opt_parameters.add_options()
		("filter-mapping-quality", boost::program_options::value< unsigned int >()->default_value(10), "Minimal phred mapping quality for a read to be considered.")
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
		cerr << "Error parsing [bamstat] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CALCULATE BASIC QC METRICS FOR A BAM FILE");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("bam")) vrb.error("Sequence data needs to be specified with --bam [file.bam]");
	if (!D.options.count("bed")) vrb.error("Annotation data needs to be specified with --bed [file.bed]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");

	//TO DO CHECK PARAMETER VALUES
	D.param_min_mapQ = D.options["filter-mapping-quality"].as < unsigned int > ();
	D.param_dup_rd = (D.options.count("filter-keep-duplicates") != 0);
	vrb.bullet("Mapping quality >= " + stb.str(D.param_min_mapQ));
	vrb.bullet("Keep duplicate reads = " + stb.str(D.param_dup_rd));

	//------------------------------------------
	// 5. READ FILES / INITIALIZE / RUN ANALYSIS
	//------------------------------------------

	D.processBasicOptions();
	D.readAnnotationsBED(D.options["bed"].as < string > ());
	D.readSequences(D.options["bam"].as < string > ());
	D.writeOutput(D.options["out"].as < string > ());
}
