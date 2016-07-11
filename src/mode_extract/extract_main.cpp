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

#include "extract_data.h"

void extract_main(vector < string > & argv) {
	extract_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF format.")
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< string >(), "Covariates in TXT format.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("region", boost::program_options::value< string >(), "Region of interest.");

	D.option_descriptions.add(opt_files).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [extract] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("DATA EXTRACTION");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if ((D.options.count("vcf") + D.options.count("bed") + D.options.count("cov")) == 0) vrb.error("At least one input file has to be specified using either --vcf [file.vcf], --bed [file.bed] or --cov [file.txt]");
	if (!D.options.count("region")) vrb.warning("Please use --region to speed up data extraction for phenotype and genotype data!");

	//--------------
	// 5. SET REGION
	//--------------
	if (D.options.count("region") && !D.regionData.parse(D.options["region"].as < string > ()))
		vrb.error("Impossible to interpret region [" + D.options["region"].as < string > () + "]");

	//--------------
	// 6. READ FILES
	//--------------
	D.processBasicOptions();
	vector < string > bed_list = D.options["bed"].as < vector < string > > ();
	for (int b = 0 ; b < bed_list.size() ; b ++) D.readSampleFromBED(bed_list[b]);
	if (D.options.count("vcf")) D.readSampleFromVCF(D.options["vcf"].as < string > ());
	if (D.options.count("cov")) D.readSampleFromCOV(D.options["cov"].as < string > ());
	D.mergeSampleLists();

	for (int b = 0 ; b < bed_list.size() ; b ++) D.readBED(bed_list[b]);
	if (D.options.count("vcf")) D.readVCF(D.options["vcf"].as < string > ());
	if (D.options.count("cov")) D.readCOV(D.options["cov"].as < string > ());

	D.imputeMissing();

	D.writeOUT(D.options["out"].as < string > ());
}
