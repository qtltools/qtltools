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

#include "rep_data.h"

void rep_main(vector < string > & argv) {
	rep_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF/BED format.")
		("bed", boost::program_options::value< string >(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< string >(), "Covariates in TXT format.")
		("qtl", boost::program_options::value< string >(), "QTLs to replicate in TXT format.")
		("out", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("normal", "Normal transform the phenotypes.");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [rep] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("REPLICATE QTLS");
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
	string outFile = D.options["out"].as < string > ();

	//--------------
	// 5. SET PARAMS
	//--------------
	if (D.options.count("cov")) vrb.bullet("Linear model: Phe ~ Var + Cov");
	else vrb.bullet("Linear model: Phe ~ Var");

	//---------------
	// 6. READ FILES
	//---------------
	D.processBasicOptions();
	D.readSampleFromBED(D.options["bed"].as < string > ());										//Read samples in BED
   	htsFile * fp = hts_open(D.options["vcf"].as < string > ().c_str(),"r");						
   	if (fp->format.format == sam) D.readSampleFromBED(D.options["vcf"].as < string > ());
	else D.readSampleFromVCF(D.options["vcf"].as < string > ());
	hts_close(fp);
	if (D.options.count("cov")) D.readSampleFromCOV(D.options["cov"].as < string > ());			//Read samples in COV
	D.mergeSampleLists();																		//Merge all sample lists

	D.readQTLs(D.options["qtl"].as < string > ());
	D.filter_phenotype.addInclusion(D.qtl_ids.first);
	D.filter_genotype.addInclusion(D.qtl_ids.second);

	D.readPhenotypes(D.options["bed"].as < string > ());										//Read data in BED
	D.readGenotypes(D.options["vcf"].as < string > ());											//Read data in VCF
	if (D.options.count("cov")) D.readCovariates(D.options["cov"].as < string > ());
	D.mapping();

	//------------------------
	// 10. INITIALIZE ANALYSIS
	//------------------------
	D.imputeGenotypes();
	D.imputePhenotypes();
	if (D.options.count("cov")) D.residualizePhenotypes();
	if (D.options.count("normal")) D.normalTransformPhenotypes();

	//-----------------
	// 11. RUN ANALYSIS
	//-----------------
	D.runNominalPass(outFile);
}
