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
		("vcf", boost::program_options::value< vector <string> >()->multitoken(), "Genotypes in VCF/BCF/BED format.")
		("bed", boost::program_options::value< vector <string> >()->multitoken(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< vector <string> >()->multitoken(), "Covariates in TXT format.")
		("qtl", boost::program_options::value< vector <string> >()->multitoken(), "QTLs to replicate in TXT format.")
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


	vector <string> vcfs = D.options["vcf"].as < vector <string> > ();
	vector <string> beds = D.options["bed"].as < vector <string> > ();
	vector <string> covs = D.options["cov"].as < vector <string> > ();
	vector <string> qtls = D.options["qtl"].as < vector <string> > ();

	if (beds.size() == 1 && vcfs.size() != 1) vrb.error("Single BED file requires a single VCF file!");
	if (vcfs.size() != 1 && beds.size() != vcfs.size()) vrb.error("Number of BED files does not match the number of VCF files!");
	if (covs.size() != 0 && covs.size() != beds.size()) vrb.error("Number of BED files does not match the number of covariate files!");
	bool skipSameIndex = qtls.size() == beds.size() && qtls.size() != 1 && beds.size() != 1;

	bool first = true;
	for (int b = 0; b < beds.size(); b++){
		for (int q = 0; q < qtls.size(); q++){
			//cerr << beds[b] << " " << qtls[q] << endl;
			if (skipSameIndex && b == q) continue;
			//---------------
			// 6. READ FILES
			//---------------
			D.processBasicOptions();
			D.readSampleFromBED(beds[b]);
			string inVCF = vcfs.size() != 1 ? vcfs[b] : vcfs[0];
			htsFile * fp = hts_open(inVCF.c_str(),"r");
			if (fp->format.format == sam) D.readSampleFromBED(inVCF);
			else D.readSampleFromVCF(inVCF);
			hts_close(fp);
			if (D.options.count("cov")) D.readSampleFromCOV(covs[b]);			//Read samples in COV
			D.mergeSampleLists();																		//Merge all sample lists

			D.readQTLs(qtls[q]);
			D.filter_phenotype.addInclusion(D.qtl_ids.first);
			D.filter_genotype.addInclusion(D.qtl_ids.second);

			D.readPhenotypes(beds[b]);										//Read data in BED
			D.readGenotypes(inVCF);											//Read data in VCF
			if (D.options.count("cov")) D.readCovariates(covs[b]);
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

			D.runNominalPass(outFile, qtls[q], beds[b],first);
			D.clear();
			first = false;
		}
	}
}
