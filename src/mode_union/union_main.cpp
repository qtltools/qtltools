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

#include "union_data.h"

void union_main(vector < string > & argv) {
	union_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();  //Mandatory

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< vector < string > >()->multitoken(), "Genotypes in VCF/BCF/BED format.")
		("bed", boost::program_options::value< vector < string > >()->multitoken(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< vector < string > >()->multitoken(), "Covariates in TXT format.")
		("hotspots", boost::program_options::value< string >(), "Hotspots in BED format.")
        ("results", boost::program_options::value< vector < string > >()->multitoken(), "QTLtools results file in TXT format.")
		("out-suffix", boost::program_options::value< string >(), "If provided output files will be suffixed with this.");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("force", "Force overwrite of union files.")
		("normal", "Normal transform the phenotypes.")
		("conditional", "Do conditional analysis.");

	boost::program_options::options_description opt_columns ("\x1B[32mColumns (1-based)\33[0m");
	opt_columns.add_options()
		("pheno-col", boost::program_options::value< unsigned int >()->default_value(1), "Phenotype column")
		("geno-col", boost::program_options::value< unsigned int >()->default_value(7), "Genotype column")
		("rank-col", boost::program_options::value< unsigned int >()->default_value(11), "Conditional analysis rank column")
		("best-col", boost::program_options::value< unsigned int >()->default_value(18), "Conditional analysis best variant column");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("chunk", boost::program_options::value< vector < int > >()->multitoken(), "Specify which chunk needs to be processed")
		("region", boost::program_options::value< string >(), "Region of interest.")
        ("window", boost::program_options::value< unsigned int >()->default_value(1000000), "Size of the cis-window.");

	D.option_descriptions.add(opt_files).add(opt_parameters).add(opt_columns).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [rtc-union] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CALCULATE UNION OF QTLS");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}


	//-----------------
	// 4. COMMON CHECKS
	//-----------------
    if (!D.options.count("vcf")) vrb.error("Genotype data needs to be specified with --vcf [file.vcf]");
    if (!D.options.count("bed")) vrb.error("Phenotype data needs to be specified with --bed [file.bed]");
    if (!D.options.count("results")) vrb.error("Results needs to be specified with --results [file.txt]");
    if (!D.options.count("hotspots")) vrb.error("Output needs to be specified with --hotspots [file.bed]");
    if (D.options["pheno-col"].as < unsigned int > () < 1) vrb.error("--pheno-col must be greater than 0");
    if (D.options["geno-col"].as < unsigned int > () < 1) vrb.error("--geno-col must be greater than 0");
    if (D.options["rank-col"].as < unsigned int > () < 1) vrb.error("--rank-col must be greater than 0");
    if (D.options["best-col"].as < unsigned int > () < 1) vrb.error("--best-col must be greater than 0");
    vector < string > bedFiles = D.options["bed"].as < vector < string > > ();
    vector < string > hitFiles = D.options["results"].as < vector < string > > ();
    vector < string > vcfFiles = D.options["vcf"].as < vector < string > > ();
    vector < string > covFiles = D.options["cov"].as < vector < string > > ();
    if (bedFiles.size() != hitFiles.size()) vrb.error("Unmatched --results and --bed files");
    if (D.options.count("cov") && covFiles.size() != hitFiles.size()) vrb.error("Unmatched --results and --cov files");
    if (vcfFiles.size() == 1) vrb.bullet("Single VCF provided, assuming common VCF file.");
    else if (vcfFiles.size() != hitFiles.size()) vrb.error("Unmatched --results and --vcf files");
    D.no_of_files = hitFiles.size();
    if(!D.options.count("force")){
    	for (int i = 0 ; i < D.no_of_files; i++){
    		string name = hitFiles[i].substr(hitFiles[i].find_last_of("/") + 1) + ".union";
    		if (D.options.count("out-suffix")) name += D.options["out-suffix"].as <string> ();
    		ifstream file(name.c_str());
    		if (file.is_open()){
    			file.close();
    			vrb.error("File [" + name + "] already exists use --force to overwrite");
    		}
    	}
    }

    //--------------
    // 5. SET PARAMS
    //--------------

	if (D.options.count("chunk")) {
		vector < int > nChunk = D.options["chunk"].as < vector < int > > ();
		if (nChunk.size() != 2 || nChunk[0] > nChunk[1]) vrb.error("Incorrect --chunk arguments!");
		vrb.bullet("Chunk = [" + stb.str(nChunk[0]) + "/" + stb.str(nChunk[1]) + "]");
	} else if (D.options.count("region")) vrb.bullet("Region = [" + D.options["region"].as < string > () +"]");
	if(D.options.count("conditional")) vrb.bullet("Doing conditional analysis.");
    D.phenotype_column = D.options["pheno-col"].as < unsigned int > () - 1;
    D.variant_column = D.options["geno-col"].as < unsigned int > () - 1;
    D.rank_column = D.options["rank-col"].as < unsigned int > () - 1;
    D.best_column = D.options["best-col"].as < unsigned int > () - 1;
    vrb.bullet("Phenotype column (0-based) " + stb.str(D.phenotype_column));
    vrb.bullet("Variant column (0-based) " + stb.str(D.variant_column));
    if (D.options.count("conditional")){
        vrb.bullet("Rank column (0-based) " + stb.str(D.rank_column));
        vrb.bullet("Best column (0-based) " + stb.str(D.best_column));
    }

    if (D.options.count("chunk") || D.options.count("region")) vrb.warning("--chunk or --region will not work for trans results");
    if (D.options.count("chunk") && !D.options.count("out-suffix")) vrb.error("--out-suffix is required when --chunk. Otherwise output files will be overwritten.");
    D.readHotspots(D.options["hotspots"].as < string > ());
    //--------------
    // 6SET REGION
    //--------------
    if (D.options.count("chunk")) {
        for (int i =0 ; i < D.no_of_files; i++) D.scanPhenotypes(bedFiles[i]);
        D.setPhenotypeRegion(D.options["chunk"].as < vector < int > > ()[0] - 1, D.options["chunk"].as < vector < int > > ()[1]);
        //outFile += "." + D.regionPhenotype.get();
        D.clearNotHotspot();
        D.deduceGenotypeRegion(D.options["window"].as < unsigned int > ());
    } else if (D.options.count("region")){
        if (!D.setPhenotypeRegion(D.options["region"].as < string > ())) vrb.error("Impossible to interpret region [" + D.options["region"].as < string > () + "]");
        D.deduceGenotypeRegion(D.options["window"].as < unsigned int > ());
    }


	D.processBasicOptions(); //Mandatory
	D.create_unions(hitFiles,vcfFiles,bedFiles,covFiles);
	D.find_unions(hitFiles,vcfFiles,bedFiles,covFiles);
	D.clear();

}
