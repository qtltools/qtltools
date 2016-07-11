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

#include "rtc_data.h"

void rtc_main(vector < string > & argv) {
	rtc_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();  //Mandatory

	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF/BED format.")
		("bed", boost::program_options::value< string >(), "Phenotypes in BED format.")
		("cov", boost::program_options::value< string >(), "Covariates in TXT format.")
		("hotspots", boost::program_options::value< string >(), "Hotspots in BED format.")
		("out", boost::program_options::value< string >(), "Output file.")
        ("stats-vcf", boost::program_options::value< string >(), "Genotypes from which D' and r2 are calculated in  VCF/BCF format. (Defaults to --vcf, MUST HAVE PHASED GENOTYPES")
        ("stats-vcf-include-samples", boost::program_options::value< string >(), "Include sample list for --stats-vcf")
        ("stats-vcf-exclude-samples", boost::program_options::value< string >(), "Exclude sample list for --stats-vcf");

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("normal", "Normal transform the phenotypes.")
		("conditional", "Do conditional analysis.")
		("debug", "Print debugging info for sampling to STDERR.")
        ("individual-Dprime", "Will calculate Dprime on an individual variant basis. If no not provided Dprime will not be calculated after first unphased genotype is encountered.")
        ("no-mem", "Don't keep results of calculations that may be used multiple times in memory (Use to save some RAM, but lose speed).")
		("window", boost::program_options::value< unsigned int >()->default_value(1000000), "Size of the cis-window.")
		("sample", boost::program_options::value< unsigned int >()->default_value(0), "Sample iterations to assess RTC significance.")
        ("R2-threshold", boost::program_options::value< double >()->default_value(0.5), "R2 threshold used in sampling")
        ("D-prime-threshold", boost::program_options::value< double >()->default_value(2,"OFF"), "If the pairs of variants have a D' greater than this the RTC calculation is extended to multiple coldspots. (Assumes D' can be calculated");

	boost::program_options::options_description opt_columns ("\x1B[32mColumns (1-based)\33[0m");
	opt_columns.add_options()
		("pheno-col", boost::program_options::value< unsigned int >()->default_value(1), "Phenotype column")
		("geno-col", boost::program_options::value< unsigned int >()->default_value(7), "Genotype column")
		("rank-col", boost::program_options::value< unsigned int >()->default_value(11), "Conditional analysis rank column")
		("best-col", boost::program_options::value< unsigned int >()->default_value(18), "Conditional analysis best variant column");

	boost::program_options::options_description opt_modes ("\x1B[32mAnalysis type\33[0m");
	opt_modes.add_options()
		("gwas-cis", boost::program_options::value< vector < string > >()->multitoken(), "MODE1: RTC for GWAS and cis-eQTL integration.")
        ("gwas-trans", boost::program_options::value< vector < string > >()->multitoken(), "MODE2: RTC for GWAS and trans-eQTL integration.")
		("mergeQTL-cis", boost::program_options::value< vector < string > >()->multitoken(), "MODE3: RTC for cis-eQTL and cis-eQTL integration.")
		("mergeQTL-trans", boost::program_options::value< vector < string > >()->multitoken(), "MODE4: RTC for trans-eQTL and trans-eQTL integration.");

	boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
	opt_parallel.add_options()
		("chunk", boost::program_options::value< vector < int > >()->multitoken(), "Specify which chunk needs to be processed")
		("region", boost::program_options::value< string >(), "Region of interest.");

	D.option_descriptions.add(opt_files).add(opt_parameters).add(opt_columns).add(opt_modes).add(opt_parallel);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [rtc] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CALCULATE RTC SCORE");
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
    if (!D.options.count("hotspots")) vrb.error("Output needs to be specified with --hotspots [file.bed]");
    int nMode = D.options.count("gwas-cis") +  D.options.count("gwas-trans") + D.options.count("mergeQTL-cis") + D.options.count("mergeQTL-trans");
    if (nMode != 1) vrb.error("Please specify only one of these options [--gwas-cis, --gwas-trans, --mergeQTL-cis, --mergeQTL-trans]");
    string outFile = D.options["out"].as < string > ();
    if (D.options["pheno-col"].as < unsigned int > () < 1) vrb.error("--pheno-col must be greater than 0");
    if (D.options["geno-col"].as < unsigned int > () < 1) vrb.error("--geno-col must be greater than 0");
    if (D.options["rank-col"].as < unsigned int > () < 1) vrb.error("--rank-col must be greater than 0");
    if (D.options["best-col"].as < unsigned int > () < 1) vrb.error("--best-col must be greater than 0");
    if (D.options["sample"].as <unsigned int> () == 0 && D.options.count("debug") ) vrb.error("--debug only applies when --sample");

	//---------
	// 5. MODES
	//---------
    vector < string > RTCfiles;
	if (D.options.count("gwas-cis")) {
        D.mode = RTC_MODE1;
		vrb.bullet("TASK: RTC on GWAS variants and cis-eQTLs");
        RTCfiles = D.options["gwas-cis"].as < vector < string > > ();
        if (RTCfiles.size() != 2) vrb.error("Please provide 2 input files for --gwas-cis");
	}

    if (D.options.count("gwas-trans")) {
        D.mode = RTC_MODE2;
        vrb.bullet("TASK: RTC on GWAS variants and trans-eQTLs");
        RTCfiles = D.options["gwas-trans"].as < vector < string > > ();
        if (RTCfiles.size() != 2) vrb.error("Please provide 2 input files for --gwas-trans");
    }

    if (D.options.count("mergeQTL-cis")) {
        D.mode = RTC_MODE3;
        vrb.bullet("TASK: RTC on pairs of cis-QTLs");
        RTCfiles = D.options["mergeQTL-cis"].as < vector < string > > ();
        if (RTCfiles.size() != 2) vrb.error("Please provide 2 input files for --mergeQTL-cis");
    }

    if (D.options.count("mergeQTL-trans")) {
        D.mode = RTC_MODE4;
        vrb.bullet("TASK: RTC on pairs of trans-QTLs");
        RTCfiles = D.options["mergeQTL-trans"].as < vector < string > > ();
        if (RTCfiles.size() != 2) vrb.error("Please provide 2 input files for --mergeQTL-trans");
    }


	//--------------
	// 6. SET PARAMS
	//--------------
    if (D.options["window"].as < unsigned int > () > 1000000000) vrb.error("Incorrect cis-window size!");
    vrb.bullet("Cis-window size is " + stb.str(D.options["window"].as < unsigned int > ()) + " bp");
    D.cis_window = D.options["window"].as < unsigned int > ();
    if (D.options.count("chunk")) {
        vector < int > nChunk = D.options["chunk"].as < vector < int > > ();
        if (nChunk.size() != 2 || nChunk[0] > nChunk[1]) vrb.error("Incorrect --chunk arguments!");
        vrb.bullet("Chunk = [" + stb.str(nChunk[0]) + "/" + stb.str(nChunk[1]) + "]");
        if (nChunk[0] == 0){
        	output_file fdo(outFile);
        	if (fdo.fail()) vrb.error("Cannot open file [" + outFile + "]");
        	fdo <<"other_variant our_variant phenotype other_variant_chr other_variant_start other_variant_rank our_variant_chr our_variant_start our_variant_rank phenotype_chr phenotype_start distance_between_variants distance_between_other_variant_and_pheno other_variant_region_index our_variant_region_index region_start region_end RTC D' r^2";
        	if (D.options.count("sample")) fdo << " better_proportion_under_H0 worse_proportion_under_H1";
        	fdo << endl;
        	fdo.close();
        	vrb.leave("Header written");
        }
    } else if (D.options.count("region")) vrb.bullet("Region = [" + D.options["region"].as < string > () +"]");
    if(D.options.count("conditional")) vrb.bullet("Doing conditional analysis.");
    D.Dprime_cutoff = D.options["D-prime-threshold"].as < double > ();
    D.R2_cutoff = D.options["R2-threshold"].as < double > ();
    D.sample_iterations = D.options["sample"].as < unsigned int > ();
    if(D.Dprime_cutoff < 0 ) vrb.error("Wrong D-prime threshold");
    if(D.R2_cutoff < 0 || D.R2_cutoff > 1) vrb.error("Wrong R2 threshold");
    if(D.sample_iterations) vrb.bullet(stb.str(D.sample_iterations) + " sample iterations.");
    if(D.sample_iterations) vrb.bullet(stb.str(D.R2_cutoff) + " R2 threshold.");
    if(D.Dprime_cutoff <= 1 )vrb.bullet(stb.str(D.Dprime_cutoff) + " D' threshold.");
    else D.calculate_Dprime_R2 = false;
    D.phenotype_column = D.options["pheno-col"].as < unsigned int > () - 1;
    D.variant_column = D.options["geno-col"].as < unsigned int > () - 1;
    D.rank_column = D.options["rank-col"].as < unsigned int > () - 1;
    D.best_column = D.options["best-col"].as < unsigned int > () - 1;
    if(D.options.count("no-mem")) D.DprimeR2inMem = false;
    vrb.bullet("Phenotype column (0-based) " + stb.str(D.phenotype_column));
    vrb.bullet("Variant column (0-based) " + stb.str(D.variant_column));
    if (D.options.count("conditional")){
        vrb.bullet("Rank column (0-based) " + stb.str(D.rank_column));
        vrb.bullet("Best column (0-based) " + stb.str(D.best_column));
    }

    //--------------
    // 7. SET REGION
    //--------------
    if (D.options.count("chunk")) {
        D.scanPhenotypes(D.options["bed"].as < string > ());
        D.setPhenotypeRegion(D.options["chunk"].as < vector < int > > ()[0] - 1, D.options["chunk"].as < vector < int > > ()[1]);
        //outFile += "." + D.regionPhenotype.get();
        D.clear();
    } else if (D.options.count("region")){
        if (!D.setPhenotypeRegion(D.options["region"].as < string > ())) vrb.error("Impossible to interpret region [" + D.options["region"].as < string > () + "]");
    }


	//---------------------------
	// 8. READ FILES & INITIALIZE
	//---------------------------
    if (D.mode == RTC_MODE1 || D.mode == RTC_MODE3){
        D.processBasicOptions(); //Mandatory
        D.readHotspots(D.options["hotspots"].as < string > ());
        D.deduceGenotypeRegion(D.options["window"].as < unsigned int > ());
        D.readSampleFromBED(D.options["bed"].as < string > ());									//Read samples in BED
        D.readSampleFromVCF(D.options["vcf"].as < string > ());									//Read samples in VCF
        if (D.options.count("cov")) D.readSampleFromCOV(D.options["cov"].as < string > ());			//Read samples in COV
        D.mergeSampleLists();																	//Merge all sample lists
        D.readPhenotypes(D.options["bed"].as < string > ());										//Read data in BED
        D.readGenotypes(D.options["vcf"].as < string > ());										//Read data in VCF
        if (D.options.count("cov")) D.readCovariates(D.options["cov"].as < string > ());
        if (D.calculate_Dprime_R2){
            if (D.options.count("stats-vcf") && D.options["stats-vcf"].as < string > () !=  D.options["vcf"].as < string > ()){
                D.setStatsVCF(D.options["stats-vcf"].as < string > ());
                if (D.options.count("stats-vcf-exclude-samples")) D.readSampleExclusionStats(D.options["stats-vcf-exclude-samples"].as <string> ());
                if (D.options.count("stats-vcf-include-samples")) D.readSampleInclusionStats(D.options["stats-vcf-include-samples"].as <string> ());
            }else{
                D.setStatsVCF(D.options["vcf"].as < string > ());
                D.copyIncludeExclude();
            }
            D.checkStatsVCF();
        }
        D.imputeGenotypes();
        D.imputePhenotypes();
        if (D.options.count("cov")) D.residualizePhenotypes();
        if (D.options.count("normal")) D.normalTransformPhenotypes();
        D.mapVariantsToColdspots();
    }else{
        //SCAN TO FIND A LIST OF SNPS TO READ INTO MEMORY
        D.processBasicOptions(); //Mandatory
        D.readHotspots(D.options["hotspots"].as < string > ());
        D.scanPhenotypes(D.options["bed"].as < string > ());
        D.scanGenotypes(D.options["vcf"].as < string > ());
        if (D.calculate_Dprime_R2){
            if (D.options.count("stats-vcf") && D.options["stats-vcf"].as < string > () !=  D.options["vcf"].as < string > ()){
                D.setStatsVCF(D.options["stats-vcf"].as < string > ());
                if (D.options.count("stats-vcf-exclude-samples")) D.readSampleExclusionStats(D.options["stats-vcf-exclude-samples"].as <string> ());
                if (D.options.count("stats-vcf-include-samples")) D.readSampleInclusionStats(D.options["stats-vcf-include-samples"].as <string> ());
            }else{
                D.setStatsVCF(D.options["vcf"].as < string > ());
                D.copyIncludeExclude();
            }
            D.checkStatsVCF();
        }
        D.mapVariantsToColdspots();
        switch (D.mode) {
            case RTC_MODE2:
                if (D.options.count("conditional")) D.gwas_trans_conditional(RTCfiles[0], RTCfiles[1]);
                else  D.gwas_trans(RTCfiles[0], RTCfiles[1]);
                break;
            case RTC_MODE4:
                if (D.options.count("conditional")) {
                	D.mergeqtl_trans_conditional(RTCfiles[0], RTCfiles[1]);
                	D.createTransLists();
                }else {
                	D.mergeqtl_trans(RTCfiles[0], RTCfiles[1]);
                	D.createTransLists();
                }
                break;
        }
        D.clear();
        //DO THE FINAL READING
        D.processBasicOptions(); //Mandatory
        D.readHotspots(D.options["hotspots"].as < string > ());
        D.readSampleFromBED(D.options["bed"].as < string > ());									//Read samples in BED
        D.readSampleFromVCF(D.options["vcf"].as < string > ());									//Read samples in VCF
        if (D.options.count("cov")) D.readSampleFromCOV(D.options["cov"].as < string > ());			//Read samples in COV
        D.mergeSampleLists();																	//Merge all sample lists
        D.readPhenotypes(D.options["bed"].as < string > ());										//Read data in BED
        D.readGenotypes(D.options["vcf"].as < string > ());										//Read data in VCF
        if (D.options.count("cov")) D.readCovariates(D.options["cov"].as < string > ());
        D.imputeGenotypes();
        D.imputePhenotypes();
        if (D.options.count("cov")) D.residualizePhenotypes();
        if (D.options.count("normal")) D.normalTransformPhenotypes();
        D.mapVariantsToColdspots();
    }

	//----------------
	// 9. RUN ANALYSIS
	//----------------

    switch (D.mode) {
        case RTC_MODE1:
            if (D.options.count("conditional")) D.gwas_cis_conditional(RTCfiles[0], RTCfiles[1]);
            else  D.gwas_cis(RTCfiles[0], RTCfiles[1]);
            break;
        case RTC_MODE2:
            if (D.options.count("conditional")) D.gwas_trans_conditional(RTCfiles[0], RTCfiles[1]);
            else  D.gwas_trans(RTCfiles[0], RTCfiles[1]);
            break;
        case RTC_MODE3:
            if (D.options.count("conditional")) D.mergeqtl_cis_conditional(RTCfiles[0], RTCfiles[1]);
            else D.mergeqtl_cis(RTCfiles[0], RTCfiles[1]);
            break;
        case RTC_MODE4:
            if (D.options.count("conditional")) D.mergeqtl_trans_conditional(RTCfiles[0], RTCfiles[1]);
            else D.mergeqtl_trans(RTCfiles[0], RTCfiles[1]);
            break;
    }
    //Print warnings
    set <string>::iterator it;
    if (D.unfound_regions.size()){
        vrb.warning(stb.str(D.unfound_regions.size()) + " genotypes were missing from [" + D.stats_vcf_file + "]");
        for (it = D.unfound_regions.begin(); it != D.unfound_regions.end(); it++) vrb.print(*it);
    }
    if (D.unphased_regions.size()){
        vrb.warning(stb.str(D.unphased_regions.size()) + " genotypes were unphased in [" + D.stats_vcf_file + "]");
        for (it = D.unphased_regions.begin(); it != D.unphased_regions.end(); it++) vrb.print(*it);
    }
    if (D.no_variance_regions.size()){
        vrb.warning(stb.str(D.no_variance_regions.size()) + " genotypes had no variance in [" + D.stats_vcf_file + "]");
        for (it = D.no_variance_regions.begin(); it != D.no_variance_regions.end(); it++) vrb.print(*it);
    }
    if (D.unmatched_alleles.size()){
        vrb.warning(stb.str(D.unmatched_alleles.size()) + " genotypes had mismatching alleles between [" + D.options["vcf"].as < string > () + "] and [" + D.stats_vcf_file + "]");
        for (it = D.unmatched_alleles.begin(); it != D.unmatched_alleles.end(); it++) vrb.print(*it);
    }
    if (D.unfound_ids.size()){
        vrb.warning(stb.str(D.unfound_ids.size()) + " genotypes could not be found in [" + D.options["vcf"].as < string > () + "]");
        for (it = D.unfound_ids.begin(); it != D.unfound_ids.end(); it++) vrb.print(*it);
    }
    if (D.unfound_phenotypes.size()){
        vrb.warning(stb.str(D.unfound_phenotypes.size()) + " phenotypes could not be found in [" + D.options["bed"].as < string > () + "]");
        for (it = D.unfound_phenotypes.begin(); it != D.unfound_phenotypes.end(); it++) vrb.print(*it);
    }
    D.calculateRTC(outFile);
}
