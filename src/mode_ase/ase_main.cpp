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

#include "ase_data.h"

void ase_main(vector < string > & argv) {
	ase_data D;

	//-------------------------
	// 1. DECLARE ALL OPTIONS
	//-------------------------
	D.declareBasicOptions();
	boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
	opt_files.add_options()
		("group-by,g", boost::program_options::value< unsigned int >()->default_value(0,"OFF"), "Minimal coverage for a genotype to be considered.")
		("max-depth,d", boost::program_options::value< int >()->default_value(16000,"16000"), "Pileup max-depth.")
		("vcf,v", boost::program_options::value< string >(), "Genotypes in VCF/BCF format.")
		("bam,b", boost::program_options::value< string >(), "Sequence data in BAM/SAM format.")
		("ind,i", boost::program_options::value< string >(), "Sample to be processed.")
		("reg,r", boost::program_options::value< string >()->default_value(""), "Genomic region(s) to be processed.")
		("filtered,f", boost::program_options::value< string >()->default_value(""), "File to output filtered variants.")
		("blacklist,B", boost::program_options::value< string >(), "BED file for blacklisted regions.")
		("fix-chr,F", "Attempt to match chromosomes to the BAM file")
		("out,o", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mFilters\33[0m");
	opt_parameters.add_options()
		("filter-mapping-quality,q", boost::program_options::value< unsigned int >()->default_value(10), "Minimal phred mapping quality for a read to be considered.")
		("filter-base-quality,Q", boost::program_options::value< unsigned int >()->default_value(13), "Minimal phred quality for a base to be considered.")
		("filter-binomial-pvalue,p", boost::program_options::value< double >()->default_value(1.0, "1.0"), "Binomial p-value threshold for ASE in output.")
		("filter-minimal-coverage,c", boost::program_options::value< unsigned int >()->default_value(16), "Minimal coverage for a genotype to be considered.")
		("filter-minimal-coverage-bias,C", boost::program_options::value< unsigned int >()->default_value(10), "Minimal coverage for a genotype to be considered.")
		("filter-minimal-sites-bias,s", boost::program_options::value< unsigned int >()->default_value(100), "Minimum number of sites to calculate a REF bias from.")
		("filter-imputation-score-info-id,I", boost::program_options::value< string >()->default_value("INFO", "INFO"), "The INFO ID of the imputation score in the VCF.")
		("filter-genotype-likelihood-id,L", boost::program_options::value< string >()->default_value("GL", "GL"), "The FORMAT ID of the genotype likelihoods for RR/RA/AA  in the VCF.")
		("filter-imputation-qual,w", boost::program_options::value< double >()->default_value(0.0, "0.0"), "Minimal imputation information score for a variant to be considered.")
		("filter-imputation-prob,z", boost::program_options::value< double >()->default_value(0.0, "0.0"), "Minimal posterior probability for a genotype to be considered.")
		("filter-sample,S", boost::program_options::value< double >()->default_value(1.0, "1.0"), "Randomly subsample sites that are greater than this percentile of all the sites in REF bias calculations.")
		("filter-both-alleles-seen,a", "Require both alleles to be observed in RNA-seq reads for a site for ASE calculations.")
		("keep-homozygote-for-bias,A", "DON'T require both alleles to be observed in RNA-seq reads for a site for REF mapping bias calculations. ")
		("filter-indel-reads,D", "Remove reads that contain indels.")
		("keep-failed-qc,e", "Keep fastq reads that fail sequencing QC (as indicated by the sequencer).")
		("keep-orphan-reads,O", "Keep paired end reads where one of mates is unmapped.")
		("check-proper-pairing,y", "If provided only properly paired reads according to the aligner.")
		("filter-remove-duplicates,u", "Remove duplicate sequencing reads in the process.");

	D.option_descriptions.add(opt_files).add(opt_parameters);

	//-------------------
	// 2. PARSE OPTIONS
	//-------------------
	boost::program_options::variables_map options;
	try {
		boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
		boost::program_options::notify(D.options);
	} catch ( const boost::program_options::error& e ) {
		cerr << "Error parsing [ase] command line :" << string(e.what()) << endl;
		exit(0);
	}

	//---------------------
	// 3. PRINT HELP/HEADER
	//---------------------
	vrb.ctitle("CALLING ALLELE SPECIFIC SITES");
	if (D.options.count("help")) {
		cout << D.option_descriptions << endl;
		exit(EXIT_SUCCESS);
	}

	//-----------------
	// 4. COMMON CHECKS
	//-----------------
	if (!D.options.count("vcf")) vrb.error("Genotype data needs to be specified with --vcf [file.vcf]");
	if (!D.options.count("bam")) vrb.error("Sequence data needs to be specified with --bam [file.bam]");
	if (!D.options.count("ind")) vrb.error("Sample ID needs to be specified with --ind [sample_id]");
	if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");

	//TO DO CHECK PARAMETER VALUES
	D.max_depth = D.options["max-depth"].as < int > ();
	if(D.max_depth == 0) {
		D.max_depth = INT_MAX;
		vrb.warning("Setting pileup max-depth to INT_MAX!");
	}
	if (D.max_depth > 1000000) vrb.warning("Pileup max-depth is above 1M. Potential memory hog!");

	D.param_min_mapQ = D.options["filter-mapping-quality"].as < unsigned int > ();
	D.param_min_baseQ = D.options["filter-base-quality"].as < unsigned int > ();
	D.param_min_cov = D.options["filter-minimal-coverage"].as < unsigned int > ();
	D.param_min_cov_for_ref_alt = D.options["filter-minimal-coverage-bias"].as < unsigned int > ();
	D.param_min_sites_for_ref_alt = D.options["filter-minimal-sites-bias"].as < unsigned int > ();
	D.param_min_pval = D.options["filter-binomial-pvalue"].as < double > ();
	D.param_min_gp = D.options["filter-imputation-prob"].as < double > ();
	D.param_min_iq = D.options["filter-imputation-qual"].as < double > ();
	D.param_sample = D.options["filter-sample"].as < double > ();
	D.param_dup_rd = D.options.count("filter-remove-duplicates");
	D.check_proper_pair = D.options.count("check-proper-pairing");
	D.param_both_alleles_seen = D.options.count("filter-both-alleles-seen");
	D.param_both_alleles_seen_bias = (D.options.count("keep-homozygote-for-bias") == 0);
	D.param_rm_indel = D.options.count("filter-indel-reads");
	D.keep_failqc = D.options.count("keep-failed-qc");
	D.keep_orphan = D.options.count("keep-orphan-reads");
	D.fix_chr = D.options.count("fix-chr");
	D.region_length = D.options["group-by"].as < unsigned int > ();
	D.param_genotype_likelihood_label = D.options["filter-genotype-likelihood-id"].as <string> ();
	D.param_imputation_score_label = D.options["filter-imputation-score-info-id"].as <string> ();
	vrb.bullet("Mapping quality >= " + stb.str(D.param_min_mapQ));
	vrb.bullet("Base quality >= " + stb.str(D.param_min_baseQ));
	vrb.bullet("Coverage >= " + stb.str(D.param_min_cov));
	vrb.bullet("Coverage for bias >= " + stb.str(D.param_min_cov_for_ref_alt));
	vrb.bullet("Binomial p-value threshold = " + stb.str(D.param_min_pval));
	vrb.bullet("Genotype probability >= " + stb.str(D.param_min_gp));
	vrb.bullet("Imputation quality >= " + stb.str(D.param_min_iq));
	vrb.bullet("Both alleles seen for ASE = " + stb.str(D.param_both_alleles_seen));
	vrb.bullet("Both alleles seen for bias = " + stb.str(D.param_both_alleles_seen_bias));
	vrb.bullet("Remove duplicate reads = " + stb.str(D.param_dup_rd));

	//if (D.param_min_cov_for_ref_alt > D.param_min_cov) vrb.error("--filter-minimal-cov-bias cannot be greater than --filter-mininal-cov")

	//------------------------------------------
	// 5. READ FILES / INITIALIZE / RUN ANALYSIS
	//------------------------------------------
	D.processBasicOptions();
	D.readSampleFromVCF(D.options["vcf"].as < string > ());
	D.readSampleFromSTR(D.options["ind"].as < string > ());
	D.mergeSampleLists();
	if (D.sample_count == 0) vrb.error("Could not find [" + D.options["ind"].as < string > () + "] in VCF/BCF file");
	else if (D.sample_count >= 2) vrb.error("More than one sample specified with --ind");
	else vrb.bullet("Target sample is [" + D.sample_id[0] + "]");
	if (D.fix_chr) D.compareChrs(D.options["vcf"].as < string > (), D.options["bam"].as < string > ());
	if (D.options.count("blacklist")) D.readBlacklist(D.options["blacklist"].as < string > ());
	D.readGenotypes2(D.options["vcf"].as < string > (), D.options["reg"].as < string > (), D.options["filtered"].as < string > ());
	if (D.region_length) D.getRegions();
	D.readSequences(D.options["bam"].as < string > ());
	D.calculateRefToAltBias(D.options["filtered"].as < string > ());
	D.calculateASE(D.options["out"].as < string > (), D.options["filtered"].as < string > ());
}
