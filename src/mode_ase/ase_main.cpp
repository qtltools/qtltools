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
		("vcf,v", boost::program_options::value< string >(), "Genotypes in VCF/BCF format.")
		("bam,b", boost::program_options::value< string >(), "Sequence data in BAM/SAM format.")
		("fasta,G", boost::program_options::value< string >(), "Genome sequence in FASTA format.")
		("ind,i", boost::program_options::value< string >(), "Sample to be processed.")
		("reg,r", boost::program_options::value< string >()->default_value(""), "Genomic region(s) to be processed.")
		("blacklist,B", boost::program_options::value< string >(), "BED file for blacklisted regions.")
		("fix-chr,F", "Attempt to match chromosome names to the BAM file")
		("auto-flip,x", "Attempt to fix reference allele mistaches. Requires a fasta file for the reference sequence.")
		("group-by,g", boost::program_options::value< unsigned int >()->default_value(0,"OFF"), "Group variants separated by this much into batches. This allows you not to stream the whole BAM file and may improve running time.")
		("max-depth,d", boost::program_options::value< int >()->default_value(16000,"16000"), "Pileup max-depth.")
		("filtered,f", boost::program_options::value< string >()->default_value(""), "File to output filtered variants.")
		("out,o", boost::program_options::value< string >(), "Output file.");

	boost::program_options::options_description opt_parameters ("\x1B[32mFilters\33[0m");
	opt_parameters.add_options()
		("filter-mapping-quality,q", boost::program_options::value< unsigned int >()->default_value(10), "Minimum phred mapping quality for a read to be considered.")
		("filter-base-quality,Q", boost::program_options::value< unsigned int >()->default_value(13), "Minimum phred quality for a base to be considered.")
		("filter-binomial-pvalue,p", boost::program_options::value< double >()->default_value(1.0, "1.0"), "Binomial p-value threshold for ASE in output.")
		("filter-coverage,c", boost::program_options::value< unsigned int >()->default_value(16), "Minimum coverage for a genotype to be considered.")
		("filter-coverage-bias,C", boost::program_options::value< unsigned int >()->default_value(10), "Minimum coverage for a genotype to be considered.")
		("filter-min-sites-bias,s", boost::program_options::value< unsigned int >()->default_value(100), "Minimum number of sites to calculate a REF bias from for a specific REF/ALT pair. The REF bias for pairs with less than this many sites will be calculated from all sites.")
		("filter-imputation-score-info-id,I", boost::program_options::value< string >()->default_value("INFO", "INFO"), "The INFO ID of the imputation score in the VCF.")
		("filter-genotype-likelihood-id,L", boost::program_options::value< string >()->default_value("GL", "GL"), "The FORMAT ID of the genotype likelihoods for RR/RA/AA  in the VCF.")
		("filter-imputation-qual,w", boost::program_options::value< double >()->default_value(0.0, "0.0"), "Minimum imputation information score for a variant to be considered.")
		("filter-imputation-prob,z", boost::program_options::value< double >()->default_value(0.0, "0.0"), "Minimum posterior probability for a genotype to be considered.")
		("filter-sample,S", boost::program_options::value< double >()->default_value(1.0, "1.0"), "Randomly subsample sites that are greater than this percentile of all the sites in REF bias calculations.")
		("filter-both-alleles-seen,a", "Require both alleles to be observed in RNA-seq reads for a site for ASE calculations.")
		("keep-homozygotes-for-bias,A", "DON'T require both alleles to be observed in RNA-seq reads for a site for REF mapping bias calculations. (NOT RECOMMENDED)")
		("filter-indel-reads,D", "Remove reads that contain indels.")
		("keep-failed-qc,e", "Keep fastq reads that fail sequencing QC (as indicated by the sequencer).")
		("keep-orphan-reads,O", "Keep paired end reads where one of mates is unmapped.")
		("check-proper-pairing,y", "If provided only properly paired reads according to the aligner.")
		("check-orientation,X", "If provided only mate pairs where both mates are on the same chromosome, first mate is on the +ve strand and the second is on the -ve strand, and the second mate does not start before the first mate will be considered.")
		("filter-remove-duplicates,u", "Remove duplicate sequencing reads in the process.")
		("legacy-options,j", "Replicate legacy options used. (DO NOT USE).");

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

	if (D.options["filter-mapping-quality"].defaulted()) vrb.warning("You used the default mapping quality threshold. Usually this should be adjusted for the aligner used, is this intentional?");

	if (D.options.count("fasta") == 0 && D.options.count("auto-flip")) vrb.warning("Ignoring --auto-flip since no --fasta is provided!");

	D.param_min_mapQ = D.options["filter-mapping-quality"].as < unsigned int > ();
	D.param_min_baseQ = D.options["filter-base-quality"].as < unsigned int > ();
	D.param_min_pval = D.options["filter-binomial-pvalue"].as < double > ();
	D.param_min_cov = D.options["filter-coverage"].as < unsigned int > ();
	D.param_min_cov_for_ref_alt = D.options["filter-coverage-bias"].as < unsigned int > ();
	D.param_min_sites_for_ref_alt = D.options["filter-min-sites-bias"].as < unsigned int > ();
	D.param_imputation_score_label = D.options["filter-imputation-score-info-id"].as <string> ();
	D.param_genotype_likelihood_label = D.options["filter-genotype-likelihood-id"].as <string> ();
	D.param_min_gp = D.options["filter-imputation-prob"].as < double > ();
	D.param_min_iq = D.options["filter-imputation-qual"].as < double > ();
	D.param_sample = D.options["filter-sample"].as < double > ();
	D.param_both_alleles_seen = D.options.count("filter-both-alleles-seen");
	D.param_both_alleles_seen_bias = (D.options.count("keep-homozygotes-for-bias") == 0);
	D.param_rm_indel = D.options.count("filter-indel-reads");
	D.keep_failqc = D.options.count("keep-failed-qc");
	D.keep_orphan = D.options.count("keep-orphan-reads");
	D.check_proper_pair = D.options.count("check-proper-pairing");
	D.check_orientation = D.options.count("check-orientation");
	D.param_dup_rd = D.options.count("filter-remove-duplicates");
	D.fix_chr = D.options.count("fix-chr");
	D.auto_flip = D.options.count("auto-flip");
	D.region_length = D.options["group-by"].as < unsigned int > ();

	if(D.options.count("legacy-options")){
		if (D.param_min_baseQ != 13) vrb.warning("--filter-base-quality is overwritten!");
		if (D.max_depth != 8000) vrb.warning("--max-depth is overwritten!");
		if (!D.param_dup_rd) vrb.warning("--filter-remove-duplicates is overwritten!");
		if (!D.check_proper_pair) vrb.warning("--check-proper-pairing is overwritten!");
		if (D.keep_orphan) vrb.warning("--keep-orphan-reads is overwritten!");
		if (D.keep_failqc) vrb.warning("--keep-failqc is overwritten!");
		if (D.param_rm_indel) vrb.warning("--filter-indel-reads is overwritten!");
		if (D.check_orientation) vrb.warning("--check-orientation is overwritten!");
		D.legacy_options = true;
		D.param_min_baseQ = 13;
		D.max_depth = 8000;
		D.param_dup_rd = true;
		D.check_proper_pair = true;
		D.keep_orphan = false;
		D.keep_failqc = false;
		D.param_rm_indel = false;
		D.check_orientation = false;
	}


	vrb.bullet("Mapping quality >= " + stb.str(D.param_min_mapQ));
	vrb.bullet("Base quality >= " + stb.str(D.param_min_baseQ));
	vrb.bullet("Coverage ASE >= " + stb.str(D.param_min_cov));
	vrb.bullet("Coverage for REF bias >= " + stb.str(D.param_min_cov_for_ref_alt));
	vrb.bullet("Binomial p-value threshold = " + stb.str(D.param_min_pval));
	vrb.bullet("Genotype probability >= " + stb.str(D.param_min_gp));
	vrb.bullet("Imputation quality >= " + stb.str(D.param_min_iq));
	vrb.bullet("Both alleles seen for ASE = " + stb.str(D.param_both_alleles_seen));
	vrb.bullet("Both alleles seen for REF bias = " + stb.str(D.param_both_alleles_seen_bias));
	vrb.bullet("Remove duplicate reads = " + stb.str(D.param_dup_rd));
	vrb.bullet("Remove indel reads = " + stb.str(D.param_rm_indel));
	vrb.bullet("Keep failed qc reads = " + stb.str(D.keep_failqc));
	vrb.bullet("Keep orphan reads = " + stb.str(D.keep_orphan));
	vrb.bullet("Check proper pairing = " + stb.str(D.check_proper_pair));
	vrb.bullet("Subsample above this percentile for REF bias = " + stb.str(D.param_sample));
	vrb.bullet("Max depth for pileup = " + stb.str(D.max_depth));

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
	if (D.options.count("fasta")) D.readGenome(D.options["fasta"].as < string > ());
	D.readGenotypes2(D.options["vcf"].as < string > (), D.options["reg"].as < string > (), D.options["filtered"].as < string > ());
	if (D.region_length) D.getRegions();
	D.readSequences(D.options["bam"].as < string > ());
	D.calculateRefToAltBias(D.options["filtered"].as < string > ());
	D.calculateASE(D.options["out"].as < string > (), D.options["filtered"].as < string > ());
}
