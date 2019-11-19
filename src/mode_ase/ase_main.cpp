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
		("vcf,v", boost::program_options::value< string >(), "Genotypes in VCF/BCF format. (REQUIRED and RECOMMENED to use a BCF file for performance)")
		("bam,b", boost::program_options::value< string >(), "Sequence data in BAM/SAM format sorted by position. (REQUIRED)")
		("fasta,f", boost::program_options::value< string >(), "Genome sequence in FASTA format. (RECOMMENED)")
		("ind,i", boost::program_options::value< string >(), "Sample to be processed. (REQUIRED)")
		("reg,r", boost::program_options::value< string >()->default_value(""), "Genomic region(s) to be processed.")
		("blacklist,B", boost::program_options::value< string >(), "BED file for blacklisted regions. RECOMMENED to use to filter out low mappability regions.")
		("merge-on-the-fly,t", "Merges the blacklisted regions on the fly. This can reduce memory usage if there are many overlapping or contiguous regions in the blacklist.")
		("gtf,g", boost::program_options::value< string >(), "Annotation in GTF format. If provided variants will be matched to exons. (RECOMMENDED)")
		("fix-chr,F", "Attempt to match chromosome names to the BAM file by adding or removing chr to chromosome names. Does not apply to --{include,exclude}-positions options. These should be in the VCF chromosome names.")
		("fix-id,R", "Convert missing VCF variant IDs to chr_pos_refalt")
		("auto-flip,x", "Attempt to fix reference allele mismatches. Requires a fasta file for the reference sequence. (NOT RECOMMENDED)")
		("print-stats,P", "Print out stats for the filtered reads for ASE sites.")
		("suppress-warnings,k", "Suppress the warnings about individual variants.")
		("illumina13,j", "Base quality is in the Illumina-1.3+ encoding")
		("group-by,G", boost::program_options::value< int >()->default_value(0,"OFF"), "Group variants separated by this much into batches. This allows you not to stream the whole BAM file and may improve running time.")
		("max-depth,d", boost::program_options::value< int >()->default_value(65536,"65536"), "Pileup max-depth. Set to 0 if you want maximum but this will be slower and use more memory.")
		("filtered,l", boost::program_options::value< string >()->default_value(""), "File to output filtered variants. RECOMMENED for troubleshooting especially if --suppress-warnings.")
		("out,o", boost::program_options::value< string >(), "Output file. (REQUIRED)");

	boost::program_options::options_description opt_parameters ("\x1B[32mFilters\33[0m");
	opt_parameters.add_options()
		("mapq,q", boost::program_options::value< int >(), "Minimum mapping quality for a read to be considered. (REQUIRED)")
		("baseq,Q", boost::program_options::value< int >()->default_value(10), "Minimum phred quality for a base to be considered.")
		("pvalue,p", boost::program_options::value< double >()->default_value(1.0, "1.0"), "Binomial p-value threshold for ASE output.")
		("cov,c", boost::program_options::value< int >()->default_value(16), "Minimum coverage for a genotype to be considered in ASE analysis.")
		("cov-bias,C", boost::program_options::value< int >()->default_value(10), "Minimum coverage for a genotype to be considered in REF bias.")
		("sites,s", boost::program_options::value< int >()->default_value(200), "Minimum number of sites to calculate a REF bias from for a specific REF/ALT pair. The REF bias for pairs with less than this many sites will be calculated from all sites.")
		("imp-qual-id,I", boost::program_options::value< string >()->default_value("INFO", "INFO"), "The INFO ID of the imputation score in the VCF.")
		("geno-prob-id,L", boost::program_options::value< string >()->default_value("GP", "GP"), "The FORMAT ID of the genotype posterior probabilities for RR/RA/AA  in the VCF.")
		("imp-qual,W", boost::program_options::value< double >()->default_value(0.0, "0.0"), "Minimum imputation score for a variant to be considered.")
		("geno-prob,V", boost::program_options::value< double >()->default_value(0.0, "0.0"), "Minimum posterior probability for a genotype to be considered.")
		("subsample,S", boost::program_options::value< double >()->default_value(0.75, "0.75"), "Randomly subsample sites that have greater coverage than this percentile of all the sites in REF bias calculations. Set to 1 to turn off which is NOT RECOMMENED.")
		("both-alleles-seen,a", "Require both alleles to be observed in RNA-seq reads for a site for ASE calculations.")
		("keep-bans-for-bias,A", "DON'T require both alleles to be observed in RNA-seq reads for a site for REF mapping bias calculations. (NOT RECOMMENDED)")
		("keep-discordant-for-bias,E", "If given sites with more discordant alleles than REF or ALT alleles will be included in the REF bias calculations. (NOT RECOMMENDED)")
		("filter-indel-reads,D", "Remove reads that contain indels.")
		("keep-failed-qc,e", "Keep fastq reads that fail sequencing QC (as indicated by the sequencer).")
		("keep-orphan-reads,O", "Keep paired end reads where one of mates is unmapped.")
		("check-proper-pairing,y", "If provided only properly paired reads according to the aligner will be considered.")
		("ignore-orientation,X", "If NOT provided only mate pairs where both mates are on the same chromosome and where the first mate is on the +ve strand and the second is on the -ve strand will be considered. (NOT RECOMMENED)")
		("filter-duplicates,u", "Remove duplicate sequencing reads in the process. (NOT RECOMMENED)")
		("legacy-options,J", "Replicate legacy options used. (NOT RECOMMENED).");

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
	if (!D.options.count("mapq")) vrb.error("Mapping quality threshold should be set to uniquely mapping reads with --mapq [quality]");

	//CHECK PARAMETER VALUES
	D.max_depth = D.options["max-depth"].as < int > ();
	if (D.max_depth < 0) vrb.error("--max-depth cannot be less than 0!");
	if(D.max_depth == 0) {
		D.max_depth = INT_MAX;
		vrb.warning("Setting pileup max-depth to INT_MAX!");
	}
	if (D.max_depth > 1000000) vrb.warning("Pileup max-depth is above 1M. Potential memory hog!");

	if (D.options.count("fasta") == 0 && D.options.count("auto-flip")) vrb.warning("Ignoring --auto-flip since no --fasta is provided!");
	int ot = D.options["mapq"].as < int > ();
	if (ot < 0) vrb.error("--mapq cannot be negative!");
	D.param_min_mapQ = ot;
	ot =  D.options["baseq"].as < int > ();
	if (ot < 0) vrb.error("--baseq cannot be negative!");
	D.param_min_baseQ = ot;
	D.param_min_pval = D.options["pvalue"].as < double > ();
	if (D.param_min_pval < 0.0 || D.param_min_pval > 1.0) vrb.error("--pvalue must be between 0 and 1!");
	ot = D.options["cov"].as < int > ();
	if (ot < 0) vrb.error("--cov cannot be negative!");
	D.param_min_cov = ot;
	ot =  D.options["cov-bias"].as < int > ();
	if (ot < 0) vrb.error("--cov-bias cannot be negative!");
	D.param_min_cov_for_ref_alt = ot;
	ot = D.options["sites"].as < int > ();
	if (ot < 0) vrb.error("--sites cannot be negative!");
	D.param_min_sites_for_ref_alt = ot;
	D.param_imputation_score_label = D.options["imp-qual-id"].as <string> ();
	D.param_genotype_likelihood_label = D.options["geno-prob-id"].as <string> ();
	D.param_min_gp = D.options["geno-prob"].as < double > ();
	if (D.param_min_gp < 0.0) vrb.error("--geno-prob cannot be negative!");
	D.param_min_iq = D.options["imp-qual"].as < double > ();
	if (D.param_min_iq < 0.0) vrb.error("--imp-qual cannot be negative!");
	D.param_sample = D.options["subsample"].as < double > ();
	if (D.param_sample < 0.0 || D.param_sample > 1.0) vrb.error("--subsample must be between 0 and 1!");
	D.param_both_alleles_seen = D.options.count("both-alleles-seen");
	D.param_both_alleles_seen_bias = (D.options.count("keep-bans-for-bias") == 0);
	D.param_rm_indel = D.options.count("filter-indel-reads");
	D.keep_failqc = D.options.count("keep-failed-qc");
	D.keep_orphan = D.options.count("keep-orphan-reads");
	D.check_proper_pair = D.options.count("check-proper-pairing");
	D.check_orientation = (D.options.count("ignore-orientation") == 0);
	D.param_dup_rd = D.options.count("filter-duplicates");
	D.fix_chr = D.options.count("fix-chr");
	D.fix_id = D.options.count("fix-id");
	D.auto_flip = D.options.count("auto-flip");
	D.print_stats = D.options.count("print-stats");
	D.print_warnings = (D.options.count("suppress-warnings") == 0);
	D.illumina13 = D.options.count("illumina13");
	D.on_the_fly = D.options.count("merge-on-the-fly");
	D.keep_discordant = D.options.count("keep-discordant-for-bias");
	ot = D.options["group-by"].as < int > ();
	if (ot < 0) vrb.error("--group-by cannot be negative!");
	D.region_length = ot;

	if(D.options.count("legacy-options")){
		if (D.param_min_baseQ != 13) vrb.warning("--baseq is overwritten!");
		if (D.max_depth != 8000) vrb.warning("--max-depth is overwritten!");
		if (!D.param_dup_rd) vrb.warning("--filter-duplicates is overwritten!");
		if (!D.check_proper_pair) vrb.warning("--check-proper-pairing is overwritten!");
		if (D.keep_orphan) vrb.warning("--keep-orphan-reads is overwritten!");
		if (D.keep_failqc) vrb.warning("--keep-failqc is overwritten!");
		if (!D.keep_discordant) vrb.warning("--keep-discordant-for-bias is overwritten!");
		if (D.param_rm_indel) vrb.warning("--filter-indel-reads is overwritten!");
		if (D.check_orientation) vrb.warning("--ignore-orientation is overwritten!");
		D.legacy_options = true;
		D.param_min_baseQ = 13;
		D.max_depth = 8000;
		D.param_dup_rd = true;
		D.check_proper_pair = true;
		D.keep_orphan = false;
		D.keep_failqc = false;
		D.param_rm_indel = false;
		D.check_orientation = false;
		D.keep_discordant = true;
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
	vrb.bullet("Keep sites with more discordant alleles than ref or alt for REF bias = " + stb.str(D.keep_discordant));
	vrb.bullet("Remove duplicate reads = " + stb.str(D.param_dup_rd));
	vrb.bullet("Remove indel reads = " + stb.str(D.param_rm_indel));
	vrb.bullet("Keep failed qc reads = " + stb.str(D.keep_failqc));
	vrb.bullet("Keep orphan reads = " + stb.str(D.keep_orphan));
	vrb.bullet("Check orientation = " + stb.str(D.check_orientation));
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
	D.compareChrs(D.options["vcf"].as < string > (), D.options["bam"].as < string > (), D.options["reg"].as < string > ());
	if (D.options.count("blacklist")) D.readBlacklist(D.options["blacklist"].as < string > ());
	if (D.options.count("fasta")) D.readGenome(D.options["fasta"].as < string > ());
	if (D.options.count("gtf")) D.readGTF(D.options["gtf"].as < string > ());
	D.readGenotypes(D.options["vcf"].as < string > (), D.options["filtered"].as < string > ());
	if (D.region_length) D.getRegions();
	else if (D.options["reg"].as < string > () == "") D.collapseRegions();
	D.readSequences(D.options["bam"].as < string > ());
	D.calculateRefToAltBias(D.options["filtered"].as < string > ());
	D.calculateASE(D.options["out"].as < string > (), D.options["filtered"].as < string > ());
}
