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

#include "quan_data.h"

void quan2_main(vector < string > & argv) {
    quan2_data D;
    //-------------------------
    // 1. DECLARE ALL OPTIONS
    //-------------------------
    D.declareBasicOptions();
    boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
    opt_files.add_options()
    	("gtf", boost::program_options::value< string >(), "Annotation in GTF format")
		("bam", boost::program_options::value< string >(), "Sequence data in BAM/SAM format.")
		("sample",boost::program_options::value< string > (), "Sample name. [Optional]")
        ("out-prefix", boost::program_options::value< string >(), "Output file prefix.");
		

	boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
	opt_parameters.add_options()
		("rpkm", "Output RPKM values.")
		("tpm", "Output TPM values.")
		("gene-types", boost::program_options::value< vector < string > > ()->multitoken(), "Gene types to quantify. (Requires gene_type attribute in GTF. It will also use transcript_type if present).")
        ("xxhash", "Rather than using GTF file name to generate unique hash for the options used, use hash of the GTF file.")
        ("no-hash", "Don't include a hash signifying the options used in the quantification in the file names.");

    boost::program_options::options_description opt_filters ("\x1B[32mFilters\33[0m");
    opt_filters.add_options()
    	("filter-mapping-quality", boost::program_options::value< unsigned int >()->default_value(10), "Minimum mapping quality for a read to be considered.")
		("filter-mismatch", boost::program_options::value< double >()->default_value(-1.0,"OFF"), "Maximum mismatches allowed in a read. If between 0 and 1 taken as the fraction of read length. (Requires NM attribute)")
		("filter-mismatch-total", boost::program_options::value< double >()->default_value(-1.0,"OFF"), "Maximum total mismatches allowed in paired reads. If between 0 and 1 taken as the fraction of combined read length. (Requires NM attribute)")
		("check-proper-pairing", "If provided only properly paired reads according to the aligner that are in correct orientation will be considered. Otherwise all pairs in correct orientation will be considered.")
        ("check-consistency", "If provided checks the consistency of split reads with annotation, rather than pure overlap of one of the blocks of the split read.")
        ("no-merge", "If provided overlapping mate pairs will not be merged. Default behaviour is to merge overlapping mate pairs based on the amount of overlap, such that each mate pair counts for less than 1 read.")
		("legacy-options", "Exactly replicate Dermitzakis lab original quantification script. (DO NOT USE UNLESS YOU HAVE A GOOD REASON). Sets --no-merge as well.")
		("filter-failed-qc", "Remove fastq reads that fail sequencing QC (as indicated by the sequencer).")
		("filter-min-exon", boost::program_options::value< unsigned int >()->default_value(0), "Minimal exon length to consider. Exons smaller than this will not be printed out in the exon quantifications, but will still count towards gene quantifications.")
		("filter-remove-duplicates", "Remove duplicate sequencing reads in the process (as indicated by the aligner.");

    boost::program_options::options_description opt_parallel ("\x1B[32mParallelization\33[0m");
    opt_parallel.add_options()
		("region", boost::program_options::value< string >(), "Region of interest.");

    D.option_descriptions.add(opt_files).add(opt_parameters).add(opt_filters).add(opt_parallel);

    //-------------------
    // 2. PARSE OPTIONS
    //-------------------
    boost::program_options::variables_map options;
    try {
        boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
        boost::program_options::notify(D.options);
    } catch ( const boost::program_options::error& e ) {
        cerr << "Error parsing [quan] command line :" << string(e.what()) << endl;
        exit(0);
    }

    //---------------------
    // 3. PRINT HELP/HEADER
    //---------------------
    vrb.ctitle("QUANTIFY GENES AND EXONS FROM BAM FILES");
    vrb.warning("OUTPUT IS NOT COMPATABLE WITH QUANTIFICATIONS GENERATED BEFORE VERSION 1.2\n");
    if (D.options.count("help")) {
        cout << D.option_descriptions << endl;
        exit(EXIT_SUCCESS);
    }

    //-----------------
    // 4. COMMON CHECKS
    //-----------------
    if (!D.options.count("gtf")) vrb.error("Gene annotation needs to be specified with --gtf [file.gtf]");
    if (!D.options.count("bam")) vrb.error("Sequence data needs to be specified with --bam [file.bam]");
    if (!D.options.count("out-prefix")) vrb.error("Output needs to be specified with --out-prefix [file.out]");
    if ((D.options.count("rpkm") || D.options.count("tpm"))  && D.options.count("region") ) vrb.error("Option --region and --rpkm or --tpm cannot be combined since we won't be parsing the whole BAM file");


    D.filter.min_mapQ = D.options["filter-mapping-quality"].as < unsigned int > ();
    vrb.bullet("Minimum mapping quality: " + stb.str(D.filter.min_mapQ));
    //D.filter.max_read_length = D.options["max-read-length"].as < unsigned int > ();
    //vrb.bullet("Maximum read length: " + stb.str(D.filter.max_read_length));
    double intpart;

    D.filter.max_mismatch_count_total = D.options["filter-mismatch-total"].as < double > ();
    if(D.filter.max_mismatch_count_total >= 0 && modf(D.filter.max_mismatch_count_total, &intpart) != 0.0) {
    	if(D.filter.max_mismatch_count_total > 1) vrb.error("--filter-mismatch-total cannot be greater than 1 when not an integer");
    	else D.filter.fraction_mmt = true;
    }
    vrb.bullet("Maximum mismatch count per mate-pair: " + stb.str(D.filter.max_mismatch_count_total));

    D.filter.max_mismatch_count = D.options["filter-mismatch"].as < double > ();
    if(D.filter.max_mismatch_count >= 0 && modf(D.filter.max_mismatch_count, &intpart) != 0.0) {
    	if(D.filter.max_mismatch_count > 1) vrb.error("--filter-mismatch cannot be greater than 1 when not an integer");
    	else D.filter.fraction_mm = true;
    }
    if (D.filter.max_mismatch_count < 0 && D.filter.max_mismatch_count_total >= 0 && !D.filter.fraction_mmt) D.filter.max_mismatch_count = D.filter.max_mismatch_count_total;
    vrb.bullet("Maximum mismatch count per read: " + stb.str(D.filter.max_mismatch_count));



    if (D.options.count("check-proper-pairing")){
        vrb.bullet("Checking properly paired flag");
        D.filter.proper_pair = true;
    }else vrb.bullet("Not checking properly paired flag");

    if (D.options.count("check-consistency")){
        vrb.bullet("Checking if all blocks of a split read are consistent with the annotation");
        D.filter.check_consistency = true;
    }else vrb.bullet("Not checking if all blocks of a split read are consistent with the annotation");

    if (D.options.count("filter-remove-duplicates")){
        vrb.bullet("Filtering reads flagged as duplicate");
        D.filter.dup_remove = true;
    }else vrb.bullet("Not filtering reads flagged as duplicate");

    if (D.options.count("filter-failed-qc")){
        vrb.bullet("Filtering reads flagged as failing QC");
        D.filter.fail_qc = true;
    }else vrb.bullet("Not filtering reads flagged as failing QC");

    if (D.options.count("no-merge")){
        vrb.bullet("Not merging overlapping mate pairs");
        D.filter.merge = false;
    }else if (!D.options.count("legacy-options")) vrb.bullet("Merging overlapping mate pairs");

    if (D.options.count("legacy-options")){
    	if (!D.options.count("no-merge")) vrb.bullet("Not merging overlapping mate pairs");
    	D.filter.min_exon = 2;
    	vrb.bullet("Excluding exons smaller than " + stb.str(D.filter.min_exon) );
        D.filter.old_wrong_split = true;
        D.filter.merge = false;
    }else{
    	D.filter.min_exon = D.options["filter-min-exon"].as < unsigned int > ();
    	vrb.bullet("Excluding exons smaller than " + stb.str(D.filter.min_exon) + " bp only in exon quantifications" );
    }
    ostringstream gts;
    if (D.options.count("gene-types")){
    	vector < string > t = D.options["gene-types"].as < vector < string > > ();
    	D.gene_types = set < string > (t.begin(),t.end());
        const char* const delim = " ";
        copy(D.gene_types.begin(), D.gene_types.end(), ostream_iterator<string>(gts, delim));
        vrb.bullet("Genes included: " + gts.str());
    }else vrb.bullet("Including all gene types");

    if(D.options.count("region")){
    	vrb.bullet("Region = [" + D.options["region"].as < string > () +"]");
    	if(!D.setRegion(D.options["region"].as < string >())) vrb.error("Cannot parse region!");
    }

    if(!D.options.count("no-hash")){
    	const char * gtff = D.options["gtf"].as < string > ().c_str();
        if (FILE *file = fopen(gtff, "r")) fclose(file);
        else vrb.error("Cannot open GTF file " + D.options["gtf"].as < string > ());
    	string r_p(realpath(gtff, NULL));
    	if (D.options.count("xxhash")){
    		vrb.bullet("Hashing: " + r_p);
    		XXH64_state_t* const state = XXH64_createState();
    		if (state==NULL) vrb.error("Cannot create XXH64 state! Try not using --xxhash");

    		size_t const bufferSize = 1024 * 64;
    		void * buffer = malloc(bufferSize);
    		if (buffer==NULL) vrb.error("Not enough memmory");

    		unsigned long long const seed = 0;   /* or any other value */
    		XXH_errorcode const resetResult = XXH64_reset(state, seed);
    		if (resetResult == XXH_ERROR) vrb.error("Cannot reset XXH64 state! Try not using --xxhash");;
    		ifstream fin(r_p, ifstream::binary);
    		int readCount;
    		while(fin.read( (char *) buffer, bufferSize) || (readCount = fin.gcount()) != 0 )  {
    			unsigned long long length = fin.gcount();
    			XXH_errorcode const addResult = XXH64_update(state, buffer, length);
    			if (addResult == XXH_ERROR) vrb.error("Cannot update XXH64 state! Try not using --xxhash");

    		}

    		unsigned long long const hash = XXH64_digest(state);
    		fin.close();
    		r_p = stb.str(hash);
    		vrb.bullet("64-bit xxhash is " + r_p);
    		free(buffer);
    		XXH64_freeState(state);
    	}
    	string opts_hash = stb.str(QUAN_VERSION) + "#" + r_p + "#" + stb.str(D.filter.min_mapQ) + "#" + stb.str(D.filter.max_mismatch_count) + "#" + stb.str(D.filter.max_mismatch_count_total) + "#" +
    					   stb.str(D.filter.proper_pair) + "#" + stb.str(D.filter.check_consistency) + "#" + stb.str(D.filter.dup_remove) + "#" + stb.str(D.filter.fail_qc) + "#" + stb.str(D.filter.min_exon) + "#" +
						   stb.str(D.filter.old_wrong_split) + "#" + gts.str() + "#" + D.region.get();
    	unsigned long long int h = D.fnv1a_hash(opts_hash);
    	D.hash = D.convertToBase(h);
    	vrb.bullet("Unique hash for this combination of options and GTF file: " + D.hash);
    }

    if (D.options.count("legacy-options")) vrb.warning("You are using --legacy-options, do you know what you are doing?");


    //------------------------------------------
    // 5. READ FILES / INITIALIZE / RUN ANALYSIS
    //------------------------------------------

    D.processBasicOptions();
    D.bam = D.options["bam"].as < string  > ();
    if (D.options.count("sample")) D.sample = D.options["sample"].as < string  > ();
    else D.sample = D.bam;
    D.readGTF(D.options["gtf"].as < string > ());
    D.readBam(D.bam);
    D.printBEDcount(D.options["out-prefix"].as < string > ());
    if (D.options.count("rpkm")) D.printBEDrpkm(D.options["out-prefix"].as < string > ());
    if (D.options.count("tpm")) D.printBEDtpm(D.options["out-prefix"].as < string > ());
    D.printStats(D.options["out-prefix"].as < string > ());
}
