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

#include "pca_data.h"

void pca_main(vector < string > & argv) {
    pca_data D;
    
    //-------------------------
    // 1. DECLARE ALL OPTIONS
    //-------------------------
    D.declareBasicOptions();  //Mandatory
    boost::program_options::options_description opt_files ("\x1B[32mI/O\33[0m");
    opt_files.add_options()
    ("vcf", boost::program_options::value< string >(), "Genotypes in VCF/BCF/BED format.")
    ("bed", boost::program_options::value< string >(), "Phenotypes in BED format.")
    ("out", boost::program_options::value< string >(), "Output file prefix.");
    
    boost::program_options::options_description opt_parameters ("\x1B[32mParameters\33[0m");
    opt_parameters.add_options()
    ("center", "Center the quantifications or genotypes before PCA.")
    ("scale", "Scale the quantifications or genotypes to unit variance before PCA.")
    //("use-cor", "Use correlation rather than SVD in PCA (Valid only when number of samples is greater than number of phenotypes or genotypes.")
    ("maf", boost::program_options::value< double >()->default_value(0.0, "0"), "Exclude sites with MAF less than this.")
    ("distance", boost::program_options::value< unsigned int >()->default_value(0,"0"), "Only include sites separated with this many bp");
    
    D.option_descriptions.add(opt_files).add(opt_parameters);
    
    //-------------------
    // 2. PARSE OPTIONS
    //-------------------
    try {
        boost::program_options::store(boost::program_options::command_line_parser(argv).options(D.option_descriptions).run(), D.options);
        boost::program_options::notify(D.options);
    } catch ( const boost::program_options::error& e ) {
        cerr << "Error parsing [pca] command line :" << string(e.what()) << endl;
        exit(0);
    }
    
    //---------------------
    // 3. PRINT HELP/HEADER
    //---------------------
    vrb.ctitle("CONDUCT PCA ON DATA");
    if (D.options.count("help")) {
        cout << D.option_descriptions << endl;
        exit(EXIT_SUCCESS);
    }
    
    //-----------------
    // 4. COMMON CHECKS
    //-----------------
    if (!D.options.count("vcf") && !D.options.count("bed") ) vrb.error("Genotypes with --vcf [file.vcf] or phenotypes with --bed [file.bed] must be specified");
    if (D.options.count("vcf") && D.options.count("bed") ) vrb.error("Provide only one of --bed or --vcf");
    if (!D.options.count("out")) vrb.error("Output needs to be specified with --out [file.out]");
    if (D.options.count("bed") && (D.options["distance"].as < unsigned int > () > 0 || D.options["maf"].as < double > () > 0)) vrb.error("--distance and --maf cannot be combined with --bed");
    
    //--------------
    // 5. SET PARAMS
    //--------------
    if (D.options.count("vcf")){
        if (D.options["maf"].as < double > () < 0 || D.options["maf"].as < double > () > 1.0) vrb.error ("Incorrect --maf");
        D.maf_cutoff = D.options["maf"].as < double > ();
        vrb.bullet("MAF greater than " + stb.str(D.maf_cutoff));
        D.distance_separator = D.options["distance"].as < unsigned int > ();
        vrb.bullet("Sites every " + stb.str(D.distance_separator) + " bp");
    }
    string outFile = D.options["out"].as < string > ();
    
    //--------------
    // 6. READ FILES
    //--------------
    D.processBasicOptions();
    if (D.options.count("bed")) D.readSampleFromBED(D.options["bed"].as < string > ());										//Read samples in BED
    else {
        htsFile * fp = hts_open(D.options["vcf"].as < string > ().c_str(),"r");
        if (fp->format.format == sam) D.readSampleFromBED(D.options["vcf"].as < string > ());
        else D.readSampleFromVCF(D.options["vcf"].as < string > ());
        hts_close(fp);
    }
    D.mergeSampleLists();																		//Merge all sample lists
    if (D.options.count("bed")) D.readDataPhenoBED(D.options["bed"].as < string > ());										//Read data in BED
    else D.readData(D.options["vcf"].as < string > ());											//Read data in VCF
    
    //-----------------
    // 7. RUN ANALYSIS
    //-----------------
    D.imputeData();
    D.PCA.Calculate(D.options.count("use-cor"), D.options.count("center") , D.options.count("scale") );
    D.printPCA(outFile);
}