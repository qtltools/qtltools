/*
 * quan2_data.h
 *
 *  Created on: Aug 21, 2017
 *      Author: halit
 */

#ifndef SRC_MODE_QUAN2_QUAN2_DATA_H_
#define SRC_MODE_QUAN2_QUAN2_DATA_H_

#define QUAN_VERSION 1

#include "../common/data.h"
#include <functional>


enum FILTER { PASS, UNMAP, SECD, FAILQC , DUP, NPP, MAPQ, STRAND, MISM};

class my_basic_overlap{
public:
	unsigned int start,end;
	my_basic_overlap(){start=0;end=0;}
	my_basic_overlap(unsigned int s , unsigned int e){start=s;end=e;}
	bool isSet(){return (start!=0 || end!=0);}
	unsigned int overlap(unsigned int s , unsigned int e) {return isSet() && start <= e && end >= s ? min(end,e) - max(start,s) + 1 : 0;}
};

class my_stats{
public:
	unsigned long int mapQ,mismatch,unmapped,unpaired,dup,notexon,good,total,exonicint,exonicint_multi,secondary,failqc,merged;
    double exonic,exonic_multi;
    my_stats(){mapQ=0;mismatch=0;unmapped=0;dup=0;unpaired=0;notexon=0;exonic=0.0;total=0;good=0;exonicint=0;secondary=0;failqc=0;exonic_multi=0.0;exonicint_multi=0;merged=0;}
};


class my_cont_blocks{
public:
	string chr;
    vector < unsigned int > starts;
    vector < unsigned int > ends;
    vector < unsigned int > lengths;
    vector < vector <my_basic_overlap> > base_count;
    unsigned int read_length;
    unsigned int mmc;
    FILTER filter;
    bam1_core_t core;
    string name;
    bool merged;
    my_cont_blocks(){read_length=0;mmc=0;filter=PASS;name="";chr="NA";merged=false;}
    void merge(my_cont_blocks &B){
        for (int i =0 ; i < starts.size(); i++){
            for (int j =i; j < B.starts.size() && ends[i] >= B.starts[j] ; j++){
                if(starts[i] <= B.ends[j] && ends[i] >= B.starts[j]){
                	B.merged = true;
                	merged=true;
                	unsigned int s = max(starts[i], B.starts[j]);
                	unsigned int e = min(ends[i], B.ends[j]);
                    base_count[i].push_back(my_basic_overlap(s,e));
                    B.base_count[j].push_back(my_basic_overlap(s,e));
#ifdef DEBUG
                	cerr<< "MERGING " << name << " " << starts[i] << "," << ends[i] << " and " << B.starts[j] << "," << B.ends[j]  << " to " << s << " " << e << endl;
#endif
                }
            }
        }
    }
    friend ostream& operator<<(ostream& out, my_cont_blocks& g){
        out<< g.name << "\t" << g.read_length << "\t" << g.mmc << "\t" << g.core.pos << "\t" << g.core.mpos << "\t" << g.base_count.size() <<"\t" << g.chr;
        for (int i =0; i < g.starts.size(); i++ ) out << "\t" << g.starts[i] << "," << g.ends[i];
        //out << endl;
        return out;
    }
};

class filters{
public:
	int max_intron_length;
    double max_mismatch_count,max_mismatch_count_total;
    unsigned int min_mapQ,max_read_length,min_exon;
    bool dup_remove,proper_pair,check_consistency,debug,merge,fraction_mm,fraction_mmt,fail_qc,old_wrong_split;
    filters(){
    	max_intron_length=0;
    	max_mismatch_count=0.0;
    	max_mismatch_count_total=0.0;
    	min_mapQ=0;max_read_length=0;
    	dup_remove=false;
    	proper_pair=false;
    	check_consistency=false;
    	debug=false;
    	merge=true;
    	fraction_mm=false;
    	fraction_mmt=false;
    	fail_qc=false;
    	old_wrong_split=false;
    	min_exon=0;
    }
};

class my_block{
public:
	string chr,gene_name,gene_id,name;
    unsigned int start;
    unsigned int end;
    unsigned int length;
    short int strand;
    double read_count;
    my_block(){start=0,end=0,length=0,strand=0;read_count=0.0;chr="NA",gene_name="NA";gene_id="NA",name="NA";}
    my_block(string c , unsigned int s, unsigned int e){start=s,end=e,length=0,strand=0;read_count=0.0;chr=c,gene_name="NA";gene_id="NA",name="NA";}

    bool operator < (my_block const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (start < p.start) return true;
        if (start > p.start) return false;
        if (end < p.end) return true;
        return false;
    }

    friend ostream& operator<<(ostream& out, my_block& g){
        out << g.chr << "\t" << g.start << "\t" << g.end;
        return out;
    }

};

class my_exon: public my_block{
public:
    string gene_type;
    bool merged;
    my_exon(){chr="NA";name="NA";start=0;end=0;gene_id="NA";gene_name="NA";gene_type="NA";strand = 0;;merged=false;length=0;read_count=0.0;}
    my_exon(string c , unsigned int s, unsigned int e, string gi, string gn, string gt, string st){chr=c;start=s;end=e;gene_id=gi;name = gi + "_" + stb.str(start) + "_" + stb.str(end);gene_name = gn; gene_type =gt; strand = st == "-" ? -1 : 1; length = end - start + 1; merged=false;read_count=0.0;}

    void resize(unsigned int s, unsigned int e){start = s; end = e;length = end - start + 1; name = gene_id + "_" + stb.str(start) + "_" + stb.str(end);merged=true;}
    void merge(vector < my_exon > &exs){
        unsigned int sr = start;
        unsigned int er = end;
        for (int i =0 ; i < exs.size(); i++){
            if (exs[i].start < sr) sr = exs[i].start;
            if (exs[i].end > er) er = exs[i].end;
        }
        resize(sr,er);
    }

    bool overlap(unsigned int s, unsigned int e){
    	if (s <= end && e >= start) return true;
    	else return false;
    }

    bool overlap(string c, unsigned int s, unsigned int e){
    	if (c == chr && s <= end && e >= start) return true;
    	else return false;
    }

    bool operator < (my_exon const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (start < p.start) return true;
        if (start > p.start) return false;
        if (end < p.end) return true;
        return false;
    }

    friend ostream& operator<<(ostream& out, const my_exon& p){
        out << p.chr << "\t" << p.start << "\t" << p.end << "\t" << p.name << "\t" << p.length << "\t" << (p.strand == 1 ? "+" : "-") << "\t"  << p.merged << endl;
        return out;
    }
};

class my_gene: public my_block{
public:
    string region;
    unsigned int tss;
    vector < my_exon > exons;
    my_gene(){chr="NA";name="NA";start=0;end=0;gene_id="NA";strand = 0;region="NA"; length = 0 ;tss=0;read_count=0.0;gene_name = name;};
    ~my_gene(){exons.clear();}

    void assign(my_exon &e){
        if (strand && e.strand != strand) vrb.error("Strand mismatch");
        if (chr != "NA" && e.chr != chr) vrb.error("Chr mismatch");
        if (gene_name != "NA" && e.gene_name != gene_name) vrb.error("Name mismatch");
        if ( gene_id != "NA" && e.gene_id != gene_id) vrb.error("ID mismatch");
        strand = e.strand;
        gene_name = e.gene_name;
        chr = e.chr;
        gene_id = e.gene_id;
        if (start == 0 || e.start < start ) start = e.start;
        if (end == 0 || e.end > end) end = e.end;
        tss = strand == -1 ? end : start;
        vector < my_exon > overlaping;
        vector < int > overlaping_idx;
        unsigned int tot = 0;
        sort(exons.begin(), exons.end());
        for (int i =0 ; i < exons.size() ; i++){
            if (exons[i].end >= e.start && exons[i].start <= e.end ){
                overlaping.push_back(exons[i]);
                overlaping_idx.push_back(i);
                tot += exons[i].length;
            }
        }
        if (overlaping.size()){
            e.merge(overlaping);
            exons.erase(exons.begin() + overlaping_idx[0], exons.begin() + (overlaping_idx.back() + 1));
            length -= tot;
        }
        exons.push_back(e);
        length += e.length;
        region = chr + ":" + stb.str(start) + "-" + stb.str(end);
    }

    bool operator < (my_gene const & p) const {
        if (chr < p.chr) return true;
        if (chr > p.chr) return false;
        if (start < p.start) return true;
        if (start > p.start) return false;
        if (end < p.end) return true;
        return false;
    }

    bool overlap(unsigned int s, unsigned int e){
    	if (s <= end && e >= start) return true;
    	else return false;
    }

    bool overlap(string c, unsigned int s, unsigned int e){
    	if (c == chr && s <= end && e >= start) return true;
    	else return false;
    }

    bool overlap(genomic_region &gr){
        if (gr.chr == chr && gr.start <= end && gr.end >= start) return true;
        else return false;
    }

    friend ostream& operator<<(ostream& out, my_gene& g){
        out << g.chr << "\t" << g.start << "\t" << g.end << "\t" << g.gene_id << "\t" << g.length << "\t" << (g.strand == 1 ? "+" : "-") << "\t" << 2<<endl;
        //sort(g.exons.begin(), g.exons.end());
        for (int i =0; i < g.exons.size(); i++ ) out << g.exons[i];
        return out;
    }

};

class quan2_data : public data {
public:
    map < string, unsigned int > genes_map;
    vector < my_gene > genes;
    string bam;
    my_stats  stats;
    string  sample;
    set < string > gene_types;
    filters filter;
    genomic_region region;
    string hash;

    map < string , map < int , vector < int > > > genome;

    static const int binsize = 100000;

    bool setRegion(string reg) {return region.parse(reg);}
    inline vector <my_gene> getOverlappingGenesQ(my_block &);
    void readGTF(string);
    void readBam(string);
    my_cont_blocks keep_read(bam1_t *b);
    void printBEDcount(string);
    void printBEDrpkm(string);
    void printBEDtpm(string);
    void printStats(string);
    unsigned long long int fnv1a_hash (string &str){
        unsigned long long int h = 14695981039346656037ull;
        for (int i = 0; i < str.length(); i++) {
            h = (h ^ (unsigned char) str[i]) * 1099511628211ull;
        }
        return h;
    }
    string convertToBase(unsigned long long int, unsigned int = 62);
};

void quan2_main(vector < string > & );

/*
   xxHash - Extremely Fast Hash algorithm
   Header File
   Copyright (C) 2012-2016, Yann Collet.
   BSD 2-Clause License (http://www.opensource.org/licenses/bsd-license.php)
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are
   met:
       * Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
       * Redistributions in binary form must reproduce the above
   copyright notice, this list of conditions and the following disclaimer
   in the documentation and/or other materials provided with the
   distribution.
   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
   OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
   SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
   LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
   DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
   THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
   (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
   You can contact the author at :
   - xxHash source repository : https://github.com/Cyan4973/xxHash
*/

/* Notice extracted from xxHash homepage :
xxHash is an extremely fast Hash algorithm, running at RAM speed limits.
It also successfully passes all tests from the SMHasher suite.
Comparison (single thread, Windows Seven 32 bits, using SMHasher on a Core 2 Duo @3GHz)
Name            Speed       Q.Score   Author
xxHash          5.4 GB/s     10
CrapWow         3.2 GB/s      2       Andrew
MumurHash 3a    2.7 GB/s     10       Austin Appleby
SpookyHash      2.0 GB/s     10       Bob Jenkins
SBox            1.4 GB/s      9       Bret Mulvey
Lookup3         1.2 GB/s      9       Bob Jenkins
SuperFastHash   1.2 GB/s      1       Paul Hsieh
CityHash64      1.05 GB/s    10       Pike & Alakuijala
FNV             0.55 GB/s     5       Fowler, Noll, Vo
CRC32           0.43 GB/s     9
MD5-32          0.33 GB/s    10       Ronald L. Rivest
SHA1-32         0.28 GB/s    10
Q.Score is a measure of quality of the hash function.
It depends on successfully passing SMHasher test set.
10 is a perfect score.
A 64-bits version, named XXH64, is available since r35.
It offers much better speed, but for 64-bits applications only.
Name     Speed on 64 bits    Speed on 32 bits
XXH64       13.8 GB/s            1.9 GB/s
XXH32        6.8 GB/s            6.0 GB/s
*/

#ifndef XXHASH_H_5627135585666179
#define XXHASH_H_5627135585666179 1

#if defined (__cplusplus)
extern "C" {
#endif


/* ****************************
*  Definitions
******************************/
#include <stddef.h>   /* size_t */
typedef enum { XXH_OK=0, XXH_ERROR } XXH_errorcode;


/* ****************************
*  API modifier
******************************/
/** XXH_PRIVATE_API
*   This is useful to include xxhash functions in `static` mode
*   in order to inline them, and remove their symbol from the public list.
*   Methodology :
*     #define XXH_PRIVATE_API
*     #include "xxhash.h"
*   `xxhash.c` is automatically included.
*   It's not useful to compile and link it as a separate module.
*/
#ifdef XXH_PRIVATE_API
#  ifndef XXH_STATIC_LINKING_ONLY
#    define XXH_STATIC_LINKING_ONLY
#  endif
#  if defined(__GNUC__)
#    define XXH_PUBLIC_API static __inline __attribute__((unused))
#  elif defined (__cplusplus) || (defined (__STDC_VERSION__) && (__STDC_VERSION__ >= 199901L) /* C99 */)
#    define XXH_PUBLIC_API static inline
#  elif defined(_MSC_VER)
#    define XXH_PUBLIC_API static __inline
#  else
#    define XXH_PUBLIC_API static   /* this version may generate warnings for unused static functions; disable the relevant warning */
#  endif
#else
#  define XXH_PUBLIC_API   /* do nothing */
#endif /* XXH_PRIVATE_API */

/*!XXH_NAMESPACE, aka Namespace Emulation :
If you want to include _and expose_ xxHash functions from within your own library,
but also want to avoid symbol collisions with other libraries which may also include xxHash,
you can use XXH_NAMESPACE, to automatically prefix any public symbol from xxhash library
with the value of XXH_NAMESPACE (therefore, avoid NULL and numeric values).
Note that no change is required within the calling program as long as it includes `xxhash.h` :
regular symbol name will be automatically translated by this header.
*/
#ifdef XXH_NAMESPACE
#  define XXH_CAT(A,B) A##B
#  define XXH_NAME2(A,B) XXH_CAT(A,B)
#  define XXH_versionNumber XXH_NAME2(XXH_NAMESPACE, XXH_versionNumber)
#  define XXH32 XXH_NAME2(XXH_NAMESPACE, XXH32)
#  define XXH32_createState XXH_NAME2(XXH_NAMESPACE, XXH32_createState)
#  define XXH32_freeState XXH_NAME2(XXH_NAMESPACE, XXH32_freeState)
#  define XXH32_reset XXH_NAME2(XXH_NAMESPACE, XXH32_reset)
#  define XXH32_update XXH_NAME2(XXH_NAMESPACE, XXH32_update)
#  define XXH32_digest XXH_NAME2(XXH_NAMESPACE, XXH32_digest)
#  define XXH32_copyState XXH_NAME2(XXH_NAMESPACE, XXH32_copyState)
#  define XXH32_canonicalFromHash XXH_NAME2(XXH_NAMESPACE, XXH32_canonicalFromHash)
#  define XXH32_hashFromCanonical XXH_NAME2(XXH_NAMESPACE, XXH32_hashFromCanonical)
#  define XXH64 XXH_NAME2(XXH_NAMESPACE, XXH64)
#  define XXH64_createState XXH_NAME2(XXH_NAMESPACE, XXH64_createState)
#  define XXH64_freeState XXH_NAME2(XXH_NAMESPACE, XXH64_freeState)
#  define XXH64_reset XXH_NAME2(XXH_NAMESPACE, XXH64_reset)
#  define XXH64_update XXH_NAME2(XXH_NAMESPACE, XXH64_update)
#  define XXH64_digest XXH_NAME2(XXH_NAMESPACE, XXH64_digest)
#  define XXH64_copyState XXH_NAME2(XXH_NAMESPACE, XXH64_copyState)
#  define XXH64_canonicalFromHash XXH_NAME2(XXH_NAMESPACE, XXH64_canonicalFromHash)
#  define XXH64_hashFromCanonical XXH_NAME2(XXH_NAMESPACE, XXH64_hashFromCanonical)
#endif


/* *************************************
*  Version
***************************************/
#define XXH_VERSION_MAJOR    0
#define XXH_VERSION_MINOR    6
#define XXH_VERSION_RELEASE  3
#define XXH_VERSION_NUMBER  (XXH_VERSION_MAJOR *100*100 + XXH_VERSION_MINOR *100 + XXH_VERSION_RELEASE)
XXH_PUBLIC_API unsigned XXH_versionNumber (void);


/*-**********************************************************************
*  32-bits hash
************************************************************************/
typedef unsigned int XXH32_hash_t;

/*! XXH32() :
    Calculate the 32-bits hash of sequence "length" bytes stored at memory address "input".
    The memory between input & input+length must be valid (allocated and read-accessible).
    "seed" can be used to alter the result predictably.
    Speed on Core 2 Duo @ 3 GHz (single thread, SMHasher benchmark) : 5.4 GB/s */
XXH_PUBLIC_API XXH32_hash_t XXH32 (const void* input, size_t length, unsigned int seed);

/*======   Streaming   ======*/
typedef struct XXH32_state_s XXH32_state_t;   /* incomplete type */
XXH_PUBLIC_API XXH32_state_t* XXH32_createState(void);
XXH_PUBLIC_API XXH_errorcode  XXH32_freeState(XXH32_state_t* statePtr);
XXH_PUBLIC_API void XXH32_copyState(XXH32_state_t* dst_state, const XXH32_state_t* src_state);

XXH_PUBLIC_API XXH_errorcode XXH32_reset  (XXH32_state_t* statePtr, unsigned int seed);
XXH_PUBLIC_API XXH_errorcode XXH32_update (XXH32_state_t* statePtr, const void* input, size_t length);
XXH_PUBLIC_API XXH32_hash_t  XXH32_digest (const XXH32_state_t* statePtr);

/*
These functions generate the xxHash of an input provided in multiple segments.
Note that, for small input, they are slower than single-call functions, due to state management.
For small input, prefer `XXH32()` and `XXH64()` .
XXH state must first be allocated, using XXH*_createState() .
Start a new hash by initializing state with a seed, using XXH*_reset().
Then, feed the hash state by calling XXH*_update() as many times as necessary.
Obviously, input must be allocated and read accessible.
The function returns an error code, with 0 meaning OK, and any other value meaning there is an error.
Finally, a hash value can be produced anytime, by using XXH*_digest().
This function returns the nn-bits hash as an int or long long.
It's still possible to continue inserting input into the hash state after a digest,
and generate some new hashes later on, by calling again XXH*_digest().
When done, free XXH state space if it was allocated dynamically.
*/

/*======   Canonical representation   ======*/

typedef struct { unsigned char digest[4]; } XXH32_canonical_t;
XXH_PUBLIC_API void XXH32_canonicalFromHash(XXH32_canonical_t* dst, XXH32_hash_t hash);
XXH_PUBLIC_API XXH32_hash_t XXH32_hashFromCanonical(const XXH32_canonical_t* src);

/* Default result type for XXH functions are primitive unsigned 32 and 64 bits.
*  The canonical representation uses human-readable write convention, aka big-endian (large digits first).
*  These functions allow transformation of hash result into and from its canonical format.
*  This way, hash values can be written into a file / memory, and remain comparable on different systems and programs.
*/


#ifndef XXH_NO_LONG_LONG
/*-**********************************************************************
*  64-bits hash
************************************************************************/
typedef unsigned long long XXH64_hash_t;

/*! XXH64() :
    Calculate the 64-bits hash of sequence of length "len" stored at memory address "input".
    "seed" can be used to alter the result predictably.
    This function runs faster on 64-bits systems, but slower on 32-bits systems (see benchmark).
*/
XXH_PUBLIC_API XXH64_hash_t XXH64 (const void* input, size_t length, unsigned long long seed);

/*======   Streaming   ======*/
typedef struct XXH64_state_s XXH64_state_t;   /* incomplete type */
XXH_PUBLIC_API XXH64_state_t* XXH64_createState(void);
XXH_PUBLIC_API XXH_errorcode  XXH64_freeState(XXH64_state_t* statePtr);
XXH_PUBLIC_API void XXH64_copyState(XXH64_state_t* dst_state, const XXH64_state_t* src_state);

XXH_PUBLIC_API XXH_errorcode XXH64_reset  (XXH64_state_t* statePtr, unsigned long long seed);
XXH_PUBLIC_API XXH_errorcode XXH64_update (XXH64_state_t* statePtr, const void* input, size_t length);
XXH_PUBLIC_API XXH64_hash_t  XXH64_digest (const XXH64_state_t* statePtr);

/*======   Canonical representation   ======*/
typedef struct { unsigned char digest[8]; } XXH64_canonical_t;
XXH_PUBLIC_API void XXH64_canonicalFromHash(XXH64_canonical_t* dst, XXH64_hash_t hash);
XXH_PUBLIC_API XXH64_hash_t XXH64_hashFromCanonical(const XXH64_canonical_t* src);
#endif  /* XXH_NO_LONG_LONG */


#ifdef XXH_STATIC_LINKING_ONLY

/* ================================================================================================
   This section contains definitions which are not guaranteed to remain stable.
   They may change in future versions, becoming incompatible with a different version of the library.
   They shall only be used with static linking.
   Never use these definitions in association with dynamic linking !
=================================================================================================== */

/* These definitions are only meant to make possible
   static allocation of XXH state, on stack or in a struct for example.
   Never use members directly. */

struct XXH32_state_s {
   unsigned total_len_32;
   unsigned large_len;
   unsigned v1;
   unsigned v2;
   unsigned v3;
   unsigned v4;
   unsigned mem32[4];   /* buffer defined as U32 for alignment */
   unsigned memsize;
   unsigned reserved;   /* never read nor write, will be removed in a future version */
};   /* typedef'd to XXH32_state_t */

#ifndef XXH_NO_LONG_LONG   /* remove 64-bits support */
struct XXH64_state_s {
   unsigned long long total_len;
   unsigned long long v1;
   unsigned long long v2;
   unsigned long long v3;
   unsigned long long v4;
   unsigned long long mem64[4];   /* buffer defined as U64 for alignment */
   unsigned memsize;
   unsigned reserved[2];          /* never read nor write, will be removed in a future version */
};   /* typedef'd to XXH64_state_t */
#endif

#ifdef XXH_PRIVATE_API
#  include "xxhash.c"   /* include xxhash function bodies as `static`, for inlining */
#endif

#endif /* XXH_STATIC_LINKING_ONLY */


#if defined (__cplusplus)
}
#endif

#endif /* XXHASH_H_5627135585666179 */


#endif /* SRC_MODE_QUAN2_QUAN2_DATA_H_ */
