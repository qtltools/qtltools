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

#include "cis_data.h"

void cis_data::writeHeader(string fout) {
	vrb.title("Writing header for output in [" + fout + "]");
	output_file fdo (fout);

	if (mode == CIS_PERM) {
		if (grp_mode == GRP_NONE) fdo << "phe_id";
		else fdo << "grp_id";
		fdo << " phe_chr phe_from phe_to phe_strd";
		switch (grp_mode) {
		case GRP_BEST: fdo << " phe_id n_phe_in_grp"; break;
		case GRP_PCA1: fdo << " ve_by_pc1 n_phe_in_grp"; break;
		case GRP_MEAN: fdo << " n_phe_in_grp"; break;
		}
		fdo << " n_var_in_cis dist_phe_var var_id var_chr var_from var_to dof1 dof2 bml1 bml2 nom_pval r_squared slope ";
		if (std_err) fdo<<"slope_se ";
		fdo << "adj_emp_pval adj_beta_pval" << endl;
	}
	if (mode == CIS_NOMI) {
		if (grp_mode == GRP_NONE) fdo << "phe_id";
		else fdo << "grp_id";
		fdo << " phe_chr phe_from phe_to phe_strd";
		switch (grp_mode) {
		case GRP_BEST: fdo << " phe_id n_phe_in_grp"; break;
		case GRP_PCA1: fdo << " ve_by_pc1 n_phe_in_grp"; break;
		case GRP_MEAN: fdo << " n_phe_in_grp"; break;
		}
		fdo << " n_var_in_cis dist_phe_var var_id var_chr var_from var_to nom_pval r_squared slope ";
		if(std_err) fdo << "slope_se ";
		fdo << "best_hit" << endl;
	}
	if (mode == CIS_COND) {
		if (grp_mode == GRP_NONE) fdo << "phe_id";
		else fdo << "grp_id";
		fdo << " phe_chr phe_from phe_to phe_strd";
		switch (grp_mode) {
		case GRP_BEST: fdo << " phe_id n_phe_in_grp"; break;
		case GRP_PCA1: fdo << " ve_by_pc1 n_phe_in_grp"; break;
		case GRP_MEAN: fdo << " n_phe_in_grp"; break;
		}
		fdo << " n_var_in_cis dist_phe_var var_id var_chr var_from var_to rank fwd_pval fwd_r_squared fwd_slope ";
		if(std_err) fdo << "fwd_slope_se ";
		fdo << "fwd_best_hit fwd_sig bwd_pval bwd_r_squared ";
		if(std_err) fdo << "bwd_slope_se ";
		fdo << "bwd_slope bwd_best_hit bwd_sig" << endl;
	}
	fdo.close();
}

bool cis_data::setPhenotypeRegion(string reg) {
	return regionPhenotype.parse(reg);
}

bool cis_data::setGenotypeRegion(string reg) {
	return regionGenotype.parse(reg);
}

class pgroup {
public:
	int start, end;
	string chr;

	pgroup(string pc, int ps, int pe) {
		chr = pc;
		start = ps;
		end = pe;
	}

	void merge(int ps, int pe) {
		if (start > ps) start = ps;
		if (end < pe) end = pe;
	}

	void merge(pgroup & p) {
		if (start > p.start) start = p.start;
		if (end < p.end) end = p.end;
	}

	bool overlap(pgroup & p) {
		if (chr != p.chr) return false;
		//cout << start << " " << end << " vs " << p.start << " " << p.end;
		if (start <= p.end && p.start <= end) {
			//cout << " Y" << endl;
			return true;
		} else {
			//cout << " N" << endl;
			return false;
		}
	}

	bool operator < (pgroup const & p) const {
		if (chr < p.chr) return true;
		if (chr > p.chr) return false;
		if (start < p.start) return true;
		if (start >= p.start) return false;
		return false;
	}
};

void cis_data::setPhenotypeRegion(int k, int K) {
	//STEP0: check input values
	if (K < 1) vrb.error("Number of chunks needs to be > 0");
	if (K > phenotype_count) vrb.error("Number of chunks (" + stb.str(K) + ") is greater than the number of phenotypes (" + stb.str(phenotype_count) + ")");
	if (k < 0) vrb.error("Chunk index needs to be > 0");
	if (k >= K) vrb.error("Chunk index needs to be smaller than the total number of chunks [=" + stb.str(K) + "]");

	//STEP1: regroup by group
	vector < pgroup > v_pgroup;
	if (phenotype_grp.size() > 0) {
		map < string, int > grp2idx;
		map < string, int > :: iterator it_grp2idx;
		for (int p = 0 ; p < phenotype_count ; p ++) {
			it_grp2idx = grp2idx.find (phenotype_grp[p]);
			if (it_grp2idx == grp2idx.end()) {
				grp2idx.insert(pair < string, int > (phenotype_grp[p], v_pgroup.size()));
				v_pgroup.push_back(pgroup(phenotype_chr[p], phenotype_start[p], phenotype_end[p]));
			} else v_pgroup[it_grp2idx->second].merge(phenotype_start[p], phenotype_end[p]);
		}
	} else {
		for (int p = 0 ; p < phenotype_count ; p ++) {
			v_pgroup.push_back(pgroup(phenotype_chr[p], phenotype_start[p], phenotype_end[p]));
		}
	}
	sort(v_pgroup.begin(), v_pgroup.end());

	//STEP2: merge overlapping groups
	stack < pgroup > s_pgroup;
	s_pgroup.push(v_pgroup[0]);
	for (int i = 1 ; i < v_pgroup.size(); i++) {
		pgroup ptop = s_pgroup.top();
		if (!ptop.overlap(v_pgroup[i])) s_pgroup.push(v_pgroup[i]);
		else {
			ptop.merge(v_pgroup[i]);
			s_pgroup.pop();
			s_pgroup.push(ptop);
		}
	}
	v_pgroup.clear();
	while (!s_pgroup.empty()) {
		v_pgroup.push_back(s_pgroup.top());
		s_pgroup.pop();
	}
	sort(v_pgroup.begin(), v_pgroup.end());

	//STEP3: build one cluster per chromosome
	vector < vector < int > > cluster_idx;
	map < string , int > chr2idx;
	for (int p = 0 ; p < v_pgroup.size() ; p ++) {
		map < string , int > :: iterator it_chr2idx = chr2idx.find(v_pgroup[p].chr);
		if (it_chr2idx == chr2idx.end()) {
			chr2idx.insert(make_pair(v_pgroup[p].chr, cluster_idx.size()));
			cluster_idx.push_back(vector < int > (1, p));
		} else cluster_idx[it_chr2idx->second].push_back(p);
	}

	//STEP4: split until number of chunks is reached
	bool done = (cluster_idx.size() >= K);
	while (!done) {

		int max_idx = -1, max_val = 1;
		for (int p = 0 ; p < cluster_idx.size() ; p ++) {
			if (cluster_idx[p].size() > max_val) {
				max_val = cluster_idx[p].size();
				max_idx = p;
			}
		}

		if (max_idx >= 0) {
			int max_mid = cluster_idx[max_idx].size() / 2;
			cluster_idx.push_back(vector < int > (cluster_idx[max_idx].begin() + max_mid, cluster_idx[max_idx].end()));
			cluster_idx[max_idx].erase(cluster_idx[max_idx].begin() + max_mid, cluster_idx[max_idx].end());
			if (cluster_idx.size() >= K) done = true;
		} else done = true;
	}

	//STEP5: extract coordinates
	if (k < cluster_idx.size()) {
		regionPhenotype.chr = v_pgroup[cluster_idx[k][0]].chr;
		regionPhenotype.start = 1000000000;
		regionPhenotype.end = 0;
		for (int c = 0 ; c < cluster_idx[k].size() ; c ++) {
			if (v_pgroup[cluster_idx[k][c]].start < regionPhenotype.start) regionPhenotype.start = v_pgroup[cluster_idx[k][c]].start;
			if (v_pgroup[cluster_idx[k][c]].end > regionPhenotype.end) regionPhenotype.end = v_pgroup[cluster_idx[k][c]].end;
		}
	} else vrb.leave("Empty chunk, no data to process!");
}
