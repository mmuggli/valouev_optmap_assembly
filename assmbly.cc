#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <list>
#include <vector>
#include <algorithm>
#include <math.h>
#include <ext/hash_map>

using namespace std;
using namespace __gnu_cxx;

#define hash_key long
#define fr_size float
#define site_ind short

#define UNDEF_HASH_KEY -1

vector<list<fr_size> > all_map_reads;
vector<string> all_map_names; //put this in a separate file later
vector<bool> skipped_maps;

string return_map_name(hash_key _key){
  assert(_key >= 0 && _key < all_map_names.size());
  return all_map_names[_key];
}
vector<fr_size> return_map_read(hash_key _key){
  assert(_key >= 0 && _key < all_map_names.size());
  assert(_key >= 0 && _key <all_map_reads.size());
  vector<fr_size> v;
  for(list<fr_size>::iterator it = all_map_reads[_key].begin();
      it != all_map_reads[_key].end(); it++){
    v.push_back(*it);
  }
  return v;
}

struct eqstr{
  bool operator()(const char* s1, const char* s2) const{
    return strcmp(s1,s2)==0;
  }
};

hash_map<const char*, hash_key, hash<const char*>, eqstr> map_hks;

void fill_hk_hash_map(){
  hash_key counter = 0;
  for(vector<string>::iterator it = all_map_names.begin();
      it != all_map_names.end(); it++){
    const char* cur_map = it->c_str();
    map_hks[cur_map] = counter;
    counter++;
  }
}

hash_key get_map_hash_key(const string& str){
  const char* cur_map = str.c_str();
  hash_key cur_hk = map_hks[cur_map];
  
  if(all_map_names[cur_hk] != str){
    cout<<str<<" "<<cur_hk<<endl;
   }
  //assert(all_map_names[cur_hk] == str);
  if(all_map_names[cur_hk] == str) return cur_hk;
  else return UNDEF_HASH_KEY;
}

void mark_skipped_map(hash_key _hk){
  assert(_hk >= 0 && _hk < skipped_maps.size());
  skipped_maps[_hk] = true;
}

#include "graph.cpp"

bool t_score_compare(const edge& e1, const edge& e2){
  if(e1.t_score>e2.t_score) return true;
  else return false;
}

bool is_in_list(const string& cur_name, const vector<string>& names){
  for(int i=0; i<names.size(); i++){
    if(cur_name == names[i]) return true;
  }
  return false;
}

vector<double> get_map(string& cur_map, vector<string>& map_names, 
		       vector<vector<double> >& map_reads, int& id){
  for(int i=0; i<map_names.size(); i++){
    if(map_names[i] == cur_map){
      id = i;
      return map_reads[i];
    }
  }
  vector<double> empty_vec;
  id = UNDEF_IND;
  return empty_vec;
}

/*
vector<double> get_map(hash_key map_hash_key, vector<vector<double> >& map_reads){
  assert(map_hash_key >=0 && map_hash_key < map_reads.size());
  vector<double> v;
  v = map_reads[map_hash_key];
  return v;
}

list<double> get_map(hash_key map_hash_key, vector<list<double> >& map_reads){
  assert(map_hash_key >=0 && map_hash_key < map_reads.size());
  list<double> l;
  l = map_reads[map_hash_key];
  return l;
}
*/

vector<double> get_map(hash_key map_hash_key, vector<list<double> >& map_reads){
  assert(map_hash_key >=0 && map_hash_key < map_reads.size());
  vector<double> v;
  for(list<double>::iterator it = map_reads[map_hash_key].begin(); 
      it != map_reads[map_hash_key].end(); it++){
    v.push_back(*it);
  }
  return v;
}

bool ridiculous(list<fr_size>& l){
  fr_size min_size = 0;
  fr_size max_size = 1000;
  for(list<fr_size>::iterator it=l.begin(); it!=l.end(); it++){
    if(*it < min_size || *it > max_size) return true;
  }
  return false;
}


vector<double> read_map_from_file(ifstream& ifs, string& map_id){
  vector<double> ret_vec;

  ifs.seekg(0, ios::beg);
  if(ifs.eof())
    ifs.clear();
  ifs.seekg(0);

  bool map_found = false;
  bool map_eof = false;
  while(!map_found && !map_eof){
    string line1, line2, line3;
    getline(ifs, line1);
    map_eof = getline(ifs,line2).eof();
    getline(ifs, line3);

    if(!map_eof){
      if(line1 == map_id){
	map_found = true;

	istringstream cur_str;
	cur_str.str(line2);

	string enz, enz_acr;
	cur_str>>enz>>enz_acr;

	int last_pos = cur_str.tellg();
	while(cur_str.good()){
	  double cur_fr;
	  cur_str>>cur_fr;
	  
	  int cur_pos = cur_str.tellg();

	  if(cur_pos > last_pos){
	    assert(cur_fr > 0);
	    ret_vec.push_back(cur_fr);
	  }
	}	
      }
    }
  }
  ifs.clear();
  ifs.seekg(0);
  return ret_vec;
}

float t_score_drop(vector<site_ind>& al_sites1, vector<site_ind>& al_sites2){
  double nu = 0.9;
  double lambda= 1.2;

  double max_t_drop = 0;//-3.0;

  double cur_max_t = 0;
  double cur_t = 0;

  for(int i=1; i<al_sites1.size(); i++){
    short cur_gap1 = al_sites1[i]-al_sites1[i-1];
    short cur_gap2 = al_sites2[i]-al_sites2[i-1];

    cur_t += nu-lambda*(cur_gap1+cur_gap2-2);
    if(cur_t > cur_max_t) cur_max_t = cur_t;

    if(cur_t - cur_max_t < max_t_drop) max_t_drop = cur_t-cur_max_t;
  }
  return max_t_drop;
}

double mm_ratio(vector<site_ind>& al_sites1, vector<site_ind>& al_sites2){
  int single_match_count = 0;
  int long_match_count = 0;
  for(int i=1; i<al_sites1.size(); i++){
    int gap1 = al_sites1[i]-al_sites1[i-1];
    int gap2 = al_sites2[i]-al_sites2[i-1];

    int max_gap = int_max(gap1,gap2);
    if(max_gap==1) single_match_count++;
    else long_match_count += max_gap;
  }
  if(long_match_count == 0) return 1000;
  return ((double) single_match_count)/((double)long_match_count);
}

bool sat_ovlp(edge& e){
  bool use_containment = true;//false;

  double t_score_drop_thresh = -3.0;
  double t_score_thresh_l = 8.0;//8.0;//8.0;
  double t_score_thresh_u = 10.0;
  double mm_ratio_thresh_l = 1.88;//1.3;//1.66;
  double mm_ratio_thresh_u = 1.5;
  
  float cur_t_score_drop = t_score_drop(e.left_al_sites,e.right_al_sites);
  double cur_mm_ratio = mm_ratio(e.left_al_sites,e.right_al_sites);
  double cur_t_score = e.t_score;

  if(cur_t_score_drop < t_score_drop_thresh) return false;
  if(cur_t_score < t_score_thresh_l) return false;
  if(cur_t_score >= t_score_thresh_l && 
     cur_t_score < t_score_thresh_u){
    if(cur_mm_ratio >= mm_ratio_thresh_u) return true;
    else return false;
  }
  if(cur_t_score >= t_score_thresh_u){
    if(cur_mm_ratio >= mm_ratio_thresh_l) return true;
    else return false;
  }
  return false;

  /*
  if(cur_t_score >= t_score_thresh &&
     cur_t_score_drop >= t_score_drop_thresh &&
     cur_mm_ratio >= mm_ratio_thresh){
    if(use_containment) return true;
    else{
      if(e.containment) return false;
      else return true;
    }
  }
  else return false;
  */
}

int main(){
  const char* map_filename = 
    //"./../datasets/human/mole_maps/mole.maps_350Kb_15fr.maps";
    //"./../datasets/human/gm/lymph_12fr.maps";
    //"./../datasets/rice/assm3/tmp/tmpmaps";
    // "./../datasets/rice/assm2/tmp/rice.contigs.133.maps";
    //"./../datasets/rice/rice_450Kb_15fr.maps";
    "./../datasets/rice/rice_all.maps";
    //"./../datasets/ecoli/StuI/Ecoli_StuI.maps";
    //"./../datasets/pseudonana/tpnana_all.maps";
    // "./../datasets/pseudonana/test.maps";
    //"./../datasets/pseudonana/Tpnana_20fr.maps";
    //"./../datasets/rice/rice_700Kb_20fr.maps";
    //"./../datasets/rice/chr02.opt.maps";
    //"./../datasets/ecoli/StuI/Ecoli_StuI_top700_maps.stm";
    //"./../datasets/rice/chr12.opt.maps";
    //"./../datasets/pestis/Ypestis_non_dupl.maps";
  //"./../datasets/rice/selected.maps";
    
  ifstream map_str;
  map_str.open(map_filename);
  if(!map_str.good()){
    cerr<<"wrong map filename: ";
    cerr<<map_filename<<endl;
  }

  const char* ovlp_filename = 
    //"./../datasets/human/mole_maps/mole_350Kb_15fr_new.ovlps";
    //"./../datasets/human/gm/lymph_12fr.ovlps";
    "./../datasets/rice/rice_all_new.ovlps";
    //"./../datasets/rice/rice_450Kb_15fr_new.ovlps";
    //"./../datasets/ecoli/StuI/Ecoli_StuI.ovlps";
    //"./../datasets/pseudonana/tpnana_all_new.ovlps";
    //"./../datasets/pseudonana/tpnana_20fr.ovlps";
    //"./../datasets/rice/rice.700Kb_20fr.ovlps";
    //"./../filtration/chr01.ovlps";
    //"./../datasets/rice/chr02.ovlps";
    //"./../datasets/ecoli/StuI/Ecoli_StuI_top700_maps.ovlps";
    //"./../datasets/rice/chr12.ovlps";
    //"./../datasets/pestis/all_ovlps";
    //"./../datasets/pestis/ovlps";
    //"./../ovlp/out";
    //"./../datasets/rice/ovlps";
  ifstream ovlp_str;
  ovlp_str.open(ovlp_filename);

  if(!ovlp_str.good()){
    cerr<<"wrong overlap filename: ";
    cerr<<ovlp_filename;
    cerr<<endl;
    assert(false);
  }

  //const char* contig_file_name = 
    
  string comp_prefix = 
    //"./../datasets/human/mole_maps/assm4/mole.contigs";
    //"./../datasets/human/gm/lymph.contigs";
    "./../datasets/rice/assm8/rice.contigs";
    //"./../datasets/rice/assm3/tmp/rice.contigs";
    //"./../datasets/ecoli/StuI/assm1/contigs";
    //"./../datasets/pseudonana/assm4/tpnana.contigs";
    //"./../datasets/pseudonana/assm2/tpnana.contigs";
    //"./../datasets/rice/assm2/tmp/rice.contigs";
    //"./../datasets/rice/chr12.contigs";
    //"./../datasets/ecoli/StuI/Ecoli_StuI_top700_maps.contigs";
    //"./../datasets/pestis/pestis.contigs";
    //"./../datasets/selected.contigs";
    //"./../datasets/ecoli/StuI/ecoli.contigs";
  //remove(contig_file_name);
  //ofstream contig_str;
  //contig_str.open(contig_file_name);

  const char* skipped_maps_file_name = 
    "./../datasets/rice/skipped_maps";
  //"./../datasets/pseudonana/chimeras";
  ifstream skipped_maps_str;
  skipped_maps_str.open(skipped_maps_file_name);

  /*
  if(!contig_str.good()){
    cerr<<"wrong contig file name: ";
    cerr<<contig_file_name;
  }
  */

  vector<string> skipped_names;
  while(skipped_maps_str.good()){
    string line1;
    getline(skipped_maps_str, line1);
    if(skipped_maps_str.good())
      skipped_names.push_back(line1);
  }
  skipped_maps_str.close();

  int end_delta = 1; //maximum number of allowed danglers

  vector<edge> edges;

  int ovlps = 0;
  int contained = 0;

  hash_key max_hash_key = UNDEF_HASH_KEY;

  all_map_reads.clear();

  map_str.seekg(0, ios::beg);
  if(map_str.eof())
    map_str.clear();
  map_str.seekg(0);
  

  bool map_eof = false;
  //long map_counter = 0;
  while(!map_eof){
    string line1, line2, line3;
    getline(map_str, line1);
    map_eof = getline(map_str,line2).eof();
    getline(map_str, line3);

    if(!map_eof){
      list<float> cur_read;
      istringstream cur_str;
      cur_str.str(line2);      

      string enz, enz_acr;
      cur_str>>enz>>enz_acr;

      int last_pos = cur_str.tellg();
      while(cur_str.good()){
	fr_size cur_fr;
	cur_str>>cur_fr;
	  
	int cur_pos = cur_str.tellg();
	
	if(cur_pos > last_pos){
	  assert(cur_fr > 0);
	  cur_read.push_back(cur_fr);
	}
      }
      all_map_reads.push_back(cur_read);
      all_map_names.push_back(line1);

      //map_counter ++;
    }
  }
  
  cerr<<"maps loaded"<<endl;
  max_hash_key = all_map_names.size();

  //fill structures begin
  //------------------------------------
  cout<<"filling the structures...";
  fill_hk_hash_map();
  skipped_maps.clear();
  for(int i=0; i<all_map_names.size(); i++){
    skipped_maps.push_back(false);
  }
  for(int i=0; i<skipped_names.size(); i++){
    cout<<i<<endl;
    hash_key cur_hk = get_map_hash_key(skipped_names[i]);
    if(cur_hk != UNDEF_HASH_KEY)
      skipped_maps[cur_hk] = true;
 //skip this map in future
  }
  cout<<" done!"<<endl;
  //------------------------------------
  //fill structures end

  while(ovlp_str.good() /*&& ovlps < 60000*/){
    bool cur_status_good = true;

    string line1, line2, line3;
    
    getline(ovlp_str, line1);
    cur_status_good = getline(ovlp_str, line2).good();
    getline(ovlp_str, line3);

    if(cur_status_good){
      ovlps++;
      //bool edge_good = true;

      istringstream cur_string_stream1;
      cur_string_stream1.str(line1);

      string map1, map2;
      hash_key map1_hash_key, map2_hash_key;
      int map1_size, map2_size;
      orientation or1, or2;
      double s_score, t_score;      

      cur_string_stream1>>map1_hash_key>>map2_hash_key;
      cur_string_stream1>>map1>>map2>>map1_size>>map2_size;
      cur_string_stream1>>or1>>or2>>s_score>>t_score;
      istringstream cur_string_stream2;
      cur_string_stream2.str(line2);

      if(!valid_orient(or1)) {cerr<<ovlps<<endl; assert(false);}
      if(!valid_orient(or2)) {cerr<<ovlps<<endl; assert(false);}

      //if(map1_hash_key > max_hash_key) max_hash_key = map1_hash_key;
      //if(map2_hash_key > max_hash_key) max_hash_key = map2_hash_key;

      string map1_id = return_map_name(map1_hash_key);
      string map2_id = return_map_name(map2_hash_key);

      //debug begin
      if(map1_id != map1 || map2_id != map2){
	cerr<<map1_hash_key<<" "<<map1<<" ("<<map1_id<<")"<<endl;
	cerr<<map2_hash_key<<" "<<map2<<" ("<<map2_id<<")"<<endl;
      }
      //else edge_good = false;
      //debug end

      assert(map1_id == map1);
      assert(map2_id == map2);
      
      vector<site_ind> al_sites1, al_sites2;
      
      int last_pos = cur_string_stream2.tellg();
      while(cur_string_stream2.good()){
	int map1_ind, map2_ind;

	cur_string_stream2>>map1_ind>>map2_ind;

	int cur_pos = cur_string_stream2.tellg();
	if(cur_pos>last_pos){
	  al_sites1.push_back(map1_ind);
	  al_sites2.push_back(map2_ind);
	}
	last_pos = cur_pos;
      }

      //cerr<<"here"<<endl;
      //cerr<<or1<<" "<<or2<<endl;
      //cerr<<map1_hash_key<<" "<<map2_hash_key<<endl;
      int al_sites_num = al_sites1.size();
      if(al_sites2[0] <= end_delta){
	if(al_sites2[al_sites_num-1]>= map2_size-end_delta){
	  //right map is contained
	  
	  contained++;
	  edge cur_edge(map1, map2, 
			map1_hash_key, map2_hash_key,
			or1, or2, map1_size, map2_size,
			s_score, t_score, al_sites1, al_sites2, true);
	  if(sat_ovlp(cur_edge)) edges.push_back(cur_edge);
	}
	else{
	  if(al_sites1[0] <= end_delta &&
	     al_sites1[al_sites_num-1] >= map1_size-end_delta){
	    //left map is contained
	    contained++;
	    edge cur_edge(map2, map1, 
			  map2_hash_key, map1_hash_key,
			  or2, or1, map2_size, map1_size,
			  s_score, t_score, al_sites2, al_sites1, true);
	    if(sat_ovlp(cur_edge)) edges.push_back(cur_edge);
	  }
	  else{
	    //no containment
	    edge cur_edge(map1, map2,
			  map1_hash_key, map2_hash_key,
			  or1, or2, map1_size, map2_size,
			  s_score, t_score, al_sites1, al_sites2, false);
	    if(sat_ovlp(cur_edge)) edges.push_back(cur_edge);
	  }
	}
      }
      else{
	//right map is not contained
	if((al_sites1[0] <= end_delta)==false){
	  cout<<al_sites1[0]<<" "<<al_sites1[al_sites1.size()-1];
	  cout<<" total: "<<map1_size<<endl;
	}
	assert(al_sites1[0] <= end_delta);
	
	if(al_sites1[0] <= end_delta){
	 
	  if(al_sites1[al_sites_num-1] >= map1_size-end_delta){
	    //left map is contained
	    edge cur_edge(map2, map1, 
			  map2_hash_key, map1_hash_key,
			  or2, or1, map2_size, map1_size,
			  s_score, t_score, al_sites2, al_sites1, true);
	    if(sat_ovlp(cur_edge)){
	      contained++;
	      edges.push_back(cur_edge);
	    }
	  }
	  else{
	    //no containment
	    edge cur_edge(map2, map1,
			  map2_hash_key, map1_hash_key,
			  or2, or1, map2_size, map1_size,
			  s_score, t_score, al_sites2, al_sites1, false);
	    
	    if(sat_ovlp(cur_edge)) edges.push_back(cur_edge);
	  }
	}
      }
    }
  }
  
  //max_hash_key = all_map_names.size();
  //cerr<<"max_hash_key: "<<max_hash_key<<endl;

  cerr<<"constructing first_graph"<<endl;
  sort(edges.begin(), edges.end(), t_score_compare);
  cerr<<"edges sorted"<<endl;

  cerr<<"constructing hash table for maps:"<<endl;
  //vector<bool> maps_used;

  //maps_used.clear();
  //for(hash_key i=0; i<max_hash_key; i++) maps_used.push_back(false);

  int ovlps_used = 0;

  //vector<string> map_names_to_use;
  //vector< vector<double> > map_reads_to_use;
  //list<hash_key> map_hks_to_use;

  //vector<edge> good_edges;
  list<int> good_edge_inds;

  cerr<<"finding maps corresponding to overlaps"<<endl;
  

  for(int i=0; i<edges.size(); i++){
    //cout<<i<<endl;
    if(i%1000 == 0){
      cerr<<"processed "<<i<<" edges of ";
      cerr<<edges.size()<<" used: "<<ovlps_used<<endl;
    }
    /*
    if(i==testi){
      cerr<<"here3"<<endl;
      int _size = edges[i].left_al_sites.size();
      cerr<<edges[i].left_map<<" "<<edges[i].right_map;
      cerr<<" "<<edges[i].left_map_hash_key<<" "<<edges[i].right_map_hash_key;
      cerr<<" "<<edges[i].left_map_size<<" "<<edges[i].right_map_size<<" "<<edges[i].left_al_sites[0];
      cerr<<" "<<edges[i].left_al_sites[_size-1];
      cerr<<" "<<edges[i].right_al_sites[0]<<" "<<edges[i].right_al_sites[_size-1]<<endl;
    }
    */

    //int bad_edge_counter = 0;
    if(sat_ovlp(edges[i])){
      bool good_edge = true;
      //if(i==testi) cerr<<"here1"<<endl;
      ovlps_used++;

      hash_key map1_hash_key = edges[i].left_map_hash_key;
      hash_key map2_hash_key = edges[i].right_map_hash_key;

      assert(map1_hash_key>=0 && map1_hash_key <= max_hash_key);
      assert(map2_hash_key>=0 && map2_hash_key <= max_hash_key);
      assert(map1_hash_key < all_map_reads.size());
      assert(map2_hash_key < all_map_reads.size());
      
      string map1 = edges[i].left_map;
      string map2 = edges[i].right_map;
      
      int map1_size = edges[i].left_map_size;
      int map2_size = edges[i].right_map_size;      
      /*
      vector<double> map_read1 = 
	get_map(map1_hash_key, orig_map_reads);
      vector<double> map_read2 = 
	get_map(map2_hash_key, orig_map_reads);
      */
      //if(i==testi) cout<<"here2"<<endl;
      
      if(!(map1_size == all_map_reads[map1_hash_key].size())) cerr<<i<<endl;

      if(map1_size != all_map_reads[map1_hash_key].size() &&
	 map2_size != all_map_reads[map2_hash_key].size()) good_edge = false;
      //assert(map1_size == all_map_reads[map1_hash_key].size());
      //assert(map2_size == all_map_reads[map2_hash_key].size());      

      /*
      if(i==testi){
	int _size = edges[i].left_al_sites.size();
	cerr<<map1<<" "<<map2<<" "<<map1_hash_key<<" "<<map2_hash_key;
	cerr<<" "<<map1_size<<" "<<map2_size<<" "<<edges[i].left_al_sites[0];
	cerr<<" "<<edges[i].left_al_sites[_size-1];
	cerr<<" "<<edges[i].right_al_sites[0]<<" "<<edges[i].right_al_sites[_size-1]<<endl;
      }
      */

      //if(!map_read1.empty() && !map_read2.empty()
      // && !ridiculous(map_read1) && !ridiculous(map_read2)
      // && !is_in_list(map1, skipped_names) && !is_in_list(map2, skipped_names)){
      if(!all_map_reads[map1_hash_key].empty() 
	 && !all_map_reads[map2_hash_key].empty()
	 && !ridiculous(all_map_reads[map1_hash_key]) 
	 && !ridiculous(all_map_reads[map2_hash_key])
	 //&& !is_in_list(all_map_names[map1_hash_key], skipped_names)
	 //&& !is_in_list(all_map_names[map2_hash_key], skipped_names)
	 && !skipped_maps[map1_hash_key]
	 && !skipped_maps[map2_hash_key]
	 && good_edge){

        //good_edges.push_back(edges[i]);
	good_edge_inds.push_back(i);

	/*
	bool map1_used_before = maps_used[map1_hash_key];
	bool map2_used_before = maps_used[map2_hash_key];

	if(!map1_used_before){
	  //map_names_to_use.push_back(map1);
	  //map_reads_to_use.push_back(map_read1);
	  maps_used[map1_hash_key] = true;
	  map_hks_to_use.push_back(map1_hash_key);
	}
	if(!map2_used_before){
	  //map_names_to_use.push_back(map2);
	  //map_reads_to_use.push_back(map_read2);
	  maps_used[map2_hash_key] = true;
	  map_hks_to_use.push_back(map2_hash_key);
	}
	*/
      }
      else{
	//map not found, do nothing
      }
    }
  }
  
  //cerr<<"using "<<maps_used.size()<<" out of "<<max_hash_key+1<<" maps"<<endl;
  cerr<<"using "<<good_edge_inds.size()<<" overlaps out of "<<edges.size()<<endl;
  cerr<<"compiled a list of good edges"<<endl;
  //sort(good_edges.begin(), good_edges.end(), t_score_compare);
  //cerr<<"edges sorted."<<endl;
  //edges.clear();

  /*
  //vector<node> nodes;
  for(list<hash_key>::iterator it = map_hks_to_use.begin();it != map_hks_to_use.end(); it++){
    hash_key cur_hk = *it;
    vector<double> map_read = get_map(cur_hk, orig_map_reads);
    string map_name = all_map_names[cur_hk];
    //node new_node(cur_hk,map_name,forward,map_read);
    //nodes.push_back(new_node);
  }
  */
  //cerr<<"created the first list of nodes"<<endl;
  //maps_used.clear();
  //map_hks_to_use.clear();

  graph g(max_hash_key);
  g.construct_graph(good_edge_inds, edges);
  /*
  g.construct_graph(good_edges, nodes);
  good_edges.clear();
  nodes.clear();
  */

  //all_map_reads.clear();
  
  cerr<<"first graph constructed"<<endl;

  //g.init_data();
  cerr<<"structures initialized"<<endl;
  //g.print();
  //g.output_graph("./orig_graph.dot");
  
  cout<<"total ovlps: "<<ovlps<<endl;
  cout<<"contained: "<<contained<<endl;
  cout<<"using: "<<ovlps_used<<endl;

  
  //int comp_ind = 0;
  //g.components[comp_ind].strongly_connected_components(4);
  //g.components[comp_ind].confirmed_edges();

  graph conf_g(max_hash_key);
  
  cerr<<"components: "<<g.components.size()<<endl;
  cerr<<"confirmed components routine"<<endl;

  for(int i=0; i<g.components.size(); i++){
    cerr<<i<<"-th component of "<<g.components.size()<<endl;;
    cout<<"processing component: "<<i<<endl;
    vector<graph_component> new_comps;
    new_comps = g.components[i].confirmed_components(2);
    for(int j=0; j<new_comps.size(); j++){
      cout<<"adding component"<<endl;
      conf_g.components.push_back(new_comps[j]);
    }
    cerr<<"finished"<<endl;
  }
  
  g.clear();

  //conf_g.print();
  //conf_g.components[0].output_graph("./conf_graph.dot");
  
  int chimera_cycle = 4;
  for(int c=0; c<chimera_cycle; c++){
    
    //eliminating chimeras
    cerr<<"elimination of chimeric nodes"<<endl;
    //cout<<"elimination of chimeric nodes"<<endl;
    //conf_g.mark_chimeras();
    for(int i=0; i<conf_g.components.size(); i++){
      cerr<<"finding chimeras in comp: "<<i<<" of "<<conf_g.components.size()-1<<endl;
      conf_g.components[i].mark_chimeras();
    }

    cerr<<"chimeras found"<<endl;
    list<int> good_edges2;
    for(list<int>::iterator it = good_edge_inds.begin(); 
	it!= good_edge_inds.end(); it++){
      int cur_edge_ind = *it;
      hash_key cur_left_map_hk = edges[cur_edge_ind].left_map_hash_key;
      hash_key cur_right_map_hk = edges[cur_edge_ind].right_map_hash_key;
      if(!skipped_maps[cur_left_map_hk] && 
	 !skipped_maps[cur_right_map_hk]) good_edges2.push_back(cur_edge_ind);
    }
    
    cerr<<"New edge set contains "<<good_edges2.size()<<" overlaps"<<endl;
    cerr<<"making a new graph"<<endl;
    conf_g.clear();
    conf_g.construct_graph(good_edges2, edges);
    
    cerr<<"confirming edges"<<endl;
    
    //graph new_g(max_hash_key);
    vector<graph_component> new_comps;
    for(int i=0; i<conf_g.components.size(); i++){
      cerr<<i<<"-th component of "<<conf_g.components.size()<<endl;;
      cout<<"processing component: "<<i<<endl;
      vector<graph_component> temp_comps 
	= conf_g.components[i].confirmed_components(2);
      for(int j=0; j<temp_comps.size(); j++){
	cout<<"adding component"<<endl;
	new_comps.push_back(temp_comps[j]);
      }
      cerr<<"finished"<<endl;
    }
    
    conf_g.clear();
    conf_g.components = new_comps;
  }

  graph new_g(max_hash_key);
  new_g.components = conf_g.components;
  conf_g.clear();

  
  /*
  int cycle = 2;
  for(int i=0; i<2; i++){
    
    vector<graph_component> new_comps;
    for(int i=0; i<conf_g.components.size(); i++){    
      vector<graph_component> cur_comps = 
	conf_g.components[i].components_free_of_chimeras();
      for(int j=0; j<cur_comps.size(); j++){
	new_comps.push_back(cur_comps[j]);
      }
    }
    
    conf_g.components.clear();
    conf_g.components = new_comps;
    new_comps.clear();
    
    vector<graph_component> conf_comps;
    for(int i=0; i<conf_g.components.size(); i++){
      vector<graph_component> cur_comps = 
	conf_g.components[i].confirmed_components(2);
      for(int j=0; j<cur_comps.size(); j++){
	cout<<"adding component"<<endl;
	conf_comps.push_back(cur_comps[j]);
      }
    }
    conf_g.components.clear();
    conf_g.components = conf_comps;
    conf_comps.clear();
  }
  
  graph new_g(max_hash_key);
  new_g.components = conf_g.components;
  conf_g.clear();
  */
  // cerr<<"constructed chimera-free graph"<<endl;

  int min_map_num = 10;
  int contig_count = 0;
  for(int i=0; i<new_g.components.size(); i++){
    if(new_g.components[i].nodes.size()>=min_map_num){
      std::stringstream contig_id;
      contig_id<<contig_count;
      
      string loc_graph_fname = comp_prefix+"."+contig_id.str()+".dot";
      new_g.components[i].output_graph(loc_graph_fname.c_str());
      
      string contig_maps_fname = comp_prefix+"."+contig_id.str()+".maps";
      remove(contig_maps_fname.c_str());
      ofstream cur_maps_str;
      cur_maps_str.open(contig_maps_fname.c_str());
      assert(cur_maps_str.good());
      
      //storing optical maps
      for(int j=0; j<new_g.components[i].nodes.size(); j++){
	if(j!=0) cur_maps_str<<endl;
	cur_maps_str<<new_g.components[i].nodes[j].read_name<<endl;
	cur_maps_str<<char(9)<<"NheI"<<char(9)<<"N";
	for(int k=0; k<new_g.components[i].nodes[j].map_read.size(); k++){
	  cur_maps_str<<char(9)<<new_g.components[i].nodes[j].map_read[k];
	}
	cur_maps_str<<endl;
      }
      
      cur_maps_str.close();
      
      //storing consensus
      string contig_fname = comp_prefix+"."+contig_id.str()+".cons";
      remove(contig_fname.c_str());
      ofstream cons_str;
      cons_str.open(contig_fname.c_str());
      assert(cons_str.good());
      
      cout<<"best path for the component: "<<i<<endl;
      
      vector<int> best_path;
      
      best_path = new_g.components[i].extract_linear_path();
      vector<double> cur_contig = new_g.components[i].extract_contig1(best_path);
      
      cons_str<<"contig_"<<contig_count<<endl;
      cons_str<<char(9)<<"NheI"<<char(9)<<"N";
      for(int k=0; k<cur_contig.size(); k++){
	cons_str<<char(9)<<cur_contig[k];
      }
      cons_str<<endl<<endl;
      
      
      cons_str.close();
      contig_count++;
    }
  }
  
  /*
  for(int i=0; i<new_g.components.size(); i++){
    cout<<"best path for the component: "<<i<<endl;

    vector<int> best_path;

    best_path = new_g.components[i].extract_linear_path();
    vector<double> cur_contig = new_g.components[i].extract_contig1(best_path);

    contig_str<<"contig_"<<i<<endl;
    contig_str<<char(9)<<"XhoI"<<char(9)<<"X";
    for(int k=0; k<cur_contig.size(); k++){
      contig_str<<char(9)<<cur_contig[k];
    }
    contig_str<<endl<<endl;
  }
  */


  //new_g.output_graph("./new_graph.dot");
  

  /*
  cout<<"discovering the cycles:"<<endl;
  for(int i=0; i<g.components[comp_ind].nodes.size(); i++){
    //if(g.components[comp_ind].cycle(i) == true) assert(false);
    cout<<"cycle extraction: processing node "<<i<<endl;
    bool cycle_exists = g.components[comp_ind].cycle(i);
    if(cycle_exists){
      cout<<"found a cycle"<<endl;
      vector<double> cur_contig = 
	g.components[comp_ind].extract_contig();
      contig_str<<"contig_"<<i<<endl;
      contig_str<<char(9)<<"XhoI"<<char(9)<<"X";
      for(int k=0; k<cur_contig.size(); k++){
	contig_str<<char(9)<<cur_contig[k];
      }
      contig_str<<endl<<endl;
    }
    else{
      cout<<"no cycle here"<<endl<<endl;
    }
  }
  
  
  for(int i=0; i<g.components[comp_ind].nodes.size(); i++){
    vector<int> path = 
      //g.components[comp_ind].linear_path(i);
      g.components[comp_ind].longest_branch(i);
    
    
    if(path.size() >= 2){
      vector<double> cur_contig = 
	g.components[comp_ind].extract_contig(path);
          
      contig_str<<"contig_";
      //contig_str<<i<<endl;
      contig_str<<path[0]<<"->"<<path[path.size()-1]<<endl;
      contig_str<<char(9)<<"XhoI"<<char(9)<<"X";
      for(int k=0; k<cur_contig.size(); k++){
	contig_str<<char(9)<<cur_contig[k];
      }
      contig_str<<endl<<endl;
    }
    
  }
  */
  //contig_str.close();
}
