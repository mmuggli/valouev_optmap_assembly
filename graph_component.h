#define _GRAPH_COMPONENT_H

#ifndef __IOSTREAM__
#include <iostream>
#endif

#ifndef __CASSERT__
#include <cassert>
#endif

#ifndef __CSTDIO__
#include <cstdio>
#endif


#ifndef __CMATH__
#include <cmath>
#endif

#ifndef __FSTREAM__
#include <fstream>
#endif

#ifndef __STRING__
#include <string>
#endif

#ifndef __SGI_STL_VECTOR
#include <vector>
#endif

#ifndef __SGI_STL_ALGORITHM
#include <algorithm>
#endif

#ifndef _NODE_H
#include "node.h"
#endif

#ifndef _EDGE_H
#include "edge.h"
#endif

#ifndef _PATH_H
#include "path.hpp"
#endif

class graph_component{
 public:
  vector<node> nodes;
  vector<edge> edges;

  hash_key max_hash_key;


  vector<int> node_ind_ht; //given a hash key will find the index of the node in O(1) time

  graph_component(hash_key _max_hash_key);
  void make_hash_table();
  void remove_hash_table();
  void clean();
  /*
  void add_edge(edge& e, vector<fr_size>& left_map_read, 
		vector<fr_size>& right_map_read);
  */

  list<int> consistent_edge_inds(list<int>& good_edge_inds,				 
				 const vector<edge>& _edges,
				 vector<vector<hash_key> >& map_hks,
				 vector<vector<orientation> >& map_orients, 
				 vector<vector<int> >& edge_comp_inds);
  /*
  vector<int> consistent_edge_inds
    (const vector<int>& orig_edge_inds, 
     vector<vector<hash_key> >& map_hks,
     vector<vector<orientation> >& map_orients, 
     vector<vector<int> >& edge_comp_inds);
  */

  void mark_chimeras();

  void print();
  bool cycle(int n); //returns true if there is a cycle
  //involving node n (like in case of circular genomes)
  //stores the path in the cycled_path

  //bool left_to_right_path(string from_node, string to_node);
  //returns true if there is left to right path
  //that goes from_node to_node 
  //the path itself is stored in vector<int> helper1
  
  //bool all_way_path(string from_node, string to_node);

  vector<int> linear_path(int node_ind);
  vector<int> longest_branch(int node_ind, double weight);
  //this is a recursive function
  vector<int> longest_branch(int node_ind);

  vector<double> extract_contig(vector<int>& map_inds);
  vector<double> extract_contig();
  vector<double> extract_contig1(vector<int>& map_inds);

  void assign_edge_weights();
  //must be called after the nodes are indexed
  double distance_between_maps(int cur_edge_ind, int map1_ind, int map2_ind, 
			       orientation map1_orient, 
			       orientation map2_orient);

  bool chimeric_map1(int ind);
  bool chimeric_map2(int ind);
  bool overlap_spans_site(int node_ind, int map_site, orientation _or,
			  int edge_ind);
  void not_left_directed_dfs(int from_node, vector<bool>& discovered_nodes,
			     int search_depth);

  vector<graph_component> components_free_of_chimeras();
  //returns components free of chimeras
  
  void index_nodes(); //does the internal indexing of all 
  //nodes and edges

  void output_overlap(edge& e);

  void output_edges(ofstream& out_str, ostringstream& prefix);
  //for the graph file. requires call from graph

  void output_graph(string fname);

  void output_component(const char* output_file);

  list<int> confirmed_edges(int dfs_depth);
  vector<graph_component> confirmed_components(int dfs_depth);
  //this function breaks a component into strongly
  //connected components with highly accurate edges
  //some nodes and some edges may be removed
  
  void possible_paths_for_consistency_check
    (vector<int>& node_ind, vector<vector<path> >& paths, 
     vector<int>& loc_node_color, vector<int>& loc_pred_list, 
     vector<int>& branch_root_list,
     int mark_color, int depth, int cur_node_ind, int prev_node_ind,
     int master_node_ind, int brach_root_ind);

  void fill_the_component_data();
  //fills out the adj node structure for all nodes

  vector<int> terminal_nodes();
  //returns the list of terminal nodes that
  //can be used as ends for contig extraction 

  vector<int> extract_linear_path();
  void extract_linear_path_dfs
    (int cur_node_ind, int prev_node_ind, double cur_score, 
     vector<int>& local_best_pred_inds,  vector<bool>& node_discovered, 
     vector<int>& node_colors, int master_color);

  
  //int get_node_index(string node_id);
  int get_node_index(hash_key hk);

  int rightmost_al_site(int ind);
  int leftmost_al_site(int ind);
  double finalizing_dist(int ind);
  double starting_dist(int ind);

 private:
  vector<int> cycled_path;
  vector<int> preds;

  vector<short> node_color; //for dfs purposes 
  //0 - white, 1 - grey, 2 - black (for dfs)

  void breadth_first_search(int from_node, vector<bool>& discovered_nodes,
			     int search_depth);

  void dfs_left_to_right_visit(int n);
  bool dfs_left_to_right_visit(int n, int to_reach);
  //n is the index of the current node
  //returns true if dfs reaches to_reach
  //only propogates from left to right to
  //discover things like circular genomes

  void dfs_all_way_visit(int n);
  bool dfs_all_way_visit(int n, int to_reach);
  //propogates in all directions at the same time:
  //left to right, right to left and through containment edges
  //discovers all-way connectivity between nodes within the component

  void dfs_all_way_visit_to_separate_components
    (int map_ind, vector<string>& map_set,
     vector<int>& map_colors, vector<int>& edge_inds, 
     vector<int>& edge_colors, int cur_color);
  //marks edges and maps by the component indicator

  /*
  int edge_ind(int map1_ind, int map2_ind);
  int edge_ind(string map1, string map2);
  */

  //returns the index of the edge corresponding to the 
  //overlap between map1 and map2; -1 if there is no edge

  //  int get_node_index(string node_id);

  vector<int> helper1;
  int int_helper1;
 
  //everyone can overwrite this vector locally.
  //so make sure there is no usage conflict
};
