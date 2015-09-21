#include "graph_component.h"

int int_max(int x, int y){
  if(x>y) return x;
  else return y;
}
int int_min(int x, int y){
  if(x<y) return x;
  else return y;
}
orientation rev(orientation _or){
  assert(_or == forward || _or == reverse);
  if(_or == reverse) return forward;
  else return reverse;
}

string suffix(string str, int suff_size){
  assert(suff_size <= str.size());
  string res;

  for(int i = str.size()-suff_size; i<str.size(); i++){
    assert(i>=0);
    if(str[i]=='_') res+="_";
    else res += str[i];
  }
  //cerr<<res.c_str()<<endl;
  return res;
}

void graph_component::make_hash_table(){
  node_ind_ht.clear();
  for(int i=0; i<=max_hash_key; i++){
    node_ind_ht.push_back(UNDEF_IND);
  }
  for(int i=0; i<nodes.size(); i++){
    assert(nodes[i].hk != UNDEF_HASH_KEY);
    assert(nodes[i].hk >= 0 && nodes[i].hk <= max_hash_key);
    node_ind_ht[nodes[i].hk] = i;
  }
}
void graph_component::remove_hash_table(){
  node_ind_ht.clear();
}
int suffix_length = 10;

int graph_component::get_node_index(hash_key hk){
  assert(!node_ind_ht.empty());
  assert(node_ind_ht.size() == max_hash_key+1);
  assert(node_ind_ht[hk] != UNDEF_IND);
  return node_ind_ht[hk];
}

void graph_component::mark_chimeras(){
  int counter = 0;
  for(int i=0; i<nodes.size(); i++){
    if(chimeric_map2(i)){
      mark_skipped_map(nodes[i].hk);
      counter++;
    }
  }
  cerr<<"marked "<<counter<<" maps as possible chimeras "<<endl;
}

vector<graph_component> graph_component::components_free_of_chimeras(){
  vector<string> marked_chimeras;
  vector<bool> chim_nodes;
  int chim_counter = 0;
  for(int i=0; i<nodes.size(); i++){
    if(chimeric_map2(i)){
      cerr<<"chimera: "<<suffix(nodes[i].read_name,suffix_length)<<endl;
      marked_chimeras.push_back(nodes[i].read_name);
      chim_nodes.push_back(true);
      chim_counter ++;
    }
    else chim_nodes.push_back(false);
  }

  if(chim_counter>0) cerr<<char(9)<<"eliminated "<<chim_counter<<" chimeras"<<endl;

  list<int> good_edge_inds;
  for(int i=0; i<edges.size(); i++){
    bool bad_edge = false;
    for(int j=0; j<marked_chimeras.size(); j++){
      if(edges[i].left_map == marked_chimeras[j] ||
	 edges[i].right_map == marked_chimeras[j])
	bad_edge = true;
    }
    if(!bad_edge) good_edge_inds.push_back(i);
  }

  //vector<vector< string > > map_names;
  vector<vector< hash_key > > map_hash_keys;
  vector<vector< orientation > > map_orients;
  vector<vector< int > > new_edges;
 
  list<int> consistent_edge_set = 
    consistent_edge_inds(good_edge_inds, edges, map_hash_keys, map_orients,
			 new_edges);
    //consistent_edge_inds(good_edge_inds, map_hash_keys, map_orients, new_edges);

  vector<graph_component> components;

  for(int i=0; i<map_hash_keys.size(); i++){
    cout<<"processing "<<i<<" component"<<endl;
    assert(map_hash_keys[i].size() == map_orients[i].size());
    
    graph_component new_comp(max_hash_key);
    assert(map_hash_keys[i].size() == map_orients[i].size());
    for(int j=0; j<map_hash_keys[i].size(); j++){
      hash_key cur_map_hash_key = map_hash_keys[i][j];      
      string cur_map_name = return_map_name(cur_map_hash_key);
      //int cur_node_index = get_node_index(cur_map_name);
      //cout<<"adding a map ("<<cur_node_index<<") "<<cur_map_name<<endl;
      //assert(cur_node_index != -1);

      vector<fr_size> cur_map_read = return_map_read(cur_map_hash_key);
      //node new_node(nodes[cur_node_index].hk, nodes[cur_node_index].read_name, 
      //		    map_orients[i][j], nodes[cur_node_index].map_read);
      //new_node.id = nodes[cur_node_index].id;
      node new_node(cur_map_hash_key, cur_map_name, map_orients[i][j], cur_map_read);

      new_comp.nodes.push_back(new_node);
    }
    
    /*
    for(list<int>::iterator it = consistent_edge_set.begin();
	it != consistent_edge_set.end(); it++){
      int cur_edge_ind = *it;
      assert(cur_edge_ind >= 0 && cur_edge_ind<edges.size());
      new_comp.edges.push_back(edges[cur_edge_ind]);    
    }
    */
    
    for(int j=0; j<new_edges[i].size(); j++){
      int cur_edge_ind = new_edges[i][j];
      cout<<"adding an edge:"<<cur_edge_ind<<endl;
      assert(cur_edge_ind >= 0 && cur_edge_ind<edges.size());
      new_comp.edges.push_back(edges[cur_edge_ind]);    
    }
    

    cout<<"filling out the component structure."<<endl;
    new_comp.fill_the_component_data();
    //new_comp.assign_edge_weights();
    cout<<"structure completed"<<endl;
    if(new_comp.nodes.size() > 0) components.push_back(new_comp);
  }  
  
  return components;
}

bool graph_component::
overlap_spans_site(int node_ind, int map_site, orientation _or, int edge_ind){
  assert(node_ind>=0 && node_ind<nodes.size());
  assert(edge_ind>=0 && edge_ind<edges.size());

  int al_start, al_end;
  int sites = edges[edge_ind].left_al_sites.size();

  if(nodes[node_ind].read_name == edges[edge_ind].left_map){
    if(edges[edge_ind].left_orient == _or){
      al_start = edges[edge_ind].left_al_sites[0];
      al_end = edges[edge_ind].left_al_sites[sites-1];
    }
    else{
      al_start = nodes[node_ind].map_read.size()-
	edges[edge_ind].left_al_sites[sites-1];
      al_end = nodes[node_ind].map_read.size()-
	edges[edge_ind].left_al_sites[0];
    }
  }
  else{
    assert(nodes[node_ind].read_name == edges[edge_ind].right_map);
    if(edges[edge_ind].right_orient == _or){
      al_start = edges[edge_ind].right_al_sites[0];
      al_end = edges[edge_ind].right_al_sites[sites-1];
    }
    else{
      al_start = nodes[node_ind].map_read.size() - 
	edges[edge_ind].right_al_sites[sites-1];
      al_end = nodes[node_ind].map_read.size() - 
	edges[edge_ind].right_al_sites[0];
    }
  }
  assert(al_start < al_end);
  if(map_site>= al_start && map_site <= al_end) return true;
  else return false;
}

bool graph_component::chimeric_map1(int ind){
  assert(ind >= 0 && ind < nodes.size());
  
  int search_depth = 3;
  int end_delta = 4;

  int cur_rightmost_al_site = rightmost_al_site(ind);
  int cur_leftmost_al_site = leftmost_al_site(ind);
  cout<<"for node: "<<ind<<" (";
  cout<<suffix(nodes[ind].read_name,suffix_length);
  cout<<") cur_leftmost_al_site: "<<cur_leftmost_al_site;
  cout<<" cur_rightmost_al_site: "<<cur_rightmost_al_site<<endl;

  int left_site = cur_leftmost_al_site + end_delta;
  int right_site = cur_rightmost_al_site - end_delta;
  
  assert(left_site <= right_site);

  vector<int> left_neighbors;
  vector<int> right_neighbors;

  vector<int> neighbors;
  vector<int> inc_edge_inds;
  for(int i=0; i<nodes[ind].right_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].right_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].right_edges_inds[i]);
  }
  for(int i=0; i<nodes[ind].left_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].left_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].left_edges_inds[i]);
  }
  for(int i=0; i<nodes[ind].contained_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].contained_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].contained_edges_inds[i]);
  }
  for(int i=0; i<nodes[ind].container_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].container_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].container_edges_inds[i]);
  }

  bool both_sites_spanned = false;
  for(int i=0; i<neighbors.size(); i++){
    int cur_edge_ind = inc_edge_inds[i];
    bool covers_left_site = 
      overlap_spans_site(ind, left_site, nodes[ind].orient, cur_edge_ind);
    bool covers_right_site = 
      overlap_spans_site(ind, right_site, nodes[ind].orient, cur_edge_ind);

    if(covers_left_site) left_neighbors.push_back(neighbors[i]);
    if(covers_right_site) right_neighbors.push_back(neighbors[i]);
    if(covers_left_site && covers_right_site) both_sites_spanned = true;
  }

  cout<<"for node: "<<suffix(nodes[ind].read_name, suffix_length)<<endl;
  cout<<"left neighbors: "<<endl;
  for(int i=0; i<left_neighbors.size(); i++){
    cout<<char(9)<<suffix(nodes[left_neighbors[i]].read_name, suffix_length)<<endl;
  }
  cout<<endl;
  cout<<"right neighbors: "<<endl;
  for(int i=0; i<right_neighbors.size(); i++){
    cout<<char(9)<<suffix(nodes[right_neighbors[i]].read_name, suffix_length)<<endl;
  }
  cout<<endl;

  if(both_sites_spanned) return false; //for full containments

  if(right_neighbors.empty() || left_neighbors.empty()) return false;

  vector<bool> discovered_nodes;
  for(int i=0; i<nodes.size(); i++){
    discovered_nodes.push_back(false);
  }
  discovered_nodes[ind] = true; //to mark that it is not searched through

  for(int i=0; i<left_neighbors.size(); i++){
    int cur_neighbor = left_neighbors[i];
    cout<<"cur_left_neighbor: "<<cur_neighbor<<endl;
    not_left_directed_dfs(cur_neighbor, discovered_nodes, search_depth);
  }
  
  //check if any right nodes were reached
  for(int i=0; i<right_neighbors.size(); i++){
    int cur_right_neighbor = right_neighbors[i];
    if(discovered_nodes[cur_right_neighbor]) return false;
  }
  return true; //nothing was discovered, it must be a chimera
}

bool graph_component::chimeric_map2(int ind){

  assert(ind >= 0 && ind < nodes.size());
  
  int search_depth = 7;//12

  vector<int> neighbors;
  vector<int> inc_edge_inds;
  for(int i=0; i<nodes[ind].right_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].right_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].right_edges_inds[i]);
  }
  for(int i=0; i<nodes[ind].left_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].left_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].left_edges_inds[i]);
  }
  for(int i=0; i<nodes[ind].contained_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].contained_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].contained_edges_inds[i]);
  }
  for(int i=0; i<nodes[ind].container_maps_inds.size(); i++){
    neighbors.push_back(nodes[ind].container_maps_inds[i]);
    inc_edge_inds.push_back(nodes[ind].container_edges_inds[i]);
  }

  assert(neighbors.size() > 0);
  
  bool map_chimeric = true;
  for(int n=0; n<neighbors.size(); n++){
    vector<bool> discovered_nodes;
    for(int i=0; i<nodes.size(); i++){
      discovered_nodes.push_back(false);
    }
    discovered_nodes[ind] = true; //to mark that it is not searched through
    
    int cur_test_neighbor = neighbors[n];

    breadth_first_search(cur_test_neighbor, discovered_nodes,search_depth);

    //check if any nodes were not reached
    bool reached_all = true;
    for(int i=0; i<neighbors.size(); i++){
      int cur_neighbor = neighbors[i];
      if(!discovered_nodes[cur_neighbor]){ 
	reached_all = false;
	if(suffix(nodes[ind].read_name,suffix_length) == "4230_0_273"){
	cerr<<"not reached: "<<suffix(nodes[cur_neighbor].read_name,suffix_length);
	cerr<<" from: "<<suffix(nodes[cur_test_neighbor].read_name,suffix_length)<<endl;
	}
      }
    }
    if(reached_all) map_chimeric = false;
  }
  return map_chimeric; //nothing was discovered, it must be a chimera
}

void graph_component::
breadth_first_search(int from_node, vector<bool>& discovered_nodes, int search_depth){
  //cout<<"from_node: "<<from_node<<" depth: "<<search_depth<<endl;
  assert(from_node >=0 && from_node < nodes.size());
  if(discovered_nodes[from_node] == true) return;
  else{
    discovered_nodes[from_node]=true; // mark discovered
    if(search_depth<=0) return;
    else{
      vector<int> neighbor_list;
      for(int i=0; i<nodes[from_node].right_maps_inds.size(); i++)
	neighbor_list.push_back(nodes[from_node].right_maps_inds[i]);
      for(int i=0; i<nodes[from_node].left_maps_inds.size(); i++)
	neighbor_list.push_back(nodes[from_node].left_maps_inds[i]);
      for(int i=0; i<nodes[from_node].contained_maps_inds.size(); i++)
	neighbor_list.push_back(nodes[from_node].contained_maps_inds[i]);
      for(int i=0; i<nodes[from_node].container_maps_inds.size(); i++)
	neighbor_list.push_back(nodes[from_node].container_maps_inds[i]);


      for(int i=0; i<neighbor_list.size(); i++){
	int cur_neighbor = neighbor_list[i];
	if(!discovered_nodes[cur_neighbor])
	  breadth_first_search(cur_neighbor, discovered_nodes, search_depth-1);
      }
      return;
    }
  }
}



void graph_component::
not_left_directed_dfs(int from_node, vector<bool>& discovered_nodes, int search_depth){
  //cout<<"from_node: "<<from_node<<" depth: "<<search_depth<<endl;
  assert(from_node >=0 && from_node < nodes.size());
  if(discovered_nodes[from_node] == true) return;
  else{
    discovered_nodes[from_node]=true; // mark discovered
    if(search_depth<=0) return;
    else{
      vector<int> neighbor_list;
      for(int i=0; i<nodes[from_node].right_maps_inds.size(); i++)
	neighbor_list.push_back(nodes[from_node].right_maps_inds[i]);
      for(int i=0; i<nodes[from_node].contained_maps_inds.size(); i++)
	neighbor_list.push_back(nodes[from_node].contained_maps_inds[i]);
      for(int i=0; i<nodes[from_node].container_maps_inds.size(); i++)
	neighbor_list.push_back(nodes[from_node].container_maps_inds[i]);

      for(int i=0; i<neighbor_list.size(); i++){
	int cur_neighbor = neighbor_list[i];
	not_left_directed_dfs(cur_neighbor, discovered_nodes, search_depth-1);
      }
      return;
    }
  }
}

vector<double> graph_component::
extract_contig1(vector<int>& map_inds){
  vector<double> contig;

  int last_al_map_site = -1;

  cout<<"contig extraction routine"<<endl;

  vector<site_ind> last_step_cur_map_al_sites;
  vector<site_ind> last_step_prev_map_al_sites;
  int last_step_cur_map_ind;
  int last_step_prev_map_ind;
  orientation last_step_cur_map_orient;
  orientation last_step_prev_map_orient;

  int last_step_last_al_site_ind=200;
  int al_ind;

  for(int i=1; i<map_inds.size(); i++){
    int cur_map_ind = map_inds[i];
    int prev_map_ind = map_inds[i-1];
    int cur_edge_ind = nodes[cur_map_ind].get_edge_ind(prev_map_ind);
    //int cur_edge_ind = edge_ind(cur_map_ind, prev_map_ind);
    assert(cur_edge_ind != UNDEF_IND);

    cout<<"processing overlap : "<<prev_map_ind<<"->"<<cur_map_ind<<endl;

    string cur_map_name = nodes[cur_map_ind].read_name;
    string prev_map_name = nodes[prev_map_ind].read_name;

    orientation cur_map_orient;
    orientation prev_map_orient;

    int cur_map_left_al_site;
    int cur_map_right_al_site;
    int prev_map_left_al_site;
    int prev_map_right_al_site;

    vector<site_ind> prev_map_al_sites;
    vector<site_ind> cur_map_al_sites;

    int al_sites = edges[cur_edge_ind].left_al_sites.size();

    prev_map_orient = nodes[prev_map_ind].orient;
    cur_map_orient = nodes[cur_map_ind].orient;

    if(prev_map_name == edges[cur_edge_ind].left_map){   
      cout<<"prev-left, cur-right"<<endl;
      assert(cur_map_name == edges[cur_edge_ind].right_map);     

      if(nodes[prev_map_ind].orient == edges[cur_edge_ind].left_orient){
	//same orientation
	assert(nodes[cur_map_ind].orient == edges[cur_edge_ind].right_orient);

	prev_map_left_al_site = edges[cur_edge_ind].left_al_sites[0];
	prev_map_right_al_site = edges[cur_edge_ind].left_al_sites[al_sites-1];
	
	cur_map_left_al_site = edges[cur_edge_ind].right_al_sites[0];
	cur_map_right_al_site = edges[cur_edge_ind].right_al_sites[al_sites-1];     
	
	prev_map_al_sites = edges[cur_edge_ind].left_al_sites;
	cur_map_al_sites = edges[cur_edge_ind].right_al_sites;
      }
      else{
	//orientation different, need to reverse

	assert(nodes[cur_map_ind].orient == rev(edges[cur_edge_ind].right_orient));

	prev_map_left_al_site = nodes[prev_map_ind].map_read.size()
	  - edges[cur_edge_ind].left_al_sites[al_sites-1];
	prev_map_right_al_site = nodes[prev_map_ind].map_read.size()
	  - edges[cur_edge_ind].left_al_sites[0];
	
	cur_map_left_al_site = nodes[cur_map_ind].map_read.size()
	  - edges[cur_edge_ind].right_al_sites[al_sites-1];
	cur_map_right_al_site = nodes[cur_map_ind].map_read.size()
	  - edges[cur_edge_ind].right_al_sites[0];     
	
	prev_map_al_sites.clear();
	cur_map_al_sites.clear();

	for(int j=edges[cur_edge_ind].left_al_sites.size()-1; j>=0; j--){
	  prev_map_al_sites.push_back(nodes[prev_map_ind].map_read.size()-
				      edges[cur_edge_ind].left_al_sites[j]);
	  cur_map_al_sites.push_back(nodes[cur_map_ind].map_read.size()-
				     edges[cur_edge_ind].right_al_sites[j]);
	}
      }
    }
    else{
      cout<<"prev-right, cur-left"<<endl;
      assert(cur_map_name == edges[cur_edge_ind].left_map);

      if(nodes[prev_map_ind].orient == edges[cur_edge_ind].right_orient){
	assert(nodes[cur_map_ind].orient == edges[cur_edge_ind].left_orient);
	assert(false);
      }
      else{
	assert(nodes[cur_map_ind].orient == rev(edges[cur_edge_ind].left_orient));
	prev_map_left_al_site = nodes[prev_map_ind].map_read.size() -
	  edges[cur_edge_ind].right_al_sites[al_sites-1];
	prev_map_right_al_site = nodes[prev_map_ind].map_read.size() - 
	  edges[cur_edge_ind].right_al_sites[0];
	
	cur_map_left_al_site = nodes[cur_map_ind].map_read.size() - 
	  edges[cur_edge_ind].left_al_sites[al_sites-1];
	cur_map_right_al_site = nodes[cur_map_ind].map_read.size() -
	  edges[cur_edge_ind].left_al_sites[0];
	
	prev_map_al_sites.clear();
	cur_map_al_sites.clear();
	
	for(int j=edges[cur_edge_ind].right_al_sites.size()-1; j>=0; j--){
	  prev_map_al_sites.push_back(nodes[prev_map_ind].map_read.size() - 
				      edges[cur_edge_ind].right_al_sites[j]);
	  cur_map_al_sites.push_back(nodes[cur_map_ind].map_read.size() - 
				     edges[cur_edge_ind].left_al_sites[j]);
	}
      }
    }

    cout<<"prev_or:"<<prev_map_orient<<" cur_or:"<<cur_map_orient<<endl;
    
    assert(cur_map_left_al_site < cur_map_right_al_site);
    assert(prev_map_left_al_site < prev_map_right_al_site);

    cout<<"al: ["<<prev_map_left_al_site<<","<<prev_map_right_al_site<<"]->[";
    cout<<cur_map_left_al_site<<","<<cur_map_right_al_site<<"]"<<endl;
    for(int j=0; j<prev_map_al_sites.size(); j++){
      cout<<" "<<prev_map_al_sites[j];
    }
    cout<<endl;
    for(int j=0; j<prev_map_al_sites.size(); j++){
      cout<<" "<<cur_map_al_sites[j];
    }
    cout<<endl<<endl;

    if(i==1){
      //add the hanging part of left map
      //for(int j=0; j<prev_map_left_al_site; j++){
      //double cur_fr;
      //if(prev_map_orient == forward) cur_fr = nodes[prev_map_ind].map_read[j];
      //else cur_fr = nodes[prev_map_ind].map_read[nodes[prev_map_ind].map_read.size()-j-1];
	//contig.push_back(cur_fr);
      //}
      
      for(int j=leftmost_al_site(prev_map_ind); j<prev_map_left_al_site; j++){
	double cur_fr;
	if(prev_map_orient == forward) cur_fr = nodes[prev_map_ind].map_read[j];
	else cur_fr = nodes[prev_map_ind].map_read[nodes[prev_map_ind].map_read.size()-j-1];
	contig.push_back(cur_fr);
      }

      last_al_map_site = cur_map_left_al_site;

      cout<<"overhang: from 0 to "<<prev_map_left_al_site<<endl;

      last_step_last_al_site_ind = 0;
      al_ind = 0;
    }
    else{
      cout<<"last_al_site:"<<last_al_map_site<<endl;

      //extending
      assert(i>1);
      
      assert(last_al_map_site >= 0 && last_al_map_site < nodes[prev_map_ind].map_read.size());

      int next_prev_map_site = -1;
      bool next_prev_map_site_found = false;
      //int al_ind;

      for(int j=0; j<prev_map_al_sites.size() && !next_prev_map_site_found; j++){
	int cur_prev_map_site = prev_map_al_sites[j];
	cout<<"cur_prev_map_site:"<<cur_prev_map_site<<endl;
	if(cur_prev_map_site >= last_al_map_site){
	  next_prev_map_site_found = true;
	  next_prev_map_site = cur_prev_map_site;
	  al_ind = j;
	}
      }
      assert(next_prev_map_site_found);


      for(int j=last_al_map_site; j<next_prev_map_site; j++){
	double cur_fr;
	if(prev_map_orient == forward) cur_fr = nodes[prev_map_ind].map_read[j];
	else cur_fr = nodes[prev_map_ind].map_read[nodes[prev_map_ind].map_read.size()-j-1];
	contig.push_back(cur_fr);
      }
      
      assert(al_ind >= 0 && al_ind <cur_map_al_sites.size());
      last_al_map_site = cur_map_al_sites[al_ind];
    }

    if(i==map_inds.size()-1){
      //first append the overlap part, then append the overhanging (non-overlapping part)

      //last overlap, append the overhang of the right map
      assert(last_al_map_site <= nodes[cur_map_ind].map_read.size());

      for(int j=last_al_map_site; j<cur_map_right_al_site; j++){
	double cur_fr;
	if(cur_map_orient == forward) cur_fr = nodes[cur_map_ind].map_read[j];
	else cur_fr = nodes[cur_map_ind].map_read[nodes[cur_map_ind].map_read.size()-j-1];
	contig.push_back(cur_fr);
      }
      
      last_al_map_site = cur_map_right_al_site;
      cout<<"last_al_map_site: "<<last_al_map_site<<" rightmost: "<<rightmost_al_site(cur_map_ind)<<endl;
      
      for(int j=last_al_map_site; j<rightmost_al_site(cur_map_ind); j++){
	double cur_fr;
	if(cur_map_orient == forward) cur_fr = nodes[cur_map_ind].map_read[j];
	else cur_fr = nodes[cur_map_ind].map_read[nodes[cur_map_ind].map_read.size()-j-1];
	contig.push_back(cur_fr);
      }

      //last_al_map_site = rightmost_al_site(cur_map_ind);
      //for(int j=last_al_map_site; j<nodes[cur_map_ind].map_read.size(); j++){
      //double cur_fr;
      //if(cur_map_orient == forward) cur_fr = nodes[cur_map_ind].map_read[j];
      //else cur_fr = nodes[cur_map_ind].map_read[nodes[cur_map_ind].map_read.size()-j-1];
      //contig.push_back(cur_fr);
      //}
    }

    last_step_cur_map_al_sites = cur_map_al_sites;
    last_step_prev_map_al_sites = prev_map_al_sites;
    last_step_cur_map_ind = cur_map_ind;
    last_step_prev_map_ind = prev_map_ind;
    last_step_cur_map_orient = cur_map_orient;
    last_step_prev_map_orient = prev_map_orient;
    last_step_last_al_site_ind = al_ind;
  }

  return contig;
}

vector<int> good_paths(vector<path>& paths){

  double std_num = 2.0;//1.0;
  double sigma = 0.5579;
  vector<int> best_path_set;
  //returns the set of inds of feasible paths
 
  for(int i=0; i<paths.size(); i++){
    double cur_path_length = paths[i].get_length();
    double upper_bound = cur_path_length + std_num*sigma*sqrt(fabs(cur_path_length));
    double lower_bound = cur_path_length - std_num*sigma*sqrt(fabs(cur_path_length));

    vector<int> cur_matching_set;
    for(int j=0; j<paths.size(); j++){
      double other_path_length = paths[i].get_length();
      if(other_path_length < upper_bound &&
	 other_path_length > lower_bound){
	cur_matching_set.push_back(j);
      }
    }
    if(cur_matching_set.size() > best_path_set.size() &&
       cur_matching_set.size() > 1){
      best_path_set = cur_matching_set;
    }
  }
  
  //debugging
  /*
  for(int i=0; i<paths.size(); i++){
    for(int j=0; j<paths[i].edge_weights.size(); j++){
      if((int)paths[i].edge_weights[j]==107){
	for(int k=0; k<paths.size(); k++){
	  cout<<"path "<<k<<": ";
	  for(int m=0; m<paths[k].edge_inds.size(); m++){
	    cout<<paths[k].edge_inds[m]<<" ";
	  }
	  cout<<" w:";
	  for(int m=0; m<paths[k].edge_weights.size(); m++){
	    cout<<" "<<paths[k].edge_weights[m];
	  }
	  cout<<endl;
	}
	cout<<endl;
	cout<<"best path set:";
	for(int l=0; l<best_path_set.size(); l++){
	  cout<<" "<<l;
	}
	cout<<endl;
	assert(false);
      }
    }
  }
  */
  //end debugging
  return best_path_set;
}

graph_component::graph_component(hash_key _max_hash_key){
  max_hash_key = _max_hash_key;
}

//vector<int> graph_component::successor_list(int cur_node_ind, int dfs_depth){
  //this function expects the node index lists 
  //be filled out in advance  
  //}

void graph_component::
possible_paths_for_consistency_check
(vector<int>& node_inds, vector<vector<path> >& paths, 
 vector<int>& loc_node_color, vector<int>& loc_pred_list, 
 vector<int>& branch_root_list,
 int mark_color, int depth, int cur_node_ind, int prev_node_ind,
 int master_node_ind, int branch_root_ind){
  assert(cur_node_ind < nodes.size());


  cout<<"processing node:"<<cur_node_ind;
  cout<<" master:"<<master_node_ind;
  cout<<" branch_root:"<<branch_root_ind;
  cout<<endl;

  if(loc_node_color[cur_node_ind] == mark_color){
    //node was discovered before
    if(branch_root_list[cur_node_ind] == branch_root_ind){
      //this is not forked right after the master node, but later
      cout<<"branch discovered before, but has the same branch root"<<endl;
      return;
    }

    cout<<"node was discovered before and has different branch root"<<endl;
    //node has been discovered before
    //add to the path list
    
    //look if the path to the current node has already been added
    bool node_was_added = false;
    int node_ind_in_the_list = -1;
    for(int i=0; i<node_inds.size(); i++){
      int cur_listed_node = node_inds[i];
      if(cur_listed_node == cur_node_ind){
	node_was_added = true;
	node_ind_in_the_list = i;
      }
    }

    if(node_was_added){
      cout<<"node was added before"<<endl;
      //no need to add a node to the list,
      //just add the current path
      assert(node_ind_in_the_list >=0 && node_ind_in_the_list < node_inds.size());
    }
    else{
      //need to add the node to the list
      //and need to add both paths to the list of paths
      
      cout<<"node not added before. Tracing back to master."<<endl;

      int sanity_check_counter = 0;
      
      vector<int> inv_path_node_inds;
      vector<double> inv_path_weights;
      
      int prev_trace_back_node_ind = loc_pred_list[cur_node_ind];
      int cur_trace_back_node_ind = cur_node_ind;

      bool master_node_reached = false;
      if(cur_trace_back_node_ind == master_node_ind){
	master_node_reached = true;
	inv_path_node_inds.push_back(master_node_ind);
      }
      else{
	assert(loc_node_color[prev_trace_back_node_ind] == mark_color);
      }
      while(!master_node_reached){
	cout<<"prev_tr_back_node:"<<prev_trace_back_node_ind<<endl;
	cout<<"cur_tr_back_node:"<<cur_trace_back_node_ind<<endl;
	cout<<endl;


	sanity_check_counter++;
	assert(sanity_check_counter < 100);
	
	//find the edge_weight
	bool corresponding_edge_found = false;
	double cur_edge_weight;
	
	for(int i=0; i<nodes[prev_trace_back_node_ind].right_maps_inds.size() 
	      && !corresponding_edge_found; i++){
	  if(nodes[prev_trace_back_node_ind].right_maps_inds[i] == 
	     cur_trace_back_node_ind){
	    corresponding_edge_found = true;
	    assert(nodes[prev_trace_back_node_ind].right_edge_weights.size() >= i);
	    cur_edge_weight = nodes[prev_trace_back_node_ind].right_edge_weights[i];
	  }
	}
	for(int i=0; i<nodes[prev_trace_back_node_ind].contained_maps_inds.size()
	      && !corresponding_edge_found; i++){
	  if(nodes[prev_trace_back_node_ind].contained_maps_inds[i] == 
	     cur_trace_back_node_ind){
	    corresponding_edge_found = true;
	    assert(nodes[prev_trace_back_node_ind].contained_edge_weights.size() >= i);
	    cur_edge_weight = 
	      nodes[prev_trace_back_node_ind].contained_edge_weights[i];
	  }
	}
	for(int i=0; i<nodes[prev_trace_back_node_ind].container_maps_inds.size()
	      && !corresponding_edge_found; i++){
	  if(nodes[prev_trace_back_node_ind].container_maps_inds[i] == 
	     cur_trace_back_node_ind){
	    corresponding_edge_found = true;
	    assert(nodes[prev_trace_back_node_ind].container_edge_weights.size() >= i);
	    cur_edge_weight = 
	      nodes[prev_trace_back_node_ind].container_edge_weights[i];
	  }
	}

	assert(corresponding_edge_found);
	
	inv_path_node_inds.push_back(cur_trace_back_node_ind);
	inv_path_weights.push_back(cur_edge_weight);
	
	if(prev_trace_back_node_ind == master_node_ind){
	  master_node_reached = true;
	  inv_path_node_inds.push_back(prev_trace_back_node_ind);
	}
	else{
	  cur_trace_back_node_ind = prev_trace_back_node_ind;
	  prev_trace_back_node_ind = loc_pred_list[prev_trace_back_node_ind];
	  assert(loc_node_color[prev_trace_back_node_ind] == mark_color);
	}
      }
      
      path orig_path;
      for(int i=inv_path_node_inds.size()-1; i>=0; i--){
	orig_path.map_inds.push_back(inv_path_node_inds[i]);
      }
      for(int i=inv_path_weights.size()-1; i>=0; i--){
	orig_path.edge_weights.push_back(inv_path_weights[i]);
      }
      orig_path.set_length();

      node_inds.push_back(cur_node_ind);
      vector<path> new_path_vec;
      new_path_vec.push_back(orig_path);
      paths.push_back(new_path_vec);
      
      node_ind_in_the_list = node_inds.size()-1;
    }
    
    int sanity_check_counter = 0;
    
    vector<int> inv_path_node_inds;
    vector<double> inv_path_weights;
    
    int prev_trace_back_node_ind = prev_node_ind;
    int cur_trace_back_node_ind = cur_node_ind;
    bool master_node_reached = false;
    
    if(cur_trace_back_node_ind == master_node_ind){
      master_node_reached = true;
      inv_path_node_inds.push_back(master_node_ind);
    }

    while(!master_node_reached){
      sanity_check_counter++;
      assert(sanity_check_counter < 100);
      
      //find the edge_weight
      bool corresponding_edge_found = false;
      double cur_edge_weight;
      
      for(int i=0; i<nodes[prev_trace_back_node_ind].right_maps_inds.size() 
	    && !corresponding_edge_found; i++){
	if(nodes[prev_trace_back_node_ind].right_maps_inds[i] == 
	   cur_trace_back_node_ind){
	  corresponding_edge_found = true;
	  assert(nodes[prev_trace_back_node_ind].right_edge_weights.size() >= i);
	  cur_edge_weight = nodes[prev_trace_back_node_ind].right_edge_weights[i];
	}
      }
      for(int i=0; i<nodes[prev_trace_back_node_ind].contained_maps_inds.size()
	    && !corresponding_edge_found; i++){
	if(nodes[prev_trace_back_node_ind].contained_maps_inds[i] == 
	   cur_trace_back_node_ind){
	  corresponding_edge_found = true;
	  assert(nodes[prev_trace_back_node_ind].contained_edge_weights.size() >= i);
	  cur_edge_weight = nodes[prev_trace_back_node_ind].contained_edge_weights[i];
	}
      }
      for(int i=0; i<nodes[prev_trace_back_node_ind].container_maps_inds.size()
	    && !corresponding_edge_found; i++){
	if(nodes[prev_trace_back_node_ind].container_maps_inds[i] == 
	   cur_trace_back_node_ind){
	  corresponding_edge_found = true;
	  assert(nodes[prev_trace_back_node_ind].container_edge_weights.size() >= i);
	  cur_edge_weight = nodes[prev_trace_back_node_ind].container_edge_weights[i];
	}
      }

      assert(corresponding_edge_found);
      
      inv_path_node_inds.push_back(cur_trace_back_node_ind);
      inv_path_weights.push_back(cur_edge_weight);
      
      if(prev_trace_back_node_ind == master_node_ind){
	master_node_reached = true;
	inv_path_node_inds.push_back(prev_trace_back_node_ind);
      }
      else{
	cur_trace_back_node_ind = prev_trace_back_node_ind;
	prev_trace_back_node_ind = loc_pred_list[prev_trace_back_node_ind];
	assert(loc_node_color[prev_trace_back_node_ind] == mark_color);
      }
    }
    
    path new_path;
    for(int i=inv_path_node_inds.size()-1; i>=0; i--){
      new_path.map_inds.push_back(inv_path_node_inds[i]);
    }
    for(int i=inv_path_weights.size()-1; i>=0; i--){
      new_path.edge_weights.push_back(inv_path_weights[i]);
    }
    
    new_path.set_length();
    assert(paths.size() > node_ind_in_the_list);
    paths[node_ind_in_the_list].push_back(new_path);
    
    //now add the new path
    return;
  }
  else{
    //recursive step
    loc_node_color[cur_node_ind] = mark_color;
    branch_root_list[cur_node_ind] = branch_root_ind;
    loc_pred_list[cur_node_ind] = prev_node_ind;
 
    if(depth == 0) return; //depth exhausted


    vector<int> succ;
    for(int i=0; i<nodes[cur_node_ind].right_maps_inds.size(); i++){
      succ.push_back(nodes[cur_node_ind].right_maps_inds[i]);
    }
    for(int i=0; i<nodes[cur_node_ind].contained_maps_inds.size(); i++){
      succ.push_back(nodes[cur_node_ind].contained_maps_inds[i]); 
    }
    //changed
    for(int i=0; i<nodes[cur_node_ind].container_maps_inds.size(); i++){
      succ.push_back(nodes[cur_node_ind].container_maps_inds[i]); 
    }

    cout<<"succ:";
    for(int i=0; i<succ.size(); i++){
      cout<<" "<<succ[i];
    }
    cout<<endl;

    if(cur_node_ind == master_node_ind){
      //the search has just started
      //all successors are roots of branches
      
      for(int i=0; i<succ.size(); i++){
	int cur_succ = succ[i];
	if(cur_succ != prev_node_ind){ //don't go back along the same edge
	  possible_paths_for_consistency_check
	    (node_inds, paths, loc_node_color, loc_pred_list, branch_root_list,
	     mark_color, depth-1, cur_succ, cur_node_ind, master_node_ind,
	     cur_succ);
	}
      }
    }
    else{
      //this is beyond the branch root
      //just keep going
      for(int i=0; i<succ.size(); i++){
	int cur_succ = succ[i];
       
	possible_paths_for_consistency_check
	  (node_inds, paths, loc_node_color, loc_pred_list, branch_root_list,
	   mark_color, depth-1, cur_succ, cur_node_ind, master_node_ind,
	   branch_root_ind);
      }
    }
  }
}
vector<graph_component> graph_component::
confirmed_components(int dfs_depth){
  vector<graph_component> new_components;

  //first make a list of highly confident edges
  //based on the distance consistency check
  
  //to do this go through all nodes and do dfs 
  //search of depth dfs_depth to locate pairs
  //of nodes connected with more than one path
  //with consistent distances. Select the 
  //edges that contribute to at least
  //two paths with consistent distances
  cerr<<"confirming edges ... ";
  list<int> new_edge_set = confirmed_edges(dfs_depth);
  cerr<<"finished"<<endl;

  new_edge_set.sort();
  //sort(new_edge_set.begin(), new_edge_set.end());
  
  //vector<vector< string > > map_names;
  vector<vector< hash_key > > map_hash_keys;
  vector<vector< orientation > > map_orients;
  vector<vector< int > > new_edges;
  list<int> consistent_edge_set = 
    consistent_edge_inds(new_edge_set, edges,
			 map_hash_keys, map_orients, new_edges);
    //  consistent_edge_inds(new_edge_set, map_hash_keys, map_orients, new_edges);

  if(consistent_edge_set.size() != new_edge_set.size()){
    cout<<"consistent edge set:"<<endl;
    for(list<int>::iterator it = consistent_edge_set.begin();
	it!=consistent_edge_set.end(); it++){
      cout<<*it<<endl;
    }
    /*
    for(int i=0; i<consistent_edge_set.size(); i++){
      cout<<consistent_edge_set[i]<<endl;
    }
    */
    cout<<endl;
    cout<<"new edge set:"<<endl;
    for(list<int>::iterator it = new_edge_set.begin(); 
	it!= new_edge_set.end(); it++){
      cout<<(*it)<<endl;
    }
    /*
    for(int i=0; i<new_edge_set.size(); i++){
      cout<<new_edge_set[i]<<endl;
    }
    */
    cout<<endl;
    assert(false);
  }

  assert(map_hash_keys.size() == map_orients.size());
  assert(map_hash_keys.size() == new_edges.size());

  for(int i=0; i<map_hash_keys.size(); i++){
    cout<<"processing "<<i<<" component"<<endl;
    assert(map_hash_keys[i].size() == map_orients[i].size());

    graph_component new_comp(max_hash_key);
    for(int j=0; j<map_hash_keys[i].size(); j++){
      hash_key cur_map_hash_key = map_hash_keys[i][j];
      string cur_map_name = return_map_name(cur_map_hash_key);
      //int cur_node_index = get_node_index(cur_map_name);
      vector<fr_size> cur_map_read = return_map_read(cur_map_hash_key);
      
      //cout<<"adding a map ("<<cur_node_index<<") "<<cur_map_name<<endl;
      //assert(cur_node_index != -1);

      //node new_node(nodes[cur_node_index].hk, nodes[cur_node_index].read_name, 
      //	    map_orients[i][j], nodes[cur_node_index].map_read);
      //new_node.id = nodes[cur_node_index].id;

      node new_node(cur_map_hash_key, cur_map_name, map_orients[i][j], cur_map_read);
      new_comp.nodes.push_back(new_node);
    }
    
    for(int j=0; j<new_edges[i].size(); j++){
      int cur_edge_ind = new_edges[i][j];
      cout<<"adding an edge:"<<cur_edge_ind<<endl;
      assert(cur_edge_ind >= 0 && cur_edge_ind<edges.size());
      new_comp.edges.push_back(edges[cur_edge_ind]);    
    }
    cout<<"filling out the component structure."<<endl;
    new_comp.fill_the_component_data();
    new_comp.assign_edge_weights();
    cout<<"structure completed"<<endl;
    new_components.push_back(new_comp);
  }
  return new_components;
}

void graph_component::fill_the_component_data(){
  make_hash_table();
  for(int i=0; i<nodes.size(); i++){

    nodes[i].right_maps.clear();
    nodes[i].right_maps_inds.clear();
    nodes[i].right_edges_inds.clear();
    nodes[i].right_edge_weights.clear();

    nodes[i].left_maps.clear();
    nodes[i].left_maps_inds.clear();
    nodes[i].left_edges_inds.clear();
    nodes[i].left_edge_weights.clear();

    nodes[i].contained_maps.clear();
    nodes[i].contained_maps_inds.clear();
    nodes[i].contained_edges_inds.clear();
    nodes[i].contained_edge_weights.clear();

    nodes[i].container_maps.clear();
    nodes[i].container_maps_inds.clear();
    nodes[i].container_edges_inds.clear();
    nodes[i].container_edge_weights.clear();
  }
  for(int i=0; i<edges.size(); i++){
    if (i%100 == 0) 
      cout<<"fill_the_component_data: processed "<<i<<" edges of "<<edges.size()<<endl;

    hash_key cur_left_map_hk = edges[i].left_map_hash_key;
    hash_key cur_right_map_hk = edges[i].right_map_hash_key;

    //string cur_left_map = edges[i].left_map;
    //string cur_right_map = edges[i].right_map;

    orientation left_orient = edges[i].left_orient;
    orientation right_orient = edges[i].right_orient;

    //int cur_left_map_ind = get_node_index(cur_left_map);
    //int cur_right_map_ind = get_node_index(cur_right_map);

    int cur_left_map_ind = get_node_index(cur_left_map_hk);
    int cur_right_map_ind = get_node_index(cur_right_map_hk);

    assert(cur_left_map_ind != UNDEF_IND);
    assert(cur_right_map_ind != UNDEF_IND);

    //debug begin    
    /*
    for(int j=0; j<nodes.size(); j++){
      cerr<<"node: "<<nodes[j].read_name<<" hk: "<<nodes[j].hk<<endl;
    }
    for(int j=0; j<edges.size(); j++){
      cerr<<edges[j].left_map<<" ("<<edges[j].left_map_hash_key<<") -> ";
      cerr<<edges[j].right_map<<" ("<<edges[j].right_map_hash_key<<")"<<endl;
    }
    */
    //cerr<<edges[0].left_map<<" -> "<<nodes[cur_left_map_ind].read_name<<endl;
    //cerr<<edges[0].right_map<<" -> "<<nodes[cur_right_map_ind].read_name<<endl;
    //cerr<<edges[1].left_map<<" -> "<<nodes[cur_left_map_ind].read_name<<endl;
    //cerr<<edges[1].right_map<<" -> "<<nodes[cur_right_map_ind].read_name<<endl;
    //debug end
    string cur_left_map = nodes[cur_left_map_ind].read_name;
    string cur_right_map = nodes[cur_right_map_ind].read_name;

    if(edges[i].left_map != cur_left_map ||
       edges[i].right_map != cur_right_map){
      cerr<<i<<" "<<cur_left_map_hk<<" "<<cur_left_map_ind;
      cerr<<" "<<cur_left_map<<" ("<<edges[i].left_map<<")"<<endl;
      cerr<<i<<" "<<cur_right_map_hk<<" "<<cur_right_map_ind;
      cerr<<" "<<cur_right_map<<" ("<<edges[i].right_map<<")"<<endl;
    }
    
    assert(edges[i].left_map == cur_left_map);//nodes[cur_left_map_ind].read_name);
    assert(edges[i].right_map == cur_right_map);//nodes[cur_right_map_ind].read_name);

    assert(cur_left_map_ind != UNDEF_IND);
    assert(cur_right_map_ind != UNDEF_IND);

    if(!edges[i].containment){
      //non-contained edge
      if(left_orient == nodes[cur_left_map_ind].orient){
	//maps in the same orientation in the component as in the overlap
	assert(right_orient == nodes[cur_right_map_ind].orient); //must be consistent

	nodes[cur_left_map_ind].right_maps.push_back(cur_right_map);
	nodes[cur_left_map_ind].right_maps_inds.push_back(cur_right_map_ind);
	nodes[cur_left_map_ind].right_edges_inds.push_back(i);

	nodes[cur_right_map_ind].left_maps.push_back(cur_left_map);
	nodes[cur_right_map_ind].left_maps_inds.push_back(cur_left_map_ind);
	nodes[cur_right_map_ind].left_edges_inds.push_back(i);
      }
      else{
	assert(right_orient == rev(nodes[cur_right_map_ind].orient)); //must be consistent
	
	nodes[cur_left_map_ind].left_maps.push_back(cur_right_map);
	nodes[cur_left_map_ind].left_maps_inds.push_back(cur_right_map_ind);
	nodes[cur_left_map_ind].left_edges_inds.push_back(i);

	nodes[cur_right_map_ind].right_maps.push_back(cur_left_map);
	nodes[cur_right_map_ind].right_maps_inds.push_back(cur_left_map_ind);
	nodes[cur_right_map_ind].right_edges_inds.push_back(i);
      }
    }
    else{
      //containment edge
      nodes[cur_left_map_ind].contained_maps.push_back(cur_right_map);
      nodes[cur_left_map_ind].contained_maps_inds.push_back(cur_right_map_ind);
      nodes[cur_left_map_ind].contained_edges_inds.push_back(i);
      
      nodes[cur_right_map_ind].container_maps.push_back(cur_left_map);
      nodes[cur_right_map_ind].container_maps_inds.push_back(cur_left_map_ind);
      nodes[cur_right_map_ind].container_edges_inds.push_back(i);
    }
  }
  index_nodes();
  assign_edge_weights();
  remove_hash_table();
}

int string_index_in_the_set(string cur_string, vector<string> strings){
  int failure_code = -1;

  for(int i=0; i<strings.size(); i++){
    if(cur_string == strings[i])
      return i;
  }
  return failure_code;//map was not found
}

void graph_component::
dfs_all_way_visit_to_separate_components(int map_ind, vector<string>& map_set,
					 vector<int>& map_colors, vector<int>& edge_inds, vector<int>& edge_colors, int cur_color){
  assert(cur_color >= 0);
  assert(map_ind >= 0 && map_ind < map_set.size());
  assert(map_colors.size() == map_set.size());
  assert(edge_colors.size() == edge_inds.size());


  if(map_colors[map_ind] != -1){
    //map was discovered before
    //has to be of the same color since has to be in the same component
    assert(map_colors[map_ind] == cur_color);

    return;
  }
  else{
    //map was never discovered before
    map_colors[map_ind] = cur_color;
    
    //now collect all the adjacent maps
    vector<int> adj_map_inds;
    for(int i=0; i<edge_inds.size(); i++){
      int cur_edge_ind = edge_inds[i];
      string cur_left_map = edges[cur_edge_ind].left_map;
      string cur_right_map = edges[cur_edge_ind].right_map;

      int cur_left_map_ind = string_index_in_the_set(cur_left_map, map_set);
      int cur_right_map_ind = string_index_in_the_set(cur_right_map, map_set);
      
      assert(cur_left_map_ind != -1); //map found
      assert(cur_right_map_ind != -1); //map_found


      if(cur_left_map_ind == map_ind){
	//edge is incident to this node (map)
	assert(cur_right_map_ind != map_ind);
	
	assert(edge_colors[i] == -1 || edge_colors[i] == cur_color);
	//edge shouldn't have been discovered before
	edge_colors[i] = cur_color;
	adj_map_inds.push_back(cur_right_map_ind);
      }
      if(cur_right_map_ind == map_ind){
	//edge is incident to this node (map)
	assert(cur_left_map_ind != map_ind);
	assert(edge_colors[i] == -1 || edge_colors[i] == cur_color);
	//edge shouldn't have been discovered before
	edge_colors[i] = cur_color;
	adj_map_inds.push_back(cur_left_map_ind);
      }
    }

    for(int i=0; i<adj_map_inds.size(); i++){
      int cur_adj_map_ind = adj_map_inds[i];
      dfs_all_way_visit_to_separate_components(cur_adj_map_ind, map_set, map_colors,
					       edge_inds, edge_colors, cur_color);
    }
    
    return;
  }
}

/*
int which_comp(string& map_name, vector<vector<string> >& maps, 
	       vector<vector<orientation> >& map_orients, orientation& map_orient){
  for(int i=0; i<maps.size(); i++){
    for(int j=0; j<maps[i].size(); j++){
      if(map_name == maps[i][j]){
	map_orient = map_orients[i][j];
	return i;
      }
    }
  }
  map_orient = undef_orient;
  return -1; 
}
*/

int which_comp(const hash_key cur_map_hash_key,
	       const vector<hash_key>& map_to_comp_ht,
	       const vector<int>& comp_to_ind_ht){
  assert(cur_map_hash_key>=0 && 
	 cur_map_hash_key < map_to_comp_ht.size() &&
	 cur_map_hash_key < comp_to_ind_ht.size());
  if(map_to_comp_ht[cur_map_hash_key] == UNDEF_HASH_KEY) return UNDEF_IND;
  else return comp_to_ind_ht[map_to_comp_ht[cur_map_hash_key]];
}

/*
struct eqstr{
  bool operator()(const char* s1, const char* s2) const{return strcmp(s1,s2) == 0;}
};
*/
/*
int index_in_comp(hash_key _key, int comp_ind, const vector<vector<hash_key> >& map_hks){
  assert(comp_ind >= 0 && comp_ind < map_hks.size());
  for(int i=0; i<map_hks[comp_ind].size(); i++){
    if(map_hks[comp_ind][i] == _key) return i;
  }
  assert(false);
  return UNDEF_IND;
}
*/
list<int> graph_component::consistent_edge_inds
  (list<int>& good_edge_inds, const vector<edge>& _edges,
   vector<vector<hash_key> >& map_hks,
   vector<vector<orientation> >& map_orients, 
   vector<vector<int> >& edge_comp_inds){
  
  //this function assumes that edges are already sorted according to some criteria
  //so that orig_edge_inds has them in the decreasing order of importance.
  //orientation inconsistent edges will not be included in constructed object
  
  //this function constructs an outline of the graph by filling out
  //the component structure using which the graph will be built on

  vector<hash_key> map_to_comp_ht; //hash table that tells the hash key
  //for the component given a hash key of the map it is in.
  //the components are hashed by the hash key of the first map added to it
  vector<int> comp_to_ind_ht;
  //gives the index of the component based on the hash key of that component
  vector<int> map_ind_in_comp_ht;

  for(int i=0; i<=max_hash_key; i++){
    map_to_comp_ht.push_back(UNDEF_HASH_KEY);
    comp_to_ind_ht.push_back(UNDEF_IND);
    map_ind_in_comp_ht.push_back(UNDEF_IND);
  }

  list<int> new_edge_inds;

  map_hks.clear();        //structure that defines which component a given map is in
  map_orients.clear();    //duplicate for orientation
  edge_comp_inds.clear(); //holds edge_inds for components

  vector<hash_key> comp_hks; //structure that contains hash_key of components

  cerr<<"consistent_edge_inds routine"<<endl;

  //for(int i=0; i<orig_edge_inds.size(); i++){
  int _counter = 0;
  for(list<int>::iterator it = good_edge_inds.begin(); 
      it != good_edge_inds.end(); it++){
    if(_counter%1000==0){
      cerr<<"processed "<<_counter<<" edges out of "<<good_edge_inds.size();
      //orig_edge_inds.size();
      cerr<<" comps: "<<map_hks.size()<<endl;
    }
    _counter++;
    int cur_edge_ind = *it; //orig_edge_inds[i];

    hash_key cur_left_map_hk = _edges[cur_edge_ind].left_map_hash_key;
    hash_key cur_right_map_hk = _edges[cur_edge_ind].right_map_hash_key;

    orientation cur_left_orient = _edges[cur_edge_ind].left_orient;
    orientation cur_right_orient = _edges[cur_edge_ind].right_orient;

    int left_map_comp_ind = which_comp(cur_left_map_hk, map_to_comp_ht, comp_to_ind_ht);
    int right_map_comp_ind = which_comp(cur_right_map_hk, map_to_comp_ht, comp_to_ind_ht);

    orientation left_map_list_orient = UNDEF_ORIENT;
    orientation right_map_list_orient = UNDEF_ORIENT;

    if(left_map_comp_ind != UNDEF_IND){
      //int left_map_ind_in_comp = index_in_comp(cur_left_map_hk, left_map_comp_ind, map_hks);
      int left_map_ind_in_comp = map_ind_in_comp_ht[cur_left_map_hk];
      assert(map_hks[left_map_comp_ind][left_map_ind_in_comp] == cur_left_map_hk);
      left_map_list_orient = map_orients[left_map_comp_ind][left_map_ind_in_comp];      
    }
    if(right_map_comp_ind != UNDEF_IND){
      //int right_map_ind_in_comp = index_in_comp(cur_right_map_hk, right_map_comp_ind, map_hks);
      int right_map_ind_in_comp = map_ind_in_comp_ht[cur_right_map_hk];
      assert(map_hks[right_map_comp_ind][right_map_ind_in_comp] == cur_right_map_hk);
      right_map_list_orient = map_orients[right_map_comp_ind][right_map_ind_in_comp];
    }
    
    if(!(cur_left_orient == forward || cur_left_orient == reverse)) 
      cerr<<_counter<<" "<<cur_left_orient<<endl;
    assert(cur_left_orient == forward || cur_left_orient == reverse);
    assert(cur_right_orient == forward || cur_right_orient == reverse);
    
    if(left_map_comp_ind != UNDEF_IND) assert(valid_orient(left_map_list_orient));
    if(right_map_comp_ind != UNDEF_IND) assert(valid_orient(right_map_list_orient));
    

    if(left_map_comp_ind == UNDEF_IND){
      //left map not present
      if(right_map_comp_ind == UNDEF_IND){
	//right map not present
	//make a new component

	//maintain hash tables begin
	hash_key new_comp_hk = cur_left_map_hk; //component is known by this first map added to it
	map_to_comp_ht[cur_left_map_hk] = new_comp_hk;
	map_to_comp_ht[cur_right_map_hk] = new_comp_hk;
	comp_to_ind_ht[new_comp_hk] = comp_hks.size(); //because will be added to the end
	map_ind_in_comp_ht[cur_left_map_hk] = 0;  
	map_ind_in_comp_ht[cur_right_map_hk] = 1;
	//maintain hash tables end
	
	//make new component begin
	vector<hash_key> new_comp_map_hks;
	vector<orientation> new_comp_map_orients;
	vector<int> new_comp_edge_inds;
	new_comp_map_hks.push_back(cur_left_map_hk);
	new_comp_map_hks.push_back(cur_right_map_hk);
	new_comp_map_orients.push_back(cur_left_orient);
	new_comp_map_orients.push_back(cur_right_orient);
	new_comp_edge_inds.push_back(cur_edge_ind);
	//make new component end
	
	//maintain component structures begin - add a new component
	map_hks.push_back(new_comp_map_hks);
	map_orients.push_back(new_comp_map_orients);
	edge_comp_inds.push_back(new_comp_edge_inds);
	comp_hks.push_back(new_comp_hk);
	//maintain component structures end

	new_edge_inds.push_back(cur_edge_ind);//save a good edge
      }
      else{
	//left map not present, right map present, add to right map comp

	//maintain hash table begin
	hash_key this_comp_hk = comp_hks[right_map_comp_ind];
	map_to_comp_ht[cur_left_map_hk] = this_comp_hk;
	map_ind_in_comp_ht[cur_left_map_hk] = map_hks[right_map_comp_ind].size(); //putting at the end
	//maintain hash table end

	//no need to make a new component because addind to the old

	//maintain component structures begin - add left map to component strucutres
	map_hks[right_map_comp_ind].push_back(cur_left_map_hk);
	if(right_map_list_orient == cur_right_orient)
	  map_orients[right_map_comp_ind].push_back(cur_left_orient);	  
	else map_orients[right_map_comp_ind].push_back(rev(cur_left_orient));
	edge_comp_inds[right_map_comp_ind].push_back(cur_edge_ind);
	//maintain component structures end

	new_edge_inds.push_back(cur_edge_ind); //save a good edge
      }
    }
    else{
      //left map is present in comp
      if(right_map_comp_ind == UNDEF_IND){
	//left map is present right map not present
	//add the right map to the left map comp

	//maintain hash tables begin
	hash_key this_comp_hk = comp_hks[left_map_comp_ind];
	map_to_comp_ht[cur_right_map_hk] = this_comp_hk; //assign right map to the new comp	
	map_ind_in_comp_ht[cur_right_map_hk] = map_hks[left_map_comp_ind].size(); //putting at the end
	//maintain hash tables end

	//no need to make a new component because adding to the old

	//maintain component structures begin - add right map to component structures
	map_hks[left_map_comp_ind].push_back(cur_right_map_hk);
	if(left_map_list_orient == cur_left_orient)
	  map_orients[left_map_comp_ind].push_back(cur_right_orient);
	else map_orients[left_map_comp_ind].push_back(rev(cur_right_orient));
	edge_comp_inds[left_map_comp_ind].push_back(cur_edge_ind);
	//maintain component structures end	

	new_edge_inds.push_back(cur_edge_ind); //save a good edge
      }
      else{
	//left map is present and right map is present
	if(left_map_comp_ind == right_map_comp_ind){
	  //both left and right in the same component, check orientation
	  assert(map_to_comp_ht[cur_left_map_hk] == map_to_comp_ht[cur_right_map_hk]); 
	  //keys must match or else there is a problem
	  
	  bool orientation_consistent = false;
	  if((left_map_list_orient == cur_left_orient && right_map_list_orient == cur_right_orient) ||
	     (left_map_list_orient == rev(cur_left_orient) && 
	      right_map_list_orient == rev(cur_right_orient))) orientation_consistent = true;
	      
	  //no need to maintain hash tables or make new components

	  //maintain component structures begin
	  if(orientation_consistent) edge_comp_inds[left_map_comp_ind].push_back(cur_edge_ind);
	  //maintain component structures end

	  if(orientation_consistent) new_edge_inds.push_back(cur_edge_ind); //save a good edge
	}
	else{
	  //both maps exist, but in different components, need to merge: 
	  //erase both components and make a new one and add it to the end

	  hash_key new_comp_hk = comp_hks[left_map_comp_ind];

	  //maintain hash tables part begin
	  for(int j=0; j<map_hks[right_map_comp_ind].size(); j++){
	    map_ind_in_comp_ht[map_hks[right_map_comp_ind][j]]=map_hks[left_map_comp_ind].size()+j;
	  }
	  //maintain hash tables part end
	  
	  //make new component begin
	  vector<hash_key> new_comp_map_hks;
	  vector<orientation> new_comp_map_orients;
	  vector<int> new_comp_edge_inds;
	  
	  new_comp_map_hks = map_hks[left_map_comp_ind];
	  new_comp_map_orients = map_orients[left_map_comp_ind];
	  new_comp_edge_inds = edge_comp_inds[left_map_comp_ind];
	  new_comp_edge_inds.push_back(cur_edge_ind);
	  
	  vector<hash_key> map_hks_to_update = map_hks[right_map_comp_ind];
	  
	  for(int j=0; j<map_hks[right_map_comp_ind].size(); j++){ //map hash keys
	    new_comp_map_hks.push_back(map_hks[right_map_comp_ind][j]);
	  }
	  for(int j=0; j<edge_comp_inds[right_map_comp_ind].size(); j++){ //edge inds
	    new_comp_edge_inds.push_back(edge_comp_inds[right_map_comp_ind][j]);
	  }	  
	  if(left_map_list_orient == cur_left_orient){ //map orients
	    if(right_map_list_orient == cur_right_orient){
	      for(int j=0; j<map_orients[right_map_comp_ind].size(); j++)
		new_comp_map_orients.push_back(map_orients[right_map_comp_ind][j]);	     
	    }
	    else{
	      for(int j=0; j<map_orients[right_map_comp_ind].size(); j++){
		new_comp_map_orients.push_back(rev(map_orients[right_map_comp_ind][j]));
	      }
	    }
	  }
	  else{	    	    
	    if(right_map_list_orient == rev(cur_right_orient)){
	      for(int j=0; j<map_orients[right_map_comp_ind].size(); j++)
		new_comp_map_orients.push_back(map_orients[right_map_comp_ind][j]);
	    }
	    else{
	      for(int j=0; j<map_orients[right_map_comp_ind].size(); j++){
		new_comp_map_orients.push_back(rev(map_orients[right_map_comp_ind][j]));
	      }
	    }
	  }
	  //make a new component end
	  
	  //erase two old components begin
	  vector<vector<hash_key> >::iterator 
	    min_hks_it = map_hks.begin()+int_min(left_map_comp_ind, right_map_comp_ind);
	  vector<vector<hash_key> >::iterator 
	    max_hks_it = map_hks.begin()+int_max(left_map_comp_ind, right_map_comp_ind);	  

	  vector<vector<orientation> >::iterator min_or_it = 
	    map_orients.begin() + int_min(left_map_comp_ind, right_map_comp_ind);
	  vector<vector<orientation> >::iterator max_or_it = 
	    map_orients.begin() + int_max(left_map_comp_ind, right_map_comp_ind);

	  vector<vector<int> >::iterator min_e_it = 
	    edge_comp_inds.begin() +  int_min(left_map_comp_ind, right_map_comp_ind);
	  vector<vector<int> >::iterator max_e_it = 
	    edge_comp_inds.begin() + int_max(left_map_comp_ind, right_map_comp_ind);

	  vector<hash_key>::iterator min_comp_hk_it = 
	    comp_hks.begin() + int_min(left_map_comp_ind, right_map_comp_ind);
	  vector<hash_key>::iterator max_comp_hk_it = 
	    comp_hks.begin() + int_max(left_map_comp_ind, right_map_comp_ind);

	  map_hks.erase(max_hks_it);
	  map_hks.erase(min_hks_it);	  
	  map_orients.erase(max_or_it);
	  map_orients.erase(min_or_it);
	  edge_comp_inds.erase(max_e_it);
	  edge_comp_inds.erase(min_e_it);
	  comp_hks.erase(max_comp_hk_it);
	  comp_hks.erase(min_comp_hk_it);
	  //erase two old components end

	  //add new component begin
	  map_hks.push_back(new_comp_map_hks);	  
	  map_orients.push_back(new_comp_map_orients);
	  edge_comp_inds.push_back(new_comp_edge_inds);
	  comp_hks.push_back(new_comp_hk);
	  //add new component end

	  //maintain hash tables begin
	  for(int j=0; j<map_hks_to_update.size(); j++) map_to_comp_ht[map_hks_to_update[j]] = new_comp_hk;
	  for(int j=0; j<comp_hks.size(); j++){
	    hash_key cur_comp_hash_key = comp_hks[j];
	    comp_to_ind_ht[cur_comp_hash_key] = j;
	  } //potentially can improve by updating only affected ones.
	  //maintain hash tables end

	  new_edge_inds.push_back(cur_edge_ind); //save a good edge
	}
      }
    }
  }
  return new_edge_inds;
}

vector<int> graph_component::terminal_nodes(){
  vector<int> term_nodes;

  for(int i=0; i<nodes.size(); i++){
    if(nodes[i].left_maps.size() == 0 &&
       nodes[i].container_maps.size()==0){
      term_nodes.push_back(i);
    }
  }

  return term_nodes;
}

int graph_component::leftmost_al_site(int ind){
  assert(ind >=0 && ind < nodes.size());
  int leftmost_site = -1;
  vector<int> adj_nodes;

  for(int i=0; i<nodes[ind].right_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].right_maps_inds[i]);
  for(int i=0; i<nodes[ind].left_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].left_maps_inds[i]);
  for(int i=0; i<nodes[ind].contained_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].contained_maps_inds[i]);
  for(int i=0; i<nodes[ind].container_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].container_maps_inds[i]);

  for(int i=0; i<adj_nodes.size(); i++){
    int cur_left_site = -1;
    //int cur_edge_ind = edge_ind(ind,adj_nodes[i]);
    int cur_edge_ind = nodes[ind].get_edge_ind(adj_nodes[i]);
    assert(cur_edge_ind != UNDEF_IND);
    int sites = edges[cur_edge_ind].left_al_sites.size();
    if(edges[cur_edge_ind].left_map == nodes[ind].read_name){
      if(edges[cur_edge_ind].left_orient == nodes[ind].orient)
	cur_left_site = edges[cur_edge_ind].left_al_sites[0];
      else cur_left_site = nodes[ind].map_read.size() - 
	     edges[cur_edge_ind].left_al_sites[sites-1];
    }
    else{
      assert(edges[cur_edge_ind].right_map == nodes[ind].read_name);
      if(edges[cur_edge_ind].right_orient != nodes[ind].orient)
	cur_left_site = nodes[ind].map_read.size() - 
	  edges[cur_edge_ind].right_al_sites[sites-1];
      else cur_left_site = edges[cur_edge_ind].right_al_sites[0];
    }

    if(i==0) leftmost_site = cur_left_site;
    else{
      if(cur_left_site < leftmost_site) leftmost_site = cur_left_site;
    }
  }
  assert(leftmost_site>=0 && leftmost_site <=nodes[ind].map_read.size());
  return leftmost_site;
}
int graph_component::rightmost_al_site(int ind){
  assert(ind >=0 && ind < nodes.size());
  int rightmost_site = -1;
  double dist = -1;
  vector<int> adj_nodes;

  for(int i=0; i<nodes[ind].right_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].right_maps_inds[i]);
  for(int i=0; i<nodes[ind].left_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].left_maps_inds[i]);
  for(int i=0; i<nodes[ind].contained_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].contained_maps_inds[i]);
  for(int i=0; i<nodes[ind].container_maps_inds.size(); i++)
    adj_nodes.push_back(nodes[ind].container_maps_inds[i]);

  for(int i=0; i<adj_nodes.size(); i++){
    int cur_right_site = -1;
    int cur_edge_ind = nodes[ind].get_edge_ind(adj_nodes[i]);
    //int cur_edge_ind = edge_ind(ind,adj_nodes[i]);
    assert(cur_edge_ind != UNDEF_IND);
    int sites = edges[cur_edge_ind].left_al_sites.size();
    if(edges[cur_edge_ind].left_map == nodes[ind].read_name){
      if(edges[cur_edge_ind].left_orient == nodes[ind].orient)
	cur_right_site = edges[cur_edge_ind].left_al_sites[sites-1];
      else cur_right_site = nodes[ind].map_read.size() - 
	     edges[cur_edge_ind].left_al_sites[0];
    }
    else{
      assert(edges[cur_edge_ind].right_map == nodes[ind].read_name);
      if(edges[cur_edge_ind].right_orient != nodes[ind].orient)
	cur_right_site = nodes[ind].map_read.size() - 
	  edges[cur_edge_ind].right_al_sites[0];
      else cur_right_site = edges[cur_edge_ind].right_al_sites[sites-1];
    }

    if(i==0) rightmost_site = cur_right_site;
    else{
      if(cur_right_site > rightmost_site) rightmost_site = cur_right_site;
    }
  }
  assert(rightmost_site>=0 && rightmost_site <=nodes[ind].map_read.size());
  return rightmost_site;
}
double graph_component::finalizing_dist(int ind){
  assert(ind >=0 && ind < nodes.size());
   
  int cur_right_site = rightmost_al_site(ind);
  
  double cur_dist = - nodes[ind].distance_to_center(cur_right_site,nodes[ind].orient);

  return cur_dist;
}
double graph_component::starting_dist(int ind){
  assert(ind >=0 && ind < nodes.size());
   
  int cur_left_site = leftmost_al_site(ind);
  
  double cur_dist = nodes[ind].distance_to_center(cur_left_site,nodes[ind].orient);

  return cur_dist;
}


vector<int> graph_component::extract_linear_path(){
  //path best_path;
  vector<int> best_path;


  //nodes in the component must be indexed. call index_nodes()
  vector<int> term_nodes = terminal_nodes();

  if(term_nodes.size()==0) return best_path; //needs to be 
  //addressed later, this may non-processed components due to cycles involving
  //potential terminal nodes

  vector<bool> node_discovered;
  vector<int> local_best_pred_inds;
  vector<int> node_colors; 
  //keeps a track of a current terminal node to avoid looping in short cycles

  for(int i=0; i<nodes.size(); i++){
    node_discovered.push_back(false);
    local_best_pred_inds.push_back(-2);
    node_colors.push_back(-1);
  }

  

  for(int i=0; i<term_nodes.size(); i++){
    int cur_term_node_ind = term_nodes[i];
    string cur_map_name = suffix(nodes[cur_term_node_ind].read_name, suffix_length);
    cout<<"dfs_path for terminal node: "<<cur_term_node_ind<<cur_map_name.c_str()<<endl;
    extract_linear_path_dfs
      (cur_term_node_ind , -1, 0, local_best_pred_inds, 
       node_discovered, node_colors, cur_term_node_ind);
  }


  int best_scoring_node_ind = -1;
  bool best_score_assigned = false;
  double top_score = 0;

  for(int i=0; i<nodes.size(); i++){
    int cur_node_ind = i;
    int best_pred_ind = local_best_pred_inds[cur_node_ind];

    string map1_id = suffix(nodes[i].read_name,suffix_length);
    string map2_id = "missing";
    if(best_pred_ind != -1 && best_pred_ind !=-2){
      cout<<"best_pred_ind: "<<best_pred_ind<<endl;
      //cout<<"here: "<<nodes[best_pred_ind].read_name.c_str()<<endl;
      assert(best_pred_ind < nodes.size());
      map2_id = suffix(nodes[best_pred_ind].read_name,suffix_length);
    }

    
    cout<<"for_node "<<i<<" ("<<map1_id.c_str()<<")";
    cout<<" best_pred: "<<best_pred_ind<<" ("<<map2_id.c_str()<<")"<<endl;
    if(node_discovered[cur_node_ind] && best_pred_ind != -1){
      double termination_score = 0;
      /*
      //int best_pred_ind = local_best_pred_inds[cur_node_ind];
      assert(best_pred_ind>=0 && best_pred_ind<nodes.size());
      int ovlp_ind = -1;
      bool ovlp_found = false;
      for(int j=0; j<nodes[best_pred_ind].right_maps_inds.size(); j++){
	if(cur_node_ind == nodes[best_pred_ind].right_maps_inds[j]){
	  ovlp_found = true;
	  ovlp_ind = nodes[best_pred_ind].right_edges_inds[j];
	}
      }
      for(int j=0; j<nodes[best_pred_ind].contained_maps_inds.size(); j++){
	if(cur_node_ind == nodes[best_pred_ind].contained_maps_inds[j]){
	  ovlp_found = true;
	  ovlp_ind = nodes[best_pred_ind].contained_edges_inds[j];
	}
      }
      assert(ovlp_ind>=0 && ovlp_ind<edges.size());
      assert(ovlp_found);
      
      //find the terminantion score
      
      int right_al_site=-1;
      int al_sites = edges[ovlp_ind].left_al_sites.size();
      assert(al_sites > 0);
      if(edges[ovlp_ind].right_map == nodes[cur_node_ind].read_name){
	if(edges[ovlp_ind].right_orient == nodes[cur_node_ind].orient)
	  right_al_site = edges[ovlp_ind].right_al_sites[al_sites-1];
	else right_al_site = nodes[cur_node_ind].map_read.size() - 
	       edges[ovlp_ind].right_al_sites[al_sites-1];
      }
      else{
	assert(edges[ovlp_ind].left_map == nodes[cur_node_ind].read_name);
	if(edges[ovlp_ind].left_orient != nodes[cur_node_ind].orient)
	  right_al_site = nodes[cur_node_ind].map_read.size() - 
	    edges[ovlp_ind].left_al_sites[0];
	else right_al_site = edges[ovlp_ind].left_al_sites[0];	  
      }
      assert(right_al_site >=0);
      assert(right_al_site <= nodes[cur_node_ind].map_read.size());
      //termination_score = 
      //- nodes[cur_node_ind].distance_to_center(right_al_site,nodes[cur_node_ind].orient);
      */

      termination_score = finalizing_dist(cur_node_ind);
      
      cout<<"term score: "<<termination_score<<" total: ";
      cout<<nodes[cur_node_ind].weight+termination_score<<endl;
      if(!best_score_assigned){
	top_score = nodes[cur_node_ind].weight+termination_score;
	best_score_assigned = true;
	best_scoring_node_ind = cur_node_ind;
      }
      else{
	if(nodes[cur_node_ind].weight+termination_score > top_score){
	  top_score = nodes[cur_node_ind].weight+termination_score;
	  best_scoring_node_ind = cur_node_ind;
	}
      }
    }
  }

  cout<<"best_terminal_node: "<<best_scoring_node_ind<<" total dist: "<<top_score<<endl;

  if(best_scoring_node_ind != -1){
    vector<int> path_inds;
    int sanity_check_counter = 0;
    bool extension_complete = false;
    
    int cur_node_ind = best_scoring_node_ind;
    
    while(!extension_complete){
      cout<<cur_node_ind<<endl;
      assert(cur_node_ind >= 0);
      assert(cur_node_ind < nodes.size());
      path_inds.push_back(cur_node_ind);
      
      if(local_best_pred_inds[cur_node_ind] != -1 && local_best_pred_inds[cur_node_ind] != -2){
	cur_node_ind = local_best_pred_inds[cur_node_ind];
      }
      else{
	extension_complete = true;
      }
      sanity_check_counter++;
      assert(sanity_check_counter < 200);
    }
    
    //cout<<"best_terminal_node: "<<best_scoring_node_ind<<" total dist: "<<top_score<<endl;
    cout<<"path: ";
    
    for(int i=path_inds.size()-1; i>=0; i--){
      cout<<path_inds[i]<<" ";
      best_path.push_back(path_inds[i]);
    }
    cout<<endl;
  }

  return best_path;
}

void graph_component::extract_linear_path_dfs
(int cur_node_ind, int prev_node_ind, double cur_score, vector<int>& local_best_pred_inds, 
 vector<bool>& node_discovered, vector<int>& node_colors, int master_color){
  assert(cur_node_ind >= 0 && cur_node_ind<nodes.size());
  
  if(node_colors[cur_node_ind] != master_color){
    //node was not reached from the same master before

    cout<<"   dfs_path: processing node: "<<suffix(nodes[cur_node_ind].read_name,suffix_length);
    cout<<" suggested score: "<<cur_score<<endl;

    //if(node_discovered[cur_node_ind]){
      //node discovered, maximize
    {
      if(node_discovered[cur_node_ind]){
	cout<<"      node discovered before. current weight: "<<cur_score;
	cout<<" prev. weight: "<<nodes[cur_node_ind].weight<<endl;
      }
      else{
	cout<<"      node never discovered before"<<endl;
      }
      //mark colors, assign discovered
      node_colors[cur_node_ind] = master_color;
      
      if(nodes[cur_node_ind].weight <= cur_score || !node_discovered[cur_node_ind]){	
	//better score, store a better score
	if(node_discovered[cur_node_ind]){
	  if(nodes[cur_node_ind].weight < cur_score){
	    nodes[cur_node_ind].weight = cur_score; //set a better score
	    local_best_pred_inds[cur_node_ind] = prev_node_ind;
	  }
	}
	else{
	  nodes[cur_node_ind].weight = cur_score;
	  local_best_pred_inds[cur_node_ind] = prev_node_ind;
	  node_discovered[cur_node_ind] = true;
	}
	
	vector<int> succs; //compile a list of successors
	vector<double> succ_scores;
	for(int i=0; i<nodes[cur_node_ind].right_maps_inds.size(); i++){
	  succs.push_back(nodes[cur_node_ind].right_maps_inds[i]);
	  if(cur_node_ind != master_color){ //current node is not the first
	    succ_scores.push_back(nodes[cur_node_ind].right_edge_weights[i]);
	    cout<<"not a master node"<<endl;
	  }
	  else{ 
	    //current node is the first (leftmost), 
	    //accomodate the overlap length
	    cout<<"master node!"<<endl;
	    int cur_ovlp_ind = nodes[cur_node_ind].right_edges_inds[i];
	    string term_map_id = nodes[cur_node_ind].read_name;
	    int term_map_left_al_site=-1;
	    if(edges[cur_ovlp_ind].left_map == term_map_id){
	      if(edges[cur_ovlp_ind].left_orient == nodes[cur_node_ind].orient)
		term_map_left_al_site = edges[cur_ovlp_ind].left_al_sites[0];
	      else term_map_left_al_site = nodes[cur_node_ind].map_read.size()-
		     edges[cur_ovlp_ind].left_al_sites[0];
	    }
	    else{
	      assert(edges[cur_ovlp_ind].right_map == term_map_id);
	      int sites = edges[cur_ovlp_ind].left_al_sites.size();
	      if(edges[cur_ovlp_ind].right_orient != nodes[cur_node_ind].orient)
		term_map_left_al_site = nodes[cur_node_ind].map_read.size()-
		  edges[cur_ovlp_ind].right_al_sites[sites-1];
	      else term_map_left_al_site = edges[cur_ovlp_ind].right_al_sites[0];
	    }
	    assert(term_map_left_al_site>=0);
	    assert(term_map_left_al_site<=nodes[cur_node_ind].map_read.size());
	    
	    //double cur_weight =nodes[cur_node_ind].right_edge_weights[i] + 
	    //nodes[cur_node_ind].distance_to_center
	    //(term_map_left_al_site,nodes[cur_node_ind].orient);	      
	    double cur_weight = nodes[cur_node_ind].right_edge_weights[i] + starting_dist(cur_node_ind);
	    succ_scores.push_back(cur_weight);
	    
	    cout<<"ncont: master: "<<master_color<<" cur_node: "<<cur_node_ind;
	    cout<<" ovlp: "<<cur_ovlp_ind;
	    cout<<" term_map_site: "<<term_map_left_al_site;
	    cout<<" total_dist: "<<cur_weight<<endl;
	  }
	}
	for(int i=0; i<nodes[cur_node_ind].contained_maps_inds.size(); i++){
	  succs.push_back(nodes[cur_node_ind].contained_maps_inds[i]);
	  if(cur_node_ind != master_color) //current node is not the first
	    succ_scores.push_back(nodes[cur_node_ind].contained_edge_weights[i]);
	  else{ //current node is the first (leftmost), accomodate the overlap length
	    int cur_ovlp_ind = nodes[cur_node_ind].contained_edges_inds[i];
	    string term_map_id = nodes[cur_node_ind].read_name;
	    int term_map_left_al_site;
	    int sites = edges[cur_ovlp_ind].left_al_sites.size();	      
	    if(edges[cur_ovlp_ind].left_map == term_map_id){
	      if(edges[cur_ovlp_ind].left_orient == nodes[cur_node_ind].orient){
		term_map_left_al_site = edges[cur_ovlp_ind].left_al_sites[0]; cout<<"1"<<endl;}
	      else{ term_map_left_al_site = nodes[cur_node_ind].map_read.size()-
		      edges[cur_ovlp_ind].left_al_sites[sites-1]; cout<<"2"<<endl;}
	    }
	    else{
	      assert(edges[cur_ovlp_ind].right_map == term_map_id);
	      if(edges[cur_ovlp_ind].right_orient != nodes[cur_node_ind].orient){
		term_map_left_al_site = nodes[cur_node_ind].map_read.size()-
		  edges[cur_ovlp_ind].right_al_sites[sites-1]; cout<<"3"<<endl;}
	      else{ term_map_left_al_site = edges[cur_ovlp_ind].right_al_sites[0]; cout<<"4"<<endl;}
	    }
	    assert(term_map_left_al_site>=0);
	    assert(term_map_left_al_site<=nodes[cur_node_ind].map_read.size());
	    double left_dist = nodes[cur_node_ind].distance_to_center
	      (term_map_left_al_site,nodes[cur_node_ind].orient);	

	    double cur_weight =nodes[cur_node_ind].contained_edge_weights[i] + 
	      left_dist;      
	    succ_scores.push_back(cur_weight);
	    
	    cout<<"cont: master: "<<master_color<<" cur_node: "<<cur_node_ind;
	    cout<<" ovlp: "<<cur_ovlp_ind;
	    cout<<" term_map_site: "<<term_map_left_al_site;
	    cout<<" edge_w: "<<nodes[cur_node_ind].contained_edge_weights[i];
	    cout<<" total_dist: "<<cur_weight;
	    cout<<" left_dist: "<<left_dist<<endl;
	  }
	}

	for(int i=0; i<succs.size(); i++){
	  int cur_succ = succs[i];
	  double cur_succ_score = nodes[cur_node_ind].weight + succ_scores[i];
	  
	  cout<<"       succ: "<<suffix(nodes[cur_succ].read_name,suffix_length)<<" dist: "<<cur_succ_score<<endl;
	}
	for(int i=0; i<succs.size(); i++){
	  int cur_succ = succs[i];
	  double cur_succ_score = nodes[cur_node_ind].weight + succ_scores[i];	
	  extract_linear_path_dfs(cur_succ, cur_node_ind, cur_succ_score, local_best_pred_inds,
				  node_discovered, node_colors, master_color);
	}  
	return;
      }
      else{
	//cout<<"node not discovered"<<endl;
	//previous score better, do nothing
	return;
      }
    }
    /*
    else{
      //node never discovered,
      //either extend or set this to be a terminal node
      
      cout<<"         node never discovered."<<endl;

      //mark colors, assign discovered
      node_colors[cur_node_ind] = master_color;
      node_discovered[cur_node_ind] = true;

      double terminal_score = -10000000;//nodes[cur_node_ind].size()/2;
      if(cur_score > terminal_score){
	//extend
	nodes[cur_node_ind].weight = cur_score; //set an extension score
	local_best_pred_inds[cur_node_ind] = prev_node_ind;	
      }
      else{
	//set terminal node
	nodes[cur_node_ind].weight = terminal_score;
	local_best_pred_inds[cur_node_ind] = -1;//terminal
      }

      vector<int> succs; //compile a list of successors
      vector<double> succ_scores;
      for(int i=0; i<nodes[cur_node_ind].right_maps_inds.size(); i++){
	succs.push_back(nodes[cur_node_ind].right_maps_inds[i]);
	succ_scores.push_back(nodes[cur_node_ind].right_edge_weights[i]);
      }
      for(int i=0; i<nodes[cur_node_ind].contained_maps_inds.size(); i++){
	succs.push_back(nodes[cur_node_ind].contained_maps_inds[i]);
	succ_scores.push_back(nodes[cur_node_ind].contained_edge_weights[i]);
      }
      
      for(int i=0; i<succs.size(); i++){
	int cur_succ = succs[i];
	double cur_succ_score = nodes[cur_node_ind].weight + succ_scores[i];
	cout<<"       succ: "<<suffix(nodes[cur_succ].read_name,suffix_length)<<" dist: "<<cur_succ_score<<endl;
		  
	extract_linear_path_dfs(cur_succ, cur_node_ind, cur_succ_score, local_best_pred_inds,
				node_discovered, node_colors, master_color);
      }
      return;
    }
*/
  }
  else{
    //node was reached from the same master before, do nothing
    return;
  }
}

list<int> graph_component::
confirmed_edges(int dfs_depth){

  list<int> confirmed_edge_set;
  //this is the list of edges that are confirmed

  vector<bool> marked_edges;
  for(int i=0; i<edges.size(); i++){
    marked_edges.push_back(false);
  }

  vector<int> loc_pred_list;
  vector<int> loc_node_color;
  vector<int> branch_root_list;

  for(int i=0; i<nodes.size(); i++){
    loc_node_color.push_back(-1);
    loc_pred_list.push_back(-1);
    branch_root_list.push_back(-1);
  }
  
  vector<short> edge_mask;
  for(int i=0; i<edges.size(); i++){
    edge_mask.push_back(0);
  }

  for(int i=0; i<nodes.size(); i++){
    int cur_node_mark_color = i;
    int cur_node_ind = i;

    vector<int> node_inds; 
    //nodes to which there are multiple paths from the 
    //current node

    vector<vector<path> > paths;
    //set of multiple paths going to each of the nodes
    cout<<endl<<"master node:"<<cur_node_ind<<endl;
    possible_paths_for_consistency_check
      (node_inds, paths, loc_node_color, loc_pred_list, branch_root_list,
       cur_node_mark_color, dfs_depth, cur_node_ind,
       -1, cur_node_ind, -1);
    
    assert(node_inds.size() == paths.size());
    for(int j=0; j<node_inds.size(); j++){
      //no common region

      vector<int> best_path_set = good_paths(paths[j]);
      if(best_path_set.size()>1){
	//update the list of confirmed edges with the edges from this
	//path set

	//this is not very efficient search, will improve later
	for(int k=0; k<best_path_set.size(); k++){
	  int cur_path_ind = best_path_set[k];
	  //current path is paths[j][cur_path_ind];
	  for(int m=1; m<paths[j][cur_path_ind].map_inds.size(); m++){
	    int map1_ind = paths[j][cur_path_ind].map_inds[m-1];
	    int map2_ind = paths[j][cur_path_ind].map_inds[m];

	    int cur_edge_ind = nodes[map1_ind].get_edge_ind(map2_ind);
	    marked_edges[cur_edge_ind] = true;
	    //int cur_edge_ind = edge_ind(map1_ind, map2_ind);
	    /*
	    assert(cur_edge_ind != UNDEF_IND);
	    marked_edges[cur_edge_ind] = true;
	    bool edge_included = false;
	    for(list<int>::iterator it = confirmed_edge_set.begin();
		it != confirmed_edge_set.end(); it++){
	      //for(int l=0; l<confirmed_edge_set.size(); l++){
	      if(*it == cur_edge_ind) edge_included = true;
	    }
	    if(!edge_included){	      
	      confirmed_edge_set.push_back(cur_edge_ind);
	    }
	    */
	  }
	}
      }
    }
  }
  //for(vector<bool>::iterator it = maked_edges.begin(); it != marked_edges.end(); it++){
  for(int i=0; i<marked_edges.size(); i++){
    if(marked_edges[i]) confirmed_edge_set.push_back(i);
  }
  
  /*
  cout<<"+++++++++++++++++++++++++++++++++"<<endl;
  cout<<"edges for current component:"<<endl;
  cout<<"originally:"<<edges.size()<<" confirmed:"<<confirmed_edge_set.size()<<endl;
  for(list<int>::iterator it = confirmed_edge_set.begin();
      it != confirmed_edge_set.end(); it++){
    //for(int i=0; i<confirmed_edge_set.size(); i++){
    cout<<*it<<endl;
  }
  */

  return confirmed_edge_set;
}

void graph_component::clean(){
  for(int i=0; i<nodes.size(); i++){
    nodes[i].right_maps_inds.clear();
    nodes[i].right_edges_inds.clear();
    nodes[i].right_edge_weights.clear();
    
    nodes[i].left_maps_inds.clear();
    nodes[i].left_edges_inds.clear();
    nodes[i].left_edge_weights.clear();

    nodes[i].contained_edges_inds.clear();
    nodes[i].contained_maps_inds.clear();
    nodes[i].contained_edge_weights.clear();

    nodes[i].container_maps_inds.clear();
    nodes[i].container_edges_inds.clear();
    nodes[i].container_edge_weights.clear();
  }
}

void graph_component::output_component(const char* output_file){
  cout<<"ouputting component as a graph."<<endl;

  remove (output_file);
  ofstream out_str;
  out_str.open(output_file);

  if(!out_str.good()){
    cout<<"bad graph output file name: "<<output_file;
    cout<<endl;
    assert(false);
  }

  out_str<<"digraph G{"<<endl;
  out_str<<"size="<<char(34)<<"8,5"<<char(34)<<";"<<endl;
  out_str<<"rankdir=LR;"<<endl;
  //out_str<<"node [shape=circle];"<<endl;

  //index_nodes();
  //assign_edge_weights();
  ostringstream prefix;
  output_edges(out_str,prefix);
  /*
  for(int i=0; i<components.size(); i++){
    ostringstream comp_prefix;
    comp_prefix<<i<<":";
    //string prefix = comp_prefix.c_str();
    components[i].index_nodes();
    components[i].output_edges(out_str, comp_prefix);
  }
  */

  out_str<<"}"<<endl;
  out_str.close();
}

void graph_component::output_edges(ofstream& out_str, ostringstream& prefix){
  int counter=0;
  for(int i=0; i<nodes.size(); i++){
    for(int j=0; j<nodes[i].right_maps_inds.size(); j++){
      
      counter++;
      out_str<<char(34)<<prefix.str();
      out_str<<nodes[i].id;
      //<<nodes[i]._index<<":";
      //out_str<<nodes[i].orient;
      //out_str<<":"<<nodes[i].read_name; //to kill
      out_str<<char(34);
      out_str<<" -> ";
      out_str<<char(34)<<prefix.str();
      out_str<<nodes[nodes[i].right_maps_inds[j]].id;
      //out_str<<nodes[i].right_maps_inds[j]<<":";
      //out_str<<nodes[nodes[i].right_maps_inds[j]].orient;
      //out_str<<((int)nodes[nodes[i].right_maps_inds[j]].weight)<<char(34);
      //out_str<<":"<<nodes[nodes[i].right_maps_inds[j]].read_name; //to kill
      out_str<<char(34); 
      out_str<<" [color="<<char(34)<<"0.002 0.999 0.999"<<char(34)<<"]";
      out_str<<" [label="<<char(34);
      //out_str<<"e:"<<prefix.str()<<nodes[i].right_edges_inds[j]<<":";
      out_str<<((int)nodes[i].right_edge_weights[j]);
      out_str<<":";
      out_str<<edges[nodes[i].right_edges_inds[j]].t_score;
      out_str<<char(34)<<"]";
      out_str<<";"<<endl;
    }
    for(int j=0; j<nodes[i].contained_maps_inds.size(); j++){
      counter++;
      out_str<<char(34)<<prefix.str();
      out_str<<nodes[i].id;
      //out_str<<nodes[i]._index<<":";
      //out_str<<nodes[i].orient;
      //out_str<<":"<<nodes[i].read_name; //to kill
      out_str<<char(34);
      out_str<<" -> ";
      out_str<<char(34)<<prefix.str();
      out_str<<nodes[nodes[i].contained_maps_inds[j]].id;
      //out_str<<nodes[i].contained_maps_inds[j]<<":";
      //out_str<<nodes[nodes[i].contained_maps_inds[j]].orient;
      //out_str<<":"<<nodes[nodes[i].contained_maps_inds[j]].read_name; // to kill
      out_str<<char(34);
      //out_str<<((int)nodes[nodes[i].contained_maps_inds[j]].weight)<<char(34);
      out_str<<" [color="<<char(34)<<"0.650 0.700 0.700"<<char(34)<<"]";
      
      out_str<<" [label="<<char(34);
      
      //out_str<<"e:"<<prefix.str()<<nodes[i].contained_edges_inds[j]<<":";
      out_str<<((int)nodes[i].contained_edge_weights[j]);
      out_str<<":";
      out_str<<edges[nodes[i].contained_edges_inds[j]].t_score;
      out_str<<char(34)<<"]";
      out_str<<";"<<endl;
    }
  }

  cout<<"stored "<<counter<<" edges "<<endl;
  //for(int i=0; i<edges.size(); i++){
  //  out_str<<edges[i].left_map<<" -> "<<edges[i].right_map<<";"<<endl;
  //}
}


void graph_component::output_graph(string fname){
  int suffix_size = suffix_length; //this many letters from the end of node_name

  remove(fname.c_str());
  ofstream graph_str;
  graph_str.open(fname.c_str());
  assert(graph_str.good());

  graph_str<<"digraph G{"<<endl;
  graph_str<<"size="<<char(34)<<"8,5"<<char(34)<<";"<<endl;
  graph_str<<"rankdir=LR;"<<endl;

  int counter=0;
  for(int i=0; i<nodes.size(); i++){
    string map1_id = suffix(nodes[i].read_name,suffix_size);

    for(int j=0; j<nodes[i].right_maps_inds.size(); j++){
      counter++;
      
      string map2_id = suffix
	(nodes[nodes[i].right_maps_inds[j]].read_name, suffix_size);           
      graph_str<<char(34)<<map1_id.c_str();

      //out_str<<char(34)<<prefix.str();
      //out_str<<nodes[i].id;
      //<<nodes[i]._index<<":";
      //out_str<<nodes[i].orient;
      //out_str<<":"<<nodes[i].read_name; //to kill
      graph_str<<char(34);
      graph_str<<" -> ";
      
      graph_str<<char(34)<<map2_id.c_str();
      //out_str<<char(34)<<prefix.str();
      //out_str<<nodes[nodes[i].right_maps_inds[j]].id;
      //out_str<<nodes[i].right_maps_inds[j]<<":";
      //out_str<<nodes[nodes[i].right_maps_inds[j]].orient;
      //out_str<<((int)nodes[nodes[i].right_maps_inds[j]].weight)<<char(34);
      //out_str<<":"<<nodes[nodes[i].right_maps_inds[j]].read_name; //to kill
      graph_str<<char(34); 
      graph_str<<" [color="<<char(34)<<"0.002 0.999 0.999"<<char(34)<<"]";
      graph_str<<" [label="<<char(34);
      //out_str<<"e:"<<prefix.str()<<nodes[i].right_edges_inds[j]<<":";
      graph_str<<((int)nodes[i].right_edge_weights[j]);
      graph_str<<":";
      graph_str<<edges[nodes[i].right_edges_inds[j]].t_score;
      graph_str<<char(34)<<"]";
      graph_str<<";"<<endl;
    }
    for(int j=0; j<nodes[i].contained_maps_inds.size(); j++){
      counter++;

      string map2_id = suffix
	(nodes[nodes[i].contained_maps_inds[j]].read_name, suffix_size);           
      graph_str<<char(34)<<map1_id.c_str();
      

      //out_str<<char(34)<<prefix.str();
      //out_str<<nodes[i].id;
      //out_str<<nodes[i]._index<<":";
      //out_str<<nodes[i].orient;
      //out_str<<":"<<nodes[i].read_name; //to kill
      graph_str<<char(34);
      graph_str<<" -> ";

      graph_str<<char(34)<<map2_id.c_str();

      //out_str<<char(34)<<prefix.str();
      //out_str<<nodes[nodes[i].contained_maps_inds[j]].id;
      //out_str<<nodes[i].contained_maps_inds[j]<<":";
      //out_str<<nodes[nodes[i].contained_maps_inds[j]].orient;
      //out_str<<":"<<nodes[nodes[i].contained_maps_inds[j]].read_name; // to kill
      graph_str<<char(34);
      //out_str<<((int)nodes[nodes[i].contained_maps_inds[j]].weight)<<char(34);
      graph_str<<" [color="<<char(34)<<"0.650 0.700 0.700"<<char(34)<<"]";
      
      graph_str<<" [label="<<char(34);
      
      //out_str<<"e:"<<prefix.str()<<nodes[i].contained_edges_inds[j]<<":";
      graph_str<<((int)nodes[i].contained_edge_weights[j]);
      graph_str<<":";
      graph_str<<edges[nodes[i].contained_edges_inds[j]].t_score;
      graph_str<<char(34)<<"]";
      graph_str<<";"<<endl;
    }
  }

  cout<<"stored "<<counter<<" edges "<<endl;

  graph_str<<"}"<<endl;
  graph_str.close();
}

/*
void graph_component::add_edge(edge& e, vector<fr_size>& left_map_read,
			       vector<fr_size>& right_map_read){
  if(edges.empty()){
    assert(nodes.empty());

    
    node new_node1(e.left_map_hash_key, e.left_map, e.left_orient, left_map_read);
    node new_node2(e.right_map_hash_key, e.right_map, e.right_orient, right_map_read);

    if(!e.containment){
      new_node1.right_maps.push_back(e.right_map);
      new_node2.left_maps.push_back(e.left_map);
    }
    else{
      new_node1.contained_maps.push_back(e.right_map);
      new_node2.container_maps.push_back(e.left_map);
    }
    
    nodes.push_back(new_node1);
    nodes.push_back(new_node2);

    edges.push_back(e);
    return;
  }
  else{
    //find the map already included in 
    //the component
    bool left_map_included = false;
    bool right_map_included = false;

    int left_map_ind = -1;
    int right_map_ind = -1;
    for(int i=0; i<nodes.size(); i++){
      if(e.left_map == nodes[i].read_name){
	left_map_included = true;
	left_map_ind = i;
      }
      if(e.right_map == nodes[i].read_name){
	right_map_included = true;
	right_map_ind = i;
      }
    }
    
    if(!left_map_included){
      assert(right_map_included);
      //right map is included in the component,
      
      //add a new node
      node new_node(e.left_map_hash_key, e.left_map, e.left_orient, left_map_read);
      if(!e.containment)
	new_node.right_maps.push_back(e.right_map);
      else
	new_node.contained_maps.push_back(e.right_map);

      assert(right_map_ind >=0 && right_map_ind < nodes.size());

      //check the orientation: connected map must be 
      //in the consistent orientation
      if(e.right_orient == nodes[right_map_ind].orient){
	cout<<"including left map: using a given overlap"<<endl;
	if(!e.containment) nodes[right_map_ind].left_maps.push_back(e.left_map);
	else nodes[right_map_ind].container_maps.push_back(e.left_map);
	nodes.push_back(new_node);
      }
      else{
	cout<<"including left map: using inverse overlap"<<endl;
	if(!e.containment) nodes[right_map_ind].right_maps.push_back(e.left_map);
	else nodes[right_map_ind].container_maps.push_back(e.left_map);
	nodes.push_back(new_node.inverse());
	//the map in the component is in another orientation,
	//so need to flip overlap and adjust appropriately
      }
    }
    else{
      //the left map is included, the right may or may not be included
      if(!right_map_included){
	//make new node for the right_map

	node new_node(e.right_map_hash_key, e.right_map, e.right_orient, right_map_read);
	if(!e.containment) new_node.left_maps.push_back(e.left_map);
	else new_node.container_maps.push_back(e.left_map);

	assert(left_map_ind >=0 && left_map_ind < nodes.size());

	if(e.left_orient == nodes[left_map_ind].orient){
	  cout<<"including right map: using a given overlap"<<endl;
	  if(!e.containment) nodes[left_map_ind].right_maps.push_back(e.right_map);
	  else nodes[left_map_ind].contained_maps.push_back(e.right_map);
	  nodes.push_back(new_node);
	}
	else{
	  cout<<"including right map: using inverse overlap"<<endl;
	  if(!e.containment) nodes[left_map_ind].left_maps.push_back(e.right_map);
	  else nodes[left_map_ind].contained_maps.push_back(e.right_map);
	  nodes.push_back(new_node.inverse());
	}
      }
      else{
	//both nodes exist, just update the the lists of left/right maps
	assert(left_map_ind >= 0 && left_map_ind < nodes.size());
	assert(right_map_ind >= 0 && right_map_ind < nodes.size());

	cout<<"both nodes exist, updating the neighbor list"<<endl;

	bool consistent_orientation = 
	  (e.left_orient == nodes[left_map_ind].orient &&
	   e.right_orient == nodes[right_map_ind].orient) ||
	  (rev(e.left_orient) == nodes[left_map_ind].orient &&
	   rev(e.right_orient) == nodes[right_map_ind].orient);

	if(!consistent_orientation){
	  output_component("./test.dot");	 

	  cout<<"inconsistent orientation of maps:"<<endl;
	  cout<<"component contains an error"<<endl;
	  cout<<"outputting exact path:"<<endl;

	  bool path_exists = all_way_path(e.left_map, e.right_map);
	  if(path_exists){
	    cout<<"path exists "<<endl;
	    for(int i=0; i<helper1.size(); i++){
	      if(i!=0) cout<<"->";
	      cout<<helper1[i];
	    }
	    cout<<endl;
	    for(int i=0; i<helper1.size(); i++){
	      if(i!=0) cout<<"->";
	      cout<<helper1[i];
	      cout<<" ("<<nodes[helper1[i]].read_name<<") ";	  
	    }
	    cout<<endl;
	  
	    //assert(false);	  
	  }
	  assert(false);
	}

	if(e.left_orient == nodes[left_map_ind].orient){
	  if(!e.containment) nodes[left_map_ind].right_maps.push_back(e.right_map);
	  else nodes[left_map_ind].contained_maps.push_back(e.right_map);
	}
	else{	  
	  if(!e.containment) nodes[left_map_ind].left_maps.push_back(e.right_map);
	  else nodes[left_map_ind].contained_maps.push_back(e.right_map);
	}
	
	if(e.right_orient == nodes[right_map_ind].orient){
	  if(!e.containment) nodes[right_map_ind].left_maps.push_back(e.left_map);
	  else nodes[right_map_ind].container_maps.push_back(e.left_map);
	}
	else{
	  if(!e.containment) nodes[right_map_ind].right_maps.push_back(e.left_map);
	  else nodes[right_map_ind].container_maps.push_back(e.left_map);
	}
      }
    }
    edges.push_back(e);
    return;
  }
}
*/

void graph_component::print(){
  cout<<"printing component:"<<endl;
  cout<<"nodes: "<<nodes.size();
  cout<<" edges: "<<edges.size();

  for(int i=0; i<nodes.size(); i++){
    cout<<endl;
    nodes[i].print();
  }

  for(int i=0; i<edges.size(); i++){
    cout<<endl;
    //edges[i].print();
    //int left_map_ind = get_node_index(edges[i].left_map);
    //int right_map_ind = get_node_index(edges[i].right_map);

    //assert(left_map_ind != -1);
    //assert(right_map_ind != -1);

    hash_key left_map_hk = edges[i].left_map_hash_key;
    hash_key right_map_hk = edges[i].right_map_hash_key;

    cout<<"overlap: "<<left_map_hk<<" => "<<right_map_hk;
    cout<<" ("<<return_map_name(left_map_hk)<<" => ";
    cout<<return_map_name(right_map_hk)<<")";
    cout<<" lo:"<<edges[i].left_orient<<" ro:"<<edges[i].right_orient;
    //cout<<" "<<nodes[left_map_ind].map_read.size();
    //cout<<" "<<nodes[right_map_ind].map_read.size();;
    cout<<endl;
    
    //int left_center = nodes[left_map_ind].map_center(edges[i].left_orient);
    //int right_center = nodes[right_map_ind].map_center(edges[i].right_orient);
    //cout<<"lc:"<<left_center<<" rc:"<<right_center<<endl;
    
    int left_map_al_begin = edges[i].left_al_sites[0];
    int left_map_al_end = 
      edges[i].left_al_sites[edges[i].left_al_sites.size()-1];

    int right_map_al_begin = edges[i].right_al_sites[0];
    int right_map_al_end = 
      edges[i].right_al_sites[edges[i].right_al_sites.size()-1];

    /*
    for(int k=0; k<left_map_al_begin; k++){
      if(edges[i].left_orient == forward){
	cout<<"["<<k<<":"<<nodes[left_map_ind].map_read[k]<<"]"<<endl;
      }
      else{
	cout<<"["<<k<<":"<<nodes[left_map_ind].map_read
				 [nodes[left_map_ind].
				  map_read.size()-k-1]<<"]"<<endl;
      }
    }
    */
    /*
    for(int k=1; k<edges[i].left_al_sites.size(); k++){
      int cur_left_map_gap_start = edges[i].left_al_sites[k-1];
      int cur_left_map_gap_end = edges[i].left_al_sites[k];
      int cur_right_map_gap_start = edges[i].right_al_sites[k-1];
      int cur_right_map_gap_end = edges[i].right_al_sites[k];
      
      cout<<"[";
      for(int m=cur_left_map_gap_start; m<cur_left_map_gap_end; m++){
	if(m!=cur_left_map_gap_start)
	  cout<<", ";	
	if(edges[i].left_orient == forward)
	  cout<<m<<":"<<nodes[left_map_ind].map_read[m];
	else
	  cout<<m<<":"<<nodes[left_map_ind].
	    map_read[nodes[left_map_ind].map_read.size()-m-1];
      }
      cout<<"] -> [";
      for(int m=cur_right_map_gap_start; m<cur_right_map_gap_end; m++){
	if(m!=cur_right_map_gap_start) cout<<", ";
	if(edges[i].right_orient == forward)
	  cout<<m<<":"<<nodes[right_map_ind].map_read[m];
	else
	  cout<<m<<":"<<nodes[right_map_ind].
	    map_read[nodes[right_map_ind].map_read.size()-m-1];
      }
      cout<<"]"<<endl;
    }
    int left_map_leftover = nodes[left_map_ind].map_read.size()-left_map_al_end;
    int right_map_leftover = nodes[right_map_ind].map_read.size()-right_map_al_end;

    for(int k=0; k<int_max(left_map_leftover,right_map_leftover); k++){
      if(k<left_map_leftover){
	int cur_left_ind = k+left_map_al_end;
	if(edges[i].left_orient == forward){
	  cout<<"["<<cur_left_ind<<":"<<nodes[left_map_ind].map_read[cur_left_ind];
	  cout<<"] ";
	}
	else{
	  cout<<"["<<cur_left_ind<<":"<<nodes[left_map_ind].
	    map_read[nodes[left_map_ind].map_read.size()-cur_left_ind-1];
	  cout<<"] ";
	}
      }
      else{
	cout<<char(9)<<char(9);
      }
      if(k<right_map_leftover){
	int cur_right_ind = k+right_map_al_end;
	if(edges[i].right_orient == forward){
	  cout<<"["<<cur_right_ind<<":"<<nodes[right_map_ind].map_read[cur_right_ind];
	  cout<<"] ";
	}
	else{
	  cout<<"["<<cur_right_ind<<":"<<nodes[right_map_ind].
	    map_read[nodes[right_map_ind].map_read.size()-cur_right_ind-1];
	  cout<<"]";
	}
      }
      cout<<endl;
    }
    */
    /*
    for(int k=right_map_al_end; k<nodes[right_map_ind].map_read.size(); k++){
      cout<<char(9)<<"["<<k<<":";
      if(edges[i].right_orient == forward)
	cout<<nodes[right_map_ind].map_read[k];
      else
	cout<<nodes[right_map_ind].
	  map_read[nodes[right_map_ind].map_read.size()-k-1];
      cout<<"]"<<endl;
    }
    */
    cout<<endl;
  }
}

/*
int map_center(vector<double>& map_read){
  double total_size = 0;
  for(int i=0; i<map_read.size(); i++){
    total_size += map_read[i];
  }
  double center = total_size/2.0;
  
  bool center_found = false;

  if(map_read.size()<1) assert(false);

  int cur_pos = 0; //0-th fragment
  double cur_cum_size = map_read[cur_pos];
  
  while(!center_found){
    assert(cur_pos < map_read.size());
    if(cur_cum_size >= center) return cur_pos;
    cur_pos++;    
  }
  assert(center_found);
}
*/
/*
double distance_between_maps(vector<double>& map1_read, vector<double>& map2_read,		  		      orientation map1_orient, orientation map2_orient,
		      vector<int>& al_sites1, vector<int>& al_sites2){
  //the distance is given from map1 to map2, not the other way around

  double distance = 0;
  int dist_sign=1;

  int abs_center1 = map_center(map1_read);
  int abs_center2 = map_center(map2_read);

  int rel_center1, rel_center2;

  if(map1_orient
  
  assert(al_sites1.size()>1);
  assert(al_sites1.size() == al_sites2.size());

  int al_size = al_sites1.size();

  int map1_common_al_begin;
  int map1_common_al_end;
  int map2_common_al_begin;
  int map2_common_al_end;

  bool common_region_exist = false;


  int left_common_al_begin;
  int left_common_al_end;
  int right_common_al_begin;
  int right_common_al_end;
  
  bool common_region_exists = false;
  bool left_end_found = false;
  for(int i=0; i<edges[cur_edge].left_al_sites.size(); i++){
    int cur_left_site_ind = edges[cur_edge].left_al_sites[i];
    int cur_right_site_ind =  edges[cur_edge].right_al_sites[i];

    if((left_center < cur_left_site_ind &&
	right_center >= cur_right_site_ind) ||
       (left_center >= cur_left_site_ind &&
	right_center < cur_right_site_ind)){
      //there is aligned site between two centers

      if(!left_end_found){
	left_end_found = true;
	common_region_exists = true;
	
	left_common_al_begin = cur_left_site_ind;
	right_common_al_begin = cur_right_site_ind;
      }
      
      left_common_al_end = cur_left_site_ind;
      right_common_al_end = cur_right_site_ind;
    }
  }

 
  if(common_region_exists){

    if(left_center <= left_common_al_begin){
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
      else{
	if(same_orient_as_overlap){
	  dist_sign = -1;
	}
	else{
	  dist_sign = 1;
	}
      }
    }
    else{
      assert(left_center >= left_common_al_end);
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = -1;
	else
	  dist_sign = 1;
      }
      else{
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
    }
    //else{
    //assert(left_center >= common_al_end);
    //if(same_orient_as_overlap) dist_sign = -1;
    //else dist_sign = 1;
    //}

    cout<<"dist_sign:"<<dist_sign<<endl;

    cout<<"common al region: (";
    cout<<left_common_al_begin<<","<<left_common_al_end<<")->(";
    cout<<right_common_al_begin<<","<<right_common_al_end<<")"<<endl;

    assert(left_common_al_begin <= left_common_al_end);
    assert(right_common_al_begin <= right_common_al_end);

    cout<<"left map region:"<<left_common_al_begin;
    cout<<" to: "<<left_common_al_end<<endl;
    double common_left_map_region = 
      nodes[left_map_ind].region_size(left_common_al_begin, left_common_al_end, left_orient);

    cout<<"right map region:"<<right_common_al_begin;
    cout<<" to: "<<right_common_al_end<<endl;
    double common_right_map_region = 
      nodes[right_map_ind].region_size(right_common_al_begin, right_common_al_end, right_orient);
    
    cout<<"left common region size:"<<common_left_map_region<<endl;
    cout<<"right common region size:"<<common_right_map_region<<endl;
    distance += 0.5*(common_left_map_region+common_right_map_region);
    //add the distance of the common region

    //now add the tails
    
    //left_map:
    distance += 0.5*
      nodes[left_map_ind].region_size(left_center, left_center+1, left_orient);
    
    double left_tail;
    if(left_center < left_common_al_begin){
      left_tail = nodes[left_map_ind].
	region_size(left_center+1, left_common_al_begin, left_orient);
    }
    else{
      assert(left_common_al_end <= left_center);
      left_tail = nodes[left_map_ind].
	region_size(left_common_al_end, left_center, left_orient);
    } 

    cout<<"left_tail:"<<left_tail<<endl;
    distance += left_tail;

    //right_map:
    distance += 0.5*nodes[right_map_ind].
      region_size(right_center, right_center+1, right_orient);
    
    double right_tail;

    if(right_center >= right_common_al_end){
      right_tail = nodes[right_map_ind].
	region_size(right_common_al_end, right_center, right_orient);
    }
    else{
      assert(right_common_al_begin >= right_center);
      right_tail = nodes[right_map_ind].
	region_size(right_center+1,right_common_al_begin,right_orient);
    }

    cout<<"right_tail:"<<right_tail<<endl;
    distance += right_tail;
  }
  else{
    cout<<"no common region"<<endl;
    cout<<"map1_ind: "<<map1_ind<<" map2_ind: "<<map2_ind<<endl;
    //debug start
    cout<<"left_map:"<<endl;


    int left_al_begin = edges[cur_edge].left_al_sites[0];
    int left_al_end = 
      edges[cur_edge].left_al_sites[edges[cur_edge].left_al_sites.size()-1];

    int right_al_begin = edges[cur_edge].right_al_sites[0];
    int right_al_end = 
      edges[cur_edge].right_al_sites[edges[cur_edge].right_al_sites.size()-1];

    if(!(left_center >= left_al_begin && left_center<=left_al_end) ||
       !(right_center>=right_al_begin && right_center<=right_al_end)) return left_al_begin*right_al_end*10000; //error

    assert(left_center>=left_al_begin && left_center<=left_al_end);
    assert(right_center>=right_al_begin && right_center<=right_al_end);

    int left_map_left_al_site;
    int left_map_right_al_site;

    int right_map_left_al_site;
    int right_map_right_al_site;

    bool left_info_found = false;
    bool right_info_found = false;

    int cur_left_al_site = 0;
    int cur_right_al_site = 0;
    int prev_left_al_site = cur_left_al_site;
    int prev_right_al_site = cur_left_al_site;

    for(int k=0; k<edges[cur_edge].left_al_sites.size(); k++){
      cur_left_al_site = edges[cur_edge].left_al_sites[k];
      cur_right_al_site = edges[cur_edge].right_al_sites[k];
      
      //cout<<"left: "<<prev_left_al_site<<" -> "<<cur_left_al_site<<endl;
      //cout<<"right: "<<prev_left_al_site<<" -> "<<cur_right_al_site<<endl;

      if(left_center>=prev_left_al_site &&
	 left_center<cur_left_al_site &&
	 !left_info_found){
	//cout<<"left_info_found"<<endl;
	left_info_found = true;
	left_map_left_al_site = prev_left_al_site;
	left_map_right_al_site = cur_left_al_site;
	assert(right_center>=prev_right_al_site &&
	       right_center<cur_right_al_site);
	right_info_found = true;
	right_map_left_al_site = prev_right_al_site;
	right_map_right_al_site = cur_right_al_site;
      }

      prev_left_al_site = cur_left_al_site;
      prev_right_al_site = cur_right_al_site;
    }

    if(!left_info_found || !right_info_found) return 100000; //error

    assert(left_info_found && right_info_found);
    assert(left_map_left_al_site < left_map_right_al_site);
    assert(right_map_left_al_site < right_map_right_al_site);
    
    cout<<"left center:"<<left_center<<" in ("<<left_map_left_al_site;
    cout<<","<<left_map_right_al_site<<")"<<endl;
    cout<<"right center:"<<right_center<<" in ("<<right_map_left_al_site;
    cout<<","<<right_map_right_al_site<<")"<<endl;
    double left_distance_to_left_al_site = 
      0.5* nodes[left_map_ind].
      region_size(left_center, left_center+1, left_orient) 
      + nodes[left_map_ind].
      region_size(left_map_left_al_site, left_center,left_orient);
    double left_distance_to_right_al_site = 
      0.5* nodes[left_map_ind].
      region_size(left_center, left_center+1, left_orient) 
      + nodes[left_map_ind].
      region_size(left_center+1, left_map_right_al_site, left_orient);
    double right_distance_to_left_al_site =
      0.5* nodes[right_map_ind].
      region_size(right_center, right_center+1, right_orient)
      + nodes[right_map_ind].
      region_size(right_map_left_al_site, right_center, right_orient);
    double right_distance_to_right_al_site = 
      0.5* nodes[right_map_ind].
      region_size(right_center, right_center+1, right_orient)
      + nodes[right_map_ind].
      region_size(right_center+1,right_map_right_al_site, right_orient);

    if(left_distance_to_left_al_site < right_distance_to_left_al_site){
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
      else{
	if(same_orient_as_overlap){
	  dist_sign = -1;
	}
	else{
	  dist_sign = 1;
	}
      }
    }
    else{
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = -1;
	else
	  dist_sign = 1;
      }
      else{
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
    }


    distance = 0.5*(fabs(left_distance_to_left_al_site-
			 right_distance_to_left_al_site)+
		    fabs(left_distance_to_right_al_site-
			 right_distance_to_right_al_site));
    
    //assert(false);
  }

  //assert(distance > 0);
  cout<<"distance: "<<dist_sign*distance<<endl<<endl;
  return dist_sign*distance;
}
*/
double graph_component::
distance_between_maps(int cur_edge_ind, int map1_ind, int map2_ind, 
		      orientation map1_orient, orientation map2_orient){
  //the distance is given from map1 to map2, not the other way around

  assert(cur_edge_ind >=0 && cur_edge_ind < edges.size());
  assert(map1_ind >= 0 && map1_ind < nodes.size());
  assert(map2_ind >= 0 && map2_ind < nodes.size());
  assert(edges[cur_edge_ind].left_map == nodes[map1_ind].read_name ||
	 edges[cur_edge_ind].right_map == nodes[map1_ind].read_name);
  assert(edges[cur_edge_ind].left_map == nodes[map2_ind].read_name ||
	 edges[cur_edge_ind].right_map == nodes[map2_ind].read_name);

  double distance = 0;
  int dist_sign=1;
  int cur_edge = cur_edge_ind;

  bool same_orient_as_overlap;
  
  if(nodes[map1_ind].read_name == edges[cur_edge].left_map){
    assert((map1_orient == edges[cur_edge].left_orient &&
	    map2_orient == edges[cur_edge].right_orient) ||
	   (rev(map1_orient) == edges[cur_edge].left_orient &&
	    rev(map2_orient) == edges[cur_edge].right_orient));
    if(nodes[map1_ind].orient == edges[cur_edge].left_orient){
      same_orient_as_overlap = true;
      cout<<"same orientation as the edge"<<endl;
    }
    else{
      same_orient_as_overlap = false;
      cout<<"reverse orientation relative to the edge"<<endl;
    }
  }
  else{
    assert((map1_orient == edges[cur_edge].right_orient &&
	    map2_orient == edges[cur_edge].left_orient) ||
	   (rev(map1_orient) == edges[cur_edge].right_orient &&
	    rev(map2_orient) == edges[cur_edge].left_orient));
    if(nodes[map1_ind].orient == edges[cur_edge].right_orient){
      same_orient_as_overlap = true;
      cout<<"same orientation as the edge"<<endl;
    }
    else{
      same_orient_as_overlap = false;
      cout<<"reverse orientation relative to the edge"<<endl;
    }
  }


  orientation left_orient = edges[cur_edge].left_orient;
  orientation right_orient = edges[cur_edge].right_orient;

  int left_map_ind;
  int right_map_ind;

  int left_center;
  int right_center;
 
  bool same_order;
  if(nodes[map1_ind].read_name == edges[cur_edge].left_map){
    //map1 is on the left side of the overlap

    assert(nodes[map2_ind].read_name == edges[cur_edge].right_map);
    //everything should be consistent
    
    left_center = nodes[map1_ind].map_center(edges[cur_edge].left_orient);
    right_center = nodes[map2_ind].map_center(edges[cur_edge].right_orient);
    
    left_map_ind = map1_ind;
    right_map_ind = map2_ind;

    same_order = true;
    
    assert(nodes[map2_ind].read_name == edges[cur_edge].right_map);
    //everything should be consistent
  }
  else{
    assert(nodes[map1_ind].read_name == edges[cur_edge].right_map);
    assert(nodes[map2_ind].read_name == edges[cur_edge].left_map);
    
    left_center = nodes[map2_ind].map_center(edges[cur_edge].left_orient);
    right_center = nodes[map1_ind].map_center(edges[cur_edge].right_orient);
    
    left_map_ind = map2_ind;
    right_map_ind = map1_ind;  
    
    same_order = false;
  }

  assert(map1_ind != -1 && map2_ind != -1);

  cout<<"overlap: "<<left_map_ind<<"->"<<right_map_ind;
  cout<<" "<<nodes[left_map_ind].map_read.size();
  cout<<" "<<nodes[right_map_ind].map_read.size();
  cout<<" lo:"<<edges[cur_edge].left_orient;
  cout<<" ro:"<<edges[cur_edge].right_orient<<endl;
  cout<<"("<<edges[cur_edge].left_al_sites[0]<<",";
  cout<<edges[cur_edge].left_al_sites[edges[cur_edge].left_al_sites.size()-1];
  cout<<") -> ("<<edges[cur_edge].right_al_sites[0]<<",";
  cout<<edges[cur_edge].right_al_sites[edges[cur_edge].right_al_sites.size()-1]<<")"<<endl;
  cout<<"left_center:"<<left_center;
  cout<<", right_center:"<<right_center<<endl;
  //cout<<endl;

  int left_common_al_begin;
  int left_common_al_end;
  int right_common_al_begin;
  int right_common_al_end;
  
  bool common_region_exists = false;
  bool left_end_found = false;
  for(int i=0; i<edges[cur_edge].left_al_sites.size(); i++){
    int cur_left_site_ind = edges[cur_edge].left_al_sites[i];
    int cur_right_site_ind =  edges[cur_edge].right_al_sites[i];

    if((left_center < cur_left_site_ind &&
	right_center >= cur_right_site_ind) ||
       (left_center >= cur_left_site_ind &&
	right_center < cur_right_site_ind)){
      //there is aligned site between two centers

      if(!left_end_found){
	left_end_found = true;
	common_region_exists = true;
	
	left_common_al_begin = cur_left_site_ind;
	right_common_al_begin = cur_right_site_ind;
      }
      
      left_common_al_end = cur_left_site_ind;
      right_common_al_end = cur_right_site_ind;
    }
  }

 
  if(common_region_exists){

    if(left_center <= left_common_al_begin){
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
      else{
	if(same_orient_as_overlap){
	  dist_sign = -1;
	}
	else{
	  dist_sign = 1;
	}
      }
    }
    else{
      assert(left_center >= left_common_al_end);
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = -1;
	else
	  dist_sign = 1;
      }
      else{
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
    }
    //else{
    //assert(left_center >= common_al_end);
    //if(same_orient_as_overlap) dist_sign = -1;
    //else dist_sign = 1;
    //}

    cout<<"dist_sign:"<<dist_sign<<endl;

    cout<<"common al region: (";
    cout<<left_common_al_begin<<","<<left_common_al_end<<")->(";
    cout<<right_common_al_begin<<","<<right_common_al_end<<")"<<endl;

    assert(left_common_al_begin <= left_common_al_end);
    assert(right_common_al_begin <= right_common_al_end);

    cout<<"left map region:"<<left_common_al_begin;
    cout<<" to: "<<left_common_al_end<<endl;
    double common_left_map_region = 
      nodes[left_map_ind].region_size(left_common_al_begin, left_common_al_end, left_orient);

    cout<<"right map region:"<<right_common_al_begin;
    cout<<" to: "<<right_common_al_end<<endl;
    double common_right_map_region = 
      nodes[right_map_ind].region_size(right_common_al_begin, right_common_al_end, right_orient);
    
    cout<<"left common region size:"<<common_left_map_region<<endl;
    cout<<"right common region size:"<<common_right_map_region<<endl;
    distance += 0.5*(common_left_map_region+common_right_map_region);
    //add the distance of the common region

    //now add the tails
    
    //left_map:
    distance += 0.5*
      nodes[left_map_ind].region_size(left_center, left_center+1, left_orient);
    
    double left_tail;
    if(left_center < left_common_al_begin){
      left_tail = nodes[left_map_ind].
	region_size(left_center+1, left_common_al_begin, left_orient);
    }
    else{
      assert(left_common_al_end <= left_center);
      left_tail = nodes[left_map_ind].
	region_size(left_common_al_end, left_center, left_orient);
    } 

    cout<<"left_tail:"<<left_tail<<endl;
    distance += left_tail;

    //right_map:
    distance += 0.5*nodes[right_map_ind].
      region_size(right_center, right_center+1, right_orient);
    
    double right_tail;

    if(right_center >= right_common_al_end){
      right_tail = nodes[right_map_ind].
	region_size(right_common_al_end, right_center, right_orient);
    }
    else{
      assert(right_common_al_begin >= right_center);
      right_tail = nodes[right_map_ind].
	region_size(right_center+1,right_common_al_begin,right_orient);
    }

    cout<<"right_tail:"<<right_tail<<endl;
    distance += right_tail;
  }
  else{
    cout<<"no common region"<<endl;
    cout<<"map1_ind: "<<map1_ind<<" map2_ind: "<<map2_ind<<endl;
    //debug start
    cout<<"left_map:"<<endl;
    /*
    for(int h=0; h<nodes[map1_ind].size(); h++){
      cout<<" "<<nodes[map1_ind].map_read[h];
    }
    cout<<endl<<endl;
    cout<<"right_map:"<<endl;
    for(int h=0; h<nodes[map2_ind].size(); h++){
      cout<<" "<<nodes[map2_ind].map_read[h];
    }
    cout<<endl<<endl;
    */
    //debug end

    int left_al_begin = edges[cur_edge].left_al_sites[0];
    int left_al_end = 
      edges[cur_edge].left_al_sites[edges[cur_edge].left_al_sites.size()-1];

    int right_al_begin = edges[cur_edge].right_al_sites[0];
    int right_al_end = 
      edges[cur_edge].right_al_sites[edges[cur_edge].right_al_sites.size()-1];

    if(!(left_center >= left_al_begin && left_center<=left_al_end) ||
       !(right_center>=right_al_begin && right_center<=right_al_end)) return left_al_begin*right_al_end*10000; //error

    assert(left_center>=left_al_begin && left_center<=left_al_end);
    assert(right_center>=right_al_begin && right_center<=right_al_end);

    int left_map_left_al_site;
    int left_map_right_al_site;

    int right_map_left_al_site;
    int right_map_right_al_site;

    bool left_info_found = false;
    bool right_info_found = false;

    int cur_left_al_site = 0;
    int cur_right_al_site = 0;
    int prev_left_al_site = cur_left_al_site;
    int prev_right_al_site = cur_left_al_site;

    for(int k=0; k<edges[cur_edge].left_al_sites.size(); k++){
      cur_left_al_site = edges[cur_edge].left_al_sites[k];
      cur_right_al_site = edges[cur_edge].right_al_sites[k];
      
      //cout<<"left: "<<prev_left_al_site<<" -> "<<cur_left_al_site<<endl;
      //cout<<"right: "<<prev_left_al_site<<" -> "<<cur_right_al_site<<endl;

      if(left_center>=prev_left_al_site &&
	 left_center<cur_left_al_site &&
	 !left_info_found){
	//cout<<"left_info_found"<<endl;
	left_info_found = true;
	left_map_left_al_site = prev_left_al_site;
	left_map_right_al_site = cur_left_al_site;
	assert(right_center>=prev_right_al_site &&
	       right_center<cur_right_al_site);
	right_info_found = true;
	right_map_left_al_site = prev_right_al_site;
	right_map_right_al_site = cur_right_al_site;
      }

      prev_left_al_site = cur_left_al_site;
      prev_right_al_site = cur_right_al_site;
    }

    if(!left_info_found || !right_info_found) return 100000; //error

    assert(left_info_found && right_info_found);
    assert(left_map_left_al_site < left_map_right_al_site);
    assert(right_map_left_al_site < right_map_right_al_site);
    
    cout<<"left center:"<<left_center<<" in ("<<left_map_left_al_site;
    cout<<","<<left_map_right_al_site<<")"<<endl;
    cout<<"right center:"<<right_center<<" in ("<<right_map_left_al_site;
    cout<<","<<right_map_right_al_site<<")"<<endl;
    double left_distance_to_left_al_site = 
      0.5* nodes[left_map_ind].
      region_size(left_center, left_center+1, left_orient) 
      + nodes[left_map_ind].
      region_size(left_map_left_al_site, left_center,left_orient);
    double left_distance_to_right_al_site = 
      0.5* nodes[left_map_ind].
      region_size(left_center, left_center+1, left_orient) 
      + nodes[left_map_ind].
      region_size(left_center+1, left_map_right_al_site, left_orient);
    double right_distance_to_left_al_site =
      0.5* nodes[right_map_ind].
      region_size(right_center, right_center+1, right_orient)
      + nodes[right_map_ind].
      region_size(right_map_left_al_site, right_center, right_orient);
    double right_distance_to_right_al_site = 
      0.5* nodes[right_map_ind].
      region_size(right_center, right_center+1, right_orient)
      + nodes[right_map_ind].
      region_size(right_center+1,right_map_right_al_site, right_orient);

    if(left_distance_to_left_al_site < right_distance_to_left_al_site){
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
      else{
	if(same_orient_as_overlap){
	  dist_sign = -1;
	}
	else{
	  dist_sign = 1;
	}
      }
    }
    else{
      if(same_order){
	if(same_orient_as_overlap)
	  dist_sign = -1;
	else
	  dist_sign = 1;
      }
      else{
	if(same_orient_as_overlap)
	  dist_sign = 1;
	else
	  dist_sign = -1;
      }
    }


    distance = 0.5*(fabs(left_distance_to_left_al_site-
			 right_distance_to_left_al_site)+
		    fabs(left_distance_to_right_al_site-
			 right_distance_to_right_al_site));
    
    //assert(false);
  }

  //assert(distance > 0);
  cout<<"distance: "<<dist_sign*distance<<endl<<endl;
  return dist_sign*distance;
}
void graph_component::assign_edge_weights(){
  //distance between maps within overlaps
  //is given by the distance between centers
  //of two maps within a given overlap
  //the center of each map
  //is given by the center of the fragment
  //within which the map center falls

  for(int i=0; i<nodes.size(); i++){
    assert(nodes[i].structures_ok());
    nodes[i].right_edge_weights.clear();
    for(int j=0; j<nodes[i].right_maps_inds.size(); j++){      
      int cur_edge_ind = nodes[i].right_edges_inds[j];
      double cur_weight = distance_between_maps
	(cur_edge_ind, i,nodes[i].right_maps_inds[j],
	 nodes[i].orient, nodes[nodes[i].right_maps_inds[j]].orient);
      nodes[i].right_edge_weights.push_back(cur_weight);
    }
    nodes[i].left_edge_weights.clear();
    for(int j=0; j<nodes[i].left_maps_inds.size(); j++){
      int cur_edge_ind = nodes[i].left_edges_inds[j];
      double cur_weight = distance_between_maps
	(cur_edge_ind, i,nodes[i].left_maps_inds[j],
	 nodes[i].orient, nodes[nodes[i].left_maps_inds[j]].orient);
      nodes[i].left_edge_weights.push_back(cur_weight);
    }
    nodes[i].contained_edge_weights.clear();
    for(int j=0; j<nodes[i].contained_maps_inds.size(); j++){
      int cur_edge_ind = nodes[i].contained_edges_inds[j];
      double cur_weight = distance_between_maps
	(cur_edge_ind, i,nodes[i].contained_maps_inds[j],
	 nodes[i].orient,nodes[nodes[i].contained_maps_inds[j]].orient);
      nodes[i].contained_edge_weights.push_back(cur_weight);
    }
    nodes[i].container_edge_weights.clear();
    for(int j=0; j<nodes[i].container_maps_inds.size(); j++){
      int cur_edge_ind = nodes[i].container_edges_inds[j];
      double cur_weight = distance_between_maps
	(cur_edge_ind, i,nodes[i].container_maps_inds[j],
	 nodes[i].orient,nodes[nodes[i].container_maps_inds[j]].orient);
      nodes[i].container_edge_weights.push_back(cur_weight);
    }

  }  
}
void graph_component::index_nodes(){
  for(int i=0; i<nodes.size(); i++){
    nodes[i]._index = i;
  }
}
/*
void graph_component::index_nodes(){
  
  for(int i=0; i<nodes.size(); i++){
    //if(i%100 == 0) cout<<"indexing nodes: "<<i<<" out of "<<nodes.size()<<endl;

    nodes[i]._index = i;
    
    nodes[i].right_maps_inds.clear();    

    for(int j=0; j<nodes[i].right_maps.size(); j++){      
      bool ind_found = false;
      int cur_ind = -1;

      for(int k=0; k<nodes.size() && !ind_found; k++){
	if(nodes[i].right_maps[j] == nodes[k].read_name){
	  ind_found = true;
	  cur_ind = k;
	  nodes[i].right_maps_inds.push_back(k);
	}
      }
      assert(ind_found);
    }

    
    nodes[i].left_maps_inds.clear();    
    for(int j=0; j<nodes[i].left_maps.size(); j++){
      bool ind_found = false;
      int cur_ind = -1;
      
      for(int k=0; k<nodes.size() && !ind_found; k++){
	if(nodes[i].left_maps[j] == nodes[k].read_name){
	  ind_found = true;
	  cur_ind = k;
	  nodes[i].left_maps_inds.push_back(k);	  
	}
      }
      assert(ind_found);
    }

    nodes[i].contained_maps_inds.clear();
    for(int j=0; j<nodes[i].contained_maps.size(); j++){
      bool ind_found = false;
      int cur_ind = -1;
      
      for(int k=0; k<nodes.size() && !ind_found; k++){
	if(nodes[i].contained_maps[j] == nodes[k].read_name){
	  ind_found = true;
	  cur_ind = k;
	  nodes[i].contained_maps_inds.push_back(k);	  
	}
      }
      assert(ind_found);
    }

    nodes[i].container_maps_inds.clear();
    for(int j=0; j<nodes[i].container_maps.size(); j++){
      bool ind_found = false;
      int cur_ind = -1;
      
      for(int k=0; k<nodes.size() && !ind_found; k++){
	if(nodes[i].container_maps[j] == nodes[k].read_name){
	  ind_found = true;
	  cur_ind = k;
	  nodes[i].container_maps_inds.push_back(k);	  
	}
      }
      assert(ind_found);
    }

    nodes[i].right_edges_inds.clear();
    for(int j=0; j<nodes[i].right_maps.size(); j++){
      int cur_edge = edge_ind(nodes[i].read_name, nodes[i].right_maps[j]);
      assert(cur_edge != -1);
      nodes[i].right_edges_inds.push_back(cur_edge);
    }

    
    nodes[i].left_edges_inds.clear();
    for(int j=0; j<nodes[i].left_maps.size(); j++){
      int cur_edge = edge_ind(nodes[i].read_name, nodes[i].left_maps[j]);
      assert(cur_edge != -1);
      nodes[i].left_edges_inds.push_back(cur_edge);
    }

    nodes[i].contained_edges_inds.clear();
    for(int j=0; j<nodes[i].contained_maps.size(); j++){
      int cur_edge = edge_ind(nodes[i].read_name, nodes[i].contained_maps[j]);
      assert(cur_edge != -1);
      nodes[i].contained_edges_inds.push_back(cur_edge);
    }
    
    nodes[i].container_edges_inds.clear();
    for(int j=0; j<nodes[i].container_maps.size(); j++){
      int cur_edge = edge_ind(nodes[i].read_name, nodes[i].container_maps[j]);
      //assert(cur_edge != -1);
      if(cur_edge == -1){
	cout<<"edge list inconsistent with node info."<<endl;
	cout<<"processing continers:"<<endl;
	cout<<"master: "<<get_node_index(nodes[i].container_maps[j]);
	cout<<" slave: "<<get_node_index(nodes[i].read_name)<<endl;
	cout<<"listing overlaps:"<<endl;
	for(int k=0; k<edges.size(); k++){
	  cout<<get_node_index(edges[k].left_map)<<" -> ";
	  cout<<get_node_index(edges[k].right_map)<<endl;
	}
	assert(false);
	//better not be here
	
      }
      nodes[i].container_edges_inds.push_back(cur_edge);
    }
  }
}
*/

bool graph_component::dfs_all_way_visit(int n, int to_reach){
  assert(to_reach < nodes.size());
  assert(n>=0 && n<nodes.size());
  assert(n<node_color.size());
  assert(node_color.size() == nodes.size());
  assert(preds.size() == nodes.size());
  
  cout<<" dfs_visit: current node "<<n;
  cout<<" ("<<nodes[n].read_name<<")"<<endl;
  cout<<"right maps:";
  for(int i=0; i<nodes[n].right_maps_inds.size(); i++){
    cout<<" "<<nodes[n].right_maps_inds[i];
  }
  cout<<endl;
  cout<<"left maps:";
  for(int i=0; i<nodes[n].left_maps_inds.size(); i++){
    cout<<" "<<nodes[n].left_maps_inds[i];
  }
  cout<<endl;
  cout<<"contained maps:";
  for(int i=0; i<nodes[n].contained_maps_inds.size(); i++){
    cout<<" "<<nodes[n].contained_maps_inds[i];   
  }
  cout<<endl;
  cout<<"container maps:";
  for(int i=0; i<nodes[n].container_maps_inds.size(); i++){
    cout<<" "<<nodes[n].container_maps_inds[i];
  }
  cout<<endl;
  

  if(n==to_reach) return true;
  else
    if(node_color[n] != 0) return false;
  //there is a cycle, but not reaching the desired destination

  assert(node_color[n] == 0);

  bool result = false;

  node_color[n] = 1; //node is discovered;

  //compiling a list of neighbors
  vector<int> neighbors;
  for(int i=0; i<nodes[n].right_maps_inds.size(); i++){
    neighbors.push_back(nodes[n].right_maps_inds[i]);
  }
  for(int i=0; i<nodes[n].left_maps_inds.size(); i++){
    neighbors.push_back(nodes[n].left_maps_inds[i]);
  }
  for(int i=0; i<nodes[n].contained_maps_inds.size(); i++){
    neighbors.push_back(nodes[n].contained_maps_inds[i]);
  }
  for(int i=0; i<nodes[n].container_maps_inds.size(); i++){
    neighbors.push_back(nodes[n].container_maps_inds[i]);
  }

  for(int i=0; i<neighbors.size(); i++){
    int cur_ind = neighbors[i];
    assert(cur_ind < preds.size());
    assert(cur_ind >= 0);
    
    if(node_color[cur_ind] == 0){
      preds[cur_ind] = n;
      bool cur_res = dfs_all_way_visit(cur_ind, to_reach);      
      
      if(cur_res == true){
	result = true;
	return true;
      }
    }
  }
  node_color[n] = 2; //all the children explored
  assert(result == false);
  return result;
}

bool graph_component::dfs_left_to_right_visit(int n, int to_reach){
  assert(to_reach < nodes.size());
  assert(n>=0 && n<nodes.size());
  assert(n<node_color.size());
  assert(node_color.size() == nodes.size());
  assert(preds.size() == nodes.size());
  
  cout<<" dfs_visit: current node "<<n;
  cout<<" ("<<nodes[n].read_name<<")"<<endl;
  cout<<"right maps:";
  for(int i=0; i<nodes[n].right_maps_inds.size(); i++){
    cout<<" "<<nodes[n].right_maps_inds[i];
  }
  cout<<endl;

  if(n==to_reach) return true;
  else
    if(node_color[n] != 0) return false;
  //there is a cycle, but not reaching the desired destination

  assert(node_color[n] == 0);

  bool result = false;

  node_color[n] = 1; //node is discovered;
  for(int i=0; i<nodes[n].right_maps_inds.size(); i++){
    int cur_ind = nodes[n].right_maps_inds[i];
    assert(cur_ind < preds.size());
    assert(cur_ind >= 0);
    
    if(node_color[cur_ind] == 0){
      preds[cur_ind] = n;
      bool cur_res = dfs_left_to_right_visit(cur_ind, to_reach);      
      
      if(cur_res == true){
	result = true;
	return true;
      }
    }
  }
  node_color[n] = 2; //all the children explored
  assert(result == false);
  return result;
}

/*
int graph_component::get_node_index(string node_id){
  bool node_found = false;
  for(int i=0; i<nodes.size(); i++){
    if(nodes[i].read_name == node_id){
      node_found = true;
      return i;
    }
  }
  
  assert(false);
  return -1;
}
*/

/*
bool graph_component::left_to_right_path(string from_node, string to_node){
  
  //index_nodes();

  int from_node_ind = get_node_index(from_node);
  int to_node_ind = get_node_index(to_node);

  preds.clear();
  node_color.clear();

  for(int i=0; i<nodes.size(); i++){
    preds.push_back(-1);
    node_color.push_back(0);
  }

  for(int i=0; i<nodes[from_node_ind].right_maps_inds.size(); i++){
    int cur_succ_ind = nodes[from_node_ind].right_maps_inds[i];
    bool path_exists = dfs_left_to_right_visit(cur_succ_ind, to_node_ind);
    
    if(path_exists){
      int it_counter = 0;
      
      bool from_node_reached = false;
      //backtracing to discover exact path
      int cur_node_ind = to_node_ind;
      while(!from_node_reached){
	vector<int>::iterator vec_begin = preds.begin();
	helper1.insert(vec_begin, cur_node_ind);
	if(cur_node_ind == cur_succ_ind) from_node_reached = true;
	
	cur_node_ind = preds[cur_node_ind];
	it_counter++;
	assert(it_counter < 1000); //make sure end is reached
      }
      vector<int>::iterator vec_begin = preds.begin();
      preds.insert(vec_begin,from_node_ind);
    }
  }
  //no path discovered
  return false;
}
*/
/*
bool graph_component::all_way_path(string from_node, string to_node){
  
  //index_nodes();

  int from_node_ind = get_node_index(from_node);
  int to_node_ind = get_node_index(to_node);

  preds.clear();
  node_color.clear();

  for(int i=0; i<nodes.size(); i++){
    preds.push_back(-1);
    node_color.push_back(0);
  }

  cout<<"looking for path from node: "<<from_node_ind;
  cout<<" to node: "<<to_node_ind<<endl;

  cout<<"current node: "<<from_node_ind;
  cout<<" ("<<nodes[from_node_ind].read_name<<")"<<endl;
  cout<<"right maps:";
  for(int i=0; i<nodes[from_node_ind].right_maps_inds.size(); i++){
    cout<<" "<<nodes[from_node_ind].right_maps_inds[i];
  }
  cout<<endl;
  cout<<"left_maps:";
  for(int i=0; i<nodes[from_node_ind].left_maps_inds.size(); i++){
    cout<<" "<<nodes[from_node_ind].left_maps_inds[i];
  }
  cout<<endl;
  cout<<"contained maps:";
  for(int i=0; i<nodes[from_node_ind].contained_maps_inds.size(); i++){
    cout<<" "<<nodes[from_node_ind].contained_maps_inds[i];
  }
  cout<<endl;
  cout<<"container maps:";
  for(int i=0; i<nodes[from_node_ind].container_maps_inds.size(); i++){
    cout<<" "<<nodes[from_node_ind].container_maps_inds[i];
  }
  cout<<endl;

  //compiling the neighbor list
  vector<int> neighbors;
  for(int i=0; i<nodes[from_node_ind].right_maps_inds.size(); i++)
    neighbors.push_back(nodes[from_node_ind].right_maps_inds[i]);
  for(int i=0; i<nodes[from_node_ind].left_maps_inds.size(); i++)
    neighbors.push_back(nodes[from_node_ind].left_maps_inds[i]);
  for(int i=0; i<nodes[from_node_ind].contained_maps_inds.size(); i++)
    neighbors.push_back(nodes[from_node_ind].contained_maps_inds[i]);
  for(int i=0; i<nodes[from_node_ind].container_maps_inds.size(); i++)
    neighbors.push_back(nodes[from_node_ind].container_maps_inds[i]);
  
  node_color[from_node_ind]=1;
  for(int i=0; i<neighbors.size(); i++){	     
    int cur_neighbor_ind = neighbors[i];

    node_color[cur_neighbor_ind]=1;
    preds[cur_neighbor_ind] = from_node_ind;
    bool path_exists = dfs_all_way_visit(cur_neighbor_ind, to_node_ind);
    node_color[cur_neighbor_ind]=2;

    if(path_exists){
      cout<<"path exists"<<endl;

      int it_counter = 0;
      
      bool from_node_reached = false;
      //backtracing to discover exact path
      int cur_node_ind = to_node_ind;
      while(!from_node_reached){
	vector<int>::iterator vec_begin = preds.begin();
	helper1.insert(vec_begin, cur_node_ind);
	if(cur_node_ind == cur_neighbor_ind) from_node_reached = true;
	
	cur_node_ind = preds[cur_node_ind];
	it_counter++;
	assert(it_counter < 1000); //make sure end is reached
      }
      vector<int>::iterator vec_begin = preds.begin();
      preds.insert(vec_begin,from_node_ind);
    }
  }
  //no path discovered
  return false;
}
*/

bool graph_component::cycle(int n){
  cout<<endl;
  cout<<"identifying cycles for node "<<n<<endl;

  assert(n<nodes.size() && n>=0);

  int n_ind = n;

  index_nodes(); //make sure everything is indexed

  preds.clear();
  //use preds to store predecessors of each node
  node_color.clear();

  for(int i=0; i<nodes.size(); i++){
    preds.push_back(-1);
    node_color.push_back(0);
  }
  
  bool parent_reached = false;
  int cycle_child = -1;
  for(int i=0; i<nodes[n_ind].right_maps_inds.size() && !parent_reached; i++){
    int cur_succ = nodes[n_ind].right_maps_inds[i];
    cout<<"processing succ: "<<cur_succ;
    //cout<<" ("<<nodes[cur_succ].read_name<<")"<<endl;
    cout<<endl;

    preds[cur_succ] = n_ind;
    bool cur_reach = dfs_left_to_right_visit(cur_succ, n_ind);
    if(cur_reach == true){
      parent_reached = true;
      cycle_child = cur_succ;
    }
    //cout<<"pred of "<<nodes[cur_succ].read_name;
    //cout<<" is "<<nodes[preds[cur_succ]].read_name<<endl;
  }

  if(parent_reached){
    cout<<"found a cycle for node: "<<n_ind;
    cout<<" through "<<cycle_child<<endl;
  }
  
  vector<int> pred_helper;

  if(parent_reached){
    //trace back and extract the cycle

    bool cycle_extracted = false;
    
    int cur_node_ind = n_ind;
    int cur_counter = 0;

    //pred_helper.push_back(cur_node_ind);
    while(!cycle_extracted){
      assert(cur_counter < 1000); //sanity check
      cur_counter++;
      int cur_pred = preds[cur_node_ind];
      assert(cur_pred >= 0);
      assert(cur_pred<nodes.size());
     
      cout<<"n:"<<n_ind;
      cout<<" for_node: "<<cur_node_ind;
      cout<<" ("<<nodes[cur_node_ind].read_name<<")";
      cout<<" pred: "<<cur_pred;
      cout<<" ("<<nodes[cur_pred].read_name<<")"<<endl;

      if(cur_pred == n_ind) cycle_extracted = true;
      
      pred_helper.push_back(cur_pred);
	

      cur_node_ind = cur_pred;
    }
   

    cout<<"done extracting the cycle"<<endl<<endl;

    cycled_path.clear();
    for(int i=pred_helper.size()-1; i>=0; i--){
      cycled_path.push_back(pred_helper[i]);
    }

    cout<<"cycle:"<<endl;
    for(int i=0; i<cycled_path.size(); i++){
      if(i!=0) cout<<" -> ";
      cout<<cycled_path[i];
    }
    cout<<endl;
  }

  return parent_reached;
}

vector<double> graph_component::extract_contig(vector<int>& map_inds){
  cout<<endl;
  cout<<"extracting the contig:"<<endl;

  vector<double> contig;

  if(map_inds.empty()) return contig;

  int cur_map_ind = -1;
  int next_map_ind = -1;
  
  int first_al_index = 0;
  assert(map_inds.size()>=2);


  for(int i=0; i<map_inds.size()-1; i++){

    cur_map_ind = map_inds[i];
    next_map_ind = map_inds[i+1];

    assert(cur_map_ind >= 0 && cur_map_ind < nodes.size());
    assert(next_map_ind >= 0 && next_map_ind < nodes.size());

    string cur_map_name = nodes[cur_map_ind].read_name;
    string next_map_name = nodes[next_map_ind].read_name;

    //locate a corresponding edge

    bool edge_found = false;
    int edge_index = -1;
    for(int j=0; j<edges.size() && !edge_found; j++){
      if((edges[j].left_map == cur_map_name && 
	  edges[j].right_map == next_map_name) ||
	 (edges[j].left_map == next_map_name &&
	  edges[j].right_map == cur_map_name)){
	edge_found = true;
	edge_index = j;
      }
    }

    assert(edge_found == true);
    assert(!edges[edge_index].containment);
    edge cur_e;

    if(edges[edge_index].left_map == cur_map_name){
      cur_e = edges[edge_index];
    }
    else{
      assert(edges[edge_index].right_map == cur_map_name);
      cur_e = edges[edge_index].inverse();
    }

    int last_al_index = cur_e.left_al_sites[cur_e.left_al_sites.size()-1];

    cout<<"overlap: "<<cur_map_ind<<" -> "<<next_map_ind;
    cout<<"   "<<cur_map_name<<" -> "<<next_map_name<<endl;
    edges[edge_index].print();
    cout<<"first_al_index: "<<first_al_index;
    cout<<" last_al_index: "<<last_al_index<<endl<<endl;

    output_overlap(cur_e);

    assert(first_al_index <= last_al_index); //take care of this later
    
    if(cur_e.left_orient == 1){
      for(int k=first_al_index; k<last_al_index; k++){
	contig.push_back(nodes[cur_map_ind].map_read[k]);
      }
    }
    else{     
      int _start = nodes[cur_map_ind].map_read.size()-first_al_index;
      int _end = nodes[cur_map_ind].map_read.size()-last_al_index;

      cout<<"in inverse orient: adding sites from ";
      cout<<_start<<" to "<<_end<<endl;;

      for(int k=_start-1; k>=_end; k--){
	cout<<"k: "<<k;
	cout<<" adding1 : "<<nodes[cur_map_ind].map_read[k]<<endl;
	contig.push_back(nodes[cur_map_ind].map_read[k]);	      
      }
    }

    first_al_index = cur_e.right_al_sites[cur_e.right_al_sites.size()-1];

    if(i==map_inds.size()-2){
      //add the fragments from the last map
      if(cur_e.right_orient == 1){
	for(int k=first_al_index; k<nodes[next_map_ind].map_read.size(); k++){
	  contig.push_back(nodes[next_map_ind].map_read[k]);
	}
      }
      else{
	for(int k=nodes[next_map_ind].map_read.size()-first_al_index-1; 
	    k>=0; k--){
	  contig.push_back(nodes[next_map_ind].map_read[k]);
	}
      }
    }
  }
  return contig;
}

vector<double> graph_component::extract_contig(){
  for(int i=0; i<cycled_path.size(); i++){
    cout<<cycled_path[i]<<endl;
  }
  vector<double> cur_contig;
  if(!cycled_path.empty()){
    cur_contig = extract_contig(cycled_path);
  }
  return cur_contig;
}

void graph_component::output_overlap(edge& e){
  bool left_map_found = false;
  bool right_map_found = false;

  int left_map_ind = -1;
  int right_map_ind = -1;

  for(int i=0; i<nodes.size() && !left_map_found; i++){
    if(e.left_map == nodes[i].read_name){
      left_map_found = true;
      left_map_ind = i;
    }
  }

  for(int i=0; i<nodes.size() && !right_map_found; i++){
    if(e.right_map == nodes[i].read_name){
      right_map_found = true;
      right_map_ind = i;
    }
  }

  cout<<"outputting alignment: "<<endl;
  cout<<e.left_map<<" -> "<<e.right_map;
  cout<<" "<<e.left_orient<<" "<<e.right_orient;
  cout<<endl;

  if(!left_map_found || !right_map_found){
    cout<<"failed to find one of the maps among the nodes"<<endl;
    return;
  }
  else{
    if(e.left_al_sites.size() > 1){
      int first_left_al_site = e.left_al_sites[0];
      for(int i=0; i<first_left_al_site; i++){
	cout<<"[ "<<i<<":";
	if(e.left_orient == 1){
	  cout<<nodes[left_map_ind].map_read[i];
	}
	else{
	  cout<<nodes[left_map_ind].map_read
	    [nodes[left_map_ind].map_read.size()-i-1];
	}
	cout<<" ]"<<endl;
      }

      for(int i=0; i<e.left_al_sites.size()-1; i++){
	int cur_left_al_site = e.left_al_sites[i];
	int next_left_al_site = e.left_al_sites[i+1];

	int cur_right_al_site = e.right_al_sites[i];
	int next_right_al_site = e.right_al_sites[i+1];

	cout<<"[ ";
	for(int j=cur_left_al_site; j<next_left_al_site; j++){
	  if(j!=cur_left_al_site){
	    cout<<", ";
	  }
	  cout<<j<<":";
	  if(e.left_orient == 1){
	    cout<<nodes[left_map_ind].map_read[j];
	  }
	  else{
	    cout<<nodes[left_map_ind].map_read[nodes[left_map_ind].
					       map_read.size() - j-1];
	  }
	}
	cout<<" ] -> [ ";
	for(int j=cur_right_al_site; j<next_right_al_site; j++){
	  if(j!= cur_right_al_site){
	    cout<<", ";
	  }
	  cout<<j<<":";
	  if(e.right_orient == 1){
	    cout<<nodes[right_map_ind].map_read[j];
	  }
	  else{
	    cout<<nodes[right_map_ind].map_read[nodes[right_map_ind].
						map_read.size()-j-1];
	  }
	}
	cout<<" ]";
	cout<<endl;
      }

      int last_left_al_site = e.right_al_sites[e.right_al_sites.size()-1];
      int last_left_site = nodes[right_map_ind].map_read.size();

      for(int i=last_left_al_site; i<last_left_site; i++){
	cout<<char(9)<<char(9)<<"[ "<<i<<":";
	if(e.right_orient == 1){
	  cout<<nodes[right_map_ind].map_read[i];
	}
	else{
	  cout<<nodes[right_map_ind].map_read
	    [nodes[right_map_ind].map_read.size()-i-1];
	}
	cout<<" ]"<<endl;
      }	
    }
    cout<<"s_score: "<<e.s_score;
    cout<<" t_score: "<<e.t_score<<endl;
    cout<<endl;
  }
}

/*
int graph_component::edge_ind(int map1_ind, int map2_ind){
  assert(map1_ind >= 0 && map1_ind < nodes.size());
  assert(map2_ind >= 0 && map2_ind < nodes.size());
  
  bool edge_found = false;
  for(int i=0; i<edges.size() && !edge_found; i++){
    if((edges[i].left_map == nodes[map1_ind].read_name &&
	edges[i].right_map == nodes[map2_ind].read_name) ||
       (edges[i].left_map == nodes[map2_ind].read_name &&
	edges[i].right_map == nodes[map1_ind].read_name)){
      edge_found = true;
      return i;
    }
  }
 
  return -1;
}

int graph_component::edge_ind(string map1, string map2){
  bool edge_found = false;
  for(int i=0; i<edges.size() && !edge_found; i++){
    if((edges[i].left_map == map1 && edges[i].right_map == map2) ||
       (edges[i].left_map == map2 && edges[i].right_map == map1)){
      edge_found = true;
      return i;
    }
  }
  return -1;
}
*/
vector<int> graph_component::linear_path(int node_ind){
  assert(node_ind >= 0 && node_ind < nodes.size());
  
  vector<int> node_visit;
  for(int i=0; i<nodes.size(); i++){
    node_visit.push_back(0); //never visited before
  }

  vector<int> res;
  
  bool search_exhausted = false;
  int cur_node = node_ind;

  while(!search_exhausted){
    assert(cur_node >=0 && cur_node<nodes.size());
    assert(node_visit[cur_node] == 0);
    
    res.push_back(cur_node);
    node_visit[cur_node] = 1;

    if(nodes[cur_node].right_maps_inds.size() == 0){
      search_exhausted = true;
    }
    else{
      //select an edge with the largest t-score
      int maximal_edge = -1;
      double maximal_t_score;
      int best_next_node = -1;

      for(int i=0; i<nodes[cur_node].right_maps_inds.size(); i++){
	int map1_ind = cur_node;
	int map2_ind = nodes[cur_node].right_maps_inds[i];

	int cur_edge_ind = nodes[map1_ind].get_edge_ind(map2_ind);
	assert(cur_edge_ind != UNDEF_IND);
	//int cur_edge_ind = edge_ind(map1_ind, map2_ind);
	if(i==0){
	  maximal_edge = i;
	  maximal_t_score = edges[cur_edge_ind].t_score;
	  best_next_node = map2_ind;
	}
	else{
	  if(maximal_t_score < edges[cur_edge_ind].t_score){
	    maximal_edge = i;
	    maximal_t_score = edges[cur_edge_ind].t_score;
	    best_next_node = map2_ind;
	  }
	}
      }

      assert(maximal_edge >= 0);
      assert(best_next_node >= 0);
      cur_node = best_next_node;
    }
  }
  
  return res;
}


vector<int> graph_component::longest_branch(int node_ind, double weight){
  cout<<"processing node: "<<node_ind;
  cout<<" next_maps: "<<nodes[node_ind].right_maps.size()<<endl;
  
  for(int i=0; i<nodes[node_ind].right_maps.size(); i++){
    cout<<" "<<nodes[node_ind].right_maps_inds[i];
  }

  cout<<endl<<endl;
  double max_weight = 0;
  vector<int> heaviest_path;
 
  for(int i=0; i<nodes[node_ind].right_maps.size(); i++){
    double cur_edge_weight = nodes[node_ind].right_edge_weights[i];
    int next_node_ind = nodes[node_ind].right_maps_inds[i];

    double cur_heaviest_w = 0;
    double next_branch_w = 0;
    vector<int> cur_heaviest_path;

    cur_heaviest_path = longest_branch(next_node_ind, next_branch_w);
    cur_heaviest_w = cur_edge_weight + next_branch_w;

    bool assign = false;
    if(i==0) assign = true;    
    else{
      if(cur_heaviest_w > max_weight) assign = true;      
    }
    if(assign){
      heaviest_path = cur_heaviest_path;
      max_weight = cur_heaviest_w;
    }
  }

  
  heaviest_path.push_back(node_ind);
  weight = max_weight;
  return heaviest_path;
}

vector<int> graph_component::longest_branch(int node_ind){
  vector<int> path;

  double max_weight;
  vector<int> heaviest_path;
  
  heaviest_path = longest_branch(node_ind, max_weight);

  for(int i=heaviest_path.size()-1; i>=0; i--){
    path.push_back(heaviest_path[i]);
  }
  
  return path;
}
