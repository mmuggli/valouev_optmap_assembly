#include "node.h"

double node::size(){
  double res = 0;
  for(int i=0; i<map_read.size(); i++){
    res += map_read[i];
  }
  return res;
}

bool node::structures_ok(){
  if(right_maps.size() != right_maps_inds.size()) return false;
  if(right_maps.size() != right_edges_inds.size()) return false;
  if(left_maps.size() != left_maps_inds.size()) return false;
  if(left_maps.size() != left_edges_inds.size()) return false;
  if(contained_maps.size() != contained_maps_inds.size()) return false;
  if(contained_maps.size() != contained_edges_inds.size()) return false;
  if(container_maps.size() != container_maps_inds.size()) return false;
  if(container_maps.size() != container_edges_inds.size()) return false;    
  return true;
}

int node::get_edge_ind(int to_node_ind){
  assert(right_maps_inds.size() == right_edges_inds.size());
  assert(left_maps_inds.size() == left_edges_inds.size());
  assert(contained_maps_inds.size() == contained_edges_inds.size());
  assert(container_maps_inds.size() == container_edges_inds.size());

  for(int i=0; i<left_maps_inds.size(); i++){
    if(left_maps_inds[i] == to_node_ind) return left_edges_inds[i];
  }
  for(int i=0; i<right_maps_inds.size(); i++){
    if(right_maps_inds[i] == to_node_ind) return right_edges_inds[i];
  }
  for(int i=0; i<contained_maps_inds.size(); i++){
    if(contained_maps_inds[i] == to_node_ind) return contained_edges_inds[i];
  }
  for(int i=0; i<container_maps_inds.size(); i++){
    if(container_maps_inds[i] == to_node_ind) return container_edges_inds[i];
  }
  return UNDEF_IND;
}

double node::distance_to_center(int from_site, orientation _orient){
  assert(from_site <= map_read.size());
  int cur_center = map_center(_orient);
  assert(cur_center < map_read.size());

  double distance = 0;
  if(from_site <= cur_center){
    for(int i=from_site; i<cur_center; i++){
      if(_orient == forward){
	assert(i>=0 && i<map_read.size());		
	distance += map_read[i];
      }
      else{
	assert(map_read.size()-i-1>=0 && 
	       map_read.size()-i-1<map_read.size());		
	distance += map_read[map_read.size()-i-1];
      }
    }
    assert(cur_center>=0 && cur_center<map_read.size());
    if(_orient == forward) distance += 0.5*map_read[cur_center];
    else distance += 0.5*map_read[map_read.size()-cur_center-1];
    return distance;
  }
  else{
    for(int i=cur_center+1; i<from_site; i++){
      if(_orient == forward){
	assert(i>=0 && i<map_read.size());
	distance -= map_read[i];
      }
      else{
	assert(map_read.size()-i-1>=0 && 
	       map_read.size()-i-1<map_read.size());
	distance -= map_read[map_read.size()-i-1];
      }
    }
    assert(cur_center>=0 && cur_center<map_read.size());
    if(_orient == forward) distance -= 0.5*map_read[cur_center];
    else distance -= 0.5*map_read[map_read.size()-cur_center-1];
    return distance;
  }
}

int node::map_center(orientation _orient){
  //cout<<"map_center: current map has "<<map_read.size()<<" frags"<<endl;
  double total_map_size = region_size(0,map_read.size(),orient);
  double _center = total_map_size/2;
  
  //cout<<"total size:"<<total_map_size<<endl;
  //now track to the middle of the map
  bool center_found = false;

  double cur_dist = 0;
  double prev_dist = 0;
  int cur_fr_ind = 0;
  while(!center_found){  
    //cout<<cur_fr_ind<<" :cur_dist:"<<cur_dist<<endl;
    if(_orient == forward)
      cur_dist += map_read[cur_fr_ind];
    else
      cur_dist += map_read[map_read.size()-cur_fr_ind-1];
    //if(_center >= prev_dist && _center < cur_dist){
    if(_center < cur_dist){
      return cur_fr_ind;
      center_found = true;
    }
    prev_dist = cur_dist;
    cur_fr_ind++;
  }
  assert(false); //better not be here
}

vector<double> node::region(int start, int end, orientation _orient){
  vector<double> chunk;
  assert(start<=end);
  assert(start>=0 && end<=map_read.size());

  double cur_size = 0;
  for(int i=start; i<end; i++){
    if(_orient == forward){
      cur_size += map_read[i];
      chunk.push_back(map_read[i]);
    }
    else{
      cur_size += map_read[map_read.size()-i-1];
      chunk.push_back(map_read[map_read.size()-i-1]);
    }
  }
  return chunk;
}

double node::region_size(int start, int end, orientation _orient){
  assert(start<=end);
  assert(start>=0 && end<=map_read.size());

  double cur_size = 0;
  for(int i=start; i<end; i++){
    if(_orient == forward){
      cur_size += map_read[i];
    }
    else{
      cur_size += map_read[map_read.size()-i-1];
    }
  }
  return cur_size;
}

node::node(){
  hk = UNDEF_HASH_KEY;
}
node::node(hash_key _hk, string _read_name, orientation _orient, 
	   vector<fr_size>& _map_read){
  hk = _hk;
  read_name = _read_name;
  map_read = _map_read;
  orient = _orient;
}

void node::print(){
  cout<<"node: ";//<<endl;
  cout<<read_name<<" or: "<<orient<<endl;

  cout<<endl;
  cout<<"left maps: "<<endl;
  for(int i=0; i<left_maps.size(); i++){
    cout<<left_maps[i]<<endl;
  }
  
  cout<<endl;
  cout<<"right_maps: "<<endl;
  for(int i=0; i<right_maps.size(); i++){
    cout<<right_maps[i]<<endl;
  }
}

node node::inverse(){
  short new_orient;
  if(orient == 1) new_orient = 0;
  else new_orient = 1;

  node new_node(hk, read_name, new_orient, map_read);

  /*
  new_node.left_maps = right_maps;
  new_node.left_maps_inds = left_maps_inds;
  new_node.right_maps = left_maps;
  new_node.contained_maps = contained_maps;
  new_node.container_maps = container_maps;
  */

  return new_node;
}
