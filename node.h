#define _NODE_H

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

//using namespace std;

class node{
 public:
  string read_name;
  orientation orient;
  vector<fr_size> map_read;

  int _index;
  int id;
  double weight;
  hash_key hk;
  
  //vector<string> left_maps;
  ///vector<string> right_maps;
  //vector<string> contained_maps;
  //vector<string> container_maps;

  node();
  node(hash_key _hk, string _read_name, orientation _orient, 
       vector<fr_size>&_map_read);
  
  double distance_to_center(int from_site, orientation _orient);
  double region_size(int start, int end, orientation _orient);
  int map_center(orientation _orient);
  vector<double> region(int start, int end, orientation _orient);
  //returns the fragment index within which the center
  //of the map is located

  void print();
  node inverse();
  int get_edge_ind(int to_node_ind);
  double size();

  bool structures_ok();

  vector<string> right_maps;
  vector<int> right_maps_inds;
  vector<int> right_edges_inds;
  vector<double> right_edge_weights;

  vector<string> left_maps;
  vector<int> left_maps_inds;
  vector<int> left_edges_inds;
  vector<double> left_edge_weights;

  vector<string> contained_maps;
  vector<int> contained_maps_inds;
  vector<int> contained_edges_inds;
  vector<double> contained_edge_weights;

  vector<string> container_maps;
  vector<int> container_maps_inds;
  vector<int> container_edges_inds;
  vector<double> container_edge_weights;
};
