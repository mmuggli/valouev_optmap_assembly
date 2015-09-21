#define _GRAPH_H

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

//#define hash_key long
#define orientation short
#define forward 1
#define reverse 0
#define UNDEF_ORIENT -1000
#define UNDEF_IND -1
#define UNDEF_HASH_KEY -1
#define UNDEF_MAP_NAME "undefined_map_name"

#ifndef _NODE_H
#include "node.cpp"
#endif

#ifndef _EDGE_H
#include "edge.cpp"
#endif

#ifndef _GRAPH_COMPONENT_H
#include "graph_component.cpp"
#endif

class graph{
 public:
  vector<graph_component> components;
  hash_key max_hash_key;
  graph(hash_key _max_hash_key);

  /*
  void add_edge(edge& e, vector<fr_size>& left_map_read, 
		vector<fr_size>& right_map_read);
  */

  void clear();
  void print();
  void output_graph(const char* output_file);
  void assign_weights();
  void init_data();
  //this function assigns weights to edges
  //and nodes based on local paths

  void construct_graph(list<int>& good_edge_inds, 
		       const vector<edge>& _edges);

  void mark_chimeras();
  //void construct_graph(const vector<edge>& _edges, vector<node>& nodes);
};
