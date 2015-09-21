#define _EDGE_H

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

class edge{
 public:
  double weight;

  string left_map;
  string right_map;
  
  hash_key left_map_hash_key;
  hash_key right_map_hash_key;

  short left_orient;
  short right_orient;

  int left_map_size;
  int right_map_size;
 
  float s_score;
  float t_score;

  vector<site_ind> left_al_sites;
  vector<site_ind> right_al_sites;

  bool containment;

  short rev(short _or);
  edge();
  edge(string& _left_map, string& _right_map, 
       hash_key _left_map_hash_key, hash_key _right_map_hash_key,
       short _left_orient, short _right_orient,
       int _left_map_size, int _right_map_size,
       float _s_score, float _t_score, 
       vector<site_ind> _left_al_sites, vector<site_ind> _right_al_sites,
       bool _containment);

  void print();
 
  edge& operator=(const edge& e);
  edge inverse();
};
