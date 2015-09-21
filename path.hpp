#define _PATH_H

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

class path{
public:
  vector<int> map_inds;
  vector<int> edge_inds;
  vector<double> edge_weights;
private:
  double length;
public:
  path(){};
  void set_length();
  double get_length();
  double get_abs_length();
};

void path::set_length(){
  length = 0;
  for(int i=0; i<edge_weights.size(); i++){
    length += edge_weights[i];
  }
}

double path::get_length(){
  return length;
}

double path::get_abs_length(){
  return fabs(get_length());
}
