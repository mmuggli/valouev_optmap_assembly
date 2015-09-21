#include "edge.h"

bool valid_orient(orientation _or){
  if(_or == forward || _or == reverse) return true;
  else return false;
}

orientation edge::rev(orientation _or){
  assert(_or == forward || _or == reverse);
  if(_or == forward) return reverse;
  else return forward;
}

edge::edge(){}

edge::edge(string& _left_map, string& _right_map, 
	   hash_key _left_map_hash_key, hash_key _right_map_hash_key,
	   orientation _left_orient, orientation _right_orient,
	   int _left_map_size, int _right_map_size,
	   float _s_score, float _t_score, 
	   vector<site_ind> _left_al_sites, vector<site_ind> _right_al_sites,
	   bool _containment){
  assert(_left_al_sites.size() == _right_al_sites.size());
  
  left_map = _left_map;
  right_map = _right_map;

  left_map_hash_key = _left_map_hash_key;
  right_map_hash_key = _right_map_hash_key;

  if(!(left_orient == forward || left_orient == reverse)){
    cerr<<_left_orient<<endl;
  }

  assert(_left_orient == forward || _left_orient == reverse);
  assert(_right_orient == forward || _right_orient == reverse);

  left_orient = _left_orient;
  right_orient = _right_orient;

  left_map_size = _left_map_size;
  right_map_size = _right_map_size;
  
  s_score = _s_score;
  t_score = _t_score;

  left_al_sites = _left_al_sites;
  right_al_sites = _right_al_sites;

  containment = _containment;
}

void edge::print(){
  cout<<"edge: ";
  if(containment) cout<<"(containment)";
  cout<<endl;
  cout<<left_map<<" ("<<left_map_size<<")  -> ";
  cout<<right_map<<" ("<<right_map_size<<")"<<endl;
  cout<<left_orient<<" "<<right_orient<<" ";
  cout<<"["<<left_al_sites[0]<<", ";
  cout<<left_al_sites[left_al_sites.size()-1]<<"] -> [";
  cout<<right_al_sites[0]<<", ";
  cout<<right_al_sites[right_al_sites.size()-1]<<"] ";
  //cout<<endl;
  cout<<"s_score: "<<s_score<<" t_score: "<<t_score;
  cout<<" ls:"<<left_map_size<<" rs:"<<right_map_size<<endl;
}

edge& edge::operator=(const edge& e){
  left_map = e.left_map;
  right_map = e.right_map;

  left_map_hash_key = e.left_map_hash_key;
  right_map_hash_key = e.right_map_hash_key;

  left_orient = e.left_orient;
  right_orient = e.right_orient;

  left_map_size = e.left_map_size;
  right_map_size = e.right_map_size;
  
  s_score = e.s_score;
  t_score = e.t_score;

  left_al_sites = e.left_al_sites;
  right_al_sites = e.right_al_sites;

  containment = e.containment;

  weight = e.weight;
  
  return *this;
}

edge edge::inverse(){
  assert(left_al_sites.size() == right_al_sites.size());

  vector<site_ind> inv_left_al_sites;
  vector<site_ind> inv_right_al_sites;

  for(int i=0; i<left_al_sites.size(); i++){
    int new_left_site = left_map_size - left_al_sites[i];
    int new_right_site = right_map_size - right_al_sites[i];

    assert(new_left_site >= 0);
    assert(new_right_site >= 0);

    inv_left_al_sites.insert(inv_left_al_sites.begin(), new_left_site);
    inv_right_al_sites.insert(inv_right_al_sites.begin(), new_right_site);
  }

  if(!containment){
    edge new_edge(right_map, left_map, 
		  right_map_hash_key, left_map_hash_key,
		  rev(right_orient), rev(left_orient),
		  right_map_size, left_map_size, s_score, t_score,
		  inv_right_al_sites, inv_left_al_sites, containment);
    new_edge.weight = weight;
    return new_edge;
  }
  else{
    edge new_edge(left_map, right_map, 
		  left_map_hash_key, right_map_hash_key,
		  rev(left_orient), rev(right_orient),
		  left_map_size, right_map_size, s_score, t_score,
		  inv_left_al_sites, inv_right_al_sites, containment);
    new_edge.weight = weight;
    return new_edge;
  }
}
