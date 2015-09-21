#include "graph.h"

void graph::clear(){
  components.clear();
}
graph::graph(hash_key _max_hash_key){
  max_hash_key = _max_hash_key;
}

void graph::construct_graph(list<int>& good_edge_inds,
			    const vector<edge>& _edges){//, 
			    //const list<int>& good_edge_inds,
			    //vector<node>& _nodes){
  graph_component starter(max_hash_key);
  
  //starter.edges = _edges;
  //starter.nodes = _nodes;

  /*
  vector<int> active_edges;
  for(int i=0; i<starter.edges.size(); i++){
    active_edges.push_back(i);
  }
  */

  vector<vector< hash_key > > map_hash_keys;
  vector<vector< orientation > > map_orients;
  vector<vector< int > > new_edges;
  list<int> consistent_edge_set = 
    starter.consistent_edge_inds(good_edge_inds, _edges, map_hash_keys, 
				 map_orients, new_edges);
  //starter.make_hash_table();

  for(int i=0; i<map_hash_keys.size(); i++){
    cerr<<"construct_graph: making "<<i;
    cerr<<" component of "<<map_hash_keys.size();
    cerr<<" nodes: "<<map_hash_keys[i].size();
    cerr<<" edges: "<<new_edges[i].size()<<endl;
    assert(map_hash_keys[i].size() == map_orients[i].size());
    
    graph_component new_comp(max_hash_key);
    assert(map_hash_keys[i].size() == map_orients[i].size());
    for(int j=0; j<map_hash_keys[i].size(); j++){
      hash_key cur_map_hk = map_hash_keys[i][j];
      //cerr<<"cur_map_hk = "<<cur_map_hk<<endl;
      //string cur_map_name = return_map_name(cur_map_hash_key);
      //int cur_node_index = starter.get_node_index(cur_map_name);
      //int cur_node_index = starter.get_node_index(cur_map_hash_key);
      //assert(cur_node_index != UNDEF_IND);

      /*
      node new_node(starter.nodes[cur_node_index].hk,
		    starter.nodes[cur_node_index].read_name, 
		    map_orients[i][j], 
		    starter.nodes[cur_node_index].map_read);
      new_node.id = starter.nodes[cur_node_index].id;
      */
      vector<fr_size> cur_map_read = return_map_read(cur_map_hk);
      node new_node(cur_map_hk, return_map_name(cur_map_hk),
		    map_orients[i][j], cur_map_read);
      new_comp.nodes.push_back(new_node);
    }
    
    for(int j=0; j<new_edges[i].size(); j++){
      int cur_edge_ind = new_edges[i][j];
      //cout<<"adding an edge:"<<cur_edge_ind<<endl;
      //assert(cur_edge_ind >= 0 && cur_edge_ind<starter.edges.size());
      assert(cur_edge_ind >= 0 && cur_edge_ind < _edges.size());
      new_comp.edges.push_back(_edges[cur_edge_ind]);    
    }
    cout<<"filling out the component structure."<<endl;
    new_comp.fill_the_component_data();
    cout<<"structure completed"<<endl;
    if(new_comp.nodes.size() > 0) components.push_back(new_comp);
  }
  //starter.remove_hash_table();
}

void graph::init_data(){
  for(int i=0; i<components.size(); i++){
    components[i].clean();
    components[i].fill_the_component_data();
    //components[i].index_nodes();
    //components[i].assign_edge_weights();
  }
}
void graph::assign_weights(){
  for(int i=0; i<components.size(); i++){
    components[i].assign_edge_weights();
  }
}
void graph::output_graph(const char* output_file){
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

  for(int i=0; i<components.size(); i++){
    if(components[i].edges.size() >= 10){
      ostringstream comp_prefix;
      comp_prefix<<i<<":";
      //string prefix = comp_prefix.c_str();
      //components[i].index_nodes();
      components[i].output_edges(out_str, comp_prefix);
    }
  }

  out_str<<"}"<<endl;
  out_str.close();
}

/*
void graph::add_edge(edge& e, vector<fr_size>& left_map_read,
		     vector<fr_size>& right_map_read){

  bool left_map_found = false;
  bool right_map_found = false;

  int left_comp_ind = -1;
  int right_comp_ind = -1;

  int left_read_ind = -1;
  int right_read_ind = -1;

  for(int i=0; i<components.size(); i++){
    //first find which components the maps belong to
    for(int j=0; j<components[i].nodes.size(); j++){
      if(e.left_map == components[i].nodes[j].read_name){
	left_map_found = true;
	left_comp_ind = i;
	left_read_ind = j;
      }
      if(e.right_map == components[i].nodes[j].read_name){
	right_map_found = true;
	right_comp_ind = i;
	right_read_ind = j;
      }
    }
  }

  if(!left_map_found && !right_map_found){
    //make a new component
    cout<<"making a new component"<<endl;
    graph_component new_component(max_hash_key);
    new_component.add_edge(e, left_map_read, right_map_read);
    components.push_back(new_component);
    return;
  }

  if(left_map_found && !right_map_found){
    assert(left_comp_ind >= 0 && left_comp_ind < components.size());
    cout<<"updating existing component with ";
    cout<<e.right_map<<endl;
    components[left_comp_ind].add_edge(e, left_map_read, right_map_read);
    return;
  }
  
  if(!left_map_found && right_map_found){
    assert(right_comp_ind >= 0 && right_comp_ind < components.size());
    cout<<"updating existing component with ";
    cout<<e.left_map<<endl;
    components[right_comp_ind].add_edge(e, left_map_read, right_map_read);
    return;
  }

  assert(left_map_found && right_map_found);
  
  if(right_comp_ind != left_comp_ind){
    cout<<"merging components "<<left_comp_ind;
    cout<<" and "<<right_comp_ind<<endl;
    assert(left_comp_ind != -1 && right_comp_ind != -1);
    //merge two components
    
    graph_component new_comp(max_hash_key);

    for(int i=0; i<components[left_comp_ind].nodes.size(); i++){
      new_comp.nodes.push_back(components[left_comp_ind].nodes[i]);
    }
    if(e.left_orient == components[left_comp_ind].nodes[left_read_ind].orient){
      if(e.right_orient == components[right_comp_ind].nodes[right_read_ind].orient){
	//both maps in the same orientation as in the overlap
	for(int i=0; i<components[right_comp_ind].nodes.size(); i++){
	  new_comp.nodes.push_back(components[right_comp_ind].nodes[i]);
	}
      }
      else{
	//left map in the same orient, right map in the opposite
	cout<<"reversing required 1"<<endl;
	for(int i=0; i<components[right_comp_ind].nodes.size(); i++){
	  new_comp.nodes.push_back
	    (components[right_comp_ind].nodes[i].inverse());
	}
      }
    }
    else{
      if(e.right_orient == components[right_comp_ind].nodes[right_read_ind].orient){ 
	cout<<"reversing required 2"<<endl;
	for(int i=0; i<components[right_comp_ind].nodes.size(); i++){
	  new_comp.nodes.push_back
	    (components[right_comp_ind].nodes[i].inverse());
	}
      }
      else{
	cout<<"both in reverse orient, no reversing required"<<endl;
	for(int i=0; i<components[right_comp_ind].nodes.size(); i++){
	  new_comp.nodes.push_back(components[right_comp_ind].nodes[i]);
	}
      }      
    }

    for(int i=0; i<components[left_comp_ind].edges.size(); i++){
      new_comp.edges.push_back(components[left_comp_ind].edges[i]);
    }
    for(int i=0; i<components[right_comp_ind].edges.size(); i++){
      new_comp.edges.push_back(components[right_comp_ind].edges[i]);
    }
    new_comp.add_edge(e, left_map_read, right_map_read);
    
    if(left_comp_ind < right_comp_ind){
      vector<graph_component>::iterator it;
      it = components.begin() + right_comp_ind;
      components.erase(it);
      
      it = components.begin() + left_comp_ind;
      components.erase(it);
    }
    else{
      vector<graph_component>::iterator it;
      it = components.begin() + left_comp_ind;
      components.erase(it);
      
      it = components.begin() + right_comp_ind;
      components.erase(it);
    }

    components.push_back(new_comp);
    return;
  }
  else{
    //both are in the same component
    assert(left_comp_ind != -1 && right_comp_ind != -1);
    components[left_comp_ind].add_edge(e, left_map_read, right_map_read);
    return;
  }
  assert(false);
}
*/

void graph::print(){
  cout<<"--------------"<<endl;
  cout<<"printing the graph data:"<<endl;
  cout<<"components: "<<components.size()<<endl;

  for(int i=0; i<components.size(); i++){
    cout<<endl;
    cout<<"+++++"<<endl;
    cout<<"component: "<<i<<endl;
    components[i].print();
    cout<<"+++++"<<endl;
  }

  cout<<"--------------"<<endl<<endl;
    
}
