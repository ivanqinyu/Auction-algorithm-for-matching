
#include <boost/config.hpp>

#include "stdafx.h"
#include "Snap.h"
#include <map>
#include <time.h>
#include <iostream>
#include <vector>

#include "gnuplot-iostream.h"

#include <fstream>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/bipartite.hpp>
#include <boost/graph/max_cardinality_matching.hpp>

#define DIRNAME "desktop/project/graphviz"

using namespace std;
using namespace boost;
using namespace TSnap;

typedef adjacency_list< vecS, vecS, undirectedS > vector_graph_t;
typedef std::pair< long, long > E;

int upper(float k){
  if (k == int(k)){
    return int(k);
  }
  else{
    return int(k)+1;
  }
}

vector<pair<int, int>> auction(vector_graph_t bipartite_vector_graph, vector<bool> bipart,float epsilon, unsigned long numSetA, unsigned long numSetB, float agressive, int maximum_matching){
  
  vector<pair<int, int>> plot_dat;
  plot_dat.push_back(make_pair(0, 0));
  
  int num_m_edges=0;
  map<int, int> a;
  
  
  bool f09 = false;
  bool f095= false;
  map<int, float> p;
  
  int ii=0,jj=0;
  for (int i=0;i<numSetA+numSetB;i++){
    if (bipart[i]==0){
      a.insert(make_pair(i, -1));
      ii++;
    }
    else{
      p.insert(make_pair(i, 0));
      jj++;
    }
  }
  float price_epsilon = epsilon;
  vector_graph_t pre_Mr((numSetA+numSetB));
//  cout<<numSetA<<","<<numSetB<<endl;
  clock_t start,end;
  start=clock();
  for (int k=1;k<upper(2/epsilon/epsilon); k++){
    if (price_epsilon<0.5 and agressive == -1){
      if (plot_dat[plot_dat.size()-1].second == plot_dat[plot_dat.size()-2].second and plot_dat[plot_dat.size()-3].second == plot_dat[plot_dat.size()-2].second){
        price_epsilon+=0.0001;
      }
    }
    else if (price_epsilon<0.5 and agressive != 0){
      if (plot_dat[plot_dat.size()-1].second == plot_dat[plot_dat.size()-2].second and plot_dat[plot_dat.size()-3].second == plot_dat[plot_dat.size()-2].second){
        price_epsilon*=agressive;
      }
      else{
        price_epsilon=epsilon;
      }
    }
//    cout<<price_epsilon<<endl;
    
    vector_graph_t Gr((numSetA+numSetB));

    map<int, int>::iterator a_i = a.begin();

    while (a_i!=a.end()){
      if (a_i->second==-1){
        float min_d = 99999;
        int s = a_i->first;
        graph_traits<vector_graph_t>::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = out_edges(s, bipartite_vector_graph); ei != ei_end; ++ei){
          int t = target(*ei, bipartite_vector_graph);
//          cout<<"s,t  in  Gr =  "<<s <<", "<<t<<"s bipart = "<<bipart[s]<<" t bipart = "<<bipart[t]<<endl;
          if (p[t]<min_d and p[t]<1){
            min_d = p[t];
          }
        }
//        cout<<"min_d"<<min_d<<endl;
        if (min_d<99999){
          graph_traits<vector_graph_t>::out_edge_iterator ei, ei_end;
          for (boost::tie(ei, ei_end) = out_edges(s, bipartite_vector_graph); ei != ei_end; ++ei){
            int t = target(*ei, bipartite_vector_graph);
            
            if (p[t]-min_d<0.000001){
              add_edge(s,t,Gr);
            }
          }
          
        }
      }
      a_i++;
    }
    if (num_edges(Gr)==0){
      end=clock();
      cout<<"iteration: "<<k<<"     number_edge: "<<num_m_edges<<" with auction algorithm.  Runtime:"<<(double) (end-start)/CLOCKS_PER_SEC << std::endl;
      return plot_dat;
    }
    vector_graph_t Mr((numSetA+numSetB));
    bool mate_Mr[numSetA+numSetB];
    for (int i=0;i<numSetA+numSetB;i++){
      mate_Mr[i]=false;
    }
    
    graph_traits<vector_graph_t>:: vertex_iterator vi,vi_end;
    for (boost::tie(vi, vi_end) = vertices(Gr); vi != vi_end; ++vi){
      int s = *vi;
      if (bipart[s]==1){
        continue;
      }
      graph_traits<vector_graph_t>::out_edge_iterator ei, ei_end;
      
      vector<int> t_vec;
      t_vec.clear();
      for (boost::tie(ei, ei_end) = out_edges(*vi, Gr); ei != ei_end; ++ei){
        int t = target(*ei, Gr);
        t_vec.push_back(t);
      }
      
     for (int kk=t_vec.size()-1;kk>=0;kk--){
        int t=t_vec[kk];
        if (mate_Mr[s]==false and mate_Mr[t]==false and s!=t){
          mate_Mr[s]=true;
          mate_Mr[t]=true;
          add_edge(s,t,Mr);
          
        }
      }
    }
    
    graph_traits<vector_graph_t>::edge_iterator ei, ei_end;
    for (boost::tie(ei, ei_end) = edges(Mr); ei != ei_end; ++ei){
      int s = source(*ei, Mr);
      int t = target(*ei, Mr);
      
      if (boost::degree(t, pre_Mr)!=0 ){
        graph_traits<vector_graph_t>::out_edge_iterator eei, eei_end;
        boost::tie (eei,eei_end) = out_edges(t,pre_Mr);
        
        int pre_val = target(*eei ,pre_Mr);

        if (pre_val != s){
          remove_edge(t,pre_val,pre_Mr);
          num_m_edges--;
          a[pre_val]=-1;
        }
      }
      
      a[s]=t;
      p[t]+=price_epsilon;
//
//      cout<<t<<"."<<p[t]<<",";
      add_edge(s,t,pre_Mr);
      num_m_edges++;
      
      
      
      
      if (num_m_edges > maximum_matching){
//        cout<<"price"<<endl;
//        for (int i=numSetA;i<numSetA+numSetB;i++){
//          cout<<p[i]<<",";
//        }
        end=clock();
        cout<<"iteration: "<<k<<"     number_edge: "<<num_m_edges<<" with auction algorithm.  Runtime:"<<(double) (end-start)/CLOCKS_PER_SEC << std::endl;
        
        plot_dat.push_back(make_pair(k, num_m_edges));
        return plot_dat;
      }
      if (num_m_edges > 0.95*maximum_matching and f095==false){
        f095=true;
        cout<<"iteration: "<<k<<"     number_edge: "<<num_m_edges<< std::endl;
      }
      else if (num_m_edges > 0.9*maximum_matching and f09==false){
        f09=true;
        cout<<"iteration: "<<k<<"     number_edge: "<<num_m_edges << std::endl;
      }
    }

//    cout<<"iteration: "<<k<<endl;
//    for (boost::tie(ei, ei_end) = edges(pre_Mr); ei != ei_end; ++ei){
//      int s = source(*ei, Mr);
//      int t = target(*ei, Mr);
//      cout<<"\\draw(f"<<s<<")--(s"<<t<<");"<<endl;
//    }
//    cout<<k<<","<<num_m_edges<<endl;
    plot_dat.push_back(make_pair(k, num_m_edges));
  }
  
  cout<<"M has "<<num_m_edges<<" edges with auction algorithm"<<endl;
  return plot_dat;
}

#define SEMICOMPLETE 0;
#define RANDOM_GRAPH 1;
#define FACEBOOK_DATA 2;
#define Google_DATA 3;

int main(int argc, char* argv[]) {
  typedef PUNGraph PGraph; // undirected graph
  printf("Creating graph:\n");
  PGraph G = PGraph::TObj::New();
  unsigned long numSetA=0,numSetB=0;
  vector<bool> bipart;
  
  // what type of graph do you want to use?
  int selection = FACEBOOK_DATA;

  if (selection == 0){
//  semicomplete DATA
    for (int n = 0; n < 2048; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...,9
      bipart.push_back(false);
    }
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      if (NI.GetId()<(G->GetNodes()/2)){
        bipart[NI.GetId()]=0;
        numSetA+=1;
      }
      else{
        bipart[NI.GetId()]=1;
        numSetB+=1;
      }
    }
  
    for (int i=0;i<numSetA;i++){
      if (i>=numSetB){
        break;
      }
      for (int j=numSetB-1;j>=i;j--){
        G->AddEdge(i, numSetA+j);
      }
    }
  }
  else if (selection == 1){
//    RANDOM GENERATE DATASET
    G = PGraph::TObj::New();
    for (int n = 0; n < 500; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...,9
      bipart.push_back(false);
    }
    for (int e = 0; e < 10000; e++) {
      const int NId1 = G->GetRndNId();
      const int NId2 = G->GetRndNId();
      if (NId1<G->GetNodes()/2 and NId2>=G->GetNodes()/2 and NId1 != NId2 and G->AddEdge(NId1, NId2) != -2) {
        continue;
//        printf("  Edge %d -- %d added\n", NId1,  NId2);
      }
      else {
        continue;
//        printf("  Edge %d -- %d already exists\n", NId1, NId2);
      }
    }
    IAssert(G->IsOk());
    //   delete
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
//      bipart[NI.GetId()] = int(rand()%2);
//      cout<<bipart[NI.GetId()]<<endl;
      if (NI.GetId()<G->GetNodes()/2){
        bipart[NI.GetId()]=0;
        numSetA+=1;
      }
      else{
        bipart[NI.GetId()]=1;
        numSetB+=1;
      }
    }
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      for (int e = 0; e < NI.GetDeg(); e++) {
//        cout<<NI.GetId()<<NI.GetNbrNId(e)<<"\n";
        if (bipart[NI.GetId()] == bipart[NI.GetNbrNId(e)] ) {
          G->DelEdge(NI.GetId(), NI.GetNbrNId(e));
  //        printf("Delete edge %d -- %d\n", NI.GetId(), NI.GetNbrNId(e));
          e--;
        }
      }
    }

  }
  else if (selection == 2){
    //FACEBOOK DATASET
    G = LoadEdgeList<PGraph>(TStr::Fmt("%s/matching_def.txt", DIRNAME));

    printf("Loaded Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
    
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      bipart.push_back(0);
    }
    //   delete
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      
      bipart[NI.GetId()] = int(rand()%2);
      if (bipart[NI.GetId()]==0){
        numSetA+=1;
      }
      else{
        numSetB+=1;
      }
    }
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      for (int e = 0; e < NI.GetDeg(); e++) {
  //      cout<<NI.GetId()<<NI.GetNbrNId(e)<<"\n";
        if (bipart[NI.GetId()] == bipart[NI.GetNbrNId(e)] ) {
          G->DelEdge(NI.GetId(), NI.GetNbrNId(e));
  //        printf("Delete edge %d -- %d\n", NI.GetId(), NI.GetNbrNId(e));
          e--;
        }
      }
    }

  }
  else if (selection == 3){
    //GooglePlus DATASET
    G = LoadEdgeList<PGraph>(TStr::Fmt("%s/transed_google.txt", DIRNAME));
    printf("Loaded Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
    
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      bipart.push_back(0);
    }
    //   delete
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      
      bipart[NI.GetId()] = int(rand()%2);
      if (bipart[NI.GetId()]==0){
        numSetA+=1;
      }
      else{
        numSetB+=1;
      }
    }
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      for (int e = 0; e < NI.GetDeg(); e++) {
  //      cout<<NI.GetId()<<NI.GetNbrNId(e)<<"\n";
        if (bipart[NI.GetId()] == bipart[NI.GetNbrNId(e)] ) {
          G->DelEdge(NI.GetId(), NI.GetNbrNId(e));
  //        printf("Delete edge %d -- %d\n", NI.GetId(), NI.GetNbrNId(e));
          e--;
        }
      }
    }
  }
  
  // transfer from SNAP to Boost Library
  
  vector<E> bipartite_edges;
  long eds = 0;
  for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    for (int e = 0; e < NI.GetDeg(); e++) {
      
      if (NI.GetId() < NI.GetNbrNId(e)){
        bipartite_edges.push_back( make_pair(NI.GetId(), NI.GetNbrNId(e)));
        eds++;
      }
    }
  }
  vector_graph_t bipartite_vector_graph;
  cout<<"edges:"<<eds<<endl;
  for (int i = 0; i < eds; i++){
    add_edge(bipartite_edges[i].first, bipartite_edges[i].second, bipartite_vector_graph);
    
    cout<<"\\draw(f"<<bipartite_edges[i].first<<")--(s"<<bipartite_edges[i].second<<");"<<endl;
  }
  if (is_bipartite(bipartite_vector_graph) == true){
    cout<<"It is a Bipartite";
  }
  else{
    cout<<"It is not a Bipartite";
    return 0;
  }
  
  cout<<endl;
  cout<<endl;
  cout<<"----------------------------"<<endl;
  cout<<"----------------------------"<<endl;
  cout<<"----------Matching----------"<<endl;
  cout<<"----------------------------"<<endl;
  cout<<"----------------------------"<<endl;
  cout<<endl;
  cout<<endl;
  
  
  clock_t start,end;
  start=clock();
  // maximum matching from boost library
  vector< graph_traits< vector_graph_t >::vertex_descriptor > mate(G->GetNodes());
  edmonds_maximum_cardinality_matching(bipartite_vector_graph, &mate[0]);
  int maximum_matiching = matching_size(bipartite_vector_graph, &mate[0]);
  end=clock();
  cout<<"M has "<< maximum_matiching <<" edges with build-in algorithm.  Runtime:"<<(double) (end-start)/CLOCKS_PER_SEC << std::endl;

  std::cout << "The matching is:" << std::endl;
  graph_traits< vector_graph_t >::vertex_iterator vi, vi_end;
  for (boost::tie(vi, vi_end) = vertices(bipartite_vector_graph); vi != vi_end; ++vi)
    if (mate[*vi] != graph_traits< vector_graph_t >::null_vertex() && *vi < mate[*vi])
      std::cout << "{" << *vi << ", " << mate[*vi] << "}" << std::endl;
  

//
  //maximum auction algorithm

  Gnuplot gp;
  gp<<"set term post color\n";
  gp<<"set output 'a.ps'\n";

  vector<pair<int, int>> plot_dat;

  gp<<"set xlabel 'Iterations'\n";
  gp<<"set ylabel 'Number of Matching'\n";
  gp<<"set title 'Get matching with different epsilon'\n";
  gp << "set yrange [0:500]\n";
  gp << "set xrange [0:500]\n";

//  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 0,maximum_matiching);
//  gp <<"plot "<<maximum_matiching<<"title 'maximum',"<<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'constant epsilon 0.001'";

  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 0,maximum_matiching);
  gp <<"plot "<<maximum_matiching<<"title 'maximum',"<<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'constant epsilon 0.001'";

//  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 0,maximum_matiching);
//  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'constant epsilon 0.001',";
//
//  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, -1,maximum_matiching);
//  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon +0.0001',";
//
//  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 1.5,maximum_matiching);
//  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon *1.1',";
//
//  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 1.001,maximum_matiching);
//  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon *1.001',";
//
//  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 1.0001,maximum_matiching);
//  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon *1.0001'";

//  gp <<endl;
  return 0;
}

