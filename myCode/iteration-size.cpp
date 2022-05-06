

#include <boost/config.hpp>

#include "stdafx.h"
#include "Snap.h"
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
typedef std::pair< int, int > E;

int upper(float k){
  if (k == int(k)){
    return int(k);
  }
  else{
    return int(k)+1;
  }
}

pair<int, int> auction(vector_graph_t bipartite_vector_graph, bool bipart[],float epsilon, int numSetA, int numSetB){
  
  pair<int, int> plot_dat;
  int num_m_edges=0;
  int a[numSetA];
  float p[numSetB];
  int index_a[numSetA],index_b[numSetB];

  for (int i=0;i<numSetA;i++){
    a[i]=-1;
  }
  for (int i=0;i<numSetB;i++){
    p[i]=0;
  }
  int ii=0,jj=0;
  for (int i=0;i<numSetA+numSetB;i++){
    if (bipart[i]==true){
      index_a[ii]=i;
      ii++;
    }
    else{
      index_b[jj]=i;
      jj++;
    }
  }

  vector_graph_t pre_Mr((numSetA+numSetB));
  for (int k=1;k<upper(2/epsilon/epsilon); k++){
//    getchar();
    vector_graph_t Gr((numSetA+numSetB));

    for (int i=0;i<numSetA;i++){

      if (a[i]==-1){
        float min_d = 99999;

        int s = index_a[i];
        graph_traits<vector_graph_t>::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = out_edges(index_a[i], bipartite_vector_graph); ei != ei_end; ++ei){
          int t = target(*ei, bipartite_vector_graph);
          int real_t = -1;
          for (int kk=0;kk<numSetB;kk++){
            if (index_b[kk]==t){
              real_t=kk;
              break;
            }
          }
          if (real_t!=-1 and p[real_t]<min_d and p[real_t]<1){
            min_d = p[real_t];
          }
        }
        
        if (min_d<99999){
          for (int kk=0;kk<numSetB;kk++){
            if (p[kk]==min_d and edge(s,index_b[kk],bipartite_vector_graph).second == true){
              add_edge(s,index_b[kk],Gr);
//              cout<<"Gr: "<<s<<" "<<index_b[kk]<<"  price ="<<min_d<<endl;
            }
          }
          
        }
      }
    }
    
//    cout<<"Gr: "<<num_edges(Gr)<<endl;
    vector_graph_t Mr((numSetA+numSetB));
    bool mate_Mr[numSetA+numSetB];
    for (int i=0;i<numSetA+numSetB;i++){
      mate_Mr[i]=false;
    }
    
    graph_traits<vector_graph_t>:: vertex_iterator vi,vi_end;
    for (boost::tie(vi, vi_end) = vertices(Gr); vi != vi_end; ++vi){
      graph_traits<vector_graph_t>::out_edge_iterator ei, ei_end;
      
      int s = *vi;
      vector<int> t_vec;
      t_vec.clear();
      for (boost::tie(ei, ei_end) = out_edges(*vi, Gr); ei != ei_end; ++ei){
        int t = target(*ei, Gr);
        t_vec.push_back(t);
      }
      
//      for (int kk=t_vec.size()-1;kk>=0;kk--){
//        cout<<" "<<s<<"--"<<t_vec[kk]<<" ";
//      }
//      cout<<endl;
      for (int kk=t_vec.size()-1;kk>=0;kk--){
        int t=t_vec[kk];
        if (mate_Mr[s]==false and mate_Mr[t]==false and s!=t){
          mate_Mr[s]=true;
          mate_Mr[t]=true;
//          cout<<"Mr: "<<s<<" "<<t<<endl;
//
          add_edge(s,t,Mr);
        }
      }
    }
    
//    cout<<"ADD Mr: "<<num_edges(Mr)<<endl;
    
    
    graph_traits<vector_graph_t>::edge_iterator ei, ei_end;
    
    for (boost::tie(ei, ei_end) = edges(Mr); ei != ei_end; ++ei){
      int s = source(*ei, Mr);
      int t = target(*ei, Mr);
      
      int real_s=-1;
      for (int kk=0;kk<numSetA;kk++){
        if (index_a[kk]==s){
          real_s=kk;
          break;
        }
      }
      a[real_s]=t;
//      cout<<"a["<<real_s<<"]"<<s<<"  = "<<t<<endl;
 
      if (boost::degree(t, pre_Mr)!=0){
        graph_traits<vector_graph_t>::out_edge_iterator eei, eei_end;
        boost::tie (eei,eei_end) = out_edges(t,pre_Mr);
        int pre_val = target(*eei ,pre_Mr);
//        cout<<boost::degree(t, pre_Mr)<<","<<t<<","<<pre_val<<endl;
        
        remove_edge(t,pre_val,pre_Mr);
        num_m_edges--;
        
        int real_pre=-1;
        for (int kk=0;kk<numSetA;kk++){
          if (index_a[kk]==pre_val){
            real_pre=kk;
            break;
          }
        }
        a[real_pre]=-1;
      }
      
      
      int real_t=-1;
      for (int kk=0;kk<numSetB;kk++){
        if (index_b[kk]==t){
          real_t=kk;
          break;
        }
      }
      p[real_t]+=epsilon;
      
      add_edge(s,t,pre_Mr);
      num_m_edges++;
//      cout<<num_m_edges<<","<<k<<endl;
      if (num_m_edges == numSetA){
        cout<<numSetA<<","<<k<<endl;
        plot_dat = make_pair(numSetA, k);

        return plot_dat;
      }
    }
  }
  
  
  return plot_dat;
}


int main(int argc, char* argv[]) {
  //// what type of graph do you want to use?
  typedef PUNGraph PGraph; // undirected graph

  // this code is independent of what particular graph implementation/type we use
  printf("Creating graph:\n");
  
//
//  //FACEBOOK DATASET
//  PGraph G;
//  G = LoadEdgeList<PGraph>(TStr::Fmt("%s/facebook.txt", DIRNAME));
//
//  printf("Loaded Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
//

  //RANDOM GENERATE DATASET
//  PGraph G = PGraph::TObj::New();
//  for (int n = 0; n < 1000; n++) {
//    G->AddNode(); // if no parameter is given, node ids are 0,1,...,9
//  }
//  for (int e = 0; e < 20000; e++) {
//    const int NId1 = G->GetRndNId();
//    const int NId2 = G->GetRndNId();
//    if (NId1 != NId2 and G->AddEdge(NId1, NId2) != -2) {
//      printf("  Edge %d -- %d added\n", NId1,  NId2); }
//    else {
//      printf("  Edge %d -- %d already exists\n", NId1, NId2); }
//  }
//  IAssert(G->IsOk());

//  printf("Create Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
//
////   delete
//
//  PGraph::TObj::TNodeI NI = G->GetNI(0);
//  printf("Delete edge %d -- %d\n", NI.GetId(), NI.GetOutNId(0));
//  G->DelEdge(NI.GetId(), NI.GetOutNId(0));
//
//
//  bool bipart[G->GetNodes()];
//  int numSetA=0,numSetB=0;
//  for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
//    bipart[NI.GetId()] = int(rand()%2);
//    if (bipart[NI.GetId()]==0){
//      numSetA+=1;
//    }
//    else{
//      numSetB+=1;
//    }
//  }
//  for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
//    for (int e = 0; e < NI.GetDeg(); e++) {
////      cout<<NI.GetId()<<NI.GetNbrNId(e)<<"\n";
//      if (bipart[NI.GetId()] == bipart[NI.GetNbrNId(e)] ) {
//        G->DelEdge(NI.GetId(), NI.GetNbrNId(e));
////        printf("Delete edge %d -- %d\n", NI.GetId(), NI.GetNbrNId(e));
//        e--;
//      }
//    }
//  }


//semicomplete DATA
  //maximum auction algorithm
  Gnuplot gp;
  gp<<"set term post color\n";
  gp<<"set output 'a.ps'\n";
  
  vector<pair<int, int>> plot_dat;
  
  vector<pair<int, int>> plot_dat_regression;
  
  
  vector<pair<int, int>> plot_dat_regression2;
  gp<<"set xlabel 'Size of maximum matching (size of graph)'\n";
  gp<<"set ylabel 'Number of Iterations to get every matching'\n";
  gp<<"set title 'The relation between  Graph Size  and  Number of Iteration  of Auction Algorithm'\n";
  gp << "set yrange [0:15000]\n";
  gp << "set xrange [0:250]\n";
  
  for (int kkk = 2; kkk<500;kkk=kkk+2){
    PGraph G = PGraph::TObj::New();
    for (int n = 0; n < kkk; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...,9
    }
    bool bipart[G->GetNodes()];
    int numSetA=0,numSetB=0;
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      if (NI.GetId()<(G->GetNodes()/2)){
        bipart[NI.GetId()]=1;
        numSetA+=1;
      }
      else{
        bipart[NI.GetId()]=0;
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

    
    
    E bipartite_edges [G->GetEdges()];
    int index=0;
    // transfer from SNAP to Boost Library
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      for (int e = 0; e < NI.GetDeg(); e++) {
        if (NI.GetId() < NI.GetNbrNId(e)){
          bipartite_edges[index++] = make_pair(NI.GetId(), NI.GetNbrNId(e));
        }
      }
    }

    vector_graph_t bipartite_vector_graph(&bipartite_edges[0], &bipartite_edges[0] + sizeof(bipartite_edges) / sizeof(E), G->GetNodes());
    

    // maximum matching from boost library
    vector< graph_traits< vector_graph_t >::vertex_descriptor > mate(G->GetNodes());
    edmonds_maximum_cardinality_matching(bipartite_vector_graph, &mate[0]);

    plot_dat.push_back(auction(bipartite_vector_graph, bipart, 0.0001, numSetA, numSetB));
    plot_dat_regression.push_back(make_pair(kkk, 87.7563 * pow(1.0231,kkk)));
    plot_dat_regression2.push_back(make_pair(kkk, 0.3136 * pow(kkk,1.874)));
    
  }

  gp << "plot"<<gp.file1d(plot_dat_regression)<<"lw 2 smooth mcsplines title 'y = 87.756 * 1.0231^x',";
  gp <<gp.file1d(plot_dat_regression2)<<"lw 2 smooth mcsplines title 'y = 0.3136 * x^{1.874}',";
  gp<< gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'number of iteration'";
  
  gp <<endl;
  return 0;
}


