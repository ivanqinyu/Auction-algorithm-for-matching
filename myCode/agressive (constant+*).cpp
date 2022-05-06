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

vector<pair<int, int>> auction(vector_graph_t bipartite_vector_graph, vector<bool> bipart,float epsilon, unsigned long numSetA, unsigned long numSetB, float agressive, int maximum_matching){
  
  vector<pair<int, int>> plot_dat;
  plot_dat.push_back(make_pair(0, 0));
  
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
    if (bipart[i]==0){
      index_a[ii]=i;
      ii++;
    }
    else{
      index_b[jj]=i;
      jj++;
    }
  }
//  cout<<ii<<",,"<<jj<<endl;
  float price_epsilon = epsilon;
  vector_graph_t pre_Mr((numSetA+numSetB));
//  cout<<numSetA<<","<<numSetB<<endl;
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

    for (int i=0;i<numSetA;i++){

      if (a[i]==-1){
        float min_d = 99999;

        int s = index_a[i];
        graph_traits<vector_graph_t>::out_edge_iterator ei, ei_end;
        for (boost::tie(ei, ei_end) = out_edges(index_a[i], bipartite_vector_graph); ei != ei_end; ++ei){
          int t = target(*ei, bipartite_vector_graph);
//          cout<<"s,t  in  Gr =  "<<s <<", "<<t<<"s bipart = "<<bipart[s]<<" t bipart = "<<bipart[t]<<endl;
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
//        cout<<"min_d"<<min_d<<endl;
        if (min_d<99999){
          for (int kk=0;kk<numSetB;kk++){
//            cout<<"before Add  Gr: "<<s<<" "<<index_b[kk]<<" "<<p[kk]<<" bool=="<< bool(p[kk]-min_d<0.000001) <<" "<< edge(s,index_b[kk],bipartite_vector_graph).first<<" +++ "<< edge(s,index_b[kk],bipartite_vector_graph).second <<endl;
            if (p[kk]-min_d<0.000001 and edge(s,index_b[kk],bipartite_vector_graph).second == true){
              add_edge(s,index_b[kk],Gr);
//              cout<<"Add  Gr: "<<s<<" "<<index_b[kk]<<"    a[i]="<<a[i]<<"  price ="<<min_d<<endl;
            }
          }
          
        }
      }
    }
    
    
    if (num_edges(Gr)==0){
      cout<<"------iteration: "<<k<<"     number_edge: "<<num_m_edges<<endl;
      return plot_dat;
    }
//    cout<<"Gr: "<<num_edges(Gr)<<endl;
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
        
//        cout<<"s,t  in  Gr2 =  "<<s <<", "<<t<<"s bipart = "<<bipart[s]<<" t bipart = "<<bipart[t]<<endl;
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
          add_edge(s,t,Mr);
          
        }
      }
    }
    
//    cout<<"ADD Mr: "<<num_edges(Mr)<<endl;
    
    
    graph_traits<vector_graph_t>::edge_iterator ei, ei_end;
    
    for (boost::tie(ei, ei_end) = edges(Mr); ei != ei_end; ++ei){
      int s = source(*ei, Mr);
      int t = target(*ei, Mr);
//      cout<<"s,t  in  Mr =  "<<s <<", "<<t<<"s bipart = "<<bipart[s]<<" t bipart = "<<bipart[t]<<endl;
      
//      cout<<boost::degree(3, pre_Mr)<<" + + +"<<boost::degree(76, pre_Mr)<<endl;
      
      if (boost::degree(t, pre_Mr)!=0 ){
        graph_traits<vector_graph_t>::out_edge_iterator eei, eei_end;
        boost::tie (eei,eei_end) = out_edges(t,pre_Mr);
        
        int pre_val = target(*eei ,pre_Mr);
//        cout<<"remove"<<boost::degree(t, pre_Mr)<<"   "<<t<<","<<pre_val<<endl;

        if (pre_val != s){
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
//          cout<<"removeeeeee a["<<real_pre<<"]"<<"  = "<<-1<<endl;
        }
      }
      for (int kk=0;kk<numSetA;kk++){
        if (index_a[kk]==s){
          a[kk]=t;
//          cout<<"adddd a["<<kk<<"]"<<"  = "<<t<<endl;
          for (int pp=0;pp<numSetB;pp++){
            if (index_b[pp]==t){
              p[pp]+=price_epsilon;
              
//              cout<<"adddd p["<<pp<<"]"<<endl;
              break;
            }
          }
          break;
        }
      }
      

    
      add_edge(s,t,pre_Mr);
//        cout<<"add edge: "<<s<<","<<t<<endl;
      num_m_edges++;
      
    
      if (num_m_edges == maximum_matching){
//        graph_traits<vector_graph_t>::edge_iterator eei, eei_end;
//
//        for (boost::tie(eei, eei_end) = edges(pre_Mr); eei != eei_end; ++eei){
//          cout<<source(*eei, pre_Mr)<<","<<target(*eei, pre_Mr)<<endl;
//        }
//
        cout<<"------iteration: "<<k<<"     number_edge: "<<num_m_edges<<endl;
        plot_dat.push_back(make_pair(k, num_m_edges));
        return plot_dat;
      }
      
    }
    
    plot_dat.push_back(make_pair(k, num_m_edges));
//    graph_traits<vector_graph_t>::edge_iterator eei, eei_end;
    
//    for (boost::tie(eei, eei_end) = edges(pre_Mr); eei != eei_end; ++eei){
//      cout<<source(*eei, pre_Mr)<<","<<target(*eei, pre_Mr)<<endl;
//    }
    
//    cout<<k<<","<<num_m_edges<<endl;
  }
  
  cout<<"M has "<<num_m_edges<<" edges with auction algorithm"<<endl;
  return plot_dat;
}

#define SEMICOMPLETE 0;
#define RANDOM_GRAPH 1;
#define FACEBOOK_DATA 2;
#define Twitter_DATA 3;

int main(int argc, char* argv[]) {
  typedef PUNGraph PGraph; // undirected graph
  printf("Creating graph:\n");
  PGraph G = PGraph::TObj::New();
  unsigned long numSetA=0,numSetB=0;
  vector<bool> bipart;
  
  // what type of graph do you want to use?
  int selection = SEMICOMPLETE;

  if (selection == 0){
//  semicomplete DATA
    for (int n = 0; n < 1000; n++) {
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
    for (int n = 0; n < 50; n++) {
      G->AddNode(); // if no parameter is given, node ids are 0,1,...,9
      bipart.push_back(false);
    }
    for (int e = 0; e < 100; e++) {
      const int NId1 = G->GetRndNId();
      const int NId2 = G->GetRndNId();
      if (NId1 != NId2 and G->AddEdge(NId1, NId2) != -2) {
        printf("  Edge %d -- %d added\n", NId1,  NId2); }
      else {
        printf("  Edge %d -- %d already exists\n", NId1, NId2); }
    }
    IAssert(G->IsOk());
    //   delete
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      bipart[NI.GetId()] = int(rand()%2);
      cout<<bipart[NI.GetId()]<<endl;
      if (bipart[NI.GetId()]==0){
        numSetA+=1;
      }
      else{
        numSetB+=1;
      }
    }
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      for (int e = 0; e < NI.GetDeg(); e++) {
        cout<<NI.GetId()<<NI.GetNbrNId(e)<<"\n";
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
    G = LoadEdgeList<PGraph>(TStr::Fmt("%s/facebook.txt", DIRNAME));

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
    G = LoadEdgeList<PGraph>(TStr::Fmt("%s/twitter.txt", DIRNAME));

    printf("Loaded Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
    
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      bipart.push_back(0);
    }
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      cout<<"pass1   "<<NI.GetId()<<endl;
      bipart[NI.GetId()] = int(rand()%2);
      cout<<"pass2"<<endl;
      if (bipart[NI.GetId()]==0){
        numSetA+=1;
      }
      else{
        numSetB+=1;
      }
      cout<<"pass3"<<endl;
    }
    cout<<"pass"<<endl;
    for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      for (int e = 0; e < NI.GetDeg(); e++) {
        if (bipart[NI.GetId()] == bipart[NI.GetNbrNId(e)] ) {
          G->DelEdge(NI.GetId(), NI.GetNbrNId(e));
          e--;
        }
      }
    }

  }
  printf("Create Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());

//   dump the graph
//  printf("Graph (%d, %d)\n", G->GetNodes(), G->GetEdges());
//  for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
//    printf("  %d: ", NI.GetId());
//    for (int e = 0; e < NI.GetDeg(); e++) {
//      printf(" %d", NI.GetNbrNId(e));
//    }
//    printf("\n");
//  }

//  Draw the "PNG" pictures
//  mkdir(DIRNAME, S_IRWXU | S_IRWXG | S_IRWXO);
//
//  TStrV LNames; //  gvlDot, gvlNeato, gvlTwopi, gvlCirco
//  LNames.Add("Dot");
//  LNames.Add("Neato");
//  LNames.Add("Twopi");
//  LNames.Add("Circo");
//
//  TStrV Exts;
//  Exts.Add("png");
//
//  for (int i = 0; i < LNames.Len(); i++) {
//    for (int e = 0; e < Exts.Len(); e++) {
//      // Baseline file is already created (use as benchmark)
//      TStr FNameDemo = TStr::Fmt("%s/demo_%s_%s.%s", DIRNAME,
//                                 "ungraph" ,
//                                 LNames[i].CStr(), Exts[e].CStr());
//
//      // Remove test graph if it already exists
//      remove(FNameDemo.CStr());
//
//      // Draw new graph
//      TSnap::DrawGViz(G, TGVizLayout(i), FNameDemo, LNames[i], true);
//
//      printf("Drawing graph '%s'\n", FNameDemo.CStr());
//    }
//  }
  
  // transfer from SNAP to Boost Library
  E bipartite_edges [G->GetEdges()];
  unsigned long index=0;
  for (PGraph::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    for (int e = 0; e < NI.GetDeg(); e++) {
      if (NI.GetId() < NI.GetNbrNId(e)){
        bipartite_edges[index++] = make_pair(NI.GetId(), NI.GetNbrNId(e));
      }
    }
  }

  vector_graph_t bipartite_vector_graph(&bipartite_edges[0], &bipartite_edges[0] + sizeof(bipartite_edges) / sizeof(E), G->GetNodes());

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
  
  // maximum matching from boost library
  vector< graph_traits< vector_graph_t >::vertex_descriptor > mate(G->GetNodes());
  edmonds_maximum_cardinality_matching(bipartite_vector_graph, &mate[0]);
  int maximum_matiching = matching_size(bipartite_vector_graph, &mate[0]);
  cout<<"M has "<< maximum_matiching <<" edges with build-in algorithm"<< std::endl;

//  std::cout << "The matching is:" << std::endl;
//  graph_traits< vector_graph_t >::vertex_iterator vi, vi_end;
//  for (boost::tie(vi, vi_end) = vertices(bipartite_vector_graph); vi != vi_end; ++vi)
//    if (mate[*vi] != graph_traits< vector_graph_t >::null_vertex() && *vi < mate[*vi])
//      std::cout << "{" << *vi << ", " << mate[*vi] << "}" << std::endl;
//
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
  gp << "set xrange [0:400]\n";

  plot_dat = auction(bipartite_vector_graph, bipart, 0.01, numSetA, numSetB, 0,maximum_matiching);
  gp <<"plot "<<maximum_matiching<<"title 'maximum',"<<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'constant epsilon 0.01',";

  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 0,maximum_matiching);
  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'constant epsilon 0.001',";

  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, -1,maximum_matiching);
  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon +0.0001',";

  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 1.1,maximum_matiching);
  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon *1.1',";
  
  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 1.001,maximum_matiching);
  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon *1.001',";
  
  plot_dat = auction(bipartite_vector_graph, bipart, 0.001, numSetA, numSetB, 1.0001,maximum_matiching);
  gp <<gp.file1d(plot_dat) << "lw 2 smooth mcsplines title 'agress epsilon *1.0001'";
  
  gp <<endl;
  return 0;
}
