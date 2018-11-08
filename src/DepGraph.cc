#include "DepGraph.h"

//using namespace boost;
//using namespace std;
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> mygraph;
#define cutint(x) static_cast<int>(x)
#define intcut(x) static_cast<CUTS>(x)

DepGraph::DepGraph() {
  
  add_edge(cutint(CUTS::eGTau), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGTop), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGElec), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGMuon), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGZ), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGW), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGHiggs), cutint(CUTS::eGen), g);
  add_edge(cutint(CUTS::eGJet), cutint(CUTS::eGen), g);
  
}

    
void DepGraph::dfs(int vertex) {
  mygraph::adjacency_iterator it, end;
  neededCuts.insert(vertex);

  for(tie(it, end) = adjacent_vertices(vertex, g); it != end; ++it) {
    if(neededCuts.find(*it) != neededCuts.end()) continue;
    dfs(*it);
  }
  return;
}


void DepGraph::loadCuts(std::vector<CUTS> cutVec) {
  for(auto cut: cutVec) {
    int icut = cutint(cut);
    if(neededCuts.find(icut) != neededCuts.end()) continue;
    
    dfs(icut);
  }
}

void DepGraph::loadCuts(CUTS cut) {
  int icut = cutint(cut);
  if(neededCuts.find(icut) != neededCuts.end()) return;

  dfs(icut);
}

bool DepGraph::isPresent(CUTS cut) {
  return (neededCuts.find(cutint(cut)) != neededCuts.end());
}

std::unordered_set<int> DepGraph::getCuts() {
  return neededCuts;
}
