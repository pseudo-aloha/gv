#ifndef ECO_MGR_H
#define ECO_MGR_H

#include "cirMgr.h"
#include "ecoNtk.h"
// #include "sat.h"
#include <unordered_set>
#include <map>
#include <queue>
#include <set>
namespace gv {
namespace eco {
// forward decleration of classes
class EcoMgr;
class EcoCir;
// class EcoNtk;
// class CirMgr;


// the main class for Andrew's ECO Approach
class EcoMgr {
public:
  // constructor
  EcoMgr () { _oldNtk = new gv::cir::EcoNtk;
              _newNtk = new gv::cir::EcoNtk; }
  ~EcoMgr () { delete _oldNtk; delete _newNtk; }
  
  // main step function
  void doEco(const string& oldDesignName, const string& newDesignName);
  void readDesigns(const string& oldDesignName, const string& newDesignName);
  void doFraig(); // conduct abc fraig on the designs
  void doMatching();


  // dfs
  void dfs(gv::cir::EcoGate* g) {
    cout << g->getGateFullName() << endl;
    if(isMerged(g)) return;
    assert(g->getGateType() != gv::cir::EcoGate::ECO_PI_GATE && g->getGateType() != gv::cir::EcoGate::ECO_CONST_0_GATE && g->getGateType() != gv::cir::EcoGate::ECO_CONST_1_GATE);
    for(size_t i=0; i<g->getNumFanins(); i++) {
      dfs(g->getFanin(i));
    }
  }

  // record the merge information
  void setMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _mergeTable[g].insert(mg); }
  void setInvMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _invMergeTable[g].insert(mg); }
  unordered_set<gv::cir::EcoGate*> getMergedGates(gv::cir::EcoGate* g) { if(!_mergeTable.count(g)) return {}; return _mergeTable.at(g); }
  unordered_set<gv::cir::EcoGate*> getInvMergedGates(gv::cir::EcoGate* g) { if(!_invMergeTable.count(g)) return {}; return _invMergeTable.at(g); }
  bool isMerged(gv::cir::EcoGate* g) { return (_mergeTable.count(g) || _invMergeTable.count(g)); }
private:
  gv::cir::EcoNtk* _oldNtk;
  gv::cir::EcoNtk* _newNtk;

  // record the merge information
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _mergeTable;
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _invMergeTable;
};

}}

#endif