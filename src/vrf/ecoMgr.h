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
class NPNHash;

// class to store the cut hashing things
class EcoNPNHash {
public:
  // constructor / destructor
  EcoNPNHash() : _cutSizeFrom(2), _cutSizeTo(4) {}
  EcoNPNHash(unsigned cutSizeFrom, unsigned cutSizeTo) : _cutSizeFrom(cutSizeFrom), _cutSizeTo(cutSizeTo) {}
  ~EcoNPNHash() {}

  void computeNpnMatchWays(int k);
  void npnHash(const vector<vector<int>>& npnMatchWays, int i, int k);
  void computeNpnHash();

private:
  unsigned _cutSizeFrom; // from which we comute the NPN-eq class
  unsigned _cutSizeTo; // to which we comute the NPN-eq class
  vector<vector<pair<string, vector<int>>>> _npnHashTable;
  vector<vector<int>> _npnHashWays;
};


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

  // matching functions
  // output side matching functions
  void doOutputSideMatching();
  // input side prepatch functions

  // recycle matching functions

  // record the merge information
  void setMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _mergeTable[g].insert(mg); }
  void setInvMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _invMergeTable[g].insert(mg); }
  unordered_set<gv::cir::EcoGate*> getMergedGates(gv::cir::EcoGate* g) { if(!_mergeTable.count(g)) return {}; return _mergeTable.at(g); }
  unordered_set<gv::cir::EcoGate*> getInvMergedGates(gv::cir::EcoGate* g) { if(!_invMergeTable.count(g)) return {}; return _invMergeTable.at(g); }
  bool isMerged(gv::cir::EcoGate* g) { return (_mergeTable.count(g) || _invMergeTable.count(g)); }

  // compute cut hash table
  void computeCutHashTable();

private:
  gv::cir::EcoNtk* _oldNtk;
  gv::cir::EcoNtk* _newNtk;

  // record the merge information
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _mergeTable;
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _invMergeTable;

  // record the NPN hash information
  EcoNPNHash* _pNpnHash;
};

}
}

#endif