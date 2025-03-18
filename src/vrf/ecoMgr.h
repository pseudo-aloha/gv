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
  
  // pre-compute things
  void computeNpnMatchWays(int k);
  void npnHash(const vector<vector<int>>& npnMatchWays, int i, int k);
  void computeNpnHash();

  // get NPN hash
  pair<string, vector<int>> getNPNHash(size_t cutTT, unsigned cutSize);
  pair<string, vector<vector<int>>> getNPNHashFull(size_t cutTT, unsigned simSize);

  // encode function (used to save memory)
  size_t encodeMatch2SizeT(int outputMatch, vector<int>& inputMatch);
  pair<int, vector<int>> decodeEncodedSizeTMatch(size_t encode, unsigned cutSize);

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
  void doMatching(unsigned kFeassible);

  // matching functions
  // general matching function
  void matchCutsAtGatePair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate); // match the cuts at the gate pair
  bool match2Cuts(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut);
  pair<int, vector<int>> getOneMatchWay(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut);
  vector<size_t> simNFindValidMatch(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, vector<int>& comb);
  bool checkMatchValid(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch); // function to check that if the matching is indeed valid
  bool checkMatchValidWithConst(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch, unordered_map<gv::cir::EcoGate*, bool>& constAssignmentOld, unordered_map<gv::cir::EcoGate*, bool>& constAssignmentNew);
  
  // output side matching functions
  void doOutputSideMatching();
  void matchOnePo(unsigned ithPo);

  // input side prepatch functions

  // recycle matching functions

  // record the merge information
  void setMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _mergeTable[g].insert(mg); }
  void setInvMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _invMergeTable[g].insert(mg); }
  unordered_set<gv::cir::EcoGate*> getMergedGates(gv::cir::EcoGate* g) { if(!_mergeTable.count(g)) return {}; return _mergeTable.at(g); }
  unordered_set<gv::cir::EcoGate*> getInvMergedGates(gv::cir::EcoGate* g) { if(!_invMergeTable.count(g)) return {}; return _invMergeTable.at(g); }
  bool isMerged(gv::cir::EcoGate* g) { return (_mergeTable.count(g) || _invMergeTable.count(g)); }

  // cut hashing things
  // compute cut hash table
  void computeCutHashTable();
  // get NPN class
  pair<string, vector<int>> getNPNHash(gv::cir::EcoCut* cut);

  // Record RP  Pair
  void addRPPair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate) { _RPPair[oldGate] = newGate; }
  gv::cir::EcoGate* getRPGate(gv::cir::EcoGate* oldGate) { if(!_RPPair.count(oldGate)) return nullptr; return _RPPair.at(oldGate);}

  // generate patch
  

private:
  gv::cir::EcoNtk* _oldNtk;
  gv::cir::EcoNtk* _newNtk;

  // record the merge information
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _mergeTable;
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _invMergeTable;

  // record the NPN hash information
  EcoNPNHash* _pNpnHash;

  // record the RP pair
  unordered_map<gv::cir::EcoGate*, gv::cir::EcoGate*> _RPPair; // record the RP pair
};

}
}

#endif