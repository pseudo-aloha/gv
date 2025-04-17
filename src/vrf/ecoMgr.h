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
class EcoRPInfo;

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

class EcoRPInfo {
public:
  EcoRPInfo(gv::cir::EcoGate* mappedGate, gv::cir::EcoGate* fixedFanout, bool mappedPole) : _mappedGate(mappedGate), _fixedFanout(fixedFanout), _mappedPole(mappedGate) {}
  gv::cir::EcoGate* getMappedGate() { return _mappedGate; }
  gv::cir::EcoGate* getFixedFanout() { return _fixedFanout; }
  bool getMappedPole() { return _mappedPole; }
  
private:
  gv::cir::EcoGate* _mappedGate;
  gv::cir::EcoGate* _fixedFanout;
  bool _mappedPole;
};


// the main class for Andrew's ECO Approach
class EcoMgr {
public:
  // constructor
  EcoMgr () { _oldNtk = new gv::cir::EcoNtk;
              _newNtk = new gv::cir::EcoNtk;
              _patchNtk = new gv::cir::EcoNtk;
              _selectorNtk = new gv::cir::EcoNtk; }
  ~EcoMgr () { delete _oldNtk; delete _newNtk; delete _patchNtk; delete _selectorNtk; }
  
  // main step function
  void doEco(const string& oldDesignName, const string& newDesignName);
  void readDesigns(const string& oldDesignName, const string& newDesignName);
  void doFraig(); // conduct abc fraig on the designs
  void doMatching(unsigned kFeassible);

  // matching functions
  // general matching function
  void matchCutsAtGatePair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, int ithPo = -1); // match the cuts at the gate pair
  bool match2Cuts(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int ithPo = -1);
  vector<size_t>  getMatchWays(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut);
  vector<size_t> simNFindValidMatch(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, vector<int>& comb);
  bool checkMatchValid(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch); // function to check that if the matching is indeed valid
  bool checkMatchValidWithConst(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch, unordered_map<gv::cir::EcoGate*, bool>& constAssignmentOld, unordered_map<gv::cir::EcoGate*, bool>& constAssignmentNew);
  
  // score computation functions
  void sortCandCutsByScore(); // collect the enumerated cuts and compute their scores

  // signature computation functions

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
  unsigned getGatesEqStatus(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate);

  // get the aig of merged aig and pole
  pair<gv::cir::CirGate*, bool> getMergedAig(gv::cir::EcoGate* g);

  // cut hashing things
  // compute cut hash table
  void computeCutHashTable();

  // get NPN class
  pair<string, vector<int>> getNPNHash(gv::cir::EcoCut* cut);
  pair<string, vector<vector<int>>> getNPNHashFull(gv::cir::EcoCut* cut);

  // simulation methods
  void doRandomSim();
  const unsigned getNumSim() { return _nSim; }

  // get Gate similarity functions
  double getCosineSimilarity(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, unsigned poId);

  // Record RP  Pair
  void addRPPair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, gv::cir::EcoGate* fixedFanout, bool inv, unsigned fixedPo);
  void reportRPPair();
  
  // generate patch
  void genPatch();
  void collectPatchGates(gv::cir::EcoGate* g, bool isEntry);
  bool applyNCheckPatch();

  // enum
  enum EQStatus {
    ECO_GATES_NEQ = 0,     // the two gates are not eq
    ECO_GATES_EQ = 1,      // the two gates are eq
    ECO_GATES_INV_EQ = 2   // the two gates are eq after inverted
  };
  

private:
  // input ntks 
  gv::cir::EcoNtk* _oldNtk;
  gv::cir::EcoNtk* _newNtk;
  
  // selector ntk things
  void buildSelector();
  void buildSelectorForIthPo(unsigned i);
  void buildSelectorForIthOldPoRec(gv::cir::EcoGate* selectorGate, gv::cir::EcoGate* origNtkGate, unsigned ithPo);
  void buildSelectorForIthNewPoRec(gv::cir::EcoGate* selectorGate, gv::cir::EcoGate* origNtkGate, unsigned ithPo);
  void addSelectorGate(gv::cir::EcoGate* );
  gv::cir::EcoNtk* _selectorNtk;
  unordered_map<gv::cir::EcoGate*, gv::cir::EcoGate*> _selectorGateMap;

  // patch ntk
  gv::cir::EcoNtk* _patchNtk;

  // record the merge information
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _mergeTable;
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _invMergeTable;

  // record the match information
  unordered_map<gv::cir::EcoGate*, vector<pair<gv::cir::EcoGate*, bool>>> _matchTable;

  // record all the cand cuts for matching
  vector<gv::cir::EcoCut*> _oldCandCuts;
  vector<gv::cir::EcoCut*> _newCandCuts;

  // record the NPN hash information
  EcoNPNHash* _pNpnHash;

  // record the RP pair
  vector<unordered_map<gv::cir::EcoGate*, EcoRPInfo*>> _rpTable; // record the ith circuit in old gate matched to jth gate in new circuit and also record the pole

  // record how many patterns has been simmed
  unsigned _nSim;
};

}
}

#endif