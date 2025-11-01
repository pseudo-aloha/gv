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
typedef vector<int> ConstInsertList;

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
  size_t encodeMatch2SizeT(int outputMatch, vector<int>& inputMatch, const ConstInsertList& ConstInsert = {});
  pair<int, vector<int>> decodeEncodedSizeTMatch(size_t encode, unsigned cutSize);

  // used to indicate const0/const1 in size_t encode
  enum ConstInsertSizetEncode {
    CONST0SIZETENCODE = 6,
    CONST1SIZETENCODE = 7
  };

private:
  unsigned _cutSizeFrom; // from which we comute the NPN-eq class
  unsigned _cutSizeTo; // to which we comute the NPN-eq class
  vector<vector<pair<string, vector<int>>>> _npnHashTable;
  vector<vector<int>> _npnHashWays;
};

class EcoRPInfo {
public:
  EcoRPInfo(gv::cir::EcoGate* mappedGate, bool mappedPole) : _mappedGate(mappedGate), _mappedPole(mappedPole) {}
  gv::cir::EcoGate* getMappedGate() { return _mappedGate; }
  bool getMappedPole() { return _mappedPole; }
  
private:
  gv::cir::EcoGate* _mappedGate;
  bool _mappedPole;
};


// the main class for Andrew's ECO Approach
class EcoMgr {
public:
  // constructor
  EcoMgr () : _patchWireCount(0) { _oldNtk = new gv::cir::EcoNtk;
              _newNtk = new gv::cir::EcoNtk;
              _patchNtk = new gv::cir::EcoNtk;
              _selectorNtk = new gv::cir::EcoNtk; }
  ~EcoMgr () { delete _oldNtk; delete _newNtk; delete _patchNtk; delete _selectorNtk; }
  
  // main step function
  void doEco(const string& oldDesignName, const string& newDesignName);
  void readDesigns(const string& oldDesignName, const string& newDesignName);
  void doFraig(); // conduct abc fraig on the designs
  void doMatching(unsigned kFeassible);
  void doRecycle();

  // set function
  void setOldDesignName(const string& name) { _oldDesignName = name; }
  void setNewDesignName(const string& name) { _newDesignName = name; }

  // get function
  string getOldDesignName() { return _oldDesignName; }
  string getNewDesignName() { return _newDesignName; }
  
  // ---------------------
  // matching functions
  // ---------------------
  // general matching function
  unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> matchCutsAtGatePair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, int ithPo = -1, bool doConstInsert = false); // match the cuts at the gate pair
  void computeCutsSignatures(unordered_map<string, vector<pair<gv::cir::EcoCut*, ConstInsertList>>>& NPNClass2Cuts, vector<gv::cir::EcoCut *>& cuts, bool doConstInsert = false);
  pair<bool, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>>> match2Cuts(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, const ConstInsertList& ConstInsert, int ithPo = -1);
  vector<size_t>  getMatchWays(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, const ConstInsertList& ConstInsert);
  vector<size_t> simNFindValidMatch(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, vector<int>& comb, const ConstInsertList& ConstInsert);
  bool checkMatchValid(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch); // function to check that if the matching is indeed valid
  bool checkMatchValidWithConst(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch, const ConstInsertList& constInsertOld, const ConstInsertList& constInsertNew);
  
  // ---------------------
  // score computation functions
  // ---------------------
  // void sortCandCutsByScore(); // collect the enumerated cuts and compute their scores
  void sortCutsByNumMergedGates(vector<gv::cir::EcoCut*>& cuts);
  void sortCutsByNumMergedGates(vector<pair<gv::cir::EcoCut*, ConstInsertList>>& cuts);
  vector<string> sortNPNClass(const unordered_map<string, vector<pair<gv::cir::EcoCut*, ConstInsertList>>>& newNPNClass2Cuts);

  // signature computation functions

  // output side matching functions
  void doOutputSideMatching();
  bool check2ConeEq(int ithPo);
  void matchOnePo(unsigned ithPo);

  // input side prepatch functions

  // recycle matching functions

  // record the merge information
  void setMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _mergeTable[g].insert(mg); }
  void setInvMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* mg) { _invMergeTable[g].insert(mg); }
  unordered_set<gv::cir::EcoGate*> getMergedGates(gv::cir::EcoGate* g) { if(!_mergeTable.count(g)) return {}; return _mergeTable.at(g); }
  unordered_set<gv::cir::EcoGate*> getInvMergedGates(gv::cir::EcoGate* g) { if(!_invMergeTable.count(g)) return {}; return _invMergeTable.at(g); }
  pair<gv::cir::EcoGate*, bool> getOneMergedGate(gv::cir::EcoGate* g, bool pole); // the first arguement is the gate we want to find merged gate, the second arguement is the preferred pole of the gate.
  bool isMerged(gv::cir::EcoGate* g) { return (_mergeTable.count(g) || _invMergeTable.count(g)); }
  unsigned getGatesEqStatus(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate);
  void addUnderMergeFrontierSet(gv::cir::EcoGate* g) { _underMergeFrontierSet.insert(g); }
  bool isUnderMergeFrontierSet(gv::cir::EcoGate* g)  { _underMergeFrontierSet.count(g);  }
  void markMergeFrontier();

  // get the aig of merged aig and pole
  pair<gv::cir::CirGate*, bool> getMergedAig(gv::cir::EcoGate* g);

  // cut hashing things
  // compute cut hash table
  void computeCutHashTable();

  // get NPN class
  pair<string, vector<int>> getNPNHash(gv::cir::EcoCut* cut);
  pair<string, vector<vector<int>>> getNPNHashFull(gv::cir::EcoCut* cut, const ConstInsertList& constInsert = {});

  // simulation methods
  void doRandomSim();
  const unsigned getNumSim() { return _nSim; }

  // get Gate similarity functions
  double getCosineSimilarity(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, unsigned poId);

  // Record RP  Pair
  void addRPPair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, bool inv, unsigned fixedPo);
  void reportRPPair();
  
  // ---------------------
  // generate patch
  // ---------------------
  void genPatch(const string& patchName);
  void decideOutputRewire();
  // write patch ntk functions
  void generateFinalRpRewire();
  void generatePatchForIthPo(unsigned i);
  void generateNewGateLogics();
  // patch check functions
  bool applyNCheckPatch(const string& patchName);
  // patch generation helper functions
  void addPatchPoName(const string& name) { _patchPoNames.insert(name); }
  bool isInPatchPoNames(const string& name) const { return _patchPoNames.count(name); }
  string getPatchGateName(gv::cir::EcoGate* g); //decide whether the old gate in the patch needs to add _in suffix
  string getPatchWireName() { return "eco_wire_" + to_string(_patchWireCount++); }
  gv::cir::EcoGate* createOrGetPatchGate(const string& gateTypeName, const string& gateName);
  void dupMergedGates();
  void addDupedMergedGate(gv::cir::EcoGate* g, gv::cir::EcoGate* dupG) { _dupGateMap[g] = dupG; }
  gv::cir::EcoGate* getDupedMergedGate(gv::cir::EcoGate* g) { if(!_dupGateMap.count(g)) return nullptr; return _dupGateMap.at(g); }
  string getDupGateName(gv::cir::EcoGate* g, unsigned ithPo);
  bool checkIfHasToDup(gv::cir::EcoGate* g, unsigned ithPo);
  bool checkIfHasToFixPo(unsigned ithPo);
  void addNeedTodupGate(gv::cir::EcoGate* g) { if(g->isConstGate()) assert(0); /* donot want to dup const gate*/ _gatesNeedToDup.insert(g); }
  bool isNeedToDup(gv::cir::EcoGate* g) { return _gatesNeedToDup.count(g); }
  bool isFixedToAnotherGate(gv::cir::EcoGate* g);
  bool isFixedToItSelf(gv::cir::EcoGate* g);

  // recycle
  void collectFloatingGates();
  void matchFloatingGates();
  
  // resynthesis circuit
  void reSynsethesis(const string& oldDir, const string& newDir);

  // ---------------------
  // enum definitions
  // ---------------------
  // Used to indicate if two gates eq status, used when we want to find the merged gates
  enum EQStatus {
    ECO_GATES_NEQ = 0,     // the two gates are not eq
    ECO_GATES_EQ = 1,      // the two gates are eq
    ECO_GATES_INV_EQ = 2   // the two gates are eq after inverted
  };
  

private:
  // input ntks 
  gv::cir::EcoNtk* _oldNtk;
  gv::cir::EcoNtk* _newNtk;

  string _oldDesignName;
  string _newDesignName;
  
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
  unsigned _patchWireCount;
  unordered_set<string> _patchPoNames; // used to record the patch po names, and add "_in" string for gate used for both po and pi in patch circuit
  unordered_set<gv::cir::EcoGate*> _usedNewGate;
  unordered_map<gv::cir::EcoGate*, unsigned> _oldGateUsedByPo; // record the old gate that is use by po i to fix itself
  unordered_map<gv::cir::EcoGate*, gv::cir::EcoGate*> _dupGateMap; // use the original gate to find the duplicated gate in the patch circuit
  unordered_set<gv::cir::EcoGate*> _gatesNeedToDup;
  unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> _reWiredPiMap;

  // record the merge information
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _mergeTable;
  unordered_map<gv::cir::EcoGate*, unordered_set<gv::cir::EcoGate*>> _invMergeTable;
  unordered_set<gv::cir::EcoGate*>                                   _underMergeFrontierSet; // record the gates that are strictly under the merge frontier

  // record the match information
  unordered_map<gv::cir::EcoGate*, vector<pair<gv::cir::EcoGate*, bool>>> _matchTable;

  // record all the cand cuts for matching
  vector<gv::cir::EcoCut*> _oldCandCuts;
  vector<gv::cir::EcoCut*> _newCandCuts;

  // record the NPN hash information
  EcoNPNHash* _pNpnHash;

  // record the RP pair
  vector<unordered_map<gv::cir::EcoGate*, EcoRPInfo*>> _rpTable; // record the ith circuit in old gate matched to jth gate in new circuit and also record the pole
  unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> _finalRpPair;

  // record how many patterns has been simmed
  unsigned _nSim;

  // recycle
  unordered_set<gv::cir::EcoGate*> _floatingGates;
};

}
}

#endif