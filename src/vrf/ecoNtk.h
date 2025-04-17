#ifndef ECO_NTK_H
#define ECO_NTK_H

#include "cirMgr.h"
#include "cirGate.h"
#include <vector>
#include <limits>

#include "abcMgr.h"
#include "yosysMgr.h"
#include "cirGate.h"
#include "cirDef.h"
// using namespace std;

// TODO: Feel free to define your own classes, variables, or functions.
#include "gvType.h"
namespace gv {
namespace cir {
class EcoNtk;
class EcoGate;
class EcoCir;
class EcoCirGate;
class EcoCut;
}}

namespace gv {
  namespace eco {
    class EcoMgr;
  }}

namespace gv {
namespace cir {

  class EcoGate {
    public:
      friend class EcoNtk;
      EcoGate(string gateType, string gateName);
      ~EcoGate();
      string getGateName() { return _gateName; }
      string getGateFullName();
      string getGateTypeName();
      unsigned getGateType() { return _gateType; }
      unsigned getNumFanins() { return _fanins.size(); }
      EcoGate* getFanin(unsigned i) { if(i>=_fanins.size()) return nullptr; return _fanins.at(i); }
      void reportGate();
      bool getAigNodeInv() { return ecoGateVComp; }
      CirGate* getAigNode() { return ecoGateV; }

      // record the gate belongs to which ntk
      unsigned getGateNtk() { return _gateNtk; }
      void setOld(unsigned gateNtk) { _gateNtk = gateNtk; }

      // add function
      void addFanin(EcoGate* g) { _fanins.push_back(g); }
  
      // traversal things
      static void setGlobalTrav() { _globalTravFlag++; }
      void setGlobalTrav(unsigned i) { _globalTravFlag += i; }
      void setToGlobalTrav() { _travFlag = _globalTravFlag; }
      bool isGlobalTrav() { return (_travFlag == _globalTravFlag);}

      // sim val
      const size_t getSimVal(unsigned i) { if(_gateType == ECO_CONST_0_GATE) return size_t(0); if(_gateType == ECO_CONST_1_GATE) return size_t(std::numeric_limits<size_t>::max());  return _simVals.at(i); }

      // gate type check function
      bool isConstGate() { return (getGateType() == ECO_CONST_0_GATE || getGateType() == ECO_CONST_1_GATE); }
      bool isPi()        { return (getGateType() == ECO_PI_GATE); }
      bool isPiOrConst() { return (isConstGate() || isPi()); }

    enum EcoGateType {
      ECO_CONST_0_GATE = 0,
      ECO_CONST_1_GATE = 1,
      ECO_AND_GATE = 2,
      ECO_OR_GATE = 3,
      ECO_NAND_GATE = 4,
      ECO_NOR_GATE = 5,
      ECO_XOR_GATE = 6,
      ECO_XNOR_GATE = 7,
      ECO_BUF_GATE = 8,
      ECO_NOT_GATE = 9,
      ECO_PI_GATE = 10,
      ECO_PO_GATE = 11
    };

    // record the gate belong to which ntk
    enum EcoGateNtk {
      ECO_OLD_NTK = 0,
      ECO_NEW_NTK = 1,
      ECO_NONE_NTK = 2 // other ntk, like selector or patch ntk etc.
    };
  private:
    // Gate attributes
    unsigned _gateType; // store the gate type e.g. and / or / not
    string _gateName; // store the gate name (the output net name)

    // Gate fanins
    vector<EcoGate*> _fanins;
    vector<string> _faninNames;

    // internal gate
    CirGate* ecoGateV; // internal AIG node that maps to the EcoGate
    bool ecoGateVComp; // record if the gate's polation when mapping to the internal AIG node
    Abc_Obj_t* _pAbcNode;

    // traverse flag
    static unsigned _globalTravFlag; // global trav flag, shared by all the gates. Used to check if the gate is traversed (increment it before traversing)
    unsigned _travFlag;
    
    
    unsigned _gateNtk; // record gate belongs to old/new circuiit or others...

    // simulation stuffs
    vector<size_t> _simVals;
};


// Cuts in the EcoNtk
class EcoCut {
public:
  // consturctors / destructors
  // EcoCut() : _root(nullptr), _leaves({}), _signature("") {}
  // EcoCut(EcoGate* root) : _leaves({}), _signature("") { _root = root; }
  EcoCut(EcoGate* root, unordered_set<EcoGate*> leaves) : _signature("") { _root = root; _leaves = leaves; }
  EcoCut(EcoGate* root, unordered_map<EcoGate*, int> leaves) : _signature("") { _root = root; for(auto&[leaf, cnt] : leaves) _leaves.insert(leaf); }
  ~EcoCut() { _root = nullptr; _leaves.clear(); _signature.clear(); }
  
  // Basic setting functions
  void setRoot(EcoGate* g) { _root = g; }
  void addLeaf(EcoGate* g) { _leaves.insert(g); }
  void setSig(const string& sig) { _signature = sig; }
  void setMgAigSig(const vector<pair<string, bool>>& sig) { _mgAigSig = sig; }
  void setNPNClass(const string& npnClass) { _npnClass = npnClass; }

  // get functions
  EcoGate* const getRoot() { return _root; }
  unordered_set<EcoGate*> const getLeaves() { return _leaves; }
  unsigned getCutSize() const { return _leaves.size(); }
  static unsigned getMaxCutsPerNode() { return _maxCutsPerNode; }
  string getSig() { return _signature; }
  const vector<pair<string, bool>> getMgAigSig() { return _mgAigSig; }
  unsigned getNumMergedLeaves() {return _numMergedLeaves;}
  const string getNPNClass() const { return _npnClass; }

  void setNumMergedLeaves(unsigned i) { _numMergedLeaves = i; };

  // collect the gates within the cut
  vector<EcoGate*> collectCurConeGate();
  void collectCurConeGateRec(EcoGate* g, vector<EcoGate*>& gateList);

  // report functions
  void reportCut();

private:
  // basib members for a cut
  EcoGate* _root;
  unordered_set<EcoGate*> _leaves;

  // used to limit the max # of cuts per node
  static unsigned _maxCutsPerNode;

  // signature things
  string _signature; // used to uniquefy cuts
  vector<pair<string, bool>> _mgAigSig; // signature that is formed by the ids of merged AIG on the cut  
  string _npnClass; // the npn class of the cut
  unsigned _numMergedLeaves; // used for sorting the score of a cut
};

// Wrap CirMgr, modified to store some extra information for ECO usage
class EcoCir {
  
  public:
    EcoCir() { _ecoCirV = new CirMgr; }
    ~EcoCir() { delete _ecoCirV; }
    void readCirFromAbcNtk(Abc_Ntk_t* pNtk);
    CirGate* getGate(unsigned gid) const { _ecoCirV->getGate(gid); }
    void printSummary() { _ecoCirV->printSummary(); }
    CirMgr* getEcoCirV() { return _ecoCirV; }
  private:
    gv::cir::CirMgr* _ecoCirV;
    unordered_map<CirGate*, vector<CirGate*>> _fanoutMap; // record the fanout mapping of the gates
  
  
  };

// used to store the primitive gate level network
class EcoNtk {
  public:
    friend class EcoMgr;
    // constructor
    EcoNtk () { cirV = new EcoCir(); }
    ~EcoNtk () { delete cirV; }
    // file parsing functions
    void readNtkFile(const string& dir);
    void rewriteDesign(const string& dir);
    void abcReadFile();
    void parsePrimitiveGates(const string& dir);
    void parsePI(const string& dir);
    void parsePO(const string& dir);
    void parseGate(const vector<string>& line);
    void genConnection();

    // cut enumeration function
    void enumerateCuts(unsigned k, gv::eco::EcoMgr* pEco); // enumerate k-feasible cuts
    vector<EcoCut*> enumerateCutsRec(const unsigned& k, EcoGate* g, gv::eco::EcoMgr* pEco); // enumerate k-feasible cuts
    void getCutCombs(unsigned faninIdx, EcoGate* root, unordered_map<EcoGate*, int>& leaves, vector<EcoCut*>& cuts, vector<vector<EcoCut*>>& faninCutVec, const unsigned& k, gv::eco::EcoMgr* pEco);
    bool checkCut(EcoCut* pCut); // rule out the invalid cuts
    bool checkCutRec(EcoCut* pCut, EcoGate* g, bool isSelf); // rule out the invalid cuts

    // cut utils
    size_t computeCutTT(EcoCut* pCut); // compute the truth table of the cut
    size_t computeCutTTWithConst(EcoCut* pCut, const vector<pair<int, bool>>& constAssignment); // compute the truth table of the cut with constant inserted

    // set functions
    void setGateByAbcNode(Abc_Obj_t* pObj, EcoGate* pEcoGate) { _abcObj2EcoGate[Abc_ObjRegular(pObj)].insert(pEcoGate); }

    // add function
    void addPo(EcoGate* g) { _POList.push_back(g); }
    void addPi(EcoGate* g) { _PIList.push_back(g); }
    void addGate(EcoGate* g) { if(!_gateName2Gate.count(g->getGateName())) _gateName2Gate[g->getGateName()] = g; GateVec.push_back(g); if(g->isPi()) addPi(g); }
    

    // get functions
    EcoGate* getGateByName(const string& name);
    unordered_set<EcoGate*> getGateByAbcNode(Abc_Obj_t* pObj) { if(!_abcObj2EcoGate.count(Abc_ObjRegular(pObj))) return {}; return _abcObj2EcoGate.at(Abc_ObjRegular(pObj)); }
    EcoGate* getConst0Gate();
    EcoGate* getConst1Gate();
    EcoGate* getPoByName(const string& name);
    EcoGate* getGate(unsigned id) { return GateVec.at(id); }
    EcoGate* getPi(unsigned id) { return _PIList.at(id); }
    EcoGate* getPo(unsigned id) { return _POList.at(id); }
    unsigned getNumGates() { return GateVec.size(); }
    unsigned getNumPis() { return _PIList.size(); }
    unsigned getNumPos() { return _POList.size(); }
    Abc_Ntk_t* getAbcNtk() { return _pAbcNtk; }
    vector<EcoCut*> getGateCuts(EcoGate* g) { if(!_gate2Cuts.count(g)) return {}; return _gate2Cuts.at(g); }
    const vector<CirGate*> getAigDfsList() const { CirMgr* pCirMgr = cirV->getEcoCirV(); return pCirMgr->_dfsList; }

    // simulation functions
    void simOnPats(size_t** pats, unsigned nPats);

    // write the ntk as a verilog file
    void writeNtkVerilog(const string& fileName);
    
  private:
    unordered_set<string> gateTypeStrings = {"and", "or", "nand", "nor", "not", "buf", "xor", "xnor"};
    // map that records gate name 2 gates
    unordered_map<string, EcoGate*> _gateName2Gate;
    unordered_map<string, EcoGate*> _poName2PoGate;

    // used for fraig
    unordered_map<Abc_Obj_t*, unordered_set<EcoGate*>> _abcObj2EcoGate;

    // cut members
    unordered_map<EcoGate*, vector<EcoCut*>> _gate2Cuts;
    unordered_set<string> _enumeratedSignatures;

    // cut signature functions
    bool computeAndInsertSigature(EcoCut* pCut);
    bool checkSignatureExists();

    // PI list
    vector<EcoGate*> _PIList;
    // PO list
    vector<EcoGate*> _POList;
    // Gate List
    vector<EcoGate*> GateVec;
    // store the AIG version of the ntk
    EcoCir* cirV;
    Abc_Ntk_t* _pAbcNtk;
};




    


// end of name space gv::cir
}}


#endif