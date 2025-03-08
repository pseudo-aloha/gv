#ifndef ECO_NTK_H
#define ECO_NTK_H

#include "cirMgr.h"
#include "cirGate.h"
#include <vector>

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
namespace cir {

// Cuts in the EcoNtk
class EcoCut {
public:
  EcoCut() : _root(nullptr), _leaves({}) {}
  EcoCut(EcoGate* root) : _leaves({}) { _root = root; }
  EcoCut(EcoGate* root, vector<EcoGate*> leaves) { _root = root; _leaves = leaves; }
  ~EcoCut() { _root = nullptr; _leaves.clear(); }
  void setRoot(EcoGate* g) { _root = g; }
  void addLeaf(EcoGate* g) { _leaves.push_back(g); }
  void appendLeaves(vector<EcoGate*> leaves) { _leaves.insert(_leaves.end(), leaves.begin(), leaves.end()); }
  void popLeaves(int n) { for(size_t i=0; i<n; i++) _leaves.pop_back(); }
  EcoGate* const getRoot() { return _root; }
  vector<EcoGate*> const getLeaves() { return _leaves; }
  unsigned getCutSize() { return _leaves.size(); }
  void reportCut();
private:
  EcoGate* _root;
  vector<EcoGate*> _leaves;
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
    void enumerateCuts(unsigned k); // enumerate k-feasible cuts
    vector<EcoCut*> enumerateCutsRec(const unsigned& k, EcoGate* g); // enumerate k-feasible cuts

    // set functions
    void setGateByAbcNode(Abc_Obj_t* pObj, EcoGate* pEcoGate) { _abcObj2EcoGate[Abc_ObjRegular(pObj)].insert(pEcoGate); }


    // get functions
    EcoGate* getGateByName(const string& name);
    unordered_set<EcoGate*> getGateByAbcNode(Abc_Obj_t* pObj) { if(!_abcObj2EcoGate.count(Abc_ObjRegular(pObj))) return {}; return _abcObj2EcoGate.at(Abc_ObjRegular(pObj)); }
    EcoGate* getConst0Gate();
    EcoGate* getPoByName(const string& name);
    EcoGate* getGate(unsigned id) { return GateVec.at(id); }
    EcoGate* getPi(unsigned id) { return _PIList.at(id); }
    EcoGate* getPo(unsigned id) { return _POList.at(id); }
    unsigned getNumGates() { return GateVec.size(); }
    unsigned getNumPis() { return _PIList.size(); }
    unsigned getNumPos() { return _POList.size(); }
    Abc_Ntk_t* getAbcNtk() { return _pAbcNtk; }
  private:
    unordered_set<string> gateTypeStrings = {"and", "or", "nand", "nor", "not", "buf", "xor", "xnor"};
    // map that records gate name 2 gates
    unordered_map<string, EcoGate*> _gateName2Gate;
    unordered_map<string, EcoGate*> _poName2PoGate;

    // used for fraig
    unordered_map<Abc_Obj_t*, unordered_set<EcoGate*>> _abcObj2EcoGate;

    // get the cuts from gate
    unordered_map<EcoGate*, vector<EcoCut*>> _gate2Cuts;

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


class EcoGate {
  public:
    friend class EcoNtk;
    EcoGate(string gateType, string gateName);
    ~EcoGate();
    string getGateName() { return _gateName; }
    string getGateFullName() { return _gateName + (_isOld ? "_O" : "_N"); }
    string getGateTypeName();
    unsigned getGateType() { return _gateType; }
    unsigned getNumFanins() { return _fanins.size(); }
    EcoGate* getFanin(unsigned i) { if(i>=_fanins.size()) return nullptr; return _fanins.at(i); }
    void reportGate();
    bool isOld() { return _isOld; }
    void setOld(bool isOld) { _isOld = isOld; }
    bool getInv() { return ecoGateVComp; }

    // traversal things
    static void setGlobalTrav() { _globalTravFlag++; }
    void setGlobalTrav(unsigned i) { _globalTravFlag += i; }
    void setToGlobalTrav() { _travFlag = _globalTravFlag; }
    bool isGlobalTrav() { return (_travFlag == _globalTravFlag);}

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
    
    
    bool _isOld; // record gate belongs to old/new circuiit
};



// end of name space gv::cir
}}


#endif