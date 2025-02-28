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
}}

namespace gv {
namespace cir {
// inherited from CirMgr, modified to store some extra information for ECO usage
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

    // get functions
    EcoGate* getGateByName(const string& name);
    EcoGate* getPoByName(const string& name);
    Abc_Ntk_t* getAbcNtk() { return _pAbcNtk; }
  private:
    unordered_set<string> gateTypeStrings = {"and", "or", "nand", "nor", "not", "buf", "xor", "xnor"};
    unordered_map<string, EcoGate*> _gateName2Gate;
    unordered_map<string, EcoGate*> _poName2PoGate;
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
    string getGateTypeName();
    unsigned getGateType() { return _gateType; }
    unsigned getNumFanins() { return _fanins.size(); }
    void reportGate();
  private:
    // Use to represent the primitive gate
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
    unsigned _gateType;
    string _gateName;
    vector<EcoGate*> _fanins;
    vector<string> _faninNames;
    CirGate* ecoGateV;
    bool ecoGateVComp;
};



// end of name space gv::cir
}}


#endif