#ifndef ECO_NTK_CPP
#define ECO_NTK_CPP

#include "ecoNtk.h"
#include "base/abc/abc.h"
#include <iostream>
#include <string>
#include <cassert>
#include <queue>

unsigned gv::cir::EcoGate::_globalTravFlag = 0;
unordered_set<string> gv::cir::EcoNtk::gateTypeStrings = {"and", "or", "nand", "nor", "not", "buf", "xor", "xnor"};

namespace gv {
namespace cir {
  extern void printBits(size_t tt);
// print the gate type name
string
EcoGate::getGateTypeName() {
  switch (_gateType) {
    case ECO_CONST_0_GATE:
      return "const0";
    case ECO_CONST_1_GATE:
      return "const1";
    case ECO_AND_GATE:
      return "and";
    case ECO_OR_GATE:
      return "or";
    case ECO_NAND_GATE:
      return "nand";
    case ECO_NOR_GATE:
      return "nor";
    case ECO_XOR_GATE:
      return "xor";
    case ECO_XNOR_GATE:
      return "xnor";
    case ECO_BUF_GATE:
      return "buf";
    case ECO_NOT_GATE:
      return "not";
    case ECO_PI_GATE:
      return "pi";
    case ECO_PO_GATE:
      return "po";
    default:
      return "none";
  }
}

EcoGate::EcoGate(string gateType, string gateName) : ecoGateVComp(false), ecoGateV(nullptr), _pAbcNode(nullptr), _gateNtk(ECO_NONE_NTK), _id(-1) {
  _travFlag = 0;
  _gateName = gateName;
  if(gateType == "const0")
    _gateType = ECO_CONST_0_GATE;
  else if(gateType == "const1")
    _gateType = ECO_CONST_1_GATE;
  else if(gateType == "and")
    _gateType = ECO_AND_GATE;
  else if(gateType == "or")
    _gateType = ECO_OR_GATE;
  else if(gateType == "nand")
    _gateType = ECO_NAND_GATE;
  else if(gateType == "nor")
    _gateType = ECO_NOR_GATE;
  else if(gateType == "xor")
    _gateType = ECO_XOR_GATE;
  else if(gateType == "xnor")
    _gateType = ECO_XNOR_GATE;
  else if(gateType == "buf")
    _gateType = ECO_BUF_GATE;
  else if(gateType == "not")
    _gateType = ECO_NOT_GATE;
  else if(gateType == "pi")
    _gateType = ECO_PI_GATE;
  else if(gateType == "po")
    _gateType = ECO_PO_GATE;
}

EcoGate::~EcoGate() {
  _simVals.clear();
}

// add the given gate as the PI of the circuit
void
EcoNtk::addPi(EcoGate* g) {
  // check if the gate is already in the PI list of the ntk (checked by gate name)
  // if it already exists, we don't add it
  for(const auto& pi : _PIList) {
    if(pi->getGateName() == g->getGateName())
      return;
  }
  _PIList.push_back(g);
}

// add the given gate as the PO of the circuit
void 
EcoNtk::addPo(EcoGate* g) {
  // check if the gate is already in the PO list of the ntk (checked by gate name)
  // if it already exists, we don't add it
  for(const auto& po : _POList) {
    if(po->getGateName() == g->getGateName())
      return;
    assert(_poName2PoGate.count(po->getGateName()));
  }
  assert(_poName2PoGate.size() == _POList.size());
  
  if(!_poName2PoGate.count(g->getGateName())) {
    _poName2PoGate[g->getGateName()] = g;
  }
  else {
    cout << "po name " << g->getGateName() << "already exists!!!" << endl;
    cout << "damn " << _poName2PoGate.at(g->getGateName()) << " " << g << endl;
    assert(0);
  }
  _POList.push_back(g);
}

// add the given gate as a gate in EcoNtk
void 
EcoNtk::addGate(EcoGate* g) {
  // add name mapping
  if(_gateName2Gate.count(g->getGateName())) {
    cout << "gate name : " << g->getGateName() << " already exists!!!" << endl;
  }
  assert(!_gateName2Gate.count(g->getGateName()));
    _gateName2Gate[g->getGateName()] = g;
  // if(g->getGateName() == )
    // cout << "ff " << endl;
  // push the gate
  g->setGateId(getNumGates());
  GateVec.push_back(g);
  
  // if the gate type is pi, add it to PI list
  if(g->isPi())
    addPi(g);
}

void
EcoGate::addFaninName(const string& name) {
  assert(!isPiOrConst());
  // if the fanin name already exists, no need to add
  if(find(_faninNames.begin(), _faninNames.end(), name) != _faninNames.end())
    return;
  if(getGateType() == ECO_BUF_GATE || getGateType() == ECO_NOT_GATE) {
    // if(getGateName() == "prim_out[6]_N") {
    //   cout << "prim_out[6]_N adding fanin name " << name << endl;
    //   assert(0)
    // }
    if(!_faninNames.empty()) {
      cout << "adding 2nd fanin for a buf/not gate " << getGateFullName() << " " << "dummy fanin " << name << endl;
      // assert(0);
    }
  }
  // add the fanin name
  _faninNames.push_back(name);
}


// report the gate name appendded with _O/_N if they are in old or new circuit
string
EcoGate::getGateFullName() {
  if(isPiOrConst() || _gateNtk > ECO_NEW_NTK)
    return getGateName();
  return getGateName() + ((_gateNtk == ECO_OLD_NTK) ? "_O" : "_N");
}

void
EcoGate::reportGate() {
  cout << getGateTypeName() << " " << getGateFullName() << " ";
  for(const auto& fanin : _fanins) {
    cout << fanin->getGateFullName() << " ";
  }
  cout << endl;
}

void
EcoNtk::setGateName2Gate(const string& newName, const string& oldName, EcoGate* g) {
  // cout << "bbb " << _gateName2Gate.at(oldName)->getGateFullName() << " " << g->getGateFullName() << endl;
  // cout << "ccc " << _gateName2Gate.at(oldName) << " " << g << endl;
  cout << "old name : " << oldName << " new name : " << newName << endl;
  // cout << "gate type : " << g->getGateTypeName() << " " << "name type : "<< _gateName2Gate.at(oldName)->getGateTypeName() << endl;
  // assert(0);
  // assert(_gateName2Gate.at(oldName) == g);
  g->setGateName(newName);
  _gateName2Gate.erase(oldName);
  _gateName2Gate[newName] = g;
}

// get the eco gate by name
EcoGate*
EcoNtk::getGateByName(const string& name) {
  if(!_gateName2Gate.count(name)) return nullptr;
  return _gateName2Gate.at(name);
}

// get the const0 gate
EcoGate*
EcoNtk::getConst0Gate() {
  return getGateByName("1'b0");
}

// get the const1 gate
EcoGate*
EcoNtk::getConst1Gate() {
  return getGateByName("1'b1");
}

// get po by po name
EcoGate*
EcoNtk::getPoByName(const string& name) {
  if(!_poName2PoGate.count(name)) return nullptr;
  return _poName2PoGate.at(name);
}

// split the line by white spaces
vector<string> splitLine(const string& line) {
  vector<string> ret;
  string buf;
  for(size_t i = 0, n = line.size(); i < n; i++) {
    if(line[i] != ' ')
      buf.push_back(line[i]);
    if(line[i] == ' ' || i == n - 1) {
      if(!buf.empty()) {
        ret.push_back(buf);
        buf.clear();
      }
    }
  }
  return ret;
}

string stripSpecialTok(const string& str) {
  string ret;
  unordered_set<char> specialTok = {'(', ')', ';', ',', '{', '}'};
  size_t i = 0;
  while(i < str.size()) {
    if(!specialTok.count(str.at(i))) break;
  }
  for( ; i < str.size(); i++) {
    if(specialTok.count(str.at(i))) break;
    ret.push_back(str.at(i));
  }
  assert(!ret.empty());
  return ret;
}

// get the net names in the primitive nets
vector<string> getGateNets(const string& line) {
  vector<string> ret;
  string buf;
  bool flag = false;
  for(size_t i = 0, n = line.size(); i < n; i++) {
    if(line[i] == '(')
      flag = true;
    
    if(flag && line[i] != ' ' && line[i] != ',' && line[i] != '(' && line[i] != ')')
      buf.push_back(line[i]);
    else if(flag) {
      if(!buf.empty()) {
        ret.push_back(buf);
        buf.clear();
      }
    }
    if(line[i] == ')')
      break;
  }
  return ret;
}

// get the net names in the assign line
vector<string> getAssignNets(const string& line) {
  vector<string> ret;
  string buf;
  // bool flag = false;
  for(size_t i = 0, n = line.size(); i < n; i++) {
    if(line[i] != ' ' && line[i] != ',' && line[i] != '(' && line[i] != ')' && line[i] != '=' && line[i] != ';')
      buf.push_back(line[i]);
    else {
      if(!buf.empty() && buf != "assign")
        ret.push_back(buf);
      buf.clear();
    }
  }
  assert(ret.size() == 2);
  return ret;
}

vector<string> getWiresPorts(const string& line) {
  vector<string> ret;
  string buf;
  int lsb = -1, msb = -1;
  bool flag = false;
  for(size_t i = 0, n = line.size(); i < n; i++) {
    if(line[i] == ' ')
      flag = true;
    
    if(flag && line[i] != ' ' && line[i] != ',' && line[i] != ';')
      buf.push_back(line[i]);
    else if(flag) {
      if(!buf.empty()) {
        if(buf[0] == '[') {
          string lsbStr, msbStr;
          bool afterColon = false;
          for(size_t j=0; j<buf.size(); j++) {
            if(buf[j] == ':') {
              afterColon = true;
              continue;
            }
            else if(buf[j] == '[' || buf[j] == ']')
              continue;
            if(!afterColon)
              lsbStr.push_back(buf[j]);
            else
              msbStr.push_back(buf[j]);
          }
          
          lsb = stoi(lsbStr);
          msb = stoi(msbStr);
          if(lsb > msb)
            swap(lsb, msb);
        }
        else if(buf != "wire" && buf != "input" && buf != "output"){
          if(lsb >= 0) {
            for(int i=lsb; i<=msb; i++) {
              ret.push_back(buf + "[" + to_string(i) + "]");
            }
          }
          else
            ret.push_back(buf);
        }
        buf.clear();
      }
    }
  }
  return ret;
}

// rewrite the design file (handle capital chars and gate name stuff)
void
EcoNtk::rewriteDesign(const string& dir) {
  ifstream file(dir);
    assert(file.is_open());
    ofstream fout("./tmp.v");
    string buf;
    unordered_set<string> wires;
    while (getline(file, buf))
    {
      // splitting the items in the line by white space
      vector<string> items = splitLine(buf);
      
      // if the line is not empty
      if(!items.empty()) {
        string firstTok = items.at(0);
        // assure that the first token is not capital
        if(firstTok.at(0) >= 'A' && firstTok.at(0) <= 'Z') {
          for(size_t i=0; i<firstTok.size(); i++)
            firstTok[i] = static_cast<char>(firstTok[i] + ('a' - 'A'));
        }
        // if it is a primitive
        if(gateTypeStrings.count(firstTok)) {
          assert(items.size() > 1);
          vector<string> nets = getGateNets(buf);
          for(const auto& net : nets) {
            if(!wires.count(net)) {
              fout << "wire " << net << ";" << endl;
              wires.insert(net);
            }
          }
          fout << firstTok << " ";
          // ignore gate name
          string gateNameItem = items.at(1);
          bool flag = false;
          for(size_t i=0; i<gateNameItem.size(); i++) {
            if(gateNameItem.at(i) == '(')
              flag = true;
            if(!flag) continue;
            fout << gateNameItem.at(i);
          }
          fout << " ";
          // print the remaining items
          for(size_t i=2; i<items.size(); i++) {
            fout << items.at(i) << " ";
          }
          fout << endl;
        }
        else {
          if(firstTok == "wire") {
            while(1) {
              vector<string> lineWires = getWiresPorts(buf);
              for(const auto& w : lineWires)
                wires.insert(w);
              while(buf[buf.size()-1] == ' ')
                buf.pop_back();
              if(buf[buf.size()-1] == ';')
                break;
              fout << buf << endl;
              getline(file, buf);
            }
          }
          fout << buf << endl;
        }
      }
    }

}



// parse the primitve information and store them
void
EcoNtk::parsePrimitiveGates(const string& dir) {
  ifstream file(dir);
  assert(file.is_open());
  string buf;
  while (getline(file, buf))
  {
    // splitting the items in the line by white space
    vector<string> items = splitLine(buf);
    
    // if the line is not empty
    if(!items.empty()) {
      string firstTok = items.at(0);
      // if it is a primitive
      if(gateTypeStrings.count(firstTok)) {
        string gateType = firstTok;
        vector<string> nets = getGateNets(buf);
        string gateName = nets.at(0);
        EcoGate* gate = new EcoGate(gateType, gateName);
        for(size_t i=1; i<nets.size(); i++)
          gate->_faninNames.push_back(nets.at(i));
        addGate(gate);
      }
      else if(firstTok == "assign") { // deal with assigns
        assert(items.size() >= 3);
        vector<string> nets = getAssignNets(buf);
        string gateName = nets.at(0);
        EcoGate* gate = new EcoGate("buf", gateName);
        gate->_faninNames.push_back(nets.at(1));
        addGate(gate);
      }
    }
  }
  file.close();
}

void
EcoNtk::abcReadFile() {
  // abc read file parameters
  int fAllNodes = 1;
  Abc_Obj_t * pObj;
  Abc_Obj_t * pNode;
  int i;

  Abc_Ntk_t* pNtk = Io_Read( "/home/yenlu_mepu/gv/tmp.v", IO_FILE_VERILOG, 0, 0 );
  assert(pNtk && Abc_NtkCheck(pNtk)); // check that the read circuit is OK
  
  Abc_Ntk_t* pNtkStrash = Abc_NtkStrash( pNtk, fAllNodes, !fAllNodes, 0 ); // strash the circuit
  _pAbcNtk = pNtkStrash;

  // new the internal EcoCir
  cirV->readCirFromAbcNtk(pNtkStrash);
  
  Abc_NtkForEachObj( pNtk, pNode, i )
  {
    string objName = Abc_ObjName( pNode );
    
    if(!pNode->pCopy)
      continue;
    
    CirGate* cirGate;
    EcoGate* ecoGate;
    if(Abc_ObjType(Abc_ObjRegular(pNode->pCopy)) != ABC_OBJ_CONST1) {
      cirGate = cirV->getGate(Abc_ObjId(Abc_ObjRegular(pNode->pCopy)));
      ecoGate = getGateByName(objName);
      ecoGate->ecoGateVComp = Abc_ObjIsComplement(pNode->pCopy) ? true : false;
    }
    else { // handle the const gate case (const node name in abc is not the same as in my data structure )
      cirGate = cirV->getEcoCirV()->_const0; // get the const 0 gate of cirV
      ecoGate = Abc_ObjIsComplement(pNode->pCopy) ? getConst0Gate() : getConst1Gate();
      ecoGate->ecoGateVComp = Abc_ObjIsComplement(pNode->pCopy) ? false : true; // since abc's const gate is const1 and ours is const0
    }
    ecoGate->ecoGateV = cirGate;
    ecoGate->_pAbcNode = pNode->pCopy;
    _abcObj2EcoGate[Abc_ObjRegular(pNode->pCopy)].insert(ecoGate);
  }
  for(size_t i=0; i<getNumPos(); i++)
    getPo(i)->ecoGateV = getPo(i)->getFanin(0)->ecoGateV;
}

void
EcoNtk::genConnection() {
  for(auto& gate : GateVec) {
    for(auto& faninName : gate->_faninNames) {
      EcoGate* fanin = getGateByName(faninName);
      if(fanin == nullptr) {
        cout << "err " << gate->getGateFullName() << " not found fanin " << faninName << endl;
      }
      assert(fanin != nullptr);
      gate->_fanins.push_back(fanin);
    }
  }
  sortGatesInTopoOrder();
}

// sort the gates in GateVec by topological order
void
EcoNtk::sortGatesInTopoOrder() {
  vector<int> inOrder(getNumGates(), 0); // record the number of fanin of the gate
  vector<vector<EcoGate*>> fanouts(getNumGates());
  vector<EcoGate*> sortedGateVec;
  queue<EcoGate*> q;

  // check the fanins of each gate
  for(unsigned i=0; i<getNumGates(); ++i) {
    auto gate = getGate(i);
    inOrder[i] = gate->getNumFanins();

    for(unsigned j=0; j<gate->getNumFanins(); ++j) {
      auto fanin = gate->getFanin(j);
      fanouts[fanin->getGateId()].push_back(gate);
    }

    if(inOrder[i] == 0) {
      q.push(gate);
    }
  }

  while(!q.empty()) {
    auto cur = q.front();
    q.pop();

    sortedGateVec.push_back(cur);
    for(auto fanout : fanouts.at(cur->getGateId())) {
      inOrder.at(fanout->getGateId())--;
      if(inOrder.at(fanout->getGateId()) == 0)
        q.push(fanout);
    }
  }
  for(unsigned i=0; i<getNumGates(); ++i) {
    auto g = getGate(i);
    if(inOrder[g->getGateId()] != 0) {
      cout << g->getGateName() << " in order not 0 but " << inOrder[g->getGateId()] << " i = " << i << endl;
    }
  }
  assert(GateVec.size() == sortedGateVec.size());

  // assign the sorted container back to GateVec
  GateVec = sortedGateVec;
}

void
EcoNtk::parsePO(const string& dir) {
  ifstream file(dir);
  assert(file.is_open());
  string buf;
  while (getline(file, buf))
  {
    // splitting the items in the line by white space
    vector<string> items = splitLine(buf);
    
    // if the line is not empty
    if(!items.empty()) {
      string firstTok = items.at(0);
      if(firstTok == "output") {
        while(1) {
          vector<string> POs = getWiresPorts(buf);
          for(const auto& PO : POs) {
            EcoGate* gate = new EcoGate("po", PO);
            gate->_faninNames.push_back(PO);
            gate->_fanins.push_back(_gateName2Gate.at(PO));
            _POList.push_back(gate);
            _poName2PoGate[PO] = gate;
          }
          while(buf[buf.size()-1] == ' ')
            buf.pop_back();
          if(buf[buf.size()-1] == ';')
            break;
          getline(file, buf);
        }
      }
    }
  }
  file.close();
}

void
EcoNtk::parsePI(const string& dir) {
  ifstream file(dir);
  assert(file.is_open());
  string buf;
  while (getline(file, buf))
  {
    // splitting the items in the line by white space
    vector<string> items = splitLine(buf);
    
    // if the line is not empty
    if(!items.empty()) {
      string firstTok = items.at(0);
      if(firstTok == "input") {
        while(1) {
          vector<string> PIs = getWiresPorts(buf);
          for(const auto& PI : PIs) {
            EcoGate* gate = new EcoGate("pi", PI);
            addGate(gate);
          }
          while(buf[buf.size()-1] == ' ')
            buf.pop_back();
          if(buf[buf.size()-1] == ';')
            break;
          getline(file, buf);
        }
      }
    }
  }
  file.close();
  if(!getGateByName("1'b0")) {
    EcoGate* const0 = new EcoGate("const0", "1'b0");
    addGate(const0);
  }
  // _PIList.push_back(const0);
  if(!getGateByName("1'b1")) {
    EcoGate* const1 = new EcoGate("const1", "1'b1");
    addGate(const1);
  }
  // _PIList.push_back(const1);
  sort(_PIList.begin(), _PIList.end(), [](EcoGate* g1, EcoGate* g2) {
    return (g1->getGateName() < g2->getGateName());
  });
}

void
EcoNtk::readNtkFile(const string& dir) {
  parsePI(dir);
  parsePrimitiveGates(dir);
  parsePO(dir);
  genConnection();
  rewriteDesign(dir); // rewrite the design format so that abc can read it
  abcReadFile(); // read the rewrited file using abc
  // for(auto& gate : GateVec) {
  //   gate->reportGate();
  // }
  // for(auto& gate : _POList) {
  //   gate->reportGate();
  // }
}

void
EcoNtk::simOnPats(size_t** pats, unsigned nPats) {
  for(unsigned i=0; i<nPats; ++i) {
    for(unsigned j=0; j<getNumPis(); ++j) {
      auto pi = getPi(j);
      auto aig = pi->getAigNode();
      assert(!pi->getAigNodeInv());
      aig->setPValue(pats[i][j]);
      // printBits(pats[i][j]);
    }
    // sim using the dfs list
    CirMgr* pCirMgr = cirV->getEcoCirV();
    for(const auto& g : pCirMgr->_dfsList) {
      g->pSim();
    }
    // add the sim value to the gate
    for(auto& g : GateVec) {
      if(g->getGateType() == EcoGate::ECO_CONST_0_GATE || g->getGateType() == EcoGate::ECO_CONST_1_GATE) continue;
      size_t val = g->getAigNode()->getPValue()();
      if(g->getAigNodeInv())
        val = ~val;
      g->_simVals.push_back(val);
    }
  }
}

void
EcoNtk::writeNtkVerilog(const string& fileName) {
  gv::cir::EcoGate::setGlobalTrav();
    // write the ntk
    ofstream f(fileName);
    f << "module top(";
    for(unsigned i=0; i<getNumPos(); ++i) {
        f << getPo(i)->getGateName();
        if(i < getNumPos() - 1 || getNumPis() > 0)
          f << ", ";
        if((i + 1) % 10 == 0)
          f << endl;
    }
    for(unsigned i=0; i<getNumPis(); ++i) {
        f << getPi(i)->getGateName();
        if(i < getNumPis() - 1)
            f << ", ";
    }
    f << ");" << endl << endl;
    f << "output ";
    for(unsigned i=0; i<getNumPos(); ++i) {
        if((i + 1) % 10 == 0) {
            f << ";" << endl;
            f << "output ";
        }
        f << getPo(i)->getGateName();
        if(i < getNumPos() - 1 && (i + 2) % 10 != 0)
            f << ", ";
    }
    f << ";" << endl << endl;
    
    if(getNumPis() > 0) {
      f << "input ";
      for(unsigned i=0; i<getNumPis(); ++i) {
        if((i + 1) % 10 == 0) {
            f << ";" << endl;
            f << "input ";
        }
        f << getPi(i)->getGateName();
        if(i < getNumPis() - 1 && (i + 2) % 10 != 0)
            f << ", ";
      }
      f << ";" << endl << endl;
    }

    f << "wire ";
    vector<string> wireVec;
    for(unsigned i=0; i<getNumGates(); ++i) {
      if(getGate(i)->isConstGate()) continue;
      wireVec.push_back(getGate(i)->getGateName());
    }
    for(unsigned i=0; i<wireVec.size(); ++i) {
        if((i + 1) % 10 == 0) {
            f << ";" << endl;
            f << "wire ";
        }
        f << wireVec.at(i);
        if(i < wireVec.size() - 1 && (i + 2) % 10 != 0)
            f << ", ";
    }
    f << ";" << endl << endl;
    for(unsigned i=0; i<getNumGates(); ++i) {
        auto g = getGate(i);
        if(g->isPiOrConst()) continue; // no need to write pi/const gates
        f << g->getGateTypeName() << " (" << g->getGateName() << ", ";
        for(unsigned j=0; j<g->getNumFanins(); ++j) {
            auto fanin = g->getFanin(j);
            f << fanin->getGateName();
            if(j < g->getNumFanins() - 1)
                f << ", ";
        }
        f << ");" << endl;
    }
    f << endl << "endmodule" << endl;
    f.close();
}


// compute the ntk cost based on cad contest 2021 problem A
int
EcoNtk::computeCadContestCost() {
  unordered_set<string> wireNames;
  int wireCost = 0, gateCost = 0;
  int totalCost = 0;

  for(const auto& g : GateVec) {
    wireNames.insert(g->getGateName());
    if(g->isPiOrConst()) continue;
    for(unsigned i=0; i<g->getNumFanins(); ++i) {
      auto fanin = g->getFanin(i);
      wireNames.insert(fanin->getGateName());
    }
    gateCost += (int)g->getNumFanins() - 2;
  }

  wireCost = wireNames.size();
  totalCost = (wireCost + gateCost);

  cout << "patch cost report :" << endl;
  cout << "wire cost : " << wireCost << endl;
  cout << "gate cost : " << gateCost << endl;
  cout << "-----------" << endl;
  cout << "total cost : " << totalCost << endl;
  cout << "-----------" << endl;

  return totalCost;
}

}
}
#endif
