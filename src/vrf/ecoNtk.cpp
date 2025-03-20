#ifndef ECO_NTK_CPP
#define ECO_NTK_CPP

#include "ecoNtk.h"
#include "base/abc/abc.h"
#include <iostream>
#include <string>
#include <cassert>

unsigned gv::cir::EcoGate::_globalTravFlag = 0;

namespace gv {
namespace cir {
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

EcoGate::EcoGate(string gateType, string gateName) : ecoGateVComp(false), ecoGateV(nullptr), _pAbcNode(nullptr) {
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

void
EcoGate::reportGate() {
  cout << getGateTypeName() << " " << getGateName() << " ";
  for(const auto& fanin : _fanins) {
    cout << fanin->getGateName() << " ";
  }
  cout << endl;
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
    ofstream fout("/home/yenlu_mepu/gv/tmp.v");
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
        _gateName2Gate[gateName] = gate;
        GateVec.push_back(gate);
      }
      else if(firstTok == "assign") { // deal with assigns
        assert(items.size() >= 3);
        vector<string> nets = getAssignNets(buf);
        string gateName = nets.at(0);
        EcoGate* gate = new EcoGate("buf", gateName);
        gate->_faninNames.push_back(nets.at(1));
        _gateName2Gate[gateName] = gate;
        GateVec.push_back(gate);
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
      cirGate = cirV->getGate(0); // get the const 0 gate of cirV
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
        cout << "err " << faninName << endl;
      }
      assert(fanin != nullptr);
      gate->_fanins.push_back(fanin);
    }
  }
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
            _PIList.push_back(gate);
            _gateName2Gate[PI] = gate;
            GateVec.push_back(gate);
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
  EcoGate* const0 = new EcoGate("const0", "1'b0");
  _gateName2Gate["1'b0"] = const0;
  _PIList.push_back(const0);
  EcoGate* const1 = new EcoGate("const1", "1'b1");
  _gateName2Gate["1'b1"] = const1;
  _PIList.push_back(const1);
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

}
}
#endif
