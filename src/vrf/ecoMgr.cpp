#ifndef ECO_MGR_CPP
#define ECO_MGR_CPP

#include "ecoMgr.h"
#include "ecoNtk.h"
#include "proof/fraig/fraig.h"

namespace gv {
namespace eco {
// top function to do ECO
void
EcoMgr::doEco(const string& oldDesignName, const string& newDesignName) {
  // read designs
  readDesigns(oldDesignName, newDesignName);
  doFraig();
}

// read input designs
void
EcoMgr::readDesigns(const string& oldDesignName, const string& newDesignName) {
  _oldNtk->readNtkFile(oldDesignName);
  cout << "------------------------------------" << endl;
  _newNtk->readNtkFile(newDesignName);
}

void Net2PO( Abc_Ntk_t* pNtk)
{
    int c;
    int fCheck, fBarBufs;
    char * pFileName;

	  fCheck = 1;
    fBarBufs = 0;

    //read verilog file
    //pull internal signals to po
    int i;
    Abc_Obj_t* 	pNode;
    Abc_Obj_t* 	pNet;
    char*		str;

    Abc_NtkForEachNode( pNtk, pNode, i )
    {
      str = Abc_ObjName( pNode );
      if ( Abc_ObjFaninNum(pNode) == 0 ) 
      {
        if ( !strcmp( Abc_ObjName(pNode), "1'b0" ) || !strcmp( Abc_ObjName(pNode), "1'b1") ) continue;
        if ( !Abc_ObjFanoutNum(pNode) ) continue;
        str = Abc_ObjName( Abc_ObjFanout0(pNode) );
      }

      pNet = Abc_NtkCreatePo( pNtk );
      Abc_ObjAddFanin( pNet, pNode );
      Abc_ObjAssignName( pNet, str, "_INT" );

      
      // printf( "create node after miter (old): %x\t%x\t(%s)\n", pNode, pNode->pCopy, str);
    }
}

// perform abc fraig on the circuits
void
EcoMgr::doFraig() {
  Abc_Ntk_t* pNtkOld = _oldNtk->getAbcNtk();
  Abc_Ntk_t* pNtkNew = _newNtk->getAbcNtk();


  // step 1 : build the miter
  Abc_Ntk_t* pNtkMiter = Abc_NtkMiter( pNtkOld, pNtkNew, 0, 0, 0, 0 );

  Abc_Obj_t * pNode;
  int i;
  Abc_NtkForEachObj( pNtkOld, pNode, i ) {
    // if(!pNode->pCopy) continue;
    auto gates = _oldNtk->getGateByAbcNode(Abc_ObjRegular(pNode));
    for(const auto& gate : gates) {
      // cout << pNode->pCopy << " name " << gate->getGateName() << endl;
      _oldNtk->setGateByAbcNode(Abc_ObjRegular(pNode->pCopy), gate);
    }
  }
// cout << "-----" << endl;
  Abc_NtkForEachObj( pNtkNew, pNode, i ) {
    // if(!pNode->pCopy) continue;
    auto gates = _newNtk->getGateByAbcNode(Abc_ObjRegular(pNode));
    for(const auto& gate : gates) {
      _newNtk->setGateByAbcNode(Abc_ObjRegular(pNode->pCopy), gate);
      // cout << pNode->pCopy << " name " << gate->getGateName() << endl;
    }
  }
  unordered_map<string, Abc_Obj_t *> name2Gate;
  Abc_NtkForEachObj( pNtkMiter, pNode, i ) {
    name2Gate[Abc_ObjName(pNode)] = pNode;
    // if(_oldNtk->getGateByAbcNode(pNode).empty() && _newNtk->getGateByAbcNode(pNode).empty()) continue;
    auto gates = _oldNtk->getGateByAbcNode(pNode);
    for(const auto& gate : gates) {
      _oldNtk->setGateByAbcNode(Abc_ObjRegular(pNode), gate);
    }
    

    gates = _newNtk->getGateByAbcNode(Abc_ObjRegular(pNode));
    for(const auto& gate : gates) {
      _newNtk->setGateByAbcNode(Abc_ObjRegular(pNode), gate);
    }
  }

  // step 2 : conduct fraig
  Fraig_Params_t Params, * pParams = &Params;
  int fAllNodes = 1;
  int fExdc = 0;
  memset( pParams, 0, sizeof(Fraig_Params_t) );
  pParams->nPatsRand  = 2048; // the number of words of random simulation info
  pParams->nPatsDyna  = 2048; // the number of words of dynamic simulation info
  pParams->nBTLimit   =  100; // the max number of backtracks to perform
  pParams->fFuncRed   =    1; // performs only one level hashing
  pParams->fFeedBack  =    1; // enables solver feedback
  pParams->fDist1Pats =    1; // enables distance-1 patterns
  pParams->fDoSparse  =    1; // performs equiv tests for sparse functions
  pParams->fChoicing  =    0; // enables recording structural choices
  pParams->fTryProve  =    0; // tries to solve the final miter
  pParams->fVerbose   =    0; // the verbosiness flag
  pParams->fVerboseP  =    0; // the verbosiness flag

  Net2PO(pNtkMiter);
  Abc_Ntk_t* pNtkMiterFraig = Abc_NtkFraig( pNtkMiter, &Params, fAllNodes, fExdc );


  // step 3 record the merge information
  unordered_map<Abc_Obj_t*, vector<pair<gv::cir::EcoGate*, bool>>> EqClass;
  Abc_NtkForEachPo( pNtkMiterFraig, pNode, i ) {
    string objName = Abc_ObjName(pNode);
    if(objName.substr(objName.length()-4, 4) != "_INT") continue;
    objName = objName.substr(0, objName.length()-4);
    assert(name2Gate.count(objName));
    auto pMiterNode = name2Gate.at(objName);
    auto oldGates = _oldNtk->getGateByAbcNode(Abc_ObjRegular(pMiterNode));
    auto newGates = _newNtk->getGateByAbcNode(Abc_ObjRegular(pMiterNode));
    bool inv = Abc_ObjFaninC0(pNode);
    for(auto g : oldGates)
      EqClass[Abc_ObjRegular(Abc_ObjFanin0(pNode))].push_back({g, inv ^ g->getInv()});

    for(auto g : newGates)
      EqClass[Abc_ObjRegular(Abc_ObjFanin0(pNode))].push_back({g, inv ^ g->getInv()});
  }
  cout << "Eq class" << endl;
  for(const auto&[_, gates] : EqClass) {
    if(gates.size() < 2) continue;
    for(const auto& [g, inv] : gates)
      cout << g->getGateName() << "(" << (inv ? "inv" : "pos") << ") ";
    cout << endl;
  }
}

// end of namespace gv::eco
}}
#endif