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

  // do fraig
  doFraig();

  // do matching
  doMatching(6); // enumerate 6-feasible cuts

  // generate patch

}

// read input designs
void
EcoMgr::readDesigns(const string& oldDesignName, const string& newDesignName) {
  _oldNtk->readNtkFile(oldDesignName);
  _newNtk->readNtkFile(newDesignName);
  for(size_t i=0; i<_oldNtk->getNumGates(); i++)
    _oldNtk->getGate(i)->setOld(true);
  for(size_t i=0; i<_oldNtk->getNumPos(); i++) {
    _oldNtk->getPo(i)->setOld(true);
  }
  for(size_t i=0; i<_newNtk->getNumGates(); i++)
    _newNtk->getGate(i)->setOld(false);
  for(size_t i=0; i<_newNtk->getNumPos(); i++)
    _newNtk->getPo(i)->setOld(false);
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
    auto gates = _oldNtk->getGateByAbcNode(Abc_ObjRegular(pNode));
    for(const auto& gate : gates) {
      _oldNtk->setGateByAbcNode(Abc_ObjRegular(pNode->pCopy), gate);
    }
  }
// cout << "-----" << endl;
  Abc_NtkForEachObj( pNtkNew, pNode, i ) {
    auto gates = _newNtk->getGateByAbcNode(Abc_ObjRegular(pNode));
    for(const auto& gate : gates) {
      _newNtk->setGateByAbcNode(Abc_ObjRegular(pNode->pCopy), gate);
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
  unordered_map<Abc_Obj_t*, vector<pair<gv::cir::EcoGate*, bool>>> oldEqClass, newEqClass;
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
      oldEqClass[Abc_ObjRegular(Abc_ObjFanin0(pNode))].push_back({g, inv ^ g->getAigNodeInv()});

    for(auto g : newGates)
      newEqClass[Abc_ObjRegular(Abc_ObjFanin0(pNode))].push_back({g, inv ^ g->getAigNodeInv()});
  }
  // cout << "Eq class" << endl;
  for(const auto&[pEqNode, oldGates] : oldEqClass) {
    if(!newEqClass.count(pEqNode)) continue;
    auto& newGates = newEqClass.at(pEqNode);
    for(const auto&[oldGate, oldGateComp] : oldGates) {
      // cout << oldGate->getGateFullName() << (oldGateComp ? "(inv)" : "(pos)") << " ";
      for(const auto&[newGate, newGateComp] : newGates) {
        if((oldGateComp ^ newGateComp) == 0) {
          _mergeTable[oldGate].insert(newGate);
          _mergeTable[newGate].insert(oldGate);
        }
        else {
          _invMergeTable[oldGate].insert(newGate);
          _invMergeTable[newGate].insert(oldGate);
        }
      }
    }
    // for(const auto&[newGate, newGateComp] : newGates)
    //   cout << newGate->getGateFullName() << (newGateComp ? "(inv)" : "(pos)") << " ";
    // cout << endl;
  }

  // merge constant and PIs
  for(size_t i=0; i<_oldNtk->getNumPis(); i++) {
    auto oldGate = _oldNtk->getPi(i);
    auto newGate = _newNtk->getPi(i);
    
    _mergeTable[oldGate].insert(newGate);
    _mergeTable[newGate].insert(oldGate);
  }
  // for(size_t i=0; i<_oldNtk->getNumPos(); i++)
  //   dfs(_oldNtk->getPo(i));
  // for(size_t i=0; i<_newNtk->getNumPos(); i++)
  //   dfs(_newNtk->getPo(i));
}

void
EcoMgr::doMatching(unsigned kFeassible) {
  
  
  // build cut hashing table
  _pNpnHash = new EcoNPNHash(2, 4); // we compute the npn cut hash from 2 <= cut size <= 4
  _pNpnHash->computeNpnHash();
  
  // enumerate k-feasible cuts
  _oldNtk->enumerateCuts(kFeassible);
  _newNtk->enumerateCuts(kFeassible);

  // output side matching
  doOutputSideMatching();
}


// end of namespace gv::eco
}}
#endif