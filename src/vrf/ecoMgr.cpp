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
  _nSim = 10000;

  setOldDesignName(oldDesignName);
  setNewDesignName(newDesignName);

  // read designs
  readDesigns(oldDesignName, newDesignName);

  // do fraig
  doFraig();

  // dp random sim
  doRandomSim();

  // do matching
  doMatching(6); // enumerate 6-feasible cuts

  // generate patch
  genPatch("patch.v");

  // do recycle
  // doRecycle();
}

// read input designs
void
EcoMgr::readDesigns(const string& oldDesignName, const string& newDesignName) {
  _oldNtk->readNtkFile(oldDesignName);
  _newNtk->readNtkFile(newDesignName);
  for(size_t i=0; i<_oldNtk->getNumGates(); i++)
    _oldNtk->getGate(i)->setOld(gv::cir::EcoGate::ECO_OLD_NTK);
  for(size_t i=0; i<_oldNtk->getNumPos(); i++) {
    _oldNtk->getPo(i)->setOld(gv::cir::EcoGate::ECO_OLD_NTK);
  }
  for(size_t i=0; i<_newNtk->getNumGates(); i++)
    _newNtk->getGate(i)->setOld(gv::cir::EcoGate::ECO_NEW_NTK);
  for(size_t i=0; i<_newNtk->getNumPos(); i++)
    _newNtk->getPo(i)->setOld(gv::cir::EcoGate::ECO_NEW_NTK);
}

void
EcoMgr::Net2PO( Abc_Ntk_t* pNtk)
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
      // cout << Abc_ObjType(pNode) << " ooo " << gate->getGateFullName() << " " << Abc_ObjId(pNode) << endl;
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

  // write the eq class in a file for debug use
  ofstream f;
  f.open("./.eq_class.txt", ios::out);
  f << "Eq class" << endl;
  for(const auto&[pEqNode, oldGates] : oldEqClass) {
    if(!newEqClass.count(pEqNode)) continue;
    auto& newGates = newEqClass.at(pEqNode);
    for(const auto&[oldGate, oldGateComp] : oldGates) {
      f << oldGate->getGateFullName() << (oldGateComp ? "(inv)" : "(pos)") << " ";
      oldGate->setIsMerged();
      
      for(const auto&[newGate, newGateComp] : newGates) {
        newGate->setIsMerged();
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
    for(const auto&[newGate, newGateComp] : newGates)
      f << newGate->getGateFullName() << (newGateComp ? "(inv)" : "(pos)") << " ";
    f << endl;
  }
  f.close();

  auto oldPiEqClass = _oldNtk->getPiEqClass();
  auto newPiEqClass = _newNtk->getPiEqClass();
  for(size_t i=0; i<_oldNtk->getNumPis(); i++) {
    auto oldGates = oldPiEqClass.at(i);
    auto newGates = newPiEqClass.at(i);
    for(auto& oldGate : oldGates) {
      oldGate->setIsMerged();
      auto oldComp = oldGate->getAigNodeInv();
      f << oldGate->getGateFullName() << (oldComp ? "(inv)" : "(pos)") << " ";
      for(auto& newGate : newGates) {
        newGate->setIsMerged();
        auto newComp = newGate->getAigNodeInv();
        if((oldComp ^ newComp) == 0) {
          _mergeTable[oldGate].insert(newGate);
          _mergeTable[newGate].insert(oldGate);
        }
        else {
          _invMergeTable[oldGate].insert(newGate);
          _invMergeTable[newGate].insert(oldGate);
        }
      }
    }
    for(auto& newGate : newGates) {
      auto newComp = newGate->getAigNodeInv();
      f << newGate->getGateFullName() << (newComp ? "(inv)" : "(pos)") << " ";
    }
    f << endl;
  }

  // merge constant and PIs
  for(size_t i=0; i<_oldNtk->getNumPis(); i++) {
    auto oldGate = _oldNtk->getPi(i);
    auto newGate = _newNtk->getPi(i);
    
    _mergeTable[oldGate].insert(newGate);
    _mergeTable[newGate].insert(oldGate);
  }

  markMergeFrontier();
}

// mark the gates that are strictly under the merge frontier
void
EcoMgr::markMergeFrontier() {
  unordered_set<gv::cir::EcoGate*> underMergedFrontierSet;
  unordered_set<gv::cir::EcoGate*> removeUnderMergedFrontierSet;
  for(unsigned i=0; i<_oldNtk->getNumGates(); ++i) {
        auto g = _oldNtk->getGate(i);
        
        if(isMerged(g))
            underMergedFrontierSet.insert(g);
        else
            continue;
        
        queue<gv::cir::EcoGate*> q;
        gv::cir::EcoGate::setGlobalTrav();
        q.push(g);

        // mark merged gate's fanin cone as under merged frontier
        while(!q.empty()) {
            auto cur = q.front();
            q.pop();
            if(cur->isGlobalTrav()) continue;
            cur->setToGlobalTrav();

            underMergedFrontierSet.insert(cur);

            // traverse its children
            for(unsigned j=0; j<cur->getNumFanins(); ++j) {
                auto fanin = cur->getFanin(j);
                if(!fanin->isGlobalTrav())
                    q.push(fanin);
            }
        }
    }

    removeUnderMergedFrontierSet = underMergedFrontierSet;
    gv::cir::EcoGate::setGlobalTrav();
    for(unsigned i=0; i<_oldNtk->getNumPos(); ++i) {
      auto po = _oldNtk->getPo(i)->getFanin(0);
      queue<gv::cir::EcoGate*> q;
      q.push(po);

      // mark merged gate's fanin cone as under merged frontier
      while(!q.empty()) {
        auto cur = q.front();
        q.pop();
        if(cur->isGlobalTrav()) continue;
        cur->setToGlobalTrav();


        // traverse its children
        for(unsigned j=0; j<cur->getNumFanins(); ++j) {
          auto fanin = cur->getFanin(j);
          if(!fanin->isGlobalTrav())
            q.push(fanin);
          if(underMergedFrontierSet.count(cur) && cur->getNumFanins() > 1) {
            removeUnderMergedFrontierSet.erase(fanin);
          }
        }
      }
    }

    for(const auto& g : underMergedFrontierSet) {
      if(!removeUnderMergedFrontierSet.count(g)) {
        addUnderMergeFrontierSet(g);
      }
    }
}

// get the aig of merged aig and pole
pair<gv::cir::CirGate*, bool>
EcoMgr::getMergedAig(gv::cir::EcoGate* g) {
  assert(isMerged(g)); // make sure that g is a merged gate
  gv::cir::CirGate* mergedAig;
  bool pole;
  auto mergedGates = getMergedGates(g);
  if(!mergedGates.empty()) {
    gv::cir::EcoGate* mergedGate = *(mergedGates.begin());
    mergedAig = mergedGate->getAigNode();
    pole = mergedGate->getAigNodeInv();
  }
  else {
    auto invMergedGates = getInvMergedGates(g);
    assert(!invMergedGates.empty());
    gv::cir::EcoGate* invMergedGate = *(invMergedGates.begin());
    mergedAig = invMergedGate->getAigNode();
    pole = (invMergedGate->getAigNodeInv() ^ true);
  }
  return {mergedAig, pole};
}

void
EcoMgr::doMatching(unsigned kFeassible) {
  
  // build cut hashing table
  _pNpnHash = new EcoNPNHash(2, 4); // we compute the npn cut hash from 2 <= cut size <= 4
  _pNpnHash->computeNpnHash();
  
  // enumerate k-feasible cuts
  _oldNtk->enumerateCuts(kFeassible, this);
  _newNtk->enumerateCuts(kFeassible, this);

  // output side matching (along with fault analysis)
  doOutputSideMatching();

  // compute the score of each cut and match the remaining cuts
  // sortCandCutsByScore();

  // do floating gate recycle
  doRecycle();
  


  // report RP pairs
  // reportRPPair();
  // for(auto[oldGate, newGate] : _RPPair) {
  //   cout << oldGate->getGateFullName() << " -> " << newGate->getGateFullName() << endl;
  // }
}

// do floating gate recycle
void
EcoMgr::doRecycle() {
  // 1. Collect the floating gates in the old circuit
  collectFloatingGates();

  // 2. Map the floating gates to the sub-circuit in the patch
  matchFloatingGates();

}


// end of namespace gv::eco
}}
#endif