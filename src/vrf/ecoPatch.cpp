#ifndef ECO_PATCH_CPP
#define ECO_PATCH_CPP

#include "ecoMgr.h"
#include "proof/fraig/fraig.h"
#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"

namespace gv {
namespace eco {
    extern "C"{
        int Abc_NtkIvyProve( Abc_Ntk_t ** ppNtk, void * pPars );
        void Abc_NtkShow( Abc_Ntk_t * pNtk, int fGateNames, int fSeq, int fUseReverse, int fKeepDot );
        void Abc_TruthNpnTest( char * pFileName, int NpnType, int nVarNum, int fDumpRes, int fBinary, int fVerbose );
        void Abc_NtkCecFraig( Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2, int nSeconds, int fVerbose );
        void Abc_NtkCecSat( Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2, int nConfLimit, int nInsLimit );
        Abc_Ntk_t * Abc_NtkMulti( Abc_Ntk_t * pNtk, int nThresh, int nFaninMax, int fCnf, int fMulti, int fSimple, int fFactor );
      }
// decide the output side rewire
void
EcoMgr::decideOutputRewire() {
    // iterate over the rp pair to see the old net is fixed to which net most frequently
    unordered_set<gv::cir::EcoGate*> oldGates;
    unordered_map<gv::cir::EcoGate*, vector<pair<gv::cir::EcoGate*, array<unsigned, 2>>>> mp;
    for(unsigned i=0; i<_rpTable.size(); ++i) {
        for(auto[oldGate, pEcoRpInfo] : _rpTable.at(i)) {
            auto mappedGate = pEcoRpInfo->getMappedGate();
            auto inv = pEcoRpInfo->getMappedPole();
            if(!mp.count(oldGate)) {
                pair<gv::cir::EcoGate*, array<unsigned, 2>> mappedCountPair = {mappedGate, {0, 0}};
                mappedCountPair.second[inv]++;
                mp[oldGate].push_back(mappedCountPair);
            }
            else {
                bool foundMappedGate = false;
                for(unsigned j=0; j<mp.at(oldGate).size(); ++j) {
                    if(mp.at(oldGate).at(j).first == mappedGate) {
                        mp.at(oldGate).at(j).second[inv]++;
                        foundMappedGate = true;
                        break;
                    }
                }
                if(!foundMappedGate) {
                    pair<gv::cir::EcoGate*, array<unsigned, 2>> mappedCountPair = {mappedGate, {0, 0}};
                    mappedCountPair.second[inv]++;
                    mp[oldGate].push_back(mappedCountPair);
                }
                
            }
            cout << oldGate->getGateFullName() << " " << mappedGate->getGateFullName() << " inv : " << inv << endl;
        }
    }
    cout << "mapped counts : " << endl;
    for(const auto&[oldGate, vec] : mp) {
        cout << oldGate->getGateFullName() << " << ";
        for(const auto&[mappedGate, poleCounts] : vec) {
            if(poleCounts[0])
                cout << mappedGate->getGateFullName() << " : " << poleCounts[0] << " | ";
            if(poleCounts[1])
                cout << "~" << mappedGate->getGateFullName() << " : " << poleCounts[1] << " | ";
        }
        cout << endl;
    }

}

// generate the patch
void
EcoMgr::genPatch(const string& patchName) {
    gv::cir::EcoGate::setGlobalTrav();
    // handle the rewire for output cut matching
    // decide which net should use rewire
    decideOutputRewire();

    // collect patch circuit from each PO of the old circuit
    // if a gate that is reached is fixed to another gate in the new circuit, replace it with the corresponding gate
    for(unsigned i=0; i<_oldNtk->getNumPos(); ++i) {
        collectPatchGates(_oldNtk->getPo(i)->getFanin(0), nullptr, i, false);
        cout << "output " << _oldNtk->getPo(i)->getFanin(0)->getGateFullName() << endl;
    }

    unordered_set<string> patchPoNameSet;
    for(unsigned i=0; i<_patchNtk->getNumPos(); ++i) {
        auto patchPo = _patchNtk->getPo(i);
        patchPoNameSet.insert(patchPo->getGateName());
    }
    for(unsigned i=0; i<_patchNtk->getNumPis(); ++i) {
        auto patchGate = _patchNtk->getPi(i);
        if(patchPoNameSet.count(patchGate->getGateName())) {
            auto oldName = patchGate->getGateName();
            auto newName = oldName + "_in";
            patchGate->setGateName(newName);
            _patchNtk->setGateName2Gate(newName, oldName, patchGate);
        }
    }

    // sort the gates in topological order
    _patchNtk->sortGatesInTopoOrder();

    // write the patch ntk
    _patchNtk->writeNtkVerilog("patch.v"); // write the patch content
    _patchNtk->computeCadContestCost();

    // apply the patch to old circuit and do equivalence checking between it and the new circuit
    if(!applyNCheckPatch(patchName))
        cout << "patched circuit NEQ to new circuit!!!!" << endl;
    else
        cout << "patched circuit eq to new circuit!!!!" << endl;
}

void
EcoMgr::collectPatchGates(gv::cir::EcoGate* g, gv::cir::EcoGate* curPatchGate, const unsigned& ithPo, bool isEntry) {
    gv::cir::EcoGate* patchGate;
    if(g->isGlobalTrav()) {
        if(g->isInNewCircuit()) {
            if(curPatchGate) {
                if(isMerged(g)) {
                    // TODO : check merge gate thing
                    auto mergedGate = *getMergedGates(g).begin();
                    curPatchGate->addFanin(mergedGate);
                }
                if(_patchNtk->getGateByName(g->getGateFullName())) {
                    patchGate = _patchNtk->getGateByName(g->getGateFullName());
                    curPatchGate->addFanin(patchGate);
                }
            }
        }
        return;
    }
    g->setToGlobalTrav();
    
    // check if the gate is in old circuit to decide how we handle it
    if(g->isInOldCircuit()) {
        // 1. check if the gate is fixed to another gate
        auto rpPairsForIthPo = _rpTable.at(ithPo);
        if(rpPairsForIthPo.count(g)) {
            auto pEcoRpInfo = rpPairsForIthPo.at(g);
            gv::cir::EcoGate* mappedGate = pEcoRpInfo->getMappedGate();
            auto inv = pEcoRpInfo->getMappedPole();
            cout << "fixed to another gate : " << g->getGateFullName() << " " << mappedGate->getGateFullName() << endl;
            // add the entry gate as the po of the patch circuit
            // TODO :  1. check if there are multiple rewire
            //         2. decide buf/inv based on the mapping pole (done)
            //         3. May need to check for constant case (?)
            gv::cir::EcoGate* patchPoGate = new gv::cir::EcoGate("po", g->getGateName());
            string rewireGateType = inv ? "not" : "buf";
            gv::cir::EcoGate* patchRewireGate = new gv::cir::EcoGate(rewireGateType, g->getGateName());
            patchPoGate->addFanin(patchRewireGate);
            _patchNtk->addPo(patchPoGate);
            _patchNtk->addGate(patchRewireGate);
            
            collectPatchGates(mappedGate, patchRewireGate, ithPo, true);
            return;
        }
    }
    // new circuit case
    else {
        assert(g->isInNewCircuit()); // check that the circuit is in the new circuit XDDD, just in case
        
        // if the gate is merged, replace it by the merged gate
        if(isMerged(g)) {
            // TODO : check merge gate thing
            auto mergedGate = *getMergedGates(g).begin();
            patchGate = new gv::cir::EcoGate("pi", mergedGate->getGateName());
            _patchNtk->addGate(patchGate);
            // add the old gate into the PI of patch circuit
            
        }
        else {
            // check if the gate is already exist in the patch circuit
            if(!_patchNtk->getGateByName(g->getGateFullName())) {
                patchGate = new gv::cir::EcoGate(g->getGateTypeName(), g->getGateFullName());
                _patchNtk->addGate(patchGate);
            }
            else {
                patchGate = _patchNtk->getGateByName(g->getGateFullName());
            }
        }

        // add the patch gate as a fanin of the fanout of patch gate
        if(curPatchGate)
            curPatchGate->addFanin(patchGate);

        // if the gate itself is merged, no further traversal is needed
        if(isMerged(g)) return;
    }

    // if the gate is a PI gate or const gate, stop traversing
    if(g->isPiOrConst()) return;
    
    for(unsigned i=0; i<g->getNumFanins(); ++i) {
        auto faninGate = g->getFanin(i);
        collectPatchGates(faninGate, patchGate, ithPo, false);
    }
}

// check if the two circuits are eqivalent
// modified from Abc_NtkCecFraig
bool
isNtkEq(Abc_Ntk_t* pNtkOld, Abc_Ntk_t* pNtkNew) {
//   extern int Abc_NtkIvyProve( Abc_Ntk_t ** ppNtk, void * pPars );

  Prove_Params_t Params, * pParams = &Params;
  // pParams->fVerbose = 1;
  // build the miter
//   assert(Abc_NtkCheck(pNtkOld));
//   assert(Abc_NtkCheck(pNtkNew));
  Abc_Ntk_t* pNtkMiter = Abc_NtkMiter( pNtkOld, pNtkNew, 1, 0, 0, 0 );
//   assert(Abc_NtkCheck(pNtkMiter));
  // handle the trivial case
  int ret = Abc_NtkMiterIsConstant( pNtkMiter );
  if(ret >= 0)
    return ret;
  Abc_Ntk_t * pCnf;
  int nConfLimit = 0;
  int  nInsLimit  = 0;
  pCnf = Abc_NtkMulti( pNtkMiter, 0, 100, 1, 0, 0, 0 );
  // ret = Abc_NtkMiterSat( pCnf, nConfLimit, nInsLimit, 0, NULL, NULL );
  Prove_ParamsSetDefault( pParams );
  pParams->nItersMax = 5;
  ret = Abc_NtkIvyProve( &pNtkMiter, pParams );
  if ( ret == -1 )
      assert(0); // undecided after running out of resources
  else if ( ret >= 0 )
    return ret;
  Abc_NtkDelete(pNtkMiter);
}

// check that after applying patch to the old circuit, old circuit will become functionally equivalent to the new circuit
// if the equilvalence holds, return true; return false otherwise
bool
EcoMgr::applyNCheckPatch(const string& patchName) {
    gv::cir::EcoNtk* patchNtk = new gv::cir::EcoNtk;
    gv::cir::EcoNtk* patchedNtk = new gv::cir::EcoNtk;
    unordered_map<gv::cir::EcoGate*, gv::cir::EcoGate*> gateMap;
    unordered_set<string> patchNtkPoNameSet;
    unordered_set<string> oldNtkPiNames;

    // 1. read the patch NTK (since I don't want to assume patch is written properly)
    patchNtk->readNtkFile(patchName);

    // 3. add patch logic
    for(unsigned i=0; i<_oldNtk->getNumPis(); ++i)
        oldNtkPiNames.insert(_oldNtk->getPi(i)->getGateName());
    // (a) add patch gate
    // add gates
    for(unsigned i=0; i<patchNtk->getNumGates(); ++i) {
        auto patchGate = patchNtk->getGate(i);
        auto gateTypeName = patchGate->getGateTypeName();
        if(gateTypeName == "pi" && !oldNtkPiNames.count(patchGate->getGateName()))
            continue;
        gv::cir::EcoGate* patchedGate = new gv::cir::EcoGate(gateTypeName, patchGate->getGateName());
        
        // do the gate mapping
        gateMap[patchGate] = patchedGate;
        
        // add the gate into patched gate
        patchedNtk->addGate(patchedGate);

        // add fanin names
        for(unsigned j=0; j<patchGate->getNumFanins(); ++j) {
            auto patchFanin = patchGate->getFanin(j);
            patchedGate->addFaninName(patchFanin->getGateName());
        }

        // If a gate is PI, then also add it to the PI of patched circuit
        // if(patchGate->isPi())
        //     patchedNtk->addPi(patchedGate);
    }

    // 2. add the old circuit function
    
    // merged the pi's
    unordered_map<string, gv::cir::EcoGate*> patchedNtkPiNameGateMap;
    for(unsigned i=0; i<patchedNtk->getNumPis(); ++i) {
        auto pi = patchedNtk->getPi(i);
        patchedNtkPiNameGateMap[pi->getGateName()] = pi;
    }
    for(unsigned i=0; i<_oldNtk->getNumPis(); ++i) {
        auto oldNtkPi = _oldNtk->getPi(i);
        if(patchedNtkPiNameGateMap.count(oldNtkPi->getGateName()))
            gateMap[oldNtkPi] = patchedNtkPiNameGateMap.at(oldNtkPi->getGateName());
    }

    // (a) add po's
    for(unsigned i=0; i<_oldNtk->getNumPos(); ++i) {
        auto oldPo = _oldNtk->getPo(i);
        gv::cir::EcoGate* patchedPo = new gv::cir::EcoGate(oldPo->getGateTypeName(), oldPo->getGateName());
        patchedNtk->addPo(patchedPo);
        gateMap[oldPo] = patchedPo;
    }

    // (b) add old circuit logics and also add patched circuit pi if it is PI in the old circuit
    
    // collect the outputs' name of patch circuit
    for(unsigned i=0; i<_patchNtk->getNumPos(); ++i) {
        auto patchPo = _patchNtk->getPo(i);
        patchNtkPoNameSet.insert(patchPo->getGateName());
    }
    
    for(unsigned i=0; i<_oldNtk->getNumGates(); ++i) {
        auto oldGate = _oldNtk->getGate(i);

        if(gateMap.count(oldGate)) continue;
        auto patchedGateName = oldGate->getGateName();

        // check if the gate exists in the output of the patch circuit
        if(patchNtkPoNameSet.count(patchedGateName))
            patchedGateName += "_in";

        gv::cir::EcoGate* patchedGate = new gv::cir::EcoGate(oldGate->getGateTypeName(), patchedGateName);
        
        // do the gate mapping
        gateMap[oldGate] = patchedGate;
        
        // add the fanins
        for(unsigned j=0; j<oldGate->getNumFanins(); ++j) {
            auto oldFanin = oldGate->getFanin(j);
            // auto mappedGate = gateMap.at(oldFanin);
            // auto mappedGate = patchedNtk->getGateByName(oldFanin->getGateName());
            patchedGate->addFaninName(oldFanin->getGateName());
        }
        
        patchedNtk->addGate(patchedGate);

        // If a gate is PI, then also add it to the PI of patched circuit
        // if(oldGate->isPi())
        //     patchedNtk->addPi(patchedGate);
    }

    

    // 4. do euivalence checking
    patchedNtk->genConnection();
    patchedNtk->writeNtkVerilog("patched.v");

    // patchedNtk->readNtkFile("./patched.v");
    gv::cir::EcoNtk::rewriteDesign(getNewDesignName());
    Abc_Ntk_t* pNtkNew = Io_Read("./tmp.v", IO_FILE_VERILOG, 0, 0 );
    Abc_Ntk_t* pNtkPatched = Io_Read("./patched.v", IO_FILE_VERILOG, 0, 0 );
    assert(Abc_NtkCheck(pNtkPatched));
    assert(Abc_NtkCheck(pNtkNew));
    pNtkNew = Abc_NtkStrash(pNtkNew, 0, 1, 0);
    pNtkPatched = Abc_NtkStrash(pNtkPatched, 0, 1, 0);

    assert(Abc_NtkCoNum(pNtkPatched) == Abc_NtkCoNum(pNtkNew));
    bool ret = true;
    for(int i=0; i < Abc_NtkCoNum(pNtkPatched); ++i) {
        
        Abc_Obj_t* pachedPo = Abc_NtkCo(pNtkPatched, i);
        Abc_Obj_t* newPo = Abc_NtkCo(pNtkNew, i);
        // cout << "patched po " << Abc_ObjName(pachedPo) << " new po " << Abc_ObjName(newPo) << endl;
        Abc_Ntk_t * pNtkPatchedCone = Abc_NtkCreateCone( pNtkPatched, Abc_ObjFanin0(pachedPo), Abc_ObjName(pachedPo), 0 );
        Abc_Ntk_t * pNtkNewCone = Abc_NtkCreateCone( pNtkNew, Abc_ObjFanin0(newPo), Abc_ObjName(newPo), 0 );
        if ( Abc_ObjFaninC0(pachedPo) ) Abc_ObjXorFaninC( Abc_NtkPo(pNtkPatchedCone, 0), 0 );
        if ( Abc_ObjFaninC0(newPo) ) Abc_ObjXorFaninC( Abc_NtkPo(pNtkNewCone, 0), 0 );
        if(!isNtkEq(pNtkPatchedCone, pNtkNewCone)) {
            cout << "patched po " << Abc_ObjName(pachedPo) << " neq to new po " << Abc_ObjName(newPo) << endl;
            ret = false;
        }
    }
    // bool ret = isNtkEq(pNtkNew, pNtkPatched);
    // extern bool isNtkEq(Abc_Ntk_t* pNtkOld, Abc_Ntk_t* pNtkNew);

    delete patchedNtk;
    return ret;
}

}
}

#endif