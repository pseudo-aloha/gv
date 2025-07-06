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


// find one merged gate of the given gate
// tend to find the gate of the specified pole
// if same pole is not found, return the opposite polfe gate and the bool will be true
// user should make sure that the gate is a merged gate
pair<gv::cir::EcoGate*, bool>
EcoMgr::getOneMergedGate(gv::cir::EcoGate* g, bool pole) {
    assert(isMerged(g)); // make sure that the gate is merged

    gv::cir::EcoGate* mergedGate = nullptr;
    auto mergedGates = getMergedGates(g);
    auto invMergedGates = getInvMergedGates(g);

    if(!pole) {
        if(!mergedGates.empty())
            return make_pair(*mergedGates.begin(), false);
        else
            return make_pair(*invMergedGates.begin(), true);
    }
    else {
        if(!invMergedGates.empty())
            return make_pair(*invMergedGates.begin(), false);
        else
            return make_pair(*mergedGates.begin(), true);
    }
}

// decide the output side rewire
void
EcoMgr::decideOutputRewire() {
    unordered_map<gv::cir::EcoGate*, vector<pair<gv::cir::EcoGate*, array<unsigned, 2>>>> mp;
    
    unordered_set<gv::cir::EcoGate*> frozenGates;

    // iterate over the rp pair to see the old net is fixed to which net most frequently
    for(unsigned i=0; i<_rpTable.size(); ++i) {
        unordered_set<gv::cir::EcoGate*> visited; // record the visited gates that are used for checking gate fix condition
        auto rpForIthPo = _rpTable.at(i);
        for(auto[oldGate, pEcoRpInfo] : rpForIthPo) {
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

            // find the gates in the fanin cone that is merged back from the new circuit
            // use BFS to traverse the circuit
            queue<gv::cir::EcoGate*> q;
            q.push(mappedGate);
            while(!q.empty()) {
                auto cur = q.front();
                q.pop();
                
                // mark visited
                if(visited.count(cur)) continue;
                visited.insert(cur);
                
                //
                if(isMerged(cur) && cur->isInNewCircuit()) {
                    auto[mergedGate, mergedPole] = getOneMergedGate(cur, false);
                    // cout << "mp " << mergedGate->getGateFullName() << endl;
                    if(!visited.count(mergedGate)) 
                        q.push(mergedGate);
                }
                else {
                    if(cur->isInOldCircuit() && cur != mappedGate) {
                        // cout << "freeze1 : " << cur->getGateFullName() << endl;
                        frozenGates.insert(cur);
                    }
                    for(unsigned j=0; j<cur->getNumFanins(); ++j) {
                        auto fanin = cur->getFanin(j);
                        if(!visited.count(fanin))
                            q.push(fanin);
                    }
                }
            }
        }

        // mark the output side cone that should not be changed
        queue<gv::cir::EcoGate*> q;
        q.push(_oldNtk->getPo(i)->getFanin(0));
        while(!q.empty()) {
            auto cur = q.front();
            q.pop();
            
            // mark visited
            if(visited.count(cur)) continue;
            visited.insert(cur);
            
            

            if(rpForIthPo.count(cur)) continue;
            
            if(cur != _oldNtk->getPo(i)->getFanin(0)) {
                // cout << "freeze2 : " << cur->getGateFullName() << endl;
                frozenGates.insert(cur);
            }

            for(unsigned j=0; j<cur->getNumFanins(); ++j) {
                auto fanin = cur->getFanin(j);
                if(!visited.count(fanin))
                    q.push(fanin);
            }
        }
    }

    


    // find the rp pairs that is used most frequently



    // write the final mapping info to a file
    ofstream f;
    f.open(".eco_map_info.txt", ios::out);
    f << "mapped counts : " << endl;
    for(const auto&[oldGate, vec] : mp) {
        unsigned maxCount = 0;
        pair<gv::cir::EcoGate*, bool> cand;
        f << oldGate->getGateFullName() << " << ";
        for(const auto&[mappedGate, poleCounts] : vec) {
            if(poleCounts[0]) {
                if(poleCounts[0] > maxCount) {
                    maxCount = poleCounts[0];
                    cand = make_pair(mappedGate, false);
                }
                f << mappedGate->getGateFullName() << " : " << poleCounts[0] << " | ";
            }
            if(poleCounts[1]) {
                f << "~" << mappedGate->getGateFullName() << " : " << poleCounts[1] << " | ";
                if(poleCounts[1] > maxCount) {
                    maxCount = poleCounts[1];
                    cand = make_pair(mappedGate, true);
                }
            }
        }
        _finalRpPair[oldGate] = cand;
        f << endl;
        if(frozenGates.count(oldGate))
            f << "......................... " << oldGate->getGateFullName() << endl;
    }
    
    f << endl << "final selection : " << endl;
    for(const auto&[oldGate, rpInfo] : _finalRpPair) {
        auto[mappedGate, mappedPole] = rpInfo;
        f << oldGate->getGateFullName() << " " << (mappedPole ? "~" : "") << mappedGate->getGateFullName() << endl;
    }
    f.close();

}

// generate the patch
void
EcoMgr::genPatch(const string& patchName) {
    gv::cir::EcoGate::setGlobalTrav();
    // handle the rewire for output cut matching
    // decide which net should use rewire
    decideOutputRewire();
    
    // record the used new gate, will generate the logic later
    unordered_set<gv::cir::EcoGate*> usedNewGate;
    // rewire the final rp
    generateFinalRpRewire(usedNewGate);

    // For each po, write the rewire info (to fix the net in the old circuit to the net in the new circuit)
    // this part takes care the gate logic duplication
    for(unsigned i=0; i<_oldNtk->getNumPos(); ++i) {
        cout << "output " << _oldNtk->getPo(i)->getFanin(0)->getGateFullName() << endl;
        // traverse the gates in the old circuit
        generatePatchForIthPo(i);
    }
    generateNewGateLogics(usedNewGate);

    // generate the connections for the patch gate
    _patchNtk->genConnection();

    // write the patch ntk
    _patchNtk->writeNtkVerilog("patch.v"); // write the patch content
    _patchNtk->computeCadContestCost();

    // apply the patch to old circuit and do equivalence checking between it and the new circuit
    if(!applyNCheckPatch(patchName))
        cout << "patched circuit NEQ to new circuit!!!!" << endl;
    else
        cout << "patched circuit eq to new circuit!!!!" << endl;
}

string
EcoMgr::getPatchGateName(gv::cir::EcoGate* g) {
    assert(g->isInNewCircuit());
    if(isMerged(g)) {
        auto[mergedGate, pole] = getOneMergedGate(g, false);
        string gateName = mergedGate->getGateName();
        if(isInPatchPoNames(gateName)) gateName += "_in";
        if(pole) gateName += "_inv"; // add _inv suffix

        return gateName;
    }
    else {
        return g->getGateFullName(); // return the full name of the new circuit
    }
}

// write the rewiring for the final selection of rp pairs
void
EcoMgr::generateFinalRpRewire(unordered_set<gv::cir::EcoGate*>& usedNewGate) {
    for(auto[oldGate, mappedGateInfo] : _finalRpPair) {
        auto mappedGate = mappedGateInfo.first;
        auto inv = mappedGateInfo.second;

        // create rewiring gate
        string gateTypeName = (inv ? "not" : "buf");
        string gateName = oldGate->getGateName();
        // create patch po for the rewiring
        gv::cir::EcoGate* po = createOrGetPatchGate("po", gateName);
        // create the gate
        gv::cir::EcoGate* g = createOrGetPatchGate(gateTypeName, gateName);
        
        // add the rewire gate name to the fanin, will create the gate later
        g->addFaninName(mappedGate->getGateFullName());

        // gv::cir::EcoGate* fanin = createOrGetPatchGate(mappedGate->getGateTypeName(), mappedGate->getGateFullName());
        // g->addFanin(fanin);
        
        usedNewGate.insert(mappedGate);
    }
}

// generate the new gates' logic in the patch circuit
void
EcoMgr::generateNewGateLogics(unordered_set<gv::cir::EcoGate*>& usedNewGate) {
    // Use BFS to generate the logics
    queue<gv::cir::EcoGate*> q;
    gv::cir::EcoGate::setGlobalTrav();
    
    for(const auto& g : usedNewGate)
        q.push(g);
    
    while(!q.empty()) {
        auto cur = q.front();
        q.pop();
        
        if(cur->isGlobalTrav()) continue;
        cur->setToGlobalTrav();
        
        // create the corresponding gate in the patch circuit
        gv::cir::EcoGate* patchGate = nullptr;

        // stop traversal when const or pi is reached
        if(cur->isPiOrConst()) {
            patchGate = createOrGetPatchGate(cur->getGateTypeName(), cur->getGateFullName());
            continue;
        }

        // if the gate is merged to a gate in the old circuit, replace it with the merged gate and stop traversal in the new circuit
        if(isMerged(cur)) {
            auto[mergedGate, pole] = getOneMergedGate(cur, false);
            string mergedPatchGateName = getPatchGateName(cur);
            string mergedPatchGateTypeStr = (pole ? "not" : "pi");
            gv::cir::EcoGate* mergedPatchGate = createOrGetPatchGate(mergedPatchGateTypeStr, mergedPatchGateName);
            patchGate = createOrGetPatchGate("buf", cur->getGateFullName());
            patchGate->addFaninName(mergedPatchGateName);

            // if the gate is in inverted pole, we need add the original gate to the inverter
            if(pole) {
                string faninGateName = mergedGate->getGateName();
                if(isInPatchPoNames(faninGateName)) faninGateName += "_in";
                gv::cir::EcoGate* patchFaninGate = createOrGetPatchGate("pi", faninGateName);
                mergedPatchGate->addFaninName(faninGateName);
            }
            continue;
        }
        
        patchGate = createOrGetPatchGate(cur->getGateTypeName(), cur->getGateFullName());

        // traverse its children
        for(unsigned i=0; i<cur->getNumFanins(); ++i) {
            auto fanin = cur->getFanin(i);
            patchGate->addFaninName(getPatchGateName(fanin));
            if(!fanin->isGlobalTrav())
                q.push(fanin);
        }
    }
}

// generate the patch logic for ith po, except the ones that are defined in final rp
// this is used to generate the duplicated circuit part
void
EcoMgr::generatePatchForIthPo(unsigned i) {
    auto rpForIthPo = _rpTable.at(i);
    const string poSuffix = "_po" + to_string(i);
    // add the rewire things
    for(auto[oldGate, pEcoRpInfo] : rpForIthPo) {
        auto mappedGate = pEcoRpInfo->getMappedGate();
        auto inv = pEcoRpInfo->getMappedPole();

        // if the rewire is handled by our final decision, no need to rewire again
        if(inv == _finalRpPair.at(oldGate).second && mappedGate == _finalRpPair.at(oldGate).first) continue;
        cout << "old gate " << oldGate->getGateFullName() << endl;
        cout << "is fixed to " << (inv ? "~" : "") << mappedGate->getGateFullName() << endl;
        
        // create rewiring gate
        string gateTypeName = (inv ? "not" : "buf");
        string gateName = oldGate->getGateName();
        gv::cir::EcoGate* g = createOrGetPatchGate(gateTypeName, gateName);

        gv::cir::EcoGate* fanin = createOrGetPatchGate(mappedGate->getGateTypeName(), mappedGate->getGateFullName() + poSuffix);
        g->addFanin(fanin);
    }

    // duplicate the output side gates
    queue<gv::cir::EcoGate*> q;
    unordered_set<gv::cir::EcoGate*> visited;
    q.push(_oldNtk->getPo(i)->getFanin(0));
    while(!q.empty()) {
        auto cur = q.front();
        q.pop();
        
        // mark visited
        if(visited.count(cur)) continue;
        visited.insert(cur);
        
        // this means that the output frontier is reached, no need to further duplicate the gates
        if(rpForIthPo.count(cur)) continue;

        // add its fanins
        bool dupGate = false;
        for(unsigned j=0; j<cur->getNumFanins(); ++j) {
            auto fanin = cur->getFanin(j);
            // if the fanin is changed by other rp point, dup the gate
            if(_finalRpPair.count(fanin) && !rpForIthPo.count(fanin))
                dupGate = true;
            if(!visited.count(fanin))
                q.push(fanin);
        }
        if(dupGate) {
            
        }
    }

}

// check if the gate's name is already in the patch gate's fanin
bool checkGateNameExists(gv::cir::EcoGate* g, const string& name) {
    for(unsigned i=0; i<g->getNumFanins(); ++i) {
        auto fanin = g->getFanin(i);
        if(fanin->getGateName() == name) return true;
    }
    return false;
}

// check if the patch gate exists, if not create it
gv::cir::EcoGate*
EcoMgr::createOrGetPatchGate(const string& gateTypeName, const string& gateName) {
    if(gateTypeName == "po") {
        auto g = _patchNtk->getPoByName(gateName);
        if(!g) {
            g = new gv::cir::EcoGate(gateTypeName, gateName);
        }
        _patchNtk->addPo(g);
        addPatchPoName(gateName); // add the name to the po name set, it is uswd to check if the pi is also po
        return g;
    }
    else {
        auto g = _patchNtk->getGateByName(gateName);
        // there are 3 possibilities
        // 1. The gate name does not exist -> create the gate along with the gate type
        // 2. the gate name exists and the gate type matches -> No need to create a new gate and return the existing gate
        // 3. the gate name exists, but the gate type does not matches -> we are using the old gate in both PI and PO, and we need to rename it with postfix "_in"
        if(g) {
            if(g->getGateTypeName() != gateTypeName) {
                if(g->getGateType() == gv::cir::EcoGate::ECO_PI_GATE) {
                    auto name = g->getGateName();
                    if(name.substr(name.size()-3, 3) == "_in") {
                        cout << "yabe 1 " << name << " " << gateName << " " << g->getGateTypeName() << " " << gateTypeName << endl;
                    }
                    _patchNtk->setGateName2Gate(g->getGateName() + "_in", g->getGateName(), g);
                    g = new gv::cir::EcoGate(gateTypeName, gateName);
                    _patchNtk->addGate(g);
                }
                else {
                    if(gateName.substr(gateName.size()-3, 3) == "_in") {
                        cout << "yabe 2 " << gateName << endl;
                    }
                    g = _patchNtk->getGateByName(gateName+"_in");
                    if(!g) {
                        g = new gv::cir::EcoGate(gateTypeName, gateName + "_in");
                        _patchNtk->addGate(g);
                    }
                }   
            }
        }
        // case 1.
        else {
            g = new gv::cir::EcoGate(gateTypeName, gateName);
            _patchNtk->addGate(g);
        }
        return g;
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
  Abc_Ntk_t* pNtkMiter = Abc_NtkMiter( pNtkOld, pNtkNew, 1, 0, 0, 0 );
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
        
        auto  patchedGate = patchedNtk->getGateByName(patchedGateName);
        if(!patchedGate) {
            patchedGate = new gv::cir::EcoGate(oldGate->getGateTypeName(), patchedGateName);
            patchedNtk->addGate(patchedGate);
        }
        // do the gate mapping
        gateMap[oldGate] = patchedGate;
        
        // add the fanins
        for(unsigned j=0; j<oldGate->getNumFanins(); ++j) {
            auto oldFanin = oldGate->getFanin(j);
            patchedGate->addFaninName(oldFanin->getGateName());
        }
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
    int neqCount = 0, eqCount = 0;
    for(int i=0; i < Abc_NtkCoNum(pNtkPatched); ++i) {
        
        Abc_Obj_t* pachedPo = Abc_NtkCo(pNtkPatched, i);
        Abc_Obj_t* newPo = Abc_NtkCo(pNtkNew, i);
        Abc_Ntk_t * pNtkPatchedCone = Abc_NtkCreateCone( pNtkPatched, Abc_ObjFanin0(pachedPo), Abc_ObjName(pachedPo), 1 );
        Abc_Ntk_t * pNtkNewCone = Abc_NtkCreateCone( pNtkNew, Abc_ObjFanin0(newPo), Abc_ObjName(newPo), 1 );
        if ( Abc_ObjFaninC0(pachedPo) ) Abc_ObjXorFaninC( Abc_NtkPo(pNtkPatchedCone, 0), 0 );
        if ( Abc_ObjFaninC0(newPo) ) Abc_ObjXorFaninC( Abc_NtkPo(pNtkNewCone, 0), 0 );
        if(!isNtkEq(pNtkPatchedCone, pNtkNewCone)) {
            cout << "patched po " << Abc_ObjName(pachedPo) << " neq to new po " << Abc_ObjName(newPo) << endl;
            ret = false;
            ++neqCount;
        }
        else {
            ++eqCount;
            cout << "eq po : " << Abc_ObjName(pachedPo) << endl;
        }
    }
    cout << "eq report :" << endl;
    cout << "# eq po's : " << eqCount << endl;
    cout << "# neq po's : " << neqCount << endl;
    // bool ret = isNtkEq(pNtkNew, pNtkPatched);
    // extern bool isNtkEq(Abc_Ntk_t* pNtkOld, Abc_Ntk_t* pNtkNew);

    delete patchedNtk;
    return ret;
}

}
}

#endif