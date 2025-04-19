#ifndef ECO_PATCH_CPP
#define ECO_PATCH_CPP

#include "ecoMgr.h"

namespace gv {
namespace eco {

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
EcoMgr::genPatch() {
    gv::cir::EcoGate::setGlobalTrav();
    // handle the rewire for output cut matching
    // decide which net should use rewire
    decideOutputRewire();

    // collect patch circuit from each PO of the old circuit
    // if a gate that is reached is fixed to another gate in the new circuit, replace it with the corresponding gate
    for(unsigned i=0; i<_newNtk->getNumPos(); ++i)
        collectPatchGates(_newNtk->getPo(i)->getFanin(0), nullptr, i, true);

    // write the patch ntk
    _patchNtk->writeNtkVerilog("patch.v"); // write the patch content

}

void
EcoMgr::collectPatchGates(gv::cir::EcoGate* g, gv::cir::EcoGate* curPatchGate, const unsigned& ithPo, bool isEntry) {
    if(g->isGlobalTrav()) return;
    g->setToGlobalTrav();

    gv::cir::EcoGate* patchGate;
    
    // check if the gate is in old circuit to decide how we handle it
    if(g->isInOldCircuit()) {
        // 1. check if the gate is fixed to another gate
        auto rpPairsForIthPo = _rpTable.at(ithPo);
        if(rpPairsForIthPo.count(g)) {
            auto pEcoRpInfo = rpPairsForIthPo.at(g);
            gv::cir::EcoGate* mappedGate = pEcoRpInfo->getMappedGate();
            auto inv = pEcoRpInfo->getMappedPole();
            collectPatchGates(mappedGate, nullptr, ithPo, true);
            return;
        }
        // if not found, further traverse in the old circuit
        // else {
        //     for(unsigned i=0; i<g->getNumFanins(); ++i) {
        //         auto faninGate = g->getFanin(i);
        //         collectPatchGates(faninGate, nullptr, ithPo, false);
        //     }
        // }
    }
    else {
        assert(g->isInNewCircuit()); // check that the circuit is in the new circuit XDDD, just in case
        cout << "new circuit : " << g->getGateFullName() << endl;
        if(!_patchNtk->getGateByName(g->getGateName())) {
            patchGate = new gv::cir::EcoGate(g->getGateTypeName(), g->getGateName());
            _patchNtk->addGate(patchGate);
            if(isEntry) {
                gv::cir::EcoGate* patchPoGate = new gv::cir::EcoGate("po", g->getGateName());
                patchPoGate->addFanin(patchGate);
                _patchNtk->addPo(patchPoGate);
            }
        }
        else {
            patchGate = _patchNtk->getGateByName(g->getGateName());
        }

        // add the patch gate as a fanin of the fanout of patch gate
        if(curPatchGate)
            curPatchGate->addFanin(patchGate);

        // if the gate itself is merged, no further traversal is needed
        if(isMerged(g)) return;

        // not merged, further traverse in the new circuit
        // g = patchGate;
        // for(unsigned i=0; i<g->getNumFanins(); ++i) {
        //     auto faninGate = g->getFanin(i);
        //     collectPatchGates(faninGate, patchGate, ithPo, true);
        // }
    }

    cout << "collecting " << g->getGateFullName() << endl;

    // if the gate is a PI gate or const gate, stop traversing
    if(g->isPiOrConst()) return;
    g->reportGate();
    for(unsigned i=0; i<g->getNumFanins(); ++i) {
        auto faninGate = g->getFanin(i);
        cout << "fanin : " << faninGate->getGateFullName() << endl;
        collectPatchGates(faninGate, patchGate, ithPo, true);
    }

    // // get fanin
    // for(unsigned i=0; i<g->getNumFanins(); ++i) {
    //     auto fanin = g->getFanin(i);
    //     collectPatchGates(fanin, ithPo, false);
    //     auto faninPatchGate = new gv::cir::EcoGate(fanin->getGateTypeName(), fanin->getGateName());
    //     patchGate->addFanin(faninPatchGate);
    //     _patchNtk->addGate(faninPatchGate);
    // }
}

}
}

#endif