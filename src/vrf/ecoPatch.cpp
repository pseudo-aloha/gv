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

    // collect patch circuit from each PO of the new circuit
    // for(unsigned i=0; i<_newNtk->getNumPos(); ++i) {
    //     collectPatchGates(_newNtk->getPo(i)->getFanin(0), true);
    // }

    // write the patch ntk
    _patchNtk->writeNtkVerilog("patch.v"); // write the patch content

}

void
EcoMgr::collectPatchGates(gv::cir::EcoGate* g, bool isEntry) {
    if(g->isGlobalTrav()) return;
    g->setToGlobalTrav();
    if(g->getGateType() == gv::cir::EcoGate::ECO_PI_GATE || g->getGateType() == gv::cir::EcoGate::ECO_CONST_0_GATE || g->getGateType() == gv::cir::EcoGate::ECO_CONST_1_GATE) return;
    if(isMerged(g)) return;
    gv::cir::EcoGate* patchGate;
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
    
    
    cout << "collecting " << g->getGateFullName() << endl;

    // get fanin
    for(unsigned i=0; i<g->getNumFanins(); ++i) {
        auto fanin = g->getFanin(i);
        collectPatchGates(fanin, false);
        auto faninPatchGate = new gv::cir::EcoGate(fanin->getGateTypeName(), fanin->getGateName());
        patchGate->addFanin(faninPatchGate);
        _patchNtk->addGate(faninPatchGate);
    }
}

}
}

#endif