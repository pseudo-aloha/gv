#ifndef ECO_SELECTOR_CPP
#define ECO_SELECTOR_CPP

#include "ecoNtk.h"
#include "ecoMgr.h"

// this file includes the functions that builds the selector circuit

namespace gv {
namespace eco {
// build the selector circuit
void
EcoMgr::buildSelector() {
    gv::cir::EcoGate::setGlobalTrav(); // increment the trav flag
    for(unsigned i=0; i<_oldNtk->getNumPos(); ++i) {
        buildSelectorForIthPo(i);
    }

    _selectorNtk->writeNtkVerilog("selector.v");
}

// build the selector circuit for ith po
void
EcoMgr::buildSelectorForIthPo(unsigned i) {
    gv::cir::EcoGate* selectorPo = new gv::cir::EcoGate("po", "selectorPo" + to_string(i));
    _selectorNtk->addPo(selectorPo);
    
    gv::cir::EcoGate* selectorPoXor = new gv::cir::EcoGate("xor", "outputXor" + to_string(i));
    _selectorNtk->addGate(selectorPoXor);

    // add po i logic of new circuit and add it to the xor gate
    auto newPo = _newNtk->getPo(i)->getFanin(0); // ignore the dummy po
    buildSelectorForIthNewPoRec(selectorPoXor, newPo, i);

    // add mux selected logic for old circuit
    auto oldPo = _oldNtk->getPo(i)->getFanin(0); // ignore the dummy po
    buildSelectorForIthOldPoRec(selectorPoXor, oldPo, i);

}

// recursively build the new circuit part of selecotor circuit
void
EcoMgr::buildSelectorForIthNewPoRec(gv::cir::EcoGate* selectorGate, gv::cir::EcoGate* origNtkGate, unsigned ithPo) {
    // if the gate is already traversed, return
    if(origNtkGate->isGlobalTrav())
        return;
    origNtkGate->setToGlobalTrav();

    gv::cir::EcoGate* g = new gv::cir::EcoGate(origNtkGate->getGateTypeName(), origNtkGate->getGateFullName());
    _selectorGateMap[g] = origNtkGate; // record the gate in the selector maps to which gate
    _selectorGateMap[origNtkGate] = g;
    selectorGate->addFanin(g);
    _selectorNtk->addGate(g);
    selectorGate = g;

    // if the gate is const gate or PI gate, return
    if(origNtkGate->isPiOrConst()) return;

    // recursive traverse
    for(unsigned i=0; i<origNtkGate->getNumFanins(); ++i)  
        buildSelectorForIthNewPoRec(g, origNtkGate->getFanin(i), ithPo);
    
    return;
}

// recursively build the old circuit part of selecotor circuit
void
EcoMgr::buildSelectorForIthOldPoRec(gv::cir::EcoGate* selectorGate, gv::cir::EcoGate* origNtkGate, unsigned ithPo) {
    // create the corresponding gates in the selector circuit
    auto ithPoRp = _rpTable.at(ithPo);
    
    // if the gate is fixed to a gate in the new circuit
    if(ithPoRp.count(origNtkGate)) {
        auto pEcoRpInfo = ithPoRp.at(origNtkGate);
        auto mappedGate = pEcoRpInfo->getMappedGate();
        
        // if the gate is mapped in inv pole
        if(pEcoRpInfo->getMappedPole()) {
            // if not the const gate case
            if(!mappedGate->isConstGate()) {
                cout << "inv " << endl;
                gv::cir::EcoGate* g = new gv::cir::EcoGate("not", mappedGate->getGateFullName() + "_INV");
                g->addFanin(mappedGate);
                _selectorNtk->addGate(g);
                selectorGate->addFanin(g);
                selectorGate = g;
            }
            // do some special trick for const gate
            else {

            }
        }
        cout << "chaing orignatkgate : " << origNtkGate->getGateFullName() << " " << pEcoRpInfo->getMappedGate()->getGateFullName() << endl;
        origNtkGate = pEcoRpInfo->getMappedGate();
        
    }

    // if the gate is already traversed, return
    if(origNtkGate->isGlobalTrav())
        return;
    // mark the gate as traversed
    origNtkGate->setToGlobalTrav();
    cout << "trav : " << origNtkGate->getGateFullName() << endl;
    gv::cir::EcoGate* g = new gv::cir::EcoGate(origNtkGate->getGateTypeName(), origNtkGate->getGateFullName());
    _selectorGateMap[g] = origNtkGate; // record the gate in the selector maps to which gate
    _selectorGateMap[origNtkGate] = g;
    selectorGate->addFanin(g);
    cout << "map " << origNtkGate->getGateFullName() << " " << g->getGateFullName() << endl;
    _selectorNtk->addGate(g);
    selectorGate = g;

    // if the gate is const gate or PI gate, return
    if(origNtkGate->isPiOrConst()) return;

    // recursive traverse
    for(unsigned i=0; i<origNtkGate->getNumFanins(); ++i) {    
        buildSelectorForIthOldPoRec(g, origNtkGate->getFanin(i), ithPo);
        // selectorGate->addFanin(selectorGateFanin); // add the newly created gate as fanin
    }
    return;
}

}}

#endif