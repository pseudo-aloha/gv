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
}
// build the selector circuit for ith po
void
EcoMgr::buildSelectorForIthPo(unsigned i) {
    gv::cir::EcoGate* g = new gv::cir::EcoGate("xor", "outputXor" + to_string(i));
    _selectorNtk->addPo(g);
    auto oldPo = _oldNtk->getPo(i);
    buildSelectorForIthPoRec(g, oldPo);
}

// recursively build the selecotor circuit
void 
EcoMgr::buildSelectorForIthPoRec(gv::cir::EcoGate* selectorGate, gv::cir::EcoGate* origNtkGate) {
    // if the gate is already traversed, return
    if(origNtkGate->isGlobalTrav()) return;
    origNtkGate->setToGlobalTrav();

    // if the gate is const gate or PI gate, return
    if(origNtkGate->getGateType() == gv::cir::EcoGate::ECO_PI_GATE || origNtkGate->getGateType() == gv::cir::EcoGate::ECO_CONST_0_GATE || origNtkGate->getGateType() == gv::cir::EcoGate::ECO_CONST_1_GATE) return;


    // create the corresponding gates in the selector circuit
    gv::cir::EcoGate* g = new gv::cir::EcoGate(origNtkGate->getGateTypeName(), origNtkGate->getGateFullName());
    _selectorNtk->addGate(g);
    selectorGate->addFanin(g); // add the newly created gate as fanin
    
    // recursive traverse
    for(unsigned i=0; i<origNtkGate->getNumFanins(); ++i) {
        buildSelectorForIthPoRec(g, origNtkGate->getFanin(i));
    }
}

}}

#endif