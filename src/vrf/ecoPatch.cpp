#ifndef ECO_PATCH_CPP
#define ECO_PATCH_CPP

#include "ecoMgr.h"

namespace gv {
namespace eco {

// generate the patch
void
EcoMgr::genPatch() {
    gv::cir::EcoGate::setGlobalTrav();
    
    // collect patch circuit from each PO of the new circuit
    for(unsigned i=0; i<_newNtk->getNumPos(); ++i) {
        collectPatchGates(_newNtk->getPo(i)->getFanin(0), true);
    }

    // write the patch ntk
    ofstream f("patch.v");
    f << "module top(";
    for(unsigned i=0; i<_patchNtk->getNumPos(); ++i) {
        f << _patchNtk->getPo(i)->getGateName();
        f << ", ";
    }
    for(unsigned i=0; i<_patchNtk->getNumPis(); ++i) {
        f << _patchNtk->getPi(i)->getGateName();
        if(i < _patchNtk->getNumPis() - 1)
            f << ", ";
    }
    f << ");" << endl;
    f << "output ";
    for(unsigned i=0; i<_patchNtk->getNumPos(); ++i) {
        if((i + 1) % 10 == 0) {
            f << ";" << endl;
            f << "output ";
        }
        f << _patchNtk->getPo(i)->getGateName();
        if(i < _patchNtk->getNumPos() - 1 && (i + 2) % 10 != 0)
            f << ", ";
    }
    f << ";" << endl;
    
    f << "input ";
    for(unsigned i=0; i<_patchNtk->getNumPis(); ++i) {
        if((i + 1) % 10 == 0) {
            f << ";" << endl;
            f << "input ";
        }
        f << _patchNtk->getPi(i)->getGateName();
        if(i < _patchNtk->getNumPis() - 1 && (i + 2) % 10 != 0)
            f << ", ";
    }
    f << ";" << endl;

    f << "wire ";
    for(unsigned i=0; i<_patchNtk->getNumGates(); ++i) {
        if((i + 1) % 10 == 0) {
            f << ";" << endl;
            f << "wire ";
        }
        f << _patchNtk->getGate(i)->getGateName();
        if(i < _patchNtk->getNumGates() - 1 && (i + 2) % 10 != 0)
            f << ", ";
    }
    f << ";" << endl;
    for(unsigned i=0; i<_patchNtk->getNumGates(); ++i) {
        auto g = _patchNtk->getGate(i);
        f << g->getGateTypeName() << " (" << g->getGateName() << ", ";
        for(unsigned j=0; j<g->getNumFanins(); ++j) {
            auto fanin = g->getFanin(j);
            f << fanin->getGateName();
            if(j < g->getNumFanins() - 1)
                f << ", ";
        }
        f << ");" << endl;
    }
    f << "endmodule" << endl;
    f.close();
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