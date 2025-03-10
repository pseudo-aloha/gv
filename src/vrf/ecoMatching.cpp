#ifndef ECO_MATCHING_CPP
#define ECO_MATCHING_CPP
#include "ecoNtk.h"
#include "ecoMgr.h"

namespace gv {
namespace eco {

// do the cut matching from output side
void
EcoMgr::doOutputSideMatching() {
    unsigned nPo = _oldNtk->getNumPos();
    
    for(unsigned i=0; i<nPo; i++) {
        if(i > 0) break;
        matchOnePo(i);
    }
}

// do the matching for one po
void
EcoMgr::matchOnePo(unsigned ithPo) {
    gv::cir::EcoGate* oldPo = _oldNtk->getPo(ithPo);
    gv::cir::EcoGate* newPo = _newNtk->getPo(ithPo);

    auto oldPoCuts = _oldNtk->getGateCuts(oldPo);
    auto newPoCuts = _newNtk->getGateCuts(newPo);
    
    cout << "old cuts : " << oldPoCuts.size() << endl;
    cout << "new cuts : " << newPoCuts.size() << endl;
    
    // sort the cuts by their NPN class
    for(auto cut : oldPoCuts) {
        _oldNtk->computeCutTT(cut);
        // break;
    }

    for(auto cut : newPoCuts) {
        _newNtk->computeCutTT(cut);
        break;
    }
}

// end of namespace gv::eco
}
// end of namespace gv
}
#endif