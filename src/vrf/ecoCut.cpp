#ifndef ECO_CUT_CPP
#define ECO_CUT_CPP

#include "ecoNtk.h"
#include "base/abc/abc.h"
#include <iostream>
#include <string>
#include <cassert>

namespace gv {
namespace cir {


void
EcoCut::reportCut() {
    cout << "root : " << _root->getGateName() << endl;
    cout << "leaves : ";
    for(const auto& leaf : _leaves) {
        cout << leaf->getGateName() << " ";
    }
    cout << endl;
}

// enumerate k feasible cuts
void
EcoNtk::enumerateCuts(unsigned k) {
    EcoGate::setGlobalTrav();
    // traverse from POs
    for(size_t i=0; i<getNumPos(); i++) {
        enumerateCutsRec(k, getPo(i));
    }
}

void
getCutCombs(unsigned faninIdx, EcoGate* root, vector<EcoGate*>& leaves, vector<EcoCut*>& cuts, vector<vector<EcoCut*>>& faninCutVec, const unsigned& k) {
    unsigned nFanins = faninCutVec.size();
    // cout << "faninidx" << faninIdx << endl;
    if(faninIdx == nFanins) {
        if(leaves.size() <= k) {
            EcoCut* cut = new EcoCut(root, leaves);
            cuts.push_back(cut);
            // cout << "root ptr : " << cut->getRoot() << " " << cut->getRoot()->getGateTypeName() << endl;
            // cut->reportCut();
        }
        return;
    }

    // get the cut from idx'th fanin
    assert(!faninCutVec.empty());
    auto faninCuts = faninCutVec.at(faninIdx);
    for(size_t i=0; i<faninCuts.size(); i++) {
        auto faninCut = faninCuts.at(i);
        vector<EcoGate*> faninCutLeaves = faninCut->getLeaves();
        leaves.insert(leaves.end(), faninCutLeaves.begin(), faninCutLeaves.end());
        // recursive call
        getCutCombs(faninIdx+1, root, leaves, cuts, faninCutVec, k);
        for(size_t j=0; j<faninCut->getLeaves().size(); j++)
            leaves.pop_back();
    }

    
}

// enumerate k feasible cuts
vector<EcoCut*>
EcoNtk::enumerateCutsRec(const unsigned& k, EcoGate* g) {
    if(g->isGlobalTrav()) return _gate2Cuts.at(g);
    g->setToGlobalTrav();
    
    // boundary case, return the PI gate when PI is reached
    if(g->getGateType() == EcoGate::ECO_PI_GATE || g->getGateType() == EcoGate::ECO_CONST_0_GATE || g->getGateType() == EcoGate::ECO_CONST_1_GATE)  {
        EcoCut* cut = new EcoCut(g, {g});
        _gate2Cuts[g].push_back(cut);
        return _gate2Cuts[g];
    }

    // otherwise recursively collect k-feasible cuts
    vector<vector<EcoCut*>> faninCutVec;
    unsigned nFanins = g->getNumFanins();
    for(size_t i=0; i<nFanins; i++) {
        auto faninCut = enumerateCutsRec(k, g->getFanin(i));
        faninCutVec.push_back(faninCut);
    }
    vector<EcoCut*> ret;
    vector<size_t> faninCutIdxs(nFanins);
    vector<EcoGate*> leaves;
    getCutCombs(0, g, leaves, ret, faninCutVec, k);

    _gate2Cuts[g].insert(_gate2Cuts[g].end(), ret.begin(), ret.end());
    return ret;
    
}

}
}
#endif