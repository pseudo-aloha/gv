#ifndef ECO_CUT_CPP
#define ECO_CUT_CPP

#include "ecoNtk.h"
#include "base/abc/abc.h"
#include <iostream>
#include <string>
#include <cassert>
#include <iomanip>

namespace gv {
namespace cir {

// return the i'th (0-indexed) bit of size_t
bool getIthBit(const size_t num, int i) {
    return ((num >> i) & 1);
}

void printBits(size_t tt) {
    for (int i = sizeof(size_t) * 8 - 1; i >= 0; --i) {
        cout << getIthBit(tt, i);
    }
    cout << endl;
}

void printBits(size_t tt, bool inv) {
    for (int i = sizeof(size_t) * 8 - 1; i >= 0; --i) {
        if((getIthBit(tt, i) ^ inv))
            cout << "0";
        else 
            cout << "1";
    }
    cout << endl;
}



// report the cut's root and leaf names
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
getCutCombs(unsigned faninIdx, EcoGate* root, unordered_map<EcoGate*, int>& leaves, vector<EcoCut*>& cuts, vector<vector<EcoCut*>>& faninCutVec, const unsigned& k) {
    unsigned nFanins = faninCutVec.size();
    if(faninIdx == nFanins || leaves.size() > k) {
        if(leaves.size() <= k) {
            EcoCut* cut = new EcoCut(root, leaves);
            cuts.push_back(cut);
            // cut->reportCut();
        }
        return;
    }

    // get the cut from idx'th fanin
    assert(!faninCutVec.empty());
    auto faninCuts = faninCutVec.at(faninIdx);
    for(size_t i=0; i<faninCuts.size(); i++) {
        auto faninCut = faninCuts.at(i);
        const unordered_set<EcoGate*> faninCutLeaves = faninCut->getLeaves();
        for(auto& leaf : faninCutLeaves)
            leaves[leaf]++;
        // leaves.insert(leaves.end(), faninCutLeaves.begin(), faninCutLeaves.end());
        // recursive call
        getCutCombs(faninIdx+1, root, leaves, cuts, faninCutVec, k);
        for(auto& leaf : faninCutLeaves) {
            leaves[leaf]--;
            if(leaves[leaf] == 0)
                leaves.erase(leaf);
        }
    }
    
}

// enumerate k feasible cuts
vector<EcoCut*>
EcoNtk::enumerateCutsRec(const unsigned& k, EcoGate* g) {
    if(g->isGlobalTrav()) return _gate2Cuts.at(g);
    g->setToGlobalTrav();
    
    EcoCut* cut = new EcoCut(g, {g});
    _gate2Cuts[g].push_back(cut);

    // boundary case, return the PI gate when PI is reached
    if(g->getGateType() == EcoGate::ECO_PI_GATE || g->getGateType() == EcoGate::ECO_CONST_0_GATE || g->getGateType() == EcoGate::ECO_CONST_1_GATE)  {
        
        return _gate2Cuts[g];
    }

    // otherwise recursively collect k-feasible cuts
    vector<vector<EcoCut*>> faninCutVec;
    unsigned nFanins = g->getNumFanins();
    // if(g->getGateName() == "y[0]")
    // cout << "fanin ";
    // for(size_t i=0; i<nFanins; i++) {
    //     auto faninCut = enumerateCutsRec(k, g->getFanin(i));
    //     if(g->getGateName() == "y[0]")
    //         cout << g->getFanin(i)->getGateFullName() << " " << faninCut.size() << " ";
    // }
    // if(g->getGateName() == "y[0]")
    // cout << endl;
    for(size_t i=0; i<nFanins; i++) {
        auto faninCut = enumerateCutsRec(k, g->getFanin(i));
        faninCutVec.push_back(faninCut);
        
    }
    
    vector<EcoCut*> faninCutCombs;
    vector<size_t> faninCutIdxs(nFanins);
    unordered_map<EcoGate*, int> leaves;
    getCutCombs(0, g, leaves, faninCutCombs, faninCutVec, k);

    _gate2Cuts[g].insert(_gate2Cuts[g].end(), faninCutCombs.begin(), faninCutCombs.end());

    // cout << "cur " << g->getGateFullName() << " " << _gate2Cuts[g].size() << endl;
    return _gate2Cuts[g];
    
}

// write the cut function as an aag file
void 
EcoNtk::writeCutAag(EcoCut* pCut) {
    for(auto g : pCut->getLeaves()) {
        
    }
}

// compute the truth table of the cut
size_t
EcoNtk::computeCutTT(EcoCut* pCut) {
    const unsigned cutSize = pCut->getCutSize();
    auto root = pCut->getRoot();
    CirGate* rootAigGate = root->getAigNode();
    ofstream cutAagFile;
    size_t patterns[cutSize];
    auto leaves = pCut->getLeaves();
    unordered_set<CirGate*> leafCirGates;

    // set the patterns
    for (int i = 0; i < cutSize; i++) {
        patterns[i] = 0;
    }
    for (int i = 0; i < pow(2, cutSize); i++) {
        unsigned pattern = i;
        for (int j = 0; j < cutSize; j++) {
            patterns[j] += ((pattern & 1) << i);
            pattern >>= 1;
        }
    }
    
    // set the sim pattern at cut leaves
    int i=0;
    for(const auto& g : leaves) {
        CirGate* gAig = g->getAigNode();
        gAig->setPValue(patterns[i++]);
        // printBits(gAig->getPValue()());
        leafCirGates.insert(gAig);
    }

    // sim using the dfs list
    CirMgr* pCirMgr = cirV->getEcoCirV();
    for(const auto& g : pCirMgr->_dfsList) {
        if(leafCirGates.count(g)) continue; // if the gate is a leaf node, don't change its value
        g->pSim();
        if(g == rootAigGate) break; // if the root aig node is reached, break
    }
    
    size_t simVal = rootAigGate->getPValue()();
    // cout << setw(5) << pCut->getRoot()->getGateName() << " ";
    // printBits(simVal, root->getAigNodeInv());
}

}
}
#endif