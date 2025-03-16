#ifndef ECO_CUT_CPP
#define ECO_CUT_CPP

#include "ecoNtk.h"
#include "ecoMgr.h"
#include "base/abc/abc.h"
#include <iostream>
#include <string>
#include <cassert>
#include <iomanip>

unsigned gv::cir::EcoCut::_maxCutsPerNode = 100;

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
    cout << "root : " << _root->getGateFullName() << endl;
    cout << "leaves : ";
    for(const auto& leaf : _leaves) {
        cout << leaf->getGateFullName() << " ";
    }
    cout << endl;
}


// to avoid redundant cuts
bool
EcoNtk::computeAndInsertSigature(EcoCut* pCut) {
    assert(pCut->getSig().empty());
    
    string sig;
    sig += pCut->getRoot()->getGateName();

    vector<string> faninNames;
    for(const auto& leaf : pCut->getLeaves())
        faninNames.push_back(leaf->getGateName());

    sort(faninNames.begin(), faninNames.end()); // sort to guarantee irredundancy

    for(const auto& name : faninNames)
        sig += name;

    if(_enumeratedSignatures.count(sig))
        return false;

    pCut->setSig(sig);
    _enumeratedSignatures.insert(sig);

    return true;
}

// enumerate k feasible cuts
void
EcoNtk::enumerateCuts(unsigned k, gv::eco::EcoMgr* pEco) {
    EcoGate::setGlobalTrav();
    // traverse from POs
    for(size_t i=0; i<getNumPos(); i++) {
        enumerateCutsRec(k, getPo(i), pEco);
    }
}

bool
EcoNtk::checkCut(EcoCut* pCut) {
    for(const auto& leaf : pCut->getLeaves()) {
        if(!checkCutRec(pCut, leaf, true))
            return false;
    }
    return true;
}

bool
EcoNtk::checkCutRec(EcoCut* pCut, EcoGate* g, bool isSelf) {
    if(pCut->getLeaves().count(g) && !isSelf)
        return false;
    if(g->getGateType() == EcoGate::ECO_PI_GATE || g->getGateType() == EcoGate::ECO_CONST_0_GATE || g->getGateType() == EcoGate::ECO_CONST_1_GATE)
        return true;

    for(size_t i=0; i<g->getNumFanins(); i++) {
        auto fanin = g->getFanin(i);
        if(!checkCutRec(pCut, fanin, false)) return false;
    }
    return true;
}

void
EcoNtk::getCutCombs(unsigned faninIdx, EcoGate* root, unordered_map<EcoGate*, int>& leaves, vector<EcoCut*>& cuts, vector<vector<EcoCut*>>& faninCutVec, const unsigned& k, gv::eco::EcoMgr* pEco) {
    unsigned nFanins = faninCutVec.size();
    if(cuts.size() > EcoCut::getMaxCutsPerNode()) return; // too many cuts for a single node
    if(faninIdx == nFanins || leaves.size() > k) {
        if(leaves.size() <= k) {
            EcoCut* cut = new EcoCut(root, leaves);
            // if(!checkCut(cut))
            //     cut->reportCut();
            if(checkCut(cut) && computeAndInsertSigature(cut)) {
                cuts.push_back(cut);
            }
            else
                delete cut;
        }
        return;
    }

    // get the cut from idx'th fanin
    assert(!faninCutVec.empty());
    auto faninCuts = faninCutVec.at(faninIdx);
    for(size_t i=0; i<faninCuts.size(); i++) {
        // cout << "fanin " << faninIdx << " " << " i "<<i << endl;
        auto faninCut = faninCuts.at(i);
        const unordered_set<EcoGate*> faninCutLeaves = faninCut->getLeaves();
        unordered_set<string> names;
        names.clear();
        for(auto& leaf : faninCutLeaves) {
            leaves[leaf]++;
            names.insert(leaf->getGateName());
        }
        // recursive call
        getCutCombs(faninIdx+1, root, leaves, cuts, faninCutVec, k, pEco);
        for(auto& leaf : faninCutLeaves) {
            leaves[leaf]--;
            if(leaves[leaf] == 0) {
                leaves.erase(leaf);
            }
        }
    }
    
}

// enumerate k feasible cuts
vector<EcoCut*>
EcoNtk::enumerateCutsRec(const unsigned& k, EcoGate* g, gv::eco::EcoMgr* pEco) {
    if(g->isGlobalTrav()) return _gate2Cuts.at(g);
    g->setToGlobalTrav();
    
    EcoCut* cut = new EcoCut(g, {g});
    _gate2Cuts[g].push_back(cut);

    // boundary case, return the PI gate when PI is reached
    if(g->getGateType() == EcoGate::ECO_PI_GATE || g->getGateType() == EcoGate::ECO_CONST_0_GATE || g->getGateType() == EcoGate::ECO_CONST_1_GATE || pEco->isMerged(g))  {
        
        return _gate2Cuts[g];
    }

    // otherwise recursively collect k-feasible cuts
    vector<vector<EcoCut*>> faninCutVec;
    unsigned nFanins = g->getNumFanins();
    unsigned combs = 1;
    for(size_t i=0; i<nFanins; i++) {
        auto faninCut = enumerateCutsRec(k, g->getFanin(i), pEco);
        faninCutVec.push_back(faninCut);
        // sort(faninCut.begin(), faninCut.end(), [](EcoCut* a, EcoCut* b) {
        //     return a->getCutSize() < b->getCutSize();
        // }); //sort the cuts from small to big
    }



    // for(size_t i=0; i<nFanins; i++)
    // cout << "fanin " << g->getFanin(i)->getGateFullName() << " " << _gate2Cuts.at(g->getFanin(i)).size() << endl;
    // cout << "cur " << g->getGateName() << endl;
    // cout << "num of combs for cur : " << combs << endl;
    
    vector<EcoCut*> faninCutCombs;
    vector<size_t> faninCutIdxs(nFanins);
    unordered_map<EcoGate*, int> leaves;
    getCutCombs(0, g, leaves, faninCutCombs, faninCutVec, k, pEco);

    _gate2Cuts[g].insert(_gate2Cuts[g].end(), faninCutCombs.begin(), faninCutCombs.end());

    // cout << "cur " << g->getGateFullName() << " " << _gate2Cuts[g].size() << endl;
    return _gate2Cuts[g];
    
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
        size_t pattern = patterns[i++];
        if(g->getAigNodeInv())
            gAig->setPValue(~pattern);
        else
            gAig->setPValue(pattern);
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

    // mask to filter out the unused bits
    size_t mask = 0;
    for (size_t i = 0; i < pow(2, cutSize); ++i)
        mask += ((size_t)1 << i);
    
    if(root->getAigNodeInv())
        simVal = (~simVal);

    simVal &= mask;

    return simVal;
}

// compute the truth table of the cut with constant inserted
size_t
EcoNtk::computeCutTTWithConst(EcoCut* pCut, const vector<pair<int, bool>>& constAssignment) {
    const unsigned cutSize = pCut->getCutSize();
    const unsigned simSize = cutSize - constAssignment.size(); // the actual number of fanin we need to sim
    auto root = pCut->getRoot();
    CirGate* rootAigGate = root->getAigNode();
    ofstream cutAagFile;
    size_t patterns[simSize];
    auto leaves = pCut->getLeaves();
    unordered_set<CirGate*> leafCirGates;
    unordered_map<int, bool> constAssignmentMap;

    for(const auto& assign : constAssignment)
        constAssignmentMap[assign.first] = assign.second;

    // set the patterns
    for (int i = 0; i < simSize; i++)
        patterns[i] = 0;
    
    for (int i = 0; i < pow(2, simSize); i++) {
        unsigned pattern = i;
        for (int j = 0; j < simSize; j++) {
            patterns[j] += ((pattern & 1) << i);
            pattern >>= 1;
        }
    }
    
    // set the sim pattern at cut leaves
    int patIdx=0;
    int leafIdx = 0;
    const size_t all0 = 0;
    const size_t all1 = 0xffffffffffffffff;

    for(const auto& g : leaves) {
        CirGate* gAig = g->getAigNode();
        if(constAssignmentMap.count(leafIdx)) {
            if(constAssignmentMap.at(leafIdx) == false)
                gAig->setPValue(all0);
            else
                gAig->setPValue(all1);
        }
        else
            gAig->setPValue(patterns[patIdx++]);

        leafCirGates.insert(gAig);
        leafIdx++;
    }

    // sim using the dfs list
    CirMgr* pCirMgr = cirV->getEcoCirV();
    for(const auto& g : pCirMgr->_dfsList) {
        if(leafCirGates.count(g)) continue; // if the gate is a leaf node, don't change its value
        g->pSim();
        if(g == rootAigGate) break; // if the root aig node is reached, break
    }
    
    size_t simVal = rootAigGate->getPValue()();

    // mask to filter out the unused bits
    size_t mask = 0;
    for (size_t i = 0; i < pow(2, simSize); ++i)
        mask += ((size_t)1 << i);
    
    simVal &= mask;

    return simVal;
}

}
}
#endif