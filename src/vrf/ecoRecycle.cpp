#ifndef ECO_RECYCLE_CPP
#define ECO_RECYCLE_CPP
#include "ecoNtk.h"
#include "ecoMgr.h"
#include <unordered_set>

namespace gv {
namespace eco {

// collect floating gates in the old circuit
// there are used to fix the
void
EcoMgr::collectFloatingGates() {
    unordered_set<gv::cir::EcoGate*> usedGates;
    gv::cir::EcoGate::setGlobalTrav();

    // do output side
    
    for(unsigned i=0; i<_rpTable.size(); ++i) {
        auto rpForIthPo = _rpTable.at(i);
        queue<gv::cir::EcoGate*> q;
        q.push(_oldNtk->getPo(i)->getFanin(0));
        while(!q.empty()) {
            auto cur = q.front();
            q.pop();

            if(cur->isGlobalTrav()) continue;
            cur->setToGlobalTrav();

            // fixed to a new circuit
            if(isFixedToAnotherGate(cur))
                cur = _finalRpPair.at(cur).first;

            // merged back to the old circuit
            if(cur->isInNewCircuit() && isMerged(cur)) {
                auto[mergedGate, mergedPole] = getOneMergedGate(cur, false);
                cur = mergedGate;
            }
            
            // if the gate is in the old cir, mark as used
            if(cur->isInOldCircuit())
                usedGates.insert(cur);

            for(unsigned j=0; j<cur->getNumFanins(); ++j) {
                auto fanin = cur->getFanin(j);
                if(!fanin->isGlobalTrav())
                    q.push(fanin);
            }
        }
    }
    for(unsigned i=0; i<_oldNtk->getNumGates(); ++i) {
        auto g = _oldNtk->getGate(i);
        if(g->isConstGate()) continue; // no need to reuse const gate

        if(!usedGates.count(g)) {
            cout << "float gate : " << g->getGateFullName() << endl;
            _floatingGates.insert(g);
        }
        else
            cout << "used " << g->getGateFullName() << endl;
    }
}

// match the floating gates to minimize cost
void
EcoMgr::matchFloatingGates() {
    vector<gv::cir::EcoCut*>  candCuts;
    vector<gv::cir::EcoCut*>  patchCuts;
    // enumerate the cuts in floating gates
    for(const auto& root : _floatingGates) {
        auto cuts = _oldNtk->getGateCuts(root);
        for(const auto& cut : cuts) {
            if(cut->getLeaves().size() < 2) continue; // no need to collect single leaf cuts
            auto cutConeGates = cut->collectCutConeGate();
            bool validCut = true;
            for(const auto& g : cutConeGates) {
                if(!_floatingGates.count(g)) {
                    validCut = false;
                    break;
                }
            }
            if(validCut)
                candCuts.push_back(cut);
        }
        // cout << "floating cuts size " << cuts.size() << endl;
    }

    cout << "cand cuts :" << endl;
    for(const auto& cut : candCuts) {
        cut->reportCut();
    }

    // match the cuts to the patch logic
    gv::cir::EcoNtk* patchNtk = new gv::cir::EcoNtk;
    patchNtk->readNtkFile("patchResyn.v");
    patchNtk->enumerateCuts(6, this); // enumerate 6 feasible cuts in patch circuit
    for(unsigned i=0; i<patchNtk->getNumGates(); ++i) {
        auto g = patchNtk->getGate(i);
        auto cuts = patchNtk->getGateCuts(g);

        for(const auto& cut : cuts) {
            unsigned numMerged = 0;
            for(const auto& leaf : cut->getLeaves()) {
                auto corresOldGate = _oldNtk->getGateByName(leaf->getGateName());
                if(corresOldGate && isMerged(corresOldGate))
                    ++numMerged;
            }
            cut->setNumMergedLeaves(numMerged);
            patchCuts.push_back(cut);
        }
    }
    sortCutsByNumMergedGates(candCuts);
    sortCutsByNumMergedGates(patchCuts);
}

}}

#endif