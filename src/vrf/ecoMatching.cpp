#ifndef ECO_MATCHING_CPP
#define ECO_MATCHING_CPP
#include "ecoNtk.h"
#include "ecoMgr.h"

namespace gv {
namespace eco {

extern bool getIthBit(const size_t num, int i);
extern void printBits(size_t tt);

// do the cut matching from output side
void
EcoMgr::doOutputSideMatching() {
    unsigned nPo = _oldNtk->getNumPos();
    
    for(unsigned i=0; i<nPo; i++) {
        // if(i > 0) break;
        matchOnePo(i);
    }
}

// do the matching for one po
void
EcoMgr::matchOnePo(unsigned ithPo) {
    gv::cir::EcoGate* oldPo = _oldNtk->getPo(ithPo);
    gv::cir::EcoGate* newPo = _newNtk->getPo(ithPo);

    matchCutsAtGatePair(oldPo->getFanin(0), newPo->getFanin(0));

}

// match the cuts at the gate pair
void
EcoMgr::matchCutsAtGatePair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate) {
    auto oldPoCuts = _oldNtk->getGateCuts(oldGate);
    auto newPoCuts = _newNtk->getGateCuts(newGate);
    
    // sort the cuts by their NPN class

    // 1. collect the cuts of the same NPN class
    // here we store the key as <cutsize>_<NPN class>, for convience of sorting by cut size
    map<string, vector<gv::cir::EcoCut*>, greater<string>> oldNPNClass2Cuts;
    map<string, vector<gv::cir::EcoCut*>, greater<string>> newNPNClass2Cuts;

    for(auto cut : oldPoCuts) {
        // if(cut->getCutSize() > 4) continue;
        auto[npnClass, match] = getNPNHash(cut);
        npnClass = to_string(cut->getCutSize()) + "_" + npnClass;
        oldNPNClass2Cuts[npnClass].push_back(cut);
    }

    for(auto cut : newPoCuts) {
        // if(cut->getCutSize() > 4) continue;
        auto[npnClass, match] = getNPNHash(cut);
        npnClass = to_string(cut->getCutSize()) + "_" + npnClass;
        newNPNClass2Cuts[npnClass].push_back(cut);
    }


    bool foundMatch = false;
    // enumerate by old npn class
    for(auto&[oldNPNClass, oldCuts] : oldNPNClass2Cuts) {
        if(!newNPNClass2Cuts.count(oldNPNClass)) continue;
        if(foundMatch) break;
        auto& newCuts = newNPNClass2Cuts.at(oldNPNClass);

        // sort the cut by the number of merged gates
        for(auto& oldCut : oldCuts) {
            unsigned numMerged = 0;
            for(const auto& leaf : oldCut->getLeaves()) {
                if(isMerged(leaf))
                    numMerged++;
            }
            oldCut->setNumMergedLeaves(numMerged);
        }

        for(auto& newCut : newCuts) {
            unsigned numMerged = 0;
            for(const auto& leaf : newCut->getLeaves()) {
                
                if(isMerged(leaf)) {
                    numMerged++;
                }
            }
            newCut->setNumMergedLeaves(numMerged);
        }
        sort(oldCuts.begin(), oldCuts.end(), [](gv::cir::EcoCut* a, gv::cir::EcoCut* b) {
            return a->getNumMergedLeaves() > b->getNumMergedLeaves();
        });
        sort(newCuts.begin(), newCuts.end(), [](gv::cir::EcoCut* a, gv::cir::EcoCut* b) {
            return a->getNumMergedLeaves() > b->getNumMergedLeaves();
        });


        
        for(auto& oldCut : oldCuts) {
            for(auto& newCut : newCuts) {
                if(oldCut->getCutSize() < 2) continue;
                foundMatch = match2Cuts(oldCut, newCut);
                if(foundMatch) break;
            }
            if(foundMatch) break;
        }
    }
}

// use the bits of the size_t to control each constant assignment bits
void
getConstAssignFromSizeT(size_t constAssignSizeT, vector<pair<int, bool>>& constAssignment) {
    for(int i=0; i<constAssignment.size(); i++)
        constAssignment.at(i).second = getIthBit(constAssignSizeT, i);
}

// use the bits of the size_t to control each constant assignment bits, consider pole
void
getConstAssignFromSizeT(size_t constAssignSizeT, size_t invAssignSizeT, vector<pair<int, bool>>& constAssignment) {
    
    size_t n = constAssignment.size();
    for(int i=0; i<n; i++) {
        bool inv = getIthBit(invAssignSizeT, i);
        bool bit = getIthBit(constAssignSizeT, i);

        bit ^= inv; 

        constAssignment.at(i).second = bit;
    }
}

// for cuts with size greater than 4, use simulation to find if there is a valid solution
// if there is no valid solution, the output matching will be returned as -1
vector<size_t>
EcoMgr::simNFindValidMatch(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, vector<int>& comb) {
    extern int factorial(int n);
    vector<size_t> ret;
    const int constAssignSize = comb.size();
    const unsigned cutSize = oldCut->getCutSize();
    const unsigned simSize = cutSize - constAssignSize;
    int numPermutation = factorial(comb.size());
    vector<pair<int, bool>> constAssignmentOld(constAssignSize);
    vector<pair<int, bool>> constAssignmentNew(constAssignSize);


    vector<gv::cir::EcoGate*> oldFreeLeaves, oldLeaves;
    vector<gv::cir::EcoGate*> newFreeLeaves, newLeaves;
    int idx = 0;
    for(const auto& leaf : oldCut->getLeaves()) {
        if(find(comb.begin(), comb.end(), idx) == comb.end())
            oldFreeLeaves.push_back(leaf);
        oldLeaves.push_back(leaf);
        idx++;
    }
    idx = 0;
    for(const auto& leaf : newCut->getLeaves()) {
        if(idx >= constAssignSize)
            newFreeLeaves.push_back(leaf);
        newLeaves.push_back(leaf);
        idx++;
    }

    for(int i=0; i<constAssignSize; i++) {
        constAssignmentOld[i] = {comb[i], false};
        constAssignmentNew[i] = {i, false};
    }

    // try different permutaion
    for(int i=0; i<numPermutation; i++) {
        // try different pole
        for(size_t invAssignSizeT = 0; invAssignSizeT < pow(2, constAssignSize); invAssignSizeT++) {
            unordered_set<size_t> matches;
            // for the old cut, we choose comb[0] as the first fanin, comb[1] as second fanin, and so on...
            // 1. set 00/01/10/11 (all the permutation) to 1/2... fanins, and sim the remaining 4 fanins to get truth table
            for(size_t constAssignSizeT = 0; constAssignSizeT < pow(2, constAssignSize); constAssignSizeT++) {
                getConstAssignFromSizeT(constAssignSizeT, invAssignSizeT, constAssignmentOld); // apply pole change on old cut
                getConstAssignFromSizeT(constAssignSizeT, constAssignmentNew);

                auto oldCutTT = _oldNtk->computeCutTTWithConst(oldCut, constAssignmentOld);
                auto newCutTT = _newNtk->computeCutTTWithConst(newCut, constAssignmentNew);

                auto[oldNpnClass, oldMatches] = _pNpnHash->getNPNHashFull(oldCutTT, simSize);
                auto[newNpnClass, newMatch] = _pNpnHash->getNPNHash(newCutTT, simSize);
                
                // if the NPN class does not match, the match is not valid
                if(oldNpnClass != newNpnClass) break; 
                
                unordered_set<size_t> caseMatches;

                for(const auto& oldMatch : oldMatches) {
                    int n = oldMatch.size();
                    vector<int> newPosVec(n);
                    vector<int> inputMatch(simSize);
                    int outputMatch = -1;

                    outputMatch = (newMatch[0] & 0b1) ^ (oldMatch[0] & 0b1);

                    for(int j=1; j<n; ++j)
                        newPosVec[newMatch[j] / 2] = j;

                    for(int i=1; i<oldMatch.size(); ++i)
                        inputMatch[i - 1] = ((newPosVec[oldMatch[i] / 2] - 1) * 2 + ((newMatch[newPosVec[oldMatch[i] / 2]] & 0b1) ^ (oldMatch.at(i) & 0b1)));
                    
                    auto encode = _pNpnHash->encodeMatch2SizeT(outputMatch, inputMatch);
                    
                    if(matches.empty() || matches.count(encode))
                        caseMatches.insert(encode);
                }
                
                if(matches.empty()) {
                    for(const auto& m : caseMatches)
                        matches.insert(m);
                }
                else {
                    unordered_set<size_t> nextMatches;
                    for(const auto& m : caseMatches) {
                        if(matches.count(m))
                            nextMatches.insert(m);
                    }
                    matches = nextMatches;
                    if(matches.empty())
                        break;
                }
            }
            for(const auto& match : matches) {
                auto[outputMatch, freeInputMatch] = _pNpnHash->decodeEncodedSizeTMatch(match, simSize);
                unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch;

                cout << "output match : " << oldCut->getRoot()->getGateFullName() << " " << newCut->getRoot()->getGateFullName() << " inv " << outputMatch << endl;
                // handle the free leaves
                cout << "free input match : " << endl;
                for(unsigned j=0; j<freeInputMatch.size(); j++) {
                    const auto& m = freeInputMatch.at(j);
                    unsigned oldPos = j;
                    unsigned newPos = m / 2;
                    bool inv = (bool)(m % 2);
                    inputMatch[oldFreeLeaves.at(oldPos)] = {newFreeLeaves.at(newPos), inv};
                    cout << oldFreeLeaves.at(oldPos)->getGateFullName() << " " << newFreeLeaves.at(newPos)->getGateFullName() << " inv " << m % 2 << endl;
                }
                // handle the const leaves
                cout << "combo :" << endl;
                for(int j=0; j<comb.size(); j++) {
                    bool inv = getIthBit(invAssignSizeT, j);
                    inputMatch[oldLeaves.at(comb.at(j))] = {newLeaves.at(j), inv};
                    cout << oldLeaves.at(comb.at(j))->getGateFullName() << " " << newLeaves.at(j)->getGateFullName() << " inv " << inv << endl; 
                }
                assert(checkMatchValid(oldCut, newCut, outputMatch, inputMatch));
                cout << "---------------" << endl;
                ret.push_back(match);
            }
        }

        next_permutation(comb.begin(), comb.end());
    }

    return ret;
}

// check if the matching is valid
bool
EcoMgr::checkMatchValid(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch) {
    assert(oldCut->getCutSize() == newCut->getCutSize());
    const unsigned cutSize = oldCut->getCutSize(); // get the cut size
    unordered_set<gv::cir::CirGate*> oldLeafCirGates, newLeafCirGates;
    gv::cir::CirGate* oldRootAigGate = oldCut->getRoot()->getAigNode();
    gv::cir::CirGate* newRootAigGate = newCut->getRoot()->getAigNode();

    // set the input pattern
    size_t patterns[cutSize];

    // set the patterns
    for (int i = 0; i < cutSize; i++)
        patterns[i] = 0;
    
    for (int i = 0; i < pow(2, cutSize); i++) {
        unsigned pattern = i;
        for (int j = 0; j < cutSize; j++) {
            patterns[j] += ((pattern & 1) << i);
            pattern >>= 1;
        }
    }

    // set the patterns to gate leaves
    unsigned patternIdx = 0;
    for(const auto&[oldGate, newGateInfo] : inputMatch) {
        auto[newGate, inputInv] = newGateInfo;
        gv::cir::CirGate* pAigOld = oldGate->getAigNode();
        gv::cir::CirGate* pAigNew = newGate->getAigNode();

        size_t pattern = patterns[patternIdx++];
        size_t invPattern = (~pattern);
        
        if(inputInv ^ oldGate->getAigNodeInv())
            pAigOld->setPValue(invPattern);
        else
            pAigOld->setPValue(pattern);
        if(newGate->getAigNodeInv())
            pAigNew->setPValue(invPattern);
        else
            pAigNew->setPValue(pattern);
        
        oldLeafCirGates.insert(pAigOld);
        newLeafCirGates.insert(pAigNew);
    }

    // sim using the dfs list
    for(const auto& g : _oldNtk->getAigDfsList()) {
        if(oldLeafCirGates.count(g)) continue; // if the gate is a leaf node, don't change its value
        g->pSim();
        if(g == oldRootAigGate) break; // if the root aig node is reached, break
    }
    for(const auto& g : _newNtk->getAigDfsList()) {
        if(newLeafCirGates.count(g)) continue; // if the gate is a leaf node, don't change its value
        g->pSim();
        if(g == newRootAigGate) break; // if the root aig node is reached, break
    }

    // get the sim value at output side
    size_t oldSimVal, newSimVal;
    oldSimVal = oldRootAigGate->getPValue()();
    newSimVal = newRootAigGate->getPValue()();

    // see if the output mapping is inv
    oldSimVal = oldCut->getRoot()->getAigNodeInv() ? (~oldSimVal) : oldSimVal;
    newSimVal = newCut->getRoot()->getAigNodeInv() ? (~newSimVal) : newSimVal;

    if(outputInv)
        oldSimVal = (~oldSimVal);

    // mask to filter out the unused bits
    size_t mask = 0;
    for (size_t i = 0; i < pow(2, cutSize); ++i)
        mask += ((size_t)1 << i);

    oldSimVal &= mask;
    newSimVal &= mask;

    cout << "output inv " << outputInv << endl;
    printBits(oldSimVal);
    printBits(newSimVal);

    if(oldSimVal != newSimVal) {
        cout << "not matched, please debug for this match" << endl;
        oldCut->reportCut();
        newCut->reportCut();
    }

    return (oldSimVal == newSimVal);
}

// get one valid matching method
pair<int, vector<int>>
EcoMgr::getOneMatchWay(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut) {
    unsigned cutSize = oldCut->getCutSize();

    int outputMatch = -1;
    
    vector<int> inputMatch(cutSize);
    
    // for 4-feasible cut, directly lookup
    if(cutSize <= 4) {
        auto[oldNpnClass, oldMatch] = getNPNHash(oldCut);
        auto[newNpnClass, newMatch] = getNPNHash(newCut);

        int n = oldMatch.size();
        vector<int> newPosVec(n);

        // if they are not of the same NPN class, We can not use it to fix the circuit
        assert(oldNpnClass == newNpnClass);

        outputMatch = (newMatch[0] & 0b1) ^ (oldMatch[0] & 0b1);
        cout << "oo : " << outputMatch << endl;

        for(int j=1; j<n; ++j)
            newPosVec[newMatch[j] / 2] = j;

        for(int i=1; i<oldMatch.size(); ++i) {
            inputMatch[i - 1] = ((newPosVec[oldMatch[i] / 2] - 1) * 2 + ((newMatch[newPosVec[oldMatch[i] / 2]] & 0b1) ^ (oldMatch.at(i) & 0b1)));
            cout << "ii " << i-1 << " " << inputMatch[i - 1] << endl;
        }

        // check that the matching is valid
        unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputGateMatch;
        vector<gv::cir::EcoGate*> oldLeaves, newLeaves;

        for(const auto& leaf : oldCut->getLeaves())
            oldLeaves.push_back(leaf);
        for(const auto& leaf : newCut->getLeaves())
            newLeaves.push_back(leaf);
        for(int i=0; i<inputMatch.size(); i++) {
            inputGateMatch[oldLeaves.at(i)] = {newLeaves.at(inputMatch.at(i) / 2), (bool)(inputMatch.at(i) % 2)};
            cout << "in matxh " << oldLeaves.at(i)->getGateFullName() << " " << newLeaves.at(inputMatch.at(i) / 2)->getGateFullName() << " inv " << (bool)(inputMatch.at(i) % 2) << endl;
        }
        assert(checkMatchValid(oldCut, newCut, outputMatch, inputGateMatch));
    }
    // for 5/6-feasible cut, sim it and compute
    else {
        vector<int> comb; // choose some signals to sim, others can be got by hashing
        int k = cutSize - 4; // since 4 feasible cuts are pre-computed
        
        for(unsigned i=0; i<k; i++) 
            comb.push_back(i);
        
        while (1)
        {
            // check if we can find a valid matching 
            simNFindValidMatch(oldCut, newCut, comb); 
            int i = k-1;
            while (i>=0 && comb[i]==i + cutSize - k)
                i--;
            if (i==-1) assert(0); // should has a solution
            comb[i]++;
            for (int j=i+1; j<k; j++)
                comb[j] = comb[j-1]+1;
            
        }
    }

    return {outputMatch, inputMatch};
}

// match the 2 cuts and record the RP pair if possible
bool
EcoMgr::match2Cuts(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut) {
    vector<gv::cir::EcoGate*> oldLeaves;
    vector<gv::cir::EcoGate*> newLeaves;

    for(const auto& g : oldCut->getLeaves())
        oldLeaves.push_back(g);
    for(const auto& g : newCut->getLeaves())
        newLeaves.push_back(g);
    
    auto[outputMatch, inputMatch] = getOneMatchWay(oldCut, newCut);
    size_t encode = _pNpnHash->encodeMatch2SizeT(outputMatch, inputMatch);
    auto[decodeOutputMatch, decodeInputMatch] = _pNpnHash->decodeEncodedSizeTMatch(encode, oldCut->getCutSize());

    cout << "output match : " << outputMatch << endl;
    vector<pair<gv::cir::EcoGate*, gv::cir::EcoGate*>> candRPPair;
    for(int i=0; i<inputMatch.size(); i++) {
        auto oldGate = oldLeaves.at(i);
        auto newGate = newLeaves.at(inputMatch[i] / 2);
        if(getRPGate(oldGate) && getRPGate(oldGate) != newGate) {
            cout << "assert(0); " << oldGate->getGateFullName() << " " << getRPGate(oldGate)->getGateFullName() << " new " << newGate->getGateFullName() << endl;
            return false;
        }
        candRPPair.push_back({oldGate, newGate});
        cout << oldGate->getGateFullName() << " " << newGate->getGateFullName() << " " << inputMatch[i] % 2 << endl;
    }
    cout << "----------------" << endl;


    // if added successfully, add the cand RP Pair into RP pairs
    for(const auto&[oldGate, newGate] : candRPPair) {
        addRPPair(oldGate, newGate);
    }


    // reorder the matching using symmetry to match the merged gates




    return true;
}

// end of namespace gv::eco
}
// end of namespace gv
}
#endif