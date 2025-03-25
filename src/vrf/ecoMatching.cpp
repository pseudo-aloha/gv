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
    
    for(unsigned i=0; i<nPo; ++i) {
        matchOnePo(i);
        break;
    }
}

// do the matching for one po
void
EcoMgr::matchOnePo(unsigned ithPo) {
    gv::cir::EcoGate* oldPo = _oldNtk->getPo(ithPo);
    gv::cir::EcoGate* newPo = _newNtk->getPo(ithPo);

    // fault analysis; for the PO
    

    // match cuts
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
        auto[npnClass, match] = getNPNHash(cut);
        npnClass = to_string(cut->getCutSize()) + "_" + npnClass;
        oldNPNClass2Cuts[npnClass].push_back(cut);
    }

    for(auto cut : newPoCuts) {
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
                    ++numMerged;
            }
            oldCut->setNumMergedLeaves(numMerged);
        }

        for(auto& newCut : newCuts) {
            unsigned numMerged = 0;
            for(const auto& leaf : newCut->getLeaves()) {
                
                if(isMerged(leaf)) {
                    ++numMerged;
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
                // if(oldCut->getCutSize() < 2) continue;
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
    for(int i=0; i<constAssignment.size(); ++i)
        constAssignment.at(i).second = getIthBit(constAssignSizeT, i);
}

// use the bits of the size_t to control each constant assignment bits, consider pole
void
getConstAssignFromSizeT(size_t constAssignSizeT, size_t invAssignSizeT, vector<pair<int, bool>>& constAssignment) {
    
    size_t n = constAssignment.size();
    for(int i=0; i<n; ++i) {
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
    unordered_map<gv::cir::EcoGate*, int> gate2IdxOld ,gate2IdxNew;
    
    vector<pair<int, bool>> constAssignmentOld(constAssignSize);
    vector<pair<int, bool>> constAssignmentNew(constAssignSize);


    vector<gv::cir::EcoGate*> oldFreeLeaves, oldLeaves;
    vector<gv::cir::EcoGate*> newFreeLeaves, newLeaves;
    int idx = 0;
    for(const auto& leaf : oldCut->getLeaves()) {
        if(find(comb.begin(), comb.end(), idx) == comb.end())
            oldFreeLeaves.push_back(leaf);
        oldLeaves.push_back(leaf);
        gate2IdxOld[leaf] = idx;
        ++idx;
    }
    idx = 0;
    for(const auto& leaf : newCut->getLeaves()) {
        if(idx >= constAssignSize)
            newFreeLeaves.push_back(leaf);
        newLeaves.push_back(leaf);
        gate2IdxNew[leaf] = idx;
        ++idx;
    }

    

    // try different permutaion
    for(int i=0; i<numPermutation; ++i) {
        for(int i=0; i<constAssignSize; ++i) {
            constAssignmentOld[i] = {comb[i], false};
            constAssignmentNew[i] = {i, false};
        }
        // try different pole
        for(size_t invAssignSizeT = 0; invAssignSizeT < pow(2, constAssignSize); ++invAssignSizeT) {
            unordered_set<size_t> matches;
            // for the old cut, we choose comb[0] as the first fanin, comb[1] as second fanin, and so on...
            // 1. set 00/01/10/11 (all the permutation) to 1/2... fanins, and sim the remaining 4 fanins to get truth table
            for(size_t constAssignSizeT = 0; constAssignSizeT < pow(2, constAssignSize); ++constAssignSizeT) {
                getConstAssignFromSizeT(constAssignSizeT, invAssignSizeT, constAssignmentOld); // apply pole change on old cut
                getConstAssignFromSizeT(constAssignSizeT, constAssignmentNew);

                auto oldCutTT = _oldNtk->computeCutTTWithConst(oldCut, constAssignmentOld);
                auto newCutTT = _newNtk->computeCutTTWithConst(newCut, constAssignmentNew);

                auto[oldNpnClass, oldMatches] = _pNpnHash->getNPNHashFull(oldCutTT, simSize);
                auto[newNpnClass, newMatch] = _pNpnHash->getNPNHash(newCutTT, simSize);
                
                // if the NPN class does not match, the match is not valid
                if(oldNpnClass != newNpnClass) {
                    matches.clear();
                    break; 
                } 
                
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
                    
                    if(constAssignSizeT == 0 || matches.count(encode))
                        caseMatches.insert(encode);
                }
                
                if(constAssignSizeT == 0) {
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
                    // if there is already no match exists, break
                    if(matches.empty())
                        break;
                            
                }
                // check each submatch is valid, debug use
                // unordered_map<gv::cir::EcoGate*, bool> constAssignmentOld;
                // unordered_map<gv::cir::EcoGate*, bool> constAssignmentNew;
                // for(const auto& match : matches) {
                //     unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch;
                //     auto[outputMatch, freeInputMatch] = _pNpnHash->decodeEncodedSizeTMatch(match, simSize);
                //     for(unsigned j=0; j<freeInputMatch.size(); ++j) {
                //         const auto& m = freeInputMatch.at(j);
                //         unsigned oldPos = j;
                //         unsigned newPos = m / 2;
                //         bool inv = (bool)(m % 2);
                //         inputMatch[oldFreeLeaves.at(oldPos)] = {newFreeLeaves.at(newPos), inv};
                //         cout << oldFreeLeaves.at(oldPos)->getGateFullName() << " " << newFreeLeaves.at(newPos)->getGateFullName() << " inv " << m % 2 << endl;
                //     }
                //     for(unsigned j=0; j<constAssignSize; ++j) {
                //         bool bit = getIthBit(constAssignSizeT, j);
                //         bool oldInv = getIthBit(invAssignSizeT, j);
                //         constAssignmentOld[oldLeaves.at(comb[j])] = (bit ^ oldInv);
                //         constAssignmentNew[newLeaves.at(j)] = bit;
                //         cout << oldLeaves.at(comb[j])->getGateFullName() << " " << newLeaves.at(j)->getGateFullName() << " " << oldInv << endl;
                        
                //     }
                //     assert(checkMatchValidWithConst(oldCut, newCut, outputMatch, inputMatch, constAssignmentOld, constAssignmentNew));
                // }
            }
            for(const auto& match : matches) {
                auto[outputMatch, freeInputMatch] = _pNpnHash->decodeEncodedSizeTMatch(match, simSize);
                unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch;
                vector<int> inputMatchIntVec;
                inputMatchIntVec.resize(cutSize);

                cout << "match " << match << endl;
                cout << "output match : " << oldCut->getRoot()->getGateFullName() << " " << newCut->getRoot()->getGateFullName() << " inv " << outputMatch << endl;
                // handle the free leaves
                cout << "free input match : " << endl;
                for(unsigned j=0; j<freeInputMatch.size(); ++j) {
                    const auto& m = freeInputMatch.at(j);
                    unsigned oldPos = j;
                    unsigned newPos = m / 2;
                    bool inv = (bool)(m % 2);
                    inputMatch[oldFreeLeaves.at(oldPos)] = {newFreeLeaves.at(newPos), inv};
                    inputMatchIntVec[gate2IdxOld.at(oldFreeLeaves.at(oldPos))] = gate2IdxNew.at(newFreeLeaves.at(newPos)) * 2 + inv;
                    cout << oldFreeLeaves.at(oldPos)->getGateFullName() << " " << newFreeLeaves.at(newPos)->getGateFullName() << " inv " << m % 2 << endl;
                }
                // handle the const leaves
                cout << "combo :" << endl;
                for(int j=0; j<comb.size(); ++j) {
                    bool inv = getIthBit(invAssignSizeT, j);
                    inputMatch[oldLeaves.at(comb.at(j))] = {newLeaves.at(j), inv};
                    inputMatchIntVec[comb.at(j)] = j * 2 + inv;
                    cout << oldLeaves.at(comb.at(j))->getGateFullName() << " " << newLeaves.at(j)->getGateFullName() << " inv " << inv << endl; 
                }
                assert(checkMatchValid(oldCut, newCut, outputMatch, inputMatch));
                cout << "---------------" << endl;
                size_t encode = _pNpnHash->encodeMatch2SizeT(outputMatch, inputMatchIntVec);
                ret.push_back(encode);
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
    for (int i = 0; i < cutSize; ++i)
        patterns[i] = 0;
    
    for (int i = 0; i < pow(2, cutSize); ++i) {
        unsigned pattern = i;
        for (int j = 0; j < cutSize; ++j) {
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

    
    if(oldSimVal != newSimVal) {
        cout << "not matched, please debug for this match" << endl;
        cout << "output inv " << outputInv << endl;
        printBits(oldSimVal);
        printBits(newSimVal);

        oldCut->reportCut();
        newCut->reportCut();
        assert(0);
    }

    return (oldSimVal == newSimVal);
}


bool
EcoMgr::checkMatchValidWithConst(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch, unordered_map<gv::cir::EcoGate*, bool>& constAssignmentOld, unordered_map<gv::cir::EcoGate*, bool>& constAssignmentNew) {
    assert(oldCut->getCutSize() == newCut->getCutSize());
    const unsigned cutSize = oldCut->getCutSize(); // get the cut size
    const unsigned simSize = cutSize - constAssignmentOld.size();
    unordered_set<gv::cir::CirGate*> oldLeafCirGates, newLeafCirGates;
    gv::cir::CirGate* oldRootAigGate = oldCut->getRoot()->getAigNode();
    gv::cir::CirGate* newRootAigGate = newCut->getRoot()->getAigNode();

    // set the input pattern (without constants)
    size_t patterns[simSize];

    // set the patterns (without constants)
    for (int i = 0; i < simSize; ++i)
        patterns[i] = 0;
    
    for (int i = 0; i < pow(2, simSize); ++i) {
        unsigned pattern = i;
        for (int j = 0; j < simSize; ++j) {
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

    // set the constant bits
    const size_t all0 = 0;
    const size_t all1 = 0xffffffffffffffff;
    for(const auto&[g, val] : constAssignmentOld) {
        gv::cir::CirGate* pAig = g->getAigNode();
        if(g->getAigNodeInv() ^ val)
            pAig->setPValue(all1);
        else
            pAig->setPValue(all0);
        oldLeafCirGates.insert(pAig);
    }
    for(const auto&[g, val] : constAssignmentNew) {
        gv::cir::CirGate* pAig = g->getAigNode();
        if(g->getAigNodeInv() ^ val)
            pAig->setPValue(all1);
        else
            pAig->setPValue(all0);
        newLeafCirGates.insert(pAig);
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

    // mask to filter out the unused bits (without countsing const bits)
    size_t mask = 0;
    for (size_t i = 0; i < pow(2, simSize); ++i)
        mask += ((size_t)1 << i);

    oldSimVal &= mask;
    newSimVal &= mask;

    return (oldSimVal == newSimVal);
}

// get one valid matching method
vector<size_t>
EcoMgr::getMatchWays(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut) {
    vector<size_t> ret; // return the different match ways in encoded form
    unsigned cutSize = oldCut->getCutSize();

    int outputMatch = -1;
    
    vector<int> inputMatch(cutSize);
    
    // for 4-feasible cut, directly lookup
    if(cutSize <= 4) {
        auto[oldNpnClass, oldMatches] = getNPNHashFull(oldCut);
        auto[newNpnClass, newMatch] = getNPNHash(newCut);

        for(const auto& oldMatch : oldMatches) {
            int n = oldMatch.size();
            vector<int> newPosVec(n);
            
            // if they are not of the same NPN class, We can not use it to fix the circuit
            assert(oldNpnClass == newNpnClass);

            outputMatch = (newMatch[0] & 0b1) ^ (oldMatch[0] & 0b1);

            for(int j=1; j<n; ++j)
                newPosVec[newMatch[j] / 2] = j;

            for(int i=1; i<oldMatch.size(); ++i) {
                inputMatch[i - 1] = ((newPosVec[oldMatch[i] / 2] - 1) * 2 + ((newMatch[newPosVec[oldMatch[i] / 2]] & 0b1) ^ (oldMatch.at(i) & 0b1)));
            }

            

            // check that the matching is valid
            unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputGateMatch;
            vector<gv::cir::EcoGate*> oldLeaves, newLeaves;

            for(const auto& leaf : oldCut->getLeaves())
                oldLeaves.push_back(leaf);
            for(const auto& leaf : newCut->getLeaves())
                newLeaves.push_back(leaf);
            for(int i=0; i<inputMatch.size(); ++i) {
                inputGateMatch[oldLeaves.at(i)] = {newLeaves.at(inputMatch.at(i) / 2), (bool)(inputMatch.at(i) % 2)};
                cout << "in match " << oldLeaves.at(i)->getGateFullName() << " " << newLeaves.at(inputMatch.at(i) / 2)->getGateFullName() << " inv " << (bool)(inputMatch.at(i) % 2) << endl;
            }
            assert(checkMatchValid(oldCut, newCut, outputMatch, inputGateMatch));

            // encode it and push it into return value
            size_t encode = _pNpnHash->encodeMatch2SizeT(outputMatch, inputMatch);
            ret.push_back(encode);
        }
    }
    // for 5/6-feasible cut, sim it and compute
    else {
        vector<int> comb; // choose some signals to sim, others can be got by hashing
        int k = cutSize - 4; // since 4 feasible cuts are pre-computed
        
        for(unsigned i=0; i<k; ++i) 
            comb.push_back(i);

        while (1)
        {
            // check if we can find a valid matching 
            auto foundMatches = simNFindValidMatch(oldCut, newCut, comb); 
            ret.insert(ret.end(), foundMatches.begin(), foundMatches.end());
            int i = k-1;
            while (i>=0 && comb[i]==i + cutSize - k)
                i--;
            if (i==-1) break; // should has a solution
            ++comb[i];
            for (int j=i+1; j<k; ++j)
                comb[j] = comb[j-1]+1;
            
        }
    }

    return ret;
}

// match the 2 cuts and record the RP pair if possible
bool
EcoMgr::match2Cuts(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut) {
    vector<gv::cir::EcoGate*> oldLeaves;
    vector<gv::cir::EcoGate*> newLeaves;

    unordered_map<gv::cir::EcoGate*, int> gate2IdxOld ,gate2IdxNew;
    
    int idx = 0;
    for(const auto& leaf : oldCut->getLeaves()) {
        oldLeaves.push_back(leaf);
        gate2IdxOld[leaf] = idx;
        ++idx;
    }
    idx = 0;
    for(const auto& leaf : newCut->getLeaves()) {
        newLeaves.push_back(leaf);
        gate2IdxNew[leaf] = idx;
        ++idx;
    }
    
    auto matchWays = getMatchWays(oldCut, newCut);
    // size_t encode = _pNpnHash->encodeMatch2SizeT(outputMatch, inputMatch);

    unordered_map<gv::cir::CirGate*, pair<gv::cir::EcoGate*, bool>> newAig2MergedOldGate;
    unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> oldGate2MergedNewGate;
    for(const auto& oldGate : oldCut->getLeaves()) {
        if(!isMerged(oldGate)) continue;
        auto[newAig, inv] = getMergedAig(oldGate);
        newAig2MergedOldGate[newAig] = {oldGate, inv};
    }
    for(const auto& newGate : newCut->getLeaves()) {
        auto newAig = newGate->getAigNode();
        if(newAig2MergedOldGate.count(newAig)) {
            auto[oldGate, inv] = newAig2MergedOldGate.at(newAig);
            inv ^= newGate->getAigNodeInv();
            oldGate2MergedNewGate[oldGate] = {newGate, inv};
            // cout << "old " << gate2IdxOld.at(oldGate) << " new " << gate2IdxNew.at(newGate) << " inv " << inv << endl;
        }
    }

    size_t candMatch;
    int bestScore = -1; // number of matched
    for(unsigned i=0; i<matchWays.size(); ++i) {
        int score = 0;
        auto[outputMatch, inputMatch] = _pNpnHash->decodeEncodedSizeTMatch(matchWays.at(i), oldCut->getCutSize());
        for(unsigned j=0; j<inputMatch.size(); ++j) {
            auto oldGate = oldLeaves.at(j);
            auto m = inputMatch.at(j);
            if(oldGate2MergedNewGate.count(oldGate)) {
                auto[newMergedGate, mergeInv] = oldGate2MergedNewGate.at(oldGate);
                auto newGate = newLeaves.at(m / 2);
                auto inv = m % 2;
                if(newMergedGate == newGate && mergeInv == inv)
                    score++;
            }
        }
        if(score > bestScore) {
            bestScore = score;
            candMatch = matchWays.at(i);
        }
    }
    
    
    // find a match way that can maximally match the merged gates
    auto[outputMatch, inputMatch] = _pNpnHash->decodeEncodedSizeTMatch(candMatch, oldCut->getCutSize());
    cout << "output match : " << outputMatch << endl;
    vector<pair<gv::cir::EcoGate*, gv::cir::EcoGate*>> candRPPair;
    for(int i=0; i<inputMatch.size(); ++i) {
        auto oldGate = oldLeaves.at(i);
        auto newGate = newLeaves.at(inputMatch[i] / 2);
        cout << "cos sim : " << getCosieSimilarity(oldGate, newGate, 0) << endl;
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

// compare function for comparing cut score
// 1. compare by # of merged gates
// 2. compare by cut size
bool cutScoreCmp(gv::cir::EcoCut* a, gv::cir::EcoCut* b) {
    if(a->getNumMergedLeaves() != b->getNumMergedLeaves())
        return a->getNumMergedLeaves() > b->getNumMergedLeaves();
    return (a->getCutSize() > b->getCutSize());
}

// the merged gate matches 1 point, the pole matches, another point
int getSigMatchedScore(const vector<pair<string, bool>>& sig1, const vector<pair<string, bool>>& sig2) {
    int ret = 0;
    size_t i=0, j=0;
    while(i<sig1.size() && j<sig2.size()) {
        auto[gid1, pole1] = sig1.at(i);
        auto[gid2, pole2] = sig2.at(j);
        if(gid1 != gid2) {
            if(gid1 < gid2) ++i;
            else ++j;
        }
        else {
            if(pole1 == pole2) ++ret;
            ++ret; ++i; ++j;
        }
    }

    return ret;
}

void
EcoMgr::sortCandCutsByScore() {
    // record NPN class 2 the cuts
    unordered_map<string, vector<gv::cir::EcoCut*>> NPNClass2CutsOld;
    unordered_map<string, vector<gv::cir::EcoCut*>> NPNClass2CutsNew;

    // collect the cuts
    for(size_t i=0; i<_oldNtk->getNumGates(); ++i) {
        const auto g = _oldNtk->getGate(i);
        const auto cuts = _oldNtk->getGateCuts(g);
        for(const auto& cut : cuts) {
            if(cut->getCutSize() <= 1) continue;
            _oldCandCuts.push_back(cut);
            vector<pair<string, bool>> mgAigSigVec;
            string mgAigSig;
            for(const auto& leaf : cut->getLeaves()) {
                if(isMerged(leaf))
                    mgAigSigVec.push_back({to_string(leaf->getAigNode()->getGid()), leaf->getAigNodeInv()});
                sort(mgAigSigVec.begin(), mgAigSigVec.end());
            }
            cut->setMgAigSig(mgAigSigVec);
        }
    }
    for(size_t i=0; i<_newNtk->getNumGates(); ++i) {
        const auto g = _newNtk->getGate(i);
        const auto cuts = _newNtk->getGateCuts(g);
        for(const auto& cut : cuts) {
            if(cut->getCutSize() <= 1) continue;
            _newCandCuts.push_back(cut);
            vector<pair<string, bool>> mgAigSigVec;
            string mgAigSig;
            for(const auto& leaf : cut->getLeaves()) {
                if(isMerged(leaf)) {
                    auto[mergedAig, pole] = getMergedAig(leaf);
                    mgAigSigVec.push_back({to_string(mergedAig->getGid()), pole});
                }
                sort(mgAigSigVec.begin(), mgAigSigVec.end());
            }
            cut->setMgAigSig(mgAigSigVec);
        }
    }
    cout << "# old cuts : " << _oldCandCuts.size() << endl;
    cout << "# new cuts : " << _newCandCuts.size() << endl;


    // sort by score
    // sort the cut by the number of merged gates
    for(auto& oldCut : _oldCandCuts) {
        unsigned numMerged = 0;
        for(const auto& leaf : oldCut->getLeaves()) {
            if(isMerged(leaf))
                ++numMerged;
        }
        oldCut->setNumMergedLeaves(numMerged);
        auto[npnClass, match] = getNPNHash(oldCut);
        oldCut->setNPNClass(npnClass);
        NPNClass2CutsOld[npnClass].push_back(oldCut);
    }

    for(auto& newCut : _newCandCuts) {
        unsigned numMerged = 0;
        for(const auto& leaf : newCut->getLeaves()) {
            
            if(isMerged(leaf)) {
                ++numMerged;
            }
        }
        newCut->setNumMergedLeaves(numMerged);
        auto[npnClass, match] = getNPNHash(newCut);
        newCut->setNPNClass(npnClass);
        NPNClass2CutsNew[npnClass].push_back(newCut);
    }

    sort(_oldCandCuts.begin(), _oldCandCuts.end(), cutScoreCmp);
    sort(_newCandCuts.begin(), _newCandCuts.end(), cutScoreCmp);

    // for(const auto& oldCut : _oldCandCuts) {
    //     cout << "#merged " << oldCut->getNumMergedLeaves() << endl;
    //     oldCut->reportCut();
    // }
    // choose a cut from new cir cuit and start matching
    for(const auto& newCut : _newCandCuts) {
        const string newNPNClass = newCut->getNPNClass();

        // find the old cuts of the same NPN class as the new cut
        int bestScore = -1;
        gv::cir::EcoCut* chosenOldCand;
        if(NPNClass2CutsOld.count(newNPNClass)) {
            newCut->reportCut();
            auto sameNPNOldCuts = NPNClass2CutsOld.at(newNPNClass);
            for(const auto& oldCut : sameNPNOldCuts) {
                int score = getSigMatchedScore(oldCut->getMgAigSig(), newCut->getMgAigSig());
                if(score > bestScore) {
                    cout << "score : " << score << endl;
                    oldCut->reportCut();
                    chosenOldCand = oldCut;
                    bestScore = score;
                }
            }
            cout << "new npn " << newCut->getNPNClass() << " old npn " << chosenOldCand->getNPNClass() << endl;
            // find the match ways
            match2Cuts(chosenOldCand, newCut);
            break;
        }
        
        

        // use signature to find a most competitive cut from candidates
        // int bestScore = -1;
        // gv::cir::EcoCut* chosenOldCand;
        // for(const auto& oldCut : _oldCandCuts) {
        //     int score = getSigMatchedScore(oldCut->getMgAigSig(), newCut->getMgAigSig());
        //     // cout << "score : " << score << endl;
        //     if(score > bestScore) {
        //         cout << "score : " << score << endl;
        //         oldCut->reportCut();
        //         chosenOldCand = oldCut;
        //         bestScore = score;
                
        //     }
        // }
        // chosenOldCand->reportCut();
        // break;
        
    }
}

// end of namespace gv::eco
}
// end of namespace gv
}
#endif