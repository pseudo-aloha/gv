#ifndef ECO_MATCHING_CPP
#define ECO_MATCHING_CPP
#include "ecoNtk.h"
#include "ecoMgr.h"

#include <iomanip>

namespace gv {
namespace eco {

extern bool getIthBit(const size_t num, int i);
extern void printBits(size_t tt);

class NPNClassSignature {
public:
    NPNClassSignature(string NPNClass, unsigned cutSize, unsigned numMergedGatesUnion) : _NPNClass(NPNClass), _cutSize(cutSize), _numMergedGatesUnion(numMergedGatesUnion) {}
    ~NPNClassSignature() {}

    string   getNPNClass()            { return _NPNClass; }
    unsigned getCutSize()             { return _cutSize; }
    unsigned getNumMergedGatesUnion() { return _numMergedGatesUnion; }

private:
    string _NPNClass;
    unsigned _cutSize;
    unsigned _numMergedGatesUnion;
};

// iterative method to get combinations
vector<vector<int>> getCombs(int n, int k) {
    vector<vector<int>> rets;
    vector<int> comb;

    for(int i=0; i<k; i++) {
        comb.push_back(i);
    }
    rets.push_back(comb);

    while(1) {
        int i = k - 1;
        while(i >= 0 && comb.at(i) >= i + n - k)
            i--;
            
        if(i < 0) break;
        comb.at(i)++;
        for(int j=i+1; j<k; j++)
            comb.at(j) = comb.at(j-1) + 1;
        rets.push_back(comb);
    }
    return rets;
}

unsigned
EcoMgr::getGatesEqStatus(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate) {
    if(!isMerged(oldGate) || !isMerged(newGate)) return ECO_GATES_NEQ;
    auto[oldMergedAig, oldMergedInv] = getMergedAig(oldGate);
    auto newAig = newGate->getAigNode();
    bool newAigInv = newGate->getAigNodeInv();

    if(newAig != oldMergedAig) return ECO_GATES_NEQ;
    if(oldMergedInv ^ newAigInv) return ECO_GATES_INV_EQ;
    return ECO_GATES_EQ;
}

// add the RP pair
// gate a fixed to gate b can fix fanout #i
void
EcoMgr::addRPPair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, bool inv, unsigned fixedPo) {
    cout << "mermermer : " << oldGate->getGateFullName() << " " << newGate->getGateFullName() << " " << inv << endl;
    // if the rp gate is merged to the old gate and the pole also matches, we don't have to add it in to rp pairs
    if(isMerged(oldGate)) {
        // cout << mergedAig << 
        auto[mergedAig, mergedPole] = getMergedAig(oldGate);
        if(mergedAig == newGate->getAigNode() && ((mergedPole ^ newGate->getAigNodeInv() ^ inv) == false)) {
            if(!_rpTable.at(fixedPo).count(oldGate)) {
                EcoRPInfo* pRPInfo = new EcoRPInfo(oldGate, false);
                assert(!_rpTable.at(fixedPo).count(oldGate));
                _rpTable.at(fixedPo)[oldGate] = pRPInfo;
            }
            return;
        }
    }
    // if the oldgate is not yet be fixed to another gate, simply add it
    if(!_rpTable.at(fixedPo).count(oldGate)) {
        EcoRPInfo* pRPInfo = new EcoRPInfo(newGate, inv);
        assert(!_rpTable.at(fixedPo).count(oldGate));
        _rpTable.at(fixedPo)[oldGate] = pRPInfo;
    }
}

// report the recorded RP pairs
void
EcoMgr::reportRPPair() {
    for(auto poRP : _rpTable) {
        for(auto[oldGate, roInfo] : poRP) {

        }
    }
}

void
EcoMgr::sortCutsByNumMergedGates(vector<gv::cir::EcoCut*>& cuts) {
    sort(cuts.begin(), cuts.end(), [](gv::cir::EcoCut* a, gv::cir::EcoCut* b) {
        return a->getNumMergedLeaves() > b->getNumMergedLeaves();
    });
}

// int getMergedGateScore(const gv::cir::EcoCutInfo*& a) {
//     auto aCut = a->getCut();
//     auto aConst = a->getConstInsertList();
//     int aMergedNum = aCut->getNumMergedLeaves();
//     if(!aConst.empty()) {
//         for(const auto& aCon : aConst) {
//             int constAssignedLeafIdx = aCon / 2;
//             auto constAssignedLeaf = aCut->getLeaf(constAssignedLeafIdx);
//             if(constAssignedLeaf->isMerged() && !constAssignedLeaf->isConstGate())
//                 --aMergedNum;
//         }
//     }
//     return aMergedNum;
// }

// sort the cuts by the number of the merged gates
void
EcoMgr::sortCutsByNumMergedGates(vector<gv::cir::EcoCutInfo*>& cuts) {
    sort(cuts.begin(), cuts.end(), [](gv::cir::EcoCutInfo*& a, gv::cir::EcoCutInfo*& b) {
        int aMergedNum = a->getNumMerged();
        int bMergedNum = b->getNumMerged();

        return (aMergedNum > bMergedNum);
    });
}

// sort the NPN class by the size of union of merged gates
vector<string>
EcoMgr::sortNPNClass(const unordered_map<string, vector<gv::cir::EcoCutInfo*>>& NPNClass2Cuts) {
    vector<string> sortedNPNClass;
    vector<NPNClassSignature*> sortedNPNClassSignatures;
    
    for(auto&[NPNClass, cutsInfo] : NPNClass2Cuts) {
        unsigned cutSize = cutsInfo.front()->getCut()->getCutSize();
        
        // sortedNPNClass.push_back(NPNClass);
        unordered_set<gv::cir::EcoGate*> mergedGateSt;
        for(const auto& cutInfo : cutsInfo) {
            auto cut = cutInfo->getCut();
            auto constInsert = cutInfo->getConstInsertList();
            for(const auto& leaf : cut->getLeaves()) {
                if(isMerged(leaf))
                    mergedGateSt.insert(leaf);
            }
        }
        NPNClassSignature* NPNSig = new NPNClassSignature(NPNClass, cutSize, mergedGateSt.size());
        sortedNPNClassSignatures.push_back(NPNSig);
    }

    sort(sortedNPNClassSignatures.begin(), sortedNPNClassSignatures.end(), [](NPNClassSignature* a, NPNClassSignature* b) {
        // if(a->getCutSize() > b->getCutSize())
        //     return true;
        return a->getNumMergedGatesUnion() > b->getNumMergedGatesUnion();
    });

    for(const auto& NPNSig : sortedNPNClassSignatures)
        sortedNPNClass.push_back(NPNSig->getNPNClass());
    
    for(const auto& NPNSig : sortedNPNClassSignatures)
        delete NPNSig;

    return sortedNPNClass;
}

// check if the 2 output cone are equivalent
bool
EcoMgr::check2ConeEq(int ithPo) {
    extern bool isNtkEq(Abc_Ntk_t* pNtkOld, Abc_Ntk_t* pNtkNew);
    gv::cir::EcoNtk::rewriteDesign(getOldDesignName());
    Abc_Ntk_t* pNtkOld = Io_Read("./tmp.v", IO_FILE_VERILOG, 0, 0 );
    gv::cir::EcoNtk::rewriteDesign(getNewDesignName());
    Abc_Ntk_t* pNtkNew = Io_Read("./tmp.v", IO_FILE_VERILOG, 0, 0 );
    Abc_Obj_t* oldPo = Abc_NtkCo(pNtkOld, ithPo);
    Abc_Obj_t* newPo = Abc_NtkCo(pNtkNew, ithPo);
    cout << "old po name : " << Abc_ObjName(oldPo) << endl;
    Abc_Ntk_t * pNtkOldCone = Abc_NtkCreateCone( pNtkOld, Abc_ObjFanin0(oldPo), Abc_ObjName(oldPo), 1 );
    Abc_Ntk_t * pNtkNewCone = Abc_NtkCreateCone( pNtkNew, Abc_ObjFanin0(newPo), Abc_ObjName(newPo), 1 );
    if ( Abc_ObjFaninC0(oldPo) ) Abc_ObjXorFaninC( Abc_NtkPo(pNtkOldCone, 0), 0 );
    if ( Abc_ObjFaninC0(newPo) ) Abc_ObjXorFaninC( Abc_NtkPo(pNtkNewCone, 0), 0 );
    
    return isNtkEq(pNtkOldCone, pNtkNewCone);
}


// do the cut matching from output side
void
EcoMgr::doOutputSideMatching() {
    unsigned nPo = _oldNtk->getNumPos();
    _rpTable.resize(nPo);
    
    for(unsigned i=0; i<nPo; ++i) {
        cout << "matching po : " << _oldNtk->getPo(i)->getGateName() << endl;
        if(check2ConeEq(i)) {
            cout << "po pair is eq" << endl;
            continue;
        }
        matchOnePo(i);
    }

    // 1. build selector
    // buildSelector();

    // 2. choose rp pair

}

// do the matching for one po
void
EcoMgr::matchOnePo(unsigned ithPo) {
    gv::cir::EcoGate* oldPo = _oldNtk->getPo(ithPo)->getFanin(0);
    gv::cir::EcoGate* newPo = _newNtk->getPo(ithPo)->getFanin(0);

    // match cuts
    // recursive matching, until merged frontier is reached
    unordered_set<gv::cir::EcoGate*> visited;
    queue<pair<gv::cir::EcoGate*, gv::cir::EcoGate*>> q;
    unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> candRp;
    
    // use po as initial rp point
    q.push(make_pair(oldPo, newPo));
    candRp[oldPo] = make_pair(newPo, false);

    while(!q.empty()) {
        auto[oldGate, newGate] = q.front();
        q.pop();
        cout << "matching_pair: " << oldGate->getGateFullName() << " " << newGate->getGateFullName() << endl;
        if(visited.count(oldGate) || isUnderMergeFrontierSet(oldGate)) break;
        visited.insert(oldGate);

        auto rp = matchCutsAtGatePair(oldGate, newGate, ithPo, true);

        if(!rp.empty()) {
            candRp.erase(oldGate);
            for(auto[og, ngInfo] : rp) {
                candRp[og] = ngInfo;
                auto[ng, inv] = ngInfo;
                q.push(make_pair(og, ng));
            }
        }
        // break;

        // if the gate is visited before or merged frontier is reached, stop matching
        // if(isMerged(oldGate)) break;
    }

    for(auto[og, ngInfo] : candRp) {
        auto[ng, inv] = ngInfo;
        addRPPair(og, ng, inv, ithPo);
    }
}

// match the cuts at the gate pair (now used for po matching)
unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>>
EcoMgr::matchCutsAtGatePair(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, int ithPo, bool doConstInsert) {
    // get the cuts rooted at gate pair
    auto oldPoCuts = _oldNtk->getGateCuts(oldGate);
    auto newPoCuts = _newNtk->getGateCuts(newGate);

    // sort the cuts by their NPN class

    // 1. collect the cuts of the same NPN class
    // here we store the key as <cutsize>_<NPN class>, for convience of sorting by cut size
    unordered_map<string, vector<gv::cir::EcoCutInfo*>> oldNPNClass2Cuts;
    unordered_map<string, vector<gv::cir::EcoCutInfo*>> newNPNClass2Cuts;

    // compute the cuts signature and sort by #merged gates
    computeCutsSignatures(oldNPNClass2Cuts, oldPoCuts, doConstInsert);
    computeCutsSignatures(newNPNClass2Cuts, newPoCuts);

    // sort the NPN classes by
    // 1. size of union of merged gates in the cuts in the NPN class
    // 2. #vars of the NPN class (i.e. the cut size in the NPN class)
    auto sortedOldNPNClass = sortNPNClass(oldNPNClass2Cuts);
    auto sortedNewNPNClass = sortNPNClass(newNPNClass2Cuts);


    bool foundMatch = false;
    // enumerate by old npn class
    for(const auto& oldNPNClass  : sortedOldNPNClass) {
        auto oldCutsInfo = oldNPNClass2Cuts.at(oldNPNClass);
        if(!newNPNClass2Cuts.count(oldNPNClass)) continue;
        if(foundMatch) break;
        auto& newCutsInfo = newNPNClass2Cuts.at(oldNPNClass);
        cout << "NPN class : " << oldNPNClass << endl;
        cout << "old cuts size : " << oldCutsInfo.size() << " new cuts size : " << newCutsInfo.size() << endl;
        // assert(oldCutsInfo.size() * newCutsInfo.size() < 10000);
        int oldCutSize = oldCutsInfo.at(0)->getCut()->getCutSize();
        int maxDist = max(oldCutsInfo.size(), newCutsInfo.size()); // pick the larger one
        int maxIt = 100;
        int numIt = 0;
        bool exceededMaxIt = false;

        // collect the cuts and sort them by number of merged gates in
        
        unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> candMatch;
        // get the corresponding merged gates in the old circuit
        unordered_set<gv::cir::EcoGate*> mergedGateSet;
        for(int newIdx = 0; newIdx < newCutsInfo.size(); ++newIdx) {
            auto newCut = newCutsInfo.at(newIdx)->getCut();
            auto newConstInsert = newCutsInfo.at(newIdx)->getConstInsertList();
            for(int i=0; i<newCut->getCutSize(); ++i) {
                auto leaf = newCut->getLeaf(i);
                auto mergedGates = getMergedGates(leaf);
                auto invMergedGates = getInvMergedGates(leaf);
                for(auto& mGate : mergedGates)
                    mergedGateSet.insert(mGate);
                for(auto& invMGate : invMergedGates)
                    mergedGateSet.insert(invMGate);
            }
        }

        for(int oldIdx = 0; oldIdx < oldCutsInfo.size(); ++oldIdx) {
            auto oldCutInfo = oldCutsInfo.at(oldIdx);
            auto oldCut = oldCutInfo->getCut();
            auto oldConstInsert = oldCutInfo->getConstInsertList();
            unordered_set<int> constInsertIdxs;
            for(auto ci : oldConstInsert) {
                int idx = ci / 2;
                constInsertIdxs.insert(idx);
            }
            unsigned numMerged = 0;
            for(int i=0; i<oldCut->getCutSize(); ++i) {
                auto leaf = oldCut->getLeaf(i);
                if(constInsertIdxs.count(i)) { 
                    continue; // we don't want to count the const inserted gate
                }
                
                if(mergedGateSet.count(leaf))
                    ++numMerged;
            }
            oldCutInfo->setNumMerged(numMerged);
        }
        // sort by the num of merged gates
        sortCutsByNumMergedGates(oldCutsInfo);

        for(int dist = 0; dist < maxDist; ++dist) {
            for(int oldIdx = 0; oldIdx <= dist && oldIdx < oldCutsInfo.size(); ++oldIdx) {
                auto oldCut = oldCutsInfo.at(oldIdx)->getCut();
                auto oldConstInsert = oldCutsInfo.at(oldIdx)->getConstInsertList();
                for(int newIdx = 0; oldIdx + newIdx <= dist && newIdx < newCutsInfo.size(); ++newIdx) {
                    auto newCut = newCutsInfo.at(newIdx)->getCut();
                    auto newConstInsert = newCutsInfo.at(newIdx)->getConstInsertList();
                    auto[matchSucess, match] = match2Cuts(oldCut, newCut, oldConstInsert, ithPo);
                    if(matchSucess)
                        return match;
                    if(++numIt >= maxIt) {
                        exceededMaxIt = true;
                        break;
                    }
                }
                if(exceededMaxIt) break;
            }
            if(exceededMaxIt) break;
        }
    }

    // if there is no rp pair found, fix it from po
    unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> emptyMatch;

    return emptyMatch;
}


// check if the truth table is simply constant
bool isTruthTableConst(size_t tt, int size) {
    bool allZero = (getIthBit(tt, 0) == false);
    for(int i=1; i<pow(2, size); ++i) {
        if(allZero && (getIthBit(tt, i) == true)) {
            return false;
        }
        if(!allZero && (getIthBit(tt, i) == false)) {
            return false;
        }
    }
    return true;
}

// compute the signature and sort them by NPN class for the provided cuts
// 
void
EcoMgr::computeCutsSignatures(unordered_map<string, vector<gv::cir::EcoCutInfo*>>& NPNClass2Cuts, vector<gv::cir::EcoCut *>& cuts, bool doConstInsert) {
    for(auto& cut : cuts) {
        // handle the original cuts
        auto[npnClass, match] = getNPNHash(cut);
        unsigned numMerged = 0;
        ConstInsertList constInsertList; // empty for the original cut
        npnClass = npnClass;
        gv::cir::EcoCutInfo* cutInfo = new gv::cir::EcoCutInfo(cut, constInsertList);
        NPNClass2Cuts[npnClass].push_back(cutInfo);

        for(const auto& leaf : cut->getLeaves()) {
            if(isMerged(leaf))
                ++numMerged;
        }
        cut->setNumMergedLeaves(numMerged);

        // handle constant insertion cuts
        // skip if not to do constant insertion
        if(doConstInsert) {
            assert(cut->getRoot()->isInOldCircuit());
            unsigned cutSize = cut->getCutSize();
            for(int nConst = 1; nConst < (cutSize+1) / 2; nConst++) {
                vector<ConstInsertList> combs = getCombs(cutSize, nConst);
                for(const auto& comb : combs) {
                    for(size_t bitMask = 0; bitMask < pow(2, nConst); bitMask++) {
                        constInsertList.clear();
                        vector<pair<int, bool>> constAssignmentOld;
                        
                        
                        for(int bitIdx = 0; bitIdx < nConst; bitIdx++) {
                            // cout << (comb.at(bitIdx) * 2 + (int)getIthBit(bitMask, bitIdx)) << " ";
                            constInsertList.push_back(comb.at(bitIdx) * 2 + (int)getIthBit(bitMask, bitIdx));
                        }
                        for(const auto& constInsert : constInsertList) {
                            // cout << "ddd " << constInsert / 2 << " " << constInsert % 2 << endl;
                            constAssignmentOld.push_back(make_pair(constInsert / 2, constInsert % 2));
                        }
                        const unsigned simSize = cutSize - nConst;
                        // cout << "size : " << constAssignmentOld.size() << " " << nConst << endl;
                        auto oldCutTT = _oldNtk->computeCutTTWithConst(cut, constAssignmentOld);
                        if(isTruthTableConst(oldCutTT, simSize)) // block the trivial cuts
                            continue;
                        for(const auto&[constIdx, constVal] : constAssignmentOld) {
                            cout << "assign " << cut->getLeaf(constIdx)->getGateName() << " " << constVal << endl;
                        }
                        printBits(oldCutTT);
                        auto[constInsertNpnClass, _] = _pNpnHash->getNPNHash(oldCutTT, simSize);
                        cutInfo = new gv::cir::EcoCutInfo(cut, constInsertList);
                        NPNClass2Cuts[constInsertNpnClass].push_back(cutInfo);
                        
                        cout << "const insert NPN class : " << constInsertNpnClass << endl;
                        // cout << endl;
                    }
                }
            }
        }
    }
    for(auto&[NPNClass, NPNCuts] : NPNClass2Cuts) {
        sortCutsByNumMergedGates(NPNCuts);
        cout << "NPNClass " << NPNClass << " " << NPNCuts.size() << endl;
    }
}

// use the bits of the size_t to control each constant assignment bits
void
getConstAssignFromSizeT(size_t constAssignSizeT, vector<pair<int, bool>>& constAssignment, const int constAssignSize) {
    size_t n = constAssignment.size();
    for(int i=n-constAssignSize; i<n; ++i)
        constAssignment.at(i).second = getIthBit(constAssignSizeT, i);
}

// use the bits of the size_t to control each constant assignment bits, consider pole
void
getConstAssignFromSizeT(size_t constAssignSizeT, size_t invAssignSizeT, vector<pair<int, bool>>& constAssignment, const int constAssignSize) {
    size_t n = constAssignment.size();
    for(int i=n-constAssignSize; i<n; ++i) {
        bool inv = getIthBit(invAssignSizeT, i);
        bool bit = getIthBit(constAssignSizeT, i);

        bit ^= inv; 

        constAssignment.at(i).second = bit;
    }
}

// for cuts with size greater than 4, use simulation to find if there is a valid solution
// if there is no valid solution, the output matching will be returned as -1
vector<size_t>
EcoMgr::simNFindValidMatch(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, vector<int>& comb, const ConstInsertList& ConstInsert) {
    extern int factorial(int n);
    vector<size_t> ret;
    const int constAssignSize = comb.size(); // this const assign size does not include const insert
    const unsigned cutSize = oldCut->getCutSize();
    const unsigned simSize = cutSize - constAssignSize - ConstInsert.size();
    int numPermutation = factorial(comb.size());
    unordered_map<gv::cir::EcoGate*, int> gate2IdxOld ,gate2IdxNew;
    unordered_set<gv::cir::EcoGate*> constAssignGateSet;
    unordered_map<gv::cir::EcoGate*, bool> constInsertGateMap; // map const inserted gate to the mapped pole

    // here we merge the const assignments of 5/6 sim and const assertion algo
    vector<pair<int, bool>> constAssignmentOld(constAssignSize);
    vector<pair<int, bool>> constAssignmentNew(constAssignSize);

    for(const auto& idx : comb)
        constAssignGateSet.insert(oldCut->getLeaf(idx));
    for(const auto& ci : ConstInsert) {
        int idx = ci / 2;
        bool val = ci % 2;
        constInsertGateMap[oldCut->getLeaf(idx)] = val;
        constAssignGateSet.insert(oldCut->getLeaf(idx));
        constAssignmentOld.push_back(make_pair(idx, val));
    }

    vector<gv::cir::EcoGate*> oldFreeLeaves, oldLeaves;
    vector<gv::cir::EcoGate*> newFreeLeaves, newLeaves;
    int idx = 0;
    for(const auto& leaf : oldCut->getLeaves()) {
        if(!constAssignGateSet.count(leaf)) // the leaf is not const inserted
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
                // Do const insert sim for the extra inputs (more than 4)
                getConstAssignFromSizeT(constAssignSizeT, invAssignSizeT, constAssignmentOld, constAssignSize); // apply pole change on old cut
                getConstAssignFromSizeT(constAssignSizeT, constAssignmentNew, constAssignSize);

                // cout << "pppp " << oldCut->getCutSize() << " " << newCut->getCutSize() << " " << constAssignmentOld.size() << endl;

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
            assert(ConstInsert.empty() || matches.empty());
            for(const auto& match : matches) {
                auto[outputMatch, freeInputMatch] = _pNpnHash->decodeEncodedSizeTMatch(match, simSize);
                unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch;
                vector<int> inputMatchIntVec;
                inputMatchIntVec.resize(cutSize);
                if(outputMatch) continue; // TODO : we only find the cuts that do not need output inv for output cut (for now), which means we do not allow invert at root (only NP cuts)
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
                    assert(oldLeaves.at(gate2IdxOld.at(oldFreeLeaves.at(oldPos))) == oldFreeLeaves.at(oldPos));
                    assert(newLeaves.at(gate2IdxNew.at(newFreeLeaves.at(newPos))) == newFreeLeaves.at(newPos));
                    cout << oldFreeLeaves.at(oldPos)->getGateFullName() << " " << newFreeLeaves.at(newPos)->getGateFullName() << " inv " << m % 2 << " " << gate2IdxOld.at(oldFreeLeaves.at(oldPos)) << " " << gate2IdxNew.at(newFreeLeaves.at(newPos)) << endl;
                }
                // handle the const leaves for 5/6 sim
                cout << "combo :" << endl;
                for(int j=0; j<comb.size(); ++j) {
                    bool inv = getIthBit(invAssignSizeT, j);
                    inputMatch[oldLeaves.at(comb.at(j))] = {newLeaves.at(j), inv};
                    inputMatchIntVec[comb.at(j)] = j * 2 + inv;
                    assert(inputMatchIntVec[comb.at(j)] < cutSize);
                    cout << oldLeaves.at(comb.at(j))->getGateFullName() << " " << newLeaves.at(j)->getGateFullName() << " inv " << inv << " " << comb.at(j) << " " << j << endl; 
                }
                // todo:handle const insertion algo
                if(!ConstInsert.empty()) {
                    cout << "const insert : " << endl;
                    for(const auto& ci : ConstInsert) {
                        int idx = ci / 2;
                        int val = ci % 2;
                        // inputMatch[oldLeaves.at(comb.at(j))] = {newLeaves.at(j), val};
                        inputMatchIntVec[idx] = 2 * EcoNPNHash::CONST0SIZETENCODE + val;
                        cout << oldCut->getLeaf(idx)->getGateFullName() << " " << val << endl;
                    }
                }

                assert(checkMatchValidWithConst(oldCut, newCut, outputMatch, inputMatch, ConstInsert, {}));
                
                size_t encode = _pNpnHash->encodeMatch2SizeT(outputMatch, inputMatchIntVec);
                auto[_outputMatch, _inputMatch] = _pNpnHash->decodeEncodedSizeTMatch(encode, cutSize);
                // for(unsigned j=0; j<_inputMatch.size(); ++j) {
                //     const auto& m = _inputMatch.at(j);
                //     unsigned oldPos = j;
                //     unsigned newPos = m / 2;
                //     bool inv = (bool)(m % 2);
                //     // cout << oldLeaves.at(oldPos)->getGateFullName() << " " << newLeaves.at(newPos)->getGateFullName() << endl;
                //     // cout << inputMatch[oldLeaves.at(oldPos)].first->getGateFullName() << " " << newLeaves.at(newPos)->getGateFullName() << endl;
                //     // assert(inputMatch[oldLeaves.at(oldPos)].first == newLeaves.at(newPos));
                //     // assert(inputMatch[oldLeaves.at(oldPos)].second == inv);
                // }
                ret.push_back(encode);
            }
        }

        next_permutation(comb.begin(), comb.end());
    }
assert(ConstInsert.empty() || ret.empty());
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
EcoMgr::checkMatchValidWithConst(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, int outputInv, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> inputMatch, const ConstInsertList& constInsertOld, const ConstInsertList& constInsertNew) {
    const unsigned oldCutSize = oldCut->getCutSize() - constInsertOld.size(); // the valid cut size of the old cut
    const unsigned newCutSize = newCut->getCutSize() - constInsertNew.size();
    assert(oldCutSize == newCutSize);
    const unsigned cutSize = oldCut->getCutSize(); // get the cut size
    const unsigned simSize = oldCutSize;
    unordered_set<gv::cir::CirGate*> oldLeafCirGates, newLeafCirGates;
    unordered_map<gv::cir::EcoGate*, bool> constAssignmentOld, constAssignmentNew;

    for(const auto& cInsert : constInsertOld) {
        unsigned idx = cInsert / 2;
        bool val = cInsert % 2;
        constAssignmentOld[oldCut->getLeaf(idx)] = val;
    }
    for(const auto& cInsert : constInsertNew) {
        unsigned idx = cInsert / 2;
        bool val = cInsert % 2;
        constAssignmentNew[newCut->getLeaf(idx)] = val;
    }
    
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
        
        if(!oldGate->isConstGate()) {
            if(inputInv ^ oldGate->getAigNodeInv())
                pAigOld->setPValue(invPattern);
            else
                pAigOld->setPValue(pattern);
        }
        if(!newGate->isConstGate()) {
            if(newGate->getAigNodeInv())
                pAigNew->setPValue(invPattern);
            else
                pAigNew->setPValue(pattern);
        }
        
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
// TODO hangle const insert for 5/6 cuts
vector<size_t>
EcoMgr::getMatchWays(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, const ConstInsertList& ConstInsert) {
    vector<size_t> ret; // return the different match ways in encoded form
    const unsigned cutSize = oldCut->getCutSize();
    const unsigned constInsertSize = ConstInsert.size();
    const unsigned freeCutSize = cutSize - constInsertSize; // cut size for free inputs of the cut
    unordered_map<gv::cir::EcoGate*, bool> constInsertGateMap; // map const inserted gate to the mapped pole

    int outputMatch = -1;
    
    vector<int> inputMatch(cutSize - constInsertSize); // this is the input matching without constant insertion
    vector<int> inputMatchWithConst(cutSize); // this one adds the constant insertion things

    for(const auto& cInsert : ConstInsert) {
        int idx = cInsert / 2;
        bool val = cInsert % 2;
        constInsertGateMap[oldCut->getLeaf(idx)] = val;
    }
    
    // for 4-feasible cut, directly lookup
    if(freeCutSize <= 4) {
        auto[oldNpnClass, oldMatches] = getNPNHashFull(oldCut, ConstInsert);
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
            int inputMatchIdx = 0;

             for(const auto& cInsert : ConstInsert) {
                int idx = cInsert / 2;
                bool val = cInsert % 2;
                cout << oldCut->getLeaf(idx)->getGateFullName() << " " << val << endl;
                // constInsertGateMap[oldCut->getLeaf(idx)] = val;
            }
            for(int oldLeafIdx = 0; oldLeafIdx < cutSize; ++oldLeafIdx) {
                if(constInsertGateMap.count(oldLeaves.at(oldLeafIdx))) {
                    bool constVal = constInsertGateMap.at(oldLeaves.at(oldLeafIdx));
                    // 6 for const 0, 7 for const 1
                    inputMatchWithConst[oldLeafIdx] = 2 * EcoNPNHash::CONST0SIZETENCODE + constVal;
                    inputGateMatch[oldLeaves.at(oldLeafIdx)] = make_pair(_newNtk->getConst0Gate(), constVal);
                    // cout << "loloha " << oldLeaves.at(oldLeafIdx)->getGateFullName() << " " << constVal << endl;
                    // oldLeafIdx++;
                    continue;
                }
                inputMatchWithConst[oldLeafIdx] = inputMatch.at(inputMatchIdx);
                inputGateMatch[oldLeaves.at(oldLeafIdx)] = {newLeaves.at(inputMatch.at(inputMatchIdx) / 2), (bool)(inputMatch.at(inputMatchIdx) % 2)};
                cout << "in match " << oldLeaves.at(oldLeafIdx)->getGateFullName() << " " << newLeaves.at(inputMatch.at(inputMatchIdx) / 2)->getGateFullName() << " inv " << (bool)(inputMatch.at(inputMatchIdx) % 2) << endl;
                // oldLeafIdx++;
                inputMatchIdx++;
            }
            cout << "m : ";
            for(auto& m : inputMatchWithConst) {
                cout << m << " ";
            }
            cout << endl;
            assert(checkMatchValidWithConst(oldCut, newCut, outputMatch, inputGateMatch, ConstInsert, {}));

            // encode it and push it into return value
            size_t encode = _pNpnHash->encodeMatch2SizeT(outputMatch, inputMatchWithConst);
            ret.push_back(encode);
        }
    }
    // for 5/6-feasible cut, sim it and compute
    else {
        vector<int> comb;    // choose some signals to sim, others can be got by hashing
        int k = cutSize - 4 - ConstInsert.size(); // since 4 feasible cuts are pre-computed, and we also need to minus const inserted part
        
        // insert 0~k-1 as the first comb
        for(unsigned i=0; i<k; ++i) 
            comb.push_back(i);

        while (1)
        {
            // check if the comb conflicts with constant insertion
            bool conflict = false;
            for(const auto& c : comb) {
                if(constInsertGateMap.count(oldCut->getLeaf(c))) {
                    conflict = true;
                    break;
                }
            }

            // check if we can find a valid matching
            if(!conflict) {
                auto foundMatches = simNFindValidMatch(oldCut, newCut, comb, ConstInsert); 
                ret.insert(ret.end(), foundMatches.begin(), foundMatches.end());
            }
            
            // compute next combination
            int i = k-1;
            while (i>=0 && comb[i]==i + cutSize - k)
                i--;
            if (i==-1) break; // should has a solution
            ++comb[i];
            for (int j=i+1; j<k; ++j)
                comb[j] = comb[j-1]+1;
            
        }
    }
    // for(unsigned i=0; i<ret.size(); ++i) {
    //     auto[outputMatch, inputMatch] = _pNpnHash->decodeEncodedSizeTMatch(ret[i], cutSize);
    //     for(unsigned j=0; j<cutSize; ++j)
    //         assert(inputMatch.at(j) < 14);
    // }
    return ret;
}


// TODO : add contant matching
// match the 2 cuts and record the RP pair if possible
// ith po is given when matching a particular po
pair<bool, unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>>>
EcoMgr::match2Cuts(gv::cir::EcoCut* oldCut, gv::cir::EcoCut* newCut, const ConstInsertList& ConstInsert, int ithPo) {
    const unsigned oldCutFreeSize = oldCut->getCutSize() - ConstInsert.size(); // the cut size without const insertion
    const unsigned newCutFreeSize = newCut->getCutSize();
    assert(oldCutFreeSize == newCutFreeSize);
    unordered_map<gv::cir::EcoGate*, pair<gv::cir::EcoGate*, bool>> candRPPair;
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
    
    auto matchWays = getMatchWays(oldCut, newCut, ConstInsert);
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
    cout << "num match ways : " << matchWays.size() << " " << ConstInsert.size() << endl;
    size_t candMatch = 0;
    int bestScore = -1; // number of matched
    for(unsigned i=0; i<matchWays.size(); ++i) {
        int score = 0;
        auto[outputMatch, inputMatch] = _pNpnHash->decodeEncodedSizeTMatch(matchWays.at(i), oldCut->getCutSize());
        for(unsigned j=0; j<inputMatch.size(); ++j) {
            auto oldGate = oldLeaves.at(j);
            auto m = inputMatch.at(j);

            // skip if it is a const inserted gate
            if(m >= 2 * EcoNPNHash::CONST0SIZETENCODE)
                continue;

            if(oldGate2MergedNewGate.count(oldGate)) {
                auto[newMergedGate, mergeInv] = oldGate2MergedNewGate.at(oldGate);
                auto newGate = newLeaves.at(m / 2);
                auto inv = m % 2;
                if(newMergedGate == newGate && mergeInv == inv)
                    score++;
            }
        }
        if(score > bestScore) {
            cout << "best score " << score << " " << i << endl;
            bestScore = score;
            candMatch = matchWays.at(i);
        }
    }

    // if no valid match is found, return false
    if(bestScore < 0) return make_pair(false, candRPPair);
    
    
    // find a match way that can maximally match the merged gates
    auto[outputMatch, inputMatch] = _pNpnHash->decodeEncodedSizeTMatch(candMatch, oldCut->getCutSize());
    cout << "output match : " << outputMatch << " const insert size " << ConstInsert.size() << endl;
    if(outputMatch) return make_pair(false, candRPPair); // we don't want invert at po, TODO : check if we have to change to NPN
    for(int i=0; i<inputMatch.size(); ++i) {
        auto oldGate = oldLeaves.at(i);
        gv::cir::EcoGate* newGate = nullptr;
        unsigned matchedIdx = inputMatch.at(i) / 2;
        bool matchedInv = inputMatch.at(i) % 2;

        if(matchedIdx >= EcoNPNHash::CONST0SIZETENCODE) {
            newGate = _newNtk->getConst0Gate();
        }
        else {
            newGate = newLeaves.at(inputMatch.at(i) / 2);
        }

        bool inv = inputMatch[i] % 2;

        if(isUnderMergeFrontierSet(oldGate)) {
            auto eqStatus = getGatesEqStatus(oldGate, newGate);
            if(!((eqStatus == ECO_GATES_INV_EQ && inv) || (eqStatus == ECO_GATES_EQ && !inv))) {
                candRPPair.clear();
                return make_pair(false, candRPPair);
            }
        }

        candRPPair[oldGate] = {newGate, inv};
        cout << "fanin match : " << oldGate->getGateFullName() << " " << newGate->getGateFullName() << " " << inv << endl;
        
        // if matching is with respect to a particular, compute the cosine similarity
        if(ithPo >= 0) {    
            cout << "cos sim : " << getCosineSimilarity(oldGate, newGate, ithPo) << endl;
        }
    }
    cout << "kkk" << endl;

    return make_pair(true, candRPPair);
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


// end of namespace gv::eco
}
// end of namespace gv
}
#endif