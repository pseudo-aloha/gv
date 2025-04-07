#ifndef ECO_SIM_CPP
#define ECO_SIM_CPP

#include "ecoMgr.h"
#include "util.h"

using namespace gv::eco;

// extern void printBits(size_t tt);

// count the bits of 1's in size_t in constant time
// modified from 32 bits version from https://web.archive.org/web/20151229003112/http://blogs.msdn.com/b/jeuge/archive/2005/06/08/hakmem-bit-count.aspx
unsigned BitCount(size_t u) {
    unsigned uCount1, uCount2;
    unsigned u1, u2;
    unsigned mask = 0xffffffff;
    u1 = (unsigned)u & mask;
    u2 = (unsigned)(u >> 32) & mask;

    uCount1 = u1 - ((u1 >> 1) & 033333333333) - ((u1 >> 2) & 011111111111);
    uCount1 = ((uCount1 + (uCount1 >> 3)) & 030707070707) % 63;

    uCount2 = u2 - ((u2 >> 1) & 033333333333) - ((u2 >> 2) & 011111111111);
    uCount2 = ((uCount2 + (uCount2 >> 3)) & 030707070707) % 63;

    return uCount1 + uCount2;
}


// conduct random simulation on old/new circuit
void
EcoMgr::doRandomSim() {
    extern void printBits(size_t tt);
    assert(_oldNtk->getNumPis() == _newNtk->getNumPis());
    unsigned nPis = _oldNtk->getNumPis();
    size_t** patterns = new size_t*[getNumSim()];

    // generate simulation patterns for PIs
    for(unsigned patId=0; patId<getNumSim(); ++patId) {
        patterns[patId] = new size_t[nPis];
        for(unsigned piId=0; piId<nPis; ++piId) {
            patterns[patId][piId] = size_t(0);
            for (int i = 0; i < 4; ++i) {
                patterns[patId][piId] += size_t(size_t(rnGen(1 << 16)) << i * 16);
            }
        }
    }

    _oldNtk->simOnPats(patterns, getNumSim());
    _newNtk->simOnPats(patterns, getNumSim());
    
    // free the pointers
    for(unsigned patId=0; patId<getNumSim(); ++patId)
        delete patterns[patId];
    delete patterns;
}



// get the cosine similarity of two gates with respect to the po id
// using O(nPats)
double
EcoMgr::getCosineSimilarity(gv::cir::EcoGate* oldGate, gv::cir::EcoGate* newGate, unsigned poId) {
    extern void printBits(size_t tt);
    gv::cir::EcoGate* oldPo = _oldNtk->getPo(poId)->getFanin(0);
    gv::cir::EcoGate* newPo = _newNtk->getPo(poId)->getFanin(0);

    // cout << "nPats : " << oldPo->getSimVal().size() << " " << oldGate->getSimVal().size() << " " << newGate->getSimVal().size() << " " << newPo->getSimVal().size() << endl;
    // assert((oldPo->getSimVal().size() == newPo->getSimVal().size() == oldGate->getSimVal().size() == newGate->getSimVal().size()));
    unsigned numGood = 0;
    unsigned numSame = 0;

    

    for(unsigned i=0; i<getNumSim(); ++i) {
        size_t oldPoVal = oldPo->getSimVal(i);
        size_t newPoVal = newPo->getSimVal(i);
        size_t goodVectorMask = ~(oldPoVal ^ newPoVal);

        size_t oldGateVal = oldGate->getSimVal(i);
        size_t newGateVal = newGate->getSimVal(i);

        size_t xnorVal = ~(oldGateVal ^ newGateVal);
        
        // use the mask to filter out the bits outside of good vector
        xnorVal &= goodVectorMask;

        numGood += BitCount(goodVectorMask);
        numSame += BitCount(xnorVal);
    }

    return numGood ? ((double)numSame / numGood) : -1;
}

#endif