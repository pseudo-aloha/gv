#ifndef ECO_CIR_CPP
#define ECO_CIR_CPP

#include "cirMgr.h"

#include <ctype.h>

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <iostream>

#include "cirDef.h"
#include "cirGate.h"
#include "util.h"
#include "yosysMgr.h"
#include "ecoNtk.h"

namespace gv {
namespace cir {

bool
newFaninComp(Abc_Obj_t *pObj, int inIdx) {
    assert(inIdx == 0 || inIdx == 1);
    bool ret = false;
    if(inIdx == 0) {
        ret = Abc_ObjFaninC0(pObj);
        if(Abc_AigNodeIsConst(Abc_ObjFanin0(pObj))) {
            ret ^= true;
        }
    }
    else {
        ret = Abc_ObjFaninC1(pObj);
        if(Abc_AigNodeIsConst(Abc_ObjFanin1(pObj))) {
            ret ^= true;
        }
    }
    return ret;
}

void 
EcoCir::readCirFromAbcNtk(Abc_Ntk_t* pNtk) {
    _ecoCirV = new gv::cir::CirMgr();
    // TODO : Convert abc ntk to gv aig ntk
    CirGateV gateV;
    // Abc_Ntk_t* pNtk = NULL;            // the gia pointer of abc
    Abc_Obj_t *pObj, *pObjRi, *pObjRo; // the obj element of gia
    unsigned iPi = 0, iPo = 0, iRi = 0, iRo = 0;
    int i;

    // initialize the size of the containers
    _ecoCirV->initCir(Abc_NtkPiNum(pNtk), Abc_NtkPoNum(pNtk), Abc_NtkLatchNum(pNtk), Abc_NtkObjNumMax(pNtk));

    // set the const node
    _ecoCirV->_totGateList[0] = _ecoCirV->_const0;
    // traverse the obj's in topological order
    Abc_NtkForEachObj( pNtk, pObj, i) {
        char *name = new char[strlen(Abc_ObjName(pObj))+1]; strcpy(name, Abc_ObjName(pObj));
        if(Abc_ObjIsPi(pObj)) {
            // cout<< "pi " << Abc_ObjId(pObj) << endl;
            
            CirPiGate* gate = new CirPiGate(Abc_ObjId(pObj), 0);
            
            _ecoCirV->_piList[iPi++] = gate;
            _ecoCirV->_totGateList[Abc_ObjId(pObj)] = gate;
            gate->setName(name);
        }
        else if(Abc_ObjIsPo(pObj)) {
        //  
        // cout<< "po " << Abc_ObjId(pObj) << endl;
        }
        else if(Abc_ObjIsNode(pObj)) {
            CirAigGate *gate = new CirAigGate(Abc_ObjId(pObj), 0);
            _ecoCirV->_totGateList[Abc_ObjId(pObj)] = gate;
            gate->setIn0(_ecoCirV->getGate(Abc_ObjId(Abc_ObjFanin0(pObj))), newFaninComp(pObj, 0));
            gate->setIn1(_ecoCirV->getGate(Abc_ObjId(Abc_ObjFanin1(pObj))), newFaninComp(pObj, 1));
            gate->setName(name);
        }
        // TODO : handle latches

        // else if(Abc_ObjIsBo(pObj)) {
        //     CirRoGate* gate = new CirRoGate(Abc_ObjId(pObj), 0);
        //     _roList[iRo++] = gate;
        //     _totGateList[Abc_ObjId(pObj)] = gate;
        // }
        // else if(Abc_ObjIsBi(pObj)) {
        //     CirRiGate *gate = new CirRiGate(Abc_ObjId(pObj), 0, Abc_ObjId(Abc_ObjFanin0(pObj)));
        //     gate->setIn0(getGate(Abc_ObjId(Abc_ObjFanin0(pObj))), newFaninComp(pObj, 0));
        //     _riList[iPo++] = gate;
        //     _totGateList[Abc_ObjId(pObj)] = gate;
        // }
        else if(Abc_AigNodeIsConst(pObj)) {
            // cout << "I am const1 " << Abc_ObjId(pObj) <<  endl;
        }
        else {
            // cout << "not defined gate type" << endl;
        }
    }
    // handle the po's
    Abc_NtkForEachPo(pNtk, pObj, i) {
        char *name = new char[strlen(Abc_ObjName(pObj))+1]; strcpy(name, Abc_ObjName(pObj));
        CirPoGate *gate = new CirPoGate(Abc_ObjId(pObj), 0, Abc_ObjId(Abc_ObjFanin0(pObj)));
        gate->setIn0(_ecoCirV->getGate(Abc_ObjId(Abc_ObjFanin0(pObj))), newFaninComp(pObj, 0));
        _ecoCirV->_poList[iPo++] = gate;
        _ecoCirV->_totGateList[Abc_ObjId( pObj)] = gate;
        gate->setName(name);
    }
    
    // add the fanout information
    Abc_NtkForEachObj( pNtk, pObj, i) {
        CirGate* gate = _ecoCirV->getGate(Abc_ObjId(pObj));
        if(gate->getType() == PI_GATE || gate->getType() == AIG_GATE) {
            for(int j=0; j<Abc_ObjFanoutNum(pObj); ++j)
                _fanoutMap[gate].push_back(_ecoCirV->getGate(Abc_ObjId(Abc_ObjFanout(pObj, j))));
        }
    }

    _ecoCirV->genDfsList();
    // _ecoCirV->printNetlist();
    // checkFloatList();
    // checkUnusedList();

}

// end of namespace gv::cir
}}
#endif