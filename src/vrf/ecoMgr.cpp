#ifndef ECO_MGR_CPP
#define ECO_MGR_CPP

#include "ecoMgr.h"

namespace gv {
namespace eco {
// top function to do ECO
void
EcoMgr::doEco(const string& oldDesignName, const string& newDesignName) {
  // read designs
  readDesigns(oldDesignName, newDesignName);
  doFraig();
}

// read input designs
void
EcoMgr::readDesigns(const string& oldDesignName, const string& newDesignName) {
  _oldNtk->readNtkFile(oldDesignName);
  cout << "------------------------------------" << endl;
  _newNtk->readNtkFile(newDesignName);
}

// perform abc fraig on the circuits
void
EcoMgr::doFraig() {
  Abc_Ntk_t* pNtkOld = _oldNtk->getAbcNtk();
  Abc_Ntk_t* pNtkNew = _newNtk->getAbcNtk();
  Abc_Ntk_t* pNtkMiter = Abc_NtkMiter( pNtkOld, pNtkNew, 0, 0, 0, 0 );

  Abc_Obj_t * pNode;
  int i;
  Abc_NtkForEachObj( pNtkOld, pNode, i ) {
    if(!pNode->pCopy) continue;
    cout << "obj " << Abc_ObjName(pNode) << " " << Abc_ObjName(pNode->pCopy) << endl;
  }
}

// end of namespace gv::eco
}}
#endif