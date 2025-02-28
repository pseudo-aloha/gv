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
}

// read input designs
void
EcoMgr::readDesigns(const string& oldDesignName, const string& newDesignName) {
  _oldNtk->readNtkFile(oldDesignName);
  cout << "------------------------------------" << endl;
  _newNtk->readNtkFile(newDesignName);
}

// end of namespace gv::eco
}}
#endif