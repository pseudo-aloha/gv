#ifndef ECO_MGR_H
#define ECO_MGR_H

#include "cirMgr.h"
#include "ecoNtk.h"
// #include "sat.h"
#include <unordered_set>
#include <map>
#include <queue>
#include <set>
namespace gv {
namespace eco {
// forward decleration of classes
class EcoMgr;
class EcoCir;
// class EcoNtk;
// class CirMgr;


// the main class for Andrew's ECO Approach
class EcoMgr {
public:
  // constructor
  EcoMgr () { _oldNtk = new gv::cir::EcoNtk;
              _newNtk = new gv::cir::EcoNtk; }
  ~EcoMgr () { delete _oldNtk; delete _newNtk; }
  void doEco(const string& oldDesignName, const string& newDesignName);
  void readDesigns(const string& oldDesignName, const string& newDesignName);
  void doFraig();
  
private:
  gv::cir::EcoNtk* _oldNtk;
  gv::cir::EcoNtk* _newNtk;
};

}}

#endif