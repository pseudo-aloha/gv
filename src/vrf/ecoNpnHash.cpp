#ifndef ECO_NPN_HASH_CPP
#define ECO_NPN_HASH_CPP

#include "ecoMgr.h"
#include "cmath"
#include "algorithm"
#include <sstream>
#include <iomanip>

namespace gv {
namespace eco {

extern "C"{
void Abc_NtkShow( Abc_Ntk_t * pNtk, int fGateNames, int fSeq, int fUseReverse, int fKeepDot );
void Abc_TruthNpnTest( char * pFileName, int NpnType, int nVarNum, int fDumpRes, int fBinary, int fVerbose );
void Abc_NtkCecFraig( Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2, int nSeconds, int fVerbose );
void Abc_NtkCecSat( Abc_Ntk_t * pNtk1, Abc_Ntk_t * pNtk2, int nConfLimit, int nInsLimit );
Abc_Ntk_t * Abc_NtkMulti( Abc_Ntk_t * pNtk, int nThresh, int nFaninMax, int fCnf, int fMulti, int fSimple, int fFactor );
}

bool getIthBit(const size_t num, int i) {
    return ((num >> i) & 1);
}

void printBits(size_t tt) {
    for (int i = sizeof(size_t) * 8 - 1; i >= 0; --i) {
        cout << getIthBit(tt, i);
    }
    cout << endl;
}

// write the hexidicimal truth table of a k-feasible cut function to a file
void
writeHexTT(size_t i, int k) {
//   extern void printBits(size_t  num);
  ofstream f;
  f.open("./.hexTT.txt", ios::out);
  k=pow(2, k) / 4 - 1;
  
  while(k>=0) {
    f << std::hex << (i>>(4*k)&(0xf));
    // cout << std::hex << (i>>(4*k)&(0xf)) << " ";
    // printBits((size_t)(i>>(4*k)&(0xf)));
    k--;
  }

  f.close();
}

void combination(vector<vector<int>>& npnMatchWays,vector<int>& arr, int start, int end, int idx) {
  for(int i=start; i<end; ++i) {
      arr[i] += 1;
      combination(npnMatchWays, arr, i+1, end, idx+1);
      arr[i] -= 1;
  }
  npnMatchWays.push_back(arr);
}

int factorial(int n) {
    if(n<=0) return 0;
    int ret=1;
    while(n) {
        ret*=n;
        n--;
    }
    return ret;
}

// get the binary truth table of a k-feasible cut function to a file
string
getBinTT(size_t i, int k) {
  string ret;
  k=pow(2, k) - 1;
  while(k>=0) {
    ret.push_back(((i>>k) & (size_t)1) + '0');
    k--;
  }
  return ret;
}

// generate the input match transform from given the NP match
int
getInputTransform(int idx, const vector<int>& match) {
  int ret=0, n=match.size();
  // match[0] is reserved for PO negation
  for(int i=1; i<n; ++i) {
    ret += ((getIthBit(idx, match[i]/2-1) ^ (match[i]%2)) << (i-1));
  }
  return ret;
}

string
toHex(string s, int k) {
  int i = 0;
  for(int j=s.length()-1; j>=0; --j) {
    i+=((s[j]-'0')<<j);
  }
  string ret;
  stringstream ss;
  ss << std::hex << i;
  ss >> ret;

  string padding;
  for(i=ret.size(); i<pow(2, k-2); ++i) 
    padding.push_back('0');
  return padding+ret;
}
  
void
EcoNPNHash::computeNpnMatchWays(int k) {
    _npnHashWays.clear();
    vector<int> IOMap;
    for(int i=0; i<=k; i++)
        IOMap.push_back(2*i);
    for(int i=0; i<factorial(k); ++i) {
        combination(_npnHashWays, IOMap, 0, k+1, 0); // k+1 for output negation
        next_permutation(IOMap.begin()+1, IOMap.end());
    }
    // cout << "npn match size " << _npnHashWays.size() << endl;
}

void
EcoNPNHash::npnHash(const vector<vector<int>>& npnMatchWays, int i, int k) {
  ifstream f("./.hexTT_out.txt", ios::in);
  string npnClass, npnClassTT, funcTT;
  size_t npnClassInt;
  int offset = 0;
  for(int j=2; j<k; ++j)
    offset+=pow(2, pow(2, j));

  f >> npnClass;
  f.close();
  funcTT = getBinTT(i,k);
  stringstream ss;
  ss << std::hex << npnClass;
  ss >> npnClassInt;
  npnClassTT = getBinTT(npnClassInt, k);
  

  // test which npn match way is for this boolean function
  int n=pow(2, k);
  for(const auto& m : npnMatchWays) {
    bool outputNeg=m[0], flag=true;
    
    for(int j=0; j<n; ++j) {
      int idx = getInputTransform(j, m);
      if(((funcTT[idx]-'0') ^ outputNeg) != (npnClassTT[j]-'0')) {
        // cout << "fail" << endl;
        // cout << npnClassInt << " " << "npn " << getBinTT(npnClassInt, k) << " tt " << funcTT << endl;
        // cout << j << " " << idx << " " << funcTT[idx] << " " << npnClassTT[j] << endl;
        flag=false;
        break;
      } 
    }
    // found the mathcing way
    if(flag) {
      // cout << "match way : " << endl;
      // for(auto obj : m) {
      //   cout << obj << " " ;
      // }
      // cout << endl;
      // for(int j=0; j<n; ++j) {
      //   int idx = getInputTransform(j, m);
      //   // cout << npnClassTT << " " << funcTT << endl;
      //   cout << j << " " << idx << " " << (npnClassTT[j]-'0')  << " " << ((funcTT[idx]-'0')) << " ftt= " << npnClassTT << endl;
      // }
      // cout << endl;
    //   cout << "funcTT : " << funcTT << endl;
    //   cout << "npnClassTT : " << npnClassTT << endl;
      _npnHashTable[i+offset].push_back({toHex(npnClassTT, k), m});
    //   cout << "matched" << endl;
    //   for(auto obj : m) {
    //     cout << obj << " " ;
    //   }
      // cout << endl;
      // break;
    }
  }
  // cout << npnClassInt << " " << "npn " << getBinTT(npnClassInt, k) << " tt " << getBinTT(i,k) << endl;
}


void
EcoNPNHash::computeNpnHash() {
    int numElement = 0;
    for(int i=_cutSizeFrom; i<=_cutSizeTo; ++i)
    numElement += pow(2, pow(2, i));
    // used to record the tt match to which npn-class and the I/O match transformation (include permutation and negation)
    _npnHashTable = vector<vector<pair<string, vector<int>>>>(numElement);

    // check if the precompute of NPN hash table is done
    ifstream fNPNHash("./.BoolTTHash.txt");

    // if the precompute is not yet available, compute it and store it in the file
    if(!fNPNHash.is_open()) {
        for(int k = _cutSizeFrom; k <= _cutSizeTo; ++k) {
            computeNpnMatchWays(k);
            // compute k feasible NPN hash
            for(int i=0; i<pow(2, pow(2, k)); ++i) {
                writeHexTT(i, k);
                Abc_TruthNpnTest("./.hexTT.txt", 9, -1, 1, 0, 1 );
                npnHash(_npnHashWays, i, k);
            }
        }

        // write out the hash
        // cout << "match table size " << _npnHashTable.size() << endl;
        ofstream fo;
        fo.open("./.BoolTTHash.txt", ios::out);
        fo << "npn class" << setw(10) << "match" << endl;
        for(const auto& it : _npnHashTable) {
        for(const auto& it2 : it) {
            auto[npnClass, match] = it2;
            fo << npnClass << setw(3);
            for(const auto& m : match) {
                fo << to_string(m) << " ";
            }
            fo << setw(6);
        }
        fo << "#" << endl;
        }
        fo.close();
    }
    // if the precompute is already available, just read it from file
    else {
        int lineIdx = 0;
        string buf;
        fNPNHash >> buf >> buf >> buf; // pop out headers
        while(!fNPNHash.eof()) {
        string npnClass;
        vector<int> match;
        match.reserve(5);
        fNPNHash >> npnClass; // npn class
        if(npnClass == "#") {
            lineIdx++;
            if(lineIdx == numElement)
            break;
            fNPNHash >> npnClass;
        }
        if(!npnClass.length())
            break;
        int k=log2(npnClass.length()*4);
        for(int i=0; i<k+1; ++i) { // add one for po
            fNPNHash >> buf;
            match.push_back(stoi(buf));
        }
        _npnHashTable[lineIdx].push_back({npnClass, match});
        }
        fNPNHash.close();
    }
}

// encode function (used to save memory)
// <unused bits, output inv bit, input1 pos bits (3 bits), input1 inv bit, input2 pos bits (3 bits), input2 inv bit, ...>
size_t
EcoNPNHash::encodeMatch2SizeT(int outputMatch, vector<int>& inputMatch) {
  size_t ret = 0;
  
  // add the output info
  if(outputMatch)
    ret += 1;

  for(const auto& m : inputMatch) {
    ret <<= 3;
    ret += m / 2;

    ret <<= 1;
    if(m % 2 == 1)
      ret += 1;
  }

  return ret;
}

pair<int, vector<int>>
EcoNPNHash::decodeEncodedSizeTMatch(size_t encode, unsigned cutSize) {
  vector<int> inputMatch(cutSize);
  int outputMatch;

  size_t posMask = 0b111;
  size_t invMask = 0b1;

  for(int i=cutSize-1; i>=0; i--) {
    int m = 0;
    if(encode & invMask == 1)
      m = 1;
    encode >>= 1;

    int pos = (encode & posMask);
    m += pos * 2;
    encode >>= 3;

    inputMatch[i] = m;
  }

  outputMatch = (encode & invMask);


  return {outputMatch, inputMatch};
}

// return the full list of encoded valid matching
// this function only support 4-feasible cut
pair<string, vector<vector<int>>>
EcoNPNHash::getNPNHashFull(size_t cutTT, unsigned simSize) {
   assert(simSize <= 4); // the cut size should be less or equal to 4
   vector<vector<int>> retMatch;
   string NPNClass;

   if(simSize <= 1) {
    if(cutTT)
      return {"0", {{0,2}, {1,3}}};
    else
      return {"0", {{0,3}, {1,2}}};
  }

   unsigned offset = 0;
  for(unsigned i=_cutSizeFrom; i<simSize; ++i)
    offset+=pow(2, pow(2, i));

    unsigned idx = offset + cutTT;
    assert(idx < _npnHashTable.size());

    auto matches = _npnHashTable.at(idx);

    NPNClass = matches.at(0).first;

  for(const auto&[_, m] : matches)
    retMatch.push_back(m);

    return {NPNClass, retMatch};
}


pair<string, vector<int>>
EcoNPNHash::getNPNHash(size_t cutTT, unsigned cutSize) {
  if(cutSize <= 1) {
    if(cutTT)
      return {"0", {0, 2}};
    else
      return {"0", {0, 3}};
  }
  else if(cutSize <= 4) {
    // compute the offset for the smaller cuts
    unsigned offset = 0;
    for(unsigned i=_cutSizeFrom; i<cutSize; ++i)
      offset+=pow(2, pow(2, i));

    unsigned idx = offset + cutTT;
    assert(idx < _npnHashTable.size());
    return _npnHashTable[idx][0];
  }

  else {
    assert(cutSize <= 6); // currently support 6-feasible cut, but should can be greater than this
    // cout << "don don  " << cutSize << endl;
    // printBits(cutTT);
    writeHexTT(cutTT, cutSize);
    Abc_TruthNpnTest("./.hexTT.txt", 9, -1, 1, 0, 0 );

    ifstream f("./.hexTT_out.txt", ios::in);
    string npnClassAbc, npnClassTT, funcTT;

    f >> npnClassAbc;
    f.close();

    // cout << "npn class : " << npnClassAbc << endl;
    size_t npnClassInt;
    // funcTT = getBinTT(i,k);
    stringstream ss;
    ss << std::hex << npnClassAbc;
    ss >> npnClassInt;
    npnClassTT = getBinTT(npnClassInt, cutSize);

    string myNpnClass = toHex(npnClassTT, cutSize); // append leading zero to avoid name conflix (though k = 2 and 3 still have same length in hex)
    


    return {myNpnClass, {}};
  }

}

// given a cut, compute its NPN class and matching
pair<string, vector<int>>
EcoMgr::getNPNHash(gv::cir::EcoCut* cut) {
  unsigned cutSize = cut->getCutSize(); // get the cut size
  size_t cutTT;
  // get the truth table of the cut
  if(cut->getRoot()->isOld())
    cutTT = _oldNtk->computeCutTT(cut);
  else
    cutTT = _newNtk->computeCutTT(cut);

  return _pNpnHash->getNPNHash(cutTT, cutSize);
}

// given a cut, compute its NPN class and matching
pair<string, vector<vector<int>>>
EcoMgr::getNPNHashFull(gv::cir::EcoCut* cut) {
  unsigned cutSize = cut->getCutSize(); // get the cut size
  size_t cutTT;
  // get the truth table of the cut
  if(cut->getRoot()->isOld())
    cutTT = _oldNtk->computeCutTT(cut);
  else
    cutTT = _newNtk->computeCutTT(cut);

  return _pNpnHash->getNPNHashFull(cutTT, cutSize);
}

} // end of namespace gv::eco
} // end of namespce gv


#endif