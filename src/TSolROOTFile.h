#ifndef __TSOLROOTFILE_H
#define __TSOLROOTFILE_H

#include "TFile.h"
#include "TTree.h"
#include "TSolEVIOFile.h"
#include "TSolGEMData.h"
#include <vector>
#include "TROOT.h"
#include "TMath.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "types.h"
#include "TSolDBManager.h"

using namespace std;


class TSolROOTFile {
 public:

  TSolROOTFile();
  TSolROOTFile( const char *name, int source );
  virtual ~TSolROOTFile();
  
  void  SetFilename( const char *name );
  void  Clear();
  Int_t Open();
  void  SetBranchAddress();
  Int_t Close();
  
  const char* GetFileName() const { return fFilename; }
  Int_t GetSource() const { return fSource; }
  Int_t ReadNextEvent();
  void  ExtractDetIDs();
  void  BuildData();
  void  BuildECData();
  void  BuildGenerated();
  void  ClearVectors();
  void  DeleteVectors();
  UInt_t CalcHitChamber(int detid);
  //    void  AddDatum(int crate, int slot, int chan, double data );
  
  UInt_t GetNData() const { return fHitData.size(); }
  UInt_t GetNGen() const { return fGenData.size(); }
  UInt_t GetNEC()  const { return fECData.size(); }

  UInt_t GetEvNum() const { return fEvNum; }
  UInt_t GetMaxEvNum() const { return fMaxEvNum; }

  hitdata *GetHitData(Int_t i) const { return fHitData[i]; }
  gendata *GetGenData(Int_t i) const { return fGenData[i]; }
  ECdata  *GetECData(Int_t i)  const { return fECData[i];  }
  TSolGEMData *GetGEMData();
  void GetGEMData(TSolGEMData* gd);

  void AddSignalInfo(Int_t pid, Int_t tid) 
  { fSignalInfo.push_back(SignalInfo(pid, tid)); }

  vector<SignalInfo> * GetSignalInfo() { return &fSignalInfo; }
  Int_t GetNSignal() const { return fSignalInfo.size(); }
  Int_t GetSigECBit(Int_t i) const { return fSignalInfo.at(i).fillBitsEC; }
  Double_t GetSigECEDep(Int_t i) const { return fSignalInfo.at(i).ECEDep; }
  Double_t GetSigMomentum(Int_t i) const { return fSignalInfo.at(i).momentum; }
  Double_t GetSigR(Int_t i) const { return fSignalInfo.at(i).R; }

  private:
  char  fFilename[255];
  TFile *fChan;
  TTree *tree_generated;
  TTree *tree_solid_gem;
  TTree *tree_flux;
  TTree *tree_header;
  Int_t fSource;   // User-defined source ID (e.g. MC run number)
  
  //generated branch
  vector<int> *gen_pid;
  vector<double> *gen_px;
  vector<double> *gen_py;
  vector<double> *gen_pz;
  vector<double> *gen_vx;
  vector<double> *gen_vy;
  vector<double> *gen_vz;
  
  //solid_gem branch
  vector<double> *solid_gem_id;
  vector<double> *solid_gem_hitn;
  vector<double> *solid_gem_pid;
  vector<double> *solid_gem_trid;
  vector<double> *solid_gem_x;
  vector<double> *solid_gem_y;
  vector<double> *solid_gem_z;
  vector<double> *solid_gem_lxin;
  vector<double> *solid_gem_lyin;
  vector<double> *solid_gem_lzin;
  vector<double> *solid_gem_tin;
  vector<double> *solid_gem_lxout;
  vector<double> *solid_gem_lyout;
  vector<double> *solid_gem_lzout;
  vector<double> *solid_gem_tout;
  vector<double> *solid_gem_px;
  vector<double> *solid_gem_py;
  vector<double> *solid_gem_pz;
  vector<double> *solid_gem_vx;
  vector<double> *solid_gem_vy;
  vector<double> *solid_gem_vz;
  vector<double> *solid_gem_ETot;
  vector<double> *solid_gem_trE;
  vector<double> *solid_gem_weight;

  //flux branch, used to record calorimeter hit position and energy
  //deposition, used later for tracking
  vector<double> *flux_id;
  vector<double> *flux_pid;
  vector<double> *flux_tid;
  vector<double> *flux_trackE;
  vector<double> *flux_avg_x;
  vector<double> *flux_avg_y;
  vector<double> *flux_avg_z;
  
  //header branch
  vector<double> *header_var8;
  

  vector<hitdata *>   fHitData;
  vector<gendata *>   fGenData;
  vector<ECdata *>    fECData;
  
  unsigned int fEvNum;
  unsigned int fMaxEvNum;

  vector<SignalInfo> fSignalInfo;  
  
  TSolDBManager* fManager;
};

#endif //__TSOLROOTFILE_H
