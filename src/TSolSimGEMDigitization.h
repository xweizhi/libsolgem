// Class handling digitization of GEMs

#ifndef __TSolSimGEMDigitization__
#define __TSolSimGEMDigitization__

#include "TRandom3.h"
#include "TVector3.h"
#include "TArrayI.h"
#include "TH1F.h"
#include "THaAnalysisObject.h"

#include <vector>
using namespace std;
class TFile;
class TTree;

class TSolGEMData;
class TSolGEMPlane;
class TSolGEMVStrip;
class TSolSpec;
class TSolSimEvent;
class TSolEVIOFile;
class TSolROOTFile;

// First an auxiliary class

// The whole strip plane; used to cumulate virtual strips charges
// and generate real strips

class TSolDigitizedPlane {
private:
  // ADC sampled value of strip array of each axis

  struct DigSubstrip
  {
    DigSubstrip (UShort_t nsample): fType(0), fTotADC(0), fMaxADC(0), fCharge(0), fTime(0) {fStripADC.assign (nsample, 0); fStripClusters.clear();}
    std::vector<Double_t> fStripADC;
    Short_t fType;
    Double_t   fTotADC;
    Double_t   fMaxADC;
    Float_t fCharge;
    Float_t fTime;
    std::vector<Short_t> fStripClusters; // Clusters associated with each strip
  };
  std::vector<DigSubstrip*> fDigPlane;

  TSolGEMPlane* fPlane;
  UShort_t fNSamples;
  UShort_t fNSubstrips;
  Int_t    fThreshold;
  UShort_t  fNOT;   // # substrips over threshold
  std::vector<Short_t> fOverThr;  // # list of substrips over threshold

  //used to simulate cross talk of APV25
  TRandom3 fRan;
public:
  // init and reset physics strips arrays
  TSolDigitizedPlane (const TSolGEMPlane* plane,
		      UShort_t nsample = 10,
		      Int_t    threshold = 0 );
  ~TSolDigitizedPlane();
  void Clear();

  // cumulate hits (strips signals)
  void Cumulate (const TSolGEMVStrip *vv, 
		 const TSolGEMPlane* plane,
		 Short_t type, Short_t clusterID );

  Short_t  GetType (Int_t n) const {return fDigPlane.at(n)->fType;}
  Double_t GetTotADC (Int_t n) const {return fDigPlane.at(n)->fTotADC;}
  Float_t  GetTime (Int_t n) const {return fDigPlane.at(n)->fTime;}
  Float_t  GetCharge (Int_t n) const {return fDigPlane.at(n)->fCharge;}
  Double_t GetADC (Int_t n, Int_t ks) const  {return fDigPlane.at(n)->fStripADC.at(ks);}
  UShort_t GetNSamples() const {return fNSamples;}
  UShort_t GetNSubstrips() const {return fNSubstrips;}

  UShort_t Threshold (Int_t thr);

  UShort_t GetNOverThr() const {return fNOT;}
  Short_t  GetIndexOverThr (Int_t n) const {return fOverThr[n];}

  const std::vector<Short_t>& GetStripClusters(UInt_t n) const { return fDigPlane.at(n)->fStripClusters; }

  const TSolGEMPlane* GetPlane() const { return fPlane; }
};

class TSolSimGEMDigitization: public THaAnalysisObject
{
 public:
  TSolSimGEMDigitization( const TSolSpec& spect,
			  const char* name = "testdigitization");
  virtual ~TSolSimGEMDigitization();

  void Initialize(const TSolSpec& spect);
  Int_t ReadDatabase (const TDatime& date);

  Int_t Digitize (const TSolGEMData& gdata, const TSolSpec& spect); // digitize event
  Int_t AdditiveDigitize (const TSolGEMData& gdata, const TSolSpec& spect); // add more digitized data, e.g. background
  void NoDigitize (const TSolGEMData& gdata, const TSolSpec& spect); // do not digitize event, just fill tree
  const TSolDigitizedPlane& GetDigitizedPlane (UInt_t ich, UInt_t ip) const {return *(fDP[ich][ip]);};
  void Print() const;
  void PrintCharges() const;
  void PrintSamples() const;

  Double_t GetGateWidth(){ return fGateWidth; }

  // Tree methods
  // To write a tree with digitization results:
  //   Call InitTree before main loop;
  //   Call SetTreeEvent in main loop (before or after Digitize)
  //   Call FillTree in main loop (after Digitize and SetTreeEvent)
  // Call WriteTree and CloseTree after main loop

  void InitTree (const TSolSpec& spect, const TString& ofile);
  void SetTreeEvent (const TSolGEMData& tsgd,
		     const TSolEVIOFile& f,
		     Int_t evnum = -1);
  void SetTreeEvent (const TSolGEMData& tsgd,
                     const TSolROOTFile& f,
                     Int_t evnum = -1); //overload for root file input
  Short_t SetTreeHit (const UInt_t ih,
		      const TSolSpec& spect,
		      TSolGEMVStrip* const *dh,
		      const TSolGEMData& tsgd,
		      Double_t t0 ); // called from Digitization
  void SetTreeStrips(const TSolSpec& spect); // called from Digitization
  void FillTree (const TSolSpec& spect);
  void WriteTree () const;
  void CloseTree () const;

  // Access to results
  Short_t GetType (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetType (n);}
  Int_t   GetTotADC (UInt_t ich, UInt_t ip, Int_t n) const {return fDP[ich][ip]->GetTotADC (n);}
  Float_t GetTime (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetTime (n);}
  Float_t GetCharge (UInt_t ich, UInt_t ip, UInt_t n) const {return fDP[ich][ip]->GetCharge (n);}
  Double_t   GetADC (UInt_t ich, UInt_t ip, Int_t n, Int_t ks) const {return fDP[ich][ip]->GetADC (n, ks);}
  UInt_t   GetNChambers() const {return fNChambers;};
  UInt_t   GetNPlanes (const UInt_t i) const {return fNPlanes[i];}
  UShort_t GetNSamples (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNSamples();}
  UShort_t GetNSubstrips (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNSubstrips();}
  UShort_t Threshold (UInt_t ich, UInt_t ip, Int_t thr) {return fDP[ich][ip]->Threshold (thr);}
  UShort_t GetNOverThr (UInt_t ich, UInt_t ip) const {return fDP[ich][ip]->GetNOverThr();}
  Short_t  GetIndexOverThr (UInt_t ich, UInt_t ip, Int_t n) const
  { return fDP[ich][ip]->GetIndexOverThr(n); }

  const std::vector<Short_t>& GetStripClusters(UInt_t ich, UInt_t ip, UInt_t n) const
  { return fDP[ich][ip]->GetStripClusters(n); }

  TSolSimEvent* GetEvent() const { return fEvent; }

  Bool_t IsMapSector() const { return fDoMapSector; }
  void SetMapSector( Bool_t b = true ) { fDoMapSector = b; }
  
  // APV cross talk parameters
  static Int_t    fDoCrossTalk;  //whether we want to do cross talk simulation
  static Int_t    fNCStripApart; // # of channels the induced signal is away from the mean signal
  static Double_t fCrossFactor;  //reduction factor for the induced signal
  static Double_t fCrossSigma;   //uncertainty of the reduction factor

 private:

  void IonModel (const TVector3& xi,
		 const TVector3& xo,
		 const Double_t elost );

  TSolGEMVStrip ** AvaModel (const Int_t ic,
			     const TSolSpec& spect,
			     const TVector3& xi,
			     const TVector3& xo,
			     const Double_t time_off);

  Double_t GetPedNoise(Double_t& phase, Double_t& amp, UInt_t& isample);
  void     FindFirstMaximum(Short_t bins, Double_t thres, Double_t& maxBin, Double_t& maxADC);
  Int_t    FindFirstBinAbove(Double_t thres, Double_t min, Double_t max);
  Int_t    FindLastBinAbove(Double_t thres, Double_t min, Double_t max);
  bool     IsMaximaExistBefore(Double_t thres, Int_t bin, Double_t time);

  // Gas parameters
  Double_t fGasWion;               // eV
  Double_t fGasDiffusion;          // mm2/s
  Double_t fGasDriftVelocity;      // mm/s
  Double_t fAvalancheFiducialBand; // number of sigma defining the band around the avalanche in readout plane
  Int_t    fAvalancheChargeStatistics;  // 0 Furry, 1 Gaussian
  Double_t fGainMean;
  Double_t fGain0;
  UInt_t   fMaxNIon;               //maximum amount of ion pairs allowed in the digitization
  
  Double_t fSNormNsigma;           //fSNormNsigma is an arbitrary multiplicative factor for the avalance radius.
  Int_t    fAvaModel;              //0 for Heavyside, 1 for Gaussian, 2 for Cauchy-Lorentz
  Double_t fAvaGain;
  Double_t fLateralUncertainty; // avalanche electrons can only pass through the holes of GEM foil
                                // which introduce additional uncertainty in the lateral direction
  // Electronics parameters
  vector<Double_t> fTriggerOffset;         // trigger offset (ns), incl latency & readout offset
  Double_t fTriggerJitter;         // trigger sigma jitter (ns)
  Double_t fAPVTimeJitter;         // time jitter associated with the APV internal clock
  Int_t    fEleSamplingPoints;
  Double_t fEleSamplingPeriod;     // ns
  Double_t fADCoffset;             // ADC offset
  Double_t fADCgain;               // ADC gain
  Int_t    fADCbits;               // ADC resolutions in bits
  Double_t fGateWidth;             // to be changed , ns - pulse shape width at ~1/10 max
  Double_t fGateWidthPost;         
                                
  Int_t    fChipMode;           //0 for APV25, 1 for SAMPA with 160ns shapping time, 2 for SAMPA
                                //with 80ns shapping time, 3 for VMM
  Double_t fVMMInteThreshold;   //integration threshold for VMM chip, only used if fChipMode is 3
  Int_t    fVMMInteTime;

  //parameter for GEM pedestal noise
  Double_t fPulseNoiseSigma;  // additional sigma term of the pedestal noise
  Double_t fPulseNoisePeriod; // period of the pedestal noise, assuming sinusoidal function
  Double_t fPulseNoiseAmpConst;  // constant term of the pedestal noise amplitude
  Double_t fPulseNoiseAmpSigma;  // sigma term of the pedestal noise amplitude

  // Pulse shaping parameters
  Double_t fPulseShapeTau0;   // [ns] GEM model; = 50. in SiD model
  Double_t fPulseShapeTau1;   // [ns] GEM model only; if negative assume SiD model
  
  // Geometry
  Double_t fRoutZ;            // z-distance hit entrance to readout plane [mm]
  Int_t    fUseTrackerFrame;       // tracker frame is used in the original version, but not so in my version
                                   // Weizhi Xiong
  Double_t fEntranceRef;           // z position of the copper layer right before the first GEM gas layer,
                             // relative to the center of the GEM chamber

  // Sector mapping
  Bool_t   fDoMapSector;
  Int_t    fSignalSector;

  //TDC measure for the trigger jitter
  Float_t  fJitterMeasure;
  
  //parameter for numerical integration
  UInt_t   fYIntegralStepsPerPitch;
  UInt_t   fXIntegralStepsPerPitch;

  TSolDigitizedPlane*** fDP; // 2D array of plane pointers indexed by chamber, plane #

  UInt_t fNChambers;  // # chambers
  UInt_t* fNPlanes;   // # planes in each chamber
  TRandom3 fTrnd;     // time randomizer
  UInt_t   fRNIon;    // number of ions
  struct IonPar_t {
    Double_t X;       // position of the point on the projection
    Double_t Y;
    Double_t Charge;  // Charge deposited by this ion
    Double_t SNorm;   // 3 x radius of ion diffusion area at readout
    Double_t R2;      // = SNorm^2 : radius of numerical integration area
    Double_t ggnorm;  // = Charge/R2/pi : charge per unit area
  };
  std::vector<IonPar_t> fRIon;
  Double_t fRSMax;
  Double_t fRTotalCharge;
  Double_t fRTime0;

  std::vector<Double_t> fSumA;   // Charge deposit in each bin, to be integrated
  std::vector<Double_t>  fDADC;

  TH1F* hADC;

  // Tree

  TFile* fOFile;          // Output ROOT file
  TTree* fOTree;          // Output tree
  TSolSimEvent* fEvent;   // Output event structure, written to tree

  Bool_t fFilledStrips;   // True if no data changed since last SetTreeStrips


  void MakePrefix() { THaAnalysisObject::MakePrefix(0); }
  void DeleteObjects();

  ClassDef (TSolSimGEMDigitization, 0)
};

#endif

