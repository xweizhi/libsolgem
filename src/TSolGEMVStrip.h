// Define the virtual strip class, used by digi
// This describes the response of multiple strips in a plane to a 
// single hit.
// N.B. The strip is virtual, not the class!
//

#ifndef __TSolGEMVStrip__
#define __TSolGEMVStrip__

#include <TArrayD.h>
#include <TArrayS.h>

#include <vector>

class TSolGEMVStrip {

 private:
  // working variable of digi method (virtual strips)

  struct VStrips
  {
    VStrips (Short_t nsample) : fIdx(0), fIss(0), fCharge(0) {fADC.assign (nsample, 0);}
    UInt_t              fIdx;     // index of the strip
    UInt_t              fIss;     // index of the substrip
    std::vector<UInt_t> fADC;     // value of the ADC sample 
    Double_t            fCharge;  // total charge in substrips (sampled)
  };
  std::vector<VStrips*> fVStrips;

  Float_t fTime; // time of first sampling
  Int_t fNsample; // number of samples

  Float_t fHitCharge;

  Short_t fSize; // effective size, maybe different from size of the array (some elements maybe not used)

 public:

  TSolGEMVStrip(Short_t n = 10, // number of strips in hit
 	       Short_t nsample = 10);
  ~TSolGEMVStrip();

  void AddStripAt(Short_t idx, Short_t n) {fVStrips.at(n)->fIdx=idx;};
  void AddSubstripAt(Short_t iss, Short_t n) {fVStrips.at(n)->fIss=iss;};
  void AddChargeAt(Double_t val, Short_t n, Short_t iss) {fVStrips.at(n)->fCharge=val;}
  void AddSampleAt(Short_t val, Short_t sample, Short_t n, Short_t iss) 
  {fVStrips.at(n)->fADC.at(sample)=val;}
  // ***** Just for testing, get rid of when done
  void AddChargeAt(Double_t val, Short_t n) {AddChargeAt (val, n, 0);}
  void AddSampleAt(Short_t val, Short_t sample, Short_t n) 
  {AddSampleAt (val, sample, n, 0);}
  // ***** end just for testing

  void SetTime(Double_t val) {fTime = val;};
  void SetHitCharge(Float_t val) { fHitCharge = val; }; // total charge of the avalanche
  void SetSize(Short_t n) { fSize = n; };

  Short_t GetSize() const { return fSize; };
  Short_t GetIdx(Short_t n) const { return fVStrips.at(n)->fIdx; };
  Short_t GetIss(Short_t n) const { return fVStrips.at(n)->fIss; };
  Short_t GetADC(Short_t n, Short_t iss, Short_t sample) const { return fVStrips.at(n)->fADC.at(sample); }
  Double_t GetCharge(Short_t n, Short_t iss) const { return fVStrips.at(n)->fCharge; }
  // ***** Just for testing, get rid of when done
  Short_t GetADC(Short_t n, Short_t sample) const { return GetADC (n, 0, sample); }
  Double_t GetCharge(Short_t n) const { return GetCharge (n, 0);}
  // ***** end just for testing
  Float_t GetTime() const { return fTime; };
  Double_t GetHitCharge() const { return fHitCharge; };

  void Print();

 //  // TODO: rewrite this using STL containers ...
 //  TArrayS *fIdx;  // index of the strip
 //  TArrayS *fADC; // value of the ADC sample
 //  TArrayD *fCharge; // total charge in strip (sampled)

 //  Float_t fTime; // time of first sampling
 //  Int_t fNsample; // number of samples

 //  Float_t fHitCharge;

 //  Short_t fSize; // effective size, maybe different from size of the array (some elements maybe not used)

 // public:

 //  TSolGEMVStrip(Short_t n = 10, // number of strips in hit
 // 	       Short_t nsample = 10);
 //  ~TSolGEMVStrip();

 //  void AddStripAt(Short_t idx, Short_t n) {fIdx->AddAt(idx,n);};
 //  void AddChargeAt(Double_t val, Short_t n) {fCharge->AddAt(val,n);};
 //  void AddSampleAt(Short_t val, Short_t sample, Short_t n) {fADC->AddAt(val,n*fNsample+sample);};
 //  void SetTime(Double_t val) {fTime = val;};
 //  void SetHitCharge(Float_t val) { fHitCharge = val; }; // total charge of the avalanche
 //  void SetSize(Short_t n) { fSize = n; };

 //  Short_t GetSize() const { return fSize; };
 //  Short_t GetIdx(Short_t n) const { return fIdx->At(n); };
 //  Short_t GetADC(Short_t n, Short_t sample) const { return fADC->At(n*fNsample+sample); };
 //  Double_t GetCharge(Short_t n) const { return fCharge->At(n); };
 //  Float_t GetTime() const { return fTime; };
 //  Double_t GetHitCharge() const { return fHitCharge; };

 //  void Print();
};

#endif
