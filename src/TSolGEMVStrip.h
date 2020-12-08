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
    VStrips (Short_t nsample) : fIst(0), fIss(0), fCharge(0) {fADC.assign (nsample, 0);}
    UInt_t              fIst;     // index of the strip
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

  void AddStripAt(Short_t ist, Short_t n) {fVStrips.at(n)->fIst=ist;};
  void AddSubstripAt(Short_t iss, Short_t n) {fVStrips.at(n)->fIss=iss;};
  void AddChargeAt(Double_t val, Short_t n) {fVStrips.at(n)->fCharge=val;}
  void AddSampleAt(Int_t val, Short_t sample, Short_t n) {fVStrips.at(n)->fADC.at(sample)=val;}

  void SetTime(Double_t val) {fTime = val;};
  void SetHitCharge(Float_t val) { fHitCharge = val; }; // total charge of the avalanche
  void SetSize(Short_t n) { fSize = n; };

  Short_t GetSize() const { return fSize; };
  Short_t GetStrip(Short_t n) const { return fVStrips.at(n)->fIst; };
  Short_t GetSubstrip(Short_t n) const { return fVStrips.at(n)->fIss; };
  Int_t GetADC(Short_t n, Short_t sample) const { return fVStrips.at(n)->fADC.at(sample); }
  Double_t GetCharge(Short_t n) const { return fVStrips.at(n)->fCharge; }
  Float_t GetTime() const { return fTime; };
  Double_t GetHitCharge() const { return fHitCharge; };

  void Print() const;
};

#endif
