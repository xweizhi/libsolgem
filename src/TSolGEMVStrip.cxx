#include <iostream>

#include "TSolGEMVStrip.h"

using namespace std;

TSolGEMVStrip::TSolGEMVStrip(Short_t n, // number of strips in hit 
			   Short_t nsample) {
  fSize = n;
  fNsample = nsample;
  fTime = -1.;
  fHitCharge = 0;
  fVStrips.resize (n);
  for (std::vector<VStrips*>::iterator i = fVStrips.begin(); i != fVStrips.end(); ++i)
    *i = new VStrips (nsample);
};

TSolGEMVStrip::~TSolGEMVStrip() {
  for (std::vector<VStrips*>::iterator i = fVStrips.begin(); i != fVStrips.end(); ++i)
    delete *i;
};

void 
TSolGEMVStrip::Print() {
  Int_t i,k;
  
  cerr << "Virtual strips sampled starting at " << GetTime() << endl;
  
  for (i=0;i<GetSize();i++) {
    
    cerr << i << ") " << GetStrip(i) << " / " << GetSubstrip(i) << " : ";
    cerr << GetCharge(i) << " = ";
    
    for (k=0;k<fNsample;k++) {
      cerr << GetADC(i,k) << " ";
    }
    
    cerr << endl;
  }
  
};
