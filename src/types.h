#ifndef __TYPES_
#define __TYPES_

enum GEMDir_t { kGEMX, kGEMY, kGEMR, kGEMPhi };

struct SignalInfo{
    Int_t pid;
    Int_t tid;
    Int_t fillBitsGEM;
    Int_t fillBitsEC;
    Double_t ECEDep;
    Double_t momentum;
    Double_t R;
    SignalInfo() {}
    SignalInfo(Int_t apid, Int_t atid):pid(apid), tid(atid),
    fillBitsGEM(0), fillBitsEC(0), ECEDep(0.) {}
    ~SignalInfo() {}
  };


#endif//__TYPES_
