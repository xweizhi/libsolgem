#include <iostream>
#include <cassert>
#include <cmath>

#include "TSolGEMChamber.h"
#include "TSolGEMPlane.h"
#include "TSolPolygon.h"
#include "THaEvData.h"
#include "THaApparatus.h"
#include "TMath.h"
#include "ha_compiledata.h"

using namespace std;

TSolGEMChamber::TSolGEMChamber( const char *name, const char *desc )
  : THaDetector (name, desc)
{
  // For now at least we just hard wire two chambers
  fNPlanes = 2;
  fPlanes = new TSolGEMPlane*[fNPlanes];
  fWedge = new TSolWedge;
  fVector3.SetXYZ(0, 0, 0);
  return;
}

TSolGEMChamber::~TSolGEMChamber()
{
  for (UInt_t i = 0; i < fNPlanes; ++i)
    delete fPlanes[i];
  delete[] fPlanes;
  delete fWedge;
}


const char* TSolGEMChamber::GetDBFileName() const {
    THaApparatus *app = GetApparatus();
    if( app )
      return Form ("%s.", app->GetName());
    else
      return fPrefix;
}

Int_t
TSolGEMChamber::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  err = InitPlane (0, TString (GetName()) + "x", TString (GetTitle()) +" x");
  if( err != kOK ) return err;
  err = InitPlane (1, TString (GetName()) + "y", TString (GetTitle()) +" y");
  if( err != kOK ) return err;

  return kOK;
}

Int_t
TSolGEMChamber::ReadGeometry (FILE* file, const TDatime& date,
			      Bool_t required)
{

  Int_t err = THaDetector::ReadGeometry (file, date, required);
  if (err)
    return err;

  Double_t r0 = -999.0;
  Double_t r1 = -999.0;
  Double_t phi0 = -999.0;
  Double_t dphi = -999.0;
  Double_t z0 = -999.0;
  Double_t depth = -999.0;
  Int_t nOff = 0;
  fHVSectorOff.clear();
  const DBRequest request[] =
    {
      {"r0",               &r0,           kDouble, 0, 1},
      {"r1",               &r1,           kDouble, 0, 1},
      {"phi0",             &phi0,         kDouble, 0, 1},
      {"dphi",             &dphi,         kDouble, 0, 1},
      {"z0",               &z0,           kDouble, 0, 1},
      {"depth",            &depth,        kDouble, 0, 1},
      {"frame_width",      &fFrameWidth,  kDouble, 0, 1},
      {"n_HV_sector_off",  &nOff,         kInt,    0, 1},
      {0}
    };
  err = LoadDB (file, date, request, fPrefix);

  if (err)
    return err;
  
  //now read the information for HV sector to be turned off
  
  if (nOff > 0){
    
    for (Int_t i=0; i<nOff; i++){
        vector<Double_t>* HVbound = 0;
        try{
            HVbound = new vector<Double_t>;
            DBRequest HVRequest[] = {
            {"bound",       HVbound,          kDoubleV, 0, 1},
            { 0 }
            };

	    // HVbound is x1, y1, x2, y2, ... xn, yn n >= 3
	    //   vertex coordinates for a polygon
	    // OR for backward compatibility
	    // x1, x2, y1, y2 (note, NOT x1, y1, x2, y2)
	    //   vertex coordinates for 2 corners of rectangle
            
            ostringstream tmp_prefix(fPrefix, ios_base::ate);
            tmp_prefix<<"HV"<<i+1<<".";
            string HV_prefix = tmp_prefix.str();
            err = LoadDB (file, date, HVRequest, HV_prefix.c_str());
            if (err == kOK)
	      {
		std::vector<Double_t> polyx;
		std::vector<Double_t> polyy;
		if (HVbound->size() == 4)
		  {
		    // rectangle parallel to x/y axes
		    assert(HVbound->at(1) >= HVbound->at(0) && 
			   HVbound->at(3) >= HVbound->at(2));
		    polyx.push_back (HVbound->at(0)); polyy.push_back (HVbound->at(2));
		    polyx.push_back (HVbound->at(0)); polyy.push_back (HVbound->at(3));
		    polyx.push_back (HVbound->at(1)); polyy.push_back (HVbound->at(3));
		    polyx.push_back (HVbound->at(1)); polyy.push_back (HVbound->at(2));
		    fHVSectorOff.push_back (TSolPolygon (4, polyx, polyy));
		  }
		else
		  {
		    assert (HVbound->size() > 5 || HVbound->size()%2 == 0);
		    UInt_t ncorner = HVbound->size() / 2;
		    for (UInt_t icorner = 0; icorner < ncorner; ++icorner)
		      {
			polyx.push_back (HVbound->at(icorner*2)); 
			polyy.push_back (HVbound->at(icorner*2+1));
		      }
		    fHVSectorOff.push_back (TSolPolygon (ncorner, polyx, polyy));
		  }
	      }
            
            delete HVbound;
        }catch(...) {
            delete HVbound;
            fclose(file);
            throw;
        }
        
    }
  
  }
  
  
  // Database specifies angles in degrees, convert to radians
  Double_t torad = atan(1) / 45.0;
  phi0 *= torad;
  dphi *= torad;

  fWedge->SetGeometry (r0, r1, phi0, dphi);

  fOrigin[0] = (fWedge->GetOrigin())[0];
  fOrigin[1] = (fWedge->GetOrigin())[1];
  fOrigin[2] = z0;
  fSize[0] = (fWedge->GetSize())[0];
  fSize[1] = (fWedge->GetSize())[1];
  fSize[2] = depth;

  return kOK;
}


Int_t
TSolGEMChamber::Decode (const THaEvData& ed )
{
  for (UInt_t i = 0; i < GetNPlanes(); ++i)
    {
      GetPlane (i).Decode (ed);
    }
  return 0;
}

Int_t
TSolGEMChamber::InitPlane (const UInt_t i, const char* name, const char* desc)
{
  fPlanes[i] = new TSolGEMPlane (name, desc, this);
  fPlanes[i]->SetName (name);
  return fPlanes[i]->Init();
}

void
TSolGEMChamber::Print (const Bool_t printplanes)
{
  cout << "I'm a GEM chamber named " << GetName() << endl;
  TVector3 o (GetOrigin());
  cout << "  Origin: " << o(0) << " " << o(1) << " " << o(2)
       << " (rho,theta,phi)=(" << o.Mag() << "," << o.Theta()*TMath::RadToDeg()
       << "," << o.Phi()*TMath::RadToDeg() << ")"
       << endl;

#if ANALYZER_VERSION_CODE >= ANALYZER_VERSION(1,6,0)
  const Double_t* s = GetSize();
#else
  const Float_t* s = GetSize();
#endif
  cout << "  Size:   " << s[0] << " " << s[1] << " " << s[2] << endl;

  cout << "  Wedge geometry: r0: " << fWedge->GetR0()
       << " r1: " << fWedge->GetR1()
       << " phi0: " << fWedge->GetPhi0()*TMath::RadToDeg()
       << " dphi: " << fWedge->GetDPhi()*TMath::RadToDeg()
       << endl;

  if (printplanes)
    for (UInt_t i = 0; i < GetNPlanes(); ++i)
      {
	fPlanes[i]->Print();
      }
}

//__________________________________________________________________________________
Bool_t TSolGEMChamber::IsInDeadArea(Double_t& x, Double_t& y)
{
    //fitst rotate into tracker frame
    fVector3.SetXYZ(x, y, 0);
    fVector3.RotateZ(-1.*fWedge->GetAngle());
    
    //first check whether it is on the GEM frame
    Double_t a = -1.*tan(fWedge->GetDPhi()/2.);
    if (fabs(a*fVector3.X() + fVector3.Y())/sqrt(1.+a*a) < fFrameWidth) return true;
    a *= -1.;
    if (fabs(a*fVector3.X() + fVector3.Y())/sqrt(1.+a*a) < fFrameWidth) return true;
    
    for (UInt_t i=0; i<fHVSectorOff.size(); i++){
        if (fHVSectorOff[i].Contains(fVector3.X(), fVector3.Y())) return true;
    }
    
    return false;
}
