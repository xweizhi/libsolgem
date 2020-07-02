#include "TSolGEMPlane.h"

#include "TClonesArray.h"

#include "TSolGEMChamber.h"
#include "TSolGEMCluster.h"
#include "TSolWedge.h"
#include "THaEvData.h"
#include "TMath.h"
#include "ha_compiledata.h"

#include <iostream>
#include <cassert>
#include <iomanip>


using namespace std;

TSolGEMPlane::TSolGEMPlane()
  : THaSubDetector()
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  fWedge = new TSolWedge;
  fYDiv = NULL;
  return;
}

TSolGEMPlane::TSolGEMPlane( const char *name, const char *desc,
			    THaDetectorBase* parent )
  : THaSubDetector (name, desc, parent)
{
  //  fClusters = new TClonesArray("TSolGEMCluster", 100);  
  fWedge = new TSolWedge;
  fYDiv = NULL;
  return;
}

TSolGEMPlane::~TSolGEMPlane()
{
  //  delete fClusters;
  delete fWedge;
  delete[] fYDiv;
}

Int_t 
TSolGEMPlane::ReadDatabase (const TDatime& date)
{
  FILE* file = OpenFile (date);
  if (!file) return kFileError;

  Int_t err = ReadGeometry (file, date, false);

  fclose(file);
  if (err)
    return err;

  return kOK;
}

Int_t 
TSolGEMPlane::ReadGeometry (FILE* file, const TDatime& date,
			    Bool_t required)
{
  // Get x/y position, size, and angles from database if and only
  // if parent is null otherwise copy from parent

  // Note that origin is in lab frame, size is in wedge frame.

  Int_t err;
  Double_t torad = atan(1) / 45.0;

  TSolGEMChamber* parent = (TSolGEMChamber*) GetParent();
  Double_t z0;
  Double_t depth;
  if (parent != NULL)
    {
      fOrigin = parent->GetOrigin();
      fSize[0] = (parent->GetSize())[0];
      fSize[1] = (parent->GetSize())[1];
      fSize[2] = (parent->GetSize())[2];
      fWedge->SetGeometry (parent->GetWedge().GetR0(),
			   parent->GetWedge().GetR1(),
			   parent->GetWedge().GetPhi0(),
			   parent->GetWedge().GetDPhi());

      z0 = fOrigin[2];
      depth = fSize[2];
    }
  else
    {
      Double_t r0 = -999.0;
      Double_t r1 = -999.0;
      Double_t phi0 = -999.0;
      Double_t dphi = -999.0;
      const DBRequest request[] = 
	{
	  {"r0",          &r0,           kDouble, 0, 1},
	  {"r1",          &r1,           kDouble, 0, 1},
	  {"phi0",        &phi0,         kDouble, 0, 1},
	  {"dphi",        &dphi,         kDouble, 0, 1},
	  {0}
	};
      err = LoadDB( file, date, request, fPrefix );
      
      if (err)
	return err;
      //cout<<"from database: "<<phi0<<" "<<dphi<<endl;
      // Database specifies angles in degrees, convert to radians
      phi0 *= torad;
      dphi *= torad;
      
      fWedge->SetGeometry (r0, r1, phi0, dphi);

      fOrigin[0] = (fWedge->GetOrigin())[0];
      fOrigin[1] = (fWedge->GetOrigin())[1];
      fSize[0] = (fWedge->GetSize())[0];
      fSize[1] = (fWedge->GetSize())[1];
      z0 = -999.0;
      depth = -999.0;
    }

  vector<Double_t> DivSegment (4);
  const DBRequest request[] = 
    {
      {"stripangle",  &fSAngle,      kDouble, 0, 1},
      {"pitch",       &fSPitch,      kDouble, 0, 1},
      {"divsegment",  &DivSegment,   kDoubleV, 0, 1},
      {"z0",          &z0,           kDouble, 0, 1},
      {"depth",       &depth,        kDouble, 0, 1},
      {0}
    };
  err = LoadDB( file, date, request, fPrefix );

  if (err)
    return err;

  fSAngle *= torad;

  SetRotations();
  fOrigin[2] = z0;
  fSize[2] = depth;
  
  Double_t rorigin = sqrt (fOrigin[0]*fOrigin[0]+fOrigin[1]*fOrigin[1]);
  fSBeg = rorigin * sin (fWedge->GetDPhi()/2);

  // Get numbers of strips

  Double_t xs0 = fWedge->GetR1() * cos (fWedge->GetDPhi()/2) - rorigin;
  Double_t ys0 = fWedge->GetR1() * sin (fWedge->GetDPhi()/2);
#define OTHERWAY
#ifdef OTHERWAY
  if (GetSAngle() < TMath::Pi()/2)
    ys0 = -ys0;
#else
  if (GetSAngle() > TMath::Pi()/2)
    ys0 = -ys0;
#endif

  PlaneToStrip (xs0, ys0);
#ifdef OTHERWAY
  if (GetSAngle() < TMath::Pi()/2)
    fNStrips = int ((fSBeg + xs0) / GetSPitch());
  else
    fNStrips = int ((-xs0 + fSBeg) / GetSPitch());
#else
  if (GetSAngle() > TMath::Pi()/2)
    fNStrips = int ((fSBeg - xs0) / GetSPitch());
  else
    fNStrips = int ((xs0 + fSBeg) / GetSPitch());
#endif

  // Make table of strip divisions

  // Get affected strips (the ones containing the segment endpoints)
  Int_t ds[2];
  Double_t xds[2];
  Double_t yds[2];
  for (UInt_t i = 0; i < 2; ++i)
    {
      xds[i] = DivSegment.at(2*i+0);
      yds[i] = DivSegment.at(2*i+1);
      PlaneToStrip (xds[i], yds[i]);      
      ds[i] = GetStripInRange (xds[i]);
    }
  if (ds[0] > ds[1])
    {
      swap (ds[0], ds[1]);
      swap (xds[0], xds[1]);
      swap (yds[0], yds[1]);
    }

  if (ds[1] > ds[0])
    {
      // cout << "Division segment in strip frame is (" << xds[0] << ", " << yds[0]
      // 	   << ") to (" << xds[1] << ", " << yds[1] << ")" << endl;
      // cout << "Division strips are " << ds[0] << " to " << ds[1] << endl;

      fNDiv = ds[1] - ds[0] + 1;
      fSDiv0 = ds[0];
      Double_t ms = (yds[1] - yds[0]) / (xds[1] - xds[0]);
      fYDiv = new Double_t[fNDiv];
      for (UInt_t i = 0; i < fNDiv; ++i)
	{
	  fYDiv[i] =  yds[0] + ms * (GetStripCenter (i+fSDiv0) - xds[0]);
	  // cout << "Strip " << i+fSDiv0 
	  //      << " center at " << GetStripCenter (i+fSDiv0) 
	  //      << " divides at " << fYDiv[i]
	  //      << endl;
	}
      // Print();

    }
  else
    {
      fNDiv = 0;
      fSDiv0 = 0;
    }

  return kOK;
}

Int_t TSolGEMPlane::Decode( const THaEvData &d ){
    // Clusters get made as so

  //    int i = 0;

    //    new ((*fClusters)[i]) TSolGEMCluster();

    return 0;
}

Double_t 
TSolGEMPlane::GetSAngle()   const
{
  return fSAngle;
}

void
TSolGEMPlane::LabToStrip (Double_t& x, Double_t& y) const
{
  x -= (GetOrigin())[0];
  y -= (GetOrigin())[1];
  Double_t temp = x;
  x = fCLS * x - fSLS * y;
  y = fSLS * temp + fCLS * y;
  return;
}

void
TSolGEMPlane::StripToPlane (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCWS * x + fSWS * y;
  y = -fSWS * temp + fCWS * y;
  return;
}

void
TSolGEMPlane::StripToLab (Double_t& x, Double_t& y) const
{
  Double_t temp = x;
  x = fCLS * x + fSLS * y;
  y = -fSLS * temp + fCLS * y;
  x += (GetOrigin())[0];
  y += (GetOrigin())[1];
  return;
}

Double_t TSolGEMPlane::StripNumtoStrip( Int_t strip )
{
    // Gives x coordinate in strip frame of a wire
  return (GetStripCenter (strip));
}


Double_t TSolGEMPlane::StriptoProj( Double_t s )
{
    // Gives coordinate in projection frame from strip frame x
    Double_t r = (GetWedge().GetR1()-GetWedge().GetR0())/2.0;
    return s + r*fCWS;
}


Double_t TSolGEMPlane::StripNumtoProj( Int_t s ){
    // Gives coordinate in projection frame from strip number
    return StriptoProj( StripNumtoStrip(s) );
}

Double_t 
TSolGEMPlane::GetStripLowerEdge (UInt_t is) const 
{
  if (GetSAngle() > TMath::Pi()/2)
    return fSBeg - (fNStrips - is) * GetSPitch();
  else
    return -fSBeg + is * GetSPitch();
}

Double_t 
TSolGEMPlane::GetStripUpperEdge (UInt_t is) const {return GetStripLowerEdge (is) + GetSPitch();}

Double_t 
TSolGEMPlane::GetStripCenter (UInt_t is) const {return GetStripLowerEdge (is) + 0.5 * GetSPitch();}

Bool_t
TSolGEMPlane::IsDivided (UInt_t is) const 
{
  return (is >= fSDiv0 && is < fSDiv0+fNDiv);
}

Double_t 
TSolGEMPlane::GetYDiv (UInt_t is) const 
{
  return IsDivided (is) ? fYDiv[is-fSDiv0] : 1e9;
}

Int_t
TSolGEMPlane::GetStripUnchecked( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame,
  // no questions asked. Caller must check return value

  if (GetSAngle() > TMath::Pi()/2)
    return  (fNStrips - 1) - int ((fSBeg - x) / GetSPitch());
  else
    return int ((x + fSBeg) / GetSPitch());
}

Int_t
TSolGEMPlane::GetStripInRange( Double_t x ) const
{
  // Get strip number for given x-coordinate in strip frame
  // and, if out of range, limit it to allowable values.

  Int_t s = GetStripUnchecked(x);
  if( s < 0 )              s = 0;
  if( s >= GetNStrips() )  s = GetNStrips()-1;
  return s;
}
    
Int_t
TSolGEMPlane::GetStrip (Double_t x, Double_t yc) const
{
  // Strip number corresponding to coordinates x, y in 
  // strip frame, or -1 if outside (2-d) bounds
  
  // Double_t xw = x;
  // Double_t yw = yc;
  // StripToPlane (xw, yw);

  Double_t xc = x;
  StripToLab (xc, yc);

  if (!fWedge->Contains (xc, yc))
    return -1;

  Int_t s = GetStripUnchecked(x);
  // cout << GetName()
  //      << " " << s << " " << x << " " << GetNStrips() 
  //      << " " << xc << " " << yc
  //      << " " << xw << " " << yw << endl;
  //      assert (!(TString(GetName()) == "gem92y"));

  assert( s >= 0 && s < GetNStrips() ); // by construction in ReadGeometry()
  return s;
}

Int_t 
TSolGEMPlane::GetDivision (UInt_t is, Double_t y) const
{
  return IsDivided (is) ? (y > fYDiv[is-fSDiv0] ? 1 : 0) : 0;  
}

void 
TSolGEMPlane::Print() const
{
  cout << "I'm a GEM plane named " << GetName() << endl;

  TVector3 o (GetOrigin());
  cout << "  Origin: " <<  setprecision(4) << o(0) << " " << o(1) << " " << o(2)
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

  cout << "  Wedge size: +-" << (fWedge->GetSize())[0] 
       << ", +-" << (fWedge->GetSize())[1] << endl;

  cout << "  " << GetNStrips() << " strips"
       << ", angle " << GetSAngle()*TMath::RadToDeg()
       << ", start " << GetStripCenter(0)
       << ", end " << GetStripCenter(GetNStrips()-1)
       << ", pitch " << GetSPitch()
       << endl;

  Int_t nd = GetNDividedStrips();
  Int_t sd0 = GetFirstDividedStrip();
  if (nd > 0)
    cout << "  " << nd << " divided strips"
	 << ", from " << sd0 
	 << " at " << GetYDiv (sd0)
	 << " to " << sd0 + nd - 1
	 << " at " << GetYDiv (sd0 + nd - 1)
	 << endl;

  Double_t xxx[4] = { -(fWedge->GetSize())[0],  (fWedge->GetSize())[0],  (fWedge->GetSize())[0],  -(fWedge->GetSize())[0] };
  Double_t yyy[4] = { -(fWedge->GetSize())[1],  -(fWedge->GetSize())[1],  (fWedge->GetSize())[1],  (fWedge->GetSize())[1] };
  Double_t xxxs[4];
  Double_t yyys[4];
  for (UInt_t ixxx = 0; ixxx < 4; ++ixxx)
    {
      xxxs[ixxx] = xxx[ixxx];
      yyys[ixxx] = yyy[ixxx];
      PlaneToStrip (xxxs[ixxx], yyys[ixxx]);
      cout << "Plane vertex " << xxx[ixxx] << " " << yyy[ixxx] 
	   << " -> Strip vertex " << xxxs[ixxx] << " " << yyys[ixxx] 
	   << endl;
    }
}

void
TSolGEMPlane::SetRotations()
{
  // Set rotation angle trig functions
  //cout<<"strip angles: "<<GetSAngle()<<" "<<GetAngle()<<endl;
#ifdef OTHERWAY
  fCWS = cos (GetSAngle());
  fCLS = cos (-GetAngle()+GetSAngle());
  fSWS = sin (GetSAngle());
  fSLS = sin (-GetAngle()+GetSAngle());
#else
  fCWS = cos (-GetSAngle());
  fCLS = cos (-GetAngle()-GetSAngle());
  fSWS = sin (-GetSAngle());
  fSLS = sin (-GetAngle()-GetSAngle());
#endif
}
