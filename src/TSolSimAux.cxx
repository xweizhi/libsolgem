#include "TSolSimAux.h"
#include <iostream>
#include <cassert>
#include <TMath.h>

// ADC Conversion
//
//

Short_t 
TSolSimAux::ADCConvert(Double_t val, Double_t off, Double_t gain, Int_t bits)
{
  // Convert analog value 'val' to integer ADC reading with 'bits' resolution

  assert( bits >= 0 && bits <= MAX_ADCBITS );

  if( val < 0. )
    val = 0.;
  Double_t vvv = (val - off)/gain;

  Double_t saturation = static_cast<Double_t>( (1<<bits)-1 );
  if( vvv > saturation )
    vvv = saturation;

  Short_t dval =
    static_cast<Short_t>( TMath::Floor( (vvv>saturation) ? saturation : vvv ));

  //  cerr << val << " dval = " << dval << endl;
  if( dval < 0 ) dval = 0;
  return dval;

}

Int_t TSolSimAux::CheckSaturation(Int_t ADC, Int_t bits)
{
    assert( bits >= 0 && bits <= MAX_ADCBITS );
    
    Double_t saturation = static_cast<Double_t>( (1<<bits)-1 );
    
    if( ADC > saturation )
    ADC = saturation;
    
    return ADC;
}

// Pulse Shape SiD model
// APV25 time function from  M. Friedl et al NIMA 572 (2007) pg 385-387 (APV25 for silicon!)
//

Double_t 
TSolSimAux::PulseShape(Double_t t, 
				  Double_t C,  // normalization factor
				  Double_t Tp) // shaping time 

{
  Double_t v;
  Double_t x;
  x = t/Tp;
  v = C/Tp * x * TMath::Exp(-x);
  
  return ( v>0. ) ? v : 0.;

}

// GEM model
// from Frank Simon Thesis 2001 on COMPASS 
//  input delta function: A=157, t0=44, tau0=22, tau1=100
//  input from real GEM: A=400+/-60.48, t0=16.31+/-1.744
//                       tau0 = 38.28+/-3.604, tau1=129+/-6.499

Double_t 
TSolSimAux::PulseShape(Double_t t,
				  Double_t A,  // normalization
				  Double_t tau0,  // first time constant
				  Double_t tau1)  // second time constant
{
  
  Double_t v;
  Double_t x0,x1;

    //cerr << __FUNCTION__ << " " << t << " " << A << " " << tau0 << " " << tau1 << endl;
  if (tau1<0) { return PulseShape(t, A, tau0); } // SiD model
  x0 = -t/tau0;
  x1 = -t/tau1;
  v = A * ((tau0+tau1)/tau1/tau1)*(1.-TMath::Exp(x0)) * TMath::Exp(x1);
  return ( v>0. ) ? v : 0.;
 	   
}

Double_t TSolSimAux::SAMPAPulseShape(Double_t t, Double_t A, Int_t mode)
{
    //hard code the parameters here for now for test 160ns shaping
    Double_t p1 = 0.;
    Double_t p2 = 0.;
    Double_t p3 = 0.;
    Double_t p4 = 0.;
    
    if (mode == 1){
        //160ns shaping time
        p1 = 1.89019e+02;
        p2 = 2.33628e+00;
        p3 = 1.29642e+02;
        p4 = 1.36690e+00;
    }else{
        //80ns shaping time
        p1 = 1.17276e+02;
        p2 = 2.33628e+00;
        p3 = 6.48168e+01;
        p4 = 1.36690e+00;
    }
    Double_t v = A* TMath::Power(t/p1, p2) * TMath::Exp(-TMath::Power(t/p3, p4));
    return ( v>0. ) ? v : 0.;
}

Double_t TSolSimAux::VMMPulseShape(Double_t t, Double_t A, Double_t tau0) {
  if (t<0) return 0;
  Double_t PeakingTime = tau0;
  Double_t a = (PeakingTime*(10^-9))/1.5;
  Double_t pole0 = 1.263/a;  // real pole
  Double_t Re_pole1 = 1.149/a;
  Double_t Im_pole1 = -0.786/a; //complex pole
  Double_t K0 = 1.584;
  Double_t Re_K1 = -0.792;
  Double_t Im_K1 = -0.115;
  Double_t pole1_square = Re_pole1*Re_pole1 + Im_pole1*Im_pole1;
  Double_t K1_abs = TMath::Sqrt(Re_K1*Re_K1 + Im_K1*Im_K1);
  Double_t argK1 = TMath::ATan2(Im_K1, Re_K1);
  t*=(10^-9);
  Double_t st = A * TMath::Power(a,3)*pole0*pole1_square*((K0*TMath::Exp(-t*pole0))+(2.*K1_abs*TMath::Exp(-t*Re_pole1)*cos(-t*Im_pole1+argK1)));
  return (st>0.) ? st : 0.;

}

// par[0] = norm
// par[1] = x_mean
// par[2] = x_sigma
// par[3] = y_mean
// par[4] = y_sigma

Double_t 
TSolSimAux::Gaus2D(Double_t *x, Double_t *par)
{
  if(par[2] > 0 && par[4] > 0)
    {
      double rx=(x[0]-par[1])/par[2];
      double ry=(x[1]-par[3])/par[4];
      return par[0]*TMath::Exp(-(rx*rx+ry*ry)/2.);
    }
  else
    {
      return 0.;
    }
}

//
// Return the sum on ngaus 2D-Gaussian Functions
//
Double_t 
TSolSimAux::MultiGaus2D(Double_t *x, Double_t *par) {

  Int_t i;

  Double_t f;

  Int_t ngaus = (Int_t) par[0];

  f=0.;

  //  cerr << ngaus << " (" << x[0] << "," << x[1] << "): " ;

  for (i=0;i<ngaus;i++) {
    Double_t gg=Gaus2D(x,&par[5*i+1]);

    f += gg;

    //   cerr << gg << "+";
  }

  //  cerr << " = " << f << endl;

  return f;
}


Double_t 
TSolSimAux::SimpleCircle(Double_t *x, Double_t *par) {

  Double_t dx = x[0]-par[1];
  Double_t dy = x[1]-par[2];

  if ((dx*dx+dy*dy) < par[3]*par[3]) {
    return par[0];
  } else {
    return 0;
  }

}

Double_t 
TSolSimAux::MultiCircle(Double_t *x, Double_t *par) {

  Int_t i;

  Double_t f;

  Int_t n = (Int_t) par[0];

  f=0.;

  for (i=0;i<n;i++) {
    Double_t gg=SimpleCircle(x,&par[4*i+1]);

    f += gg;

  }

  return f;

}
