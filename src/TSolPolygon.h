#ifndef __TSOLPOLYGON_H
#define __TSOLPOLYGON_H

// Class for a polygon that can tell if a point is inside it
// Implemented entirely in header

#include <vector>
#include <iostream>

class TSolPolygon
{
 public:
  TSolPolygon() {};
  TSolPolygon (UInt_t polyCorners, std::vector<Double_t>& polyX, std::vector<Double_t>& polyY)
    {
      // polyCorners is number of vertices of polygon
      // polyX and polyY are vectors of vertex coordinates in plane 
      // coordinate frame
      
      // We precompute quantities to speed up Contains method
      // see http://alienryderflex.com/polygon/
      
      fPc = polyCorners;
      fPX = new Double_t[fPc];
      fPY = new Double_t[fPc];
      fConstant = new Double_t[fPc];
      fMultiple = new Double_t[fPc];
      UInt_t j = polyCorners - 1;
      for (UInt_t i = 0; i < polyCorners; i++) 
	{
	  fPX[i] = polyX[i];
	  fPY[i] = polyY[i];
	  if (polyY[j] == polyY[i]) 
	    {
	      fConstant[i] = polyX[i];
	      fMultiple[i] = 0; 
	    }
	  else 
	    {
	      fMultiple[i] = (polyX[j] - polyX[i]) / (polyY[j] - polyY[i]);
	      fConstant[i] = polyX[i] - polyY[i] * fMultiple[i];
	    }
	  j = i; 
	}
    }
  
  Bool_t Contains (const Double_t& x, const Double_t& y)
  {
    // Check whether x and y (in plane coordinate frame) are in dead
    // region. The function will return true if the point x,y is
    // inside the polygon, or false if it is not. If the point is
    // exactly on the edge of the polygon, then the function may
    // return true or false.
    
    UInt_t i;
    UInt_t j = fPc - 1;
    
    Bool_t  oddNodes = false;
    
    for (i = 0; i < fPc; i++) 
      {
	if ((fPY[i] < y && fPY[j] >= y)
	    ||   (fPY[j] < y && fPY[i] >= y))
	  oddNodes ^= (y * fMultiple[i] + fConstant[i] < x); 
	j = i; 
      }

    return oddNodes; 
  }

  void Print() const
  {
    for (UInt_t i = 0; i < fPc; ++i)
      std::cout << i << ": (" << fPX[i] << ", " << fPY[i] << ")" << std::endl;
  }

 private:  
  UInt_t fPc;
  Double_t* fPX;  // x coordinates not actually used directly but stored for printing
  Double_t* fPY;
  Double_t* fConstant;
  Double_t* fMultiple;
};

#endif//__TSOLPOLYGON_H
