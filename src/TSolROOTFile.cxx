#include "TSolROOTFile.h"
#include "TMath.h"
#include "assert.h"
#include <algorithm>
#include <cstdlib>

#define NGENERATED 7
#define NSOLID_GEM 23
#define NECDATA 8
//#ifndef __CINT__
//_______________________________________________________________________
TSolROOTFile::TSolROOTFile() : fChan(0), fSource(0) {
  fFilename[0] = '\0';
  fManager = TSolDBManager::GetInstance();
}
//_______________________________________________________________________
TSolROOTFile::TSolROOTFile(const char *f, int source) 
: fChan(0), fSource(source) {
    SetFilename(f);
    gen_pid = 0;
    gen_px  = 0;
    gen_py  = 0;
    gen_pz  = 0;
    gen_vx  = 0;
    gen_vy  = 0;
    gen_vz  = 0;

    solid_gem_id = 0;
    solid_gem_hitn = 0;
    solid_gem_pid = 0;
    solid_gem_trid = 0;
    solid_gem_x = 0;
    solid_gem_y = 0;
    solid_gem_z = 0;
    solid_gem_lxin = 0;
    solid_gem_lyin = 0;
    solid_gem_lzin = 0;
    solid_gem_tin = 0;
    solid_gem_lxout = 0;
    solid_gem_lyout = 0;
    solid_gem_lzout = 0;
    solid_gem_tout = 0;
    solid_gem_px = 0;
    solid_gem_py = 0;
    solid_gem_pz = 0;
    solid_gem_vx = 0;
    solid_gem_vy = 0;
    solid_gem_vz = 0;
    solid_gem_ETot = 0;
    solid_gem_trE = 0;
    solid_gem_weight = 0;

    flux_id = 0;
    flux_pid = 0;
    flux_tid = 0;
    flux_trackE = 0;
    flux_avg_x = 0;
    flux_avg_y = 0;
    flux_avg_z = 0;
    
    header_var8 = 0;
    
    fSignalInfo.clear();
    fManager = TSolDBManager::GetInstance();
    
    if (fSource == 0){
        for (int i=0; i<fManager->GetNSigParticle(); i++){
            fSignalInfo.push_back(SignalInfo(fManager->GetSigPID(i), 
                                             fManager->GetSigTID(i)));
        }
    }
}
//_______________________________________________________________________
TSolROOTFile::~TSolROOTFile() {
  Clear();
  DeleteVectors();
  delete fChan;
}
//_______________________________________________________________________
void TSolROOTFile::SetFilename( const char *f ){
  if( !f ) return;
  strcpy( fFilename, f );
}
//_______________________________________________________________________
Int_t TSolROOTFile::Open(){
    // Return 0 on fail, 1 on success
    if( fFilename[0] == '\0' ){ return 0; }
    fChan = new TFile(fFilename, "READ");

    if (fChan->IsZombie()) {
      std::cout <<"Error opening file "<<fFilename<<std::endl;
      return 0;
    }

    fEvNum = -1;
    
    SetBranchAddress();
    
    return 1;
}
//_______________________________________________________________________
void TSolROOTFile::SetBranchAddress()
{
  //generated branch
  tree_generated = (TTree*)fChan->Get("generated");
  tree_solid_gem = (TTree*)fChan->Get("solid_gem");
  tree_flux      = (TTree*)fChan->Get("flux");
  tree_header    = (TTree*)fChan->Get("header");
  
  tree_generated->SetBranchAddress("pid",&gen_pid);
  tree_generated->SetBranchAddress("px",&gen_px);
  tree_generated->SetBranchAddress("py",&gen_py);
  tree_generated->SetBranchAddress("pz",&gen_pz);
  tree_generated->SetBranchAddress("vx",&gen_vx);
  tree_generated->SetBranchAddress("vy",&gen_vy);
  tree_generated->SetBranchAddress("vz",&gen_vz);
  
  //solid_gem branch
  tree_solid_gem->SetBranchAddress("hitn",&solid_gem_hitn);
  tree_solid_gem->SetBranchAddress("id",&solid_gem_id);
  tree_solid_gem->SetBranchAddress("pid",&solid_gem_pid);
  tree_solid_gem->SetBranchAddress("trid",&solid_gem_trid);
  tree_solid_gem->SetBranchAddress("x",&solid_gem_x);
  tree_solid_gem->SetBranchAddress("y",&solid_gem_y);
  tree_solid_gem->SetBranchAddress("z",&solid_gem_z);
  tree_solid_gem->SetBranchAddress("lxin",&solid_gem_lxin);
  tree_solid_gem->SetBranchAddress("lyin",&solid_gem_lyin);
  tree_solid_gem->SetBranchAddress("lzin",&solid_gem_lzin);
  tree_solid_gem->SetBranchAddress("tin",&solid_gem_tin);
  tree_solid_gem->SetBranchAddress("lxout",&solid_gem_lxout);
  tree_solid_gem->SetBranchAddress("lyout",&solid_gem_lyout);
  tree_solid_gem->SetBranchAddress("lzout",&solid_gem_lzout);
  tree_solid_gem->SetBranchAddress("tout",&solid_gem_tout);
  tree_solid_gem->SetBranchAddress("px",&solid_gem_px);
  tree_solid_gem->SetBranchAddress("py",&solid_gem_py);
  tree_solid_gem->SetBranchAddress("pz",&solid_gem_pz);
  tree_solid_gem->SetBranchAddress("vx",&solid_gem_vx);
  tree_solid_gem->SetBranchAddress("vy",&solid_gem_vy);
  tree_solid_gem->SetBranchAddress("vz",&solid_gem_vz);
  tree_solid_gem->SetBranchAddress("ETot",&solid_gem_ETot);
  tree_solid_gem->SetBranchAddress("trE",&solid_gem_trE);
  tree_solid_gem->SetBranchAddress("weight",&solid_gem_weight);
  //flux branch
  tree_flux->SetBranchAddress("id",&flux_id);
  tree_flux->SetBranchAddress("pid",&flux_pid);
  tree_flux->SetBranchAddress("tid",&flux_tid);
  tree_flux->SetBranchAddress("trackE",&flux_trackE);
  tree_flux->SetBranchAddress("avg_x",&flux_avg_x);
  tree_flux->SetBranchAddress("avg_y",&flux_avg_y);
  tree_flux->SetBranchAddress("avg_z",&flux_avg_z);
  
  //header branch
  tree_header->SetBranchAddress("var8", &header_var8);

  fMaxEvNum = tree_flux->GetEntries();  
}
//_______________________________________________________________________
Int_t TSolROOTFile::Close(){
    // Return 0 on fail, 1 on success
    Int_t ret = 1;
    fChan->Close();
    DeleteVectors();
    delete fChan; fChan = 0;
    return ret;
}
//_______________________________________________________________________
void TSolROOTFile::Clear(){
    // Clear out hit and generated data

    unsigned int i;
    for( i = 0; i < fHitData.size(); i++ ){
	delete fHitData[i];
    }

    for( i = 0; i < fGenData.size(); i++ ){
	delete fGenData[i];
    }
    
    for( i=0; i < fECData.size(); i++){
    	delete fECData[i];
    }

    fHitData.clear();
    fGenData.clear();
    fECData.clear();
	ClearVectors();
    return;
}

//_______________________________________________________________________
Int_t TSolROOTFile::ReadNextEvent(){
  Clear();
  if (fChan->IsZombie()) {
    std::cout <<"Error opening file "<<fFilename<<std::endl;
    return 0;
  }
  //cout<<"what's the event number again? "<<fEvNum<<endl;
  fEvNum++;
  //cout<<"event number after"<<fEvNum<<endl;
  if (fEvNum >= fMaxEvNum || fEvNum < 0 ){
    return 0;
  }
   
  //clear signal info array
  for (unsigned int i = 0; i < fSignalInfo.size(); i++){
    fSignalInfo.at(i).fillBitsGEM = 0;
    fSignalInfo.at(i).fillBitsEC = 0;
    fSignalInfo.at(i).ECEDep = 0.;
    fSignalInfo.at(i).momentum = 0.;
    fSignalInfo.at(i).R = 0.;
    fSignalInfo.at(i).signalSector = -1;
  }

  tree_generated->GetEntry(fEvNum);
  tree_solid_gem->GetEntry(fEvNum);
  tree_flux->GetEntry(fEvNum);
  tree_header->GetEntry(fEvNum);
  BuildGenerated();
  BuildECData();
  ExtractDetIDs();
  BuildData();
  return 1;
}
//_______________________________________________________________________
void TSolROOTFile::BuildGenerated()
{
  unsigned int length = gen_pid->size();
  unsigned int i = 0;

  for (i=0; i<length; i++){
    fGenData.push_back( new gendata() );
  }
  
  for (i=0; i<length; i++){
    if (fSource == 0){
        fGenData[i]->SetData(0, gen_pid->at(i));
        fGenData[i]->SetData(1, gen_px->at(i));
        fGenData[i]->SetData(2, gen_py->at(i));
        fGenData[i]->SetData(3, gen_pz->at(i));
        fGenData[i]->SetData(4, gen_vx->at(i));
        fGenData[i]->SetData(5, gen_vy->at(i));
        fGenData[i]->SetData(6, gen_vz->at(i));
    //assuming fSignalInfo is arranged according to tid
    }else{
        fGenData[i]->SetData(0, gen_pid->at(i));
        fGenData[i]->SetData(1, 0);
        fGenData[i]->SetData(2, 0);
        fGenData[i]->SetData(3, 0);
        fGenData[i]->SetData(4, 0);
        fGenData[i]->SetData(5, 0);
        fGenData[i]->SetData(6, 0);
    }
    if ( fSource != 0) continue;
    if ( gen_pid->at(i) == fSignalInfo.at(i).pid)
    fSignalInfo.at(i).momentum = sqrt(pow(gen_px->at(i),2) + pow(gen_py->at(i),2)
                                    + pow(gen_pz->at(i),2))*1.e-3; //To GeV

  }
  
  if ( fSource != 0 ) return;
  
  length = header_var8->size() > gen_pid->size() ? gen_pid->size() : header_var8->size();
  i=0;
  
  //sometimes there is empty entry in the generated branch but not in the header branch
  //the following assert will then result in a crush
  //assert(length == fGenData.size());
  
  for (i=0; i<length; i++){
    fGenData[i]->SetData(7, header_var8->at(i));
  }
}
//_______________________________________________________________________
void TSolROOTFile::ExtractDetIDs()
{
  unsigned int length = solid_gem_id->size();
  unsigned int i=0;
  for (i=0; i<length; i++){
    int detID = TMath::Nint(solid_gem_id->at(i));
    fHitData.push_back( new hitdata(detID, NSOLID_GEM) );
    //calculate the hit pattern, particle is considered hitting the GEM chamber if it enters the 3mm gas layer
    /*if ( (TMath::Nint(solid_gem_id->at(i)) % 100) == fManager->GetGEMCopperFrontID() && fSource == 0 ){
        int signalID = -1;
        for (unsigned int j=0; j<fSignalInfo.size(); j++){
          if (TMath::Nint(solid_gem_pid->at(i)) == fSignalInfo.at(j).pid &&
              TMath::Nint(solid_gem_trid->at(i)) == fSignalInfo.at(j).tid ) signalID = j;
        }

        if (signalID < 0) continue; //not found in the signal info array
        //at the moment the signal sector is defined as the first GEM sector that the primary particle hits
        //Maybe using the EC to define this will be better?

    	int chamberID = CalcHitChamber(detID) - 1;
    	int planeID = chamberID / fManager->GetNSector();
    	assert(planeID >= 0); // should never happen, planeID should start from 0
    	fSignalInfo.at(signalID).fillBitsGEM |= 1<<planeID;
        int n = 0;
        //find the first fill bit, should has one already filled by now by the above
        assert(fSignalInfo.at(signalID).fillBitsGEM != 0);
        while(true){
          if ((fSignalInfo.at(signalID).fillBitsGEM>>n & 1) == 1) break;
          n++;
        }
        //if planeID is by far the most upstream GEM being hit then:
        if (n == planeID) fSignalInfo.at(signalID).signalSector = chamberID % fManager->GetNSector();
    }*/

  }
}
//_______________________________________________________________________
void TSolROOTFile::BuildECData()
{
	unsigned int length = flux_id->size();
	unsigned int i=0;
	int nhit = 0;
	for (i=0; i<length; i++){
	  if (TMath::Nint(flux_id->at(i)) == fManager->GetFAECID() || TMath::Nint(flux_id->at(i)) == fManager->GetLAECID()){
	    if ( fSource != 0 ) continue; //not from signal data file

            int signalID = -1;
            for (unsigned int j=0; j<fSignalInfo.size(); j++){
              if (TMath::Nint(flux_pid->at(i)) == fSignalInfo.at(j).pid &&
                  TMath::Nint(flux_tid->at(i)) == fSignalInfo.at(j).tid ) signalID = j;
            }
            if (signalID < 0) continue; //not found in the signal info array


	    if ( TMath::Nint(flux_id->at(i)) == fManager->GetFAECID() ) fSignalInfo.at(signalID).fillBitsEC |= 1;
	    if ( TMath::Nint(flux_id->at(i)) == fManager->GetLAECID() ) fSignalInfo.at(signalID).fillBitsEC |= (1<<1);

	    fECData.push_back( new ECdata( TMath::Nint(flux_id->at(i)) ) );
	    fECData[nhit]->SetData(0, TMath::Nint(flux_pid->at(i)) );
	    fECData[nhit]->SetData(1, TMath::Nint(flux_tid->at(i)) );
	    fECData[nhit]->SetData(2, flux_trackE->at(i)*1e-3);
	    fECData[nhit]->SetData(3, flux_avg_x->at(i));
	    fECData[nhit]->SetData(4, flux_avg_y->at(i));
	    fECData[nhit]->SetData(5, flux_avg_z->at(i));
        fSignalInfo.at(signalID).ECEDep = flux_trackE->at(i)*1e-3;
        fSignalInfo.at(signalID).R = sqrt(pow(flux_avg_x->at(i) ,2) + pow(flux_avg_y->at(i) ,2))*1.e-3;
	    nhit++;
	  }
        }
}
//_______________________________________________________________________
void TSolROOTFile::BuildData()
{
  unsigned int length = solid_gem_id->size();
  unsigned int i=0;

  for (i=0; i<length; i++){
    fHitData[i]->SetData(1, solid_gem_ETot->at(i)); 
    fHitData[i]->SetData(2, solid_gem_x->at(i));
    fHitData[i]->SetData(3, solid_gem_y->at(i));
    fHitData[i]->SetData(4, solid_gem_z->at(i));
    fHitData[i]->SetData(5, solid_gem_lxin->at(i));
    fHitData[i]->SetData(6, solid_gem_lyin->at(i));
    fHitData[i]->SetData(7, solid_gem_lzin->at(i));
    fHitData[i]->SetData(8, solid_gem_tin->at(i));
    fHitData[i]->SetData(9, solid_gem_lxout->at(i));
    fHitData[i]->SetData(10, solid_gem_lyout->at(i));
    fHitData[i]->SetData(11, solid_gem_lzout->at(i));
    fHitData[i]->SetData(12, solid_gem_tout->at(i));
    fHitData[i]->SetData(13, solid_gem_pid->at(i));
    fHitData[i]->SetData(14, solid_gem_vx->at(i));
    fHitData[i]->SetData(15, solid_gem_vy->at(i));
    fHitData[i]->SetData(16, solid_gem_vz->at(i));
    fHitData[i]->SetData(17, solid_gem_trE->at(i));
    fHitData[i]->SetData(18, solid_gem_trid->at(i));
    fHitData[i]->SetData(19, solid_gem_weight->at(i));
    fHitData[i]->SetData(20, solid_gem_px->at(i));
    fHitData[i]->SetData(21, solid_gem_py->at(i));
    fHitData[i]->SetData(22, solid_gem_pz->at(i));
  }

  //save the weight factor into fGenData
  /*for( i = 0; i < fGenData.size() && length > 0; i++ ){
    fGenData[i]->SetData(7, solid_gem_weight->at(0));
  }*/
}

//_______________________________________________________________________
TSolGEMData* TSolROOTFile::GetGEMData()
{
  // Return TSolGEMData object filled with GEM data of present event.
  // The returned object pointer must be deleted by the caller!

  TSolGEMData* gd = new TSolGEMData();

  GetGEMData(gd);
  return gd;
}
//_______________________________________________________________________
void TSolROOTFile::GetGEMData(TSolGEMData* gd)
{
  // Pack data into TSolGEMData

//    printf("NEXT EVENT ---------------------------\n");

    if( !gd ) return;
    gd->ClearEvent();
    gd->SetSource(fSource);
    gd->SetEvent(fEvNum);

    if (GetNData() == 0) {
      return;
    }
    gd->InitEvent(GetNData());

    hitdata *h, *hs;
    bool matchedstrip;
    unsigned int i, j, ngdata = 0;
    for( i = 0; i < GetNData(); i++ ){
      h = GetHitData(i);

      UInt_t tmpChamberID = CalcHitChamber(h->GetDetID()) - 1;
      int    tmpTrackerID = tmpChamberID / fManager->GetNSector();
      if ( tmpTrackerID < 0 || tmpTrackerID >= fManager->GetNTracker() ) continue;


      //detector id as defined in GEMC2 should be
      //$id=1000000+$n*100000+$sec*1000+$i where n is the plane #, sec is the
      //sector # and i is the layer #
      //this way of finding GEM_DRIFT is still good
      if( h->GetDetID()%100 == fManager->GetGEMDriftID() &&  h->GetData(1)>0.0 ){
	// Vector information
	TVector3 p(h->GetData(20), h->GetData(21), h->GetData(22));
	
	gd->SetMomentum(ngdata, p);

	TVector3 li(h->GetData(5), h->GetData(6), h->GetData(7));
	
	gd->SetHitEntrance(ngdata, li);

	TVector3 lo(h->GetData(9), h->GetData(10), h->GetData(11));
	
	gd->SetHitExit(ngdata, lo);
	// Average over entrance and exit time
	gd->SetHitTime(ngdata, (h->GetData(8)+h->GetData(12))/2.0);

	TVector3 vert(h->GetData(14), h->GetData(15), h->GetData(16));
	
	gd->SetVertex(ngdata, vert);

	gd->SetHitEnergy(ngdata, h->GetData(1)*1e6 ); // Gives eV
	gd->SetParticleID(ngdata, (UInt_t) h->GetData(18) );

	gd->SetParticleType(ngdata, (UInt_t) h->GetData(13) );

	////////////////////////////////////////////
	// Search for entrance/exit hits in surrounding Cu plane

	for( j = 0; j < GetNData(); j++ ){
	  hs = GetHitData(j);

	  if( hs->GetDetID()%100 == fManager->GetGEMCopperFrontID() &&    // is prior Cu plane
	      hs->GetDetID()/100 == h->GetDetID()/100 && // same detector
	      ((UInt_t) hs->GetData(18)) == ((UInt_t) h->GetData(18))  // same particle
	      ){
	    // Found matching hit, replace entrance data
	    li = TVector3(hs->GetData(5), hs->GetData(6), hs->GetData(7));
	    
	    gd->SetHitEntrance(ngdata, li);
	    break;
	  }
	}

	for( j = 0; j < GetNData(); j++ ){
	  hs = GetHitData(j);

	  if( hs->GetDetID()%100 == fManager->GetGEMCopperBackID() &&    // is subsequent Cu plane
	      hs->GetDetID()/100 == h->GetDetID()/100 && // same detector
	      ((UInt_t) hs->GetData(18)) == ((UInt_t) h->GetData(18))  // same particle
	      ){
	    // Found matching hit, replace exit data
	    lo = TVector3(hs->GetData(9), hs->GetData(10), hs->GetData(11));
	    
	    gd->SetHitExit(ngdata, lo);
	    break;
	  }
	}

	////////////////////////////////////////////
	// Search other hits for the corresponding 
	// hit on the strip


	matchedstrip = false;
	for( j = 0; j < GetNData(); j++ ){
	  hs = GetHitData(j);

	  if( hs->GetDetID()%100 == fManager->GetGEMStripID() &&    // is strip plane
	      hs->GetDetID()/100 == h->GetDetID()/100 && // same detector
	      ((UInt_t) hs->GetData(18)) == ((UInt_t) h->GetData(18))  // same particle
	      ){
	    if( !matchedstrip ){
	      // This is the truth information
	      TVector3 lr(hs->GetData(2), hs->GetData(3), hs->GetData(4));
	      
	      gd->SetHitReadout(ngdata, lr);
	      matchedstrip = true;
	    } else {
	      std::cout<<"Found multiple readout plane hits matching drift hit.  Truth information may be inaccurate"<<std::endl;
	    }
	  } 
	}

	if( !matchedstrip && (UInt_t) h->GetData(18) == 1 ){
        // FIXME
        // SPR 12/2/2011
        // This isn't the greatest way to do this but we're usually intersted in 
        // Particle 1 when we're looking at doing tracking.  Not all things depositing
        // energy leave a corresponding hit in the cathode plane.  Maybe we can look at
        // the mother IDs or something later.  

         // fprintf(stderr, "%s %s line %d: Did not find readout plane hit corresponding to drift hits.  No truth information\n",
        // 	  __FILE__, __FUNCTION__, __LINE__);
        TVector3 lr(-1e9, -1e9, -1e9);
        gd->SetHitReadout(ngdata, lr);
	}
	
	// Chamber ID starts indexing a 0 whereas we start conventionally at 1
    // hit chamber need to be determined differently for the new detector
    // convention
        
    //in the following there are two ways, one use the chamber id definition in the
    //gemc simulation, one use self-defined (through db_gemc.dat and db_generalinfo.dat)
       
    UInt_t hitChamberID = CalcHitChamber(h->GetDetID()) - 1;
    int    hitTrackerID = hitChamberID / fManager->GetNSector();
    assert( hitTrackerID >= 0 && hitTrackerID < fManager->GetNTracker() );
    
    if (fManager->DoSelfDefineSector()){
        double x = (li.X() + lo.X()) / 2. *1.e-3; // mm to m
        double y = (li.Y() + lo.Y()) / 2. *1.e-3; // mm to m
        int thisID = fManager->GetSectorIDFromPos(x, y, hitTrackerID);
        //if the particle doesn't hit any self-defined sector, skip it
        
        if (thisID < 0) continue; 
        if (thisID < 0) cout<<hitChamberID<<" "<<hitTrackerID<<endl;
        hitChamberID = hitTrackerID*fManager->GetNSector() + thisID;
        
    }
        
    gd->SetHitChamber(ngdata, hitChamberID );
    
    gd->SetPrimary(ngdata, false);
    if (fSource == 0){
        for (unsigned int idx = 0; idx < fSignalInfo.size(); idx++){
            if ( TMath::Nint(h->GetData(13)) == fSignalInfo[idx].pid && 
                 TMath::Nint(h->GetData(18)) == fSignalInfo[idx].tid ){
                gd->SetPrimary(ngdata, true);
      
                //update GEM hit pattern if the particle is a signal partical
                fSignalInfo.at(idx).fillBitsGEM |= 1<<hitTrackerID;
                int n = 0;
                //find the first fill bit, should has one already filled by now by the above
                assert(fSignalInfo.at(idx).fillBitsGEM != 0);
                while(true){
                    if ((fSignalInfo.at(idx).fillBitsGEM>>n & 1) == 1) break;
                    n++;
                }
                //if planeID is by far the most upstream GEM being hit then:
                if (n == hitTrackerID) 
                fSignalInfo.at(idx).signalSector = hitChamberID % fManager->GetNSector();
            }
        }
    }
	ngdata++;
      }
    }
    gd->SetNHit(ngdata);
    
    //if there is more then one signal particle, it doesn't make sense to map
    //sector, and so signal sector is not useful anyway
    if (fSignalInfo.size() > 0) gd->SetSigSector(fSignalInfo[0].signalSector);
}

//_______________________________________________________________________
UInt_t TSolROOTFile::CalcHitChamber(int detid)
{
  vector<int> digi;
  UInt_t id = detid;
  while(id>0){
    digi.push_back(id%10);
    id = id/10;
  }
  std::reverse(digi.begin(),digi.end());
  assert(digi.at(0) == 1 && "GEM detector id shoud be like 1xxxxxx");
  
  return (digi.at(1)-1)*fManager->GetNSector()+digi.at(2)*10+digi.at(3);
}
//_______________________________________________________________________
Double_t TSolROOTFile::GetGenVz(Int_t i)
{
    return (fGenData[i]->GetV()).Z()*1e-3;
}
//_______________________________________________________________________
Double_t TSolROOTFile::GetGenTheta(Int_t i)
{
    return (fGenData[i]->GetP()).Theta()/TMath::Pi()*180.;
}
//_______________________________________________________________________
void TSolROOTFile::ClearVectors()
{
  //generated branch
  gen_pid->clear();
  gen_px->clear();
  gen_py->clear();
  gen_pz->clear();
  gen_vx->clear();
  gen_vy->clear();
  gen_vz->clear();
  
  //solid_gem branch
  solid_gem_id->clear();
  solid_gem_hitn->clear();
  solid_gem_pid->clear();
  solid_gem_trid->clear();
  solid_gem_x->clear();
  solid_gem_y->clear();
  solid_gem_z->clear();
  solid_gem_lxin->clear();
  solid_gem_lyin->clear();
  solid_gem_lzin->clear();
  solid_gem_tin->clear();
  solid_gem_lxout->clear();
  solid_gem_lyout->clear();
  solid_gem_lzout->clear();
  solid_gem_tout->clear();
  solid_gem_px->clear();
  solid_gem_py->clear();
  solid_gem_pz->clear();
  solid_gem_vx->clear();
  solid_gem_vy->clear();
  solid_gem_vz->clear();
  solid_gem_ETot->clear();
  solid_gem_trE->clear();
  solid_gem_weight->clear();

  //flux branch, used to record calorimeter hit position and energy
  //deposition, used later for tracking
  flux_id->clear();
  flux_pid->clear();
  flux_tid->clear();
  flux_trackE->clear();
  flux_avg_x->clear();
  flux_avg_y->clear();
  flux_avg_z->clear();
  
  header_var8->clear();
}
//______________________________________________________________________
void TSolROOTFile::DeleteVectors()
{
	
	vector<int>().swap(*gen_pid);
	
	vector<double>().swap(*gen_px);
  vector<double>().swap(*gen_py);
  vector<double>().swap(*gen_pz);
  vector<double>().swap(*gen_vx);
  vector<double>().swap(*gen_vy);
  vector<double>().swap(*gen_vz);
  
  //solid_gem branch
  vector<double>().swap(*solid_gem_id);
  vector<double>().swap(*solid_gem_hitn);
  vector<double>().swap(*solid_gem_pid);
  vector<double>().swap(*solid_gem_trid);
  vector<double>().swap(*solid_gem_x);
  vector<double>().swap(*solid_gem_y);
  vector<double>().swap(*solid_gem_z);
  vector<double>().swap(*solid_gem_lxin);
  vector<double>().swap(*solid_gem_lyin);
  vector<double>().swap(*solid_gem_lzin);
  vector<double>().swap(*solid_gem_tin);
  vector<double>().swap(*solid_gem_lxout);
  vector<double>().swap(*solid_gem_lyout);
  vector<double>().swap(*solid_gem_lzout);
  vector<double>().swap(*solid_gem_tout);
  vector<double>().swap(*solid_gem_px);
  vector<double>().swap(*solid_gem_py);
  vector<double>().swap(*solid_gem_pz);
  vector<double>().swap(*solid_gem_vx);
  vector<double>().swap(*solid_gem_vy);
  vector<double>().swap(*solid_gem_vz);
  vector<double>().swap(*solid_gem_ETot);
  vector<double>().swap(*solid_gem_trE);
  vector<double>().swap(*solid_gem_weight);

  //flux branch, used to record calorimeter hit position and energy
  //deposition, used later for tracking
  vector<double>().swap(*flux_id);
  vector<double>().swap(*flux_pid);
  vector<double>().swap(*flux_tid);
  vector<double>().swap(*flux_trackE);
  vector<double>().swap(*flux_avg_x);
  vector<double>().swap(*flux_avg_y);
  vector<double>().swap(*flux_avg_z);
  
  vector<double>().swap(*header_var8);
}
//#endif









