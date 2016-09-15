// SoLID tracking simulations: digitize (simulate detector response of) MC raw data
#include <vector>
void Digitize_sidis_6gem_qgsp( const char *filename = NULL,
      const char* outfile = NULL, double bgratio = 0.1,
      int nevents = 10000, int startevent = 0, bool interactive = false )
{
  cout << endl << "** This gets called with 'analyzer' and not 'root' **" << endl;
  cout << "** If you're getting missing symbol errors, this is likely the cause **" << endl << endl;

  const char* data_dir = "/home/xiongw/Research/digitization/sidis_workspace/signal_run/sidis_signal";
  const char* backgr_dir = "/home/xiongw/Research/digitization/sidis_workspace/background_run/sidis_bg";
  const char* name_suffix = "digitized.root"; //"s10k_b21M_(25%).root";

  const int mindata = startevent;
  const int maxdata = startevent+nevents;
  // 25.74M for 275 ns * 15 uA, (73000.0/1e8) is a normalization factor because I remove all the
  // empty entries in the 1e8 background file

  //    const int nbacktoadd = 100000;
  const bool map_sectors = false; // we don't map sectors for SIDIS
  const int mark_interval = 100;
  int nbacktoadd = (73000.0/1e8) * 25740000 * bgratio;//275ns
  gSystem->Load("libsolgem.so");

  if(filename == NULL){
    filename = "disrate";
  }

  string fnam = data_dir;
  fnam += "/";
  fnam += filename;
  if( fnam.substr(fnam.length()-5,5) != ".root" )//using root file for the input, for now
    fnam += ".root";
  if( gSystem->AccessPathName(fnam.c_str(),kReadPermission) ) {
    cout << "Error: Can't read input file " << fnam << endl;
    return;
  }
  cout << "Using file " << fnam << endl;

  ////////////////////////////////////////////////////////////////

  TSolGEMChamber *ddy;
  TSolSpec *dds;
  TSolSimGEMDigitization *ddd;
  const int NSECTORS = 30, NPLANES = 6;

  dds = new TSolSpec ("gemc", "SOLiD spectrometer");
  if( dds->Init() )
    return;

  for( int i = 0; i < NSECTORS*NPLANES; i++ ){
    ddy = new TSolGEMChamber (Form("gem%d", i+1),"GEM chamber");
    ddy->SetApparatus(dds);
    if( ddy->Init() )
      return;
    dds->AddGEM (ddy);
  }

  //dds->Print();
  
  ddd = new TSolSimGEMDigitization (*dds, "ratedig");
  ddd->SetMapSector(map_sectors);
  
  ddd->Print();
  
  ////////////////////////////////////////////////////////////////

  // Signal data file
  TSolROOTFile* f = new TSolROOTFile(fnam.c_str());
  cout<<fnam.c_str()<<endl;
  cout << "ROOT file name = " << f->GetFileName() << endl;

  int res = f->Open();

  if( res != 1 ){
    cout << "Error " << res << " when opening ROOT file" << endl;
    return;
  }

  //setting signal particle information, need to match with simulation file
  f->AddSignalInfo(11, 1); //first track is a electron track
  //f->AddSignalInfo(-211, 2); //second track is pion track 

  // Background data, will modify this later
  vector<string> backgr_files;
  void* dirp = gSystem->OpenDirectory(backgr_dir);
  if( !dirp ) {
    cerr << "Error: Can't open background run directory " << backgr_dir << endl;
    return;
  }
  cout << "Scanning for background files in " << backgr_dir << " " << flush;
  const char* dir_item;
  while( (dir_item = gSystem->GetDirEntry(dirp)) ) {
    if( !strncmp("reduced_background_solid_SIDIS_He3_03122016",dir_item,43) &&
        !strncmp(".root",dir_item+strlen(dir_item)-5,5) ) {
      fnam = backgr_dir;
      fnam.append("/");
      fnam.append(dir_item);
      backgr_files.push_back(fnam);
      cout << "." << flush;
    }
  }
  cout << endl << "Found " << backgr_files.size() << " background runs" << endl;
  gSystem->FreeDirectory(dirp);
  
  TSolROOTFile* fback = 0;
  ////////////////////////////////////////////////////////////////

  TSolGEMData *gd, *gb;
  gd = new TSolGEMData();
  gb = new TSolGEMData();

  const char* ofname = Form("%s_%s", (outfile != 0) ? outfile : filename, name_suffix);
  cout << "Writing digitized output to " << ofname << endl;
  ddd->InitTree (*dds, ofname);
  
  cout << "Digitizing events " << mindata << "-" << maxdata
       << " (nevents = " << maxdata-mindata << ")" << endl;
  cout << "Background ratio = " << 100*bgratio << "% ("
       << nbacktoadd/1e3 << "k bg/signal events)" << endl;
  
  int ndata = 0, nwritten = 0, c = 0;
  int bg_mark_interval = 1;
  if( nbacktoadd > 0 ) {
    if( map_sectors )//we don't map sectors for SIDIS
      nbacktoadd /= NSECTORS;
    int round_level = TMath::Power(10,TMath::FloorNint(TMath::Log10(nbacktoadd/10)));
    bg_mark_interval = int(nbacktoadd/10/round_level) * round_level;
  }
  while( ndata < maxdata && f->ReadNextEvent() ){
    ndata++;
    if( ndata < mindata )
      continue;
    if( interactive || ndata % mark_interval == 0 ){
      cout << "Event " << ndata << endl;
    }
    f->GetGEMData(gd);
    
    //if( gd->GetNHit() <= 4 )//don't do anything if the signal particle doesn't hit any GEM
      //continue;
    //if ( f->GetFillBitsEC() != 1<<1) continue; //LA
    //if ( f->GetFillBitsEC() != 1) continue; //FA
    //if ( f->GetECEDep() < 0.9) continue;
    bool pass = true;
    for (unsigned int ii=0; ii<f->GetNSignal(); ii++){
      //condition for event selection
      if (f->GetSigECBit(ii) != 1) pass = false;
      if (f->GetSigECEDep(ii) < 0.9) pass = false;
    }

    if (!pass) continue;
    // particle did not hit neither LAEC nor FAEC 
    ddd->SetTreeEvent(*gd,*f,nwritten+1);
    ddd->Digitize(*gd, *dds);
    // Add background events
    int nbackadded = 0;
    bool success = false;
    while( nbackadded < nbacktoadd ) {
      if( !fback ) {
        int backidx = rand() % backgr_files.size();
        fback = new TSolROOTFile(backgr_files[backidx].c_str());
        fback->EnableRotateBg();
        fback->Open();
        fback->SetSource( backidx+1 );
        success = false;
      }
      if( fback->ReadNextEvent() == 0 ) {
        fback->Close();
        delete fback; fback = 0;
	//success = false;
        if( success )  // Prevent infinite loop with bad input
          continue;
        cerr << "Error reading background ROOT file " << backgr_files[backidx] << endl;
        return;
      }
      success = true;
      fback->GetGEMData(gb);

      ddd->AdditiveDigitize(*gb,*dds);
      nbackadded++;
      if( nbackadded % bg_mark_interval == 0 ){
        //cout << "Background event ";
        //if( bg_mark_interval >= 1000 )
          //cout << nbackadded/1000 << "k";
        //else
          //cout << nbackadded;
	  //cout << endl;
      }
    }

    //delete the background file
    if (fback != NULL){
     fback->Close();
     delete fback; fback = 0;
    }

    if( ddd->GetEvent()->GetNclust() > 0 ) {
      nwritten++;
      ddd->FillTree();
    }

    if( interactive ) {
      ddd->GetEvent()->Print("all");
      // Enter -1 to quit or # of events to process until next stop
      if( c == 0 )
        cin >> c;
      if( c < 0 )
        break;
      --c;
    }

    
  }

  cout << "Read " << ndata << " events, " << nwritten << " written" << endl;

  ddd->WriteTree();
  ddd->CloseTree();
  
  if( fback ) fback->Close();
  f->Close();
  
  delete gb;
  delete gd;
  
  exit(0);
}
