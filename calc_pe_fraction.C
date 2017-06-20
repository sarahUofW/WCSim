// Simple example of reading a generated Root file

void calc_pe_fraction(char * infilename){

  // Clear global scope
  gROOT->Reset();

  // Load the library with class dictionary info
  // (create with "gmake shared")
  char* wcsimdirenv = getenv("WCSIMDIR");
  if(wcsimdirenv !=  NULL){
    gSystem->Load("${WCSIMDIR}/libWCSimRoot.so");
  }else{
    gSystem->Load("../WCSim/libWCSimRoot.so");
  }

  // Open the input file
  TFile file(infilename);
  
  // Get the a pointer to the tree from the file
  TTree *etree = (TTree*)file.Get("wcsimT");
  TTree *gtree = (TTree*)file.Get("wcsimGeoT");


  WCSimRootEvent * fEv = new WCSimRootEvent();
  WCSimRootGeom * fGeo = new WCSimRootGeom();

  // Set the branch address for reading from the tree
  TBranch *evbranch = etree->GetBranch("wcsimrootevent");
  evbranch->SetAddress(&fEv);
  TBranch *geobranch = gtree->GetBranch("wcsimrootgeom");
  geobranch->SetAddress(&fGeo);
  
  // load geometry
  gtree->GetEntry(0);

  // Open an output file, and make some histograms
  TFile * fout = new TFile("calc_pe_fraction.root","recreate");
  TH1D* htubeid = new TH1D("htubeid","tube id; tube id; count/tube",171,-0.5,170.5);
  htubeid->Sumw2();
  
  int ev;
  int nevent = etree->GetEntries();

  double nckov = 0.0;
  const double nthrown = 1.0;
  
  for (ev=0;ev<nevent; ev++)  {
    if (ev%100==0)std::cout<<ev<<"/"<<nevent<<std::endl;
    etree->GetEntry(ev);

    WCSimRootTrigger * fTrig=fEv->GetTrigger(0);
    std::cout<<"ev "<<ev<<" vtx="<<fTrig->GetVtx(0)<<", "<<fTrig->GetVtx(1)<<", "<<fTrig->GetVtx(2)<<std::endl;
    std::cout<<"  mode="<<fTrig->GetMode()
     	     <<"  NTubes="<<fTrig->GetNumTubesHit()
    	     <<"  Ntrack="<<fTrig->GetNtrack()
     	     <<"  sumQ="<<fTrig->GetSumQ()
    	     <<"  Ncherenkov="<<fTrig->GetNcherenkovhits()
    	     <<"  NDigiHits="<<fTrig->GetNcherenkovdigihits()
    	     <<std::endl;

    if ( fTrig->GetNcherenkovhits() > 0 ){
      //      nckov += fTrig->GetNcherenkovhits();

    
      WCSimRootTrigger * fSub = fEv->GetTrigger(0);
      TClonesArray* hits = fSub->GetCherenkovHits();
      if ( fSub == NULL ) cout<<"No subevent found"<<endl;
      else{
	cout<<"Ncherenkov="<<fSub->GetNcherenkovhits()<<endl;
	for (int idig=0; idig<fSub->GetNcherenkovhits(); idig++){
	  WCSimRootCherenkovHit * hit = (WCSimRootCherenkovHit*)((*hits)[idig]);
	  int tubeid=hit->GetTubeID();
	  WCSimRootPMT pmt = fGeo->GetPMT(tubeid-1);
	  TVector3 pmtpos( pmt.GetPosition(0), 
			   pmt.GetPosition(1),
			   pmt.GetPosition(2) );
	  cout<<"tubeid="<<tubeid<<" pmtpos="<<pmtpos.x()<<", "<<pmtpos.y()<<", "<<pmtpos.z()
	      <<endl;

	  htubeid->Fill( tubeid );
	  
	  if ( tubeid==38 )  nckov += fTrig->GetNcherenkovhits();

	}
      }
    }
    
  }

  std::cout<<"Fraction of photons reaching photocathode tube 38 is "
	   << nckov/nthrown/nevent << " +- " << sqrt( double(nckov) )/nevent/nthrown
	   << std::endl;

  htubeid->Scale( 1.0 / nthrown / nevent );

  fout->Write();
  fout->Close();
}
