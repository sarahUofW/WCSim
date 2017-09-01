// Simple example of reading a generated Root file

void checkwcsim(char * infilename){

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

  // default units are not too obvious!
  // guessed below based on what was put into simulation
  printf("Cyl radius %f cm\n", fGeo->GetWCCylRadius());
  printf("Cyl length %f cm\n", fGeo->GetWCCylLength());
  printf("PMT radius %f m\n", fGeo->GetWCPMTRadius());
  printf("Offset x y z %f %f %f cm\n", fGeo->GetWCOffset(0),
	 fGeo->GetWCOffset(1),fGeo->GetWCOffset(2));

  TVector3 detcenter( fGeo->GetWCOffset(0), fGeo->GetWCOffset(1), fGeo->GetWCOffset(2) );
  
  // Open an output file, and make some histograms
  TFile * fout = new TFile("checkwcsim.root","recreate");
  fout->cd();
  TH2D*hvertex = new TH2D("hvertex","Vertex z vs r; r (cm); z (cm)",100,0.,1700.,100,-3620.,3620.);
  TH1D*hnsubev = new TH1D("hnsubev","Number Of Sub Events",20,-0.5,19.5);
  TH1D*hntrack = new TH1D("hntrack","Ntrack",20,-0.5,19.5);
  TH1D*hntubes = new TH1D("hntubes","Num Tubes Hit (SubEv 0 only)",100,0.,10000.);
  TH1D*hnchdig = new TH1D("hnchdig","Ncherenkovdigihits (SubEv0 only)",100,0.,10000.);
  TH2D*hnchdig2 = new TH2D("hnchdig2","Ncherenkovdigihits SubEv1 vs 0",100,0.,10000.,100,0.,10000.);
  TH1D*hsumq   = new TH1D("hsumq","SumQ (SubEv 0 only)",100,0.,20000.);
  TH2D*hsumq2  = new TH2D("hsumq2","SumQ SubEv 1 vs 0",100,0.,20000.,100,0.,20000.);

  TH2D*hpmtlocxy = new TH2D("hpmtlocxy","PMT locations xy; X (cm); Y (cm)",
			    100,-350.,350.,
			    100,-1600.,-300.);
  TH2D*hpmtlocxz = new TH2D("hpmtlocxz","PMT locations xz; X (cm); Z (cm)",
			    100,-400.,400.,
			    100,-400.,-400.);
  TH2D*hpmtlocrz = new TH2D("hpmtlocrz","PMT locations yz; Y (cm); Z (cm)",
			    100, -1600., -300.,
			    100, -400. ,  400.);
  TH2D*hpmthitxy = new TH2D("hpmthitxy","PMT hitcount xy; X (cm); Y (cm)",
			    100,-350., 350.,
			    100,-1600., -300.);
  TH2D*hpmthitxz = new TH2D("hpmthitxz","PMT hitcount xy; X (cm); Z (cm)",
			    100,-400., 400.,
			    100,-400., -400.);
  TH2D*hpmthitrz = new TH2D("hpmthitrz","PMT hitcount yz; Y (cm); Z (cm)",
			    100, -1600., -300.,
			    100,-400., 400.);

  //loop over pmts
  for (int ipmt=0; ipmt<fGeo->GetWCNumPMT(); ipmt++){
    WCSimRootPMT pmt = fGeo->GetPMT(ipmt);
    TVector3 pmtpos( pmt.GetPosition(0), 
		     pmt.GetPosition(1),
		     pmt.GetPosition(2) );
    hpmtlocxy->Fill( pmtpos.X(), pmtpos.Y() );
    hpmtlocxz->Fill( pmtpos.X(), pmtpos.Z() );
    hpmtlocrz->Fill( pmtpos.Y(), pmtpos.Z() );
  }
    
  TH2D*hpmthitsum = new TH2D("hpmthitsum","PMT Total Hit Count;Theta(rad);Phi(rad)", 75,0.,TMath::Pi(),150,-TMath::Pi(),TMath::Pi());
  TH2D*hpmtqsum = new TH2D("hpmtqsum","PMT Total Charge Deposit;Theta(rad);Phi(rad)", 75,0.,TMath::Pi(),150,-TMath::Pi(),TMath::Pi());

  TH2D**hthetaphi;
  TDirectory* adirevts = fout->mkdir("ThetaPhiByEvent");
  adirevts->cd();

  // Now loop over "events" 
  int ev;
  int nevent = etree->GetEntries();

  hthetaphi = new TH2D*[nevent];
  TH2D** hthetaphit = new TH2D*[nevent];
  TH2D** hzxt = new TH2D*[nevent];
 
  char hhname[100];
  for (ev=0;ev<nevent; ev++)  {
    if (ev%100==0)std::cout<<ev<<"/"<<nevent<<std::endl;
    sprintf(hhname,"hthetaphi_q_%05d;Theta(rad);Phi(rad)",ev);
    hthetaphi[ev] = new TH2D(hhname,hhname,75,0.,TMath::Pi(),150,-TMath::Pi(),TMath::Pi());
    sprintf(hhname,"hthetaphi_time_%05d;Theta(rad);Phi(rad)",ev);
    hthetaphit[ev] = new TH2D(hhname,hhname,75,0.,TMath::Pi(),150,-TMath::Pi(),TMath::Pi());
    sprintf(hhname,"hzxq_%05d;X(cm);Z(cm)",ev);
    hzxt[ev] = new TH2D( hhname, hhname,
			    100,-350.,350.,
			    100,-350.,350.);    
    etree->GetEntry(ev);
    // fEv
    WCSimRootTrigger * fTrig=fEv->GetTrigger(0);
    std::cout<<"ev "<<ev<<" vtx="<<fTrig->GetVtx(0)<<", "<<fTrig->GetVtx(1)<<", "<<fTrig->GetVtx(2)<<std::endl;
    std::cout<<"  mode="<<fTrig->GetMode()
    	     <<"  NTubes="<<fTrig->GetNumTubesHit()
    	     <<"  Ntrack="<<fTrig->GetNtrack()
    	     <<"  sumQ="<<fTrig->GetSumQ()
	     <<"  Ncherenkov="<<fTrig->GetNcherenkovhits()
    	     <<"  NDigiHits="<<fTrig->GetNcherenkovdigihits()
    	     <<std::endl;
    hvertex->Fill( sqrt( fTrig->GetVtx(0)*fTrig->GetVtx(0) + 
			 fTrig->GetVtx(1)*fTrig->GetVtx(1) ), 
		   fTrig->GetVtx(2) );

    int nsubev = fEv->GetNumberOfSubEvents(); 
    hnsubev->Fill( float( nsubev ) );
    if (fTrig->GetNtrack()>0 || nsubev==0){
      hntrack->Fill( float( fTrig->GetNtrack() ) );
      hntubes->Fill( float(fTrig->GetNumTubesHit()) );
      hnchdig->Fill( float(fTrig->GetNcherenkovdigihits() ) );
      hsumq->Fill( float(fTrig->GetSumQ() ) );
    }else {
      WCSimRootTrigger * fTrig1 = fEv->GetTrigger(1);
      hntrack->Fill( float( fTrig1->GetNtrack() ) );
      hntubes->Fill( float(fTrig1->GetNumTubesHit()) );
      hnchdig->Fill( float(fTrig1->GetNcherenkovdigihits() ) );
      hsumq->Fill( float(fTrig1->GetSumQ() ) );
    }
      
    if ( nsubev >= 1 ){
      WCSimRootTrigger * fTrig1 = fEv->GetTrigger(1);
      hnchdig2->Fill( float(fTrig->GetNcherenkovdigihits()),float( fTrig1->GetNcherenkovdigihits() ) );
      hsumq2->Fill( fTrig->GetSumQ(), fTrig1->GetSumQ() );
    }

    for (int isub=0; isub<=TMath::Max(0,fEv->GetNumberOfSubEvents()); isub++){
      WCSimRootTrigger * fSub = fEv->GetTrigger(isub);
      //std::cout<<"      subev="<<isub
      //	       <<"  mode="<<fSub->GetMode()
      //	       <<"  NTubes="<<fSub->GetNumTubesHit()
      //	       <<"  Ntrack="<<fSub->GetNtrack()
      //	       <<"  sumQ="<<fSub->GetSumQ()
      //	       <<"  NDigiHits="<<fSub->GetNcherenkovdigihits()
      //	       <<std::endl;
      
      TClonesArray* digihits = fSub->GetCherenkovDigiHits();
      for (int idig=0; idig<fSub->GetNcherenkovdigihits(); idig++){
	WCSimRootCherenkovDigiHit * hit = (WCSimRootCherenkovDigiHit*)((*digihits)[idig]);
	int tubeid=hit->GetTubeId();
	float tubeq = hit->GetQ();
	WCSimRootPMT pmt = fGeo->GetPMT(tubeid-1);
	//std::cout<<"idig="<<idig<<" tubeid="<<tubeid<<" tubeid="<<pmt.GetTubeNo()<<std::endl;
	TVector3 pmtpos( pmt.GetPosition(0), 
			 pmt.GetPosition(1),
			 pmt.GetPosition(2) );
	TVector3 ppos = pmtpos - detcenter;
	cout<<"pmtpos="<<pmtpos.x()<<", "<<pmtpos.y()<<", "<<pmtpos.z()<<" Rxz="<<sqrt(pmtpos.x()*pmtpos.x()+pmtpos.z()*pmtpos.z())<<std::endl;
	cout<<" detce="<<detcenter.x()<<", "<<detcenter.y()<<", "<<detcenter.z()<<std::endl;
	cout<<" ppos ="<<ppos.x()<<", "<<ppos.y()<<", "<<ppos.z()<<std::endl;
	// define theta in zx plane from z axis:
	double phi = atan2( ppos.X(), ppos.Z() );
	double theta = acos( ppos.Y() / ppos.Mag() ); 

	float tubet = hit->GetT();
	hthetaphi[ev]->Fill( theta, phi, tubeq );
	hthetaphit[ev]->Fill( theta, phi, tubet );
	hzxt[ev]->Fill( pmtpos.X(), pmtpos.Z(), tubet );
	
	hpmthitsum->Fill( theta, phi );
	hpmtqsum->Fill( theta, phi, tubeq );

	hpmthitxy->Fill( pmtpos.X(), pmtpos.Y() );
	hpmthitxz->Fill( pmtpos.X(), pmtpos.Z() );
	hpmthitrz->Fill( pmtpos.Y(), pmtpos.Z() );
	//hthetaphi[ev]->Fill( pmtpos.Theta(), pmtpos.Phi(), 1.0 );
	
      }
    }


  }

  fout->Write();
  fout->Close();

  

}
