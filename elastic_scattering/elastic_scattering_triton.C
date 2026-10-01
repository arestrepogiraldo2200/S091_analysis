
char fname[200];
int tpat;

TH1D* h_CALIFATHETASUM;
TH1D* h_CALIFAESUM;
TH2D* h_CALIFATHETA;
TH2D* h_CALIFAE;
TH2D* h_CALIFAETHETASUMS;

TH2D* h_NEULANDET;
TH2D* h_NEULANDPADDLET;

TH2D* h_NEULANDET_FILTER;
TH2D* h_NEULANDPADDLET_FILTER;

TH2D* h_NEULANDET_MISS;
TH2D* h_NEULANDPADDLET_MISS;

void elastic_scattering_triton(){

	h_CALIFATHETASUM = new TH1D("Opening angle","Opening angle",100,0,2);
    h_CALIFATHETASUM->GetXaxis()->SetTitle("Opening angle");

	h_CALIFAESUM = new TH1D("Total energy","Total energy",100,0,600000);
    h_CALIFAESUM->GetXaxis()->SetTitle("Total energy");

    h_CALIFATHETA = new TH2D("Angular correlation","Angular correlation",100,0,2,100,0,2);
    h_CALIFATHETA->GetXaxis()->SetTitle("Angle 1");
    h_CALIFATHETA->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAE = new TH2D("Energy correlation","Energy correlation",100,0,600000,100,0,600000);
    h_CALIFAE->GetXaxis()->SetTitle("E1");
    h_CALIFAE->GetYaxis()->SetTitle("E2");

    h_CALIFAETHETASUMS = new TH2D("Energy-Angle correlation","Energy-Angle correlation",100,0,2,100,0,600000);
    h_CALIFAETHETASUMS->GetXaxis()->SetTitle("Opening angle");
    h_CALIFAETHETASUMS->GetYaxis()->SetTitle("E1+E2");

	h_NEULANDET = new TH2D("Time-Energy","Time-energy",500,0,200,200,-6000,6000);
    h_NEULANDET->GetXaxis()->SetTitle("Energy");
    h_NEULANDET->GetYaxis()->SetTitle("Time");

	h_NEULANDPADDLET = new TH2D("Time-paddle","Time-paddle",1301,0,1300,200,-4000,4000);
    h_NEULANDPADDLET->GetXaxis()->SetTitle("Paddle");
    h_NEULANDPADDLET->GetYaxis()->SetTitle("Time");

	h_NEULANDET_FILTER = new TH2D("Time-Energy two hits in califa","Time-energy two hits in califa",500,0,200,200,-6000,6000);
    h_NEULANDET_FILTER->GetXaxis()->SetTitle("Energy");
    h_NEULANDET_FILTER->GetYaxis()->SetTitle("Time");

	h_NEULANDPADDLET_FILTER = new TH2D("Time-paddle two hits in califa","Time-paddle two hits in califa",1301,0,1300,200,-4000,4000);
    h_NEULANDPADDLET_FILTER->GetXaxis()->SetTitle("Paddle");
    h_NEULANDPADDLET_FILTER->GetYaxis()->SetTitle("Time");

	h_NEULANDET_MISS = new TH2D("Time-Energy zero hits in califa - Miss","Time-energy zero hits in califa - Miss",500,0,200,200,-6000,6000);
    h_NEULANDET_MISS->GetXaxis()->SetTitle("Energy");
    h_NEULANDET_MISS->GetYaxis()->SetTitle("Time");

	h_NEULANDPADDLET_MISS = new TH2D("Time-paddle zero hits in califa - Miss","Time-paddle zero hits in califa - Miss",1301,0,1300,200,-4000,4000);
    h_NEULANDPADDLET_MISS->GetXaxis()->SetTitle("Paddle");
    h_NEULANDPADDLET_MISS->GetYaxis()->SetTitle("Time");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 
    TChain* chain = new TChain("evt");   
    sprintf(fname, "./main0067_0001_unpacked_califa_ToFD_LOS_NEULAND.root");
    //sprintf(fname, "./main0083_0001_unpacked_califa_ToFD_LOS_NEULAND.root");

    chain->Reset();
    chain->Add(fname); 
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 

	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("EventHeader.");
    branchData->SetAddress(&header);    
    
	TClonesArray*  hitLOS = new TClonesArray("R3BLosHitData",5); 
	R3BLosHitData** hitLos;                                      
	TBranch *branchHitLOS = chain->GetBranch("LosHit");         
	branchHitLOS->SetAddress(&hitLOS);                           
	Long64_t LosHitPerEvent = 0;                               

	TClonesArray*  hitTOFD = new TClonesArray("R3BTofdHitData",5);
	R3BTofdHitData** hitTofd;
	TBranch *branchHitTOFD = chain->GetBranch("TofdHit");
	branchHitTOFD->SetAddress(&hitTOFD);
	Long64_t TofdHitPerEvent = 0;  	

    TClonesArray*  hitCALIFA = new TClonesArray("R3BCalifaClusterData",5);
	R3BCalifaClusterData** clusterCalifa;
	TBranch *branchHitCALIFA = chain->GetBranch("CalifaClusterData");
	branchHitCALIFA->SetAddress(&hitCALIFA);
	Long64_t CalifaHitPerEvent = 0;  	

	TClonesArray* hitNEULAND = new TClonesArray("R3BNeulandHit",5);
	R3BNeulandHit** hitNeuland;
	TBranch *branchHitNEULAND = chain->GetBranch("NeulandHits");
	branchHitNEULAND->SetAddress(&hitNEULAND);
	Long64_t NeulandHitPerEvent = 0;  	

    /////////////////////////////////////////////////////////////////////////////////////////

	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in " << fname << " : " << nevent << endl;

	for(Long64_t e=0; e<nevent; e++){   
    
		hitLOS->Clear();
	    hitTOFD->Clear();
        hitCALIFA->Clear();
		hitNEULAND->Clear();
	    chain->GetEvent(e);
	    
        // Entries of LOS and TOFD per event in the chain of trees
        tpat            = header->GetTpat();
	    LosHitPerEvent  = hitLOS->GetEntriesFast();
	    TofdHitPerEvent = hitTOFD->GetEntriesFast();
        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
		NeulandHitPerEvent = hitNEULAND->GetEntriesFast();

		if( TofdHitPerEvent == 0 ){

            bool califaL=false;
            bool califaR=false;
            double thetaTL, thetaTR, phiTR, phiTL, EL, ER;

			/////////////////////////////////////////////////////////////////////////////////////////
			// Two hits in califa and two hits in neuland

			if( CalifaHitPerEvent == 2 && NeulandHitPerEvent > 0){

				//cout << "Two hits in califa, none in tofd" << endl;
                clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];
                double Esum     = 0;
				double thetasum_plane = 0;


				for(int p = 0; p < CalifaHitPerEvent; p++){

					clusterCalifa[p] = (R3BCalifaClusterData*)hitCALIFA->At(p);
					double E     = clusterCalifa[p]->GetEnergy();
					double theta = clusterCalifa[p]->GetTheta();
					double phi   = clusterCalifa[p]->GetPhi();

					if ( E > 20000 ){

						if ( (phi < 1.5 && phi > -1.5) && califaL == false){

							califaL = true;
							thetaTL = theta;
							EL = E;
							thetasum_plane += thetaTL;
							phiTL = phi;
							Esum += E;

							if (califaR == true && califaL == true){
								if ( abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
									h_CALIFATHETA->Fill(thetaTL,thetaTR);
									h_CALIFATHETASUM->Fill(thetasum_plane);
									h_CALIFAESUM->Fill(Esum);
									//cout << EL << "  " << ER << "  " << Esum << endl;
									h_CALIFAE->Fill(EL,ER);
									h_CALIFAETHETASUMS->Fill(thetasum_plane, Esum);
								}
								califaL = false;
								califaR = false;
								thetaTL = 0;
								thetaTR = 0;
								EL = 0;
								ER = 0;
							}
							continue;
						}

						if ( (phi>1.5 || phi<-1.5) && califaR == false) { 
													
							califaR = true;
							thetaTR = theta;
							ER = E;
							thetasum_plane += thetaTR;
							phiTR = phi;
							Esum += E;

							if (califaR == true && califaL == true){
								if ( abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
									h_CALIFATHETA->Fill(thetaTL,thetaTR);
									h_CALIFAESUM->Fill(Esum);
									//cout << EL << "  " << ER << "  " << Esum << endl;
									h_CALIFATHETASUM->Fill(thetasum_plane);
									h_CALIFAE->Fill(EL,ER);
									h_CALIFAETHETASUMS->Fill(thetasum_plane, Esum);
								}
								califaL = false;
								califaR = false;
								thetaTL = 0;
								thetaTR = 0;
								EL = 0;
								ER = 0;
							}
							continue;
						}

						if ((califaR == true && califaL == false) || (califaR == false && califaL == true)){
							califaL = false;
							califaR = false;
						}
					}
				}


				/////////////////////////////////////////////////////////////////////////////////////////
				// Neuland detection with two hits in califa
				if( NeulandHitPerEvent > 0 ){

					hitNeuland = new R3BNeulandHit*[NeulandHitPerEvent];

					for(int p = 0; p < NeulandHitPerEvent; p++){
						
						hitNeuland[p] = (R3BNeulandHit*)hitNEULAND->At(p);

						int paddle = hitNeuland[p]->GetPaddle();
						double T   = hitNeuland[p]->GetT();         // This is the time between the LOS and Neuland
						double E   = hitNeuland[p]->GetE();

						// Histograms without any filter 
						h_NEULANDET_FILTER->Fill(E,T);
						h_NEULANDPADDLET_FILTER->Fill(paddle,T);	
					
					}
				}
			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// Whole signal in Neuland
		if( NeulandHitPerEvent > 0 ){

            hitNeuland = new R3BNeulandHit*[NeulandHitPerEvent];

			for(int p = 0; p < NeulandHitPerEvent; p++){
				
				hitNeuland[p] = (R3BNeulandHit*)hitNEULAND->At(p);

				int paddle = hitNeuland[p]->GetPaddle();
				double T   = hitNeuland[p]->GetT();         // This is the time between the LOS and Neuland
				double E   = hitNeuland[p]->GetE();

				// Histograms without any filter 
				h_NEULANDET->Fill(E,T);
				h_NEULANDPADDLET->Fill(paddle,T);	

			}
		}

		/////////////////////////////////////////////////////////////////////////////////////////
		// Neuland events of miss
		if( TofdHitPerEvent == 1 ){       // All charges are 1

			if( CalifaHitPerEvent == 0 ){

				if( NeulandHitPerEvent > 0 ){

					hitNeuland = new R3BNeulandHit*[NeulandHitPerEvent];

					for(int p = 0; p < NeulandHitPerEvent; p++){
						
						hitNeuland[p] = (R3BNeulandHit*)hitNEULAND->At(p);

						int paddle = hitNeuland[p]->GetPaddle();
						double T   = hitNeuland[p]->GetT();         // This is the time between the LOS and Neuland
						double E   = hitNeuland[p]->GetE();

						// Histograms without any filter 
						h_NEULANDET_MISS->Fill(E,T);
						h_NEULANDPADDLET_MISS->Fill(paddle,T);	

					}
				}			
			}
		}
	}   

	TCanvas* c_CALIFATHETASUM = new TCanvas();
	c_CALIFATHETASUM->cd(1);
	h_CALIFATHETASUM->Draw("hist");

	TCanvas* c_CALIFATHETA = new TCanvas();
	c_CALIFATHETA->cd(1);
	h_CALIFATHETA->Draw("colz");

	TCanvas* c_CALIFAESUM = new TCanvas();
	c_CALIFAESUM->cd(1);
	h_CALIFAESUM->Draw("hist");

	TCanvas* c_CALIFAE = new TCanvas();
	c_CALIFAE->cd(1);
	h_CALIFAE->Draw("colz");

	TCanvas* c_CALIFAETHETASUMS = new TCanvas();
	c_CALIFAETHETASUMS->cd(1);
	h_CALIFAETHETASUMS->Draw("colz");

	TCanvas* c_NEULANDET = new TCanvas();
	c_NEULANDET->cd(1);
	h_NEULANDET->Draw("colz");

	TCanvas* c_NEULANDPADDLET = new TCanvas();
	c_NEULANDPADDLET->cd(1);
	h_NEULANDPADDLET->Draw("colz");

	TCanvas* c_NEULANDET_FILTER = new TCanvas();
	c_NEULANDET_FILTER->cd(1);
	h_NEULANDET_FILTER->Draw("colz");

	TCanvas* c_NEULANDPADDLET_FILTER = new TCanvas();
	c_NEULANDPADDLET_FILTER->cd(1);
	h_NEULANDPADDLET_FILTER->Draw("colz");

	TCanvas* c_NEULANDET_MISS = new TCanvas();
	c_NEULANDET_MISS->cd(1);
	h_NEULANDET_MISS->Draw("colz");

	TCanvas* c_NEULANDPADDLET_MISS = new TCanvas();
	c_NEULANDPADDLET_MISS->cd(1);
	h_NEULANDPADDLET_MISS->Draw("colz");

}