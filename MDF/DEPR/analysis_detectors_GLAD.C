/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname[200];
char histname[200];
int tpat;

TH2D* h_CALIFA;
TH1D* h_CALIFAE;
TH2D* h_CALIFAECORR;

TH1D* h_CALIFAA;
TH2D* h_CALIFATHETA;
TH1D* h_CALIFATHETASUM;
TH2D* h_CALIFAPHI;

TH1D* h_MWPC;
TH1D* h_FIBER30;
TH2D* h_MWPC_FIBERS_X;
TH2D* h_MWPC_TOFD_X;
TH1D* h_TOFDX;

vector<double> Eloss1;
vector<double> Eloss2;
vector<double> CalifaR;
vector<double> CalifaL;

void analysis_detectors_GLAD(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Set histograms canvas and axis titles here  

    h_CALIFAE = new TH1D("Energy","Energy",1000,0,900000);
    h_CALIFAE->GetXaxis()->SetTitle("Energy");

    h_CALIFAA = new TH1D("Azimutal angle","Azimutal angle",100,-4,4);
    h_CALIFAA->GetXaxis()->SetTitle("Polar angle");

    h_CALIFATHETA = new TH2D("Angular correlation","Angular correlation",100,0,2,100,0,2);
    h_CALIFATHETA->GetXaxis()->SetTitle("Angle 1");
    h_CALIFATHETA->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAPHI = new TH2D("Azimutal angular correlation","Azimutal angular correlation",100,-3.2,3.2,100,-3.2,3.2);
    h_CALIFAPHI->GetXaxis()->SetTitle("Angle 1");
    h_CALIFAPHI->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAECORR = new TH2D("Energy correlations","Energy correlations",1000,0,900000,1000,0,900000);
    h_CALIFAECORR->GetXaxis()->SetTitle("Energy 1");
    h_CALIFAECORR->GetYaxis()->SetTitle("Energy 2");

    h_CALIFATHETASUM = new TH1D("Opening angle","Opening angle",100,0,2);
    h_CALIFATHETASUM->GetXaxis()->SetTitle("Opening angle");

    h_MWPC = new TH1D("MWPC 1 x position","MWPC 1 x position",200,-200,100);
    h_MWPC->GetXaxis()->SetTitle("MWPC 1 x position");

    h_TOFDX = new TH1D("TOFD x position","TOFD x position",200,-80,80);
    h_TOFDX->GetXaxis()->SetTitle("TOFD x position");

    h_FIBER30 = new TH1D("Fibers x position","Fibers x position",512,-25,25);
    h_FIBER30->GetXaxis()->SetTitle("Fibers x position");

    h_MWPC_FIBERS_X = new TH2D("X position correlations MWPC-FIBERS","X position correlations MWPC-FIBERS",200,-200,100,512,-25,25);
    h_MWPC_FIBERS_X->GetXaxis()->SetTitle("MWPC X-position");
    h_MWPC_FIBERS_X->GetYaxis()->SetTitle("FIBERS X-position");

    h_MWPC_TOFD_X = new TH2D("X position correlations MWPC-TOFD","X position correlations MWPC-TOFD",200,-200,100,200,-100,100);
    h_MWPC_TOFD_X->GetXaxis()->SetTitle("MWPC X-position");
    h_MWPC_TOFD_X->GetYaxis()->SetTitle("TOFD X-position");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt");  
    sprintf(fname, "./main0252_0001_unpacked_MDF.root");
    //sprintf(fname, "./main0321_0001_unpacked_MDF.root"); 
    //sprintf(fname, "./main0327_0001_unpacked_MDF.root");

    chain->Reset();
    chain->Add(fname); 
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 

	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("EventHeader.");
    branchData->SetAddress(&header);    

	TClonesArray* hitTOFD = new TClonesArray("R3BTofdHitData",5);
	R3BTofdHitData** hitTofd;
	TBranch *branchHitTOFD = chain->GetBranch("TofdHit");
	branchHitTOFD->SetAddress(&hitTOFD);
	Long64_t TofdHitPerEvent = 0;  	

    TClonesArray* hitCALIFA = new TClonesArray("R3BCalifaClusterData",5);
	R3BCalifaClusterData** clusterCalifa;
	TBranch *branchHitCALIFA = chain->GetBranch("CalifaClusterData");
	branchHitCALIFA->SetAddress(&hitCALIFA);
	Long64_t CalifaHitPerEvent = 0;  	

    TClonesArray* hitMWPC = new TClonesArray("R3BMwpcHitData",5);
	R3BMwpcHitData** hitMwpc;
	TBranch *branchHitMWPC = chain->GetBranch("Mwpc1HitData");
	branchHitMWPC->SetAddress(&hitMWPC);
	Long64_t MwpcHitPerEvent = 0;  	

    TClonesArray* hitFI30 = new TClonesArray("R3BFiberMAPMTHitData",5);
	R3BFiberMAPMTHitData** hitFi30;
	TBranch *branchHitFI30 = chain->GetBranch("Fi32Hit");
	branchHitFI30->SetAddress(&hitFI30);
	Long64_t Fi30HitPerEvent = 0;  	

    /////////////////////////////////////////////////////////////////////////////////////////
	
	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in " << fname << " : " << nevent << endl;

	for(Long64_t e=0; e<nevent; e++)
	{   
	    hitTOFD->Clear();
        hitCALIFA->Clear();
        hitMWPC->Clear();
        hitFI30->Clear();
	    chain->GetEvent(e);
	    
        // Entries of LOS and TOFD per event in the chain of trees
        tpat            = header->GetTpat();
	    TofdHitPerEvent = hitTOFD->GetEntriesFast();
        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        MwpcHitPerEvent = hitMWPC->GetEntriesFast();
        Fi30HitPerEvent = hitFI30->GetEntriesFast();

        if ((tpat & 8) == 8){
            
            /////////////////////////////////////////////////////////////////////////////////////////

            if( TofdHitPerEvent > 0 ){

                hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
                
                for(int p = 0; p < TofdHitPerEvent; p++){

                    hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
                    double TOFD_X = hitTofd[p]->GetX();
                    h_TOFDX->Fill(TOFD_X);
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////

            bool califaL=false;
            bool califaR=false;
            double thetaTL, thetaTR, phiTR, phiTL, EL, ER;
            double phi1 = 0., phi2 =0., theta1=0., theta2=0.;

            if( CalifaHitPerEvent == 2 ){

                clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];
                hitMwpc       = new R3BMwpcHitData*[MwpcHitPerEvent];
                hitFi30       = new R3BFiberMAPMTHitData*[Fi30HitPerEvent];
                hitTofd       = new R3BTofdHitData*[TofdHitPerEvent];

                double thetasum_plane = 0;

                // ######################### CALIFA #######################################

                for(int p = 0; p < CalifaHitPerEvent; p++){

                    clusterCalifa[p] =(R3BCalifaClusterData*)hitCALIFA->At(p);

                    double E = clusterCalifa[p]->GetEnergy();
                    double theta = clusterCalifa[p]->GetTheta();
                    double phi = clusterCalifa[p]->GetPhi();

                    if(p==0){
                        phi1 = phi;
                        theta1 = theta; 
                    }
                    if(p==1){
                        phi2 = phi;
                        theta2 = theta;
                    }

                    if ( E > 20000 ){
                        
                        h_CALIFAA->Fill(phi);

                        if ( (phi < 1.5 && phi > -1.5) && califaL == false){

                            califaL = true;
                            thetaTL = theta;
                            thetasum_plane += thetaTL;
                            phiTL = phi;
                            EL = E;

                            if (califaR == true && califaL == true){
                                if (  abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
                                    //h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                    //h_CALIFAPHI->Fill(phiTL,phiTR);
                                    h_CALIFATHETASUM->Fill(thetasum_plane);
                                    h_CALIFAECORR->Fill(EL,ER);
                                    h_CALIFAE->Fill(EL+ER);
                                }
                                califaL = false;
                                califaR = false;
                                thetaTL = 0;
                                thetaTR = 0;
                            }
                            continue;
                        }

                        if ( (phi>1.5 || phi<-1.5) && califaR == false) { 
                                                            
                            califaR = true;
                            thetaTR = theta;
                            thetasum_plane += thetaTR;
                            phiTR = phi;
                            ER = E;

                            if (califaR == true && califaL == true){
                                if (  abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
                                    //h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                    //h_CALIFAPHI->Fill(phiTL,phiTR);
                                    h_CALIFATHETASUM->Fill(thetasum_plane);
                                    h_CALIFAECORR->Fill(EL,ER);
                                    h_CALIFAE->Fill(EL+ER);
                                }
                                califaL = false;
                                califaR = false;
                                thetaTL = 0;
                                thetaTR = 0;
                            }
                            continue;
                        }

                        if ((califaR == true && califaL == false) || (califaR == false && califaL == true)){
                            califaL = false;
                            califaR = false;
                        }

                    }
                }

                if (abs(abs(phi1)+abs(phi2))- 3.1415926535 < 0.2){
                    h_CALIFAPHI->Fill(phi1,phi2);
                    h_CALIFATHETA->Fill(theta1,theta2);
                }

                // ######################### CALIFA #######################################

                // ############################# MWPC #####################################
                
                for(int p = 0; p < MwpcHitPerEvent; p++){

                    hitMwpc[p] =(R3BMwpcHitData*)hitMWPC->At(p);
                    double x = hitMwpc[p]->GetX();
                    h_MWPC->Fill(x);

                    for(int q = 0; q < Fi30HitPerEvent; q++){
                        hitFi30[q] =(R3BFiberMAPMTHitData*)hitFI30->At(q);
                        double xfib = hitFi30[q]->GetX();
                        h_MWPC_FIBERS_X->Fill(x, xfib); 
                    }
                }

                // ############################# MWPC #####################################

                // ############################# TOFD #####################################

                for(int p = 0; p < MwpcHitPerEvent; p++){

                    hitMwpc[p] =(R3BMwpcHitData*)hitMWPC->At(p);
                    double x = hitMwpc[p]->GetX();
                    h_MWPC->Fill(x);

                    for(int q = 0; q < TofdHitPerEvent; q++){

                        hitTofd[q] =(R3BTofdHitData*)hitTOFD->At(q);
                        double TOFD_X = hitTofd[p]->GetX();
                        h_MWPC_TOFD_X->Fill(x, TOFD_X); 
                    }
                }

                // ############################# TOFD #####################################
                
                // ############################# FIBERS ###################################
                
                for(int p = 0; p < Fi30HitPerEvent; p++){
                    hitFi30[p] =(R3BFiberMAPMTHitData*)hitFI30->At(p);
                    double x = hitFi30[p]->GetX();
                    h_FIBER30->Fill(x);
                }
                
                // ############################# FIBERS ###################################
            }  
        }  
	}

TCanvas* c_CALIFAE = new TCanvas();
c_CALIFAE->cd(1);
h_CALIFAE->Draw("hist");

TCanvas* c_CALIFAECORR = new TCanvas();
c_CALIFAECORR->cd(1);
h_CALIFAECORR->Draw("colz");

TCanvas* c_CALIFATHETASUM = new TCanvas();
c_CALIFATHETASUM->cd(1);
h_CALIFATHETASUM->Draw("hist");

TCanvas* c_CALIFAA = new TCanvas();
c_CALIFAA->cd(1);
h_CALIFAA->Draw("hist");

TCanvas* c_CALIFATHETA = new TCanvas();
c_CALIFATHETA->cd(1);
h_CALIFATHETA->Draw("colz");

TCanvas* c_CALIFAPHI = new TCanvas();
c_CALIFAPHI->cd(1);
h_CALIFAPHI->Draw("colz");

TCanvas* c_MWPC = new TCanvas();
c_MWPC->cd(1);
h_MWPC->Draw("hist");

TCanvas* c_TOFDX = new TCanvas();
c_TOFDX->cd(1);
h_TOFDX->Draw("hist");

TCanvas* c_FIBER30 = new TCanvas();
c_FIBER30->cd(1);
h_FIBER30->Draw("hist");

TCanvas* c_MWPC_FIBERS_X = new TCanvas();
c_MWPC_FIBERS_X->cd(1);
h_MWPC_FIBERS_X->Draw("colz");

TCanvas* c_MWPC_TOFD_X = new TCanvas();
c_MWPC_TOFD_X->cd(1);
h_MWPC_TOFD_X->Draw("colz");

}