/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname[200];
char histname[200];
int tpat;

TH1D* h_LOSQ;
TH2D* h_LOSBeam;
TH1D* h_TOFDQ[4];
TH2D* h_TOFDCC;
TH2D* h_CALIFA;
TH1D* h_CALIFAE;
TH1D* h_CALIFAA;
TH2D* h_CALIFATHETA;
TH1D* h_CALIFATHETASUM;
TH2D* h_CALIFAPHI;

vector<double> Eloss1;
vector<double> Eloss2;

vector<double> CalifaR;
vector<double> CalifaL;

void sim_analysis_CALIFA(){

    h_CALIFAE = new TH1D("Energy","Energy",100,0,300000);
    h_CALIFAE->GetXaxis()->SetTitle("Energy");

    h_CALIFAA = new TH1D("Azimutal angle","Azimutal angle",100,-4,4);
    h_CALIFAA->GetXaxis()->SetTitle("Polar angle");

    h_CALIFATHETA = new TH2D("Angular correlation","Angular correlation",100,0,2,100,0,2);
    h_CALIFATHETA->GetXaxis()->SetTitle("Angle 1");
    h_CALIFATHETA->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAPHI = new TH2D("Azimutal angular correlation","Azimutal angular correlation",100,-3.2,3.2,100,-3.2,3.2);
    h_CALIFAPHI->GetXaxis()->SetTitle("Angle 1");
    h_CALIFAPHI->GetYaxis()->SetTitle("Angle 2");

    h_CALIFATHETASUM = new TH1D("Opening angle","Opening angle",100,0,2);
    h_CALIFATHETASUM->GetXaxis()->SetTitle("Opening angle");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt");   // ??? Not sure why is this needed
    sprintf(fname, "./sim.root");
    //sprintf(fname, "./main0252_0001_unpacked_025_final.root");
    //sprintf(fname, "./sim.root");

    chain->Reset();
    chain->Add(fname); 
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 

	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("MCEventHeader.");
    branchData->SetAddress(&header);    

    TClonesArray*  hitCALIFA = new TClonesArray("R3BCalifaClusterData",5);
	R3BCalifaClusterData** clusterCalifa;
	TBranch *branchHitCALIFA = chain->GetBranch("CalifaClusterData");
	branchHitCALIFA->SetAddress(&hitCALIFA);
	Long64_t CalifaHitPerEvent = 0;  	

    /////////////////////////////////////////////////////////////////////////////////////////
	
	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in " << fname << " : " << nevent << endl;

	for(Long64_t e=0; e<nevent; e++)
	{   
    
	    // hitLOS->Clear();
	    // hitTOFD->Clear();
        hitCALIFA->Clear();
	    chain->GetEvent(e);
	    
        // Entries of LOS and TOFD per event in the chain of trees
        tpat            = header->GetTpat();
	    // LosHitPerEvent  = hitLOS->GetEntriesFast();
	    // TofdHitPerEvent = hitTOFD->GetEntriesFast();
        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        // cout << CalifaHitPerEvent << "\n";

            /////////////////////////////////////////////////////////////////////////////////////////

            bool califaL=false;
            bool califaR=false;
            double thetaTL, thetaTR, phiTR, phiTL;

                if( CalifaHitPerEvent == 2 ){

                    clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];

                    double thetasum_plane = 0;

                    for(int p = 0; p < CalifaHitPerEvent; p++){

                        clusterCalifa[p] =(R3BCalifaClusterData*)hitCALIFA->At(p);

                        double E = clusterCalifa[p]->GetEnergy();
                        double theta = clusterCalifa[p]->GetTheta();
                        double phi = clusterCalifa[p]->GetPhi();

                        // cout << "Entered califa loop" << endl;

                        if ( E > 0 ){

                            h_CALIFAA->Fill(phi);

                            // cout << "Enterd energy conditional" << "\n";


                            if ( (phi < 1.5 && phi > -1.5) && califaL == false){

                                // cout << "Enterd phi left conditional" << "\n";

                                califaL = true;
                                thetaTL = theta;
                                thetasum_plane += thetaTL;
                                phiTL = phi;

                                if (califaR == true && califaL == true){
                                    if ( abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
                                        h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                        h_CALIFAPHI->Fill(phiTL,phiTR);

                                        h_CALIFATHETASUM->Fill(thetasum_plane);

                                    }
                                    califaL = false;
                                    califaR = false;
                                    thetaTL = 0;
                                    thetaTR = 0;
                                }
                                continue;
                            }

                            if ( (phi>1.5 || phi<-1.5) && califaR == false) { 
                                
                                // cout << "Enterd phi right conditional" << "\n";
                                
                                califaR = true;
                                thetaTR = theta;
                                thetasum_plane += thetaTR;
                                phiTR = phi;

                                if (califaR == true && califaL == true){
                                    if ( abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
                                        h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                        h_CALIFAPHI->Fill(phiTL,phiTR);

                                        h_CALIFATHETASUM->Fill(thetasum_plane);

                                    }
                                    califaL = false;
                                    califaR = false;
                                    thetaTL = 0;
                                    thetaTR = 0;
                                }
                                continue;
                            }

                            if ((califaR == true && califaL == false) || (califaR == false && califaL == true)){
                                // cout << "Cleaned" << "\n";
                                califaL = false;
                                califaR = false;
                            }

                        }
                    }
                }
            
	}



	

TCanvas* c_CALIFAE = new TCanvas();
c_CALIFAE->cd(1);
h_CALIFAE->Draw("hist");

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

}
