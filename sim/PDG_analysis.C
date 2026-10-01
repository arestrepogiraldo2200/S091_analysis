
char fname[200];

TH1D* h_CALIFAE;
TH1D* h_CALIFATHETASUM;

TH2D* h_CALIFATHETA;
TH2D* h_CALIFAPHI;
TH2D* h_CALIFAECORR;


void PDG_analysis(){

    h_CALIFAE = new TH1D("Energy","Energy",1000,0,900);
    h_CALIFAE->GetXaxis()->SetTitle("Energy");

    h_CALIFATHETASUM = new TH1D("Opening angle","Opening angle",100,0,2);
    h_CALIFATHETASUM->GetXaxis()->SetTitle("Opening angle");

    h_CALIFATHETA = new TH2D("Polar angular correlation","Polar angular correlation",60,0,2,60,0,2);
    h_CALIFATHETA->GetXaxis()->SetTitle("Angle 1");
    h_CALIFATHETA->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAPHI = new TH2D("Azimutal angular correlation","Azimutal angular correlation",60,-3.2,3.2,60,-3.2,3.2);
    h_CALIFAPHI->GetXaxis()->SetTitle("Angle 1");
    h_CALIFAPHI->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAECORR = new TH2D("Energy correlations","Energy correlations",1000,0,900,1000,0,900);
    h_CALIFAECORR->GetXaxis()->SetTitle("Energy 1");
    h_CALIFAECORR->GetYaxis()->SetTitle("Energy 2");

    /////////////////////////////////////////////////////////////////////////////////////////

    TChain* chain = new TChain("evt");  
    sprintf(fname, "./sim_16C.root");

    chain->Reset();
    chain->Add(fname); 

    /////////////////////////////////////////////////////////////////////////////////////////

	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("MCEventHeader.");
    branchData->SetAddress(&header); 

    TClonesArray* mcTrack = new TClonesArray("R3BMCTrack", 5);
	R3BMCTrack** track;
	TBranch *branchMCTrack = chain->GetBranch("MCTrack");
	branchMCTrack->SetAddress(&mcTrack);
	Long64_t mcTrackPerEvent = 0;           

    TClonesArray*  hitCALIFA = new TClonesArray("R3BCalifaClusterData",5);
	R3BCalifaClusterData** clusterCalifa;
	TBranch *branchHitCALIFA = chain->GetBranch("CalifaClusterData");
	branchHitCALIFA->SetAddress(&hitCALIFA);
	Long64_t CalifaHitPerEvent = 0;  


    Long64_t nevent = chain->GetEntries();
	cout << "Number of events in " << fname << " : " << nevent << endl;

    int twotracks_event = 0;

    /////////////////////////////////////////////////////////////////////////////////////////

	for(Long64_t e=0; e<nevent; e++){

        mcTrack->Clear();
        chain->GetEvent(e);

	    mcTrackPerEvent = mcTrack->GetEntriesFast();                       // Number of tracks in an event
        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();

        track = new R3BMCTrack*[mcTrackPerEvent];
        clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];

        bool protons_bool = true;                                          // Variable to guarantee the tracks are only protons

        if (mcTrackPerEvent == 2){

            twotracks_event++;

            ///////////////////////////////////////////////////////////////////////////
            for(int p = 0; p < mcTrackPerEvent; p++){

                track[p] =(R3BMCTrack*)mcTrack->At(p);
                int pdgcode = track[p]->GetPdgCode();

                if (pdgcode != 2212){
                    protons_bool = false;
                    break;
                }
            }

            ///////////////////////////////////////////////////////////////////////////
            if (protons_bool == true && CalifaHitPerEvent == 2){                         // Number of hits in califa

                double thetasum = 0;
                double theta1, theta2, phi1, phi2, E1, E2;

                for(int p = 0; p < CalifaHitPerEvent; p++){

                    clusterCalifa[p] = (R3BCalifaClusterData*)hitCALIFA->At(p);
                    double E     = clusterCalifa[p]->GetEnergy();
                    double theta = clusterCalifa[p]->GetTheta();
                    double phi   = clusterCalifa[p]->GetPhi();

                    h_CALIFAE->Fill(E);
                    thetasum = thetasum + theta;  
                    
                    if (p == 0){
                        phi1 = phi;
                        theta1 = theta;
                        E1 = E;
                    }

                    if (p == 1){
                        phi2 = phi;
                        theta2 = theta; 
                        E2 = E;   
                    }          
                }

                h_CALIFATHETASUM->Fill(thetasum);
                h_CALIFATHETA->Fill(theta1, theta2);
                h_CALIFAPHI->Fill(phi1, phi2);
                h_CALIFAECORR->Fill(E1, E2);
            }
        }
    }

    cout << "Total events with two tracks: " << twotracks_event << endl;

    TCanvas* c_CALIFAE = new TCanvas();
    c_CALIFAE->cd(1);
    h_CALIFAE->Draw("hist");

    TCanvas* c_CALIFATHETASUM = new TCanvas();
    c_CALIFATHETASUM->cd(1);
    h_CALIFATHETASUM->Draw("hist");

    TCanvas* c_CALIFATHETA = new TCanvas();
    c_CALIFATHETA->cd(1);
    h_CALIFATHETA->Draw("colz");

    TCanvas* c_CALIFAPHI = new TCanvas();
    c_CALIFAPHI->cd(1);
    h_CALIFAPHI->Draw("colz");

    TCanvas* c_CALIFAECORR = new TCanvas();
    c_CALIFAECORR->cd(1);
    h_CALIFAECORR->Draw("colz");
    
}