/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname[200];
char histname[200];
int tpat;

// Before GLAD
TH1D* h_TTT10_X;
TH1D* h_MWPC_X;
TH2D* h_TTT10_MWPC_X;

// After GLAD
TH1D* h_FI32_X;
TH1D* h_FI33_X;
TH2D* h_FI32_FI33_X;
TH1D* h_TOFD_X;

// Correlation
TH2D* h_TTT10_MWPC_X_corr;
TH2D* h_FI32_FI33_X_corr;

// Correlated area
TH2D* h_TTT10_MWPC_X_area;
TH2D* h_FI32_FI33_X_area;

// vector<double> Eloss1;
// vector<double> Eloss2;
// vector<double> CalifaR;
// vector<double> CalifaL;

void analysis_detectors_GLAD_xpos_corr(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Set histograms canvas and axis titles here  

    h_TTT10_X = new TH1D("TTT10 x position","TTT10 x position",32,-50,50);
    h_TTT10_X->GetXaxis()->SetTitle("TTT10 x position");

    h_MWPC_X = new TH1D("MWPC 1 x position","MWPC 1 x position",200,-150,0);
    h_MWPC_X->GetXaxis()->SetTitle("MWPC 1 x position");

    h_TTT10_MWPC_X = new TH2D("X position correlations MWPC-TTT10","X position correlations MWPC-TTT10",200,-150,0,32,-50,50);
    h_TTT10_MWPC_X->GetXaxis()->SetTitle("MWPC X-position");
    h_TTT10_MWPC_X->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI32_X = new TH1D("FI32 position","FI32 x position",200,-30,30);
    h_FI32_X->GetXaxis()->SetTitle("FI32 x position");

    h_FI33_X = new TH1D("FI33 position","FI33 x position",200,-30,30);
    h_FI33_X->GetXaxis()->SetTitle("FI33 x position");

    h_FI32_FI33_X = new TH2D("X position correlations FI32-FI33","X position correlations FI32-FI33",200,-30,30,200,-30,30);
    h_FI32_FI33_X->GetXaxis()->SetTitle("FI32 X-position");
    h_FI32_FI33_X->GetYaxis()->SetTitle("FI33 X-position");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_MWPC_X_corr = new TH2D("X position correlations MWPC-TTT10 correlation","X position correlations MWPC-TTT10 correlation",200,-150,0,32,-50,50);
    h_TTT10_MWPC_X_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_TTT10_MWPC_X_corr->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI32_FI33_X_corr = new TH2D("X position correlations FI32-FI33 correlation","X position correlations FI32-FI33 correlation",200,-30,30,200,-30,30);
    h_FI32_FI33_X_corr->GetXaxis()->SetTitle("FI32 X-position");
    h_FI32_FI33_X_corr->GetYaxis()->SetTitle("FI33 X-position");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_MWPC_X_area = new TH2D("X position correlations MWPC-TTT10 area","X position correlations MWPC-TTT10 area",200,-150,0,32,-50,50);
    h_TTT10_MWPC_X_area->GetXaxis()->SetTitle("MWPC X-position");
    h_TTT10_MWPC_X_area->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI32_FI33_X_area = new TH2D("X position correlations FI32-FI33 area","X position correlations FI32-FI33 area",200,-30,30,200,-30,30);
    h_FI32_FI33_X_area->GetXaxis()->SetTitle("FI32 X-position");
    h_FI32_FI33_X_area->GetYaxis()->SetTitle("FI33 X-position");

    // h_CALIFAE = new TH1D("Energy","Energy",1000,0,900000);
    // h_CALIFAE->GetXaxis()->SetTitle("Energy");

    // h_CALIFAA = new TH1D("Azimutal angle","Azimutal angle",100,-4,4);
    // h_CALIFAA->GetXaxis()->SetTitle("Polar angle");

    // h_CALIFATHETA = new TH2D("Angular correlation","Angular correlation",100,0,2,100,0,2);
    // h_CALIFATHETA->GetXaxis()->SetTitle("Angle 1");
    // h_CALIFATHETA->GetYaxis()->SetTitle("Angle 2");

    // h_CALIFAPHI = new TH2D("Azimutal angular correlation","Azimutal angular correlation",100,-3.2,3.2,100,-3.2,3.2);
    // h_CALIFAPHI->GetXaxis()->SetTitle("Angle 1");
    // h_CALIFAPHI->GetYaxis()->SetTitle("Angle 2");

    // h_CALIFAECORR = new TH2D("Energy correlations","Energy correlations",1000,0,900000,1000,0,900000);
    // h_CALIFAECORR->GetXaxis()->SetTitle("Energy 1");
    // h_CALIFAECORR->GetYaxis()->SetTitle("Energy 2");

    // h_CALIFATHETASUM = new TH1D("Opening angle","Opening angle",100,0,2);
    // h_CALIFATHETASUM->GetXaxis()->SetTitle("Opening angle");

    // h_TOFDX = new TH1D("TOFD x position","TOFD x position",200,-80,80);
    // h_TOFDX->GetXaxis()->SetTitle("TOFD x position");

    // h_FIBER30 = new TH1D("Fibers x position","Fibers x position",512,-25,25);
    // h_FIBER30->GetXaxis()->SetTitle("Fibers x position");

    // h_MWPC_FIBERS_X = new TH2D("X position correlations MWPC-FIBERS","X position correlations MWPC-FIBERS",200,-200,100,512,-25,25);
    // h_MWPC_FIBERS_X->GetXaxis()->SetTitle("MWPC X-position");
    // h_MWPC_FIBERS_X->GetYaxis()->SetTitle("FIBERS X-position");

    // h_MWPC_TOFD_X = new TH2D("X position correlations MWPC-TOFD","X position correlations MWPC-TOFD",200,-200,100,200,-100,100);
    // h_MWPC_TOFD_X->GetXaxis()->SetTitle("MWPC X-position");
    // h_MWPC_TOFD_X->GetYaxis()->SetTitle("TOFD X-position");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt");  
    //sprintf(fname, "./main0252_0001_unpacked_MDF.root");
    //sprintf(fname, "./main0321_0001_unpacked_MDF.root"); 
    sprintf(fname, "./main0327_0001_unpacked_MDF.root");

    chain->Reset();
    chain->Add(fname); 
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 

	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("EventHeader.");
    branchData->SetAddress(&header);    


    TClonesArray* hitTTT10 = new TClonesArray("R3BTttxHitData",5);
	R3BTttxHitData** hitTtt10;
	TBranch *branchHitTTT10 = chain->GetBranch("tttxHitData");
	branchHitTTT10->SetAddress(&hitTTT10);
	Long64_t Ttt10HitPerEvent = 0;  	


    TClonesArray* hitMWPC = new TClonesArray("R3BMwpcHitData",5);
	R3BMwpcHitData** hitMwpc;
	TBranch *branchHitMWPC = chain->GetBranch("Mwpc1HitData");
	branchHitMWPC->SetAddress(&hitMWPC);
	Long64_t MwpcHitPerEvent = 0;  	


    TClonesArray* hitFI32 = new TClonesArray("R3BFiberMAPMTHitData",5);
	R3BFiberMAPMTHitData** hitFi32;
	TBranch *branchHitFI32 = chain->GetBranch("Fi32Hit");
	branchHitFI32->SetAddress(&hitFI32);
	Long64_t Fi32HitPerEvent = 0;  	


    TClonesArray* hitFI33 = new TClonesArray("R3BFiberMAPMTHitData",5);
	R3BFiberMAPMTHitData** hitFi33;
	TBranch *branchHitFI33 = chain->GetBranch("Fi33Hit");
	branchHitFI33->SetAddress(&hitFI33);
	Long64_t Fi33HitPerEvent = 0;  	


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

    /////////////////////////////////////////////////////////////////////////////////////////
	
	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in " << fname << " : " << nevent << endl;

	for(Long64_t e=0; e<nevent; e++)
	{   
        hitTTT10->Clear();
        hitMWPC->Clear();

        hitFI32->Clear();
        hitFI33->Clear();
	    hitTOFD->Clear();

        hitCALIFA->Clear();
	    chain->GetEvent(e);
	    
        // Entries of LOS and TOFD per event in the chain of trees
        tpat             = header->GetTpat();
        Ttt10HitPerEvent = hitTTT10->GetEntriesFast();
        MwpcHitPerEvent  = hitMWPC->GetEntriesFast();

        Fi32HitPerEvent = hitFI32->GetEntriesFast();
        Fi33HitPerEvent = hitFI33->GetEntriesFast();
	    TofdHitPerEvent = hitTOFD->GetEntriesFast();

        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        MwpcHitPerEvent   = hitMWPC->GetEntriesFast();

        /////////////////////////////////////////////////////////////////////////////////////////

        if( Ttt10HitPerEvent > 0 ){

            hitTtt10 = new R3BTttxHitData*[Ttt10HitPerEvent];
            
            for(int p = 0; p < Ttt10HitPerEvent; p++){

                hitTtt10[p] =(R3BTttxHitData*)hitTTT10->At(p);
                double TTT10_X = hitTtt10[p]->GetX();
                h_TTT10_X->Fill(TTT10_X);
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////

        double MWPC_X;

        if( MwpcHitPerEvent > 0 ){

            hitMwpc = new R3BMwpcHitData*[MwpcHitPerEvent];
            
            for(int p = 0; p < MwpcHitPerEvent; p++){

                hitMwpc[p] =(R3BMwpcHitData*)hitMWPC->At(p);
                MWPC_X = hitMwpc[p]->GetX();
                h_MWPC_X->Fill(MWPC_X);

                /////////////////////////////////////////////////////////////////////////////////////////
                
                hitTtt10 = new R3BTttxHitData*[Ttt10HitPerEvent];
                
                for(int p = 0; p < Ttt10HitPerEvent; p++){

                    hitTtt10[p] =(R3BTttxHitData*)hitTTT10->At(p);
                    double TTT10_X = hitTtt10[p]->GetX();
                    h_TTT10_MWPC_X->Fill(MWPC_X, TTT10_X);
                    h_TTT10_MWPC_X_corr->Fill(MWPC_X, TTT10_X);

                    if ( TTT10_X < 0.6*MWPC_X + 45 && TTT10_X > 0.6*MWPC_X + 39 ){
                        h_TTT10_MWPC_X_area->Fill(MWPC_X, TTT10_X);
                    }
                } 
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////

        if( Fi32HitPerEvent > 0 ){

            hitFi32 = new R3BFiberMAPMTHitData*[Fi32HitPerEvent];
            
            for(int p = 0; p < Fi32HitPerEvent; p++){

                hitFi32[p] =(R3BFiberMAPMTHitData*)hitFI32->At(p);
                double hitFI32_X = hitFi32[p]->GetX();
                h_FI32_X->Fill(hitFI32_X);
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////

        double hitFI33_X;

        if( Fi33HitPerEvent > 0 ){

            hitFi33 = new R3BFiberMAPMTHitData*[Fi33HitPerEvent];
            
            for(int p = 0; p < Fi33HitPerEvent; p++){

                hitFi33[p] =(R3BFiberMAPMTHitData*)hitFI33->At(p);
                hitFI33_X = hitFi33[p]->GetX();
                h_FI33_X->Fill(hitFI33_X);

                /////////////////////////////////////////////////////////////////////////////////////////

                hitFi32 = new R3BFiberMAPMTHitData*[Fi32HitPerEvent];
                
                for(int p = 0; p < Fi32HitPerEvent; p++){

                    hitFi32[p] =(R3BFiberMAPMTHitData*)hitFI32->At(p);
                    double hitFI32_X = hitFi32[p]->GetX();
                    h_FI32_FI33_X->Fill(hitFI32_X,hitFI33_X);
                    h_FI32_FI33_X_corr->Fill(hitFI32_X,hitFI33_X);

                    if ( hitFI33_X < 1.08*hitFI32_X -20.5  && hitFI33_X > 1.08*hitFI32_X -26.5  ){
                        h_TTT10_MWPC_X_area->Fill(hitFI32_X, hitFI33_X);
                    }
                }
            }
        }

    //     if ((tpat & 8) == 8){
            
    //         /////////////////////////////////////////////////////////////////////////////////////////

    //         if( TofdHitPerEvent > 0 ){

    //             hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
                
    //             for(int p = 0; p < TofdHitPerEvent; p++){

    //                 hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
    //                 double TOFD_X = hitTofd[p]->GetX();
    //                 h_TOFDX->Fill(TOFD_X);
    //             }
    //         }

    //         /////////////////////////////////////////////////////////////////////////////////////////

    //         bool califaL=false;
    //         bool califaR=false;
    //         double thetaTL, thetaTR, phiTR, phiTL, EL, ER;
    //         double phi1 = 0., phi2 =0., theta1=0., theta2=0.;

    //         if( CalifaHitPerEvent == 2 ){

    //             clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];
    //             hitMwpc       = new R3BMwpcHitData*[MwpcHitPerEvent];
    //             hitFi30       = new R3BFiberMAPMTHitData*[Fi30HitPerEvent];
    //             hitTofd       = new R3BTofdHitData*[TofdHitPerEvent];

    //             double thetasum_plane = 0;

    //             // ######################### CALIFA #######################################

    //             for(int p = 0; p < CalifaHitPerEvent; p++){

    //                 clusterCalifa[p] =(R3BCalifaClusterData*)hitCALIFA->At(p);

    //                 double E = clusterCalifa[p]->GetEnergy();
    //                 double theta = clusterCalifa[p]->GetTheta();
    //                 double phi = clusterCalifa[p]->GetPhi();

    //                 if(p==0){
    //                     phi1 = phi;
    //                     theta1 = theta; 
    //                 }
    //                 if(p==1){
    //                     phi2 = phi;
    //                     theta2 = theta;
    //                 }

    //                 if ( E > 20000 ){
                        
    //                     h_CALIFAA->Fill(phi);

    //                     if ( (phi < 1.5 && phi > -1.5) && califaL == false){

    //                         califaL = true;
    //                         thetaTL = theta;
    //                         thetasum_plane += thetaTL;
    //                         phiTL = phi;
    //                         EL = E;

    //                         if (califaR == true && califaL == true){
    //                             if (  abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
    //                                 //h_CALIFATHETA->Fill(thetaTL,thetaTR);
    //                                 //h_CALIFAPHI->Fill(phiTL,phiTR);
    //                                 h_CALIFATHETASUM->Fill(thetasum_plane);
    //                                 h_CALIFAECORR->Fill(EL,ER);
    //                                 h_CALIFAE->Fill(EL+ER);
    //                             }
    //                             califaL = false;
    //                             califaR = false;
    //                             thetaTL = 0;
    //                             thetaTR = 0;
    //                         }
    //                         continue;
    //                     }

    //                     if ( (phi>1.5 || phi<-1.5) && califaR == false) { 
                                                            
    //                         califaR = true;
    //                         thetaTR = theta;
    //                         thetasum_plane += thetaTR;
    //                         phiTR = phi;
    //                         ER = E;

    //                         if (califaR == true && califaL == true){
    //                             if (  abs(phiTL) + abs(phiTR) - 3.1415926535 < 0.2) {
    //                                 //h_CALIFATHETA->Fill(thetaTL,thetaTR);
    //                                 //h_CALIFAPHI->Fill(phiTL,phiTR);
    //                                 h_CALIFATHETASUM->Fill(thetasum_plane);
    //                                 h_CALIFAECORR->Fill(EL,ER);
    //                                 h_CALIFAE->Fill(EL+ER);
    //                             }
    //                             califaL = false;
    //                             califaR = false;
    //                             thetaTL = 0;
    //                             thetaTR = 0;
    //                         }
    //                         continue;
    //                     }

    //                     if ((califaR == true && califaL == false) || (califaR == false && califaL == true)){
    //                         califaL = false;
    //                         califaR = false;
    //                     }

    //                 }
    //             }

    //             if (abs(abs(phi1)+abs(phi2))- 3.1415926535 < 0.2){
    //                 h_CALIFAPHI->Fill(phi1,phi2);
    //                 h_CALIFATHETA->Fill(theta1,theta2);
    //             }

    //             // ######################### CALIFA #######################################

    //             // ############################# MWPC #####################################
                
    //             for(int p = 0; p < MwpcHitPerEvent; p++){

    //                 hitMwpc[p] =(R3BMwpcHitData*)hitMWPC->At(p);
    //                 double x = hitMwpc[p]->GetX();
    //                 h_MWPC->Fill(x);

    //                 for(int q = 0; q < Fi30HitPerEvent; q++){
    //                     hitFi30[q] =(R3BFiberMAPMTHitData*)hitFI30->At(q);
    //                     double xfib = hitFi30[q]->GetX();
    //                     h_MWPC_FIBERS_X->Fill(x, xfib); 
    //                 }
    //             }

    //             // ############################# MWPC #####################################

    //             // ############################# TOFD #####################################

    //             for(int p = 0; p < MwpcHitPerEvent; p++){

    //                 hitMwpc[p] =(R3BMwpcHitData*)hitMWPC->At(p);
    //                 double x = hitMwpc[p]->GetX();
    //                 h_MWPC->Fill(x);

    //                 for(int q = 0; q < TofdHitPerEvent; q++){

    //                     hitTofd[q] =(R3BTofdHitData*)hitTOFD->At(q);
    //                     double TOFD_X = hitTofd[p]->GetX();
    //                     h_MWPC_TOFD_X->Fill(x, TOFD_X); 
    //                 }
    //             }

    //             // ############################# TOFD #####################################
                
    //             // ############################# FIBERS ###################################
                
    //             for(int p = 0; p < Fi30HitPerEvent; p++){
    //                 hitFi30[p] =(R3BFiberMAPMTHitData*)hitFI30->At(p);
    //                 double x = hitFi30[p]->GetX();
    //                 h_FIBER30->Fill(x);
    //             }
                
    //             // ############################# FIBERS ###################################
    //         }  
    //     }  
	}

///////////////////////////////////////////////////////////////////

TCanvas* c_TTT10_X  = new TCanvas();
c_TTT10_X ->cd(1);
h_TTT10_X ->Draw("hist");

TCanvas* c_MWPC_X  = new TCanvas();
c_MWPC_X ->cd(1);
h_MWPC_X ->Draw("hist");

TCanvas* c_TTT10_MWPC_X = new TCanvas();
c_TTT10_MWPC_X->cd(1);
h_TTT10_MWPC_X->Draw("colz");

TCanvas* c_FI32_X  = new TCanvas();
c_FI32_X->cd(1);
h_FI32_X->Draw("hist");

TCanvas* c_FI33_X  = new TCanvas();
c_FI33_X->cd(1);
h_FI33_X->Draw("hist");

TCanvas* c_FI32_FI33_X = new TCanvas();
c_FI32_FI33_X->cd(1);
h_FI32_FI33_X->Draw("colz");

///////////////////////////////////////////////////////////////////

double m1 = 0.6;
double b1 = 42;
double x1_1 = h_TTT10_MWPC_X_corr->GetXaxis()->GetXmin();
double x2_1 = h_TTT10_MWPC_X_corr->GetXaxis()->GetXmax();
double y1_1 = m1*x1_1 + b1;
double y2_1 = m1*x2_1 + b1;

TCanvas *c1 = new TCanvas();
h_TTT10_MWPC_X_corr->Draw("colz");

TLine *line1 = new TLine(x1_1, y1_1, x2_1, y2_1);
line1->SetLineColor(kRed);
line1->SetLineWidth(2);
line1->Draw("same");

TF1 *f1 = new TF1("line1", "[0]*x + [1]", x1_1, x2_1);
f1->SetParameters(m1, b1);
f1->SetLineColor(kRed);
f1->SetLineWidth(1);

h_TTT10_MWPC_X_corr->Draw("colz");
f1->Draw("same");

///////////////////////////////////////////////////////////////////

double m2 = 1.08;
double b2 = -23.5;
double x1_2 = h_FI32_FI33_X->GetXaxis()->GetXmin();
double x2_2 = h_FI32_FI33_X->GetXaxis()->GetXmax();
double y1_2 = m2*x1_2 + b2;
double y2_2 = m2*x2_2 + b2;

TCanvas *c2 = new TCanvas();
h_FI32_FI33_X_corr->Draw("colz");

TLine *line2 = new TLine(x1_2, y1_2, x2_2, y2_2);
line2->SetLineColor(kRed);
line2->SetLineWidth(2);
line2->Draw("same");

TF1 *f2 = new TF1("line2", "[0]*x + [1]", x1_2, x2_2);
f2->SetParameters(m2, b2);
f2->SetLineColor(kRed);
f2->SetLineWidth(1);

h_FI32_FI33_X_corr->Draw("colz");
f2->Draw("same");

///////////////////////////////////////////////////////////////////

TCanvas* c_TTT10_MWPC_X_area = new TCanvas();
c_TTT10_MWPC_X_area->cd(1);
h_TTT10_MWPC_X_area->Draw("colz");

TCanvas* c_FI32_FI33_X_area = new TCanvas();
c_FI32_FI33_X_area->cd(1);
h_FI32_FI33_X_area->Draw("colz");

}