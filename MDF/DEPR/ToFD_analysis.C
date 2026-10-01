/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname1[200];
char histname[200];
int tpat;

// Before GLAD
TH1D* h_TTT10_X;
TH1D* h_MWPC_X;
TH2D* h_TTT10_MWPC_X;

// After GLAD
TH1D* h_FI31_X;
TH1D* h_FI32_X;
TH1D* h_FI33_X;
TH1D* h_FICOMB_X;
TH1D* h_TOFD_X;

// Hits per event histograms
TH1D* h_hits_TTT10;
TH1D* h_hits_MWPC;
TH1D* h_hits_FI31;
TH1D* h_hits_FI32;
TH1D* h_hits_FI33;
TH1D* h_hits_TOFD;
TH1D* h_hits_LOS;

TH1D* h_TOFD_1;
TH1D* h_TOFD_2;
TH1D* h_TOFD_3;
TH1D* h_TOFD_4;


// Correlation
TH2D* h_TTT10_MWPC_X_corr;
TH2D* h_FI33_FI32_X;
TH2D* h_FI31_FI32_X;
TH2D* h_FI33_FI31_combined_FI32_X;
TH2D* h_Q1_Q2;

// Correlated area
TH2D* h_TTT10_MWPC_X_area;
TH2D* h_FIB31_FIB33_combined_FIB32_area;

// Correlations between detectors before and after
TH2D* h_TTT10_FI32_X_corr;
TH2D* h_TTT10_FI33_X_corr;
TH2D* h_MWPC_FI32_X_corr;
TH2D* h_MWPC_FI33_X_corr;
TH2D* h_MWPC_FI31_FI33_X_single_corr;
TH2D* h_TTT10_FI32_X_single_corr;
TH2D* h_TTT10_FI33_X_single_corr;
TH2D* h_MWPC_FI32_X_single_corr;
TH2D* h_MWPC_FI33_X_single_corr;

vector<R3BTttxHitData*> TTT10_corr_hit;
vector<R3BMwpcHitData*>  MWPC_corr_hit;

vector<R3BTttxHitData*> TTT10_corr_single_hit;
vector<R3BFiberMAPMTHitData> FI32_corr_hit;
vector<R3BFiberMAPMTHitData> FI33_corr_hit;
vector<R3BFiberMAPMTHitData> FI31_corr_hit;

// Calculate the area of the triangle between hits in the fiber detectors
double TriangleArea(double x1, double y1, double x2, double y2,  double x3, double y3){
    return 0.5*std::abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)));
}

void ToFD_analysis(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Set histograms canvas and axis titles here

    int fiber_bins = 150; 
    int mwpc1_bins = 200; 

    h_TTT10_X = new TH1D("TTT10 x position","TTT10 x position",32,-50,50);
    h_TTT10_X->GetXaxis()->SetTitle("TTT10 x position");

    h_MWPC_X = new TH1D("MWPC1 x position","MWPC1 x position",mwpc1_bins,-150,0);
    h_MWPC_X->GetXaxis()->SetTitle("MWPC1 x position");

    h_TTT10_MWPC_X = new TH2D("X position correlations MWPC-TTT10","X position correlations MWPC-TTT10",mwpc1_bins,-150,0,32,-50,50);
    h_TTT10_MWPC_X->GetXaxis()->SetTitle("MWPC X-position");
    h_TTT10_MWPC_X->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI31_X = new TH1D("FI31 x position","FI31 x position",fiber_bins,-26.4,26.4);
    h_FI31_X->GetXaxis()->SetTitle("FI31 x position [cm]");

    h_FI32_X = new TH1D("FI32 x position","FI32 x position",fiber_bins,-26.4,26.4);
    h_FI32_X->GetXaxis()->SetTitle("FI32 x position [cm]");

    h_FI33_X = new TH1D("FI33 x position","FI33 x position",fiber_bins,-26.4,26.4);
    h_FI33_X->GetXaxis()->SetTitle("FI33 x position [cm]");



    h_hits_FI31 = new TH1D("Hits per event 31", "Hits per event 31",20, 0, 20);
    h_hits_FI31->GetXaxis()->SetTitle("Hits per event 31");

    h_hits_FI32 = new TH1D("Hits per event 32", "Hits per event 32",20, 0, 20);
    h_hits_FI32->GetXaxis()->SetTitle("Hits per event 32");

    h_hits_FI33 = new TH1D("Hits per event 33", "Hits per event 33",20, 0, 20);
    h_hits_FI33->GetXaxis()->SetTitle("Hits per event 33");

    h_hits_TTT10 = new TH1D("Hits per event TTT10", "Hits per event TTT10",20, 0, 20);
    h_hits_TTT10->GetXaxis()->SetTitle("Hits per event TTT10");

    h_hits_MWPC = new TH1D("Hits per event MWPC", "Hits per event MWPC",20, 0, 20);
    h_hits_MWPC->GetXaxis()->SetTitle("Hits per event MWPC");

    h_hits_TOFD = new TH1D("Hits per event ToFD", "Hits per event ToFD",20, 0, 20);
    h_hits_TOFD->GetXaxis()->SetTitle("Hits per event ToFD");

    h_hits_LOS = new TH1D("Hits per event LOS", "Hits per event LOS",20, 0, 20);
    h_hits_LOS->GetXaxis()->SetTitle("Hits per event LOS");




    h_FICOMB_X = new TH1D("FI31 & FI33 x position","FI31 & FI33 x position", 2*fiber_bins,-52.8, 52.8);
    h_FICOMB_X->GetXaxis()->SetTitle("FI31 & FI33 x position [cm]");

    h_TOFD_X = new TH1D("TOFD x position","TOFD x position",100,-60,60);
    h_TOFD_X->GetXaxis()->SetTitle("TOFD x position [cm]");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_MWPC_X_corr = new TH2D("X position MWPC-TTT10 correlation","X position MWPC-TTT10 correlation",mwpc1_bins,-150,0,32,-50,50);
    h_TTT10_MWPC_X_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_TTT10_MWPC_X_corr->GetYaxis()->SetTitle("TTT10 X-position");

    h_TTT10_MWPC_X_area = new TH2D("X position correlations MWPC-TTT10 area","X position correlations MWPC-TTT10 area",mwpc1_bins,-150,0,32,-50,50);
    h_TTT10_MWPC_X_area->GetXaxis()->SetTitle("MWPC X-position" );
    h_TTT10_MWPC_X_area->GetYaxis()->SetTitle("TTT10 X-position");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_FI32_X_corr = new TH2D("X position TTT10-FIB32 correlation","X position TTT10-FIB32 correlation",32,-50,50,fiber_bins,-26.4,26.4);
    h_TTT10_FI32_X_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI32_X_corr->GetYaxis()->SetTitle("FIB32 X-position [cm]");

    h_TTT10_FI33_X_corr = new TH2D("X position TTT10-FIB33 correlation","X position TTT10-FIB33 correlation",32,-50,50,fiber_bins,-26.4,26.4);
    h_TTT10_FI33_X_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI33_X_corr->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    h_MWPC_FI32_X_corr = new TH2D("X position MWPC-FIB32 correlation","X position MWPC-FIB32 correlation",mwpc1_bins,-150,0,fiber_bins,-26.4,26.4);
    h_MWPC_FI32_X_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_MWPC_FI32_X_corr->GetYaxis()->SetTitle("FIB32 X-position [cm]");

    h_MWPC_FI33_X_corr = new TH2D("X position MWPC-FIB33 correlation","X position MWPC-FIB33 correlation",mwpc1_bins,-150,0,fiber_bins,-26.4,26.4);
    h_MWPC_FI33_X_corr->GetXaxis()->SetTitle("MWPC X-position" );
    h_MWPC_FI33_X_corr->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_FI32_X_single_corr = new TH2D("X position TTT10-FIB32 correlation single hit","X position TTT10-FIB32 correlation single hit",32,-50,50, fiber_bins,-26.4,26.4);
    h_TTT10_FI32_X_single_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI32_X_single_corr->GetYaxis()->SetTitle("FIB32 X-position [cm]");

    h_TTT10_FI33_X_single_corr = new TH2D("X position TTT10-FIB33 correlation single hit","X position TTT10-FIB33 correlation single hit",32,-50,50, fiber_bins,-26.4,26.4);
    h_TTT10_FI33_X_single_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI33_X_single_corr->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    h_MWPC_FI32_X_single_corr = new TH2D("X position MWPC-FIB32 correlation single hit","X position MWPC-FIB32 correlation single hit", mwpc1_bins,-150,0,fiber_bins,-26.4,26.4);
    h_MWPC_FI32_X_single_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_MWPC_FI32_X_single_corr->GetYaxis()->SetTitle("FIB32 X-position [cm]");

    h_MWPC_FI33_X_single_corr = new TH2D("X position MWPC-FIB33 correlation single hit","X position MWPC-FIB33 correlation single hit", mwpc1_bins,-150,0,fiber_bins,-26.4,26.4);
    h_MWPC_FI33_X_single_corr->GetXaxis()->SetTitle("MWPC X-position" );
    h_MWPC_FI33_X_single_corr->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    h_FI33_FI31_combined_FI32_X = new TH2D("X position FIB31&FIB33-FIB32 correlation","X position FIB31&FIB33-FIB32 correlation", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FI33_FI31_combined_FI32_X->GetXaxis()->SetTitle("FIB32 X-position [cm]");
    h_FI33_FI31_combined_FI32_X->GetYaxis()->SetTitle("FIB31&FIB33 X-position [cm]");

    h_FIB31_FIB33_combined_FIB32_area = new TH2D("X position FIB31&FIB33-FIB32 combined correlation area","X position FIB31&FIB33-FIB32 combined correlation area", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FIB31_FIB33_combined_FIB32_area->GetXaxis()->SetTitle("FIB32 X-position [cm]" );
    h_FIB31_FIB33_combined_FIB32_area->GetYaxis()->SetTitle("FIB31&FIB31 X-position [cm]");

    h_FI33_FI32_X = new TH2D("X position FIB33-FIB32 correlation area","X position FIB33-FIB32 correlation area", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FI33_FI32_X->GetXaxis()->SetTitle("FIB32 X-position [cm]" );
    h_FI33_FI32_X->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    h_FI31_FI32_X = new TH2D("X position FIB31-FIB32 correlation area","X position FIB31-FIB32 correlation area", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FI31_FI32_X->GetXaxis()->SetTitle("FIB32 X-position [cm]" );
    h_FI31_FI32_X->GetYaxis()->SetTitle("FIB31 X-position [cm]");

    h_MWPC_FI31_FI33_X_single_corr = new TH2D("X position MWPC and FIB31-FIB33 correlation single hit","X position MWPC and FIB31-FIB33 correlation single hit", mwpc1_bins, -150, 0, 300, -52.8, 52.8);
    h_MWPC_FI31_FI33_X_single_corr->GetXaxis()->SetTitle("MWPC X-position [mm]" );
    h_MWPC_FI31_FI33_X_single_corr->GetYaxis()->SetTitle("FIB31&FIB33 X-position [cm]");







    h_Q1_Q2 = new TH2D("Charge correlations","Charge correlations", 40, 0, 8, 40, 0, 8);
    h_Q1_Q2->GetXaxis()->SetTitle("Charge plane 1");
    h_Q1_Q2->GetYaxis()->SetTitle("Charge plane 2");

    h_TOFD_1 = new TH1D("Energy loss for plane 1","Energy loss for plane 1",40, 0, 8);
    h_TOFD_1->GetXaxis()->SetTitle("Eloss plane 1 [Q]");

    h_TOFD_2 = new TH1D("Energy loss for plane 2","Energy loss for plane 2",40, 0, 8);
    h_TOFD_2->GetXaxis()->SetTitle("Eloss plane 2 [Q]");

    h_TOFD_3 = new TH1D("Energy loss for plane 3","Energy loss for plane 3",40, 0, 8);
    h_TOFD_3->GetXaxis()->SetTitle("Eloss plane 3 [Q]");

    h_TOFD_4 = new TH1D("Energy loss for plane 4","Energy loss for plane 4",40, 0, 8);
    h_TOFD_4->GetXaxis()->SetTitle("Eloss plane 4 [Q]");






    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt"); 
    //sprintf(fname1,"./merged_MDF_runs.root"); 
    sprintf(fname1,"./merged_MDF_runs_v2.root");
    //sprintf(fname1,"./merged_partial.root"); 
    //sprintf(fname1,"./main0252_0001_unpacked_MDF.root");

    chain->Reset();
    chain->Add(fname1); 

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


    TClonesArray* hitFI31 = new TClonesArray("R3BFiberMAPMTHitData",5);
	R3BFiberMAPMTHitData** hitFi31;
	TBranch *branchHitFI31 = chain->GetBranch("Fi31Hit");
	branchHitFI31->SetAddress(&hitFI31);
	Long64_t Fi31HitPerEvent = 0;  	


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


    TClonesArray*  hitLOS = new TClonesArray("R3BLosHitData",5); 
	R3BLosHitData** hitLos;                                      
	TBranch *branchHitLOS = chain->GetBranch("LosHit");           
	branchHitLOS->SetAddress(&hitLOS);                            
	Long64_t LosHitPerEvent = 0;                                 

	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in the merged tree: " << nevent << endl;

    int good_events = 0;
    int collinear = 0, noncollinear = 0;
    int evts1 = 0;
    int evts2 = 0;
    int evts3 = 0;
    int evts4 = 0;
    int evts5 = 0;

	for(Long64_t e=0; e < nevent; e++)
    {

        hitTTT10->Clear();
        hitMWPC->Clear();

        hitFI31->Clear();
        hitFI32->Clear();
        hitFI33->Clear();
	    hitTOFD->Clear();

        hitCALIFA->Clear();
	    chain->GetEvent(e);

        // Cleaning the vectors 
        TTT10_corr_hit.clear();
        MWPC_corr_hit.clear();
        TTT10_corr_single_hit.clear();

        FI32_corr_hit.clear();
        FI33_corr_hit.clear();
        FI31_corr_hit.clear();

        // Entries of LOS and TOFD per event in the chain of trees
        tpat             = header->GetTpat();
        Ttt10HitPerEvent = hitTTT10->GetEntriesFast();
        MwpcHitPerEvent  = hitMWPC->GetEntriesFast();
	    LosHitPerEvent  = hitLOS->GetEntriesFast();

        Fi31HitPerEvent = hitFI31->GetEntriesFast();
        Fi32HitPerEvent = hitFI32->GetEntriesFast();
        Fi33HitPerEvent = hitFI33->GetEntriesFast();
	    TofdHitPerEvent = hitTOFD->GetEntriesFast();

        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        MwpcHitPerEvent   = hitMWPC->GetEntriesFast();

        h_hits_FI31->Fill(Fi31HitPerEvent);
        h_hits_FI32->Fill(Fi32HitPerEvent);
        h_hits_FI33->Fill(Fi33HitPerEvent);
        h_hits_TTT10->Fill(Ttt10HitPerEvent);
        h_hits_MWPC->Fill(MwpcHitPerEvent);
        h_hits_TOFD->Fill(TofdHitPerEvent);
        h_hits_LOS->Fill(LosHitPerEvent);

        /////////////////////////////////////////////////////////////////////////////////////////

        if ( TofdHitPerEvent == 4 ){

            if( TofdHitPerEvent > 0 ){                         // This condition is overrun by the previous conditional 

                hitTofd = new R3BTofdHitData*[TofdHitPerEvent];

                double Eloss1=0, Eloss2=0, Eloss3=0, Eloss4=0;
                
                for(int p = 0; p < TofdHitPerEvent; p++){

                    hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);

                    double Eloss = hitTofd[p]->GetEloss();
                    int PlaneID = hitTofd[p]->GetDetId();

                    if (PlaneID==1){
                        Eloss1 = Eloss;
                    }

                    if (PlaneID==2){
                        Eloss2 = Eloss;
                    }

                    if (PlaneID==3){
                        Eloss3 = Eloss;
                    }

                    if (PlaneID==4){
                        Eloss4 = Eloss;
                    }
                }

                if (Eloss1 != 0 && Eloss2 != 0){
                    h_Q1_Q2->Fill(Eloss1, Eloss2);
                }

                if (Eloss1 != 0 ){
                    h_TOFD_1->Fill(Eloss1);
                }
 
                if (Eloss2 != 0 ){
                    h_TOFD_2->Fill(Eloss2);
                }

                if (Eloss3 != 0 ){
                    h_TOFD_3->Fill(Eloss3);
                }

                if (Eloss4 != 0 ){
                    h_TOFD_4->Fill(Eloss4);
                }

            }
        }         
    }          

gSystem->mkdir("figures", true);

TFile* outFile = new TFile("figures/ToFD_analysis.root", "RECREATE");

TCanvas* c_Q1_Q2 = new TCanvas();
c_Q1_Q2->cd(1);
h_Q1_Q2->Draw("colz");
c_Q1_Q2->SaveAs("figures/Q1_Q2.png");
c_Q1_Q2->Write();

TCanvas* c_TOFD_1  = new TCanvas();
c_TOFD_1->cd(1);
h_TOFD_1->Draw("hist");
c_TOFD_1->SaveAs("figures/TOFD_1.png");
c_TOFD_1->Write();

TCanvas* c_TOFD_2  = new TCanvas();
c_TOFD_2->cd(1);
h_TOFD_2->Draw("hist");
c_TOFD_2->SaveAs("figures/TOFD_2.png");
c_TOFD_2->Write();

TCanvas* c_TOFD_3  = new TCanvas();
c_TOFD_3->cd(1);
h_TOFD_3->Draw("hist");
c_TOFD_3->SaveAs("figures/TOFD_3.png");
c_TOFD_3->Write();

TCanvas* c_TOFD_4  = new TCanvas();
c_TOFD_4->cd(1);
h_TOFD_4->Draw("hist");
c_TOFD_4->SaveAs("figures/TOFD_4.png");
c_TOFD_4->Write();

h_Q1_Q2->Write();
h_TOFD_1->Write();
h_TOFD_2->Write();
h_TOFD_3->Write();
h_TOFD_4->Write();

outFile->Close();

}
