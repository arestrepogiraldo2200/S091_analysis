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

// Correlation
TH2D* h_TTT10_MWPC_X_corr;
TH2D* h_FI33_FI32_X;
TH2D* h_FI31_FI32_X;
TH2D* h_FI33_FI31_combined_FI32_X;
TH2D* h_FI33_FI31_combined_FI32_E;

// Correlated area
TH2D* h_TTT10_MWPC_X_area;
TH2D* h_FIB31_FIB33_combined_FIB32_area;
TH2D* h_FIB31_FIB33_combined_FIB32_energy;

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
vector<R3BFiberMAPMTHitData> FIB32_corr_single_hit;
vector<R3BFiberMAPMTHitData> FIB33_corr_single_hit;
vector<R3BFiberMAPMTHitData> FIB31_corr_single_hit;

vector<R3BFiberMAPMTHitData> FI31_corr_hit;
vector<R3BFiberMAPMTHitData> FI32_corr_hit;
vector<R3BFiberMAPMTHitData> FI33_corr_hit;

void Fib31_Fib33_merging_X_POS_ENERGY(){

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

    h_FI33_FI31_combined_FI32_E = new TH2D("Energy FIB31&FIB33-FIB32 correlation","Energy FIB31&FIB33-FIB32 correlation", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FI33_FI31_combined_FI32_E->GetXaxis()->SetTitle("FIB32 energy [cm]");
    h_FI33_FI31_combined_FI32_E->GetYaxis()->SetTitle("FIB31&FIB33 energy [cm]");

    h_FIB31_FIB33_combined_FIB32_area = new TH2D("X position FIB31&FIB33-FIB32 combined correlation area","X position FIB31&FIB33-FIB32 combined correlation area", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FIB31_FIB33_combined_FIB32_area->GetXaxis()->SetTitle("FIB32 X-position [cm]" );
    h_FIB31_FIB33_combined_FIB32_area->GetYaxis()->SetTitle("FIB31&FIB31 X-position [cm]");

    h_FI33_FI32_X = new TH2D("X position FIB33-FIB32 correlation area","X position FIB33-FIB32 correlation area", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FI33_FI32_X->GetXaxis()->SetTitle("FIB32 X-position [cm]" );
    h_FI33_FI32_X->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    h_FI31_FI32_X = new TH2D("X position FIB31-FIB32 correlation area","X position FIB33-FIB32 correlation area", fiber_bins, -26.4, 26.4, 300, -52.8, 52.8);
    h_FI31_FI32_X->GetXaxis()->SetTitle("FIB32 X-position [cm]" );
    h_FI31_FI32_X->GetYaxis()->SetTitle("FIB31 X-position [cm]");

    h_MWPC_FI31_FI33_X_single_corr = new TH2D("X position MWPC and FIB31-FIB33 correlation single hit","X position MWPC and FIB31-FIB33 correlation single hit", mwpc1_bins, -150, 0, 300, -52.8, 52.8);
    h_MWPC_FI31_FI33_X_single_corr->GetXaxis()->SetTitle("MWPC X-position [mm]" );
    h_MWPC_FI31_FI33_X_single_corr->GetYaxis()->SetTitle("FIB31&FIB33 X-position [cm]");


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

    /////////////////////////////////////////////////////////////////////////////////////////

    TGraph* gFiberTrajectory = new TGraph();
    gFiberTrajectory->SetTitle("Fiber trajectory X-Z;X position [cm];Z position [cm]");
    gFiberTrajectory->SetMarkerStyle(20);
    gFiberTrajectory->SetMarkerSize(1.0);
    gFiberTrajectory->SetLineColor(kBlue);
    gFiberTrajectory->SetLineWidth(2);

    TH2D* hFiberTrajectoryFrame = new TH2D("hFiberTrajectoryFrame", "Fiber trajectory X-Z;X position [cm];Z position [cm]", 10, -20, 20, 10, -10, 130);
    hFiberTrajectoryFrame->SetStats(false);

    const double minFiberXRange = 20.0;
    const double minFiberZRange = 130.0;

    double fiberXMin = 1e9;
    double fiberXMax = -1e9;
    double fiberZMin = 1e9;
    double fiberZMax = -1e9;

    int fiberColors[] = {kRed, kGreen+2, kMagenta, kCyan, kOrange+7, kBlack};
    const int nFiberColors = sizeof(fiberColors)/sizeof(fiberColors[0]);
    int fiberTrajectoryColorIndex = 0;
    bool fiberTrajectoryDrawn = false;

    TCanvas* cFiberTrajectory = new TCanvas("cFiberTrajectory", "FIB fiber trajectory", 800, 600);

    /////////////////////////////////////////////////////////////////////////////////////////

	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in the merged tree: " << nevent << endl;

    int good_events = 0;

	for(Long64_t e=0; e < 10000; e++)
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

        FIB32_corr_single_hit.clear();
        FIB33_corr_single_hit.clear();
        FIB31_corr_single_hit.clear();

        FI31_corr_hit.clear();
        FI32_corr_hit.clear();
        FI33_corr_hit.clear();

        // Entries of LOS and TOFD per event in the chain of trees
        tpat             = header->GetTpat();
        Ttt10HitPerEvent = hitTTT10->GetEntriesFast();
        MwpcHitPerEvent  = hitMWPC->GetEntriesFast();

        Fi31HitPerEvent = hitFI31->GetEntriesFast();
        Fi32HitPerEvent = hitFI32->GetEntriesFast();
        Fi33HitPerEvent = hitFI33->GetEntriesFast();
	    TofdHitPerEvent = hitTOFD->GetEntriesFast();

        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        MwpcHitPerEvent   = hitMWPC->GetEntriesFast();


        /////////////////////////////////////////////////////////////////////////////////////////

        bool b_tofdq = false;

        if( TofdHitPerEvent > 0 ){

            hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
            
            for(int p = 0; p < TofdHitPerEvent; p++){

                hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
                double Eloss = hitTofd[p]->GetEloss();
                double TOFD_X = hitTofd[p]->GetX();
                h_TOFD_X->Fill(TOFD_X);

                if(Eloss<5 && Eloss>7) {
                    b_tofdq=true;                          // If at least one of the hits is charge 6 
                }
            }
        }
        
        if (b_tofdq == true || b_tofdq == false){

            /////////////////////////////////////////////////////////////////////////////////////////
            if( Ttt10HitPerEvent > 0 ){

                hitTtt10 = new R3BTttxHitData*[Ttt10HitPerEvent];

                for(int p = 0; p < Ttt10HitPerEvent; p++){

                    hitTtt10[p] =(R3BTttxHitData*)hitTTT10->At(p);
                    double TTT10_X = hitTtt10[p]->GetX();
                    h_TTT10_X->Fill(TTT10_X);                                               // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////
            double MWPC_X;

            if( MwpcHitPerEvent > 0 ){

                hitMwpc = new R3BMwpcHitData*[MwpcHitPerEvent];
                
                for(int p = 0; p < MwpcHitPerEvent; p++){

                    hitMwpc[p] =(R3BMwpcHitData*)hitMWPC->At(p);
                    MWPC_X = hitMwpc[p]->GetX();
                    h_MWPC_X->Fill(MWPC_X);                                               // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                    /////////////////////////////////////////////////////////////////////////////////////////

                    if( Ttt10HitPerEvent > 0 ){

                        hitTtt10 = new R3BTttxHitData*[Ttt10HitPerEvent];
                        
                        for(int q = 0; q < Ttt10HitPerEvent; q++){

                            hitTtt10[q] =(R3BTttxHitData*)hitTTT10->At(q);
                            double TTT10_X = hitTtt10[q]->GetX();

                            h_TTT10_MWPC_X->Fill(MWPC_X, TTT10_X);                         // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_TTT10_MWPC_X_corr->Fill(MWPC_X, TTT10_X);                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            if ( (TTT10_X<0.6*MWPC_X+47)  &&  (TTT10_X>0.6*MWPC_X+37) ){

                                h_TTT10_MWPC_X_area->Fill(MWPC_X, TTT10_X);                // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                TTT10_corr_hit.push_back(hitTtt10[q]);
                                MWPC_corr_hit.push_back(hitMwpc[p]);
                            }
                        }
                    } 
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////
            if( Fi32HitPerEvent > 0 ){

                hitFi32 = new R3BFiberMAPMTHitData*[Fi32HitPerEvent];
                double hitFI32_X, hitFI32_E;

                for(int p = 0; p < Fi32HitPerEvent; p++){

                    hitFi32[p] =(R3BFiberMAPMTHitData*)hitFI32->At(p);
                    auto* hit32 = static_cast<R3BFiberMAPMTHitData*>(hitFI32->At(p));

                    hitFI32_X = hitFi32[p]->GetX();
                    hitFI32_E = hitFi32[p]->GetEloss();
                    h_FI32_X->Fill(hitFI32_X);                                                            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Correlation with FIB31
                    if( Fi31HitPerEvent > 0 ){

                        hitFi31 = new R3BFiberMAPMTHitData*[Fi31HitPerEvent];
                        double hitFI31_X, hitFI31_E;

                        for(int q = 0; q < Fi31HitPerEvent; q++){

                            hitFi31[q] =(R3BFiberMAPMTHitData*)hitFI31->At(q);
                            auto* hit31 = static_cast<R3BFiberMAPMTHitData*>(hitFI31->At(q));

                            hitFI31_X = hitFi31[q]->GetX()-26.4+5;      
                            hitFI31_E = hitFi31[q]->GetEloss();      

                            h_FI31_X->Fill(hitFI31_X);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI33_FI31_combined_FI32_X->Fill(hitFI32_X, hitFI31_X);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI33_FI31_combined_FI32_E->Fill(hitFI32_E, hitFI31_E);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI31_FI32_X->Fill(hitFI32_X, hitFI31_X);                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FICOMB_X->Fill(hitFI31_X);                                                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            // Correlation condition
                            if ( (hitFI31_X>1*hitFI32_X-3) && (hitFI31_X<1*hitFI32_X+5) ){  

                                h_FIB31_FIB33_combined_FIB32_area->Fill( hitFI32_X , hitFI31_X);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                FI31_corr_hit.push_back(*hit31);
                                FI32_corr_hit.push_back(*hit32);
                            }
                        }
                    }

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Correlation with FIB33
                    if( Fi33HitPerEvent > 0 ){

                        hitFi33 = new R3BFiberMAPMTHitData*[Fi33HitPerEvent];
                        double hitFI33_X, hitFI33_E;

                        for(int q = 0; q < Fi33HitPerEvent; q++){

                            hitFi33[q] =(R3BFiberMAPMTHitData*)hitFI33->At(q);
                            auto* hit33 = static_cast<R3BFiberMAPMTHitData*>(hitFI33->At(q));

                            hitFI33_X = hitFi33[q]->GetX()+26.4-5;      
                            hitFI33_E = hitFi33[q]->GetEloss();      

                            h_FI33_X->Fill(hitFI33_X);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI33_FI31_combined_FI32_X->Fill(hitFI32_X, hitFI33_X);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI33_FI31_combined_FI32_E->Fill(hitFI32_E, hitFI33_E);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI33_FI32_X->Fill(hitFI32_X, hitFI33_X);                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FICOMB_X->Fill(hitFI33_X);                                                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            // Correlation condition
                            if ( (hitFI33_X>1*hitFI32_X-3) && (hitFI33_X<1*hitFI32_X+5) ){  

                                h_FIB31_FIB33_combined_FIB32_area->Fill( hitFI32_X , hitFI33_X);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                FI33_corr_hit.push_back(*hit33);
                                FI32_corr_hit.push_back(*hit32);
                            }
                        }
                    }  
                }
            }

            /////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////
            /////////////////////////////////////////////////////////////////////////////////////////

            if (TTT10_corr_hit.size()>0 && MWPC_corr_hit.size()>0 && FI32_corr_hit.size()>0 && (FI33_corr_hit.size()>0 || FI31_corr_hit.size()>0) ){   // THERE COULD BE HIT IN FI31 AND NOT IN FI33

                ////////////////////////////////////////////////////////////////////////////////////////////////////
                // THIS PART FILLS HISTOGRAMS USING THE CORRELATED REGIONS ONLY AND TAKING THE SMALLEST MULTIPLICITY
                ////////////////////////////////////////////////////////////////////////////////////////////////////

                double TTT10_X_CORR;
                double  MWPC_X_CORR;
                double FIB32_X_CORR;
                double FIB33_X_CORR;
                double FIB31_X_CORR;

                int smallest = std::min((int)FI31_corr_hit.size(), std::min(std::min((int)TTT10_corr_hit.size(), (int)MWPC_corr_hit.size()), std::min((int)FI32_corr_hit.size(), (int)FI33_corr_hit.size())));

                // TTT10 vs FIB33 ///////////////////////////////////////////////////////////////////////
                if (FI33_corr_hit.size()>0){
                    for(int i = 0; i < smallest; i++){

                        TTT10_X_CORR = TTT10_corr_hit[i]->GetX();
                        FIB33_X_CORR = FI33_corr_hit[i].GetX();
                        h_TTT10_FI33_X_corr->Fill(TTT10_X_CORR, FIB33_X_CORR);                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    }
                }

                // TTT10 vs FIB32 ///////////////////////////////////////////////////////////////////////
                for(int i = 0; i < smallest; i++){

                    TTT10_X_CORR = TTT10_corr_hit[i]->GetX();
                    FIB32_X_CORR = FI32_corr_hit[i].GetX();
                    h_TTT10_FI32_X_corr->Fill(TTT10_X_CORR, FIB32_X_CORR);                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                }

                // MWPC vs FIB33 ///////////////////////////////////////////////////////////////////////
                if (FI33_corr_hit.size()>0){
                    for(int i = 0; i < smallest; i++){

                        MWPC_X_CORR = MWPC_corr_hit[i]->GetX();
                        FIB33_X_CORR = FI33_corr_hit[i].GetX();
                        h_MWPC_FI33_X_corr->Fill(MWPC_X_CORR, FIB33_X_CORR);                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    }
                }

                // MWPC vs FIB32 ///////////////////////////////////////////////////////////////////////
                for(int i = 0; i < smallest; i++){

                    MWPC_X_CORR = MWPC_corr_hit[i]->GetX();
                    FIB32_X_CORR = FI32_corr_hit[i].GetX();
                    h_MWPC_FI32_X_corr->Fill(MWPC_X_CORR, FIB32_X_CORR);                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                }

                /////////////////////////////////////////////////////////////////////////////////////
                // THIS PART FILLS HISTOGRAMS USING THE CORRELATED REGIONS ONLY AND SINGLE HIT FOR IT
                /////////////////////////////////////////////////////////////////////////////////////

                double largestE_TTT10 = -1;
                double largestE_MWPC  = -1;
                double largestE_FIB32 = -1;
                double largestE_FIB33 = -1;
                double largestE_FIB31 = -1;
                double largestE_FIB33_only = -1;
                double largestE_FIB31_only = -1;

                double largestE_index_TTT10;
                double largestE_index_MWPC ;
                double largestE_index_FIB32;
                double largestE_index_FIB33;
                double largestE_index_FIB31;

                /////////////////////////////////////////////////////////////////////////////////////////
                for (int i=0; i < TTT10_corr_hit.size(); i++){

                    double energyttt10 = TTT10_corr_hit[i]->GetEnergy();

                    if (energyttt10 > largestE_TTT10){
                        largestE_TTT10 = energyttt10;
                        largestE_index_TTT10 = i;
                    }
                }
                TTT10_corr_single_hit.push_back(TTT10_corr_hit[largestE_index_TTT10]);

                /////////////////////////////////////////////////////////////////////////////////////////
                for (int i=0; i < FI32_corr_hit.size(); i++){

                    double energyfib32 = FI32_corr_hit[i].GetEloss();

                    if (energyfib32 > largestE_FIB32){
                        largestE_FIB32 = energyfib32;
                        largestE_index_FIB32 = i;
                    }
                }
                FIB32_corr_single_hit.push_back(FI32_corr_hit[largestE_index_FIB32]);

                /////////////////////////////////////////////////////////////////////////////////////////
                if (FI33_corr_hit.size()>0){
                    for (int i=0; i < FI33_corr_hit.size(); i++){

                        double energyfib33 = FI33_corr_hit[i].GetEloss();

                        if (energyfib33 > largestE_FIB33){
                            largestE_FIB33 = energyfib33;
                            largestE_index_FIB33 = i;
                        }
                    }
                    FIB33_corr_single_hit.push_back(FI33_corr_hit[largestE_index_FIB33]);
                }

                /////////////////////////////////////////////////////////////////////////////////////////
                if (FI31_corr_hit.size()>0){
                    for (int i=0; i < FI31_corr_hit.size(); i++){

                        double energyfib31 = FI31_corr_hit[i].GetEloss();

                        if (energyfib31 > largestE_FIB31){
                            largestE_FIB31 = energyfib31;
                            largestE_index_FIB31 = i;
                        }
                    }
                    FIB31_corr_single_hit.push_back(FI31_corr_hit[largestE_index_FIB31]);
                }

                /////////////////////////////////////////////////////////////////////////////////////////

                //cout << "TTT10: " << TTT10_corr_single_hit.size() << " FIB33: " << FIB33_corr_single_hit.size() << " FIB32: " << FIB32_corr_single_hit.size() << " MWPC: " << MWPC_corr_hit.size() << " FIB31: " << FIB31_corr_single_hit.size() << endl;

                if (FI33_corr_hit.size()>0){
                    h_TTT10_FI33_X_single_corr->Fill(TTT10_corr_single_hit[0]->GetX(), FIB33_corr_single_hit[0].GetX());
                }
                h_TTT10_FI32_X_single_corr->Fill(TTT10_corr_single_hit[0]->GetX(), FIB32_corr_single_hit[0].GetX());
                h_MWPC_FI32_X_single_corr->Fill(MWPC_corr_hit[0]->GetX(), FIB32_corr_single_hit[0].GetX());                                 // !!!!  JUST TAKING FIRST ENTRY OF MWPC 

                if (FI33_corr_hit.size()>0){
                    h_MWPC_FI31_FI33_X_single_corr->Fill(MWPC_corr_hit[0]->GetX(), FIB33_corr_single_hit[0].GetX()+26.4-5);                 // !!!!  JUST TAKING FIRST ENTRY OF MWPC 
                }
                if (FI31_corr_hit.size()>0){
                    h_MWPC_FI31_FI33_X_single_corr->Fill(MWPC_corr_hit[0]->GetX(), FIB31_corr_single_hit[0].GetX()-26.4+5);                 // !!!!  JUST TAKING FIRST ENTRY OF MWPC 
                }



                // /////////////////////////////////////////////////////////////////////////////////////////
                // // Case of hits in all three fiber detectors

                // if ( FI31_corr_hit.size()>0 && FI32_corr_hit.size()>0 && FI33_corr_hit.size()>0 ){

                //     //if ( (FIB33_corr_single_hit[0].GetX() <= 5) && (FIB31_corr_single_hit[0].GetX() >= -5) && (FIB32_corr_single_hit[0].GetX() >= -5) && (FIB32_corr_single_hit[0].GetX() <= 5) ){

                //         cout << "FIB32 X: " << FIB32_corr_single_hit[0].GetX() << "   FIB33 X: " << FIB33_corr_single_hit[0].GetX() << "   FIB31 X: " << FIB31_corr_single_hit[0].GetX() << endl;

                //         int nPoints = gFiberTrajectory->GetN();
                //         double x0 = FIB32_corr_single_hit[0].GetX();
                //         double z0 = 0.0;
                //         double x1 = FIB33_corr_single_hit[0].GetX()+26.4-5;
                //         double z1 = 100.0;
                //         double x2 = FIB31_corr_single_hit[0].GetX()-26.4+5;
                //         double z2 = 117.5;

                //         gFiberTrajectory->SetPoint(nPoints++, x0, z0);
                //         gFiberTrajectory->SetPoint(nPoints++, x1, z1);
                //         gFiberTrajectory->SetPoint(nPoints++, x2, z2);

                //         fiberXMin = std::min(fiberXMin, std::min(x0, std::min(x1, x2)));
                //         fiberXMax = std::max(fiberXMax, std::max(x0, std::max(x1, x2)));
                //         fiberZMin = std::min(fiberZMin, std::min(z0, std::min(z1, z2)));
                //         fiberZMax = std::max(fiberZMax, std::max(z0, std::max(z1, z2)));

                //         TGraph* gCurrentLine = new TGraph(3);
                //         gCurrentLine->SetPoint(0, x0, z0);
                //         gCurrentLine->SetPoint(1, x1, z1);
                //         gCurrentLine->SetPoint(2, x2, z2);
                //         //gCurrentLine->SetLineColor(fiberColors[fiberTrajectoryColorIndex % nFiberColors]);
                                            
                //         gCurrentLine->SetLineColor(kBlue);
                //         gCurrentLine->SetMarkerColor(kBlack);
                //         gCurrentLine->SetMarkerStyle(20);
                //         gCurrentLine->SetMarkerSize(2.5);
                //         gCurrentLine->SetLineWidth(2);

                //         cFiberTrajectory->cd();
                //         if (fiberXMin < fiberXMax) {
                //             double xMargin = 5.0;
                //             double zMargin = 10.0;
                //             double xMin = fiberXMin - xMargin;
                //             double xMax = fiberXMax + xMargin;
                //             double xSpan = xMax - xMin;
                //             if (xSpan < minFiberXRange) {
                //                 double xCenter = 0.5*(xMin + xMax);
                //                 xMin = xCenter - 0.5*minFiberXRange;
                //                 xMax = xCenter + 0.5*minFiberXRange;
                //             }
                //             double zMin = fiberZMin - zMargin;
                //             double zMax = fiberZMax + zMargin;
                //             double zSpan = zMax - zMin;
                //             if (zSpan < minFiberZRange) {
                //                 double zCenter = 0.5*(zMin + zMax);
                //                 zMin = zCenter - 0.5*minFiberZRange;
                //                 zMax = zCenter + 0.5*minFiberZRange;
                //             }
                //             hFiberTrajectoryFrame->GetXaxis()->SetLimits(xMin, xMax);
                //             hFiberTrajectoryFrame->GetYaxis()->SetRangeUser(zMin, zMax);
                //         }
                //         if (!fiberTrajectoryDrawn) {
                //             hFiberTrajectoryFrame->Draw();
                //             gFiberTrajectory->Draw("P same");
                //             fiberTrajectoryDrawn = true;
                //         } else {
                //             hFiberTrajectoryFrame->Draw("axis same");
                //             gFiberTrajectory->Draw("P same");
                //         }
                //         gCurrentLine->Draw("PL same");
                //         cFiberTrajectory->Modified();
                //         cFiberTrajectory->Update();

                //         fiberTrajectoryColorIndex++;
                //     //}
                // }
                // /////////////////////////////////////////////////////////////////////////////////////////
                // // Case of hits in FI33 and FI32 detectors only

                // if ( FI31_corr_hit.size()==0 && FI32_corr_hit.size()>0 && FI33_corr_hit.size()>0 ){

                //     //if ( (FIB33_corr_single_hit[0].GetX() <= 5) && (FIB31_corr_single_hit[0].GetX() >= -5) && (FIB32_corr_single_hit[0].GetX() >= -5) && (FIB32_corr_single_hit[0].GetX() <= 5) ){

                //         cout << "FIB32 X: " << FIB32_corr_single_hit[0].GetX() << "   FIB33 X: " << FIB33_corr_single_hit[0].GetX() << endl;

                //         double x0 = FIB32_corr_single_hit[0].GetX();
                //         double z0 = 0.0;
                //         double x1 = FIB33_corr_single_hit[0].GetX()+26.4-5;
                //         double z1 = 100.0;

                //         TGraph* gCurrentLine = new TGraph(2);
                //         gCurrentLine->SetPoint(0, x0, z0);
                //         gCurrentLine->SetPoint(1, x1, z1);
                //         //gCurrentLine->SetLineColor(fiberColors[fiberTrajectoryColorIndex % nFiberColors]);

                //         fiberXMin = std::min(fiberXMin, std::min(x0, x1));
                //         fiberXMax = std::max(fiberXMax, std::max(x0, x1));
                //         fiberZMin = std::min(fiberZMin, std::min(z0, z1));
                //         fiberZMax = std::max(fiberZMax, std::max(z0, z1));

                //         gCurrentLine->SetLineColor(kRed);
                //         gCurrentLine->SetMarkerColor(kBlack);
                //         gCurrentLine->SetMarkerStyle(20);
                //         gCurrentLine->SetMarkerSize(2.5);
                //         gCurrentLine->SetLineWidth(2);

                //         cFiberTrajectory->cd();
                //         if (fiberXMin < fiberXMax) {
                //             double xMargin = 5.0;
                //             double zMargin = 10.0;
                //             double xMin = fiberXMin - xMargin;
                //             double xMax = fiberXMax + xMargin;
                //             double xSpan = xMax - xMin;
                //             if (xSpan < minFiberXRange) {
                //                 double xCenter = 0.5*(xMin + xMax);
                //                 xMin = xCenter - 0.5*minFiberXRange;
                //                 xMax = xCenter + 0.5*minFiberXRange;
                //             }
                //             double zMin = fiberZMin - zMargin;
                //             double zMax = fiberZMax + zMargin;
                //             double zSpan = zMax - zMin;
                //             if (zSpan < minFiberZRange) {
                //                 double zCenter = 0.5*(zMin + zMax);
                //                 zMin = zCenter - 0.5*minFiberZRange;
                //                 zMax = zCenter + 0.5*minFiberZRange;
                //             }
                //             hFiberTrajectoryFrame->GetXaxis()->SetLimits(xMin, xMax);
                //             hFiberTrajectoryFrame->GetYaxis()->SetRangeUser(zMin, zMax);
                //         }
                //         if (!fiberTrajectoryDrawn) {
                //             hFiberTrajectoryFrame->Draw();
                //             fiberTrajectoryDrawn = true;
                //         } else {
                //             hFiberTrajectoryFrame->Draw("axis same");
                //         }
                //         gCurrentLine->Draw("PL same");
                //         cFiberTrajectory->Modified();
                //         cFiberTrajectory->Update();

                //         fiberTrajectoryColorIndex++;
                //     //}
                // }

                // /////////////////////////////////////////////////////////////////////////////////////////
                // // Case of hits in FI31 and FI32 detectors only


                // if ( FI31_corr_hit.size()>0 && FI32_corr_hit.size()>0 && FI33_corr_hit.size()==0 ){

                //     //if ( (FIB33_corr_single_hit[0].GetX() <= 5) && (FIB31_corr_single_hit[0].GetX() >= -5) && (FIB32_corr_single_hit[0].GetX() >= -5) && (FIB32_corr_single_hit[0].GetX() <= 5) ){

                //         cout << "FIB32 X:" << FIB32_corr_single_hit[0].GetX() <<  "   FIB31 X:" << FIB31_corr_single_hit[0].GetX() << endl;

                //         double x0 = FIB32_corr_single_hit[0].GetX();
                //         double z0 = 0.0;
                //         double x1 = FIB31_corr_single_hit[0].GetX()-26.4+5;
                //         double z1 = 117.5;

                //         TGraph* gCurrentLine = new TGraph(2);
                //         gCurrentLine->SetPoint(0, x0, z0);
                //         gCurrentLine->SetPoint(1, x1, z1);
                //         //gCurrentLine->SetLineColor(fiberColors[fiberTrajectoryColorIndex % nFiberColors]);        

                //         fiberXMin = std::min(fiberXMin, std::min(x0, x1));
                //         fiberXMax = std::max(fiberXMax, std::max(x0, x1));
                //         fiberZMin = std::min(fiberZMin, std::min(z0, z1));
                //         fiberZMax = std::max(fiberZMax, std::max(z0, z1));

                //         gCurrentLine->SetLineColor(kGreen);
                //         gCurrentLine->SetMarkerColor(kBlack);
                //         gCurrentLine->SetMarkerStyle(20);
                //         gCurrentLine->SetMarkerSize(2.5);
                //         gCurrentLine->SetLineWidth(2);

                //         cFiberTrajectory->cd();
                //         if (fiberXMin < fiberXMax) {
                //             double xMargin = 5.0;
                //             double zMargin = 10.0;
                //             double xMin = fiberXMin - xMargin;
                //             double xMax = fiberXMax + xMargin;
                //             double xSpan = xMax - xMin;
                //             if (xSpan < minFiberXRange) {
                //                 double xCenter = 0.5*(xMin + xMax);
                //                 xMin = xCenter - 0.5*minFiberXRange;
                //                 xMax = xCenter + 0.5*minFiberXRange;
                //             }
                //             double zMin = fiberZMin - zMargin;
                //             double zMax = fiberZMax + zMargin;
                //             double zSpan = zMax - zMin;
                //             if (zSpan < minFiberZRange) {
                //                 double zCenter = 0.5*(zMin + zMax);
                //                 zMin = zCenter - 0.5*minFiberZRange;
                //                 zMax = zCenter + 0.5*minFiberZRange;
                //             }
                //             hFiberTrajectoryFrame->GetXaxis()->SetLimits(xMin, xMax);
                //             hFiberTrajectoryFrame->GetYaxis()->SetRangeUser(zMin, zMax);
                //         }
                //         if (!fiberTrajectoryDrawn) {
                //             hFiberTrajectoryFrame->Draw();
                //             fiberTrajectoryDrawn = true;
                //         } else {
                //             hFiberTrajectoryFrame->Draw("axis same");
                //         }
                //         gCurrentLine->Draw("PL same");
                //         cFiberTrajectory->Modified();
                //         cFiberTrajectory->Update();

                //         fiberTrajectoryColorIndex++;
                //     //}
                // }

                good_events++;
            }
        }
	}

cout << "Good events: " << good_events << endl;

TCanvas* c_TTT10_X  = new TCanvas();
c_TTT10_X ->cd(1);
h_TTT10_X ->Draw("hist");

TCanvas* c_MWPC_X  = new TCanvas();
c_MWPC_X ->cd(1);
h_MWPC_X ->Draw("hist");

TCanvas* c_TTT10_MWPC_X = new TCanvas();
c_TTT10_MWPC_X->cd(1);
h_TTT10_MWPC_X->Draw("colz");

TCanvas* c_FI31_X  = new TCanvas();
c_FI31_X->cd(1);
h_FI31_X->Draw("hist");

TCanvas* c_FI32_X  = new TCanvas();
c_FI32_X->cd(1);
h_FI32_X->Draw("hist");

TCanvas* c_FI33_X  = new TCanvas();
c_FI33_X->cd(1);
h_FI33_X->Draw("hist");

TCanvas* c_FICOMB_X  = new TCanvas();
c_FICOMB_X->cd(1);
h_FICOMB_X->Draw("hist");

TCanvas* c_TOFD_X  = new TCanvas();
c_TOFD_X->cd(1);
h_TOFD_X->Draw("hist");

// ///////////////////////////////////////////////////////////////////

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

double m3 = 1.0;
double b3 = 2;
double x1_3 = h_FI33_FI31_combined_FI32_X->GetXaxis()->GetXmin();
double x2_3 = h_FI33_FI31_combined_FI32_X->GetXaxis()->GetXmax();
double y1_3 = m3*x1_3 + b3;
double y2_3 = m3*x2_3 + b3;

TCanvas *c3 = new TCanvas();
h_FI33_FI31_combined_FI32_X->Draw("colz");

TLine *line3 = new TLine(x1_3, y1_3, x2_3, y2_3);
line3->SetLineColor(kRed);
line3->SetLineWidth(2);
line3->Draw("same");

TF1 *f3 = new TF1("line3", "[0]*x + [1]", x1_3, x2_3);
f3->SetParameters(m3, b3);
f3->SetLineColor(kRed);
f3->SetLineWidth(1);

h_FI33_FI31_combined_FI32_X->Draw("colz");
f3->Draw("same");

// ///////////////////////////////////////////////////////////////////

TCanvas* c_FI33_FI31_combined_FI32_E = new TCanvas();
c_FI33_FI31_combined_FI32_E->cd(1);
h_FI33_FI31_combined_FI32_E->Draw("colz");

TCanvas* c_TTT10_MWPC_X_area = new TCanvas();
c_TTT10_MWPC_X_area->cd(1);
h_TTT10_MWPC_X_area->Draw("colz");

// ///////////////////////////////////////////////////////////////////

TCanvas* c_TTT10_FI33_X_corr = new TCanvas();
c_TTT10_FI33_X_corr->cd(1);
h_TTT10_FI33_X_corr->Draw("colz");

TCanvas* c_TTT10_FI32_X_corr = new TCanvas();
c_TTT10_FI32_X_corr->cd(1);
h_TTT10_FI32_X_corr->Draw("colz");

TCanvas* c_MWPC_FI33_X_corr = new TCanvas();
c_MWPC_FI33_X_corr->cd(1);
h_MWPC_FI33_X_corr->Draw("colz");

TCanvas* c_MWPC_FI32_X_corr = new TCanvas();
c_MWPC_FI32_X_corr->cd(1);
h_MWPC_FI32_X_corr->Draw("colz");

// ///////////////////////////////////////////////////////////////////

TCanvas* c_TTT10_FI33_X_single_corr = new TCanvas();
c_TTT10_FI33_X_single_corr->cd(1);
h_TTT10_FI33_X_single_corr->Draw("colz");

TCanvas* c_TTT10_FI32_X_single_corr = new TCanvas();
c_TTT10_FI32_X_single_corr->cd(1);
h_TTT10_FI32_X_single_corr->Draw("colz");

TCanvas* c_MWPC_FI33_X_single_corr = new TCanvas();
c_MWPC_FI33_X_single_corr->cd(1);
h_MWPC_FI33_X_single_corr->Draw("colz");

TCanvas* c_MWPC_FI32_X_single_corr = new TCanvas();
c_MWPC_FI32_X_single_corr->cd(1);
h_MWPC_FI32_X_single_corr->Draw("colz");

TCanvas* c_FI33_FI32_X = new TCanvas();
c_FI33_FI32_X->cd(1);
h_FI33_FI32_X->Draw("colz");

TCanvas* c_FI31_FI32_X = new TCanvas();
c_FI31_FI32_X->cd(1);
h_FI31_FI32_X->Draw("colz");

TCanvas* c_FIB31_FIB33_combined_FIB32_area = new TCanvas();
c_FIB31_FIB33_combined_FIB32_area->cd(1);
h_FIB31_FIB33_combined_FIB32_area->Draw("colz");

TCanvas* c_MWPC_FI31_FI33_X_single_corr = new TCanvas();
c_MWPC_FI31_FI33_X_single_corr->cd(1);
h_MWPC_FI31_FI33_X_single_corr->Draw("colz");

}