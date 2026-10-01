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

// Correlation
TH2D* h_TTT10_MWPC_X_corr;
TH2D* h_FI33_FI32_X;
TH2D* h_FI31_FI32_X;
TH2D* h_FI33_FI31_combined_FI32_X;

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

void Fib31_Fib33_merging_X_HITS_PER_EVENT(){

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

    /////////////////////////////////////////////////////////////////////////////////////////
    // TGraph* gFiberTrajectory = new TGraph();
    // gFiberTrajectory->SetTitle("Fiber trajectory X-Z;X position [cm];Z position [cm]");
    // gFiberTrajectory->SetMarkerStyle(20);
    // gFiberTrajectory->SetMarkerSize(1.0);
    // gFiberTrajectory->SetLineColor(kBlue);
    // gFiberTrajectory->SetLineWidth(2);

    // TH2D* hFiberTrajectoryFrame = new TH2D("hFiberTrajectoryFrame", "Fiber trajectory X-Z;X position [cm];Z position [cm]", 10, -20, 20, 10, -10, 130);
    // hFiberTrajectoryFrame->SetStats(false);

    // const double minFiberXRange = 20.0;
    // const double minFiberZRange = 130.0;

    // double fiberXMin = 1e9;
    // double fiberXMax = -1e9;
    // double fiberZMin = 1e9;
    // double fiberZMax = -1e9;

    // int fiberColors[] = {kRed, kGreen+2, kMagenta, kCyan, kOrange+7, kBlack};
    // const int nFiberColors = sizeof(fiberColors)/sizeof(fiberColors[0]);
    // int fiberTrajectoryColorIndex = 0;
    // bool fiberTrajectoryDrawn = false;

    // TCanvas* cFiberTrajectory = new TCanvas("cFiberTrajectory", "FIB fiber trajectory", 800, 600);

    /////////////////////////////////////////////////////////////////////////////////////////

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

        if (Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts1++; }
        if (MwpcHitPerEvent == 1 && Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts2++; }
        if (LosHitPerEvent  == 1 && MwpcHitPerEvent == 1 && Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts3++; }
        if (TofdHitPerEvent == 1 && LosHitPerEvent  == 1 && MwpcHitPerEvent == 1 && Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts4++; }
        if ((TofdHitPerEvent == 4 || TofdHitPerEvent == 5) && LosHitPerEvent  == 1 && MwpcHitPerEvent == 1 && Fi32HitPerEvent == 1 && (Fi31HitPerEvent == 1 || Fi33HitPerEvent == 1) ){ evts5++; }

        /////////////////////////////////////////////////////////////////////////////////////////

        if ((TofdHitPerEvent == 4 || TofdHitPerEvent == 5) && LosHitPerEvent  == 1 && MwpcHitPerEvent == 1 && Fi32HitPerEvent == 1 && (Fi31HitPerEvent == 1 || Fi33HitPerEvent == 1) ){
        //if (Fi32HitPerEvent >= 1 && Fi31HitPerEvent >= 1 && Fi33HitPerEvent >= 1 ){

            bool b_tofdq = false;

            if( TofdHitPerEvent > 0 ){                         // This condition is overrun by the previous conditional 

                hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
                
                for(int p = 0; p < TofdHitPerEvent; p++){

                    hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
                    double Eloss = hitTofd[p]->GetEloss();
                    double TOFD_X = hitTofd[p]->GetX();
                    h_TOFD_X->Fill(TOFD_X);
                    
                    if(Eloss>5.5 && Eloss<6.5) {                   
                        b_tofdq=true;                          // If at least one of the hits is charge 6 
                    }
                }
            }
        
            if ( b_tofdq == true ){

                good_events++;

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
                    double hitFI32_X;
                    double hitFI32_T;

                    for(int p = 0; p < Fi32HitPerEvent; p++){

                        hitFi32[p] =(R3BFiberMAPMTHitData*)hitFI32->At(p);
                        auto* hit32 = static_cast<R3BFiberMAPMTHitData*>(hitFI32->At(p));

                        hitFI32_X = hitFi32[p]->GetX();
                        hitFI32_T = hitFi32[p]->GetTime_ns();

                        h_FI32_X->Fill(hitFI32_X);                                                            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                        //h_FI32_T->Fill(hitFI32_T);                                                            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                        /////////////////////////////////////////////////////////////////////////////////////////
                        // Correlation with FIB31
                        if( Fi31HitPerEvent > 0 ){

                            hitFi31 = new R3BFiberMAPMTHitData*[Fi31HitPerEvent];
                            double hitFI31_X, hitFI31_T;

                            for(int q = 0; q < Fi31HitPerEvent; q++){

                                hitFi31[q] =(R3BFiberMAPMTHitData*)hitFI31->At(q);
                                auto* hit31 = static_cast<R3BFiberMAPMTHitData*>(hitFI31->At(q));

                                hitFI31_X = hitFi31[q]->GetX()-26.4+5;    
                                //hitFI31_T = hitFi31[p]->GetTime_ns();   

                                h_FI31_X->Fill(hitFI31_X);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                //h_FI31_T->Fill(hitFI31_T);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                h_FI33_FI31_combined_FI32_X->Fill(hitFI32_X, hitFI31_X);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                h_FI31_FI32_X->Fill(hitFI32_X, hitFI31_X);                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                h_FICOMB_X->Fill(hitFI31_X);                                                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                // Correlation condition
                                if ( (hitFI31_X>1*hitFI32_X-4) && (hitFI31_X<1*hitFI32_X+6) ){  

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
                            double hitFI33_X, hitFI33_T;

                            for(int q = 0; q < Fi33HitPerEvent; q++){

                                hitFi33[q] =(R3BFiberMAPMTHitData*)hitFI33->At(q);
                                auto* hit33 = static_cast<R3BFiberMAPMTHitData*>(hitFI33->At(q));

                                hitFI33_X = hitFi33[q]->GetX()+26.4-5;      
                                //hitFI33_T = hitFi33[p]->GetTime_ns();   

                                h_FI33_X->Fill(hitFI33_X);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                //h_FI33_T->Fill(hitFI33_T);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                h_FI33_FI31_combined_FI32_X->Fill(hitFI32_X, hitFI33_X);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                h_FI33_FI32_X->Fill(hitFI32_X, hitFI33_X);                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                                h_FICOMB_X->Fill(hitFI33_X);                                                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                // Correlation condition
                                if ( (hitFI33_X>1*hitFI32_X-4) && (hitFI33_X<1*hitFI32_X+6) ){  

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
                /////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////


                
                // CONDITIONAL FOR CORRELATED HITS IN BEFORE AND AFTER DETECTORS
                if (TTT10_corr_hit.size()>0 && MWPC_corr_hit.size()>0 && FI32_corr_hit.size()>0 && (FI33_corr_hit.size()>0 || FI31_corr_hit.size()>0) ){  

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Case of hits in all three fiber detectors

                    if ( FI31_corr_hit.size()>0 && FI32_corr_hit.size()>0 && FI33_corr_hit.size()>0 ){

                        //if ( (FI33_corr_hit[0].GetX() <= 5) && (FI31_corr_hit[0].GetX() >= -5) && (FI32_corr_hit[0].GetX() >= -5) && (FI32_corr_hit[0].GetX() <= 5) ){
                           
                            int nPoints = gFiberTrajectory->GetN();
                            double x0 = FI32_corr_hit[0].GetX();
                            double z0 = 0.0;
                            double x1 = FI33_corr_hit[0].GetX()+26.4-5;
                            double z1 = 100.0;
                            double x2 = FI31_corr_hit[0].GetX()-26.4+5;
                            double z2 = 117.5;

                            double triangle_area = TriangleArea(x0, z0, x1, z1, x2, z2);

                            cout << "FIB32 X: " << FI32_corr_hit[0].GetX() << "   FIB33 X: " << FI33_corr_hit[0].GetX() << "   FIB31 X: " << FI31_corr_hit[0].GetX() << " Area: "  <<  triangle_area << endl;

                            gFiberTrajectory->SetPoint(nPoints++, x0, z0);
                            gFiberTrajectory->SetPoint(nPoints++, x1, z1);
                            gFiberTrajectory->SetPoint(nPoints++, x2, z2);

                            fiberXMin = std::min(fiberXMin, std::min(x0, std::min(x1, x2)));
                            fiberXMax = std::max(fiberXMax, std::max(x0, std::max(x1, x2)));
                            fiberZMin = std::min(fiberZMin, std::min(z0, std::min(z1, z2)));
                            fiberZMax = std::max(fiberZMax, std::max(z0, std::max(z1, z2)));

                            TGraph* gCurrentLine = new TGraph(3);
                            gCurrentLine->SetPoint(0, x0, z0);
                            gCurrentLine->SetPoint(1, x1, z1);
                            gCurrentLine->SetPoint(2, x2, z2);
                            //gCurrentLine->SetLineColor(fiberColors[fiberTrajectoryColorIndex % nFiberColors]);
                            
                            if (triangle_area < 100){
                                gCurrentLine->SetLineColor(kRed);
                                collinear++;
                            } else {
                                gCurrentLine->SetLineColor(kBlue);
                                noncollinear++;
                            }

                            gCurrentLine->SetMarkerColor(kBlack);
                            gCurrentLine->SetMarkerStyle(20);
                            gCurrentLine->SetMarkerSize(2.5);
                            gCurrentLine->SetLineWidth(2);

                            cFiberTrajectory->cd();
                            if (fiberXMin < fiberXMax) {
                                double xMargin = 5.0;
                                double zMargin = 10.0;
                                double xMin = fiberXMin - xMargin;
                                double xMax = fiberXMax + xMargin;
                                double xSpan = xMax - xMin;
                                if (xSpan < minFiberXRange) {
                                    double xCenter = 0.5*(xMin + xMax);
                                    xMin = xCenter - 0.5*minFiberXRange;
                                    xMax = xCenter + 0.5*minFiberXRange;
                                }
                                double zMin = fiberZMin - zMargin;
                                double zMax = fiberZMax + zMargin;
                                double zSpan = zMax - zMin;
                                if (zSpan < minFiberZRange) {
                                    double zCenter = 0.5*(zMin + zMax);
                                    zMin = zCenter - 0.5*minFiberZRange;
                                    zMax = zCenter + 0.5*minFiberZRange;
                                }
                                hFiberTrajectoryFrame->GetXaxis()->SetLimits(xMin, xMax);
                                hFiberTrajectoryFrame->GetYaxis()->SetRangeUser(zMin, zMax);
                            }
                            if (!fiberTrajectoryDrawn) {
                                hFiberTrajectoryFrame->Draw();
                                gFiberTrajectory->Draw("P same");
                                fiberTrajectoryDrawn = true;
                            } else {
                                hFiberTrajectoryFrame->Draw("axis same");
                                gFiberTrajectory->Draw("P same");
                            }
                            gCurrentLine->Draw("PL same");
                            cFiberTrajectory->Modified();
                            cFiberTrajectory->Update();

                            fiberTrajectoryColorIndex++;
                        //}
                    }
                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Case of hits in FI33 and FI32 detectors only

                    if ( FI31_corr_hit.size()==0 && FI32_corr_hit.size()>0 && FI33_corr_hit.size()>0 ){

                        //if ( (FI33_corr_hit[0].GetX() <= 5) && (FI31_corr_hit[0].GetX() >= -5) && (FI32_corr_hit[0].GetX() >= -5) && (FI32_corr_hit[0].GetX() <= 5) ){

                            cout << "FIB32 X: " << FI32_corr_hit[0].GetX() << "   FIB33 X: " << FI33_corr_hit[0].GetX() << endl;

                            double x0 = FI32_corr_hit[0].GetX();
                            double z0 = 0.0;
                            double x1 = FI33_corr_hit[0].GetX()+26.4-5;
                            double z1 = 100.0;

                            TGraph* gCurrentLine = new TGraph(2);
                            gCurrentLine->SetPoint(0, x0, z0);
                            gCurrentLine->SetPoint(1, x1, z1);
                            //gCurrentLine->SetLineColor(fiberColors[fiberTrajectoryColorIndex % nFiberColors]);

                            fiberXMin = std::min(fiberXMin, std::min(x0, x1));
                            fiberXMax = std::max(fiberXMax, std::max(x0, x1));
                            fiberZMin = std::min(fiberZMin, std::min(z0, z1));
                            fiberZMax = std::max(fiberZMax, std::max(z0, z1));

                            gCurrentLine->SetLineColor(kCyan);
                            gCurrentLine->SetMarkerColor(kBlack);
                            gCurrentLine->SetMarkerStyle(20);
                            gCurrentLine->SetMarkerSize(2.5);
                            gCurrentLine->SetLineWidth(2);

                            cFiberTrajectory->cd();
                            if (fiberXMin < fiberXMax) {
                                double xMargin = 5.0;
                                double zMargin = 10.0;
                                double xMin = fiberXMin - xMargin;
                                double xMax = fiberXMax + xMargin;
                                double xSpan = xMax - xMin;
                                if (xSpan < minFiberXRange) {
                                    double xCenter = 0.5*(xMin + xMax);
                                    xMin = xCenter - 0.5*minFiberXRange;
                                    xMax = xCenter + 0.5*minFiberXRange;
                                }
                                double zMin = fiberZMin - zMargin;
                                double zMax = fiberZMax + zMargin;
                                double zSpan = zMax - zMin;
                                if (zSpan < minFiberZRange) {
                                    double zCenter = 0.5*(zMin + zMax);
                                    zMin = zCenter - 0.5*minFiberZRange;
                                    zMax = zCenter + 0.5*minFiberZRange;
                                }
                                hFiberTrajectoryFrame->GetXaxis()->SetLimits(xMin, xMax);
                                hFiberTrajectoryFrame->GetYaxis()->SetRangeUser(zMin, zMax);
                            }
                            if (!fiberTrajectoryDrawn) {
                                hFiberTrajectoryFrame->Draw();
                                fiberTrajectoryDrawn = true;
                            } else {
                                hFiberTrajectoryFrame->Draw("axis same");
                            }
                            gCurrentLine->Draw("PL same");
                            cFiberTrajectory->Modified();
                            cFiberTrajectory->Update();

                            fiberTrajectoryColorIndex++;
                        //}1
                    }

                    ///////////////////////////////////////////////////////////////////////////////////////
                    // Case of hits in FI31 and FI32 detectors only

                    if ( FI31_corr_hit.size()>0 && FI32_corr_hit.size()>0 && FI33_corr_hit.size()==0 ){

                        //if ( (FI33_corr_hit[0].GetX() <= 5) && (FI31_corr_hit[0].GetX() >= -5) && (FI32_corr_hit[0].GetX() >= -5) && (FI32_corr_hit[0].GetX() <= 5) ){

                            cout << "FIB32 X:" << FI32_corr_hit[0].GetX() <<  "   FIB31 X:" << FI31_corr_hit[0].GetX() << endl;

                            double x0 = FI32_corr_hit[0].GetX();
                            double z0 = 0.0;
                            double x1 = FI31_corr_hit[0].GetX()-26.4+5;
                            double z1 = 117.5;

                            TGraph* gCurrentLine = new TGraph(2);
                            gCurrentLine->SetPoint(0, x0, z0);
                            gCurrentLine->SetPoint(1, x1, z1);
                            //gCurrentLine->SetLineColor(fiberColors[fiberTrajectoryColorIndex % nFiberColors]);

                            fiberXMin = std::min(fiberXMin, std::min(x0, x1));
                            fiberXMax = std::max(fiberXMax, std::max(x0, x1));
                            fiberZMin = std::min(fiberZMin, std::min(z0, z1));
                            fiberZMax = std::max(fiberZMax, std::max(z0, z1));

                            gCurrentLine->SetLineColor(kGreen);
                            gCurrentLine->SetMarkerColor(kBlack);
                            gCurrentLine->SetMarkerStyle(20);
                            gCurrentLine->SetMarkerSize(2.5);
                            gCurrentLine->SetLineWidth(2);

                            cFiberTrajectory->cd();
                            if (fiberXMin < fiberXMax) {
                                double xMargin = 5.0;
                                double zMargin = 10.0;
                                double xMin = fiberXMin - xMargin;
                                double xMax = fiberXMax + xMargin;
                                double xSpan = xMax - xMin;
                                if (xSpan < minFiberXRange) {
                                    double xCenter = 0.5*(xMin + xMax);
                                    xMin = xCenter - 0.5*minFiberXRange;
                                    xMax = xCenter + 0.5*minFiberXRange;
                                }
                                double zMin = fiberZMin - zMargin;
                                double zMax = fiberZMax + zMargin;
                                double zSpan = zMax - zMin;
                                if (zSpan < minFiberZRange) {
                                    double zCenter = 0.5*(zMin + zMax);
                                    zMin = zCenter - 0.5*minFiberZRange;
                                    zMax = zCenter + 0.5*minFiberZRange;
                                }
                                hFiberTrajectoryFrame->GetXaxis()->SetLimits(xMin, xMax);
                                hFiberTrajectoryFrame->GetYaxis()->SetRangeUser(zMin, zMax);
                            }
                            if (!fiberTrajectoryDrawn) {
                                hFiberTrajectoryFrame->Draw();
                                fiberTrajectoryDrawn = true;
                            } else {
                                hFiberTrajectoryFrame->Draw("axis same");
                            }
                            gCurrentLine->Draw("PL same");
                            cFiberTrajectory->Modified();
                            cFiberTrajectory->Update();

                            fiberTrajectoryColorIndex++;
                        //}
                    }

                    good_events++;
                }

                /////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////
                /////////////////////////////////////////////////////////////////////////////////////////

            }      // ----> ToFD condition
        }         // ----> Single hit events condition
    }             // ----> Loop over events

cout << "Good events: " << good_events << endl;
cout << "Collinear events: " << collinear << endl;
cout << "Non-collinear events: " << noncollinear << endl;
cout << "Events with one hit at all fibers: " << evts1 << endl;
cout << "Events with one hit at all fibers and MWPC: " << evts2 << endl;
cout << "Events with one hit at all fibers, MWPC and LOS: " << evts3 << endl;
cout << "Events with one hit at all fibers, MWPC, LOS and ToFD: " << evts4 << endl;
cout << "Events with one hit at all fibers, MWPC, LOS and 4 or 5 hits in ToFD: " << evts5 << endl;


// TCanvas* c_TTT10_X  = new TCanvas();
// c_TTT10_X ->cd(1);
// h_TTT10_X ->Draw("hist");

// TCanvas* c_MWPC_X  = new TCanvas();
// c_MWPC_X ->cd(1);
// h_MWPC_X ->Draw("hist");

// TCanvas* c_TTT10_MWPC_X = new TCanvas();
// c_TTT10_MWPC_X->cd(1);
// h_TTT10_MWPC_X->Draw("colz");

// TCanvas* c_FI31_X  = new TCanvas();
// c_FI31_X->cd(1);
// h_FI31_X->Draw("hist");

// TCanvas* c_FI32_X  = new TCanvas();
// c_FI32_X->cd(1);
// h_FI32_X->Draw("hist");

// TCanvas* c_FI33_X  = new TCanvas();
// c_FI33_X->cd(1);
// h_FI33_X->Draw("hist");


TCanvas* c_hits_FI31  = new TCanvas();
c_hits_FI31->cd(1);
h_hits_FI31->Draw("hist");

TCanvas* c_hits_FI32  = new TCanvas();
c_hits_FI32->cd(1);
h_hits_FI32->Draw("hist");

TCanvas* c_hits_FI33  = new TCanvas();
c_hits_FI33->cd(1);
h_hits_FI33->Draw("hist");

TCanvas* c_hits_TTT10  = new TCanvas();
c_hits_TTT10->cd(1);
h_hits_TTT10->Draw("hist");

TCanvas* c_hits_MWPC  = new TCanvas();
c_hits_MWPC->cd(1);
h_hits_MWPC->Draw("hist");

TCanvas* c_hits_TOFD  = new TCanvas();
c_hits_TOFD->cd(1);
h_hits_TOFD->Draw("hist");

TCanvas* c_hits_LOS  = new TCanvas();
c_hits_LOS->cd(1);
h_hits_LOS->Draw("hist");


// TCanvas* c_FICOMB_X  = new TCanvas();
// c_FICOMB_X->cd(1);
// h_FICOMB_X->Draw("hist");

// TCanvas* c_TOFD_X  = new TCanvas();
// c_TOFD_X->cd(1);
// h_TOFD_X->Draw("hist");

// // ///////////////////////////////////////////////////////////////////

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

// Save all figures_long to disk when the macro finishes
if (gSystem) {
    gSystem->mkdir("figures_long", kTRUE);
}

// Create ROOT file for interactive browsing
TFile* rootFile = new TFile("figures_long/Fib31_Fib33_analysis.root", "RECREATE");

// Write canvases
// cFiberTrajectory->Write();
c_hits_FI31->Write();
c_hits_FI32->Write();
c_hits_FI33->Write();
c_hits_TTT10->Write();
c_hits_MWPC->Write();
c_hits_TOFD->Write();
c_hits_LOS->Write();
c1->Write();
c3->Write();
c_TTT10_MWPC_X_area->Write();
c_TTT10_FI33_X_corr->Write();
c_TTT10_FI32_X_corr->Write();
c_MWPC_FI33_X_corr->Write();
c_MWPC_FI32_X_corr->Write();
c_TTT10_FI33_X_single_corr->Write();
c_TTT10_FI32_X_single_corr->Write();
c_MWPC_FI33_X_single_corr->Write();
c_MWPC_FI32_X_single_corr->Write();
c_FI33_FI32_X->Write();
c_FI31_FI32_X->Write();
c_FIB31_FIB33_combined_FIB32_area->Write();
c_MWPC_FI31_FI33_X_single_corr->Write();

// Write histograms
h_TTT10_X->Write();
h_MWPC_X->Write();
h_TTT10_MWPC_X->Write();
h_FI31_X->Write();
h_FI32_X->Write();
h_FI33_X->Write();
h_FICOMB_X->Write();
h_TOFD_X->Write();
h_hits_TTT10->Write();
h_hits_MWPC->Write();
h_hits_FI31->Write();
h_hits_FI32->Write();
h_hits_FI33->Write();
h_hits_TOFD->Write();
h_hits_LOS->Write();
h_TTT10_MWPC_X_corr->Write();
h_FI33_FI32_X->Write();
h_FI31_FI32_X->Write();
h_FI33_FI31_combined_FI32_X->Write();
h_TTT10_MWPC_X_area->Write();
h_FIB31_FIB33_combined_FIB32_area->Write();
h_TTT10_FI32_X_corr->Write();
h_TTT10_FI33_X_corr->Write();
h_MWPC_FI32_X_corr->Write();
h_MWPC_FI33_X_corr->Write();
h_MWPC_FI31_FI33_X_single_corr->Write();
h_TTT10_FI32_X_single_corr->Write();
h_TTT10_FI33_X_single_corr->Write();
h_MWPC_FI32_X_single_corr->Write();
h_MWPC_FI33_X_single_corr->Write();

rootFile->Close();

// Save PNG figures_long
// cFiberTrajectory->SaveAs("figures_long/FiberTrajectory.png");
c_hits_FI31->SaveAs("figures_long/hits_FI31.png");
c_hits_FI32->SaveAs("figures_long/hits_FI32.png");
c_hits_FI33->SaveAs("figures_long/hits_FI33.png");
c_hits_TTT10->SaveAs("figures_long/hits_TTT10.png");
c_hits_MWPC->SaveAs("figures_long/hits_MWPC.png");
c_hits_TOFD->SaveAs("figures_long/hits_TOFD.png");
c_hits_LOS->SaveAs("figures_long/hits_LOS.png");
c1->SaveAs("figures_long/TTT10_MWPC_correlation.png");
c3->SaveAs("figures_long/FI33_FI31_combined_FI32_correlation.png");
c_TTT10_MWPC_X_area->SaveAs("figures_long/TTT10_MWPC_X_area.png");
c_TTT10_FI33_X_corr->SaveAs("figures_long/TTT10_FI33_X_corr.png");
c_TTT10_FI32_X_corr->SaveAs("figures_long/TTT10_FI32_X_corr.png");
c_MWPC_FI33_X_corr->SaveAs("figures_long/MWPC_FI33_X_corr.png");
c_MWPC_FI32_X_corr->SaveAs("figures_long/MWPC_FI32_X_corr.png");
c_TTT10_FI33_X_single_corr->SaveAs("figures_long/TTT10_FI33_X_single_corr.png");
c_TTT10_FI32_X_single_corr->SaveAs("figures_long/TTT10_FI32_X_single_corr.png");
c_MWPC_FI33_X_single_corr->SaveAs("figures_long/MWPC_FI33_X_single_corr.png");
c_MWPC_FI32_X_single_corr->SaveAs("figures_long/MWPC_FI32_X_single_corr.png");
c_FI33_FI32_X->SaveAs("figures_long/FI33_FI32_X.png");
c_FI31_FI32_X->SaveAs("figures_long/FI31_FI32_X.png");
c_FIB31_FIB33_combined_FIB32_area->SaveAs("figures_long/FIB31_FIB33_combined_FIB32_area.png");
c_MWPC_FI31_FI33_X_single_corr->SaveAs("figures_long/MWPC_FI31_FI33_X_single_corr.png");

}