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

TH2D* h_TTT10_FI32_X;
TH2D* h_TTT10_FI33_X;
TH2D* h_MWPC_FI32_X;
TH2D* h_MWPC_FI33_X;

TH2D* h_TTT10_FI32_X_corr;
TH2D* h_TTT10_FI33_X_corr;
TH2D* h_MWPC_FI32_X_corr;
TH2D* h_MWPC_FI33_X_corr;

//&&&&&&&&&&&&&&&&&&
TH2D* h_TEST;
//&&&&&&&&&&&&&&&&&&

vector<R3BTttxHitData*> TTT10_corr_hit;
vector<R3BMwpcHitData*>  MWPC_corr_hit;
vector<R3BFiberMAPMTHitData*> FIB32_corr_hit;
vector<R3BFiberMAPMTHitData*> FIB33_corr_hit;

void analysis_detectors_GLAD_xpos_corr(){


    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // CHECK THE BINING ACCORDIG TO THE MAP TO 
    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    
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

    h_TOFD_X = new TH1D("TOFD x position","TOFD x position",100,-60,60);
    h_TOFD_X->GetXaxis()->SetTitle("TOFD x position");

    h_FI32_FI33_X = new TH2D("X position correlations FI32-FI33","X position correlations FI32-FI33",200,-30,30,200,-30,30);
    h_FI32_FI33_X->GetXaxis()->SetTitle("FI32 X-position");
    h_FI32_FI33_X->GetYaxis()->SetTitle("FI33 X-position");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_MWPC_X_corr = new TH2D("X position MWPC-TTT10 correlation","X position MWPC-TTT10 correlation",200,-150,0,32,-50,50);
    h_TTT10_MWPC_X_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_TTT10_MWPC_X_corr->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI32_FI33_X_corr = new TH2D("X position FI32-FI33 correlation","X position FI32-FI33 correlation",200,-30,30,200,-30,30);
    h_FI32_FI33_X_corr->GetXaxis()->SetTitle("FI32 X-position");
    h_FI32_FI33_X_corr->GetYaxis()->SetTitle("FI33 X-position");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_MWPC_X_area = new TH2D("X position correlations MWPC-TTT10 area","X position correlations MWPC-TTT10 area",200,-150,0,32,-50,50);
    h_TTT10_MWPC_X_area->GetXaxis()->SetTitle("MWPC X-position" );
    h_TTT10_MWPC_X_area->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI32_FI33_X_area = new TH2D("X position correlations FI32-FI33 area","X position correlations FI32-FI33 area",200,-30,30,200,-30,30);
    h_FI32_FI33_X_area->GetXaxis()->SetTitle("FI32 X-position");
    h_FI32_FI33_X_area->GetYaxis()->SetTitle("FI33 X-position");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_FI32_X_corr = new TH2D("X position TTT10-FIB32 correlation","X position TTT10-FIB32 correlation",32,-50,50,200,-30,30);
    h_TTT10_FI32_X_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI32_X_corr->GetYaxis()->SetTitle("FIB32 X-position");

    h_TTT10_FI33_X_corr = new TH2D("X position TTT10-FIB33 correlation","X position TTT10-FIB33 correlation",32,-50,50,200,-30,30);
    h_TTT10_FI33_X_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI33_X_corr->GetYaxis()->SetTitle("FIB33 X-position");

    h_MWPC_FI32_X_corr = new TH2D("X position MWPC-FIB32 correlation","X position MWPC-FIB32 correlation",200,-150,0,200,-30,30);
    h_MWPC_FI32_X_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_MWPC_FI32_X_corr->GetYaxis()->SetTitle("FIB32 X-position");

    h_MWPC_FI33_X_corr = new TH2D("X position MWPC-FIB33 correlation","X position MWPC-FIB33 correlation",200,-150,0,200,-30,30);
    h_MWPC_FI33_X_corr->GetXaxis()->SetTitle("MWPC X-position" );
    h_MWPC_FI33_X_corr->GetYaxis()->SetTitle("FIB33 X-position");

    h_TEST = new TH2D("TEST","TEST",200,-30,30,200,-30,30);
    h_TEST->GetXaxis()->SetTitle("MWPC X-position");
    h_TEST->GetYaxis()->SetTitle("TTT10 X-position");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt");  
    //sprintf(fname, "./main0252_0001_unpacked_MDF.root");
    // sprintf(fname, "./main0320_0001_unpacked_MDF.root");        // not good
    // sprintf(fname, "./main0321_0001_unpacked_MDF.root");        // FI33 seems displaced
    // sprintf(fname, "./main0322_0001_unpacked_MDF.root");        // Too few events
    // sprintf(fname, "./main0323_0001_unpacked_MDF.root");        // Good! 
    // sprintf(fname, "./main0324_0001_unpacked_MDF.root");        // Good! 
    // sprintf(fname, "./main0325_0001_unpacked_MDF.root");        // Good, need to adjust slopes of linear fit
    // sprintf(fname, "./main0326_0001_unpacked_MDF.root");        // Good, but some data missing
    sprintf(fname, "./main0327_0001_unpacked_MDF.root");        // Good
    // sprintf(fname, "./main0328_0001_unpacked_MDF.root");        // Good, but some data missing

    chain->Reset();
    chain->Add(fname); 

    // /////////////////////////////////////////////////////////////////////////////////////////
    // // Save the merged tree to a new file
    // TFile* outfile = new TFile("merged_MDF_runs.root", "RECREATE");
    // TTree* newtree = chain->CloneTree(-1, "fast");
    // newtree->Write();
    // outfile->Close();
    // cout << "Merged tree saved to merged_MDF_runs.root" << endl;
    
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

        // Cleaning the vectors 
        TTT10_corr_hit.clear();
        MWPC_corr_hit.clear();
        FIB32_corr_hit.clear();
        FIB33_corr_hit.clear();

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

                if( Ttt10HitPerEvent > 0 ){

                    hitTtt10 = new R3BTttxHitData*[Ttt10HitPerEvent];
                    
                    for(int q = 0; q < Ttt10HitPerEvent; q++){

                        hitTtt10[q] =(R3BTttxHitData*)hitTTT10->At(q);
                        double TTT10_X = hitTtt10[q]->GetX();
                        h_TTT10_MWPC_X->Fill(MWPC_X, TTT10_X);
                        h_TTT10_MWPC_X_corr->Fill(MWPC_X, TTT10_X);
                        int corr1 = 0;

                        if ( (TTT10_X<0.6*MWPC_X+47)  &&  (TTT10_X>0.6*MWPC_X+37) ){//  && MWPC_X<(TTT10_X-37)/0.6  && MWPC_X<(TTT10_X-47)/0.6 ){ 

                            h_TTT10_MWPC_X_area->Fill(MWPC_X, TTT10_X);
                            TTT10_corr_hit.push_back(hitTtt10[q]);
                            corr1 = 1;
                        }

                        if (corr1 == 1){
                            MWPC_corr_hit.push_back(hitMwpc[p]);
                        }
                    }
                } 
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
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////

        if( Fi32HitPerEvent > 0 ){

            hitFi32 = new R3BFiberMAPMTHitData*[Fi32HitPerEvent];
            
            for(int p = 0; p < Fi32HitPerEvent; p++){

                hitFi32[p] =(R3BFiberMAPMTHitData*)hitFI32->At(p);
                double hitFI32_X = hitFi32[p]->GetX();
                h_FI32_X->Fill(hitFI32_X);

                /////////////////////////////////////////////////////////////////////////////////////////

                if( Fi33HitPerEvent > 0 ){

                    hitFi33 = new R3BFiberMAPMTHitData*[Fi33HitPerEvent];
                    
                    for(int q = 0; q < Fi33HitPerEvent; q++){

                        hitFi33[q] =(R3BFiberMAPMTHitData*)hitFI33->At(q);
                        double hitFI33_X = hitFi33[q]->GetX();
                        h_FI32_FI33_X->Fill(hitFI32_X, hitFI33_X);
                        h_FI32_FI33_X_corr->Fill(hitFI32_X,hitFI33_X);
                        int corr2 = 0;

                        // Area cut condition between regressions 
                        if ( (hitFI33_X<1.08*hitFI32_X-22.5)  && (hitFI33_X>1.08*hitFI32_X-24.5) ){//  && hitFI32_X<(hitFI33_X+24.5)/1.08  && hitFI32_X>(hitFI33_X+22.5)/1.08 ){

                            h_FI32_FI33_X_area->Fill(hitFI32_X, hitFI33_X);
                            FIB33_corr_hit.push_back(hitFi33[q]);
                            corr2 = 1;
                        }

                        if (corr2 == 1){
                            FIB32_corr_hit.push_back(hitFi32[p]);
                        }
                    }
                }
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////

        if( TofdHitPerEvent > 0 ){

            hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
            
            for(int p = 0; p < TofdHitPerEvent; p++){

                hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
                double TOFD_X = hitTofd[p]->GetX();
                h_TOFD_X->Fill(TOFD_X);
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////

        if (TTT10_corr_hit.size()>0 && MWPC_corr_hit.size()>0 && FIB32_corr_hit.size()>0 && FIB33_corr_hit.size()>0){

            // cout << "Hits in event " << e << " of TTT10 "<< TTT10_corr_hit.size() << endl;
            // cout << "Hits in event " << e << " of MWPC " << MWPC_corr_hit.size()  << endl;
            // cout << "Hits in event " << e << " of FI32 " << FIB32_corr_hit.size()  << endl;
            // cout << "Hits in event " << e << " of FI33 " << FIB32_corr_hit.size()  << endl;
            // cout << endl;

            double largestE_TTT10 = -1;
            double largestE_MWPC  = -1;
            double largestE_FIB32 = -1;
            double largestE_FIB33 = -1;

            double largestE_index_TTT10;
            double largestE_index_MWPC ;
            double largestE_index_FIB32;
            double largestE_index_FIB33;

            // for (int i=0; i < TTT10_corr_hit.size(); i++){

            //     energyttt10 = TTT10_corr_hit[i]->GetEnergy();
            //     if (energyttt10 > largestE_TTT10){
            //         largestE_TTT10 = energyttt10;
            //         largestE_index_TTT10 = i;
            //     }
            // }

            // for (int i=0; i < MWPC_corr_hit.size(); i++){
            //     double energy = MWPC_corr_hit[i]->GetEnergy();
            //     if (energy > largestE_MWPC){
            //         largestE_MWPC = energy;
            //         largestE_index_MWPC = i;
            //     }
            // }

            //cout << TTT10_corr_hit.size() << "  " << MWPC_corr_hit.size() << "  " << FIB32_corr_hit.size() << "  " << FIB33_corr_hit.size() << endl;

            // double energyfib32;
            // for (int i=0; i < FIB32_corr_hit.size(); i++){

            //     try {
            //         energyfib32 = FIB32_corr_hit[i]->GetEloss();
            //         cout << FIB32_corr_hit[i]->GetEloss() << endl;
            //         cout << FIB32_corr_hit.size() << endl;
            //         int a = 1;
            //         throw(a);
            //     }
            //     catch(int a){
            //         cout << "Error occurred" << endl;
            //         continue;
            //     }


            //     if (energyfib32 > largestE_FIB32){
            //         largestE_FIB32 = energyfib32;
            //         largestE_index_FIB32 = i;
            //     }
            // }

            // for (int i=0; i < FIB33_corr_hit.size(); i++){
            //     double energy = FIB33_corr_hit[i]->GetEloss();
            //     if (energy > largestE_FIB33){
            //         largestE_FIB33 = energy;
            //         largestE_index_FIB33 = i;
            //     }
            // }
            
            double TTT10_X_CORR;
            double MWPC_X_CORR;
            double FIB32_X_CORR;
            double FIB33_X_CORR;

            // TTT10 vs FIB33 ///////////////////////////////////////////////////////////////////////7
            for(int i = 0; i < TTT10_corr_hit.size(); i++){

                TTT10_X_CORR = TTT10_corr_hit[i]->GetX();

                for(int j = 0; j < FIB33_corr_hit.size(); j++){

                    FIB33_X_CORR = FIB33_corr_hit[j]->GetX();
                    h_TTT10_FI33_X_corr->Fill(TTT10_X_CORR, FIB33_X_CORR);
                }
            }

            // TTT10 vs FIB32 ///////////////////////////////////////////////////////////////////////7
            for(int i = 0; i < TTT10_corr_hit.size(); i++){

                TTT10_X_CORR = TTT10_corr_hit[i]->GetX();

                for(int j = 0; j < FIB32_corr_hit.size(); j++){

                    FIB32_X_CORR = FIB32_corr_hit[j]->GetX();
                    h_TTT10_FI32_X_corr->Fill(TTT10_X_CORR, FIB32_X_CORR);
                }
            }

            // MWPC vs FIB33 ///////////////////////////////////////////////////////////////////////7
            for(int i = 0; i < MWPC_corr_hit.size(); i++){

                MWPC_X_CORR = MWPC_corr_hit[i]->GetX();

                for(int j = 0; j < FIB33_corr_hit.size(); j++){

                    FIB33_X_CORR = FIB33_corr_hit[j]->GetX();
                    h_MWPC_FI33_X_corr->Fill(MWPC_X_CORR, FIB33_X_CORR);
                }
            }

            // MWPC vs FIB32 ///////////////////////////////////////////////////////////////////////7
            for(int i = 0; i < MWPC_corr_hit.size(); i++){

                MWPC_X_CORR = MWPC_corr_hit[i]->GetX();

                for(int j = 0; j < FIB32_corr_hit.size(); j++){

                    FIB32_X_CORR = FIB32_corr_hit[j]->GetX();
                    h_MWPC_FI32_X_corr->Fill(MWPC_X_CORR, FIB32_X_CORR);
                }
            }

            // TEST //////////////////////////////////////////////////////////////////////////////
            for(int i = 0; i < FIB32_corr_hit.size(); i++){

                TTT10_X_CORR = FIB32_corr_hit[i]->GetX();

                for(int j = 0; j < FIB33_corr_hit.size(); j++){

                    MWPC_X_CORR = FIB33_corr_hit[j]->GetX();
                    h_TEST->Fill(MWPC_X_CORR, TTT10_X_CORR);
                }
            }



        }

    //     if ((tpat & 8) == 8){
            
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

TCanvas* c_TOFD_X  = new TCanvas();
c_TOFD_X->cd(1);
h_TOFD_X->Draw("hist");

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

///////////////////////////////////////////////////////////////////

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

///////////////////////////////////////////////////////////////////

TCanvas* c_TEST = new TCanvas();
c_TEST->cd(1);
h_TEST->Draw("colz");

}