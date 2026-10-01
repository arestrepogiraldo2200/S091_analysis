/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname1[200],fname2[200],fname3[200],fname4[200],fname5[200],fname6[200];
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

TH2D* h_TTT10_FI32_X_single_corr;
TH2D* h_TTT10_FI33_X_single_corr;
TH2D* h_MWPC_FI32_X_single_corr;
TH2D* h_MWPC_FI33_X_single_corr;

vector<R3BTttxHitData*> TTT10_corr_hit;
vector<R3BMwpcHitData*>  MWPC_corr_hit;
vector<R3BFiberMAPMTHitData*> FIB32_corr_hit;
vector<R3BFiberMAPMTHitData*> FIB33_corr_hit;

vector<R3BTttxHitData*> TTT10_corr_single_hit;
vector<R3BMwpcHitData*>  MWPC_corr_single_hit;
vector<R3BFiberMAPMTHitData*> FIB32_corr_single_hit;
vector<R3BFiberMAPMTHitData*> FIB33_corr_single_hit;

void analysis_detectors_GLAD_xpos_corr_runs(){

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

    h_TOFD_X = new TH1D("TOFD x position","TOFD x position",100,-60,60);
    h_TOFD_X->GetXaxis()->SetTitle("TOFD x position [cm]");

    h_FI32_FI33_X = new TH2D("X position correlations FI32-FI33","X position correlations FI32-FI33",fiber_bins,-26.4,26.4,fiber_bins,-26.4,26.4);
    h_FI32_FI33_X->GetXaxis()->SetTitle("FI32 X-position [cm]");
    h_FI32_FI33_X->GetYaxis()->SetTitle("FI33 X-position [cm]");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_MWPC_X_corr = new TH2D("X position MWPC-TTT10 correlation","X position MWPC-TTT10 correlation",mwpc1_bins,-150,0,32,-50,50);
    h_TTT10_MWPC_X_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_TTT10_MWPC_X_corr->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI32_FI33_X_corr = new TH2D("X position FI32-FI33 correlation","X position FI32-FI33 correlation",fiber_bins,-26.4,26.4,fiber_bins,-26.4,26.4);
    h_FI32_FI33_X_corr->GetXaxis()->SetTitle("FI32 X-position [cm]");
    h_FI32_FI33_X_corr->GetYaxis()->SetTitle("FI33 X-position [cm]");

    /////////////////////////////////////////////////////////////////////////////////////////

    h_TTT10_MWPC_X_area = new TH2D("X position correlations MWPC-TTT10 area","X position correlations MWPC-TTT10 area",mwpc1_bins,-150,0,32,-50,50);
    h_TTT10_MWPC_X_area->GetXaxis()->SetTitle("MWPC X-position" );
    h_TTT10_MWPC_X_area->GetYaxis()->SetTitle("TTT10 X-position");

    h_FI32_FI33_X_area = new TH2D("X position correlations FI32-FI33 area","X position correlations FI32-FI33 area",fiber_bins,-26.4,26.4,fiber_bins,-26.4,26.4);
    h_FI32_FI33_X_area->GetXaxis()->SetTitle("FI32 X-position [cm]");
    h_FI32_FI33_X_area->GetYaxis()->SetTitle("FI33 X-position [cm]");

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

    h_TTT10_FI32_X_single_corr = new TH2D("X position TTT10-FIB32 correlation single hit","X position TTT10-FIB32 correlation single hit",32,-50,50,fiber_bins,-26.4,26.4);
    h_TTT10_FI32_X_single_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI32_X_single_corr->GetYaxis()->SetTitle("FIB32 X-position [cm]");

    h_TTT10_FI33_X_single_corr = new TH2D("X position TTT10-FIB33 correlation single hit","X position TTT10-FIB33 correlation single hit",32,-50,50,fiber_bins,-26.4,26.4);
    h_TTT10_FI33_X_single_corr->GetXaxis()->SetTitle("TTT10 X-position");
    h_TTT10_FI33_X_single_corr->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    h_MWPC_FI32_X_single_corr = new TH2D("X position MWPC-FIB32 correlation single hit","X position MWPC-FIB32 correlation single hit",mwpc1_bins,-150,0,fiber_bins,-26.4,26.4);
    h_MWPC_FI32_X_single_corr->GetXaxis()->SetTitle("MWPC X-position");
    h_MWPC_FI32_X_single_corr->GetYaxis()->SetTitle("FIB32 X-position [cm]");

    h_MWPC_FI33_X_single_corr = new TH2D("X position MWPC-FIB33 correlation single hit","X position MWPC-FIB33 correlation single hit",mwpc1_bins,-150,0,fiber_bins,-26.4,26.4);
    h_MWPC_FI33_X_single_corr->GetXaxis()->SetTitle("MWPC X-position" );
    h_MWPC_FI33_X_single_corr->GetYaxis()->SetTitle("FIB33 X-position [cm]");

    /////////////////////////////////////////////////////////////////////////////////////////

    // Load the root file 

    TChain* chain = new TChain("evt"); 
    sprintf(fname1,"./merged_MDF_runs.root");        
    //sprintf(fname1,"./main0252_0001_unpacked_MDF.root");

    chain->Reset();
    chain->Add(fname1); 

    // /////////////////////////////////////////////////////////////////////////////////////////

    // MERGED FILES 11C (runs 320 and 321), 13C (runs 323 and 324), 14C (run 325), 10C (run 326) and 12C (run 327)

    // sprintf(fname1,"./main0321_0001_unpacked_MDF.root");        // FI33 seems displaced
    // sprintf(fname2,"./main0323_0001_unpacked_MDF.root");        // Good! 
    // sprintf(fname3,"./main0324_0001_unpacked_MDF.root");        // Good! 
    // sprintf(fname4,"./main0325_0001_unpacked_MDF.root");        // Good, need to adjust slopes of linear fit
    // sprintf(fname5,"./main0326_0001_unpacked_MDF.root");        // Good, but some data missing
    // sprintf(fname6,"./main0327_0001_unpacked_MDF.root");        // Good 

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
	
	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in the merged tree: " << nevent << endl;

    int good_events = 0;

	for(Long64_t e=0; e<nevent; e++)
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
        FIB32_corr_hit.clear();
        FIB33_corr_hit.clear();
        TTT10_corr_single_hit.clear();
        MWPC_corr_single_hit.clear();
        FIB32_corr_single_hit.clear();
        FIB33_corr_single_hit.clear();

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

                        if ( (TTT10_X<0.6*MWPC_X+47)  &&  (TTT10_X>0.6*MWPC_X+37) ){      //  && MWPC_X<(TTT10_X-37)/0.6  && MWPC_X<(TTT10_X-47)/0.6 ){ 

                            h_TTT10_MWPC_X_area->Fill(MWPC_X, TTT10_X);
                            TTT10_corr_hit.push_back(hitTtt10[q]);
                            MWPC_corr_hit.push_back(hitMwpc[p]);
                        }
                    }
                } 
            }
        }

        /////////////////////////////////////////////////////////////////////////////////////////

        double hitFI31_X;

        if( Fi31HitPerEvent > 0 ){

            hitFi31 = new R3BFiberMAPMTHitData*[Fi31HitPerEvent];
            
            for(int p = 0; p < Fi31HitPerEvent; p++){

                hitFi31[p] =(R3BFiberMAPMTHitData*)hitFI31->At(p);
                hitFI31_X = hitFi31[p]->GetX();
                h_FI31_X->Fill(hitFI31_X);
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

                        // Area cut condition between regressions 
                        if ( (hitFI33_X<1.08*hitFI32_X-21.5)  && (hitFI33_X>1.08*hitFI32_X-25.5) ){   //  && hitFI32_X<(hitFI33_X+24.5)/1.08  && hitFI32_X>(hitFI33_X+22.5)/1.08 ){

                            h_FI32_FI33_X_area->Fill(hitFI32_X, hitFI33_X);
                            FIB33_corr_hit.push_back(hitFi33[q]);
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
            // cout << "Hits in event " << e << " of MWPC ".q << MWPC_corr_hit.size()  << endl;
            // cout << "Hits in event " << e << " of FI32 " << FIB32_corr_hit.size()  << endl;
            // cout << "Hits in event " << e << " of FI33 " << FIB32_corr_hit.size()  << endl;
            // cout << endl;

            ////////////////////////////////////////////////////////////////////////////////////////////
            // THIS PART FILLS HISTOGRAMS USING THE CORRELATED REGIONS ONLY AND TAKING A SUBSET OF THEM
            ////////////////////////////////////////////////////////////////////////////////////////////

            double TTT10_X_CORR;
            double MWPC_X_CORR ;
            double FIB32_X_CORR;
            double FIB33_X_CORR;

            int smallest = std::min(std::min((int)TTT10_corr_hit.size(), (int)MWPC_corr_hit.size()), std::min((int)FIB32_corr_hit.size(), (int)FIB33_corr_hit.size()));

            // TTT10 vs FIB33 ///////////////////////////////////////////////////////////////////////
            for(int i = 0; i < smallest; i++){

                TTT10_X_CORR = TTT10_corr_hit[i]->GetX();
                FIB33_X_CORR = FIB33_corr_hit[i]->GetX();
                h_TTT10_FI33_X_corr->Fill(TTT10_X_CORR, FIB33_X_CORR);
            }

            // TTT10 vs FIB32 ///////////////////////////////////////////////////////////////////////
            for(int i = 0; i < smallest; i++){

                TTT10_X_CORR = TTT10_corr_hit[i]->GetX();
                FIB32_X_CORR = FIB32_corr_hit[i]->GetX();
                h_TTT10_FI32_X_corr->Fill(TTT10_X_CORR, FIB32_X_CORR);
            }

            // MWPC vs FIB33 ///////////////////////////////////////////////////////////////////////
            for(int i = 0; i < smallest; i++){

                MWPC_X_CORR = MWPC_corr_hit[i]->GetX();
                FIB33_X_CORR = FIB33_corr_hit[i]->GetX();
                h_MWPC_FI33_X_corr->Fill(MWPC_X_CORR, FIB33_X_CORR);
            }

            // MWPC vs FIB32 ///////////////////////////////////////////////////////////////////////
            for(int i = 0; i < smallest; i++){

                MWPC_X_CORR = MWPC_corr_hit[i]->GetX();
                FIB32_X_CORR = FIB32_corr_hit[i]->GetX();
                h_MWPC_FI32_X_corr->Fill(MWPC_X_CORR, FIB32_X_CORR);
            }

            ////////////////////////////////////////////////////////////////////////////////////////////
            // THIS PART FILLS HISTOGRAMS USING THE CORRELATED REGIONS ONLY AND SINGLE HIT FOR IT
            ////////////////////////////////////////////////////////////////////////////////////////////

            double largestE_TTT10 = -1;
            double largestE_MWPC  = -1;
            double largestE_FIB32 = -1;
            double largestE_FIB33 = -1;

            double largestE_index_TTT10;
            double largestE_index_MWPC ;
            double largestE_index_FIB32;
            double largestE_index_FIB33;

            for (int i=0; i < TTT10_corr_hit.size(); i++){

                double energyttt10 = TTT10_corr_hit[i]->GetEnergy();
                if (energyttt10 > largestE_TTT10){
                    largestE_TTT10 = energyttt10;
                    largestE_index_TTT10 = i;
                }
            }
            TTT10_corr_single_hit.push_back(TTT10_corr_hit[largestE_index_TTT10]);


            /////////////////////////////////////////////////////////////////////////////////////////

            // for (int i=0; i < MWPC_corr_hit.size(); i++){

            //     double energy_mwpc = MWPC_corr_hit[i]->GetEnergy();

            //     if (energy_mwpc > largestE_MWPC){

            //         largestE_MWPC = energy_mwpc;
            //         largestE_index_MWPC = i;
            //     }
            // }
            // MWPC_corr_single_hit.push_back(MWPC_corr_hit[largestE_index_MWPC]);

            /////////////////////////////////////////////////////////////////////////////////////////

            for (int i=0; i < FIB32_corr_hit.size(); i++){

                double energyfib32 = FIB32_corr_hit[i]->GetEloss();

                if (energyfib32 > largestE_FIB32){
                    largestE_FIB32 = energyfib32;
                    largestE_index_FIB32 = i;
                }
            }
            FIB32_corr_single_hit.push_back(FIB32_corr_hit[largestE_index_FIB32]);


            for (int i=0; i < FIB33_corr_hit.size(); i++){

                double energyfib33 = FIB33_corr_hit[i]->GetEloss();

                if (energyfib33 > largestE_FIB33){
                    largestE_FIB33 = energyfib33;
                    largestE_index_FIB33 = i;
                }
            }
            FIB33_corr_single_hit.push_back(FIB33_corr_hit[largestE_index_FIB33]);

            //cout << TTT10_corr_single_hit.size() << FIB33_corr_single_hit.size() << FIB32_corr_single_hit.size() << MWPC_corr_hit.size() << endl;

            h_TTT10_FI33_X_single_corr->Fill(TTT10_corr_single_hit[0]->GetX(), FIB33_corr_single_hit[0]->GetX());
            h_TTT10_FI32_X_single_corr->Fill(TTT10_corr_single_hit[0]->GetX(), FIB32_corr_single_hit[0]->GetX());
            h_MWPC_FI33_X_single_corr->Fill(MWPC_corr_hit[0]->GetX(), FIB33_corr_single_hit[0]->GetX());                 // !!!!  JUST TAKING FIRST ENTRY OF MWPC 
            h_MWPC_FI32_X_single_corr->Fill(MWPC_corr_hit[0]->GetX(), FIB32_corr_single_hit[0]->GetX());

            good_events++;

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

double m2 = 1.2;
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

}