// // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF
#include "MDF_s118.h"
// // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF

/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname1[200], fname2[200], fname3[200], fname4[200], fname5[200];
int tpat;
// MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF
// TString MDF_PoQ_file = "./10C_MDF_params.txt";
// TString MDF_PoQ_file = "./11C_MDF_params.txt";
// TString MDF_PoQ_file = "./12C_MDF_params.txt";
// TString MDF_PoQ_file = "./13C_MDF_params.txt";
// TString MDF_PoQ_file = "./14C_MDF_params.txt";
TString MDF_PoQ_file = "./MDF_params.txt";
// TString MDF_PoQ_file = "./MDF_PoQ_1456A_20240221.txt";

const TString isotopeTag = "AllC_FREEZE";
const TString outputRootDir = Form("%s_mdf", isotopeTag.Data());
const TString outputPngDir = Form("%s_mdf", isotopeTag.Data());
const TString detectorLabFile = Form("%s_lab.root", isotopeTag.Data());
const TString outputRootFile = Form("%s/%s_mdf.root", outputRootDir.Data(), isotopeTag.Data());

Double_t mdf_data[7]={0};
// MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF

// Before GLAD
TH1D* h_MWPC1_X;

// After GLAD
TH1D* h_FI31_X;
TH1D* h_FI32_X;
TH1D* h_FI33_X;
TH1D* h_FICOMB_X;
TH1D* h_TOFD_X;

// Area histogram
TH1D* h_area;

// Hits per event histograms
TH1D* h_hits_LOS;
TH1D* h_hits_MWPC1;
TH1D* h_hits_FI31;
TH1D* h_hits_FI32;
TH1D* h_hits_FI33;
TH1D* h_hits_TOFD;

// Correlation
TH2D* h_FI33_FI32_X;
TH2D* h_FI31_FI32_X;
TH2D* h_FI33_FI31_combined_FI32_X;

// P/Q histogram
TH1D* h_PoQ;

// Incidence angle histogram
TH1D* h_incidence_angle;

// Correlated area
TH2D* h_FIB31_FIB33_combined_FIB32_area;

// Charge and momentum space
TH2D* h_PoQvsQ;

// Correlations between detectors before and after
TH2D* h_MWPC1_FI31_FI33_X_single_corr;
TH2D* h_MWPC1_FI32_X_single_corr;

vector<R3BMwpcHitData*> MWPC1_hit;
vector<R3BFiberMAPMTHitData> FI32_corr_hit;
vector<R3BFiberMAPMTHitData> FI33_corr_hit;
vector<R3BFiberMAPMTHitData> FI31_corr_hit;

// Calculate the area of the triangle between hits in the fiber detectors
double TriangleArea(double x1, double y1, double x2, double y2,  double x3, double y3){
    return 0.5*std::abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)));
}

// Function to correct the spread of the beam
double FIB3133_CORRECTION(double x){
    double a = std::abs(x);
    double y = -2.004 - 0.150 * a - 0.00418 * a * a;
    return x < 0 ? -1*y-1 : y;
}

double averageMWPC1(std::string filename){

    TChain chain("evt");
    
    chain.Add(filename.c_str());

    TClonesArray* hitMWPC1 = new TClonesArray("R3BMwpcHitData",5);
    chain.SetBranchAddress("Mwpc1HitData", &hitMWPC1);

    TClonesArray* hitLOS = new TClonesArray("R3BLosHitData",5);
    chain.SetBranchAddress("LosHit", &hitLOS);

    TClonesArray* hitFI31 = new TClonesArray("R3BFiberMAPMTHitData",5);
    chain.SetBranchAddress("Fi31Hit", &hitFI31);

    TClonesArray* hitFI32 = new TClonesArray("R3BFiberMAPMTHitData",5);
    chain.SetBranchAddress("Fi32Hit", &hitFI32);

    TClonesArray* hitFI33 = new TClonesArray("R3BFiberMAPMTHitData",5);
    chain.SetBranchAddress("Fi33Hit", &hitFI33);

    TClonesArray* hitTOFD = new TClonesArray("R3BTofdHitData",5);
    chain.SetBranchAddress("TofdHit", &hitTOFD);

    const Long64_t nevent = chain.GetEntries();
    double sumMWPC1 = 0.0;
    Long64_t nMWPC1Entries = 0;

    for (Long64_t e = 0; e < nevent; ++e) {
        hitMWPC1->Clear();
        hitLOS->Clear();
        hitFI31->Clear();
        hitFI32->Clear();
        hitFI33->Clear();
        hitTOFD->Clear();
        chain.GetEntry(e);

        const Long64_t LosHitPerEvent = hitLOS->GetEntriesFast();
        const Long64_t mwpc1Hits = hitMWPC1->GetEntriesFast();
        const Long64_t Fi31HitPerEvent = hitFI31->GetEntriesFast();
        const Long64_t Fi32HitPerEvent = hitFI32->GetEntriesFast();
        const Long64_t Fi33HitPerEvent = hitFI33->GetEntriesFast();
        const Long64_t TofdHitPerEvent = hitTOFD->GetEntriesFast();

        bool b_tofdq = false;
        for (Long64_t p = 0; p < TofdHitPerEvent; ++p) {
            auto* hitTofd = static_cast<R3BTofdHitData*>(hitTOFD->At(p));
            if (!hitTofd) {
                continue;
            }

            const double eloss = hitTofd->GetEloss();
            if (eloss > 5.5 && eloss < 6.5) {
                b_tofdq = true;
                break;
            }
        }

        if (!((TofdHitPerEvent == 4 || TofdHitPerEvent == 5) &&
              LosHitPerEvent == 1 &&
              mwpc1Hits == 1 &&
              Fi32HitPerEvent == 1 &&
              (Fi31HitPerEvent == 1 || Fi33HitPerEvent == 1) &&
              b_tofdq)) {
            continue;
        }

        for (Long64_t p = 0; p < mwpc1Hits; ++p) {
            auto* hit = static_cast<R3BMwpcHitData*>(hitMWPC1->At(p));
            if (!hit) {
                continue;
            }

            const double hitX = hit->GetX();
            if (!std::isfinite(hitX)) {
                continue;
            }

            sumMWPC1 += hitX;
            ++nMWPC1Entries;
        }
    }

    if (nMWPC1Entries == 0) {
        std::cout << "No valid MWPC1 hits found in file: " << filename << std::endl;
        return 0.0;
    }

    const double average = sumMWPC1 / static_cast<double>(nMWPC1Entries);
    std::cout << "Average MWPC1 X over " << nMWPC1Entries << " hits in " << filename << ": " << average << std::endl;

    return average;
}

void MDF_s118_FREEZE(){

    // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF
    MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_file);
    // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF

    /////////////////////////////////////////////////////////////////////////////////////////
    // Set histograms canvas and axis titles here

    int fiber_bins = 150; 
    int mwpc1_bins = 200; 

    // X-position histograms
    h_MWPC1_X = new TH1D("MWPC1 x position","MWPC1 x position",mwpc1_bins,-150,0);
    h_MWPC1_X->GetXaxis()->SetTitle("MWPC1 x position");

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

    // Hits histograms
    h_hits_LOS = new TH1D("Hits per event LOS", "Hits per event LOS",20, 0, 20);
    h_hits_LOS->GetXaxis()->SetTitle("Hits per event LOS");

    h_hits_MWPC1 = new TH1D("Hits per event MWPC1", "Hits per event MWPC1",20, 0, 20);
    h_hits_MWPC1->GetXaxis()->SetTitle("Hits per event MWPC1");

    h_hits_FI31 = new TH1D("Hits per event 31", "Hits per event 31",20, 0, 20);
    h_hits_FI31->GetXaxis()->SetTitle("Hits per event 31");

    h_hits_FI32 = new TH1D("Hits per event 32", "Hits per event 32",20, 0, 20);
    h_hits_FI32->GetXaxis()->SetTitle("Hits per event 32");

    h_hits_FI33 = new TH1D("Hits per event 33", "Hits per event 33",20, 0, 20);
    h_hits_FI33->GetXaxis()->SetTitle("Hits per event 33");

    h_hits_TOFD = new TH1D("Hits per event ToFD", "Hits per event ToFD",20, 0, 20);
    h_hits_TOFD->GetXaxis()->SetTitle("Hits per event ToFD");

    // Area histograms
    h_area = new TH1D("Area of hits in fibers","Area of hits in fibers", 100, 0, 500);
    h_area->GetXaxis()->SetTitle("Area of hits in fibers");

    /////////////////////////////////////////////////////////////////////////////////////////

    // Correlations in X-positions
    h_MWPC1_FI32_X_single_corr = new TH2D("X position MWPC1-FIB32 correlation single hit","X position MWPC1-FIB32 correlation single hit", mwpc1_bins,-150,0,fiber_bins,-26.4,26.4);
    h_MWPC1_FI32_X_single_corr->GetXaxis()->SetTitle("MWPC1 X-position");
    h_MWPC1_FI32_X_single_corr->GetYaxis()->SetTitle("FIB32 X-position [cm]");

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

    h_MWPC1_FI31_FI33_X_single_corr = new TH2D("X position MWPC1 and FIB31-FIB33 correlation single hit","X position MWPC1 and FIB31-FIB33 correlation single hit", mwpc1_bins, -150, 0, 300, -52.8, 52.8);
    h_MWPC1_FI31_FI33_X_single_corr->GetXaxis()->SetTitle("MWPC1 X-position [mm]" );
    h_MWPC1_FI31_FI33_X_single_corr->GetYaxis()->SetTitle("FIB31&FIB33 X-position [cm]");

    // Momentum P/Q
    h_PoQ = new TH1D("PoQ from MDF","PoQ from MDF", 5000, 0, 4);
    h_PoQ->GetXaxis()->SetTitle("PoQ [GeV/c]");

    h_PoQvsQ = new TH2D("P/Q vs Q from MDF","P/Q vs Q from MDF", 5000, 0, 4, 3000, 0, 8);
    h_PoQvsQ->GetXaxis()->SetTitle("PoQ [GeV/c]");
    h_PoQvsQ->GetYaxis()->SetTitle("Q");

    // Incidence angle 
    h_incidence_angle = new TH1D("Incidence angle between plane of FIB32 and FIB31-fib33","Incidence angle between plane of FIB32 and FIB31-fib33", 1000, -0.05, 0.05);
    h_incidence_angle->GetXaxis()->SetTitle("θ [rad]");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt"); 
    //sprintf(fname1,"./merged_MDF_runs.root"); 
    sprintf(fname1,"./merged_MDF_runs_v3.root");

    // sprintf(fname1,"./main0326_0001_unpacked_MDF_v3.root");        // 10C
    // sprintf(fname1,"./main0321_0001_unpacked_MDF_v3.root");        // 11C
    // sprintf(fname1,"./main0327_0001_unpacked_MDF_v3.root");        // 12C
    // sprintf(fname1,"./main0324_0001_unpacked_MDF_v3.root");        // 13C
    // sprintf(fname1,"./main0325_0001_unpacked_MDF_v3.root");        // 14C

    // sprintf(fname1,"./merged_partial.root"); 
    // sprintf(fname1,"./main0252_0001_unpacked_MDF_2.root");

    chain->Reset();
    chain->Add(fname1); 
    // chain->Add(fname2); 
    // chain->Add(fname3); 
    // chain->Add(fname4); 
    // chain->Add(fname5); 

    double MWPC_avrg = averageMWPC1(fname1);

    // Load the branches with data as TClonesArrays 
	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("EventHeader.");
    branchData->SetAddress(&header);  
    
    TClonesArray*  hitLOS = new TClonesArray("R3BLosHitData",5); 
	R3BLosHitData** hitLos;                                      
	TBranch *branchHitLOS = chain->GetBranch("LosHit");           
	branchHitLOS->SetAddress(&hitLOS);                            
	Long64_t LosHitPerEvent = 0;  


    TClonesArray* hitMWPC1 = new TClonesArray("R3BMwpcHitData",5);
    R3BMwpcHitData** hitMwpc1;
    TBranch *branchHitMWPC1 = chain->GetBranch("Mwpc1HitData");
    branchHitMWPC1->SetAddress(&hitMWPC1);
    Long64_t Mwpc1HitPerEvent = 0;  	


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
    // Struct to store the hits in the lab frame
    struct DetectorPosition {
        Double_t x;
        Double_t z;
    };

    DetectorPosition mwpc1{};
    DetectorPosition fib32{};
    DetectorPosition fib3133{};

    TFile* detectorFile = new TFile(detectorLabFile, "RECREATE");
    detectorFile->cd();

    TTree* detectorTree = new TTree("detectors", "Detector x and z positions in lab frame");
    detectorTree->Branch("MWPC1", &mwpc1, "x/D:z/D");
    detectorTree->Branch("Fib32", &fib32, "x/D:z/D");
    detectorTree->Branch("Fib31_Fib33", &fib3133, "x/D:z/D");

    /////////////////////////////////////////////////////////////////////////////////////////

	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in the tree: " << nevent << endl;

    int good_events = 0;
    int collinear = 0, noncollinear = 0;
    int evts1 = 0;
    int evts2 = 0;
    int evts3 = 0;
    int evts4 = 0;
    int evts5 = 0;

    // Constants for the transformation between reference frames of FIB32 and FIB31&33
    double FIB33frame_to_FIB32frame = +26.4;
    double FIB31frame_to_FIB32frame = -26.4;

	for(Long64_t e=0; e < nevent; e++){

        // Cleaning the vectors 
        hitLOS->Clear();
        hitMWPC1->Clear();
        hitFI31->Clear();
        hitFI32->Clear();
        hitFI33->Clear();
        hitTOFD->Clear();
        hitCALIFA->Clear();

        chain->GetEvent(e);

        MWPC1_hit.clear();
        FI32_corr_hit.clear();
        FI33_corr_hit.clear();
        FI31_corr_hit.clear();

        // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF
        mdf_tracks.clear();
        // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF

        tpat             = header->GetTpat();
        LosHitPerEvent   = hitLOS->GetEntriesFast();
        Mwpc1HitPerEvent = hitMWPC1->GetEntriesFast();
        Fi31HitPerEvent  = hitFI31->GetEntriesFast();
        Fi32HitPerEvent  = hitFI32->GetEntriesFast();
        Fi33HitPerEvent  = hitFI33->GetEntriesFast();
        TofdHitPerEvent  = hitTOFD->GetEntriesFast();
        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();

        h_hits_LOS->Fill(LosHitPerEvent);
        h_hits_MWPC1->Fill(Mwpc1HitPerEvent);
        h_hits_FI31->Fill(Fi31HitPerEvent);
        h_hits_FI32->Fill(Fi32HitPerEvent);
        h_hits_FI33->Fill(Fi33HitPerEvent);
        h_hits_TOFD->Fill(TofdHitPerEvent);

        if (Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts1++; }
        if (Mwpc1HitPerEvent == 1 && Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts2++; }
        if (LosHitPerEvent  == 1 && Mwpc1HitPerEvent == 1 && Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts3++; }
        if (TofdHitPerEvent == 1 && LosHitPerEvent  == 1 && Mwpc1HitPerEvent == 1 && Fi32HitPerEvent == 1 && Fi31HitPerEvent == 1 && Fi33HitPerEvent == 1 ){ evts4++; }
        if ((TofdHitPerEvent == 4 || TofdHitPerEvent == 5) && LosHitPerEvent  == 1 && Mwpc1HitPerEvent == 1 && Fi32HitPerEvent == 1 && (Fi31HitPerEvent == 1 || Fi33HitPerEvent == 1) ){ evts5++; }

        // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF
        SetMDFOffsets();
        // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF

        /////////////////////////////////////////////////////////////////////////////////////////
        if ((TofdHitPerEvent == 4 || TofdHitPerEvent == 5) && LosHitPerEvent  == 1 && Mwpc1HitPerEvent == 1 && Fi32HitPerEvent == 1 && (Fi31HitPerEvent == 1 || Fi33HitPerEvent == 1) ){

            /////////////////////////////////////////////////////////////////////////////////////////
            // TOFD charge analysis section
            bool b_tofdq = false;
            double averageEloss = 0.;

            hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
            double totalEloss = 0.;

            for(int p = 0; p < TofdHitPerEvent; p++){

                hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
                double Eloss = hitTofd[p]->GetEloss();
                totalEloss += Eloss;
                double TOFD_X = hitTofd[p]->GetX();
                h_TOFD_X->Fill(TOFD_X);
                
                if(Eloss>5.5 && Eloss<6.5) {                   
                    b_tofdq=true;                          // If at least one of the hits is charge 6 
                }
            }
            
            averageEloss = totalEloss / TofdHitPerEvent;
            /////////////////////////////////////////////////////////////////////////////////////////

            if ( b_tofdq == true ){

                good_events++;

                /////////////////////////////////////////////////////////////////////////////////////////
                // MWPC hit x-position analysis

                double MWPC1_X;
                
                hitMwpc1 = new R3BMwpcHitData*[Mwpc1HitPerEvent];
                
                for(int p = 0; p < Mwpc1HitPerEvent; p++){

                    hitMwpc1[p] =(R3BMwpcHitData*)hitMWPC1->At(p);
                    MWPC1_X = hitMwpc1[p]->GetX();
                    h_MWPC1_X->Fill(MWPC1_X);                                               // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    MWPC1_hit.push_back(hitMwpc1[p]); 
                }
        
                /////////////////////////////////////////////////////////////////////////////////////////
                // Fibers analysis and correlations

                hitFi32 = new R3BFiberMAPMTHitData*[Fi32HitPerEvent];
                double hitFI32_X;

                for(int p = 0; p < Fi32HitPerEvent; p++){

                    hitFi32[p] =(R3BFiberMAPMTHitData*)hitFI32->At(p);
                    auto* hit32 = static_cast<R3BFiberMAPMTHitData*>(hitFI32->At(p));
                    hitFI32_X = hitFi32[p]->GetX();
                    h_FI32_X->Fill(hitFI32_X);                                                            // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Correlation with FIB31
                    if( Fi31HitPerEvent == 1 ){

                        hitFi31 = new R3BFiberMAPMTHitData*[Fi31HitPerEvent];
                        double hitFI31_X;

                        for(int q = 0; q < Fi31HitPerEvent; q++){

                            hitFi31[q] =(R3BFiberMAPMTHitData*)hitFI31->At(q);
                            auto* hit31 = static_cast<R3BFiberMAPMTHitData*>(hitFI31->At(q));

                            hitFI31_X = hitFi31[q]->GetX() + FIB31frame_to_FIB32frame + FIB3133_CORRECTION(hitFI32_X); 

                            h_FI31_X->Fill(hitFI31_X);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            h_FI33_FI31_combined_FI32_X->Fill(hitFI32_X, hitFI31_X);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI31_FI32_X->Fill(hitFI32_X, hitFI31_X);                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FICOMB_X->Fill(hitFI31_X);                                                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            // Correlation condition
                            if ( (hitFI31_X>1.25*hitFI32_X-4) && (hitFI31_X<1.25*hitFI32_X+6) ){  

                                h_FIB31_FIB33_combined_FIB32_area->Fill( hitFI32_X , hitFI31_X);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                FI31_corr_hit.push_back(*hit31);
                                FI32_corr_hit.push_back(*hit32);
                            }
                        }
                    }

                    /////////////////////////////////////////////////////////////////////////////////////////
                    // Correlation with FIB33
                    if( Fi33HitPerEvent == 1 ){

                        hitFi33 = new R3BFiberMAPMTHitData*[Fi33HitPerEvent];
                        double hitFI33_X;

                        for(int q = 0; q < Fi33HitPerEvent; q++){

                            hitFi33[q] =(R3BFiberMAPMTHitData*)hitFI33->At(q);
                            auto* hit33 = static_cast<R3BFiberMAPMTHitData*>(hitFI33->At(q));

                            hitFI33_X = hitFi33[q]->GetX() + FIB33frame_to_FIB32frame + FIB3133_CORRECTION(hitFI32_X);  

                            h_FI33_X->Fill(hitFI33_X);                                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            h_FI33_FI31_combined_FI32_X->Fill(hitFI32_X, hitFI33_X);                      // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FI33_FI32_X->Fill(hitFI32_X, hitFI33_X);                                    // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                            h_FICOMB_X->Fill(hitFI33_X);                                                  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                            // Correlation condition
                            if ( (hitFI33_X>1*hitFI32_X-4) && (hitFI33_X<1*hitFI32_X+4) ){  

                                h_FIB31_FIB33_combined_FIB32_area->Fill( hitFI32_X , hitFI33_X);          // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                                FI33_corr_hit.push_back(*hit33);
                                FI32_corr_hit.push_back(*hit32);
                            }
                        }
                    }  
                }

                // Conditional of correlated hits
                if ( MWPC1_hit.size()>0 && (FI33_corr_hit.size()>0 || FI31_corr_hit.size()>0) && FI32_corr_hit.size()>0 ){

                    if (FI33_corr_hit.size()>0 && FI31_corr_hit.size()== 0){
                        h_MWPC1_FI31_FI33_X_single_corr->Fill( MWPC1_hit[0]->GetX() , FI33_corr_hit[0].GetX() + FIB33frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[0].GetX()) );
                    }
                    if (FI31_corr_hit.size()>0 && FI33_corr_hit.size()== 0){
                        h_MWPC1_FI31_FI33_X_single_corr->Fill( MWPC1_hit[0]->GetX() , FI31_corr_hit[0].GetX() + FIB31frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[0].GetX()));
                    }

                    h_MWPC1_FI32_X_single_corr->Fill( MWPC1_hit[0]->GetX() , FI32_corr_hit[0].GetX());


                    if (FI33_corr_hit.size()>0 && FI31_corr_hit.size()>0 && FI32_corr_hit.size()>0){

                        // Calculation of area
                        double x0 = FI32_corr_hit[0].GetX();
                        double z0 = 0.0;
                        double x1 = FI33_corr_hit[0].GetX()+FIB33frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[0].GetX());
                        double z1 = 100.0;
                        double x2 = FI31_corr_hit[0].GetX()+FIB31frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[0].GetX());
                        double z2 = 117.5;
                        double triangle_area = TriangleArea(x0, z0, x1, z1, x2, z2);

                        h_area->Fill( triangle_area );
                        
                        if (triangle_area < 150){
                            h_MWPC1_FI31_FI33_X_single_corr->Fill( MWPC1_hit[0]->GetX() , FI31_corr_hit[0].GetX()-26.4+5);
                        }
                    }
                }
            }        // ----> ToFD condition

            // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF

            // Setting the detector hit positions 
            // if( MWPC1_hit.size()>0 && MWPC2_corr_hit.size()>0 && FI32_corr_hit.size()>0 && (FI33_corr_hit.size()>0 || FI31_corr_hit.size()>0) ){
            if( MWPC1_hit.size()>0 && FI32_corr_hit.size()>0 && (FI33_corr_hit.size()>0 || FI31_corr_hit.size()>0) ){

                // This loop appends a single hit in FI32 and FI31 and then considers FI32 and FI33 separately
                
                for (int mult = 1; mult <= FI32_corr_hit.size(); mult++){

                    // Corrections for parallel transport of the trajectory so that MWPC1 peak is at zero

                    const Double_t Angle = (72 ) * TMath::Pi() / 180.;
                    const Double_t compl_Angle = (18 ) * TMath::Pi() / 180.;

                    //                                                       MWPC_avrg < 0 so this term is positive as it should be
                    Candidate.mwpc1_dhit.SetXYZ( MWPC1_hit[0]->GetX()/10.0 - MWPC_avrg/10.0  , 0., 0.);          // MWPC1_hit[0]->GetX())/10 +6.8, (MWPC1_hit[0]->GetY())/10.
                    // Candidate.mwpc1_dhit.SetXYZ(-0.44, -0.6, 0.0);                                            // MWPC1_hit[0]->GetX())/10 +6.8, (MWPC1_hit[0]->GetY())/10.
    

                    if (FI33_corr_hit.size()>0 && FI31_corr_hit.size()==0){
                        
                        double X_3133_cframe = FI33_corr_hit[0].GetX()+FIB33frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[mult-1].GetX());
                        double X_32_cframe = (FI32_corr_hit[mult-1].GetX());
                        double FIB32_FIB33_distance = 100;
                        double incidence_angle = TMath::ATan( (X_32_cframe - X_3133_cframe)/(FIB32_FIB33_distance));                        
                        // cout << "Incidence angle: " << incidence_angle << endl;

                        h_incidence_angle->Fill(incidence_angle);

                        double parallel_translation  = MWPC_avrg/10.0*(TMath::Cos(compl_Angle) + TMath::Sin(compl_Angle)*TMath::Tan(incidence_angle)); 

                        Candidate.fi31_33_dhit.SetXYZ( FI33_corr_hit[0].GetX()+FIB33frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[mult-1].GetX()) - parallel_translation, 0., 0.);            // -1*FI33_corr_hit[0].GetX()+26.4-5  ??  FI33_corr_hit[0].GetY()        
                        Candidate.fi32_dhit.SetXYZ((FI32_corr_hit[mult-1].GetX()) - parallel_translation, 0., 0.);                // -1*FI32_corr_hit[mult-1].GetX() ?? FI32_corr_hit[mult-1].GetY()
                    }

                    if (FI31_corr_hit.size()>0 && FI33_corr_hit.size()==0){

                        // Convention of projecting on FIB33
                        double X_3133_cframe = FI31_corr_hit[0].GetX() + FIB31frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[mult-1].GetX());
                        double X_32_cframe = (FI32_corr_hit[mult-1].GetX());
                        double FIB32_FIB31_distance = 100+17.5;
                        double incidence_angle = TMath::ATan( (X_32_cframe - X_3133_cframe)/(FIB32_FIB31_distance));
                        // cout << "Incidence angle: " << incidence_angle << endl;

                        h_incidence_angle->Fill(incidence_angle);

                        double correction_projection = 17.5*TMath::Tan(incidence_angle);
                        double parallel_translation  = MWPC_avrg/10.0*(TMath::Cos(compl_Angle) + TMath::Sin(compl_Angle)*TMath::Tan(incidence_angle)); 

                        Candidate.fi31_33_dhit.SetXYZ((FI31_corr_hit[0].GetX() + FIB31frame_to_FIB32frame + FIB3133_CORRECTION(FI32_corr_hit[mult-1].GetX())) - parallel_translation - correction_projection -10, 0., 0.);            // -1*FI31_corr_hit[0].GetX()-26.4+5  ?? FI31_corr_hit[0].GetY()
                        Candidate.fi32_dhit.SetXYZ((FI32_corr_hit[mult-1].GetX()) - parallel_translation, 0., 0.);                // -1*FI32_corr_hit[mult-1].GetX() ?? FI32_corr_hit[mult-1].GetY()
                    }

                    Candidate.tofd_dhit.SetXYZ(0. , 0. , 0.);    // Ignore tofd for now
                    mdf_tracks.push_back(Candidate);
                }

                if(mdf_tracks.size()>0){ 
                    for(auto & t : mdf_tracks){    
                            
                        // t.ToF      = t.tofd_t;                                     
                        t.mwpc1_lab = t.mwpc1_dhit;
                        t.fi32_lab = t.fi32_dhit;
                        t.fi31_33_lab = t.fi31_33_dhit;

                        // This will transform "det_point" vectors into lab frame
                        
                        TransformPoint(t.mwpc1_lab,    mwpc1_angles + mwpc1_ang_offset ,
                                                    mwpc1_position + mwpc1_pos_offset);

                        TransformPoint(t.fi32_lab,     fi32_angles + fi32_ang_offset ,
                                                    fi32_position + fi32_pos_offset );

                        TransformPoint(t.fi31_33_lab,  fi31_33_angles + fi31_33_ang_offset ,
                                                    fi31_33_position + fi31_33_pos_offset );
                                                                
                        t.PoQ_mdf = CalculatePoQ(t.mwpc1_lab, t.fi31_33_lab, t.fi32_lab, Actual_Current, Reference_Current, MDF_PoQ); 

                        mwpc1.x = t.mwpc1_lab.X();
                        mwpc1.z = t.mwpc1_lab.Z();

                        fib32.x = t.fi32_lab.X();                      // -1*t.fi32_lab.X()
                        fib32.z = t.fi32_lab.Z();

                        fib3133.x = t.fi31_33_lab.X();                 //  -1*t.fi31_33_lab.X()
                        fib3133.z = t.fi31_33_lab.Z();

                        detectorTree->Fill();

                        // cout << t.PoQ_mdf << endl;

                        h_PoQ->Fill(t.PoQ_mdf);
                        h_PoQvsQ->Fill(t.PoQ_mdf, averageEloss);
                    }    
                }               
            }

            // MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF-MDF

        }          // ----> Hits events condition
    }              // ----> Loop over events

cout << "Good events: " << good_events << endl;
cout << "Collinear events: " << collinear << endl;
cout << "Non-collinear events: " << noncollinear << endl;
cout << "Events with one hit at all fibers: " << evts1 << endl;
cout << "Events with one hit at all fibers and MWPC1: " << evts2 << endl;
cout << "Events with one hit at all fibers, MWPC1 and LOS: " << evts3 << endl;
cout << "Events with one hit at all fibers, MWPC1, LOS and ToFD: " << evts4 << endl;
cout << "Events with one hit at all fibers, MWPC1, LOS and 4 or 5 hits in ToFD: " << evts5 << endl;
cout << "Average of MWPC X-positions:" << MWPC_avrg << endl;

TCanvas* c_hits_FI31  = new TCanvas();
c_hits_FI31->cd(1);
h_hits_FI31->Draw("hist");

TCanvas* c_hits_FI32  = new TCanvas();
c_hits_FI32->cd(1);
h_hits_FI32->Draw("hist");

TCanvas* c_hits_FI33  = new TCanvas();
c_hits_FI33->cd(1);
h_hits_FI33->Draw("hist");

TCanvas* c_hits_MWPC1  = new TCanvas();
c_hits_MWPC1->cd(1);
h_hits_MWPC1->Draw("hist");

TCanvas* c_hits_TOFD  = new TCanvas();
c_hits_TOFD->cd(1);
h_hits_TOFD->Draw("hist");

TCanvas* c_hits_LOS  = new TCanvas();
c_hits_LOS->cd(1);
h_hits_LOS->Draw("hist");

TCanvas* c_PoQ  = new TCanvas();
c_PoQ->cd(1);
h_PoQ->Draw("hist");

TCanvas* c_MWPC1_X  = new TCanvas();
c_MWPC1_X->cd(1);
h_MWPC1_X->Draw("hist");

TCanvas* c_FI32_X  = new TCanvas();
c_FI32_X->cd(1);
h_FI32_X->Draw("hist");

TCanvas* c_FI31_X  = new TCanvas();
c_FI31_X->cd(1);
h_FI31_X->Draw("hist");

TCanvas* c_FI33_X  = new TCanvas();
c_FI33_X->cd(1);
h_FI33_X->Draw("hist");

TCanvas* c_FICOMB_X  = new TCanvas();
c_FICOMB_X->cd(1);
h_FICOMB_X->Draw("hist");

TCanvas* c_incidence_angle  = new TCanvas();
c_incidence_angle->cd(1);
h_incidence_angle->Draw("hist");

///////////////////////////////////////////////////////////////////

double m3 = 1.;
double b3 = 0;
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

TCanvas* c_MWPC1_FI32_X_single_corr = new TCanvas();
c_MWPC1_FI32_X_single_corr->cd(1);
h_MWPC1_FI32_X_single_corr->Draw("colz");

TCanvas* c_FI33_FI32_X = new TCanvas();
c_FI33_FI32_X->cd(1);
h_FI33_FI32_X->Draw("colz");

TCanvas* c_FI31_FI32_X = new TCanvas();
c_FI31_FI32_X->cd(1);
h_FI31_FI32_X->Draw("colz");

TCanvas* c_FIB31_FIB33_combined_FIB32_area = new TCanvas();
c_FIB31_FIB33_combined_FIB32_area->cd(1);
h_FIB31_FIB33_combined_FIB32_area->Draw("colz");

TCanvas* c_MWPC1_FI31_FI33_X_single_corr = new TCanvas();
c_MWPC1_FI31_FI33_X_single_corr->cd(1);
h_MWPC1_FI31_FI33_X_single_corr->Draw("colz");

TCanvas* c_area = new TCanvas();
c_area->cd(1);
h_area->Draw("hist");

// Save all isotope-specific MDF output to disk when the macro finishes
if (gSystem) {
    gSystem->mkdir(outputRootDir, kTRUE);
    gSystem->mkdir(outputPngDir, kTRUE);
}

TCanvas* C_PoQvsQ = new TCanvas();
C_PoQvsQ->cd(1);
h_PoQvsQ->Draw("colz");

// Create ROOT file for interactive browsing
TFile* rootFile = new TFile(outputRootFile, "RECREATE");

// Write canvases
c_hits_FI31->Write();
c_hits_FI32->Write();
c_hits_FI33->Write();
c_hits_MWPC1->Write();
c_hits_TOFD->Write();
c_hits_LOS->Write();
// c1->Write();
c3->Write();
c_area->Write();
c_MWPC1_FI32_X_single_corr->Write();
c_FI33_FI32_X->Write();
c_FI31_FI32_X->Write();
c_FIB31_FIB33_combined_FIB32_area->Write();
c_MWPC1_FI31_FI33_X_single_corr->Write();

// Write histograms
h_MWPC1_X->Write();
h_FI31_X->Write();
h_FI32_X->Write();
h_FI33_X->Write();
h_FICOMB_X->Write();
h_TOFD_X->Write();
h_hits_MWPC1->Write();
h_hits_FI31->Write();
h_hits_FI32->Write();
h_hits_FI33->Write();
h_hits_TOFD->Write();
h_hits_LOS->Write();
h_FI33_FI32_X->Write();
h_FI31_FI32_X->Write();
h_FI33_FI31_combined_FI32_X->Write();
h_FIB31_FIB33_combined_FIB32_area->Write();
h_MWPC1_FI31_FI33_X_single_corr->Write();
h_MWPC1_FI32_X_single_corr->Write();
h_PoQ->Write();
h_PoQvsQ->Write();
h_incidence_angle->Write();

rootFile->Close();

// Save PNG output
c_hits_FI31->SaveAs(Form("%s/hits_FI31.png", outputPngDir.Data()));
c_hits_FI32->SaveAs(Form("%s/hits_FI32.png", outputPngDir.Data()));
c_hits_FI33->SaveAs(Form("%s/hits_FI33.png", outputPngDir.Data()));
c_hits_MWPC1->SaveAs(Form("%s/hits_MWPC1.png", outputPngDir.Data()));
c_hits_TOFD->SaveAs(Form("%s/hits_TOFD.png", outputPngDir.Data()));
c_hits_LOS->SaveAs(Form("%s/hits_LOS.png", outputPngDir.Data()));
c3->SaveAs(Form("%s/FI33_FI31_combined_FI32_correlation.png", outputPngDir.Data()));
c_FI31_FI32_X->SaveAs(Form("%s/FI31_FI32_X.png", outputPngDir.Data()));
c_FIB31_FIB33_combined_FIB32_area->SaveAs(Form("%s/FIB31_FIB33_combined_FIB32_area.png", outputPngDir.Data()));
c_MWPC1_FI31_FI33_X_single_corr->SaveAs(Form("%s/MWPC1_FI31_FI33_X_single_corr.png", outputPngDir.Data()));
c_area->SaveAs(Form("%s/area.png", outputPngDir.Data()));

detectorFile->cd();
detectorTree->Write();
detectorFile->Close();

}
