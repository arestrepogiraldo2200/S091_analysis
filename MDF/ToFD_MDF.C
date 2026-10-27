
#include <fstream>
#include <algorithm>
#include <cmath>

/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms
char fname1[200], fname2[200], fname3[200], fname4[200], fname5[200], fname6[200], fname7[200], fname8[200], fname9[200], fname10[200], fname11[200], fname12[200], fname13[200];
int tpat;
double recalib_tofd_c0[4][44] = {{0.0}};
double recalib_tofd_c1[4][44] = {{0.0}};
double recalib_tofd_c2[4][44] = {{0.0}};

TH2F* bar_Eloss[4];
TH2F* bar_Eloss_rec[4];

// Loads recalibration parameters of TOFD
bool LoadCalibrationMatrix(const TString& fileName, double values[4][44]){
    std::ifstream input(fileName.Data());
    if (!input.is_open()) {
        std::cout << "Failed to open " << fileName << " for reading" << std::endl;
        return false;
    }

    for (int plane = 0; plane < 4; ++plane) {
        for (int bar = 0; bar < 44; ++bar) {
            if (!(input >> values[plane][bar])) {
                std::cout << "Failed to read " << fileName << " at row " << plane << " column " << bar << std::endl;
                return false;
            }
        }
    }

    std::cout << "Loaded " << fileName << std::endl;
    return true;
}

///////////////////////////////////////////////////////////////
void ToFD_MDF(){

    int max_charge =8;
    int nbins_Q = 2000;

    /////////////////////////////////////////////////////////////////////////////////////////
    // Setting titles 
    bar_Eloss[0] = new TH2F("Eloss plane 1", "Eloss plane 1", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss[1] = new TH2F("Eloss plane 2", "Eloss plane 2", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss[2] = new TH2F("Eloss plane 3", "Eloss plane 3", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss[3] = new TH2F("Eloss plane 4", "Eloss plane 4", 44, 0, 45, nbins_Q, 0, max_charge);

    bar_Eloss_rec[0] = new TH2F("Eloss plane 1 recalibrated", "Eloss plane 1 recalibrated", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss_rec[1] = new TH2F("Eloss plane 2 recalibrated", "Eloss plane 2 recalibrated", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss_rec[2] = new TH2F("Eloss plane 3 recalibrated", "Eloss plane 3 recalibrated", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss_rec[3] = new TH2F("Eloss plane 4 recalibrated", "Eloss plane 4 recalibrated", 44, 0, 45, nbins_Q, 0, max_charge);

    for(int a=0; a<4; a++) {
        bar_Eloss[a]->GetXaxis()->SetTitle("Bar ID");
        bar_Eloss[a]->GetYaxis()->SetTitle("Q"); 

        bar_Eloss_rec[a]->GetXaxis()->SetTitle("Bar ID");
        bar_Eloss_rec[a]->GetYaxis()->SetTitle("Q"); 
    }   
    /////////////////////////////////////////////////////////////////////////////////////////

    // TOFD recalibration
    LoadCalibrationMatrix("./recalib_tofd_c0.txt", recalib_tofd_c0);
    LoadCalibrationMatrix("./recalib_tofd_c1.txt", recalib_tofd_c1);
    LoadCalibrationMatrix("./recalib_tofd_c2.txt", recalib_tofd_c2);

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt"); 
    sprintf(fname1,"./main0274_0001_unpacked_los.root");
    sprintf(fname2,"./main0275_0001_unpacked_los.root");
    sprintf(fname3,"./main0276_0001_unpacked_los.root");
    sprintf(fname4,"./main0277_0001_unpacked_los.root");
    sprintf(fname5,"./main0278_0001_unpacked_los.root");
    sprintf(fname6,"./main0279_0001_unpacked_los.root");
    sprintf(fname7,"./main0280_0001_unpacked_los.root");
    sprintf(fname8,"./main0281_0001_unpacked_los.root");
    sprintf(fname9,"./main0282_0001_unpacked_los.root");
    sprintf(fname10,"./main0283_0001_unpacked_los.root");
    sprintf(fname11,"./main0284_0001_unpacked_los.root");
    sprintf(fname12,"./main0285_0001_unpacked_los.root");
    sprintf(fname13,"./main0286_0001_unpacked_los.root");

    chain->Reset();
    chain->Add(fname1); 
    chain->Add(fname2); 
    chain->Add(fname3); 
    chain->Add(fname4); 
    chain->Add(fname5); 
    chain->Add(fname6); 
    chain->Add(fname7); 
    chain->Add(fname8); 
    chain->Add(fname9); 
    chain->Add(fname10); 
    chain->Add(fname11); 
    chain->Add(fname12); 
    chain->Add(fname13); 

    // Load the branches with data as TClonesArrays 
	R3BEventHeader* header = new R3BEventHeader;    
    chain->SetBranchAddress("EventHeader.", &header);

	TClonesArray* hitTOFD = new TClonesArray("R3BTofdHitData",5);
	R3BTofdHitData** hitTofd;
	chain->SetBranchAddress("TofdHit", &hitTOFD);
	Long64_t TofdHitPerEvent = 0;  	

	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in the tree: " << nevent << endl;

	for(Long64_t e=0; e < 100000; e++){

        // Cleaning the vectors 
        hitTOFD->Clear();
        chain->GetEvent(e);

        tpat              = header->GetTpat();
        TofdHitPerEvent   = hitTOFD->GetEntriesFast();

        if (TofdHitPerEvent == 4){

            /////////////////////////////////////////////////////////////////////////////////////////
            // TOFD charge gate selection
            bool b_tofdq = false;
            double averageEloss = 0.;

            hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
            double totalEloss = 0.;

            for(int p = 0; p < TofdHitPerEvent; p++){

                hitTofd[p] = (R3BTofdHitData*)hitTOFD->At(p);

                int TOFD_Id  = hitTofd[p]->GetDetId();     
                int bar      = hitTofd[p]->GetBarId();   
                double Eloss = hitTofd[p]->GetEloss();
                double Eloss_new = recalib_tofd_c0[TOFD_Id-1][bar-1] + recalib_tofd_c1[TOFD_Id-1][bar-1]*Eloss + recalib_tofd_c2[TOFD_Id-1][bar-1]*Eloss*Eloss;

                switch (TOFD_Id) {
                    case 1:
                        bar_Eloss[0]->Fill(bar, Eloss);
                        bar_Eloss_rec[0]->Fill(bar, Eloss_new);
                    case 2:
                        bar_Eloss[1]->Fill(bar, Eloss);
                        bar_Eloss_rec[1]->Fill(bar, Eloss_new);
                    case 3:
                        bar_Eloss[2]->Fill(bar, Eloss);
                        bar_Eloss_rec[2]->Fill(bar, Eloss_new);
                    case 4:
                        bar_Eloss[3]->Fill(bar, Eloss);
                        bar_Eloss_rec[3]->Fill(bar, Eloss_new);
                }
            }
        }         
    }            

    // Original histograms
    TCanvas* c_bar_Eloss0 = new TCanvas();
    c_bar_Eloss0->cd(1);
    bar_Eloss[0]->Draw("colz");

    TCanvas* c_bar_Eloss1 = new TCanvas();
    c_bar_Eloss1->cd(1);
    bar_Eloss[1]->Draw("colz");

    TCanvas* c_bar_Eloss2 = new TCanvas();
    c_bar_Eloss2->cd(1);
    bar_Eloss[2]->Draw("colz");

    TCanvas* c_bar_Eloss3 = new TCanvas();
    c_bar_Eloss3->cd(1);
    bar_Eloss[3]->Draw("colz");

    // Recalibrated histograms
    TCanvas* c_bar_Eloss_rec0 = new TCanvas();
    c_bar_Eloss_rec0->cd(1);
    bar_Eloss_rec[0]->Draw("colz");

    TCanvas* c_bar_Eloss_rec1 = new TCanvas();
    c_bar_Eloss_rec1->cd(1);
    bar_Eloss_rec[1]->Draw("colz");

    TCanvas* c_bar_Eloss_rec2 = new TCanvas();
    c_bar_Eloss_rec2->cd(1);
    bar_Eloss_rec[2]->Draw("colz");

    TCanvas* c_bar_Eloss_rec3 = new TCanvas();
    c_bar_Eloss_rec3->cd(1);
    bar_Eloss_rec[3]->Draw("colz");
}
