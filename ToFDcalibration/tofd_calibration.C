
char fname[200];
char histname[200];
TH1F* plane_bar[4][44];

void tofd_calibration() {

    // Load the root file 
    TChain* chain = new TChain("evt");
    //sprintf(fname, "./main0252_0001_unpacked_025_final.root");
    sprintf(fname, "./main0252_0001_unpacked_califa_ToFD_LOS.root");

    chain->Reset();
    chain->Add(fname); 

    auto bar_Eloss1  = new TH2F("Eloss plane 1","Eloss plane 1", 44, 0, 45, 200, 1, 10 );
    auto bar_Eloss2  = new TH2F("Eloss plane 2","Eloss plane 2", 44, 0, 45, 200, 1, 10 );
    auto bar_Eloss3  = new TH2F("Eloss plane 3","Eloss plane 3", 44, 0, 45, 200, 1, 10 );
    auto bar_Eloss4  = new TH2F("Eloss plane 4","Eloss plane 4", 44, 0, 45, 200, 1, 10 );
    auto bar_Eloss   = new TH1F("Eloss ToFD","Eloss ToFD", 200, 0, 8);
    auto bar_Eloss_a = new TH2F("Eloss all planes","Eloss all planes", 44, 0, 45, 200, 1, 10 );
    auto h_npeaks    = new TH1F("Number of peaks","Number of peaks", 8, 0.5, 8.5);
    auto value_Eloss = new TH2F("Value peak vs number peaks","Value peak vs number peaks", 100, 0, 8, 8, 0.5, 8.5 );

    bar_Eloss1->GetXaxis()->SetTitle("Bar ID");
    bar_Eloss1->GetYaxis()->SetTitle("Q"); 

    bar_Eloss2->GetXaxis()->SetTitle("Bar ID");
    bar_Eloss2->GetYaxis()->SetTitle("Q"); 

    bar_Eloss3->GetXaxis()->SetTitle("Bar ID");
    bar_Eloss3->GetYaxis()->SetTitle("Q");

    bar_Eloss4->GetXaxis()->SetTitle("Bar ID");
    bar_Eloss4->GetYaxis()->SetTitle("Q"); 

    bar_Eloss->GetXaxis()->SetTitle("Eloss");

    bar_Eloss_a->GetXaxis()->SetTitle("Bar ID");
    bar_Eloss_a->GetYaxis()->SetTitle("Q"); 

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {
            sprintf(histname,"ToFD Plane: %i Bar: %i",a+1, b+1);
            plane_bar[a][b] = new TH1F(histname,histname,100,0,8);
            plane_bar[a][b]->GetXaxis()->SetTitle("Q");
            plane_bar[a][b]->GetYaxis()->SetTitle("Counts"); 
        }
    }   

    h_npeaks->GetXaxis()->SetTitle("Eloss");
    h_npeaks->GetYaxis()->SetTitle("Counts"); 

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 

	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("EventHeader.");
    branchData->SetAddress(&header);    

	TClonesArray*  hitTOFD = new TClonesArray("R3BTofdHitData",5);
	R3BTofdHitData** hitTofd;
	TBranch *branchHitTOFD = chain->GetBranch("TofdHit");
	branchHitTOFD->SetAddress(&hitTOFD);
	Long64_t TofdHitPerEvent = 0;  	

    Long64_t nevent = chain->GetEntries();
	cout << "Number of events in " << fname << " : " << nevent << endl;

    ///////////////////////////////////////////////////////////////////////////////////////

	for(Long64_t e=0; e<nevent; e++){  

	    hitTOFD->Clear();
	    chain->GetEvent(e);
	    TofdHitPerEvent = hitTOFD->GetEntriesFast();

        if( TofdHitPerEvent > 0 ){

            hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
                
            for(int p = 0; p < TofdHitPerEvent; p++){

                hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
                
                int TOFD_Id = hitTofd[p]->GetDetId();     
                double Eloss = hitTofd[p]->GetEloss();
                double TOFD_X = hitTofd[p]->GetX();
                double TOFD_Y = hitTofd[p]->GetY();
                int bar = hitTofd[p]->GetBarId();                

                switch (TOFD_Id) {
                    case 1:
                        bar_Eloss1->Fill(bar, Eloss);
                    case 2:
                        bar_Eloss2->Fill(bar, Eloss);       
                    case 3:
                        bar_Eloss3->Fill(bar, Eloss);  
                    case 4:
                        bar_Eloss4->Fill(bar, Eloss);  
                }

                bar_Eloss->Fill(Eloss);
                bar_Eloss_a->Fill(bar,Eloss);
                plane_bar[TOFD_Id-1][bar-1]->Fill(Eloss);
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            auto spectrum = new TSpectrum();
            spectrum->Search(plane_bar[a][b]);
            int npeaks = spectrum->GetNPeaks();
            Double_t* xpos = spectrum->GetPositionX();
            h_npeaks->Fill(npeaks);
            if (npeaks == 4) cout<< "Plane: " << a+1 << " Bar: " << b+1 << " :" << spectrum->GetNPeaks() << "\n";

            for (int j=0; j < npeaks; j++){
                value_Eloss->Fill(xpos[j], npeaks);
            }
        }
    }   
    
    /////////////////////////////////////////////////////////////////////////////////////////

    auto spectrum2 = new TSpectrum();
    spectrum2->Search(bar_Eloss);
    int npeaks2 = spectrum2->GetNPeaks();
    cout<< "All bars peaks: " << npeaks2 << "\n";

    ///////////////////////////////////////////////////////////////////////////////////////

    // Background
    //auto ideal_spec = new TF1("f6", "[2]*exp(-0.5*((x-3)/0.15)^2)+(x-3)**2*([0]*exp(-0.5*((x-1)/0.15)^2) +  [1]*exp(-0.5*((x-2)/0.15)^2) +  [2]*exp(-0.5*((x-3)/0.15)^2) +  [3]*exp(-0.5*((x-4)/0.15)^2) +  [4]*exp(-0.5*((x-5)/0.15)^2) +  [5]*exp(-0.5*((x-6)/0.15)^2))", 0, 7);
    // No background
    auto ideal_spec = new TF1("f6", "[0]*exp(-0.5*((x-1)/0.15)^2) + [1]*exp(-0.5*((x-2)/0.15)^2) +  [2]*exp(-0.5*((x-3)/0.15)^2) +  [3]*exp(-0.5*((x-4)/0.15)^2) +  [4]*exp(-0.5*((x-5)/0.15)^2) +  [5]*exp(-0.5*((x-6)/0.15)^2)", 0, 7);
    ideal_spec->SetParameters(500, 800, 1000, 1000, 800, 500); 
    ideal_spec->Draw();

    int cry_1 = 1;
    int cry_2 = 44;

    // 1 - 42 

    TCanvas* c_npeaksone = new TCanvas();
    c_npeaksone->cd(1);
    plane_bar[cry_1][cry_2]->Draw("");

    TCanvas* c_npeaksone_prod = new TCanvas();
    c_npeaksone_prod->cd(1);

    auto plane_bar_copy = (TH1F*)plane_bar[cry_1][cry_2]->Clone("plane_bar1_42_copy");
    plane_bar_copy->Multiply(ideal_spec);
    double max_bin = plane_bar_copy->GetMaximum();
    plane_bar_copy->Scale(1.0 / max_bin);

    plane_bar_copy->Multiply(ideal_spec);

    plane_bar_copy->Multiply(ideal_spec);
    //double max_bin2 = plane_bar_copy->GetMaximum();
    //plane_bar_copy->Scale(1.0 / max_bin2);

    auto spectrum3 = new TSpectrum();
    spectrum3->Search(plane_bar_copy);
    int npeaks3 = spectrum3->GetNPeaks();
    cout << "New peaks" << npeaks3 << "\n";

    c_npeaksone_prod->Clear();
    plane_bar_copy->Draw("");

    ///////////////////////////////////////////////////////////////////////////////////////

    // TCanvas* c_bar_Eloss1 = new TCanvas();
    // c_bar_Eloss1->cd(1);
    // bar_Eloss1->Draw("colz");

    // TCanvas* c_bar_Eloss2 = new TCanvas();
    // c_bar_Eloss2->cd(1);
    // bar_Eloss2->Draw("colz");

    // TCanvas* c_bar_Eloss3 = new TCanvas();
    // c_bar_Eloss3->cd(1);
    // bar_Eloss3->Draw("colz");

    // TCanvas* c_bar_Eloss4 = new TCanvas();
    // c_bar_Eloss4->cd(1);
    // bar_Eloss4->Draw("colz");

    // TCanvas* c_bar_Eloss = new TCanvas();
    // c_bar_Eloss->cd(1);
    // bar_Eloss->Draw("");

    // TCanvas* c_bar_Eloss_a = new TCanvas();
    // c_bar_Eloss_a->cd(1);
    // bar_Eloss_a->Draw("colz");

    // TCanvas* c_npeaks = new TCanvas();
    // c_npeaks->cd(1);
    // h_npeaks->Draw("");

    // TCanvas* c_npeaksone = new TCanvas();
    // c_npeaksone->cd(1);
    // plane_bar[1][42]->Draw("");

    // TCanvas* c_value_Eloss = new TCanvas();
    // c_value_Eloss->cd(1);
    // value_Eloss->Draw("colz");

}


