
char histname[200];
char histname_copy[200];
char fname[200];

TH1F* plane_bar[4][44];
TH1F* plane_bar_copy[4][44];

TH1F* plane_bar_rec[4][44];

int plane_npeaks[4][44];
int plane_npeaks_copy[4][44];

Double_t* plane_peak_positions[4][44];
Double_t* plane_peak_positions_copy[4][44];
Double_t* plane_peak_positions_copy_integers[4][44];

int nbins_Q = 200;
int max_charge = 10;

TH2F* bar_Eloss[4];
TH2F* bar_Eloss_rec[4];
TH2F* bar_Eloss_rec_quad[4];

double intercept[4][44];
double slope[4][44];

double quad_c0[4][44];
double quad_c1[4][44];
double quad_c2[4][44];

void ToFD_crystal_calibration_v2(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 
    TChain* chain = new TChain("evt");
    sprintf(fname, "./main0252_0001_unpacked_califa_ToFD_LOS.root");

    chain->Reset();
    chain->Add(fname); 
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

    bar_Eloss_rec_quad[0] = new TH2F("Eloss plane 1 recalibrated quadratic", "Eloss plane 1 recalibrated quadratic", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss_rec_quad[1] = new TH2F("Eloss plane 2 recalibrated quadratic", "Eloss plane 2 recalibrated quadratic", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss_rec_quad[2] = new TH2F("Eloss plane 3 recalibrated quadratic", "Eloss plane 3 recalibrated quadratic", 44, 0, 45, nbins_Q, 0, max_charge);
    bar_Eloss_rec_quad[3] = new TH2F("Eloss plane 4 recalibrated quadratic", "Eloss plane 4 recalibrated quadratic", 44, 0, 45, nbins_Q, 0, max_charge);


    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            sprintf(histname,"ToFD Plane: %i Bar: %i", a+1, b+1);
            plane_bar[a][b] = new TH1F(histname, histname, nbins_Q, 0, 8);
            plane_bar[a][b]->GetXaxis()->SetTitle("Q");
            plane_bar[a][b]->GetYaxis()->SetTitle("Counts"); 
        }

        bar_Eloss[a]->GetXaxis()->SetTitle("Bar ID");
        bar_Eloss[a]->GetYaxis()->SetTitle("Q"); 

        bar_Eloss_rec[a]->GetXaxis()->SetTitle("Bar ID");
        bar_Eloss_rec[a]->GetYaxis()->SetTitle("Q"); 

        bar_Eloss_rec_quad[a]->GetXaxis()->SetTitle("Bar ID");
        bar_Eloss_rec_quad[a]->GetYaxis()->SetTitle("Q"); 
    }   

    /////////////////////////////////////////////////////////////////////////////////////////
    // Getting the data and filling the histograms

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
                
                if (Eloss >= 0.4){
                
                    plane_bar[TOFD_Id-1][bar-1]->Fill(Eloss); 

                    switch (TOFD_Id) {
                        case 1:
                            bar_Eloss[0]->Fill(bar, Eloss);
                        case 2:
                            bar_Eloss[1]->Fill(bar, Eloss);       
                        case 3:
                            bar_Eloss[2]->Fill(bar, Eloss);  
                        case 4:
                            bar_Eloss[3]->Fill(bar, Eloss);  
                    }
                }
            }
        }
    }
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Finding peaks of original data

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            auto spectrum = new TSpectrum();
            spectrum->Search(plane_bar[a][b]);

            int npeaks = spectrum->GetNPeaks();
            Double_t* xpos = spectrum->GetPositionX();
            sort(xpos, xpos + npeaks);

            plane_npeaks[a][b] = npeaks;
            plane_peak_positions[a][b] = xpos;
        }
    }   

    /////////////////////////////////////////////////////////////////////////////////////////
    // Finding peaks of convoluted  

    // Convolution function
    auto ideal_spec = new TF1("f6", "[0]*exp(-0.5*((x-1)/0.26)^2) + [1]*exp(-0.5*((x-2)/0.26)^2) +  [2]*exp(-0.5*((x-3)/0.26)^2) +  [3]*exp(-0.5*((x-4)/0.26)^2) +  [4]*exp(-0.5*((x-5)/0.26)^2) +  [5]*exp(-0.5*((x-6)/0.26)^2)", 0, 7);
    ideal_spec->SetParameters(0.5, 0.7, 1, 1, 0.8, 0.5); 

    // Number of convolutions
    int n_conv = 1;

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            // Creates a copy
            sprintf(histname_copy,"Copy of ToFD Plane: %i Bar: %i copy", a+1, b+1);
            plane_bar_copy[a][b] = new TH1F(histname_copy, histname_copy, nbins_Q, 0, 8);

            // Copies all entries
            for (int bin = 1; bin <= plane_bar[a][b]->GetNbinsX(); ++bin) {
                double content = plane_bar[a][b]->GetBinContent(bin);
                double error = plane_bar[a][b]->GetBinError(bin);
                plane_bar_copy[a][b]->SetBinContent(bin, content);
            }

            // Convolutes the original histogram
            for (int i=0; i < n_conv; i++) {
                plane_bar_copy[a][b]->Multiply(ideal_spec);
                double max_bin = plane_bar_copy[a][b]->GetMaximum();
                plane_bar_copy[a][b]->Scale(1.0 / max_bin);
            }

            // Peak finding
            plane_bar_copy[a][b]->Smooth(3);
            auto spectrum = new TSpectrum();
            spectrum->Search(plane_bar_copy[a][b]);

            int npeaks = spectrum->GetNPeaks();
            Double_t* xpos = spectrum->GetPositionX();
            sort(xpos, xpos + npeaks);

            plane_npeaks_copy[a][b] = npeaks;
            plane_peak_positions_copy[a][b] = xpos;
        }
    }   

    /////////////////////////////////////////////////////////////////////////////////////////
    // Finding discrepancies 

    int total_discr = 0;

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            if (plane_npeaks[a][b]!= plane_npeaks_copy[a][b]){
                cout << "Different number of peaks at plane " << a+1 << " bar " << b+1 << "\n";
                total_discr++;
            }
        }
    }
    
    cout << "Total amount of discrepancies: " << total_discr << "\n";

    /////////////////////////////////////////////////////////////////////////////////////////
    // Partial calibration with the non-discrepant peaks

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            plane_peak_positions_copy_integers[a][b] = new Double_t[plane_npeaks_copy[a][b]];

            // No discrepancy between the number of peaks of the original and convoluted
            if (plane_npeaks[a][b] == plane_npeaks_copy[a][b]){

                // Integer values of the convoluted peaks
                for (int i = 0; i < plane_npeaks_copy[a][b]; ++i) {
                    plane_peak_positions_copy_integers[a][b][i] = static_cast<Double_t>(std::round(plane_peak_positions_copy[a][b][i]));
                }

                /////////////////////////////////////////////////////////////////////////////////////////
                // Prints of peaks of original, convoluted and ideal
                cout << "Plane " << a+1 << " bar " << b+1 << "," << " number of peaks scaled: " << plane_npeaks_copy[a][b] << " number of original peaks: " << plane_npeaks[a][b] << "\n";

                for (int i = 0; i < plane_npeaks[a][b]; ++i) {
                    cout << plane_peak_positions[a][b][i] << "  ";
                }
                cout << "\n";

                for (int i = 0; i < plane_npeaks_copy[a][b]; ++i) {
                    cout << plane_peak_positions_copy[a][b][i] << "  ";
                }
                cout << "\n";

                for (int i = 0; i < plane_npeaks_copy[a][b]; ++i) {
                    cout << plane_peak_positions_copy_integers[a][b][i] << "  ";
                }
                cout << "\n";
                cout << "\n";

                /////////////////////////////////////////////////////////////////////////////////////////
                // Check if peaks are repeated or missassigned

                int missasigned = 0;
                int repeated = 0;

                for (int i = 0; i < plane_npeaks[a][b]; ++i) {
                    if ((plane_peak_positions[a][b][i] - plane_peak_positions_copy_integers[a][b][i]) > 1.0) {
                        cout << "Missasignment at plane " << a+1 << " bar " << b+1 << "\n";
                        missasigned = 1;
                    }
                }

                for (int i = 0; i < plane_npeaks_copy[a][b]; ++i) {
                    for (int j = i + 1; j < plane_npeaks_copy[a][b]; ++j) {
                        if (plane_peak_positions_copy_integers[a][b][i] == plane_peak_positions_copy_integers[a][b][j]) {
                            cout << "Repeated value " << plane_peak_positions_copy_integers[a][b][i] << " at plane " << a+1 << " bar " << b+1 << " entry " << plane_peak_positions_copy_integers[a][b][i] << endl;
                            repeated=1;
                        }
                    }
                }

                /////////////////////////////////////////////////////////////////////////////////////////


                if (missasigned == 0 && repeated == 0){
                    // Points to fit original peak position to 
                    TGraph graph(plane_npeaks[a][b], plane_peak_positions[a][b], plane_peak_positions_copy_integers[a][b]);

                    // Linear fit
                    graph.Fit("pol1", "Q");
                    TF1* fit = graph.GetFunction("pol1");
                    intercept[a][b] = fit->GetParameter(0); 
                    slope[a][b]     = fit->GetParameter(1); 

                    // Quadratic fit
                    graph.Fit("pol2", "Q");
                    TF1* fit2 = graph.GetFunction("pol2");
                    quad_c0[a][b] = fit2->GetParameter(0);
                    quad_c1[a][b] = fit2->GetParameter(1);
                    quad_c2[a][b] = fit2->GetParameter(2);
                }               
                else{
                    quad_c0[a][b] = 0.000;
                    quad_c1[a][b] = 0.000;
                    quad_c2[a][b] = 0.000;
                }   
            }      
        }
    }
    /////////////////////////////////////////////////////////////////////////////////////////
    // Recalibration 

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            sprintf(histname,"ToFD Plane recalibrated: %i Bar: %i", a+1, b+1);
            plane_bar_rec[a][b] = new TH1F(histname, histname, 100, 0, 8);
            plane_bar_rec[a][b]->GetXaxis()->SetTitle("Q");
            plane_bar_rec[a][b]->GetYaxis()->SetTitle("Counts"); 

            for (int bin = 1; bin <= plane_bar[a][b]->GetNbinsX(); ++bin) {
                
                double x_bin = plane_bar[a][b]->GetBinCenter(bin);
                double x_new = quad_c0[a][b] + quad_c1[a][b]*x_bin + quad_c2[a][b]*x_bin*x_bin;
                double content = plane_bar[a][b]->GetBinContent(bin);
                plane_bar_rec[a][b]->Fill(x_new, content);
                plane_bar_rec[a][b]->SetBinError(bin, 0.0);
            }
        }
    } 

    for (int plane = 0; plane <= 3; ++plane) {
        for (int bar = 1; bar <= 44; ++bar) {
            for (int qbin = 1; qbin <= nbins_Q; ++qbin) {

                double content = bar_Eloss[plane]->GetBinContent(bar, qbin);
                double Q = bar_Eloss[plane]->GetYaxis()->GetBinCenter(qbin);
                int a = plane;
                int b = bar - 1;
                double Q_new = quad_c0[a][b] + quad_c1[a][b]*Q + quad_c2[a][b]*Q*Q;

                if (quad_c0[a][b] != 0 && quad_c1[a][b] != 0 && quad_c2[a][b] != 0){
                    bar_Eloss_rec[plane]->Fill(bar, Q_new, content);
                } 
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    int cry_1 = 0;
    int cry_2 = 20;

    // 1 - 42 
    TCanvas* c_npeaksone = new TCanvas();
    c_npeaksone->cd(1);
    plane_bar[cry_1][cry_2]->Draw("");

    TCanvas* c_npeaksone_copy = new TCanvas();
    c_npeaksone_copy->cd(1);
    plane_bar_copy[cry_1][cry_2]->Draw("");

    TCanvas* c_npeaksone_copy_rec = new TCanvas();
    c_npeaksone_copy_rec->cd(1);
    plane_bar_rec[cry_1][cry_2]->Draw("");


    // Original histograms
    TCanvas* c_bar_Eloss0 = new TCanvas();
    c_bar_Eloss0->cd(1);
    gPad->SetLogz();
    bar_Eloss[0]->Draw("colz");

    TCanvas* c_bar_Eloss1 = new TCanvas();
    c_bar_Eloss1->cd(1);
    gPad->SetLogz();
    bar_Eloss[1]->Draw("colz");

    TCanvas* c_bar_Eloss2 = new TCanvas();
    c_bar_Eloss2->cd(1);
    gPad->SetLogz();
    bar_Eloss[2]->Draw("colz");

    TCanvas* c_bar_Eloss3 = new TCanvas();
    c_bar_Eloss3->cd(1);
    gPad->SetLogz();
    bar_Eloss[3]->Draw("colz");


    // Recalibrated histograms
    TCanvas* c_bar_Eloss_rec0 = new TCanvas();
    c_bar_Eloss_rec0->cd(1);
    gPad->SetLogz();
    bar_Eloss_rec[0]->Draw("colz");

    TCanvas* c_bar_Eloss_rec1 = new TCanvas();
    c_bar_Eloss_rec1->cd(1);
    gPad->SetLogz();
    bar_Eloss_rec[1]->Draw("colz");

    TCanvas* c_bar_Eloss_rec2 = new TCanvas();
    c_bar_Eloss_rec2->cd(1);
    gPad->SetLogz();
    bar_Eloss_rec[2]->Draw("colz");

    TCanvas* c_bar_Eloss_rec3 = new TCanvas();
    c_bar_Eloss_rec3->cd(1);
    gPad->SetLogz();
    bar_Eloss_rec[3]->Draw("colz");




}


// double quad_constraint(double* x, double* par) {
//     double x1 = ...; // first x value
//     double x2 = ...; // second x value
//     double y1 = 1;   // y at x1
//     double y2 = 6;   // y at x2
//     double c2 = par[0];
//     double c1 = (y1 - y2 - c2*(x1*x1 - x2*x2)) / (x1 - x2);
//     double c0 = y1 - c1*x1 - c2*x1*x1;
//     return c0 + c1*x[0] + c2*x[0]*x[0];
// }