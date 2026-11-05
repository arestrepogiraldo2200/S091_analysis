// Components speed
#include <cmath>
#include <TF1.h>

TH2D* h_Califa_Components;
TH1D* h_Califa_hits;
TH1D* h_Califa_Ns;
TH1D* h_Califa_Nf;
TH2D* h_Califa_Ns_Nf;
TH1D* h_Califa_CC;

TH1D* h_CALIFAE;
TH2D* h_CALIFAECORR;
TH2D* h_CALIFATHETA;
TH1D* h_CALIFATHETASUM;
TH2D* h_CALIFAPHI;

#include "TH2D.h"
#include "TH1D.h"
#include "TSpectrum.h"
#include <vector>
#include <algorithm>

// Function to slice the histogram
void SlicePeaks(TH2D* h2, std::vector<double>& xs, std::vector<double>& ys, int binsPerSlice = 100, double xmin = 0, double xmax = 900000, double minEntries = 1){
    
    xs.clear(); ys.clear();
    int b0 = h2->GetXaxis()->FindBin(xmin), b1 = h2->GetXaxis()->FindBin(xmax);

    for (int ix = b0; ix + binsPerSlice - 1 <= b1; ix += binsPerSlice){
        int ixe = ix + binsPerSlice - 1;
        double xc = 0.5 * (h2->GetXaxis()->GetBinLowEdge(ix) + h2->GetXaxis()->GetBinUpEdge(ixe));

        TH1D* p = h2->ProjectionY(Form("slice_%d", ix), ix, ixe);   // Vertical slice
        if (p->GetEntries() < minEntries) { delete p; continue; }
        p->Rebin(2);

        TSpectrum s(5);
        int n = s.Search(p, 2, "nodraw", 0.1);                     // sigma, option, threshold (fraction of max)
        if (n > 0) {
            double* px = s.GetPositionX();
            xs.push_back(xc);
            ys.push_back(*std::max_element(px, px + n));          // peak with highest Y

            //cout << "Peak found" << endl;
            cout << xc << " " << *std::max_element(px, px + n) << endl; 

        } else {

            //cout << "No peak found!!!!!!!!!!!!!" << endl;
            double pk = p->GetBinCenter(p->GetMaximumBin());
            double rms = p->GetRMS();

            if (rms <= 0) rms = p->GetBinWidth(1) * 5;
            TF1 g("g_slice", "gaus", pk - 2 * rms, pk + 2 * rms);
            g.SetParameters(p->GetMaximum(), pk, rms);
            int status = p->Fit(&g, "QN0R");
            double mean = g.GetParameter(1);

            xs.push_back(xc);
            // Last resort: if the fit failed, use the mean of all entries in the slice
            double yval = (status == 0 && std::isfinite(mean)) ? mean : p->GetMean();
            ys.push_back(yval);

            cout << xc << " " << yval << endl; 
        }
        delete p;  
    }
}

// When getting a specific crystal
int target_crystal = 4272; 
// int target_crystal = 4300; 
// int target_crystal = -1; 

char fname1[200], fname2[200], fname3[200], fname4[200], fname5[200], fname6[200], fname7[200], fname8[200], fname9[200], fname10[200], fname11[200], fname12[200], fname13[200];

void CalifaPID(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Califa 
    h_Califa_Components = new TH2D("Califa fast and slow components","Califa fast and slow components",5000,-50000,900000, 200,-0.5,1);
    h_Califa_Components->GetXaxis()->SetTitle("Nf+Ns");
    h_Califa_Components->GetYaxis()->SetTitle("Nf/(Nf+Ns)");

    h_Califa_hits = new TH1D("Hits per event CALIFA", "Hits per event CALIFA",20, 0, 20);
    h_Califa_hits->GetXaxis()->SetTitle("Hits per event CALIFA");

    h_Califa_Ns = new TH1D("CALIFA Ns", "CALIFA Ns",2000, -300000, 700000);
    h_Califa_Ns->GetXaxis()->SetTitle("Hits per event CALIFA");

    h_Califa_Nf = new TH1D("CALIFA Nf", "CALIFA Nf",2000, -100000, 400000);
    h_Califa_Nf->GetXaxis()->SetTitle("Hits per event CALIFA");

    h_Califa_Ns_Nf = new TH2D("Califa fast vs slow components","Califa fast vs slow components",2000, -300000, 700000, 2000, -100000, 400000);
    h_Califa_Ns_Nf->GetXaxis()->SetTitle("Nf");
    h_Califa_Ns_Nf->GetYaxis()->SetTitle("Ns");

    h_Califa_CC = new TH1D("Califa mother crystal","Califa mother crystal",2000, 1500, 6000);
    h_Califa_CC->GetXaxis()->SetTitle("Crystal ID");

    h_CALIFAE = new TH1D("Energy","Energy",100,0,700000);
    h_CALIFAE->GetXaxis()->SetTitle("Energy");

    h_CALIFATHETA = new TH2D("Angular correlation","Angular correlation",100,0,2,100,0,2);
    h_CALIFATHETA->GetXaxis()->SetTitle("Angle 1");
    h_CALIFATHETA->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAPHI = new TH2D("Azimutal angular correlation","Azimutal angular correlation",100,-3.2,3.2,100,-3.2,3.2);
    h_CALIFAPHI->GetXaxis()->SetTitle("Angle 1");
    h_CALIFAPHI->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAECORR = new TH2D("Energy correlations","Energy correlations",100,0,700000,100,0,700000);
    h_CALIFAECORR->GetXaxis()->SetTitle("Energy 1");
    h_CALIFAECORR->GetYaxis()->SetTitle("Energy 2");

    h_CALIFATHETASUM = new TH1D("Opening angle","Opening angle",100,0,2);
    h_CALIFATHETASUM->GetXaxis()->SetTitle("Opening angle");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 
    TChain* chain = new TChain("evt"); 
    // sprintf(fname1,"./main0280_0001_unpacked_los.root");
    // sprintf(fname1,"./main0252_0001_unpacked_MDF_2.root");
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

    // sprintf(fname1,"./main0126_0001_unpacked_lo.root");

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

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 
	R3BEventHeader* header = new R3BEventHeader;    
    chain->SetBranchAddress("EventHeader.", &header);


    TClonesArray* hitCALIFA = new TClonesArray("R3BCalifaClusterData",5);
	R3BCalifaClusterData** clusterCalifa;
	chain->SetBranchAddress("CalifaClusterData", &hitCALIFA);
	Long64_t CalifaHitPerEvent = 0;  


    TClonesArray* calCALIFA = new TClonesArray("R3BCalifaCrystalCalData",5);
	R3BCalifaCrystalCalData** caldataCalifa;
	chain->SetBranchAddress("CalifaCrystalCalData", &calCALIFA);
	Long64_t CalifaCalPerEvent = 0;  

    
	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in the tree: " << nevent << endl;

    /////////////////////////////////////////////////////////////////////////////////////////
	for(Long64_t e=0; e < nevent; e++){

        hitCALIFA->Clear();
        calCALIFA->Clear();
        chain->GetEvent(e);

        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        h_Califa_hits->Fill(CalifaHitPerEvent);

        CalifaCalPerEvent = calCALIFA->GetEntriesFast();

        /////////////////////////////////////////////////////////////////////////////////////////
        clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];
        caldataCalifa = new R3BCalifaCrystalCalData*[CalifaCalPerEvent];

        // Tolerances for kinematical conditions
        double tolerance_theta = 0.2;
        double tolerance_phi = 0.2;

        // ----------------------------------------------------------------------------------------------------
        if (CalifaHitPerEvent == 2){                                         // Condition on two hits in Califa

            bool califaL = false;
            bool califaR = false;
            double thetaTL = 0, thetaTR = 0, phiTR = 0, phiTL = 0, EL = 0, ER = 0, NfL = 0, NfR = 0, NsL = 0, NsR = 0;
            int CCL = 0, CCR = 0;
            double thetasum_plane = 0;

            for(int p = 0; p < CalifaHitPerEvent; p++){

                clusterCalifa[p] = (R3BCalifaClusterData*)hitCALIFA->At(p);

                double E  = clusterCalifa[p]->GetEnergy();
                double Nf = clusterCalifa[p]->GetNf();
                double Ns = clusterCalifa[p]->GetNs();
                int CC = clusterCalifa[p]->GetMotherCrystal();
                double theta = clusterCalifa[p]->GetTheta();
                double phi   = clusterCalifa[p]->GetPhi();

                // ----------------------------------------------------------------------------------------------------
                if ( (phi < 1.5 && phi > -1.5) ){

                    califaL = true;
                    thetaTL = theta;
                    thetasum_plane += thetaTL;
                    phiTL = phi;
                    EL = E;
                    CCL = CC;
                    NfL = Nf;
                    NsL = Ns;
                }

                // ----------------------------------------------------------------------------------------------------
                if ( (phi>1.5 || phi<-1.5) ){ 
                                                    
                    califaR = true;
                    thetaTR = theta;
                    thetasum_plane += thetaTR;
                    phiTR = phi;
                    ER = E;
                    CCR = CC;
                    NfR = Nf;
                    NsR = Ns;
                }
            }

            // ----------------------------------------------------------------------------------------------------
            // Conditionals on energy and kinematics

            if ( califaR == true && califaL == true){
                if ( EL > 20000 && ER > 20000){
                    if (  abs( abs(phiTL) + abs(phiTR) - 3.1415926535)  < tolerance_phi) {
                        if (  abs( abs(thetaTL) + abs(thetaTR) - 1.41372)  < tolerance_theta) {

                            // Search over each crystal
                            for(int q = 0; q < CalifaCalPerEvent; q++){

                                caldataCalifa[q] = (R3BCalifaCrystalCalData*)calCALIFA->At(q);
                                int CCcry = caldataCalifa[q]->GetCrystalId();
                                double Ecry  = caldataCalifa[q]->GetEnergy();
                                double Nfcry = caldataCalifa[q]->GetNf();
                                double Nscry = caldataCalifa[q]->GetNs();

                                // ----------------------------------------------------------------------------------------------------
                                if ( (CCcry == CCR && Ecry >= 0.9*ER) || (CCcry == CCL && Ecry >= 0.9*EL)){        // Condition on more than 90% of the energy deposited in the mother crystal of both hits

                                    h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                    h_CALIFAPHI->Fill(phiTL,phiTR);
                                    h_CALIFATHETASUM->Fill(thetasum_plane);
                                    h_CALIFAECORR->Fill(EL,ER);
                                    h_CALIFAE->Fill(EL+ER);

                                    // If the target crystal is in the left .........................................................
                                    if (target_crystal == CCL){

                                        h_Califa_CC->Fill(CCL);
                                        h_Califa_Ns->Fill(NsL);
                                        h_Califa_Nf->Fill(NfL);
                                        h_Califa_Ns_Nf->Fill(NfL,NsL);

                                        if ( (NfL+NsL)>0 ){
                                            h_Califa_Components->Fill(NfL+NsL, NfL/(NfL+NsL));
                                        }
                                    }

                                    // If the target crystal is the right .........................................................
                                    else if (target_crystal == CCR){

                                        h_Califa_CC->Fill(CCR);
                                        h_Califa_Ns->Fill(NsR);
                                        h_Califa_Nf->Fill(NfR);
                                        h_Califa_Ns_Nf->Fill(NfR,NsR);

                                        if ( (NfR+NsR)>0 ){
                                            h_Califa_Components->Fill(NfR+NsR, NfR/(NfR+NsR));
                                        }
                                    }

                                    // All crystals ...............................................................................
                                    else if (target_crystal ==-1) {

                                        h_Califa_CC->Fill(CCR);
                                        h_Califa_Ns->Fill(NsR);
                                        h_Califa_Nf->Fill(NfR);
                                        h_Califa_Ns_Nf->Fill(NfR,NsR);

                                        h_Califa_CC->Fill(CCL);
                                        h_Califa_Ns->Fill(NsL);
                                        h_Califa_Nf->Fill(NfL);
                                        h_Califa_Ns_Nf->Fill(NfL,NsL);

                                        if ( (NfR+NsR)>0 ){
                                            h_Califa_Components->Fill(NfR+NsR, NfR/(NfR+NsR));
                                        }

                                        if ( (NfL+NsL)>0 ){
                                            h_Califa_Components->Fill(NfL+NsL, NfL/(NfL+NsL));
                                        }
                                    }
                                }
                            }   
                        }
                    }               
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    // Fitting of the Nf component
    double xMin = 30000;
    double xMax = 200000;
    double yMin = 0.4;
    double yMax = 0.7;

    TGraph* gr = new TGraph();

    // Vector to store the peak coordinates
    std::vector<double> peakX, peakY;
    SlicePeaks(h_Califa_Components, peakX, peakY, 50, xMin, xMax);

    cout << "Peaks array length: " <<  peakX.size() << endl << endl; 

    // Set the values in the X-Y plane
    for (size_t i = 0; i < peakX.size(); ++i) {
        if (peakY[i] < yMin || peakY[i] > yMax) continue;   
        gr->SetPoint(gr->GetN(), peakX[i], peakY[i]);
    }

    // Function to fit the peaks
    TF1 *f = new TF1("f",
                 "[0] + [1]/sqrt(x+[2]) + [3]/(x + [4]) + [5]/((x + [6])*(x + [6]))",
                 xMin, xMax);

    gr->Fit(f, "R");
    TCanvas* c_gr  = new TCanvas();
    gr->Draw();

    /////////////////////////////////////////////////////////////////////////////////////////

    TCanvas* c_Califa_hits  = new TCanvas();
    c_Califa_hits->cd(1);
    h_Califa_hits->Draw("hist");

    TCanvas* c_Califa_Ns  = new TCanvas();
    c_Califa_Ns->cd(1);
    h_Califa_Ns->Draw("hist");

    TCanvas* c_Califa_Nf  = new TCanvas();
    c_Califa_Nf->cd(1);
    h_Califa_Nf->Draw("hist");

    TCanvas* c_Califa_Ns_Nf = new TCanvas();
    c_Califa_Ns_Nf->cd(1);
    h_Califa_Ns_Nf->Draw("colz");

    TF1 *line1 = new TF1("line1", "[0]*x+3500", 0, 100000);
    line1->SetParameter(0,1.4);   
    line1->SetLineWidth(1);
    line1->Draw("same");

    TF1 *line2 = new TF1("line2", "[0]*x-13000", 0, 100000);
    line2->SetParameter(0, 1.4);   
    line2->SetLineWidth(1);
    line2->Draw("same");


    TCanvas* c_Califa_CC  = new TCanvas();
    c_Califa_CC->cd(1);
    h_Califa_CC->Draw("hist");

    // ===========================================================================================
    TCanvas* c_Califa_Components = new TCanvas();
    c_Califa_Components->cd(1);
    h_Califa_Components->Draw("colz");

    double drawMin = 1000;
    double drawMax = 700000;

    TGraph* g_fit = new TGraph();
    for (int i = 0; i <= 1000; ++i) {
        double x = drawMin + (drawMax - drawMin) * i / 1000.;
        double y = f->Eval(x);
        if (std::isfinite(y)) g_fit->SetPoint(g_fit->GetN(), x, y);
    }
    g_fit->SetLineColor(kRed);
    g_fit->SetLineWidth(1);
    g_fit->Draw("L same");

    TCanvas* c_CALIFAE = new TCanvas();
    c_CALIFAE->cd(1);
    h_CALIFAE->Draw("hist");

    TCanvas* c_CALIFAECORR = new TCanvas();
    c_CALIFAECORR->cd(1);
    h_CALIFAECORR->Draw("colz");

    TCanvas* c_CALIFATHETASUM = new TCanvas();
    c_CALIFATHETASUM->cd(1);
    h_CALIFATHETASUM->Draw("hist");


    TCanvas* c_CALIFATHETA = new TCanvas();
    c_CALIFATHETA->cd(1);
    h_CALIFATHETA->Draw("colz");

    TCanvas* c_CALIFAPHI = new TCanvas();
    c_CALIFAPHI->cd(1);
    h_CALIFAPHI->Draw("colz");

    TFile* fout = TFile::Open("CalifaPID.root", "RECREATE");
    h_Califa_Components->Write();
    h_Califa_hits->Write();
    h_Califa_Ns->Write();
    h_Califa_Nf->Write();
    h_Califa_Ns_Nf->Write();
    h_Califa_CC->Write();
    h_CALIFAE->Write();
    h_CALIFATHETA->Write();
    h_CALIFAPHI->Write();
    h_CALIFAECORR->Write();
    h_CALIFATHETASUM->Write();
    gr->Write("gr");
    fout->Close();
}