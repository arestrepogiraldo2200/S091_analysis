// Components speed
#include <cmath>
#include <TF1.h>

TH2D* h_Califa_Components;
TH2D* h_Califa_Components_cut;
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

// When getting a specific crystal
int target_crystal = 4272; 
// int target_crystal = -1; 

char fname1[200], fname2[200], fname3[200], fname4[200], fname5[200], fname6[200], fname7[200], fname8[200], fname9[200], fname10[200], fname11[200], fname12[200], fname13[200];

void CalifaPID(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Califa 
    h_Califa_Components = new TH2D("Califa fast and slow components","Califa fast and slow components",5000,-50000,900000, 200,-0.5,1);
    h_Califa_Components->GetXaxis()->SetTitle("Nf+Ns");
    h_Califa_Components->GetYaxis()->SetTitle("Nf/(Nf+Ns)");

    h_Califa_Components_cut = new TH2D("Califa fast and slow components cut","Califa fast and slow components cut",5000,-50000,900000, 200,-0.5,1);
    h_Califa_Components_cut->GetXaxis()->SetTitle("Nf+Ns");
    h_Califa_Components_cut->GetYaxis()->SetTitle("Nf/(Nf+Ns)");

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

    chain->Reset();
    chain->Add(fname1); 
    // chain->Add(fname2); 
    // chain->Add(fname3); 
    // chain->Add(fname4); 
    // chain->Add(fname5); 
    // chain->Add(fname6); 
    // chain->Add(fname7); 
    // chain->Add(fname8); 
    // chain->Add(fname9); 
    // chain->Add(fname10); 
    // chain->Add(fname11); 
    // chain->Add(fname12); 
    // chain->Add(fname13); 

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

        // ----------------------------------------------------------------------------------------------------
        if (CalifaHitPerEvent == 2){                                         // Condition on two hits in Califa

            bool califaL=false;
            bool califaR=false;
            double thetaTL, thetaTR, phiTR, phiTL, EL, ER;
            double phi1 = 0., phi2 =0., theta1=0., theta2=0.;
            double thetasum_plane = 0;

            for(int p = 0; p < CalifaHitPerEvent; p++){

                clusterCalifa[p] = (R3BCalifaClusterData*)hitCALIFA->At(p);

                double E  = clusterCalifa[p]->GetEnergy();
                double Nf = clusterCalifa[p]->GetNf();
                double Ns = clusterCalifa[p]->GetNs();
                int CC = clusterCalifa[p]->GetMotherCrystal();

                double theta = clusterCalifa[p]->GetTheta();
                double phi   = clusterCalifa[p]->GetPhi();

                if(p==0){
                    phi1 = phi;
                    theta1 = theta; 
                }
                if(p==1){
                    phi2 = phi;
                    theta2 = theta;
                }

                // ----------------------------------------------------------------------------------------------------
                if ( E > 20000 ){

                    double tolerance_theta = 0.2;
                    double tolerance_phi = 0.2;

                    // ----------------------------------------------------------------------------------------------------
                    if ( (phi < 1.5 && phi > -1.5) && califaL == false){

                        califaL = true;
                        thetaTL = theta;
                        thetasum_plane += thetaTL;
                        phiTL = phi;
                        EL = E;

                        if (califaR == true && califaL == true){

                            // ----------------------------------------------------------------------------------------------------
                            // Kinematical conditions
                            if (  abs(phiTL) + abs(phiTR) - 3.1415926535 < tolerance_phi) {
                                if (  abs(thetaTL) + abs(thetaTR) - 1.41372 < tolerance_theta) {

                                    // Search over the crystal
                                    for(int q = 0; q < CalifaCalPerEvent; q++){

                                        caldataCalifa[q] = (R3BCalifaCrystalCalData*)calCALIFA->At(q);
                                        int CCcry = caldataCalifa[q]->GetCrystalId();
                                        double Ecry  = caldataCalifa[q]->GetEnergy();
                                        double Nfcry = caldataCalifa[q]->GetNf();
                                        double Nscry = caldataCalifa[q]->GetNs();

                                        // ----------------------------------------------------------------------------------------------------
                                        if (CCcry == CC && Ecry >= 0.9*E ){                    // Condition on more than 90% of the energy deposited in the mother crystal

                                            h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                            h_CALIFAPHI->Fill(phiTL,phiTR);
                                            h_CALIFATHETASUM->Fill(thetasum_plane);
                                            h_CALIFAECORR->Fill(EL,ER);
                                            h_CALIFAE->Fill(EL+ER);
                                            h_Califa_CC->Fill(CC);

                                            if (CC == target_crystal || target_crystal ==-1){ 

                                                h_Califa_Ns->Fill(Ns);
                                                h_Califa_Nf->Fill(Nf);
                                                h_Califa_Ns_Nf->Fill(Nf,Ns);

                                                // Make sure they are not nan
                                                if (!std::isfinite(Nf) || !std::isfinite(Ns)) {
                                                    continue;
                                                }

                                                if ((Nf+Ns)>0 ){

                                                    h_Califa_Components->Fill(Nf+Ns, Nf/(Nf+Ns));

                                                    if ( Ns<1.4*Nf+3500 && Ns>1.4*Nf-13000 && Ns>0){
                                                        h_Califa_Components_cut->Fill(Nf+Ns, Nf/(Nf+Ns));
                                                    }
                                                }
                                            } 
                                        }
                                    }   
                                }
                            }

                            califaL = false;
                            califaR = false;
                            thetaTL = 0;
                            thetaTR = 0;
                        }
                        continue;
                    }

                    // ----------------------------------------------------------------------------------------------------
                    if ( (phi>1.5 || phi<-1.5) && califaR == false){ 
                                                        
                        califaR = true;
                        thetaTR = theta;
                        thetasum_plane += thetaTR;
                        phiTR = phi;
                        ER = E;

                        if (califaR == true && califaL == true){

                            // ----------------------------------------------------------------------------------------------------
                            // Kinematical conditions
                            if (  abs(phiTL) + abs(phiTR) - 3.1415926535 < tolerance_phi) {
                                if (  abs(thetaTL) + abs(thetaTR) - 1.41372 < tolerance_theta) {

                                    // Search over the crystal
                                    for(int q = 0; q < CalifaCalPerEvent; q++){

                                        caldataCalifa[q] = (R3BCalifaCrystalCalData*)calCALIFA->At(q);
                                        int CCcry = caldataCalifa[q]->GetCrystalId();
                                        double Ecry  = caldataCalifa[q]->GetEnergy();
                                        double Nfcry = caldataCalifa[q]->GetNf();
                                        double Nscry = caldataCalifa[q]->GetNs();

                                        // ----------------------------------------------------------------------------------------------------
                                        if (CCcry == CC && Ecry >= 0.9*E ){                    // Condition on more than 90% of the energy deposited in the mother crystal

                                            h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                            h_CALIFAPHI->Fill(phiTL,phiTR);
                                            h_CALIFATHETASUM->Fill(thetasum_plane);
                                            h_CALIFAECORR->Fill(EL,ER);
                                            h_CALIFAE->Fill(EL+ER);
                                            h_Califa_CC->Fill(CC);

                                            if (CC == target_crystal || target_crystal ==-1){ 

                                                h_Califa_Ns->Fill(Ns);
                                                h_Califa_Nf->Fill(Nf);
                                                h_Califa_Ns_Nf->Fill(Nf,Ns);

                                                // Make sure they are not nan
                                                if (!std::isfinite(Nf) || !std::isfinite(Ns)) {
                                                    continue;
                                                }

                                                if ((Nf+Ns)>0 ){

                                                    h_Califa_Components->Fill(Nf+Ns, Nf/(Nf+Ns));

                                                    if ( Ns<1.4*Nf+3500 && Ns>1.4*Nf-13000 && Ns>0){
                                                        h_Califa_Components_cut->Fill(Nf+Ns, Nf/(Nf+Ns));
                                                    }
                                                }
                                            } 
                                        }
                                    } 
                                }
                            }
                            califaL = false;
                            califaR = false;
                            thetaTL = 0;
                            thetaTR = 0;
                        }
                        continue;
                    }

                    // Case  of two hits  in the same branch of califa
                    if ((califaR == true && califaL == false) || (califaR == false && califaL == true)){
                        califaL = false;
                        califaR = false;
                    }
                }
            }
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    // Fitting of the Nf component

    double xMin = 25000;
    double xMax = 600000;
    double yMin = 0.4;
    double yMax = 0.60;

    TGraph* gr = new TGraph();
    int ip = 0;

    for (int ix = 1; ix <= h_Califa_Components_cut->GetNbinsX(); ++ix) {

        double x = h_Califa_Components_cut->GetXaxis()->GetBinCenter(ix);

        if (x < xMin || x > xMax) continue;

        double sumY = 0.0;
        double sumW = 0.0;

        for (int iy = 1; iy <= h_Califa_Components_cut->GetNbinsY(); ++iy) {

            double y = h_Califa_Components_cut->GetYaxis()->GetBinCenter(iy);

            if (y < yMin || y > yMax) continue;

            double w = h_Califa_Components_cut->GetBinContent(ix, iy);

            if (w <= 0) continue;

            // Takes a summation over the Y components inside the range wheigted by the multiplicities
            sumY += y * w*w;
            sumW += w*w;
        }

        if (sumW > 0) {
            gr->SetPoint(ip++, x, sumY / sumW);
        }
    }

    TF1 *f = new TF1("f",
                 "[0] + [1]/sqrt(x+[2]) + [3]/(x + [4]) + [5]/((x + [6])*(x + [6]))",
                 0, 900000);

    f->SetParLimits(1, 0.0, 1e12);
    f->SetParLimits(3, 0.0, 1e12);
    f->SetParLimits(5, 0.0, 1e12);

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
    // c_Califa_Components->SetLogz();
    h_Califa_Components->Draw("colz");

    f->SetLineWidth(1);
    f->Draw("same");

    TCanvas* c_Califa_Components_cut = new TCanvas();
    c_Califa_Components_cut->cd(1);
    // c_Califa_Components->SetLogz();
    h_Califa_Components_cut->Draw("colz");

    f->SetLineWidth(1);
    f->Draw("same");

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

    // TF1* hyperbola = new TF1("hyperbola", "1600.0/x + 0.47", 1000, 900000);
    // hyperbola->SetLineColor(kRed);
    // hyperbola->SetLineWidth(2);
    // hyperbola->Draw("same");

    // TF1* trendCurve = new TF1("trendCurve", "0.37 + 0.485/pow(x/1000.0 + 5.0, 0.38)", 0, 700000);
    // trendCurve->SetLineColor(kRed);
    // trendCurve->SetLineWidth(2);
    // trendCurve->Draw("same");

    gSystem->mkdir("pdf", kTRUE);
    TCanvas* canvases[] = {c_gr, c_Califa_hits, c_Califa_Ns, c_Califa_Nf, c_Califa_Ns_Nf, c_Califa_CC,
                           c_Califa_Components, c_Califa_Components_cut, c_CALIFAE, c_CALIFAECORR,
                           c_CALIFATHETASUM, c_CALIFATHETA, c_CALIFAPHI};
    const char* names[] = {"gr", "Califa_hits", "Califa_Ns", "Califa_Nf", "Califa_Ns_Nf", "Califa_CC",
                           "Califa_Components", "Califa_Components_cut", "CALIFAE", "CALIFAECORR",
                           "CALIFATHETASUM", "CALIFATHETA", "CALIFAPHI"};
    for (int i = 0; i < 13; i++) {
        canvases[i]->SaveAs(Form("pdf/%s.pdf", names[i]));
    }
}