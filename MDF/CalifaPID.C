
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

// When getting a specific crystal
int target_crystal = -1; 

char fname1[200];

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

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 
    TChain* chain = new TChain("evt"); 
    // sprintf(fname1,"./main0280_0001_unpacked_los.root");
    sprintf(fname1,"./main0252_0001_unpacked_MDF_2.root");

    chain->Reset();
    chain->Add(fname1); 

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 
	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("EventHeader.");
    branchData->SetAddress(&header);  


    TClonesArray* hitCALIFA = new TClonesArray("R3BCalifaClusterData",5);
	R3BCalifaClusterData** clusterCalifa;
	TBranch *branchHitCALIFA = chain->GetBranch("CalifaClusterData");
	branchHitCALIFA->SetAddress(&hitCALIFA);
	Long64_t CalifaHitPerEvent = 0;  

    
	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in the tree: " << nevent << endl;

    /////////////////////////////////////////////////////////////////////////////////////////
	for(Long64_t e=0; e < nevent; e++){

        hitCALIFA->Clear();
        chain->GetEvent(e);

        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        h_Califa_hits->Fill(CalifaHitPerEvent);

        /////////////////////////////////////////////////////////////////////////////////////////
        clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];

        for(int p = 0; p < CalifaHitPerEvent; p++){

            clusterCalifa[p] = (R3BCalifaClusterData*)hitCALIFA->At(p);
            double E  = clusterCalifa[p]->GetEnergy();
            double Nf = clusterCalifa[p]->GetNf();
            double Ns = clusterCalifa[p]->GetNs();
            int CC = clusterCalifa[p]->GetMotherCrystal();

            if (E > 10000){

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

                        if ( Ns<1.4*Nf-4000 && Ns>1.4*Nf-15000 && Ns>0){
                            h_Califa_Components_cut->Fill(Nf+Ns, Nf/(Nf+Ns));
                        }
                    }
                } 
            }    
        }
    }

    /////////////////////////////////////////////////////////////////////////////////////////

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
                 "[0] + [1]/(x + [2]) +[3]/((x + [4])*(x + [4]))",
                 0, 900000);

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

    TF1 *line1 = new TF1("line1", "[0]*x-5000", 0, 100000);
    line1->SetParameter(0,1.4);   
    line1->SetLineWidth(1);
    line1->Draw("same");

    TF1 *line2 = new TF1("line2", "[0]*x-15000", 0, 100000);
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

    f->SetLineWidth(2);
    f->Draw("same");

    TCanvas* c_Califa_Components_cut = new TCanvas();
    c_Califa_Components_cut->cd(1);
    // c_Califa_Components->SetLogz();
    h_Califa_Components_cut->Draw("colz");

    f->SetLineWidth(2);
    f->Draw("same");

    // TF1* hyperbola = new TF1("hyperbola", "1600.0/x + 0.47", 1000, 900000);
    // hyperbola->SetLineColor(kRed);
    // hyperbola->SetLineWidth(2);
    // hyperbola->Draw("same");

    // TF1* trendCurve = new TF1("trendCurve", "0.37 + 0.485/pow(x/1000.0 + 5.0, 0.38)", 0, 700000);
    // trendCurve->SetLineColor(kRed);
    // trendCurve->SetLineWidth(2);
    // trendCurve->Draw("same");
}