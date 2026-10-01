
// Components speed
#include <cmath>
#include <TF1.h>

TH2D* h_Califa_Components;
TH1D* h_Califa_hits;
TH1D* h_Califa_Ns;
TH1D* h_Califa_Nf;
TH2D* h_Califa_Ns_Nf;

char fname1[200];

void CalifaPID(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Califa 
    h_Califa_Components = new TH2D("Califa fast and slow components","Califa fast and slow components",500,-50000,900000, 200,0,1);
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
            double Nf = clusterCalifa[p]->GetNf();
            double Ns = clusterCalifa[p]->GetNs();

            h_Califa_Ns->Fill(Ns);
            h_Califa_Nf->Fill(Nf);
            h_Califa_Ns_Nf->Fill(Ns,Nf);

            // Make sure they are not nan
            if (!std::isfinite(Nf) || !std::isfinite(Ns)) {
                continue;
            }

            if ((Nf+Ns)>0){
                // cout << Nf/(Nf+Ns) << endl;
                h_Califa_Components->Fill(Nf+Ns, Nf/(Nf+Ns));
            }
        }
    }

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

    TCanvas* c_Califa_Components = new TCanvas();
    c_Califa_Components->cd(1);
    c_Califa_Components->SetLogz();
    h_Califa_Components->Draw("colz");

    // TF1* hyperbola = new TF1("hyperbola", "1600.0/x + 0.47", 1000, 900000);
    // hyperbola->SetLineColor(kRed);
    // hyperbola->SetLineWidth(2);
    // hyperbola->Draw("same");

    TF1* trendCurve = new TF1("trendCurve", "0.37 + 0.485/pow(x/1000.0 + 5.0, 0.38)", 0, 700000);
    trendCurve->SetLineColor(kRed);
    trendCurve->SetLineWidth(2);
    trendCurve->Draw("same");

}
