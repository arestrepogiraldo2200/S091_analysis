/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname[200];
char histname[200];

TH1D* h_E;
TH2D* h_THETA;
TH2D* h_ECORR;
TH2D* h_PHI;
TH1D* h_THETASUM;

vector<double> CalifaR;
vector<double> CalifaL;

void sim_analysis_correlations(){

    h_E = new TH1D("Energy","Energy",1000,0,900);
    h_E->GetXaxis()->SetTitle("Energy");

    h_THETA = new TH2D("Angular correlation","Angular correlation",1000,0,2,100,0,2);
    h_THETA->GetXaxis()->SetTitle("Angle 1");
    h_THETA->GetYaxis()->SetTitle("Angle 2");

    h_PHI = new TH2D("Azimutal angular correlation","Azimutal angular correlation",1000,-3.2,3.2,1000,-3.2,3.2);
    h_PHI->GetXaxis()->SetTitle("Angle 1");
    h_PHI->GetYaxis()->SetTitle("Angle 2");

    h_ECORR = new TH2D("Energy correlations","Energy correlations",1000,0,900,1000,0,900);
    h_ECORR->GetXaxis()->SetTitle("Energy 1");
    h_ECORR->GetYaxis()->SetTitle("Energy 2");

    h_THETASUM = new TH1D("Opening angle","Opening angle",1000,0,2);
    h_THETASUM->GetXaxis()->SetTitle("Opening angle");

    /////////////////////////////////////////////////////////////////////////////////////////
    
    // Load the root file 
    TFile *rootfile = TFile::Open("quasi_C16_beam.root");
    //TFile *rootfile = TFile::Open("quasi_H1_beam.root");

    TTree *tree = (TTree*) rootfile->Get("tree");    
    int n = tree->GetEntries();

    // List of trees
    TList *tree_list = new TList();
    tree_list->Add(tree);

    ////////////////////////////////////////////////////////////////////////////////////////////
    // Case of direct leaves in tree. Case of branches and leaves on it tree->SetBranchAddress("myBranch.x", &x);
    Double_t theta1, theta2;
    tree->SetBranchAddress("theta_1", &theta1);
    tree->SetBranchAddress("theta_2", &theta2);

    Double_t phi1, phi2;
    tree->SetBranchAddress("phi_1", &phi1);
    tree->SetBranchAddress("phi_2", &phi2);

    Double_t E1, E2;
    tree->SetBranchAddress("E1", &E1);
    tree->SetBranchAddress("E2", &E2);

    for (int i = 0; i < n; i++) {
        tree->GetEntry(i);
        h_THETA->Fill(theta1, theta2);
        h_THETASUM->Fill(theta1 + theta2);
        h_PHI->Fill(phi1, phi2);
        h_E->Fill(E1 + E2);
        h_ECORR->Fill(E1, E2);

    }
	
TCanvas* c_E = new TCanvas();
c_E->cd(1);
h_E->Draw("hist");

TCanvas* c_THETASUM = new TCanvas();
c_THETASUM->cd(1);
h_THETASUM->Draw("hist");

TCanvas* c_THETA = new TCanvas();
c_THETA->cd(1);
h_THETA->Draw("colz");

TCanvas* c_ECORR = new TCanvas();
c_ECORR->cd(1);
h_ECORR->Draw("colz");

TCanvas* c_PHI = new TCanvas();
c_PHI->cd(1);
h_PHI->Draw("colz");

}
