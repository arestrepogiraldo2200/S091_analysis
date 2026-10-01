
char histname[200];
char histname_copy[200];
char fname[200];

TH1F* plane_bar[4][44];
TH1F* plane_bar_copy[4][44];
int plane_npeaks[4][44];
int plane_npeaks_copy[4][44];


void ToFD_crystal_calibration(){

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

    for(int a=0; a<4; a++) {
        for(int b=0; b<44; b++) {

            sprintf(histname,"ToFD Plane: %i Bar: %i", a+1, b+1);
            plane_bar[a][b] = new TH1F(histname, histname, 100, 0, 8);
            plane_bar[a][b]->GetXaxis()->SetTitle("Q");
            plane_bar[a][b]->GetYaxis()->SetTitle("Counts"); 

        }
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

                plane_bar[TOFD_Id-1][bar-1]->Fill(Eloss);

            }
        }
    }
    









    
    int cry_1 = 2;
    int cry_2 = 24;

    // 1 - 42 
    TCanvas* c_npeaksone = new TCanvas();
    c_npeaksone->cd(1);
    plane_bar[cry_1][cry_2]->Draw("");

    TCanvas* c_npeaksone_copy = new TCanvas();
    c_npeaksone_copy->cd(1);
    plane_bar_copy[cry_1][cry_2]->Draw("");


}