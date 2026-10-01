/////////////////////////////////////////////////////////////////////////////////////////
// Definitions of strings and histograms 
char fname[200];
char histname[200];
int tpat;

TH1D* h_LOSQ;
TH2D* h_LOSBeam;
TH1D* h_TOFDQ[4];
TH2D* h_TOFDCC;
TH2D* h_CALIFA;
TH1D* h_CALIFAE;
TH1D* h_CALIFAA;
TH2D* h_CALIFATHETA;
TH1D* h_CALIFATHETASUM;
TH2D* h_CALIFAPHI;
TH1D* h_CALIFAA_NOCOPL;
TH2D* h_CALIFAPHI_NOCOPL;

TH2D* h_CALIFAPHIS_BUMP;
TH1D* h_CALIFAPHI_BUMP;
TH2D* h_EVENTS_OPENNING;

vector<double> Eloss1;
vector<double> Eloss2;

vector<double> CalifaR;
vector<double> CalifaL;

void analysis_secondary_peak(){

    /////////////////////////////////////////////////////////////////////////////////////////
    // Set histograms canvas and axis titles here

    h_LOSQ = new TH1D("LOSQ", "LOS Q", 100,0,100);
    h_LOSQ->GetXaxis()->SetTitle("LOS Q");
    h_LOSQ->GetYaxis()->SetTitle("Counts");
    
    h_LOSBeam = new TH2D("LOSBeam", "LOS Beam", 100,-5,5,100,-5,5);
    h_LOSBeam->GetXaxis()->SetTitle("LOS X");
    h_LOSBeam->GetYaxis()->SetTitle("LOS Yy");  

    for(int a=0; a<4; a++) {
        sprintf(histname,"TOFD Plane %i - Q",a+1);
        h_TOFDQ[a] = new TH1D(histname,histname,70,0,7);
        h_TOFDQ[a]->GetXaxis()->SetTitle("Q");
        h_TOFDQ[a]->GetYaxis()->SetTitle("Counts"); 
    }     

    h_TOFDCC = new TH2D("Charge correlation","Charge correlation",100,0,7,100,0,7);
    h_TOFDCC->GetXaxis()->SetTitle("TODF Eloss plane 1");
    h_TOFDCC->GetYaxis()->SetTitle("TOFD Eloss plane 2");  


    h_CALIFAE = new TH1D("Energy","Energy",100,0,300000);
    h_CALIFAE->GetXaxis()->SetTitle("Energy");

    h_CALIFAA = new TH1D("Azimutal angle","Azimutal angle",100,-4,4);
    h_CALIFAA->GetXaxis()->SetTitle("Polar angle");

    h_CALIFATHETA = new TH2D("Angular correlation","Angular correlation",100,0,2,100,0,2);
    h_CALIFATHETA->GetXaxis()->SetTitle("Angle 1");
    h_CALIFATHETA->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAPHI = new TH2D("Azimutal angular correlation","Azimutal angular correlation",100,-3.2,3.2,100,-3.2,3.2);
    h_CALIFAPHI->GetXaxis()->SetTitle("Angle 1");
    h_CALIFAPHI->GetYaxis()->SetTitle("Angle 2");

    h_CALIFATHETASUM = new TH1D("Opening angle","Opening angle",100,0,2);
    h_CALIFATHETASUM->GetXaxis()->SetTitle("Opening angle");

    h_CALIFAA_NOCOPL = new TH1D("Azimutal angle no coplanarity","Azimutal angle no coplanarity",100,0,2);
    h_CALIFAA_NOCOPL->GetXaxis()->SetTitle("Polar angle");

    h_CALIFAPHI_NOCOPL = new TH2D("Azimutal angles","Azimutal angles",100,-3.2,3.2,100,-3.2,3.2);
    h_CALIFAPHI_NOCOPL->GetXaxis()->SetTitle("Azimutal angle 1");
    h_CALIFAPHI_NOCOPL->GetYaxis()->SetTitle("Azimutal angle 2");

    h_CALIFAPHIS_BUMP = new TH2D("Azimutal angles bump","Azimutal angles bump",100,-3.2,3.2,100,-3.2,3.2);
    h_CALIFAPHIS_BUMP->GetXaxis()->SetTitle("Angle 1");
    h_CALIFAPHIS_BUMP->GetYaxis()->SetTitle("Angle 2");

    h_CALIFAPHI_BUMP = new TH1D("Azimutal bump","Azimutal bump",100,-3.2,3.2);
    h_CALIFAPHI_BUMP->GetXaxis()->SetTitle("Azimutal bump");

    h_EVENTS_OPENNING = new TH2D("Califa events per openning","Califa events per openning",10,1,10, 100,0,2);
    h_EVENTS_OPENNING->GetXaxis()->SetTitle("Events");
    h_EVENTS_OPENNING->GetYaxis()->SetTitle("Opennig angle");

    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the root file 

    TChain* chain = new TChain("evt");   
    sprintf(fname, "./main0252_0001_unpacked_califa_ToFD_LOS.root");

    chain->Reset();
    chain->Add(fname); 
    
    /////////////////////////////////////////////////////////////////////////////////////////
    // Load the branches with data as TClonesArrays 

	R3BEventHeader* header = new R3BEventHeader;    
    TBranch* branchData = chain->GetBranch("EventHeader.");
    branchData->SetAddress(&header);    
    
	TClonesArray*  hitLOS = new TClonesArray("R3BLosHitData",5);  // Create array of objects of type R3BLosHitData. What is the use of the second argument ????
	R3BLosHitData** hitLos;                                       // Creates an array of pointers of type R3BLosHitData
	TBranch *branchHitLOS = chain->GetBranch("LosHit");           // Create an object of type TBranch and save on it the branch LosHit of the tree. GetBranch returns a pointer
	branchHitLOS->SetAddress(&hitLOS);                            // Sets address of the clones to the branch  ???? Do not understand this
	Long64_t LosHitPerEvent = 0;                                  // Hits per event in LOS detector

	TClonesArray*  hitTOFD = new TClonesArray("R3BTofdHitData",5);
	R3BTofdHitData** hitTofd;
	TBranch *branchHitTOFD = chain->GetBranch("TofdHit");
	branchHitTOFD->SetAddress(&hitTOFD);
	Long64_t TofdHitPerEvent = 0;  	

    TClonesArray*  hitCALIFA = new TClonesArray("R3BCalifaClusterData",5);
	R3BCalifaClusterData** clusterCalifa;
	TBranch *branchHitCALIFA = chain->GetBranch("CalifaClusterData");
	branchHitCALIFA->SetAddress(&hitCALIFA);
	Long64_t CalifaHitPerEvent = 0;  	

    /////////////////////////////////////////////////////////////////////////////////////////
	
	Long64_t nevent = chain->GetEntries();
	cout << "Number of events in " << fname << " : " << nevent << endl;

	for(Long64_t e=0; e<nevent; e++){   
    
	    hitLOS->Clear();
	    hitTOFD->Clear();
        hitCALIFA->Clear();
	    chain->GetEvent(e);
	    
        // Entries of LOS and TOFD per event in the chain of trees
        tpat            = header->GetTpat();
	    LosHitPerEvent  = hitLOS->GetEntriesFast();
	    TofdHitPerEvent = hitTOFD->GetEntriesFast();
        CalifaHitPerEvent = hitCALIFA->GetEntriesFast();
        // cout << CalifaHitPerEvent << "\n";

        if ((tpat & 8) == 8){
            
            /////////////////////////////////////////////////////////////////////////////////////////

            if(LosHitPerEvent>0){                                     // If there are hits in event 

                hitLos = new R3BLosHitData*[LosHitPerEvent];          // Initialization of hitLos as an array of pointers. Contains all hits in event
                
                for(int p = 0; p < LosHitPerEvent; p++){

                    hitLos[p] =(R3BLosHitData*)hitLOS->At(p);

                    double x_LOS = (hitLos[p]->fX_cm);
                    double y_LOS = (hitLos[p]->fY_cm);
                    double T_LOS = (hitLos[p]->fTime_ns);
                    double z_LOS = (hitLos[p]->fZ);

                    h_LOSBeam->Fill(x_LOS,y_LOS);
                    h_LOSQ->Fill(z_LOS);
                    
                }

                if(LosHitPerEvent) delete[] hitLos;
            }
            
            /////////////////////////////////////////////////////////////////////////////////////////

            bool b_tofdq5 = false;

            if( TofdHitPerEvent > 0 ){

                hitTofd = new R3BTofdHitData*[TofdHitPerEvent];
                
                for(int p = 0; p < TofdHitPerEvent; p++){

                    hitTofd[p] =(R3BTofdHitData*)hitTOFD->At(p);
                    
                    int TOFD_Id = hitTofd[p]->GetDetId();
                    double Eloss = hitTofd[p]->GetEloss();
                    double TOFD_X = hitTofd[p]->GetX();
                    double TOFD_Y = hitTofd[p]->GetY();
                    int bar = hitTofd[p]->GetBarId();
                    
                    ///////////////////////////////////////////
                    // Charge correlation
                    if (TOFD_Id == 1) Eloss1.push_back(Eloss);
                    if (TOFD_Id == 3) Eloss2.push_back(Eloss);
                    ///////////////////////////////////////////

                    if(Eloss<5.5 && Eloss>4.8) {
                        b_tofdq5=true;
                        h_TOFDQ[TOFD_Id - 1]->Fill(Eloss);
                    }
                }

                if(TofdHitPerEvent) delete[] hitTofd;
            }

            /////////////////////////////////////////////////////////////////////////////////////////

            bool califaL=false;
            bool califaR=false;
            double thetaTL, thetaTR, phiTR, phiTL;

            if(b_tofdq5==true){

                if( CalifaHitPerEvent == 2 ){             // The condition >0 displays the secondary bump

                    clusterCalifa = new R3BCalifaClusterData*[CalifaHitPerEvent];

                    double thetasum_plane = 0;
                    double thetasum_plane_no_copl = 0;

                    for(int p = 0; p < CalifaHitPerEvent; p++){

                        clusterCalifa[p] =(R3BCalifaClusterData*)hitCALIFA->At(p);

                        double E = clusterCalifa[p]->GetEnergy();
                        double theta = clusterCalifa[p]->GetTheta();
                        double phi = clusterCalifa[p]->GetPhi();

                        if ( E > 20000 ){

                            /////////////////////////////////////////////////////////////////////////////////////////
                            // NO-COPLANARITY CONDITIONS AND LOGIC
                            thetasum_plane_no_copl += theta;

                            /////////////////////////////////////////////////////////////////////////////////////////
                            // CONDITIONS AND LOGIC

                            h_CALIFAA->Fill(phi);
                            // cout << "Enterd energy conditional" << "\n";

                            if ( (phi < 1.5 && phi > -1.5) && califaL == false){

                                // cout << "Enterd phi left conditional" << "\n";

                                califaL = true;
                                thetaTL = theta;
                                thetasum_plane += thetaTL;
                                phiTL = phi;

                                if (califaR == true && califaL == true){

                                    h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                    h_CALIFAPHI->Fill(phiTL,phiTR);
                                    h_CALIFATHETASUM->Fill(thetasum_plane);

                                    califaL = false;
                                    califaR = false;
                                    thetaTL = 0;
                                    thetaTR = 0;

                                    h_EVENTS_OPENNING->Fill(CalifaHitPerEvent, thetasum_plane);

                                    if ( 0 < thetasum_plane && thetasum_plane < 0.8 ){

                                        h_CALIFAPHIS_BUMP->Fill(phiTL,phiTR);
                                        h_CALIFAPHI_BUMP->Fill(phiTL);
                                        h_CALIFAPHI_BUMP->Fill(phiTR);
                                    }
                                
                                }
                                continue;
                            }

                            if ( (phi>1.5 || phi<-1.5) && califaR == false) { 
                                
                                // cout << "Enterd phi right conditional" << "\n";
                                
                                califaR = true;
                                thetaTR = theta;
                                thetasum_plane += thetaTR;
                                phiTR = phi;

                                if (califaR == true && califaL == true){

                                    h_CALIFATHETA->Fill(thetaTL,thetaTR);
                                    h_CALIFAPHI->Fill(phiTL,phiTR);
                                    h_CALIFATHETASUM->Fill(thetasum_plane);
                                
                                    califaL = false;
                                    califaR = false;
                                    thetaTL = 0;
                                    thetaTR = 0;

                                    h_EVENTS_OPENNING->Fill(CalifaHitPerEvent, thetasum_plane);

                                    if ( 0 < thetasum_plane && thetasum_plane < 0.8 ){

                                        h_CALIFAPHIS_BUMP->Fill(phiTL,phiTR);
                                        h_CALIFAPHI_BUMP->Fill(phiTL);
                                        h_CALIFAPHI_BUMP->Fill(phiTR);
                                    }
                                }

                                continue;

                            }

                            if ((califaR == true && califaL == false) || (califaR == false && califaL == true)){
                                // cout << "Cleaned" << "\n";
                                califaL = false;
                                califaR = false;
                            }
                        }

                        if (thetasum_plane_no_copl != 0){
                            h_CALIFAA_NOCOPL->Fill(thetasum_plane_no_copl);
                        }
                    }
                }
            }
        }  
	}

// cout << Eloss1.size() << "   " << Eloss2.size();
for (int i=1; i<=Eloss1.size(); i++){
    h_TOFDCC->Fill(Eloss1[i],Eloss2[i]);
}
	
// TCanvas* c_LOS = new TCanvas();
// c_LOS->Divide(2,1);
// c_LOS->cd(1);
// h_LOSBeam->Draw("colz");
// c_LOS->cd(2);
// h_LOSQ->Draw("hist");

// TCanvas* c_TOFD = new TCanvas();
// c_TOFD->Divide(2,2);
// c_TOFD->cd(1);
// h_TOFDQ[0]->Draw("hist");
// c_TOFD->cd(2);
// h_TOFDQ[1]->Draw("hist");
// c_TOFD->cd(3);
// h_TOFDQ[2]->Draw("hist");
// c_TOFD->cd(4);
// h_TOFDQ[3]->Draw("hist");

// TCanvas* c_TOFDCC = new TCanvas();
// c_TOFDCC->cd(1);
// h_TOFDCC->Draw("colz");

// TCanvas* c_CALIFAE = new TCanvas();
// c_CALIFAE->cd(1);
// h_CALIFAE->Draw("hist");

TCanvas* c_CALIFATHETASUM = new TCanvas();
c_CALIFATHETASUM->cd(1);
h_CALIFATHETASUM->Draw("hist");

TCanvas* c_CALIFAA = new TCanvas();
c_CALIFAA->cd(1);
h_CALIFAA->Draw("hist");

TCanvas* c_CALIFATHETA = new TCanvas();
c_CALIFATHETA->cd(1);
h_CALIFATHETA->Draw("colz");

TCanvas* c_CALIFAPHI = new TCanvas();
c_CALIFAPHI->cd(1);
h_CALIFAPHI->Draw("colz");

TCanvas* c_CALIFAPHIS_BUMP = new TCanvas();
c_CALIFAPHIS_BUMP->cd(1);
h_CALIFAPHIS_BUMP->Draw("colz");

TCanvas* c_CALIFAPHI_BUMP = new TCanvas();
c_CALIFAPHI_BUMP->cd(1);
h_CALIFAPHI_BUMP->Draw("hist");

TCanvas* c_EVENTS_OPENNING = new TCanvas();
c_EVENTS_OPENNING->cd(1);
//c_EVENTS_OPENNING->SetLogz();
h_EVENTS_OPENNING->Draw("colz");

}