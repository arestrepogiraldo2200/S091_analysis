
void p2p(){

   TFile *p2pfile = new TFile("main0252_0001_unpacked_025_final.root");
   TTree *p2ptree = (TTree*)p2pfile->Get("evt");
   //p2ptree->GetListOfBranches()->Print();
   //p2ptree->Draw("TofdHit.fEloss:TofdHit.fBarId>>h(44,1,45,100,0,10)","TofdHit.fDetId==1 && TofdHit.fEloss>4.7 && TofdHit.fEloss<5.7","colz", 50000, 0);

   TTree* p2ptree_charge = p2ptree->CopyTree("TofdHit.fEloss > 4.7");

   //p2ptree_charge->Draw("TofdHit.fEloss","","colz", 50000, 0);
   //Double_t fEloss = ;
}
