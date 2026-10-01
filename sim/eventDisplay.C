void eventDisplay() {
 // gSystem->Load("libR3BBase");
 // gSystem->Load("libR3BEventDisplay");

  FairRunAna *fRun = new FairRunAna();
  fRun->SetSource(new FairFileSource("sim_alejo.root"));
  fRun->SetSink(new FairRootFileSink("test.root"));

  FairRuntimeDb *rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo *parIo1 = new FairParRootFileIo();
  parIo1->open("par_alejo.root");
  rtdb->setFirstInput(parIo1);
  rtdb->print();

  R3BEventManager *fMan = new R3BEventManager();
  auto Track = new R3BMCTracks ("Monte-Carlo Tracks");
  //R3BMCTracks *Track = new R3BMCTracks("Monte-Carlo Tracks");

  fMan->AddTask(Track);
  fMan->Init();
//  gEve->GetDefaultGLViewer()->SetClearColor(kOrange - 4);
}
