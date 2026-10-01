void eventDisplay_bea() { FairRunAna *fRun= new FairRunAna();
 
  FairRuntimeDb* rtdb = fRun->GetRuntimeDb();
  FairParRootFileIo*  parIo1 = new FairParRootFileIo();
  parIo1->open("par_16C.root");
  rtdb->setFirstInput(parIo1);
  rtdb->print();
 
  fRun->SetSource(new FairFileSource("sim_16C.root"));
 // fRun->SetOutputFile("test.root");
  fRun->SetSink(new FairRootFileSink("test.root"));
  //FairEventManager *fMan= new FairEventManager();
  R3BEventManager *fMan= new R3BEventManager();
  //FairMCTracksDraw *Track =  new FairMCTracksDraw ("Monte-Carlo Tracks");
  //FairMCPointDraw *LandPoints =   new FairMCPointDraw ("LandPoint",kOrange,  kFullSquare);
 
  //fMan->AddTask(Track);
  //
  //fMan->AddTask(LandPoints);
  //
  auto Track = new R3BMCTracks ("Monte-Carlo Tracks");
  fMan->AddTask(Track);
 
  //gEve->GetDefaultGLViewer()->SetClearColor(kOrange - 4);
  fMan->Init();

  auto fLogger = FairLogger::GetLogger();
  fLogger->SetLogVerbosityLevel("low");
  fLogger->SetLogScreenLevel("info");
  fLogger ->SetColoredLog(true);
 
  gEve->GetDefaultGLViewer()->SetClearColor(kOrange - 4);

 // gEve->GetDefaultGLViewer()->SetClearColor(0);
}
