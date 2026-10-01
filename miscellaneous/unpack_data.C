/**
 **  Macro to run the offline for all the detectors simultaneously
 **
 **  This macro generates a root file with all the data at mapped level using
 **  a lmd file as input
 **
 **  Author: Jose Luis <j.l.rodriguez.sanchez@udc.es>
 **  @since Feb 5th, 2024
 **
 **
 **  const Int_t nev = -1; number of events to read, -1 - until CTRL+C
 **  Select experiment ID: 91 or 118
 **
 **  After defining the input file (filename), execute the macro
 **  1) if all the parameters are right by default
 **     root -l unpack_data.C
 **  2) if one wants to select a RunId, for instance 'RunId = 173'
 **     root -l 'unpack_data.C(173)'
 **  3) if one wants to select a RunId and max number of events,
 **     for instance 'RunId = 173' and 'nev = 200'
 **     root -l 'unpack_data.C(173,200)'
 **  4) if one wants to select a RunId, max number of events and ExpId,
 **     for instance 'RunId = 173', 'nev = 200', and 'ExpId = 118'
 **     root -l 'unpack_data.C(173,200,118)'
 **
 **/

typedef struct EXT_STR_h101_t {
  EXT_STR_h101_unpack_t unpack;
  EXT_STR_h101_TPAT_t tpat;
  EXT_STR_h101_TRLO_onion_t trloscaler;

  EXT_STR_h101_FRSSCI_onion_t frssci;
  
  EXT_STR_h101_ROLU_onion_t rolu;
  EXT_STR_h101_SOFMWPC_onion_t mwpc;
  EXT_STR_h101_MUSLI_onion_t musli;
  EXT_STR_h101_LOS_onion_t los;

  EXT_STR_h101_FOOT_onion_t foot;
  EXT_STR_h101_CALIFA_t califa;
  EXT_STR_h101_MOSAIC_onion_t mosaic;
  EXT_STR_h101_TTTX_onion_t tttx;

  EXT_STR_h101_FIBEO_onion_t fib30; // Fib 30
  EXT_STR_h101_FIBEI_onion_t fib31; // Fib 31
  EXT_STR_h101_FIBEZ_onion_t fib32; // Fib 32
  EXT_STR_h101_FIBEE_onion_t fib33; // Fib 33
  EXT_STR_h101_TOFD_onion_t tofd;
  EXT_STR_h101_RPC_onion_t rpc;
  EXT_STR_h101_raw_nnp_tamex_onion_t neuland;

  EXT_STR_h101_WRMASTER_t wrmaster;
  EXT_STR_h101_WRNEULAND_t wrneuland;
  EXT_STR_h101_WRFOOT_onion_t wrfoot;
} EXT_STR_h101;

void unpack_data(const Int_t fRunId = 252, const Int_t nev = -1,
                 const Int_t fExpId = 91) {
  TString cRunId = Form("%04d", fRunId);
  TString cExpId = Form("%03d", fExpId);

  FairLogger::GetLogger()->SetLogScreenLevel("info");
  FairLogger::GetLogger()->SetColoredLog(true);

  TStopwatch timer;

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
  TString time_suffix = oss.str();

  // file names and paths
  // -------------------------------------------------------
  TString ntuple_options = "RAW";
//  TString ucesb_dir = getenv("UCESB_DIR");

  TString ucesb_dir = getenv("UCESB_DIR");
  TString filename, outputFilename, upexps_dir, ucesb_path;
  std::vector<UInt_t> WRsFoot;

  if (fExpId == 118 || fExpId == 91) {
    // Input file
    // filename = "/d/land7/202202_para/lmd/main" + cRunId + "_0001.lmd";
//    filename = "/lustre/r3b/202402_s118/lmd_stitched/main" + cRunId + "_*.lmd";
//    TString outputpath = "/lustre/r3b/202402_s118/rootfiles/offline/";
//    outputFilename = outputpath + "main" + cRunId + "_map_data.root";

// main0252_0001.lmd

    filename = "./main" + cRunId + "_*.lmd";
    TString outputpath = "./";
    outputFilename = outputpath + "main" + cRunId + "_map_data.root";

    // outputFilename = "s" + cExpId + "_map_data_" + time_suffix + ".root";

    upexps_dir = ucesb_dir + "/../upexps";
    ucesb_path = upexps_dir + "/202402_s091_s118/202402_s091_s118 --allow-errors --input-buffer=200Mi";
    ucesb_path.ReplaceAll("//", "/");
    WRsFoot.resize(11);
    // FIXME
    WRsFoot[1] = 0x0400;
    WRsFoot[3] = 0x0200;
    WRsFoot[6] = 0x0500;
    WRsFoot[8] = 0x0700;
    WRsFoot[9] = 0x0800;
    WRsFoot[10] = 0x0600;
  } else {
    std::cout << "Experiment was not selected" << std::endl;
    gApplication->Terminate();
  }

  // Setup: Selection of detectors
  UShort_t fNumSci = 2;
  Bool_t fFrsSci = false;   // Start: Plastic scintillators at FRS
  // --- R3B standard -------------------------------------------
  Bool_t fMusli = false;    // Ionization chamber for charge-Z of incoming projectiles
  Bool_t fMwpc2 = false;    // MWPC2 for tracking at entrance of Cave-C
  Bool_t fRolu = false;     // Rolu detector
  Bool_t fLos = false;      // LOS detectors
  Bool_t fFoot = false;     // FOOT tracking detectors
  Bool_t fTttx = false;     // TTTX tracking detectors
  Bool_t fAlpide = false;   // ALIPDE tracking detectors
  Bool_t fCalifa = true;   // Califa calorimeter
  Bool_t fMwpc1 = false;    // MWPC1 for tracking of fragments after the target
  Bool_t fMwpc0 = false;    // MWPC0 for tracking of fragments after the target
  Bool_t fRpc = false;      // RPC detector
  Bool_t fFib30 = false;    // Fiber 30 for fragment tracking behind GLAD
  Bool_t fFib31 = false;    // Fiber 31 for fragment tracking behind GLAD
  Bool_t fFib32 = false;    // Fiber 32 for fragment tracking behind GLAD
  Bool_t fFib33 = false;    // Fiber 33 for fragment tracking behind GLAD
  Bool_t fTofd = false;     // TofD detector
  Bool_t fNeuland = false; // Neuland detector

  // Create source using ucesb for input
  // ----------------------------------------
  EXT_STR_h101 ucesb_struct;

  // Create online run
  // ----------------------------------------------------------
  auto *run = new FairRunOnline();
  auto *EvntHeader = new R3BEventHeader();
  EvntHeader->SetExpId(fExpId);
  run->SetEventHeader(EvntHeader);
  run->SetRunId(fRunId);
  run->SetSink(new FairRootFileSink(outputFilename));

  auto *source = new R3BUcesbSource(filename, ntuple_options, ucesb_path,
                                    &ucesb_struct, sizeof(ucesb_struct));
  source->SetMaxEvents(nev);

  // Add readers
  // ----------------------------------------------------------------
  source->AddReader(new R3BUnpackReader(&ucesb_struct.unpack,
                                        offsetof(EXT_STR_h101, unpack)));

  auto *trloii =
      new R3BTrloiiTpatReader(&ucesb_struct.tpat, offsetof(EXT_STR_h101, tpat));
//  trloii->SetTpatRange(1, 10);
  source->AddReader(trloii);

  auto unpackWRMaster = new R3BWhiterabbitMasterReader(
      (EXT_STR_h101_WRMASTER *)&ucesb_struct.wrmaster,
      offsetof(EXT_STR_h101, wrmaster), 0x1000);
  source->AddReader(unpackWRMaster);

  source->AddReader(new R3BTrloiiScalerReader(
			              (EXT_STR_h101_TRLO_onion*)&ucesb_struct.trloscaler, offsetof(EXT_STR_h101, trloscaler)));

  if (fFrsSci) {
    auto unpackfrssci =
        new R3BFrsSciReader((EXT_STR_h101_FRSSCI *)&ucesb_struct.frssci,
                            offsetof(EXT_STR_h101_t, frssci), fNumSci);
    source->AddReader(unpackfrssci);
  }

  if (fLos) {
    source->AddReader(new R3BLosReader((EXT_STR_h101_LOS *)&ucesb_struct.los,
                                       offsetof(EXT_STR_h101, los)));
  }

  if (fFoot) {
    source->AddReader(
        new R3BFootSiReader((EXT_STR_h101_FOOT_onion *)&ucesb_struct.foot,
                            offsetof(EXT_STR_h101, foot)));
    // source->AddReader(new
    // R3BWhiterabbitFootReader((EXT_STR_h101_WRFOOT_onion*)&ucesb_struct.wrfoot,
    //                                                      offsetof(EXT_STR_h101,
    //                                                      wrfoot), WRsFoot));
  }

  if (fTttx) {
    source->AddReader(
        new R3BTttxReader((EXT_STR_h101_TTTX*)&ucesb_struct.tttx,
                          offsetof(EXT_STR_h101, tttx)));
  }

  if (fCalifa) {
    source->AddReader(
        new R3BCalifaFebexReader((EXT_STR_h101_CALIFA *)&ucesb_struct.califa,
                                 offsetof(EXT_STR_h101, califa)));
  }
  
  if (fAlpide) {
    source->AddReader(
        new R3BMosaicReader((EXT_STR_h101_MOSAIC_onion*)&ucesb_struct.mosaic,
                                                        offsetof(EXT_STR_h101, mosaic)));
  }

  if (fMwpc0 || fMwpc1 || fMwpc2) {
    auto *mwpcreader =
        new R3BMwpcReader((EXT_STR_h101_SOFMWPC *)&ucesb_struct.mwpc,
                          offsetof(EXT_STR_h101, mwpc));
    mwpcreader->SetMaxNbDet(3);
    source->AddReader(mwpcreader);
  }

  if (fMusli) {
    source->AddReader(
        new R3BMusliReader((EXT_STR_h101_MUSLI*)&ucesb_struct.musli,
                           offsetof(EXT_STR_h101, musli)));
  }

  if (fFib30) {
    source->AddReader(new R3BFiberReader(
        "Fi30", 512, (EXT_STR_h101_FIBEO_onion *)&ucesb_struct.fib30,
        offsetof(EXT_STR_h101_t, fib30)));
  }

  if (fFib31) {
    source->AddReader(new R3BFiberReader(
        "Fi31", 512, (EXT_STR_h101_FIBEI_onion *)&ucesb_struct.fib31,
        offsetof(EXT_STR_h101_t, fib31)));
  }

  if (fFib32) {
    source->AddReader(new R3BFiberReader(
        "Fi32", 512, (EXT_STR_h101_FIBEZ_onion *)&ucesb_struct.fib32,
        offsetof(EXT_STR_h101_t, fib32)));
  }

  if (fFib33) {
    source->AddReader(new R3BFiberReader(
        "Fi33", 512, (EXT_STR_h101_FIBEE_onion *)&ucesb_struct.fib33,
        offsetof(EXT_STR_h101_t, fib33)));
  }

  if (fTofd) {
    source->AddReader(
        new R3BTofdReader((EXT_STR_h101_TOFD_onion *)&ucesb_struct.tofd,
                          offsetof(EXT_STR_h101, tofd)));
  }

  if (fRpc) {
    source->AddReader(
        new R3BRpcReader((EXT_STR_h101_RPC*)&ucesb_struct.rpc,
                         offsetof(EXT_STR_h101, rpc)));
  }

  if (fNeuland) {
    source->AddReader(new R3BNeulandTamexReader(
        (EXT_STR_h101_raw_nnp_tamex_onion*)&ucesb_struct.neuland,
        offsetof(EXT_STR_h101, neuland)));
    // source->AddReader(new
    // R3BWhiterabbitNeulandReader((EXT_STR_h101_WRNEULAND*)&ucesb_struct.wrneuland,
    //                                                               offsetof(EXT_STR_h101,
    //                                                               wrneuland),
    //                                                               0x900));
  }

  run->SetSource(source);

  // Initialize
  // -----------------------------------------------------------------
  timer.Start();
  run->Init();

  // Run
  // ------------------------------------------------------------------------
  run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);

  // Finish
  // ---------------------------------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime() / 60.;
  Double_t ctime = timer.CpuTime() / 60.;
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl;
  std::cout << "Output file is " << outputFilename << std::endl;
  std::cout << "Real time " << rtime << " min, CPU time " << ctime << " min"
            << std::endl
            << std::endl;
  gApplication->Terminate();
}
