/*
 *  Macro to run the online for all the detectors simultaneously
 *
 *  Author: Jose Luis <j.l.rodriguez.sanchez@udc.es>
 *  @since Jan 30th, 2024
 *
 */
#include <iostream>
#include <fstream>
#include <sstream>

typedef struct EXT_STR_h101_t {
  EXT_STR_h101_unpack_t unpack;
  EXT_STR_h101_TPAT_t unpacktpat;

  EXT_STR_h101_FRSSCI_onion_t frssci;
  EXT_STR_h101_SOFMWPC_onion_t mwpc;
  EXT_STR_h101_MUSLI_onion_t musli;
  EXT_STR_h101_ROLU_onion_t rolu;
  EXT_STR_h101_LOS_t los;
  EXT_STR_h101_MOSAIC202402_onion_t mosaic;
  EXT_STR_h101_FOOT_onion_t foot;
  EXT_STR_h101_CALIFA202402_onion_t califa;
  EXT_STR_h101_TTTX_onion_t tttx;
  EXT_STR_h101_RPC rpc;
  EXT_STR_h101_TOFD_onion_t tofd;
  EXT_STR_h101_raw_nnp_tamex_onion_t raw_nnp;
  
  EXT_STR_h101_FIBEO_onion_t fib30; // Fib 30
  EXT_STR_h101_FIBEI_onion_t fib31; // Fib 31
  EXT_STR_h101_FIBEZ_onion_t fib32; // Fib 32
  EXT_STR_h101_FIBEE_onion_t fib33; // Fib 33

  EXT_STR_h101_WRMASTER_t wrmaster;
  EXT_STR_h101_WRLOS_t wrlos;
  EXT_STR_h101_WRCALIFA_t wrcalifa;
} EXT_STR_h101;

void unpack_runs_s091_MDF_updated(const TString fRun = "12C_2", const Int_t fRunId = 999, const Int_t runnum = 274 ,const Int_t nev = -1, Int_t nb = 1) {
  TStopwatch timer;
  timer.Start();

  FairLogger::GetLogger()->SetLogScreenLevel("info");
  FairLogger::GetLogger()->SetColoredLog(true);

  auto t = std::time(nullptr);
  auto tm = *std::localtime(&t);
  std::ostringstream oss;
  oss << std::put_time(&tm, "%Y%m%d_%H%M%S");
  TString time_suffix = oss.str();

  UShort_t fNumSci = 2;
  const Int_t fExpId = 118;
  
	// OLD MDF Settings
	const Double_t Angle = -18.0 * TMath::Pi() / 180.;//Central turning angle of GLAD
	const Double_t Target_to_GLAD_flange = 127.9;

	const Double_t Target_to_MWPC1  = Target_to_GLAD_flange - 27.5 ;
	const Double_t Target_to_MWPC0  = Target_to_GLAD_flange - 9.5 ;

	const Double_t TP_to_flange  = 385.0 - 125.0;
	const Double_t Target_to_TP = Target_to_GLAD_flange + 165.;

	const Double_t TP_to_F32 = TP_to_flange + 125.0; 
	const Double_t TP_to_F30 = TP_to_flange + 153.5; 
	const Double_t TP_to_F33 = TP_to_flange + 225.0; 
	const Double_t TP_to_F31 = TP_to_flange + 242.5; 

	const Double_t F31_offset = -26.9;
	const Double_t F33_offset = 24.78;


  	// MDF Settings
	 Double_t GLAD_CURRENT;
         //Double_t GLAD_CURRENT = 1490.0; //currnet for MDF 
	 //Double_t GLAD_CURRENT = 1150.0; //current for 10C 1st
	 //Double_t GLAD_CURRENT = 1320.0; //current for 10C 2nd
	 //Double_t GLAD_CURRENT = 1456.0; //current for 12C
	 //Double_t GLAD_CURRENT = 1680.0; //current for 14C 1679.8
	 //Double_t GLAD_CURRENT = 1960.0; //current for 16C
	 const Double_t GLAD_CURRENT_REF = 1456.0; //trained current
	 //Double_t GLAD_CURRENT_REF = 1695.0; //trained current

	/*
	const Double_t Angle = -18.0 * TMath::Pi() / 180.;//Central turning angle of GLAD
	const Double_t Target_to_GLAD_flange = 128.6; //127.9

	const Double_t Target_to_MWPC1  = Target_to_GLAD_flange - 27.5 ;
	const Double_t Target_to_MWPC0  = Target_to_GLAD_flange - 9.5 ;

	const Double_t TP_to_flange  = 385.0 - 125.0;
	const Double_t Target_to_TP = Target_to_GLAD_flange + 165.;

	//const Double_t TP_to_F32 = TP_to_flange + 125.0; 
	//const Double_t TP_to_F30 = TP_to_flange + 153.5; 
	//const Double_t TP_to_F33 = TP_to_flange + 225.5; //225.0
	//const Double_t F33_to_F31 = 16.5;
	//const Double_t TP_to_F31 = TP_to_flange + 242.5; 
	//const Double_t TP_to_F31 = TP_to_flange + TP_to_F33 + F33_to_F31; 
	const Double_t TP_to_F32 = TP_to_flange + 125.6; 
	const Double_t TP_to_F30 = TP_to_flange + 154.1; 
	const Double_t TP_to_F33 = TP_to_flange + 225.5; //225.0
	const Double_t TP_to_F31 = TP_to_flange + 242.0; 
	//const Double_t F33_to_F31 = 16.5;
	//const Double_t TP_to_F31 = TP_to_flange + TP_to_F33 + F33_to_F31; 

	const Double_t F31_offset = -25.6;
	const Double_t F33_offset = 25.0;
*/
	const Double_t ToF_offset = 0.8;//ns added to ToF in TOFD, gives correct 124Sn mass for run 497

  //// Neuland
  const Int_t nBarsPerPlane = 50; // number of scintillator bars per plane
  const Int_t nPlanes = 26;       // number of planes (for TCAL calibration)
  const double distanceToTarget = 1564.15; // 1520 //rough estimation of |proton start-NeuLAND| distance


  // File names and paths -----------------------------
  TString dir = gSystem->Getenv("VMCWORKDIR");
  // TString ntuple_options = "RAW,time-stitch=2000";
  TString ntuple_options = "RAW";
  //TString ucesb_dir = getenv("UCESB_DIR");
  //TString ucesb_dir = "/lustre/land/mwhiteh/2020402/unpacker/";
  TString ucesb_dir = "/u/arestrep/2020402/unpacker/";
  TString filename, outputFilename, upexps_dir, ucesb_path, pardir, mydir;
  TString calibparfile;
  filename = TString::Format("/lustre/r3b/202402_s118/lmd_stitched/main%04d_%04d_stitched.lmd",runnum,nb);
  //filename = "/lustre/land/mwhiteh/s118_files/main0126_stitched.lmd";
  //outputFilename = "/lustre/land/mwhiteh/s118_files/main0126_out.root";


    if (fRun == "10C_1"){
	    cout << "************ 10C(1) Runs Selected *************" << endl;
	    calibparfile = "parameters/CalibPar_10C1.par";
	    GLAD_CURRENT = 1150.0;
//  	    outputFilename = TString::Format("/lustre/land/mwhiteh/unpacked_files/10C/main%04d_%04d_unpacked_tpats.root",runnum,nb);
    }
    if (fRun == "10C_2"){
	    cout << "************ 10C(2) Runs Selected *************" << endl;
	    calibparfile = "parameters/CalibPar_10C2.par";
	    GLAD_CURRENT = 1320.0;
//  	    outputFilename = TString::Format("/lustre/land/mwhiteh/unpacked_files/10C/main%04d_%04d_unpacked_tpats.root",runnum,nb);
    }
    if (fRun == "12C_1"){
	    cout << "************ 12C Runs Selected *************" << endl;
	    calibparfile = "parameters/CalibPar_12C1.par";
	    GLAD_CURRENT = 1456.0;
  	    outputFilename = TString::Format("/u/land/gmw537/main%04d_%04d_unpacked_los.root",runnum,nb);
    }
    if (fRun == "12C_2"){
	    cout << "************ 12C Runs Selected *************" << endl;
	    calibparfile = "parameters/CalibPar_12C2.par";
	    GLAD_CURRENT = 1456.0;
            outputFilename = TString::Format("/u/arestrep/main%04d_%04d_unpacked_los.root",runnum,nb);
//  	    outputFilename = TString::Format("/lustre/land/mwhiteh/unpacked_files/12C/main%04d_%04d_unpacked_tpats.root",runnum,nb);
    }
    if (fRun == "14C"){
	    cout << "************ 14C Runs Selected *************" << endl;
	    calibparfile = "parameters/CalibPar_14C.par";
	    GLAD_CURRENT = 1680.0;
//  	    outputFilename = TString::Format("/lustre/land/mwhiteh/unpacked_files/14C/main%04d_%04d_unpacked_tpats.root",runnum,nb);
    }
    if (fRun == "16C"){
	    cout << "************ 16C Runs Selected *************" << endl;
	    calibparfile = "parameters/CalibPar_16C.par";
	    GLAD_CURRENT = 1960.0;
            outputFilename = TString::Format("/u/land/gmw537/main%04d_%04d_unpacked_califa_ToFD_LOS_NEULAND.root",runnum,nb);
//  	    outputFilename = TString::Format("/lustre/land/mwhiteh/unpacked_files/16C/main%04d_%04d_unpacked_tpats.root",runnum,nb);
    }
    if (fRun == "s118"){
	    cout << "************ s118 Runs Selected *************" << endl;
	    calibparfile = "parameters/CalibPar_s118.par";
//  	    outputFilename = TString::Format("/lustre/land/mwhiteh/unpacked_files/s118/main%04d_%04d_unpacked_03.root",runnum,nb);
            outputFilename = TString::Format("/u/land/gmw537/main%04d_%04d_unpacked_MDF.root",runnum,nb);
	    // GLAD_CURRENT = 980.0;
            GLAD_CURRENT = 2193.0; // Current in the elog
    }
    cout << "GLAD CURRENT: " << GLAD_CURRENT << endl;

/*
  
  #define PREFIX "/lustre/r3b/202402_s118/lmd_stitched/main"
  #define POSTFIX "_*_stitched.lmd "
  filename="";
    auto addRun=[&filename](const char* name) 
    {
    filename+=PREFIX;
    filename+=name;
    filename+=POSTFIX;
    };
  
  if (fExpId == 118) {
    if (fRun == "2H_200"){
	std::cout << "200MeV deuteron runs" << std::endl;
    	addRun("0126");
    	addRun("0127");
    	addRun("0128");
    }
    
    if (fRun == "12C_1")
    {
	std::cout << "First 12C runs" << std::endl;
	calibparfile = "parameters/CalibPar_12C1.par";
    	addRun("0173");
    	addRun("0174");
    	addRun("0175");
    	addRun("0176");
    	addRun("0177");
    	addRun("0178");
    	addRun("0179");
    	addRun("0180");
    	addRun("0181");
    	addRun("0182");
    	addRun("0183");
    	addRun("0184");
    	addRun("0185");
    	addRun("0186");
    	addRun("0187");
    	addRun("0188");
    	addRun("0189");
    	addRun("0190");
    	addRun("0191");
    	addRun("0192");
    	addRun("0193");
    	addRun("0194");
    }

    if (fRun == "10C_1")
    {
	std::cout << "First 10C runs" << std::endl;
	calibparfile = "parameters/CalibPar_10C1.par";
    	addRun("0197");
    	addRun("0198");
    	addRun("0199");
    	addRun("0200");
    	addRun("0201");
    	addRun("0202");
    	addRun("0203");
    	addRun("0204");
    	addRun("0205");
    	addRun("0206");
    	addRun("0207");
    	addRun("0208");
    	addRun("0209");
    	addRun("0210");
    	addRun("0211");
    	addRun("0212");
    	addRun("0213");
    	addRun("0214");
    	addRun("0215");
    	addRun("0216");
    	addRun("0217");
    	addRun("0218");
    	addRun("0219");
    	addRun("0220");
    	addRun("0221");
    	addRun("0222");
    	addRun("0223");
    }
    
    if (fRun == "16C")
    {
	    
	calibparfile = "parameters/CalibPar_16C.par";
	std::cout << "16C runs" << std::endl;
    	//addRun("0231"); //can be included?
    	//addRun("0232"); //can be included?
    	addRun("0233");
    	addRun("0234");
    	addRun("0235");
    	addRun("0236");
    	addRun("0237");
    	addRun("0238");
    	addRun("0239");
    	addRun("0240");
    	addRun("0241");
    	addRun("0242");
    	addRun("0243");
    	addRun("0244");
    	addRun("0245");
    	addRun("0246");
    	addRun("0247");
    	addRun("0248");
    	addRun("0249");
    	addRun("0250");
    	addRun("0251");
    	addRun("0252");
    }
    
    if (fRun == "14C")
    {
	calibparfile = "parameters/CalibPar_14C.par";
	std::cout << "14C runs" << std::endl;
    	addRun("0253");
    	addRun("0254");
    	addRun("0255");
    	addRun("0256");
    	addRun("0257");
    	addRun("0258");
    	addRun("0259");
    	addRun("0260");
    	addRun("0261");
    	addRun("0262");
    	addRun("0263");
    	addRun("0264");
    	addRun("0265");
    	addRun("0266");
    	addRun("0267");
    	addRun("0268");
    	addRun("0269");
    	addRun("0270");
    	addRun("0271");
    }

    if (fRun == "12C_2")
    {
	calibparfile = "parameters/CalibPar_12C2.par";
	std::cout << "Second 12C runs" << std::endl;
    	addRun("0274");
    	addRun("0275");
    	addRun("0276");
    	addRun("0277");
    	addRun("0278");
    	addRun("0279");
    	addRun("0280");
    	addRun("0281");
    	addRun("0282");
    	addRun("0283");
    	addRun("0284");
    	addRun("0285");
    }
 
    if (fRun == "10C_2")
    {
	calibparfile = "parameters/CalibPar_10C2.par";
	std::cout << "Second 10C runs" << std::endl;
    	addRun("0290");
    	addRun("0291");
    	addRun("0292");
    	addRun("0293");
    	addRun("0294");
    	addRun("0295");
    	addRun("0296");
    	addRun("0297");
    	addRun("0298");
    	addRun("0299");
    	addRun("0300");
    	addRun("0301");
    	addRun("0302");
    	addRun("0303");
    	addRun("0304");
    	addRun("0305");
    	addRun("0306");
    	addRun("0307");
    	addRun("0308");
    	addRun("0309");
    	addRun("0310");
    	addRun("0311");
    	addRun("0312");
    	addRun("0313");
    	addRun("0314");
    	addRun("0315");
    }
  */
  /*
  #define PREFIX "/lustre/r3b/202402_s118/lmd_stitched/main"
  #define POSTFIX "_*_stitched.lmd "
  filename="";
    auto addRun=[&filename](const char* name) 
    {
    filename+=PREFIX;
    filename+=name;
    filename+=POSTFIX;
    };
  
  if (fExpId == 118) {
    //filename = "/lustre/r3b/202402_s118/lmd_stitched/main0126_*_stitched.lmd";
    if (fRun == "2H_200MeV"){
	std::cout << "200MeV deuteron runs" << std::endl;
    	addRun("0126");
    	addRun("0127");
    	addRun("0128");
    }

    if (fRun == "10C"){

	std::cout << "10C runs" << std::endl;
    	addRun("0197");
    	addRun("0198");
    	addRun("0199");
    	addRun("0200");
    	addRun("0201");
    	addRun("0202");
    	addRun("0203");
    	addRun("0204");
    	addRun("0205");
    	addRun("0206");
    	addRun("0207");
    	addRun("0208");
    	addRun("0209");
    	addRun("0210");
    	addRun("0211");
    	addRun("0212");
    	addRun("0213");
    	addRun("0214");
    	addRun("0215");
    	addRun("0216");
    	addRun("0217");
    	addRun("0218");
    	addRun("0219");
    	addRun("0220");
    	addRun("0221");
    	addRun("0222");
    	addRun("0223");
    }

    if (fRun == "12C_2"){
	std::cout << "Second 12C runs" << std::endl;
    	addRun("0274");
    	addRun("0275");
    	addRun("0276");
    	addRun("0277");
    	addRun("0278");
    	addRun("0279");
EXT_STR_h101 ucesb_struct;

    // Create online run ------------------------------------
    auto run = new FairRunOnline();
    auto EvntHeader = new R3BEventHeader();
    run->SetEventHeader(EvntHeader);
    EvntHeader->SetExpId(expId);
    run->SetRunId(5);
    run->SetSink(new FairRootFileSink(outputFileName));
    run->ActivateHttpServer(refresh, online_port);
    	addRun("0280");
    	addRun("0281");
    	addRun("0282");
    	addRun("0283");
    	addRun("0284");
    	addRun("0285");
    }
    */
   /* 
    TString outputpath = "/u/mwhiteh/s091/MDF/MDF_S091/rootfiles/";
    //outputFilename = outputpath + "online_tracker_" + time_suffix + ".root";
    outputFilename = outputpath + "mdf_12C_test_" + time_suffix + ".root";
    */
    
    upexps_dir = ucesb_dir;
    ucesb_path = upexps_dir + "/202402_s091_s118 --allow-errors --input-buffer=600Mi";
    ucesb_path.ReplaceAll("//", "/");
    pardir =  "./R3BParams_s091_s118/";
    //mydir = "/u/mwhiteh/s091/R3BParams_s091_s118/";
    //mydir = "/lustre/land/mwhiteh/s091/R3BParams_s091_s118/";
    mydir = "./R3BParams_s091_s118_Desa/";

 //} else {
  //  std::cout << "Experiment was not selected" << std::endl;
  //  gApplication->Terminate();
  //}
  

    // Create online run ------------------------------------

  // store data or not ------------------------------------
  Bool_t NOTstoremappeddata = false; // if true, don't store mapped data in the root file
  Bool_t NOTstorecaldata = false; // if true, don't store cal data in the root file
  Bool_t NOTstorehitdata = false; // if true, don't store hit data in the root file

  // Setup: Selection of detectors ------------------------
  // --- R3B standard -------------------------------------
  Bool_t fFrsSci = true;   // Start: Plastic scintillators at S2
  Bool_t fMwpc2 = true;    // MWPC2 for tracking at entrance of Cave-C
  Bool_t fMusli = true;    // Ionization chamber for charge-Z of incoming projectiles
  Bool_t fRolu = true;    // Rolu detector
  Bool_t fLos = true;     // LOS detectors
  Bool_t fAlpide = false;  // Alpide pixel detectors
  Bool_t fAlpidecal = false;  // Alpide pixel detectors
  Bool_t fFoot = false;    // FOOT tracking detectors
  Bool_t fCalifa = true;  // Califa calorimeter
  Bool_t fTttx = true;    // TTT10 silicon for charge after target
  Bool_t fMwpc1 = true;    // MWPC1 for tracking of fragments after the target
  Bool_t fMwpc0 = true;    // MWPC0 for tracking of fragments after the target
  Bool_t fRpc =   true;     // RPC detector
  Bool_t fFib30 = true;   // Fiber 30 for fragment tracking behind GLAD
  Bool_t fFib31 = true;   // Fiber 31 for fragment tracking behind GLAD
  Bool_t fFib32 = true;   // Fiber 32 for fragment tracking behind GLAD
  Bool_t fFib33 = true;   // Fiber 33 for fragment tracking behind GLAD
  Bool_t fTofd = true;    // Tofd
  Bool_t fNeuland = false; // Neuland
  
  
  // --- Tracking   ------------------------------------------------------------
  Bool_t fTrackingFrs = true; // Tracking of incoming projectiles
  Bool_t fTracking = true; // Tracking of fragments before GLAD
  Bool_t fMDF = false; //MDF tracking of fragments
  if (fMDF){
	fMwpc2 = true;
	fMwpc1 = true;
	fMwpc0 = true;
	fFib30 = true;
	fFib31 = true;
	fFib32 = true;
	fFib33 = true;
	fTofd  = true;
  }

  // Calibration files for SOFIA ----------------------------------------------
  //TString calfilename = pardir + "parameters/CalibPar.par";
  TString calfilename = mydir + calibparfile;
  calfilename.ReplaceAll("//", "/");
  
  // Parameters for LOS
  TString loscalfilename = mydir + "los/loscalpar_v2.root";
  loscalfilename.ReplaceAll("//", "/");
  // Parameters for FOOT   ----------------------------------------------------
 
 TString footparfilename = "./R3BParams_s091_s118_Desa/foot/foot_mapping_v1.par";
 TString footcalfilename = "footpar_in.par";
 TString tempFootParFile = "./R3BParams_s091_s118_Desa/foot/footpar_v5_temp.par";
 footparfilename.ReplaceAll("//", "/");


  // Parameters for CALIFA   ----------------------------------------------------
  TString califamapfilename = pardir + "califa/califamapping_v1.par";

  califamapfilename.ReplaceAll("//", "/");
  
  //TString califacalfilename = pardir + "califa/califacalpar_v3.par";
  TString califacalfilename = "./R3BParams_s091_s118/califa/califacalpar_v3.par";
  std::cout<<"./R3BParams_s091_s118/califa/califacalpar_v3.par"<<std::endl;
  if(!califacalfilename){std::cout<<"AIDSAIDAIDIADISisda"<<std::endl;}
  califacalfilename.ReplaceAll("//", "/");

  // Create source using ucesb for input ------------------
  EXT_STR_h101 ucesb_struct;

  // Create online run ------------------------------------
  auto *run = new FairRunOnline();
  auto *EvntHeader = new R3BEventHeader();
  EvntHeader->SetExpId(fExpId);
  run->SetEventHeader(EvntHeader);
  run->SetRunId(fRunId);
  run->SetSink(new FairRootFileSink(outputFilename));
  //run->ActivateHttpServer(refresh, online_port);

  auto *source = new R3BUcesbSource(filename, ntuple_options, ucesb_path, &ucesb_struct,
                                                                            sizeof(ucesb_struct));
  source->SetMaxEvents(nev);

  // Definition of reader ---------------------------------
  // Add readers ------------------------------------------
  source->AddReader(new R3BUnpackReader(&ucesb_struct.unpack,
                                        offsetof(EXT_STR_h101, unpack)));
  auto unpacktrloii = new R3BTrloiiTpatReader(&ucesb_struct.unpacktpat, offsetof(EXT_STR_h101, unpacktpat));
  //unpacktrloii->SetTrigger(1);
  //unpacktrloii->SetTpatRange(1,12);
  source->AddReader(unpacktrloii);

  source->AddReader(new R3BWhiterabbitMasterReader((EXT_STR_h101_WRMASTER*)&ucesb_struct.wrmaster,
                                                                  offsetof(EXT_STR_h101, wrmaster), 0x1000));
  
  if (fFrsSci)
  {        
    auto unpackfrssci = new R3BFrsSciReader((EXT_STR_h101_FRSSCI*)&ucesb_struct.frssci, offsetof(EXT_STR_h101_t, frssci),fNumSci);
    unpackfrssci->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackfrssci);
  }

  if (fMwpc0 || fMwpc1 || fMwpc2) 
  {
    auto unpackmwpc = new R3BMwpcReader((EXT_STR_h101_SOFMWPC*)&ucesb_struct.mwpc,
                                                      offsetof(EXT_STR_h101, mwpc));
    unpackmwpc->SetMaxNbDet(3);
    unpackmwpc->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackmwpc);
  }

  if (fMusli) 
  {
    auto unpackmusli = new R3BMusliReader((EXT_STR_h101_MUSLI_t *)&ucesb_struct.musli,
                                                       offsetof(EXT_STR_h101, musli));
    unpackmusli->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackmusli);
  }
  
//if(fRolu)
//  {
//    auto unpackrolu = new R3BRoluReader((EXT_STR_h101_ROLU_onion*)&ucesb_struct.rolu, offsetof(EXT_STR_h101, rolu));
//    unpackrolu->SetNbDet(1);
//    unpackrolu->SetOnline(NOTstoremappeddata);
//    source->AddReader(unpackrolu);
//  }

  if (fLos) 
  {
    auto unpackWRlos = new R3BWhiterabbitLosReader((EXT_STR_h101_WRLOS*)&ucesb_struct.wrlos,
                                                      offsetof(EXT_STR_h101, wrlos), 0x1100); // FIXME
    unpackWRlos->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackWRlos);

    auto unpacklos = new R3BLosReader((EXT_STR_h101_LOS*)&ucesb_struct.los,
                                                offsetof(EXT_STR_h101, los));
    unpacklos->SetOnline(NOTstoremappeddata);
    source->AddReader(unpacklos);
  }

  if (fAlpide)
  {
	  auto unpackalpide = new R3BMosaicReader((EXT_STR_h101_MOSAIC202402_onion*)&ucesb_struct.mosaic, offsetof(EXT_STR_h101, mosaic));
    std::vector<int> mosaic_map = { 0, 0, 0, 0, 0, 1, 2, 3, 0};
    unpackalpide->SetMosaicMapping(mosaic_map);
	  unpackalpide->SetOnline(false);
	  source->AddReader(unpackalpide);
  }

  if (fFoot){
    auto unpackfoot = new R3BFootSiReader((EXT_STR_h101_FOOT_onion*)&ucesb_struct.foot, offsetof(EXT_STR_h101, foot));
    // std::vector<int> footIds = { 2, 13, 3, 11, 7, 6, 9, 12, 1, 3, 5, 8, 10, 14, 15, 16};
    std::vector<int> footIds = { 2, 13, 4, 11, 6, 7, 12, 9, 1, 3, 5, 8, 10, 14, 15, 16};
    unpackfoot->SetMapping(footIds);

    unpackfoot->SetOnline(false);
    source->AddReader(unpackfoot);
  }

  if (fCalifa)
  {
	 auto unpackWRC = new R3BWhiterabbitCalifaReader((EXT_STR_h101_WRCALIFA*)&ucesb_struct.wrcalifa,
			 offsetof(EXT_STR_h101, wrcalifa), 0x1a00, 0x1b00, 0x1c00, 0x1d00);
	 unpackWRC->SetOnline(false);
	 source->AddReader(unpackWRC);
	 
	 
	  auto unpackcalifa = new R3BCalifaFebexReader((EXT_STR_h101_CALIFA202402_onion*)&ucesb_struct.califa,
			  offsetof(EXT_STR_h101, califa));
	  //unpackcalifa->SetOnline(NOTstoremappeddata);
	  unpackcalifa->SetOnline(false);
	  source->AddReader(unpackcalifa);
  }

  if (fTttx){
	R3BTttxReader *ttt10_reader = new R3BTttxReader((EXT_STR_h101_TTTX_t*)&ucesb_struct.tttx,
		       	offsetof(EXT_STR_h101, tttx));
	ttt10_reader->SetOnline(NOTstoremappeddata);
	source->AddReader(ttt10_reader);	
  }

//  if (fRpc) 
//  { 
    auto unpackRPC = new R3BRpcReader((EXT_STR_h101_RPC*)&ucesb_struct.rpc, offsetof(EXT_STR_h101, rpc));
    unpackRPC->SetOnline(NOTstoremappeddata);
    source->AddReader(unpackRPC);
//  }
  
  if (fFib30) 
  {
    auto fiber30reader =new R3BFiberReader("Fi30",512,(EXT_STR_h101_FIBEO_onion*)&ucesb_struct.fib30,
                                                                      offsetof(EXT_STR_h101_t, fib30));
    fiber30reader->SetOnline(NOTstoremappeddata);
    source->AddReader(fiber30reader);
  }
  
  if (fFib31) 
  {
    auto fiber31reader = new R3BFiberReader("Fi31",512,(EXT_STR_h101_FIBEI_onion*)&ucesb_struct.fib31,
                                                                      offsetof(EXT_STR_h101_t, fib31));
    fiber31reader->SetOnline(NOTstoremappeddata);
    source->AddReader(fiber31reader);
  }
  
  if (fFib32) 
  {
    auto fiber32reader = new R3BFiberReader("Fi32",512,(EXT_STR_h101_FIBEZ_onion*)&ucesb_struct.fib32,
                                                                      offsetof(EXT_STR_h101_t, fib32));
    fiber32reader->SetOnline(NOTstoremappeddata);
    source->AddReader(fiber32reader);
  }
  
  if (fFib33) 
  {
    auto fiber33reader = new R3BFiberReader("Fi33",512,(EXT_STR_h101_FIBEE_onion*)&ucesb_struct.fib33,
                                                                      offsetof(EXT_STR_h101_t, fib33));
    fiber33reader->SetOnline(NOTstoremappeddata);
    source->AddReader(fiber33reader);
  }
  
  if (fTofd) 
  {
    auto tofdreader = new R3BTofdReader((EXT_STR_h101_TOFD_onion*)&ucesb_struct.tofd, offsetof(EXT_STR_h101_t, tofd));
    tofdreader->SetOnline(NOTstoremappeddata);
    source->AddReader(tofdreader);
  }

  if (fNeuland){
     auto unpackneuland = new R3BNeulandTamexReader(&ucesb_struct.raw_nnp, offsetof(EXT_STR_h101, raw_nnp));
     //nl_reader->SetSkipTriggerTimes();
     source->AddReader(unpackneuland);
  }

  run->SetSource(source);

  // Runtime data base ------------------------------------
  FairRuntimeDb *rtdb = run->GetRuntimeDb();

  // Load parameter files -----------------------------------------------------
  // Root file
  Bool_t kParameterMerged = kFALSE;
  auto parIo1 = new FairParRootFileIo(kParameterMerged);
  TList *parList1 = new TList();
  if (fLos)
    parList1->Add(new TObjString(loscalfilename)); 
  if (fFrsSci)
    parList1->Add(new TObjString(mydir+"/sci2/sci2calpar_v2.root")); 


//  if (fRolu)
//    parList1->Add(new TObjString(pardir + "rolu/rolucalpar_v1.root"));
  if (fTofd) {
    //parList1->Add(new TObjString(pardir + "tofd/tofdcalpar_v1.root"));
    parList1->Add(new TObjString(mydir + "tofd/tofdhitpar_v5.root"));
  }
  if (fNeuland){
    parList1->Add(new TObjString("./params/neuland_mapping.root"));
    parList1->Add(new TObjString("./params/params_neuland_tcal_160124.root"));
    parList1->Add(new TObjString("./params/params_ig_fine2.root"));
  }
   if (fFib30 || fFib31 || fFib32 || fFib33){ 
    parList1->Add(new TObjString(mydir + "fibers/fibercalpar_v1.root"));  
    parList1->Add(new TObjString(mydir + "fibers/fiberhitpar_v2.root"));  
  }
  if (fRpc){
    parList1->Add(new TObjString(mydir + "rpc/rpcmergedparams_v2.root"));
  }
    
  parIo1->open(parList1, "in");
  rtdb->setFirstInput(parIo1);
  rtdb->print();

  // Ascii file
  auto parIo2 = new FairParAsciiFileIo();
  TList *parList2 = new TList();
//  parList2->Add(new TObjString(califacalfilename));
  parList2->Add(new TObjString(calfilename));
//  parList2->Add(new TObjString(footparfilename));
  //parList2->Add(new TObjString(footcalfilename));
//  parList2->Add(new TObjString(pardir + "alpide/alpide_mapping_v1.par"));
  parList2->Add(new TObjString(mydir + "tofd/tofd_mapping_v1.par"));
  parList2->Add(new TObjString("/lustre/land/mwhiteh/s091/R3BParams_s091_s118/califa/califacalpar_vNew.par"));

//  parList2->Add(new TObjString(mydir + "fibers/fibermappingpar_v1.par"));
  parList2->Add(new TObjString(mydir + "/tttx/tttxCalibParams_V2.par"));
  parList2->Add(new TObjString(mydir + "/parameters/CalibPar.par"));
  parIo2->open(parList2, "in");

  rtdb->setSecondInput(parIo2);

  FairParAsciiFileIo* parIoOut = new FairParAsciiFileIo();
  parIoOut->open(tempFootParFile, "out");
  rtdb->setOutput(parIoOut);

  rtdb->print();
  rtdb->addRun(fRunId);
  

  if (fFrsSci) {
    rtdb->getContainer("Sci2TCalPar");
    rtdb->setInputVersion(fRunId, (char *)"Sci2TCalPar", 1, 1);
  }
  if (fLos) {
     rtdb->getContainer("LosTCalPar");
     rtdb->setInputVersion(fRunId, (char *)"LosTCalPar", 1, 1);
  }
  if(fRolu){
    rtdb->getContainer("RoluTCalPar");
    rtdb->setInputVersion(fRunId, (char*)"RoluTCalPar", 1, 1);
  }
  if (fFib30){
    rtdb->getContainer("Fi30MAPMTTCalPar");
    rtdb->setInputVersion(fRunId, (char*)"Fi30MAPMTTCalPar", 1, 1);
    rtdb->getContainer("Fi30MAPMTTrigTCalPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi30MAPMTTrigTCalPar", 1, 1);
    rtdb->getContainer("Fi30HitPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi30HitPar", 1, 1);
  }
  if (fFib31){
    rtdb->getContainer("Fi31MAPMTTCalPar");
    rtdb->setInputVersion(fRunId, (char*)"Fi31MAPMTTCalPar", 1, 1);
    rtdb->getContainer("Fi31MAPMTTrigTCalPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi31MAPMTTrigTCalPar", 1, 1);
    rtdb->getContainer("Fi31HitPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi31HitPar", 1, 1);
  }
  if (fFib32){
    rtdb->getContainer("Fi32MAPMTTCalPar");
    rtdb->setInputVersion(fRunId, (char*)"Fi32MAPMTTCalPar", 1, 1);
    rtdb->getContainer("Fi32MAPMTTrigTCalPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi32MAPMTTrigTCalPar", 1, 1);
    rtdb->getContainer("Fi32HitPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi32HitPar", 1, 1);
  }
  if (fFib33){
    rtdb->getContainer("Fi33MAPMTTCalPar");
    rtdb->setInputVersion(fRunId, (char*)"Fi33MAPMTTCalPar", 1, 1);
    rtdb->getContainer("Fi33MAPMTTrigTCalPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi33MAPMTTrigTCalPar", 1, 1);
    rtdb->getContainer("Fi33HitPar");
    rtdb->setInputVersion(fRunId, (char *)"Fi33HitPar", 1, 1);
  }
  if (fRpc){
    rtdb->getContainer("RpcTCalPar");
    rtdb->setInputVersion(fRunId, (char*)"RpcTCalPar", 1, 1);
    rtdb->getContainer("RpcTotCalPar");
    rtdb->setInputVersion(fRunId, (char*)"RpcTotCalPar", 1, 1);
    rtdb->getContainer("RpcHitPar");
    rtdb->setInputVersion(fRunId, (char*)"RpcHitPar", 1, 1);
  }
  if (fTofd){
      rtdb->getContainer("tofdMappingPar");
      rtdb->setInputVersion(fRunId, (char *)"tofdMappingPar", 1, 1);
      rtdb->getContainer("TofdTCalPar");
      rtdb->setInputVersion(fRunId, (char *)"TofdTCalPar", 1, 1);
      rtdb->getContainer("tofdHitPar");
      rtdb->setInputVersion(fRunId, (char *)"tofdHitPar", 1, 1);
  }
  if (fNeuland){
    rtdb->getContainer("neulandMappingPar");
    rtdb->setInputVersion(999, (char*)"neulandMappingPar", 1, 1);
    rtdb->getContainer("LandTCalPar");
    rtdb->setInputVersion(999, (char*)"LandTCalPar", 1, 1);
    rtdb->getContainer("NeulandHitPar");
    rtdb->setInputVersion(999, (char*)"NeulandHitPar", 1, 1); 
  }
  rtdb->print();

  // Add analysis task ------------------------------------
  if (fFrsSci) {
    auto *frsSciMapped2Tcal = new R3BFrsSciMapped2Tcal();
    frsSciMapped2Tcal->SetOnline(NOTstorecaldata);
    run->AddTask(frsSciMapped2Tcal);
    
    auto frsSciTcal2Cal = new R3BFrsSciTcal2Cal();
    frsSciTcal2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(frsSciTcal2Cal);
  }
  
  // MWPC2: 1st mwpc located at the entrance of Cave-C
  if (fMwpc2) {
    auto MW2Map2Cal = new R3BMwpc2Mapped2Cal();
    MW2Map2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(MW2Map2Cal);

    auto MW2Cal2Hit = new R3BMwpc2Cal2Hit();
    MW2Cal2Hit->SetOnline(NOTstorehitdata);
    run->AddTask(MW2Cal2Hit);
  }

  // MWPC1: 2nd mwpc close to the vacuum chamber
  if (fMwpc1) {
    auto MW1Map2Cal = new R3BMwpc1Mapped2Cal();
    MW1Map2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(MW1Map2Cal);

    auto MW1Cal2Hit = new R3BMwpc1Cal2Hit();
    MW1Cal2Hit->SetOnline(NOTstorehitdata);
    run->AddTask(MW1Cal2Hit);
  }

  // MWPC0: 3rd mwpc close to GLAD
  if (fMwpc0) {
    auto MW0Map2Cal = new R3BMwpc0Mapped2Cal();
    MW0Map2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(MW0Map2Cal);

    auto MW0Cal2Hit = new R3BMwpc0Cal2Hit();
    MW0Cal2Hit->SetOnline(NOTstorehitdata);
    run->AddTask(MW0Cal2Hit);
  }
  
  // LOS
  if (fLos) {
    auto losMapped2Cal = new R3BLosMapped2Cal();
    losMapped2Cal->SetNofModules(1, 8);
    losMapped2Cal->SetTrigger(1);
    losMapped2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(losMapped2Cal);
    
    auto losstart = new R3BLosProvideTStart();
    run->AddTask(losstart);

    auto losCal2Hit = new R3BLosCal2Hit();
    losCal2Hit->SetTrigger(1);
    losCal2Hit->SetOnline(NOTstorehitdata);
    run->AddTask(losCal2Hit);
  }
  
  // MUSLI
  if (fMusli) {
	 
    auto MusliMap2Cal = new R3BMusliMapped2Cal();
    MusliMap2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(MusliMap2Cal);
		
    auto MusliCal2Hit = new R3BMusliCal2Hit();
    MusliCal2Hit->SetDirectBeta(0.904);
    MusliCal2Hit->SetOnline(NOTstorehitdata);
    run->AddTask(MusliCal2Hit);
  }

  if (fLos && fFrsSci){
    //NICHT GUT with CALIFA cal!!    
    auto IncomingBeta = new R3BIncomingBeta();
    //IncomingBeta->SetUseTref();
    IncomingBeta->SetLosRefCh(7);
    IncomingBeta->UseFrsSci(true);
    IncomingBeta->SetNumDets(1);
    IncomingBeta->SetStaId(1);
    IncomingBeta->SetStoId(0);
    IncomingBeta->SetOnline(NOTstorehitdata);
    run->AddTask(IncomingBeta);
  
   // auto IncomingBeta = new R3BIncomingBeta();
   // IncomingBeta->SetUseTref();
   // IncomingBeta->SetOnline(!NOTstorehitdata);
   // run->AddTask(IncomingBeta);
  }
 
  if (fTrackingFrs && fLos && fMusli && fFrsSci)
  {
        auto AnalysisIncomingID = new R3BAnalysisIncomingID();
        AnalysisIncomingID->SetOnline(NOTstorehitdata);
        AnalysisIncomingID->SetMusicForPID();
        AnalysisIncomingID->SetNumDet(1);
        run->AddTask(AnalysisIncomingID);
  }

/*  
  //// ROLU
  //if (fRolu) {
  //  auto roluMapped2Cal = new R3BRoluMapped2Cal();
  //  roluMapped2Cal->SetNofModules(1,4);
  //  roluMapped2Cal->SetTrigger(-1);
  //  roluMapped2Cal->SetOnline(NOTstorecaldata);
  //  run->AddTask( roluMapped2Cal );
  //}
 */

  // ALPIDE
  if (fAlpidecal){
  //if (false){
    auto alpidemap2cal = new R3BAlpideMapped2Cal();
    alpidemap2cal->SetOnline(NOTstorecaldata);
    run->AddTask(alpidemap2cal);
    
    auto alpidecal2hit = new R3BAlpideCal2Hit();
    alpidecal2hit->SetOnline(NOTstorehitdata);
    alpidecal2hit->SetAlgorithm("Default");
    run->AddTask(alpidecal2hit);
  }

  // FOOT
  
  if (fFoot) {  
    auto footMap2Cal = new R3BFootMapped2StripCal();
    footMap2Cal->SetThresholdSigma(5.);
    footMap2Cal->SetThresholds(10.,10.);
    footMap2Cal->SetRefRate(10000,10000);
    footMap2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(footMap2Cal);

    auto footCal2Hit = new R3BFootStripCal2Hit();
    footCal2Hit->SetOnline(NOTstorehitdata);
    footCal2Hit->SetClusterEnergy(0.);
    footCal2Hit->SetTimesSigmas(5.);
    footCal2Hit->SetMaxNumClusters(200);
    run->AddTask(footCal2Hit);

    //auto footonline = new R3BFootOnlineSpectra();

  //  footonline->SetBinParams(-100, 5000, 1000);
//    run->AddTask(footonline);
//    run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);

//    rtdb->saveOutput();

  }

  // CALIFA
  if (fCalifa){
    
	  // R3BCalifaMapped2CrystalCal ---
    auto califaMap2Cal = new R3BCalifaMapped2CrystalCal();
    califaMap2Cal->SetOnline(NOTstorehitdata);
    run->AddTask(califaMap2Cal);
    
    // R3BCalifaCrystalCal2Cluster ---
    auto Cal2Clus = new R3BCalifaCrystalCal2Cluster();
    Cal2Clus->SetCrystalThreshold(100.);//100keV
    Cal2Clus->SetProtonClusterThreshold(12000.);//12MeV
    Cal2Clus->SetGammaClusterThreshold(500.);//2MeV
    Cal2Clus->SetRoundWindow(0.15);//0.3
    //Cal2Clus->SelectGeometryVersion(2024);
    
    Cal2Clus->SetRandomization(kTRUE);
    //Cal2Clus->SetRandomizationFile("/u/mwhiteh/s091/R3BParams_s091_s118/califa/angular_distributions.root");
    //Cal2Clus->SetRandomizationFile("/u/mwhiteh/s091/R3BParams_s091_s118/califa/angular_histograms_shifted.root");
    //Cal2Clus->SetRandomizationFile("/u/land/gmw537/angular_histograms_shifted.root");
    Cal2Clus->SetRandomizationFile("/lustre/land/mwhiteh/angular_histograms_v2.root");
    Cal2Clus->SetOnline(NOTstorehitdata);
    auto califa_geo = new R3BCalifaGeometry();
    califa_geo->Init(2024);
    run->AddTask(Cal2Clus);

    }
  
  // TTT10
  if (fTttx){
	R3BTttxMapped2Cal* tttxMap2Cal = new R3BTttxMapped2Cal();
	tttxMap2Cal->SetOnline(NOTstorecaldata);
	run->AddTask(tttxMap2Cal);
	
	R3BTttxCal2Hit* tttxCal2Hit = new R3BTttxCal2Hit();
	tttxCal2Hit->SetOnline(NOTstorehitdata);
	tttxCal2Hit->SetFixedBeta(0.714);
	tttxCal2Hit->SetTimeWindow(100);
	tttxCal2Hit->SetNumNeighbours(0);
	run->AddTask(tttxCal2Hit);
  }
  
  // RPC
  if (fRpc) { 
    auto rpcMapped2PreCal=new R3BRpcMapped2PreCal();
    rpcMapped2PreCal->SetFpgaCorrelationFile("/u/mwhiteh/s091/params/RPC.txt");
    rpcMapped2PreCal->SetOnline(NOTstorecaldata);
    run->AddTask(rpcMapped2PreCal);

    auto rpcPreCal2Cal = new R3BRpcPreCal2Cal();
    rpcPreCal2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(rpcPreCal2Cal);

    auto rpcCal2Hit = new R3BRpcCal2Hit();
    rpcCal2Hit->SetOnline(NOTstorehitdata);
    run->AddTask(rpcCal2Hit);
  }
  
  // Fibers
  if (fFib30){
    auto map2calFi30 = new R3BFiberMAPMTMapped2Cal("Fi30");
    map2calFi30->SetOnline(NOTstorecaldata);
    run->AddTask(map2calFi30);
    auto cal2hitFi30 = new R3BFiberMAPMTCal2Hit("Fi30", 1, R3BFiberMAPMTCal2Hit::HORIZONTAL);
    cal2hitFi30->SetOrientation(R3BFiberMAPMTCal2Hit::INVERTED);
    cal2hitFi30->SetOnline(NOTstorehitdata);
    run->AddTask(cal2hitFi30);
  }
  
  if (fFib31){
    auto map2calFi31 = new R3BFiberMAPMTMapped2Cal("Fi31");
    map2calFi31->SetOnline(NOTstorecaldata);
    run->AddTask(map2calFi31);
    auto cal2hitFi31 = new R3BFiberMAPMTCal2Hit("Fi31", 1, R3BFiberMAPMTCal2Hit::VERTICAL);
    cal2hitFi31->SetOrientation(R3BFiberMAPMTCal2Hit::INVERTED);
    cal2hitFi31->SetOnline(NOTstorehitdata);
    run->AddTask(cal2hitFi31);
  }
  
  if (fFib32){
    auto map2calFi32 = new R3BFiberMAPMTMapped2Cal("Fi32");
    map2calFi32->SetOnline(NOTstorecaldata);
    run->AddTask(map2calFi32);
    auto cal2hitFi32 = new R3BFiberMAPMTCal2Hit("Fi32", 1, R3BFiberMAPMTCal2Hit::VERTICAL);
    cal2hitFi32->SetOnline(NOTstorehitdata);
    run->AddTask(cal2hitFi32);
  }
  
  if (fFib33){
    auto map2calFi33 = new R3BFiberMAPMTMapped2Cal("Fi33");
    map2calFi33->SetOnline(NOTstorecaldata);
    run->AddTask(map2calFi33);
    auto cal2hitFi33 = new R3BFiberMAPMTCal2Hit("Fi33", 1, R3BFiberMAPMTCal2Hit::VERTICAL);
    cal2hitFi33->SetOnline(NOTstorehitdata);
    run->AddTask(cal2hitFi33);
  }
  
  //TofD
  if (fTofd){
  /*  auto tofdMapped2Cal = new R3BTofDMapped2Cal();
    tofdMapped2Cal->SetOnline(NOTstorecaldata);
    run->AddTask(tofdMapped2Cal);

    auto tofdCal2Hit = new R3BTofDCal2Hit();
    tofdCal2Hit->SetOnline(NOTstorehitdata);
    tofdCal2Hit->SetTofdQ(6);          // Set nuclear charge Z of main beam
    tofdCal2Hit->SetTrigger(1);
    //tofdCal2Hit->SetTpat(1, 16); // Tpat -1, Tpat 1 = mult2, Tpat 2 = or12
    //tofdCal2Hit->SetTofdHisto(false); // create control histograms
    tofdCal2Hit->SetTofdTotPos(true); // use y-position from ToT instead of Tdiff
    run->AddTask(tofdCal2Hit);
*/  
    // ToFD Mapped2Cal
	auto tofdMapped2Cal = new R3BTofDMapped2Cal();
	//tofdMapped2Cal->SetOnline(true);
	tofdMapped2Cal->SetNofModules(4, 44);
	run->AddTask(tofdMapped2Cal);

	auto tofdCal2Hit = new R3BTofDCal2Hit();
	//tofdCal2Hit->SetOnline(true);
	//  tofdCal2Hit->SetTofdQ(6);          // Set nuclear charge Z of main beam
	tofdCal2Hit->SetNofModules(4, 44); // planes, bars
	tofdCal2Hit->SetTrigger(-1);
	tofdCal2Hit->SetTpat(-1, -1); // Tpat -1, Tpat 1 = mult2, Tpat 2 = or12
	//tofdCal2Hit->SetTofdHisto(false); // create control histograms
	tofdCal2Hit->SetTofdTotPos(false); // use y-position from ToT instead of Tdiff
	run->AddTask(tofdCal2Hit);


  }
  
  //Neuland
  if(fNeuland)
  {
    auto tcal = new R3BNeulandMapped2Cal();
    tcal->SetTrigger(-1);
    tcal->SetNofModules(nPlanes, nBarsPerPlane);
    tcal->SetNhitmin(1);
    tcal->EnableWalk(true);
    //tcal->SubtractTriggerTime(false);
    run->AddTask(tcal);
    
    auto nlhit = new R3BNeulandCal2Hit();
    nlhit->SetDistanceToTarget(distanceToTarget);
    //nlhit->SetGlobalTimeOffset(-367.966); // until run 120 default -400
    //nlhit->SetGlobalTimeOffset(-392.81600); // from run 121
    nlhit->SetGlobalTimeOffset(-337.79); // from run 167
    nlhit->SetEnergyCutoff(0.0);
    run->AddTask(nlhit);
    
    auto r3bNeulandClusterFinder = new R3BNeulandClusterFinder(2. * 7.5, 2. * 7.5, 2. * 7.5, 1.);
    run->AddTask(r3bNeulandClusterFinder);
  }/*
  if (fMDF)
  {
  	auto tracking = new R3BTrackingWZ();

	tracking->SetTpat(-1);

	tracking->Set_MDF_PoQ("/u/mwhiteh/s091/MDF/MDF_S091/MDF_PoQ_1456A_20240221.txt");
	tracking->Set_MDF_TX0("/u/mwhiteh/s091/MDF/MDF_S091/MDF_TX0_1456A_20240221.txt");
	tracking->Set_MDF_TY0("/u/mwhiteh/s091/MDF/MDF_S091/MDF_TY0_1456A_20240221.txt");
	tracking->Set_MDF_TX1("/u/mwhiteh/s091/MDF/MDF_S091/MDF_TX1_1456A_20240221.txt");
	tracking->Set_MDF_TY1("/u/mwhiteh/s091/MDF/MDF_S091/MDF_TY1_1456A_20240221.txt");
	tracking->Set_MDF_FlightPath("/u/mwhiteh/s091/MDF/MDF_S091/MDF_FlightPath_1456A_20240221.txt");

	tracking->SetGladCurrent(GLAD_CURRENT);//value from this run
	tracking->SetGladReferenceCurrent(GLAD_CURRENT_REF);//same as for MDF training

	//this reference is needed only for the alignment mode
	tracking->SetReferencePoQ(3.9441667);//from LISE++: 12C @ 1250MeV/u

	// Expeirmental cuts
	tracking->SetFiberEnergyMinMax(0,10000);//min, max energy of Fiber hits
	//tracking->SetFiberTimeMinMax(-10000,10000);//min, max energy of Fiber hits
	tracking->SetFiberTimeMinMax(-300,-200);//min, max energy of Fiber hits

	//-------- Set nominal lab position of the detectors
	tracking->SetPositionMwpc0(0, 0,  Target_to_MWPC0);
	tracking->SetPositionMwpc1(2.81479, -3.14525 , Target_to_MWPC1);

	tracking->SetPositionFiber32(
			TP_to_F32 * TMath::Sin(Angle) + 2.12681,
			0,
			Target_to_TP + TP_to_F32 * TMath::Cos(Angle));

	tracking->SetPositionFiber30(
			TP_to_F30 * TMath::Sin(Angle),
			3.60796,
			Target_to_TP + TP_to_F30 * TMath::Cos(Angle));

	tracking->SetPositionFiber33(
			TP_to_F33 * TMath::Sin(Angle) + F33_offset*TMath::Cos(18. * TMath::Pi() / 180.) + 1.90317,
			0,
			Target_to_TP + TP_to_F33 * TMath::Cos(Angle) + F33_offset*TMath::Sin(18. * TMath::Pi() / 180.));

	tracking->SetPositionFiber31(
			TP_to_F31 * TMath::Sin(Angle) + F31_offset*TMath::Cos(18. * TMath::Pi() / 180.) + 1.83582,
			0,
			Target_to_TP + TP_to_F31 * TMath::Cos(Angle) + F31_offset*TMath::Sin(18. * TMath::Pi() / 180.));

	//---------- Set nominal lab (Euler) rotations of the detectors
	tracking->SetEulerAnglesMwpc0(0, 0, 0);
	tracking->SetEulerAnglesMwpc1(0, 0, 0);
	tracking->SetEulerAnglesFiber32(0, Angle, 0);
	tracking->SetEulerAnglesFiber30(0, Angle, 0);
	tracking->SetEulerAnglesFiber33(0, Angle, 0);
	tracking->SetEulerAnglesFiber31(0, Angle, 0);

	run->AddTask(tracking);
  }
*/
  
  // Initialize -------------------------------------------
  run->Init();

  // Run --------------------------------------------------
  run->Run((nev < 0) ? nev : 0, (nev < 0) ? 0 : nev);

  //rtdb->saveOutput();

  // Finish -----------------------------------------------
  timer.Stop();
  Double_t rtime = timer.RealTime() / 60.;
  Double_t ctime = timer.CpuTime() / 60.;
  std::cout << std::endl << std::endl;
  std::cout << "Macro finished succesfully." << std::endl;
  std::cout << "Output file is " << outputFilename << std::endl;
  std::cout << "Real time " << rtime << " min, CPU time " << ctime << " min"
            << std::endl
            << std::endl;
  gSystem->Exec(Form("mv %s %s", tempFootParFile.Data(), footcalfilename.Data()));
   
  gApplication->Terminate();
}

