const UInt_t NVarsFunctor = 8;

//Positions in [cm] w.r.t. Turning Point (TP) of GLAD along 18 deg line
const Double_t SPEED_OF_LIGHT = 29.9792458;                                            //cm/ns
const Double_t AMU = 0.9314940038;                                                     //GeV/c2
const Double_t Angle = (18.0 ) * TMath::Pi() / 180.;                                   //Central turning angle of GLAD
const Double_t TOFD_TIME_OFFSET = 11.74 +  73;                                         //0 offset + 17.4 m - before psps and after Fi12

//Nominal detector positions (starting values)
const Double_t Reference_Current = 1490.0;        //1456.0;   // 3583.81;     //3583.81;// 2641.0;//Current used to calculated MDFs (800 MeV/u Sn)
const Double_t Actual_Current = 1490.0;      //run 326                                                   // <-------------------

// Distances in cm 
const Double_t Target_to_GLAD_flange  = 127.9;                                         // <-------------------
const Double_t targ_to_TP   = 293.1;                                                   // <-------------------
const Double_t TP_to_GLAD_exit_flange = 385.0 - 125.0;                                 // <-------------------
const Double_t Target_to_RP   = Target_to_GLAD_flange +  54.05 - 0.55580;              // <-------------------

const Double_t targ_to_mwpc1 = 95.7;                              // target to MWPC1                 // <-------------------
const Double_t fi_32_to_TP   = 385.0;                            // turning point to fib32           // <-------------------
const Double_t fi_31_to_TP   = TP_to_GLAD_exit_flange + 242.5;   // turning point to fib31           // <-------------------
const Double_t fi_33_to_TP   = TP_to_GLAD_exit_flange + 225.0;   // turning point to fib33           // <-------------------
const Double_t tofd_to_TP    = TP_to_GLAD_exit_flange + 661.8;   // turning point to fib33           // <-------------------



const Double_t AngleAl = (-18. ) * TMath::Pi() / 180.;
const Double_t AngleAl11 = (-18.05  )* TMath::Pi() / 180.;

TVector3 mwpc1_position; 
TVector3 fi32_position; 
TVector3 fi31_33_position; 


TVector3 mwpc1_angles; 
TVector3 fi32_angles; 
TVector3 fi31_33_angles; 


TVector3 mwpc1_pos_offset; 
TVector3 fi32_pos_offset; 
TVector3 fi31_33_pos_offset;


TVector3 mwpc1_ang_offset; 
TVector3 fi32_ang_offset; 
TVector3 fi31_33_ang_offset;


// Data structure to track the particles across the detectors

struct MDF_Track_Data 
{
    TVector3 mwpc1_lab; 
    TVector3 fi32_lab; 
    TVector3 fi31_33_lab; 
    TVector3 tofd_lab; 

    TVector3 mwpc1_dhit; 
    TVector3 fi32_dhit; 
    TVector3 fi31_33_dhit; 
    TVector3 tofd_dhit; 

    TVector3 mwpc1_dE_x; 
    TVector3 mwpc1_dE_y; 
    TVector3 fi32_dE; 
    TVector3 fi31_33_dE; 
    TVector3 tofd_dE; 
    TVector3 los_dE; 

    TVector3 mwpc1_t; 
    TVector3 fi32_t; 
    TVector3 fi31_33_t; 
    TVector3 tofd_t;    
    TVector3 los_t;    

    Int_t tofd_plane;
    Double_t ToF;
    Double_t Beta;
    Double_t Gamma;
    Double_t AoZ;

    Double_t PoQ_mdf;
    Double_t FlightPath_mdf;
};

// MDF training variables
R3BMDFWrapper* MDF_PoQ;
R3BMDFWrapper* MDF_TX0;
R3BMDFWrapper* MDF_TY0;
R3BMDFWrapper* MDF_FlightPath;

MDF_Track_Data Candidate;
std::vector<MDF_Track_Data> mdf_tracks;         //Storing candidates of Global mdf_tracks

R3BMDFWrapper* MDF_PoQ_glob;
R3BMDFWrapper* MDF_TX0_glob;
R3BMDFWrapper* MDF_FlightPath_glob;

void SetMDFOffsets(){

    mwpc1_position.SetXYZ(0,0,targ_to_mwpc1);
    fi31_33_position.SetXYZ(fi_33_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_33_to_TP * TMath::Cos(AngleAl));       // By convention, we project all 31 hits in 33 plane 
    fi32_position.SetXYZ(fi_32_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_32_to_TP * TMath::Cos(AngleAl));

    mwpc1_angles.SetXYZ(0, 0, 0);
    fi31_33_angles.SetXYZ(0, AngleAl, 0);
    fi32_angles.SetXYZ(0, AngleAl, 0);

    mwpc1_pos_offset.SetXYZ(0, 0, 0);
    fi31_33_pos_offset.SetXYZ(0, 0, 0);
    fi32_pos_offset.SetXYZ(0, 0, 0);

    mwpc1_ang_offset.SetXYZ(0, 0, 0);
    fi31_33_ang_offset.SetXYZ(0, 0, 0);
    fi32_ang_offset.SetXYZ(0, 0, 0);
}


void TransformPoint(TVector3& point, TVector3 rot, TVector3 trans)
{
    TRotation r;
    // First Euler rotation around Y axis
    r.RotateY(rot.Y());

    // get local X axis after first rotation
    TVector3 v3_localX(1, 0, 0);
    v3_localX.SetTheta(r.ThetaX());
    v3_localX.SetPhi(r.PhiX());

    // Second Euler rotation around local X axis
    r.Rotate(rot.X(), v3_localX);

    // get local Z axis after second rotation
    TVector3 v3_localZ(0, 0, 1);
    v3_localZ.SetTheta(r.ThetaZ());
    v3_localZ.SetPhi(r.PhiZ());

    // final rotation around local Z axis
    r.Rotate(rot.Z(), v3_localZ);

    // Applying rotation
    point.Transform(r);

    // Applying translation
    point += trans;

    return;
}

double CalculatePoQ(const TVector3& mwpc1_lab,
                    const TVector3& fi31_33_lab, const TVector3& fi32_lab,
                    double current, double reference, R3BMDFWrapper* mdf){
    
    // VARIABLES BASED ON MY OWN TRAINING OF THE MDF
    Double_t edata[8];
    edata[0] = mwpc1_lab.X();
    edata[1] = mwpc1_lab.Z();
    edata[2] = fi32_lab.X();
    edata[3] = fi32_lab.Z();
    edata[4] = fi31_33_lab.X();
    edata[5] = fi31_33_lab.Z();
    edata[6] = (fi31_33_lab.X()-fi32_lab.X())/(fi31_33_lab.Z()-fi32_lab.Z());

    // // VARIABLES BASED ON PREVIOUS TRAINING OF THE MDF 
    // Double_t edata[7];
    // edata[0] = mwpc1_lab.X();
    // edata[1] = mwpc1_lab.Y();
    // edata[2] = mwpc1_lab.Z();
    // edata[3] = fi32_lab.X();
    // edata[4] = fi32_lab.Z();
    // edata[5] = (fi31_33_lab.X()-fi32_lab.X())/(fi31_33_lab.Z()-fi32_lab.Z());
    // edata[6] = 0;

    return mdf->MDF(edata) * (current / reference);
}
