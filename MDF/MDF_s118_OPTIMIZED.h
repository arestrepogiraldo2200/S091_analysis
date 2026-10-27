const UInt_t NVarsFunctor = 8;

//Positions in [cm] w.r.t. Turning Point (TP) of GLAD along 18 deg line
const Double_t SPEED_OF_LIGHT = 29.9792458;                                            //cm/ns
const Double_t AMU = 0.9314940038;                                                     //GeV/c2
const Double_t Angle = (18.0 ) * TMath::Pi() / 180.;                                   //Central turning angle of GLAD
const Double_t TOFD_TIME_OFFSET = 11.74 +  73;                                         //0 offset + 17.4 m - before psps and after Fi12

//Nominal detector positions (starting values)
const Double_t Reference_Current = 1490.0;        //1456.0;   // 3583.81;     //3583.81;// 2641.0;//Current used to calculated MDFs (800 MeV/u Sn)
const Double_t Actual_Current = 1456;      //run 326                                                   // <-------------------

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
TVector3 fi33_position; 
TVector3 fi31_position; 

TVector3 mwpc1_angles; 
TVector3 fi32_angles; 
TVector3 fi33_angles; 
TVector3 fi31_angles; 

TVector3 mwpc1_pos_offset; 
TVector3 fi32_pos_offset; 
TVector3 fi33_pos_offset;
TVector3 fi31_pos_offset;

TVector3 mwpc1_ang_offset; 
TVector3 fi32_ang_offset; 
TVector3 fi33_ang_offset;
TVector3 fi31_ang_offset;


// Data structure to track the particles across the detectors

struct MDF_Track_Data 
{
    TVector3 mwpc1_lab; 
    TVector3 fi32_lab;
    TVector3 fi33_lab; 
    TVector3 fi31_lab; 
    TVector3 tofd_lab; 

    TVector3 mwpc1_dhit; 
    TVector3 fi32_dhit; 
    TVector3 fi33_dhit; 
    TVector3 fi31_dhit; 
    TVector3 tofd_dhit; 

    TVector3 mwpc1_dE_x; 
    TVector3 mwpc1_dE_y; 
    TVector3 fi32_dE; 
    TVector3 fi33_dE; 
    TVector3 fi31_dE; 
    TVector3 tofd_dE; 
    TVector3 los_dE; 

    TVector3 mwpc1_t; 
    TVector3 fi32_t; 
    TVector3 fi33_t; 
    TVector3 fi31_t; 
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

    // Optimum parameters from self-consistent alignement
    // double par[14] = {-3.38931, 4.81549, 
    //                   -0.907715, -4.98417, -0.488445, -0.498736,
    //                   -3.21726, 0.369934, -0.49921, 0.370774, 
    //                   -3.05182, 1.6021, -0.213565, -0.48408};

    // double par[14] = {-2.56985, 0.815912, 
    //                   0.991082, 2.46896, 0.406464, 0.259915, 
    //                   0.78244, 3.42896, 0.433834, 0.124068, 
    //                   1.15385, 2.11223, 0.337565, 0.282816};

    // double par[14] = {3.89911, 0.77072,
    //                 -1.48251, 3.20584, -0.368437, 0.204058,
    //                 -0.62791, -3.59123, 0.390708,0.333552,
    //                 -1.2314,-3.21318,-0.258049,-0.321246};

    // double par[14] = {-2.31821, -3.91973, 
    //                   -1.61989, 4.61088, 0.130884, 0.346447, 
    //                    1.30293, -4.99543, 0.363009, 0.241384, 
    //                   -1.19909, 2.91063, 0.27135, 0.0637986};

    // Constrained 
    double par[14] = {0.62705,  4.91053, 
                     -4.73033,  2.14528, -0.144427,  -0.10782,
                     -4.93043, -1.1072,   0.195477,  -0.193539,
                     -4.99979, -1.43749, -0.0373261, -0.0264633};


// double par1 = 4.5357;
// double par2 = -2.62659;

// double par3 = -0.127484;
// double par4 = 2.99808;
// double par5 = 0.0664501;
// double par6 = 0.160543;

// double par7 = 0.110671;
// double par8 = -1.05328;
// double par9 = 0.153397;
// double par10 = 0.160712;

// double par11 = -0.0697367;
// double par12 = -1.22957;
// double par13 = -0.0194854;
// double par14 = -0.182135;
// Best delta: 0.000503649


// double par[14] = {4.5357, -2.62659,
//                   -0.127484, 2.99808, 0.0664501, 0.160543,
//                    0.110671, -1.05328, 0.153397, 0.160712,
//                   -0.0697367, -1.22957, -0.0194854, -0.182135};
// // Best delta: 0.000594228


// double par[14] = { -4.99592, 4.92425, 
//                    -0.974042, 6.11133, 0.146641, -0.0208257, 
//                     3.30711, -6.5677, 0.199675 , 0.19999, 
//                    -1.15886, 7.33394, -0.0734445, 0.105279};
// Best delta: 0.000594228


    double par1 = par[0];      // mwpc x
    double par2 = par[1];      // mwpc z

    double par3 = par[2];      // fib32 x
    double par4 = par[3];      // fib32 z
    double par5 = par[4];      // fib32 y'
    double par6 = par[5];      // fib32 z'

    double par7 = par[6];      // fib33 x
    double par8 = par[7];      // fib33 z
    double par9 = par[8];      // fib33 y'
    double par10 = par[9];     // fib33 z'

    double par11 = par[10];    // fib31 x 
    double par12 = par[11];    // fib31 z 
    double par13 = par[12];    // fib31 y'
    double par14 = par[13];    // fib31 z'

    // Set the positions of central values of the planes of the detectors. Ideal case and orientations
    mwpc1_position.SetXYZ(0,0,targ_to_mwpc1);
    fi32_position.SetXYZ(fi_32_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_32_to_TP * TMath::Cos(AngleAl));
    fi33_position.SetXYZ(fi_33_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_33_to_TP * TMath::Cos(AngleAl));
    fi31_position.SetXYZ(fi_31_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_31_to_TP * TMath::Cos(AngleAl));        // Z position of FIB31 is the Z position of FIB33. Projection of 31 in 33.

    mwpc1_angles.SetXYZ(0, 0, 0);
    fi32_angles.SetXYZ(0, AngleAl, 0);
    fi33_angles.SetXYZ(0, AngleAl, 0);
    fi31_angles.SetXYZ(0, AngleAl, 0);

    mwpc1_pos_offset.SetXYZ(par1, 0,  par2);
    fi32_pos_offset.SetXYZ( par3, 0,  par4);
    fi33_pos_offset.SetXYZ( par7, 0,  par8);
    fi31_pos_offset.SetXYZ(par11, 0,  par12);    // The z-position of FIB31 is projected into the plane of FIB33

    mwpc1_ang_offset.SetXYZ(0, 0, 0);
    fi32_ang_offset.SetXYZ(0,  par5,  par6 );
    fi33_ang_offset.SetXYZ(0,  par9,  par10);
    fi31_ang_offset.SetXYZ(0, par13,  par14);
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
