extern TTree* tree;
extern TTree* tree227;
extern TTree* tree223;
extern TTree* tree229;
R3BMDFWrapper* MDF_PoQ;
R3BMDFWrapper* MDF_PoQ10;
R3BMDFWrapper* MDF_PoQ1012;
R3BMDFWrapper* MDF_FP1012;

// TVector3 fi11_lab_223, fi12_lab_223, fi11_lab_227, fi12_lab_227, fi11_lab_225, fi12_lab_225, fi11_lab_229, fi12_lab_229, fi11_lab_212, fi12_lab_212, fi11_lab_220, fi12_lab_220, fi11_lab_201, fi12_lab_201, fi11_lab_202, fi12_lab_202;
// TVector3 fi10_lab_223, fi10_lab_227, fi10_lab_225, fi10_lab_229, fi10_lab_212, fi10_lab_220, fi10_lab_201, fi10_lab_202;
// TVector3 psp2_lab, psp3_lab, fi11_lab, fi12_lab, fi13_lab, fi10_lab;
    
TVector3 fi32_lab, fi31_lab, fi33_lab, mwpc1_lab;

// TVector3 psp2_position, psp3_position, fib11_position, fib10_position, fib12_position, fib13_position;
// TVector3 psp2_angles, psp3_angles, fib11_angles, fib10_angles, fib12_angles, fib13_angles; 

TVector3 fi32_position, fi31_position, fi33_position, mwpc1_position;
TVector3 fi32_angles  , fi31_angles  , fi33_angles  , mwpc1_angles  ;

// TVector3 psp2_ang_offset, psp3_ang_offset, fib11_ang_offset, fib10_ang_offset, fib12_ang_offset, fib13_ang_offset;
// TVector3 psp2_pos_offset, psp3_pos_offset, fib11_pos_offset, fib10_pos_offset, fib12_pos_offset, fib13_pos_offset;

TVector3 fi32_pos_offset, fi31_pos_offset, fi33_pos_offset, mwpc1_pos_offset;
TVector3 fi32_ang_offset, fi31_ang_offset, fi33_ang_offset, mwpc1_ang_offset;

// const Double_t targ_to_psp2 = 81.77-10.+1.86121;//target to PSP2
// const Double_t targ_to_psp3 = 105.77-10.+1.876;//target to PSP3
// const Double_t targ_to_TP   = 116.84 + 165.;//(Target to GLAD entrance flange?)+(flange to TP form daniel)
// const Double_t fi13_to_TP   = 1301.17 + 47.;//(TP to exit flange of vac. pipe) + (flange to det)

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

// //0.5mm LOS - C density 1.851 - 100um stainless steel
// double Reference_PoQ_229 = 3.4560406;
// double Reference_PoQ_202 = 2.2147714;
// double Reference_PoQ_227 = 2.8963644;
// double Reference_PoQ_223 = 3.3233172;
// double Reference_PoQ_201 = 1.994957;
// double Reference_PoQ_225 = 3.1909179;
// double Reference_PoQ_212 = 2.7089167;
// double Reference_PoQ_220 = 3.0160662;

double Reference_PoQ_10 = 1.527872;
double Reference_PoQ_11 = 1.7092484;
double Reference_PoQ_12 = 1.9902324;
double Reference_PoQ_13 = 2.0717421;
double Reference_PoQ_14 = 2.2532898;

//0.5mm LOS - C density 1.851 - 200um stainless steel - just 229,202,227
/*double Reference_PoQ_229 = 3.4474683;
double Reference_PoQ_202 = 2.2014897;
double Reference_PoQ_227 = 2.8864796;
double Reference_PoQ_223 = 3.3153014;
double Reference_PoQ_201 = 1.9809929;
double Reference_PoQ_225 = 3.181675;
double Reference_PoQ_212 = 2.6983719;
double Reference_PoQ_220 = 3.0065355;*/
//0.5mm LOS - C density 1.851 - 150um stainless steel - just 229,202,227
/*double Reference_PoQ_229 = 3.4517581;
double Reference_PoQ_202 = 2.2081441;
double Reference_PoQ_227 = 2.8914251;
double Reference_PoQ_223 = 3.3145758;
double Reference_PoQ_201 = 1.9797101;
double Reference_PoQ_225 = 3.180928;*/
//0.5mm LOS - C density 1.851 - 120um stainless steel - just 229,202,227
/*double Reference_PoQ_229 = 3.4543278;
double Reference_PoQ_202 = 2.2121217;
double Reference_PoQ_227 = 2.8947233;
double Reference_PoQ_223 = 3.3145758;
double Reference_PoQ_201 = 1.9797101;
double Reference_PoQ_225 = 3.180928;*/
//0.5mm LOS - C density 1.851 - 90um stainless steel - just 229,202,227
/*double Reference_PoQ_229 = 3.456897;
double Reference_PoQ_202 = 2.2160922;
double Reference_PoQ_227 = 2.8973513;
double Reference_PoQ_223 = 3.3145758;
double Reference_PoQ_201 = 1.9797101;
double Reference_PoQ_225 = 3.180928;*/

// double Reference_Current = 2641.0;
//double Reference_PoQ_212 = 2.6983719;
//double Reference_PoQ_220 = 3.0065355;

double Reference_Current = 1490.0;
double Reference_PoQ, Actual_Current;

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

// double CalculatePoQ(const TVector3& fi11_lab, const TVector3& fi12_lab,
//                     const TVector3& psp2_lab, const TVector3& psp3_lab,
//                     const TVector3& fi13_lab, double current, double reference, R3BMDFWrapper* mdf)
// {
//     Double_t edata[9];
//     edata[0] = psp3_lab.X();
//     edata[1] = psp3_lab.Y();
//     edata[2] = psp3_lab.Z();
//     edata[3] = (psp3_lab.X() - psp2_lab.X()) / (psp3_lab.Z() - psp2_lab.Z());
//     edata[4] = fi12_lab.X();
//     edata[5] = fi12_lab.Z();
//     edata[6] = (fi12_lab.X() - fi11_lab.X()) / (fi12_lab.Z() - fi11_lab.Z());
//     edata[7] = fi13_lab.Y();
//     edata[8] = fi13_lab.Z();

//     return mdf->MDF(edata) * (current / reference);
// }

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

// double CalculatePoQ1012(const TVector3& fi10_lab, const TVector3& fi12_lab,
//                     const TVector3& psp2_lab, const TVector3& psp3_lab,
//                     const TVector3& fi13_lab, double current, double reference, R3BMDFWrapper* mdf)
// {
//     /*Double_t edata[9];
//     edata[0] = psp3_lab.X();
//     edata[1] = psp3_lab.Y();
//     edata[2] = psp3_lab.Z();
//     edata[3] = (psp3_lab.X() - psp2_lab.X()) / (psp3_lab.Z() - psp2_lab.Z());
//     edata[4] = fi10_lab.X();
//     edata[5] = fi10_lab.Z();
//     edata[6] = (fi12_lab.X() - fi10_lab.X()) / (fi12_lab.Z() - fi10_lab.Z());
//     edata[7] = fi13_lab.Y();
//     edata[8] = fi13_lab.Z();*/
//     Double_t edata[10];    
//     edata[0] = psp3_lab.X();
//     edata[1] = psp3_lab.Y();
//     edata[2] = psp3_lab.Z();
//     edata[3] = (psp3_lab.X() - psp2_lab.X()) / (psp3_lab.Z() - psp2_lab.Z());
//     edata[4] = fi10_lab.X();
//     edata[5] = fi10_lab.Z();
//     edata[6] = fi12_lab.X();
//     edata[7] = fi12_lab.Z();
//     edata[8] = fi13_lab.Y();
//     edata[9] = fi13_lab.Z();    

//     return mdf->MDF(edata) * (current / reference);
// }

// double CalculatePoQ1012_11(const TVector3& fi10_lab, const TVector3& fi12_lab,
//                     const TVector3& psp2_lab, const TVector3& psp3_lab,
//                     const TVector3& fi13_lab, double current, double reference, R3BMDFWrapper* mdf)
// {
//     /*
//     Double_t edata[11];    
//     edata[0] = psp2_lab.X();
//     edata[1] = psp2_lab.Z();
//     edata[2] = psp3_lab.X();
//     edata[3] = psp3_lab.Y();
//     edata[4] = psp3_lab.Z();
//     edata[5] = fi10_lab.X();
//     edata[6] = fi10_lab.Z();
//     edata[7] = fi12_lab.X();
//     edata[8] = fi12_lab.Z();
//     edata[9] = fi13_lab.Y();
//     edata[10] = fi13_lab.Z();  */
    
//     Double_t edata[9];
//     edata[0] = psp3_lab.X();
//     edata[1] = psp3_lab.Y();
//     edata[2] = psp3_lab.Z();
//     edata[3] = (psp3_lab.X() - psp2_lab.X()) / (psp3_lab.Z() - psp2_lab.Z());
//     edata[4] = fi10_lab.X();
//     edata[5] = fi10_lab.Z();
//     edata[6] = (fi12_lab.X() - fi10_lab.X()) / (fi12_lab.Z() - fi10_lab.Z());
//     edata[7] = fi13_lab.Y();
//     edata[8] = fi13_lab.Z();

//     return mdf->MDF(edata) * (current / reference);
// }

// double CalculatePoQ1012_09(const TVector3& fi10_lab, const TVector3& fi12_lab,
//                     const TVector3& psp2_lab, const TVector3& psp3_lab,
//                     const TVector3& fi13_lab, double current, double reference, R3BMDFWrapper* mdf)
// {
//     Double_t edata[9];    
//     edata[0] = psp3_lab.X();
//     edata[1] = psp3_lab.Y();
//     edata[2] = psp3_lab.Z();
//     edata[3] = fi10_lab.X();
//     edata[4] = fi10_lab.Z();
//     edata[5] = fi12_lab.X();
//     edata[6] = fi12_lab.Z();
//     edata[7] = fi13_lab.Y();
//     edata[8] = fi13_lab.Z();    

//     return mdf->MDF(edata) * (current / reference);
// }

// double CalculateFP1012(const TVector3& fi10_lab, const TVector3& fi12_lab,
//                     const TVector3& psp2_lab, const TVector3& psp3_lab,
//                     const TVector3& fi13_lab, R3BMDFWrapper* FP)
// {
//     Double_t edata[11];    
//     edata[0] = psp2_lab.X();
//     edata[1] = psp2_lab.Z();
//     edata[2] = psp3_lab.X();
//     edata[3] = psp3_lab.Y();
//     edata[4] = psp3_lab.Z();
//     edata[5] = fi10_lab.X();
//     edata[6] = fi10_lab.Z();
//     edata[7] = fi12_lab.X();
//     edata[8] = fi12_lab.Z();
//     edata[9] = fi13_lab.Y();
//     edata[10] = fi13_lab.Z();    

//     return FP->MDF(edata) + targ_to_psp2 + 63. +200. -151.;
// }


// double MinimizeDelta(const double* par) {

//     //const Double_t AngleAl = (-18.) * TMath::Pi() / 180.;
//     // double Actual_Current_227 = 2166.0; //Run227
//     // double Actual_Current_223 = 2527.0; //Run223 
//     // double Actual_Current_225 = 2413.0; //Run225
//     // double Actual_Current_229 = 2641.0; //Run229
//     // double Actual_Current_212 = 2071.0; //Run207
//     // double Actual_Current_220 = 2299.0; //Run219
//     // double Actual_Current_201 = 1520.0; //Run201
//     // double Actual_Current_202 = 1691.0; //Run202


//     double Actual_Current_227 = 2166.0; //Run227

//     double fval;

//     // Example decryption
//     double par1 = par[0]; double par2 = par[1];
//     double par3 = par[2]; double par4 = par[3];
//     double par5 = par[4]; double par6 = par[5];
//     double par7 = par[6];

//     Double_t fi10_to_TP   = 255. + 201.3 + par2;//(TP to vac.cham. internal wall) + (wall to det)
//     Double_t fi11_to_TP   = 255. + 65.3 + par1 ;//(TP to vac.cham. internal wall) + (wall to det)
//     Double_t fi12_to_TP   = 1301.17 + 18. ;//(TP to exit flange of vac. pipe) + (flange to det)

//     Double_t AngleAl = (-18.0 +par7 ) * TMath::Pi() / 180.;
    
//     fi11_lab_227.SetXYZ(-1.*(-4.60759 + par3) , 0 , 0 );
//     fi12_lab_227.SetXYZ(-1.*(-16.8725 + par4) , 0 , 0 );
//     fi10_lab_227.SetXYZ(-1.*(-5.20989 + par4) , 0 , 0 );  

//     fi11_lab_223.SetXYZ(-1.*(-3.20910 + par3) , 0 , 0 );
//     fi12_lab_223.SetXYZ(-1.*(-1.12588 + par4) , 0 , 0 );         
//     fi10_lab_223.SetXYZ(-1.*(-3.20786 + par4) , 0 , 0 ); 

//     fi11_lab_225.SetXYZ(-1.*(-3.64903 + par3) , 0 , 0 );
//     fi12_lab_225.SetXYZ(-1.*(-13.0462 + par4) , 0 , 0 );
//     fi10_lab_225.SetXYZ(-1.*(-3.84497 + par4) , 0 , 0 ); 
    
//     fi11_lab_229.SetXYZ(-1.*(-2.73996 + par3) , 0 , 0 );
//     fi12_lab_229.SetXYZ(-1.*(-9.26699 + par4) , 0 , 0 );
//     fi10_lab_229.SetXYZ(-1.*(-2.53834 + par4) , 0 , 0 ); 
    
//     fi11_lab_212.SetXYZ(-1.*(-2.83415 + par3) , 0 , 0 );
//     fi12_lab_212.SetXYZ(-1.*(-9.74823 + par4) , 0 , 0 ); 
//     fi10_lab_212.SetXYZ(-1.*(-2.68217 + par4) , 0 , 0 ); 
                    
//     fi11_lab_220.SetXYZ(-1.*(-2.83545 + par3) , 0 , 0 );
//     fi12_lab_220.SetXYZ(-1.*(-10.0403 + par4) , 0 , 0 );
//     fi10_lab_220.SetXYZ(-1.*(-2.71734 + par4) , 0 , 0 );
    
//     fi11_lab_201.SetXYZ(-1.*(-2.93167 + par3) , 0 , 0 );
//     fi12_lab_201.SetXYZ(-1.*(-10.3565 + par4) , 0 , 0 );
//     fi10_lab_201.SetXYZ(-1.*(-2.85664 + par4) , 0 , 0 );
    
//     fi11_lab_202.SetXYZ(-1.*(-2.80716 + par3) , 0 , 0 );
//     fi12_lab_202.SetXYZ(-1.*(-9.61590 + par4) , 0 , 0 );                                         
//     fi10_lab_202.SetXYZ(-1.*(-2.63749 + par4) , 0 , 0 );
    
//     psp2_lab.SetXYZ(0 , 0 , 0);
//     psp3_lab.SetXYZ(0 , 0 , 0);
//     fi11_lab.SetXYZ(0 , 0 , 0);
//     fi12_lab.SetXYZ(0 , 0 , 0);
//     fi13_lab.SetXYZ(0 , 0 , 0); 

//     fib11_ang_offset.SetXYZ(0, par5, 0);
//     fib11_pos_offset.SetXYZ(0, 0 ,0);  

//     fib12_ang_offset.SetXYZ(0, par6, 0);
//     fib12_pos_offset.SetXYZ(0, 0, 0);  

//     fib10_ang_offset.SetXYZ(0, par6, 0);
//     fib10_pos_offset.SetXYZ(0, 0, 0);          
    
//     psp2_ang_offset.SetXYZ(0 , 0 , 0);
//     psp3_ang_offset.SetXYZ(0 , 0 , 0);
//     fib13_ang_offset.SetXYZ(0 , 0 , 0);
//     psp2_pos_offset.SetXYZ(0 , 0 , 0);
//     psp3_pos_offset.SetXYZ(0 , 0 , 0);
//     fib13_pos_offset.SetXYZ(0 , 0 , 0);    


//     psp2_position.SetXYZ(0,0,targ_to_psp2);
//     psp3_position.SetXYZ(0,0,targ_to_psp3);    
//     fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
//     fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
//     fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl));
//     fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl));      
        
//     fib10_angles.SetXYZ(0, AngleAl, 0);
//     fib11_angles.SetXYZ(0, AngleAl, 0);
//     fib12_angles.SetXYZ(0, AngleAl, 0);
//     fib13_angles.SetXYZ(0, AngleAl, 0); 
    
//     TransformPoint(fi11_lab_223, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_223, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//     TransformPoint(fi11_lab_227, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_227, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//     TransformPoint(fi11_lab_225, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_225, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//     TransformPoint(fi11_lab_229, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_229, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
//     TransformPoint(fi11_lab_212, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_212, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//     TransformPoint(fi11_lab_220, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_220, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//     TransformPoint(fi11_lab_201, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_201, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
//     TransformPoint(fi11_lab_202, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//     TransformPoint(fi12_lab_202, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
//     TransformPoint(fi10_lab_202, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//     TransformPoint(fi10_lab_201, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//     TransformPoint(fi10_lab_212, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//     TransformPoint(fi10_lab_220, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//     TransformPoint(fi10_lab_223, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//     TransformPoint(fi10_lab_225, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//     TransformPoint(fi10_lab_227, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//     TransformPoint(fi10_lab_229, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    
//     TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
//     TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
//     TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset ); 

//     double delta = 0.0;
//     delta += fabs(Reference_PoQ_227 - CalculatePoQ(fi11_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ));


//     //delta += fabs(Reference_PoQ_227 - CalculatePoQ(fi11_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_223 - CalculatePoQ(fi11_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_225 - CalculatePoQ(fi11_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_229 - CalculatePoQ(fi11_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_201 - CalculatePoQ(fi11_lab_201, fi12_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_202 - CalculatePoQ(fi11_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_212 - CalculatePoQ(fi11_lab_212, fi12_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ));    
//     //delta += fabs(Reference_PoQ_227 - CalculatePoQ(fi11_lab_227, fi10_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_223 - CalculatePoQ(fi11_lab_223, fi10_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_225 - CalculatePoQ(fi11_lab_225, fi10_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_229 - CalculatePoQ(fi11_lab_229, fi10_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ10)); 
//     //delta += fabs(Reference_PoQ_201 - CalculatePoQ(fi11_lab_201, fi10_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_202 - CalculatePoQ(fi11_lab_202, fi10_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_212 - CalculatePoQ(fi11_lab_212, fi10_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ10));           


//    return fval = delta;
// }

// double MinimizeDelta2(const double* par) {

//     double fval;
//     // Beispielhafte Entschlüsselung:
//     double par1 = par[0]; double par2 = par[1];
//     double par3 = par[2]; double par4 = par[3];
//     double par5 = par[4]; double par6 = par[5];
//     double par7 = par[6]; double par8 = par[7];
//     double par9 = par[8]; double par10 = par[9];
//     double par11 = par[10]; 
//     Double_t fi10_to_TP   = 255. + 201.3 + par1 + par2;//(TP to vac.cham. internal wall) + (wall to det)
//     Double_t fi11_to_TP   = 255. + 65.3 + par1 ;//(TP to vac.cham. internal wall) + (wall to det)
//     Double_t fi12_to_TP   = 1301.17 + 18. + par7;//(TP to exit flange of vac. pipe) + (flange to det)

//     double Actual_Current_227 = 2166.0; //Run227
//     double Actual_Current_223 = 2527.0; //Run223 
//     double Actual_Current_225 = 2413.0; //Run225
//     double Actual_Current_229 = 2641.0; //Run229
//     double Actual_Current_212 = 2071.0; //Run207
//     double Actual_Current_220 = 2299.0; //Run219
//     double Actual_Current_201 = 1520.0; //Run201
//     double Actual_Current_202 = 1691.0; //Run202

//         Double_t AngleAl = (-18.0 + par10 ) * TMath::Pi() / 180.;
        
//         fi11_lab_227.SetXYZ(-1.*(-4.60759 + par3) , 0 , 0 );
//         fi12_lab_227.SetXYZ(-1.*(-16.8725 + par8) , 0 , 0 );
//         fi10_lab_227.SetXYZ(-1.*(-5.20989 + par4) , 0 , 0 );  

//         fi11_lab_223.SetXYZ(-1.*(-3.20910 + par3) , 0 , 0 );
//         fi12_lab_223.SetXYZ(-1.*(-11.2588 + par8) , 0 , 0 );         
//         fi10_lab_223.SetXYZ(-1.*(-3.20786 + par4) , 0 , 0 ); 

//         fi11_lab_225.SetXYZ(-1.*(-3.64903 + par3) , 0 , 0 );
//         fi12_lab_225.SetXYZ(-1.*(-13.0462 + par8) , 0 , 0 );
//         fi10_lab_225.SetXYZ(-1.*(-3.84497 + par4) , 0 , 0 ); 
        
//         fi11_lab_229.SetXYZ(-1.*(-2.73996 + par3) , 0 , 0 );
//         fi12_lab_229.SetXYZ(-1.*(-9.26699 + par8) , 0 , 0 );
//         fi10_lab_229.SetXYZ(-1.*(-2.53834 + par4) , 0 , 0 ); 
        
//         fi11_lab_212.SetXYZ(-1.*(-2.83415 + par3) , 0 , 0 );
//         fi12_lab_212.SetXYZ(-1.*(-9.74823 + par8) , 0 , 0 ); 
//         fi10_lab_212.SetXYZ(-1.*(-2.68217 + par4) , 0 , 0 ); 
                       
//         fi11_lab_220.SetXYZ(-1.*(-2.83545 + par3) , 0 , 0 );
//         fi12_lab_220.SetXYZ(-1.*(-10.0403 + par8) , 0 , 0 );
//         fi10_lab_220.SetXYZ(-1.*(-2.71734 + par4) , 0 , 0 );
        
//         fi11_lab_201.SetXYZ(-1.*(-2.93167 + par3) , 0 , 0 );
//         fi12_lab_201.SetXYZ(-1.*(-10.3565 + par8) , 0 , 0 );
//         fi10_lab_201.SetXYZ(-1.*(-2.85664 + par4) , 0 , 0 );
        
//         fi11_lab_202.SetXYZ(-1.*(-2.80716 + par3) , 0 , 0 );
//         fi12_lab_202.SetXYZ(-1.*(-9.61590 + par8) , 0 , 0 );                                         
//         fi10_lab_202.SetXYZ(-1.*(-2.63749 + par4) , 0 , 0 );
        
//         psp2_lab.SetXYZ(0 , 0 , 0);
//         psp3_lab.SetXYZ(0 , 0 , 0);
//         fi11_lab.SetXYZ(0 , 0 , 0);
//         fi12_lab.SetXYZ(0 , 0 , 0);
//         fi13_lab.SetXYZ(0 , 0 , 0); 

//         fib11_ang_offset.SetXYZ(0, par5, 0);
//         fib11_pos_offset.SetXYZ(0, 0 ,0);  

//         fib12_ang_offset.SetXYZ(0, par9, 0);
//         fib12_pos_offset.SetXYZ(0, 0, 0);  

//         fib10_ang_offset.SetXYZ(0, par6, 0);
//         fib10_pos_offset.SetXYZ(0, 0, 0);
        
//         psp2_ang_offset.SetXYZ(0 , 0 , 0);
//         psp3_ang_offset.SetXYZ(0 , 0 , 0);
//         fib13_ang_offset.SetXYZ(0 , 0 , 0);
//         psp2_pos_offset.SetXYZ(0 , 0 , 0);
//         psp3_pos_offset.SetXYZ(0 , 0 , 0);
//         fib13_pos_offset.SetXYZ(0 , 0 , 0);    
   
    
//         psp2_position.SetXYZ(0,0,targ_to_psp2);
//         psp3_position.SetXYZ(0,0,targ_to_psp3);    
//         fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
//         fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
//         fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl));
//         fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl));      
            
//         fib10_angles.SetXYZ(0, AngleAl, 0);
//         fib11_angles.SetXYZ(0, AngleAl, 0);
//         fib12_angles.SetXYZ(0, AngleAl, 0);
//         fib13_angles.SetXYZ(0, AngleAl, 0); 
        
//         TransformPoint(fi11_lab_223, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_223, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//         TransformPoint(fi11_lab_227, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_227, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//         TransformPoint(fi11_lab_225, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_225, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//         TransformPoint(fi11_lab_229, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_229, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
//         TransformPoint(fi11_lab_212, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_212, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//         TransformPoint(fi11_lab_220, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_220, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
//         TransformPoint(fi11_lab_201, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_201, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
//         TransformPoint(fi11_lab_202, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
//         TransformPoint(fi12_lab_202, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
//         TransformPoint(fi10_lab_202, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//         TransformPoint(fi10_lab_201, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//         TransformPoint(fi10_lab_212, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//         TransformPoint(fi10_lab_220, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//         TransformPoint(fi10_lab_223, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//         TransformPoint(fi10_lab_225, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//         TransformPoint(fi10_lab_227, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
//         TransformPoint(fi10_lab_229, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        
//         TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
//         TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
//         TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset ); 

//     double delta = 0.0;
//     //delta += fabs(Reference_PoQ_227 - CalculatePoQ(fi11_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_223 - CalculatePoQ(fi11_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_225 - CalculatePoQ(fi11_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_229 - CalculatePoQ(fi11_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_201 - CalculatePoQ(fi11_lab_201, fi12_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_202 - CalculatePoQ(fi11_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ));
//     //delta += fabs(Reference_PoQ_212 - CalculatePoQ(fi11_lab_212, fi12_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ));    
//     //delta += fabs(Reference_PoQ_227 - CalculatePoQ(fi11_lab_227, fi10_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_223 - CalculatePoQ(fi11_lab_223, fi10_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_225 - CalculatePoQ(fi11_lab_225, fi10_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_229 - CalculatePoQ(fi11_lab_229, fi10_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ10)); 
//     //delta += fabs(Reference_PoQ_201 - CalculatePoQ(fi11_lab_201, fi10_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_202 - CalculatePoQ(fi11_lab_202, fi10_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ10));
//     //delta += fabs(Reference_PoQ_212 - CalculatePoQ(fi11_lab_212, fi10_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ10));   
    
//     //delta += fabs(CalculatePoQ(fi11_lab_227, fi10_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ));
//     //delta += fabs(CalculatePoQ(fi11_lab_225, fi10_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ));
//     //delta += fabs(CalculatePoQ(fi11_lab_223, fi10_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ));        
//     //delta += fabs(CalculatePoQ(fi11_lab_229, fi10_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ));
//     //delta += fabs(CalculatePoQ(fi11_lab_202, fi10_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ));
//     //delta += fabs(CalculatePoQ(fi11_lab_201, fi10_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_201, fi12_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ));

//    //return fval = delta;
   
//    double area_227 = (0.5 * fabs( fi11_lab_227.X() * (fi10_lab_227.Z() - fi12_lab_227.Z()) + fi10_lab_227.X() * (fi12_lab_227.Z() - fi11_lab_227.Z()) + fi12_lab_227.X() * (fi11_lab_227.Z() - fi10_lab_227.Z()) ));
//    std::cout << "Area 227 = " << area_227 << std::endl;
//    double area_229 = (0.5 * fabs( fi11_lab_229.X() * (fi10_lab_229.Z() - fi12_lab_229.Z()) + fi10_lab_229.X() * (fi12_lab_229.Z() - fi11_lab_229.Z()) + fi12_lab_229.X() * (fi11_lab_229.Z() - fi10_lab_229.Z()) ));
//    std::cout << "Area 229 = " << area_229 << std::endl;
//    double area_202 = (0.5 * fabs( fi11_lab_202.X() * (fi10_lab_202.Z() - fi12_lab_202.Z()) + fi10_lab_202.X() * (fi12_lab_202.Z() - fi11_lab_202.Z()) + fi12_lab_202.X() * (fi11_lab_202.Z() - fi10_lab_202.Z()) ));
//    std::cout << "Area 202 = " << area_202 << std::endl;
//    double area_223 = (0.5 * fabs( fi11_lab_223.X() * (fi10_lab_223.Z() - fi12_lab_223.Z()) + fi10_lab_223.X() * (fi12_lab_223.Z() - fi11_lab_223.Z()) + fi12_lab_223.X() * (fi11_lab_223.Z() - fi10_lab_223.Z()) ));
//    std::cout << "Area 223 = " << area_223 << std::endl;
//    double area_225 = (0.5 * fabs( fi11_lab_225.X() * (fi10_lab_225.Z() - fi12_lab_225.Z()) + fi10_lab_225.X() * (fi12_lab_225.Z() - fi11_lab_225.Z()) + fi12_lab_225.X() * (fi11_lab_225.Z() - fi10_lab_225.Z()) ) );
//    std::cout << "Area 225 = " << area_225 << std::endl;
   
//     double errors[9];
//     errors[0] = fabs(CalculatePoQ(fi11_lab_227, fi10_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ));
//     errors[1] = fabs(CalculatePoQ(fi11_lab_229, fi10_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ));
//     errors[2] = fabs(CalculatePoQ(fi11_lab_202, fi10_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ));
//     errors[3] = fabs(CalculatePoQ(fi11_lab_225, fi10_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ));
//     errors[4] = fabs(CalculatePoQ(fi11_lab_223, fi10_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ));
//     errors[5] = fabs(Reference_PoQ_227 - CalculatePoQ(fi11_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ)) /5.;
//     errors[6] = fabs(Reference_PoQ_229 - CalculatePoQ(fi11_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ));
//     errors[7] = fabs(Reference_PoQ_202 - CalculatePoQ(fi11_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ)); 
//     errors[8] = fabs(Reference_PoQ_223 - CalculatePoQ(fi11_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ)) /6.; 

//     double max_error = std::max({errors[2],errors[7], area_225/10000., area_223/10000., area_229/10000., area_227/10000., area_202/10000.}); //, errors[5], errors[6] , errors[7], errors[8]

//     return fval = max_error;
// }

double MinimizeDelta4(const double* par, double Reference_PoQ, TString MDF_PoQ_file) {

    // Load the MDF
    MDF_PoQ = new R3BMDFWrapper(MDF_PoQ_file);

    // Paramaters loaded
    double fval;
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

    // Current of the run
    double Actual_Current = 1490.0;
    double Reference_Current = 1490.0;

    // Geometries
    const Double_t Target_to_GLAD_flange  = 127.9;                                         // <-------------------
    const Double_t targ_to_TP   = 293.1;                                                   // <-------------------
    const Double_t TP_to_GLAD_exit_flange = 385.0 - 125.0;                                 // <-------------------
    const Double_t Target_to_RP   = Target_to_GLAD_flange +  54.05 - 0.55580;              // <-------------------
    const Double_t targ_to_mwpc1 = 95.7;                              // target to MWPC1                 // <-------------------
    const Double_t tofd_to_TP    = TP_to_GLAD_exit_flange + 661.8;    // turning point to fib33           // <-------------------
    
    // Distance of fibers to TP
    const Double_t fi_32_to_TP   = 385.0;                              // turning point to fib32           // <-------------------
    const Double_t fi_33_to_TP   = TP_to_GLAD_exit_flange + 225.0;     // turning point to fib33           // <-------------------
    const Double_t fi_31_to_TP   = TP_to_GLAD_exit_flange + 242.5;     // turning point to fib31           // <-------------------

    // Double_t fi10_to_TP   = 255. + 201.3 + par1;   //(TP to vac.cham. internal wall) + (wall to det)
    // Double_t fi11_to_TP   = 255. + 65.3 ;          //(TP to vac.cham. internal wall) + (wall to det)
    // Double_t fi12_to_TP   = 1301.17 + 18. + par4;  //(TP to exit flange of vac. pipe) + (flange to det)
    // double Actual_Current_227 = 2166.0;//*par7; //Run227
    // double Actual_Current_223 = 2527.0;//*par8; //Run223 
    // double Actual_Current_225 = 2413.0;//*par13; //Run225
    // double Actual_Current_229 = 2641.0;//*par9; //Run229
    // double Actual_Current_212 = 2071.0;//*par10; //Run207
    // double Actual_Current_220 = 2299.0;//*par11; //Run219
    // double Actual_Current_201 = 1520.0;//*par14; //Run201
    // double Actual_Current_202 = 1691.0;//*par12; //Run202
    //Reference_PoQ_229 = 3.4560406+((0.653079312*par8));
    //Reference_PoQ_202 = 2.2147714+par8;
    //Reference_PoQ_223 = 3.324256 +(par8*0.214490469);
    //Reference_PoQ_201 = 1.996527 +par8;

    Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
    
    // fi11_lab_227.SetXYZ(-1.*(-4.60759 ) , 0 , 0 );
    // fi12_lab_227.SetXYZ(-1.*(-16.8725 + par5) , 0 , 0 );
    // fi10_lab_227.SetXYZ(-1.*(-5.20989 + par2) , 0 , 0 );  

    // fi11_lab_223.SetXYZ(-1.*(-3.20910 ) , 0 , 0 );
    // fi12_lab_223.SetXYZ(-1.*(-11.2588 + par5) , 0 , 0 );         
    // fi10_lab_223.SetXYZ(-1.*(-3.20786 + par2) , 0 , 0 ); 

    // fi11_lab_225.SetXYZ(-1.*(-3.64903 ) , 0 , 0 );
    // fi12_lab_225.SetXYZ(-1.*(-13.0462 + par5 ) , 0 , 0 );
    // fi10_lab_225.SetXYZ(-1.*(-3.84497 + par2) , 0 , 0 ); 
    
    // fi11_lab_229.SetXYZ(-1.*(-2.73996 ) , 0 , 0 );
    // fi12_lab_229.SetXYZ(-1.*(-9.26699 + par5 ) , 0 , 0 );
    // fi10_lab_229.SetXYZ(-1.*(-2.53834 + par2) , 0 , 0 ); 
    
    // fi11_lab_212.SetXYZ(-1.*(-2.83415 ) , 0 , 0 );
    // fi12_lab_212.SetXYZ(-1.*(-9.74823 + par5) , 0 , 0 ); 
    // fi10_lab_212.SetXYZ(-1.*(-2.68217 + par2) , 0 , 0 );  
                    
    // fi11_lab_220.SetXYZ(-1.*(-2.83545 ) , 0 , 0 );
    // fi12_lab_220.SetXYZ(-1.*(-10.0403 + par5) , 0 , 0 );
    // fi10_lab_220.SetXYZ(-1.*(-2.71734 + par2) , 0 , 0 );
    
    // fi11_lab_201.SetXYZ(-1.*(-2.93167 ) , 0 , 0 );
    // fi12_lab_201.SetXYZ(-1.*(-10.3565 + par5 ) , 0 , 0 );
    // fi10_lab_201.SetXYZ(-1.*(-2.85664 + par2 ) , 0 , 0 );
    
    // fi11_lab_202.SetXYZ(-1.*(-2.80716 ) , 0 , 0 );
    // fi12_lab_202.SetXYZ(-1.*(-9.61590 + par5) , 0 , 0 );                                         
    // fi10_lab_202.SetXYZ(-1.*(-2.63749 + par2) , 0 , 0 );
    
    // psp2_lab.SetXYZ(-0.5+par7 , 0 , 0);
    // psp3_lab.SetXYZ(-0.5+par9 , 0 , 0);
    // fi11_lab.SetXYZ(0 , 0 , 0);
    // fi12_lab.SetXYZ(0 , 0 , 0);
    // fi13_lab.SetXYZ(0 , 0 , 0); 

    // fib11_ang_offset.SetXYZ(0, 0, 0);
    // fib11_pos_offset.SetXYZ(0, 0 ,0);  
    // fib12_ang_offset.SetXYZ(0, par6, 0);
    // fib12_pos_offset.SetXYZ(0, 0, 0);  
    // fib10_ang_offset.SetXYZ(0, par3, 0);
    // fib10_pos_offset.SetXYZ(0, 0, 0);          
    // psp2_ang_offset.SetXYZ(0 , 0 , 0);
    // psp3_ang_offset.SetXYZ(0 , 0 , 0);
    // fib13_ang_offset.SetXYZ(0 , 0 , 0);
    // psp2_pos_offset.SetXYZ(0 , 0 , 0);
    // psp3_pos_offset.SetXYZ(0 , 0 , 0);
    // fib13_pos_offset.SetXYZ(0 , 0 , 0);  
    
    // Set the positions of central values of the planes of the detectors. Ideal case and orientations
    mwpc1_position.SetXYZ(0,0,targ_to_mwpc1);
    fi32_position.SetXYZ(fi_32_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_32_to_TP * TMath::Cos(AngleAl));
    fi33_position.SetXYZ(fi_33_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_33_to_TP * TMath::Cos(AngleAl));
    fi31_position.SetXYZ(fi_31_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_31_to_TP * TMath::Cos(AngleAl)); 
    
    mwpc1_angles.SetXYZ(0, 0, 0);
    fi32_angles.SetXYZ(0, AngleAl, 0);
    fi33_angles.SetXYZ(0, AngleAl, 0);
    fi31_angles.SetXYZ(0, AngleAl, 0);

    // Set the variations of positions and angles as the parameters
    mwpc1_pos_offset.SetXYZ(par1, 0,  par2);
    fi32_pos_offset.SetXYZ( par3  0,  par4);
    fi33_pos_offset.SetXYZ( par7, 0,  par8);
    fi31_pos_offset.SetXYZ(par11, 0, par12);

    mwpc1_ang_offset.SetXYZ(0, 0, 0);
    fi32_ang_offset.SetXYZ(0,  par5,  par6 );
    fi33_ang_offset.SetXYZ(0,  par9,  par10);
    fi31_ang_offset.SetXYZ(0, par13,  par14);

    // Set the positions of hits 
    mwpc1_lab.SetXYZ(0, 0, 0);
    fi32_lab.SetXYZ(0, 0, 0);
    fi33_lab.SetXYZ(0, 0, 0);
    fi31_lab.SetXYZ(0, 0, 0);

    // Up to this point the x and z positions of the fibers and mwpc are being placed in the lab frame 
    TransformPoint(mwpc1_lab, mwpc1_angles + mwpc1_ang_offset, mwpc1_position + mwpc1_pos_offset);
    TransformPoint( fi32_lab,   fi32_angles + fi32_ang_offset,   fi32_position + fi32_pos_offset);
    TransformPoint( fi33_lab,   fi33_angles + fi33_ang_offset,   fi33_position + fi33_pos_offset);
    TransformPoint( fi31_lab,   fi31_angles + fi31_ang_offset,   fi31_position + fi31_pos_offset);

    double delta = 0.0;

    // psp2_position.SetXYZ(0,0,targ_to_psp2+par8);
    // psp3_position.SetXYZ(0,0,targ_to_psp3+par10);
    // fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
    // fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
    // fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl));
    // fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl));
        
    // fib10_angles.SetXYZ(0, AngleAl, 0);
    // fib11_angles.SetXYZ(0, AngleAl, 0);
    // fib12_angles.SetXYZ(0, AngleAl, 0);
    // fib13_angles.SetXYZ(0, AngleAl, 0); 
    
    // TransformPoint(fi11_lab_223, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_223, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
    // TransformPoint(fi11_lab_227, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_227, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
    // TransformPoint(fi11_lab_225, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_225, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
    // TransformPoint(fi11_lab_229, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_229, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
    // TransformPoint(fi11_lab_212, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_212, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
    // TransformPoint(fi11_lab_220, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_220, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
    // TransformPoint(fi11_lab_201, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_201, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
    // TransformPoint(fi11_lab_202, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
    // TransformPoint(fi12_lab_202, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset ); 
    // TransformPoint(fi10_lab_202, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    // TransformPoint(fi10_lab_201, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    // TransformPoint(fi10_lab_212, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    // TransformPoint(fi10_lab_220, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    // TransformPoint(fi10_lab_223, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    // TransformPoint(fi10_lab_225, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    // TransformPoint(fi10_lab_227, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    // TransformPoint(fi10_lab_229, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
    
    // TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
    // TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
    // TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset ); 






    delta = fabs(Reference_PoQ - CalculatePoQ(mwpc1_lab, NEED TO COMBINE THE FIBERS HERE , fi32_lab, Actual_Current, Reference_Current, MDF_PoQ));
    

    //delta += fabs(Reference_PoQ_223 - CalculatePoQ(fi11_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ));
    //delta += fabs(Reference_PoQ_225 - CalculatePoQ(fi11_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ));
    //delta += fabs(Reference_PoQ_229 - CalculatePoQ(fi11_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ));
    //delta += fabs(Reference_PoQ_201 - CalculatePoQ(fi11_lab_201, fi12_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ));
    //delta += fabs(Reference_PoQ_202 - CalculatePoQ(fi11_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ));
    //delta += fabs(Reference_PoQ_212 - CalculatePoQ(fi11_lab_212, fi12_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ));    
    //delta += fabs(Reference_PoQ_227 - CalculatePoQ(fi11_lab_227, fi10_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ10));
    //delta += fabs(Reference_PoQ_223 - CalculatePoQ(fi11_lab_223, fi10_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ10));
    //delta += fabs(Reference_PoQ_225 - CalculatePoQ(fi11_lab_225, fi10_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ10));
    //delta += fabs(Reference_PoQ_229 - CalculatePoQ(fi11_lab_229, fi10_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ10)); 
    //delta += fabs(Reference_PoQ_201 - CalculatePoQ(fi11_lab_201, fi10_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ10));
    //delta += fabs(Reference_PoQ_202 - CalculatePoQ(fi11_lab_202, fi10_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ10));
    //delta += fabs(Reference_PoQ_212 - CalculatePoQ(fi11_lab_212, fi10_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ10));   
    
    //delta += fabs(CalculatePoQ(fi11_lab_227, fi10_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ));
    //delta += fabs(CalculatePoQ(fi11_lab_225, fi10_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ));
    //delta += fabs(CalculatePoQ(fi11_lab_223, fi10_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ));        
    //delta += fabs(CalculatePoQ(fi11_lab_229, fi10_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ));
    //delta += fabs(CalculatePoQ(fi11_lab_202, fi10_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ));
    //delta += fabs(CalculatePoQ(fi11_lab_201, fi10_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ10) - CalculatePoQ(fi11_lab_201, fi12_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ));

   //return fval = delta;
   
    // double errors[9];

    // errors[0] = fabs(Reference_PoQ_227 - CalculatePoQ1012(fi10_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ1012));
    // errors[1] = fabs(Reference_PoQ_229 - CalculatePoQ1012(fi10_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ1012));
    // errors[2] = fabs(Reference_PoQ_202 - CalculatePoQ1012(fi10_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ1012)); 
    // errors[3] = fabs(Reference_PoQ_223 - CalculatePoQ1012(fi10_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ1012)); 
    // errors[4] = fabs(Reference_PoQ_225 - CalculatePoQ1012(fi10_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ1012)); 
    // errors[5] = fabs(Reference_PoQ_201 - CalculatePoQ1012(fi10_lab_201, fi12_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ1012)); 
    // errors[6] = fabs(Reference_PoQ_212 - CalculatePoQ1012(fi10_lab_212, fi12_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ1012)); 
    // errors[7] = fabs(Reference_PoQ_220 - CalculatePoQ1012(fi10_lab_220, fi12_lab_220, psp2_lab, psp3_lab, fi13_lab, Actual_Current_220, Reference_Current, MDF_PoQ1012));    

    /*errors[0] = fabs(CalculatePoQ1012(fi10_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ1012) - CalculatePoQ1012(fi10_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ1012)) - fabs(Reference_PoQ_229-Reference_PoQ_227);
    errors[1] = fabs(CalculatePoQ1012(fi10_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ1012) - CalculatePoQ1012(fi10_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ1012)) - fabs(Reference_PoQ_202-Reference_PoQ_227);
    errors[2] = fabs(CalculatePoQ1012(fi10_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ1012) - CalculatePoQ1012(fi10_lab_220, fi12_lab_220, psp2_lab, psp3_lab, fi13_lab, Actual_Current_220, Reference_Current, MDF_PoQ1012)) - fabs(Reference_PoQ_229-Reference_PoQ_220);
    errors[3] = fabs(Reference_PoQ_227 - CalculatePoQ1012(fi10_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ1012));*/

    // double max_error = std::max({errors[2],errors[1],errors[4],errors[3]}); //, errors[5], errors[6] , errors[7], errors[8]
    // std::cout << "Error = " << max_error << std::endl;
    // return fval = max_error;

    return delta;
}

























double MinimizeGaus(const double* par) {

    double d_PSP2X, d_PSP2Y, d_PSP3X, d_PSP3Y, d_Fib10X, d_Fib11X, d_Fib12X, d_Fib13Y, d_TOF;
    int i_TOFD_Bar;

    // Branches verbinden

    tree->SetBranchAddress("PSP2X", &d_PSP2X);
    tree->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree->SetBranchAddress("PSP3X", &d_PSP3X);
    tree->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree->SetBranchAddress("Fib10X", &d_Fib10X);
    tree->SetBranchAddress("Fib11X", &d_Fib11X);
    tree->SetBranchAddress("Fib12X", &d_Fib12X);
    tree->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree->SetBranchAddress("TOF", &d_TOF);
    tree->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);

    double par1 = par[0]; double par2 = par[1];
    double par3 = par[2]; double par4 = par[3];
    double par5 = par[4]; double par6 = par[5];
    double par7 = par[6]; double par8 = par[7];
    double par9 = par[8]; double par10 = par[9];


    TH1D* dist;
    dist = new TH1D("dist","dist",1500,2.,3.5);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        Double_t fi10_to_TP   = 255. + 201.3 + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2166.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18. ) * TMath::Pi() / 180.;
        psp2_lab.SetXYZ(d_PSP2X-0.5+par8, d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X-0.5+par10, d_PSP3Y , 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9+ par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2+par7);
        psp3_position.SetXYZ(0,0,targ_to_psp3+par9);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
        
        double pq_227 = CalculatePoQ1012(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist->Fill(pq_227);
        }
    
    std::vector<double> expected_peaks = {2.892, 2.843, 2.951};
    double total_penalty = 0.0;

    for (double peak_pos : expected_peaks) {
        TF1 gauss("gauss", "gaus", peak_pos - 0.02, peak_pos + 0.02);
        int fit_status = dist->Fit(&gauss, "QNR");

        if (fit_status != 0) {
            std::cerr << "Fit failed at position: " << peak_pos << std::endl;
            total_penalty += 1000.0;
            continue;
        }
        double chi2_ndf = gauss.GetChisquare() / std::max(1.0, (double)gauss.GetNDF());
        double mean_error = gauss.GetParError(1);
        double amp = gauss.GetParameter(0);
        double amp_err = gauss.GetParError(0);
        
        if (gauss.GetChisquare() <= 0 || gauss.GetParameter(1) == 0) {
            total_penalty += 1000.0; // Große Strafe für fehlgeschlagenen Fit
            continue;
        }

        // Gewichtung kannst du anpassen je nach Stabilität
        total_penalty += std::abs(chi2_ndf - 1.0);
        total_penalty += mean_error * 10.0;
        if (amp_err > 0) total_penalty += 1.0 / (amp / amp_err);
    }

    delete dist;
    return total_penalty;

}

double MinimizeGaus2(const double* par) {

    double d_PSP2X, d_PSP2Y, d_PSP3X, d_PSP3Y, d_Fib10X, d_Fib11X, d_Fib12X, d_Fib13Y, d_TOF;
    int i_TOFD_Bar;

    // Branches verbinden

    tree->SetBranchAddress("PSP2X", &d_PSP2X);
    tree->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree->SetBranchAddress("PSP3X", &d_PSP3X);
    tree->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree->SetBranchAddress("Fib10X", &d_Fib10X);
    tree->SetBranchAddress("Fib11X", &d_Fib11X);
    tree->SetBranchAddress("Fib12X", &d_Fib12X);
    tree->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree->SetBranchAddress("TOF", &d_TOF);
    tree->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);
    
    tree229->SetBranchAddress("PSP2X", &d_PSP2X);
    tree229->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree229->SetBranchAddress("PSP3X", &d_PSP3X);
    tree229->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree229->SetBranchAddress("Fib10X", &d_Fib10X);
    tree229->SetBranchAddress("Fib11X", &d_Fib11X);
    tree229->SetBranchAddress("Fib12X", &d_Fib12X);
    tree229->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree229->SetBranchAddress("TOF", &d_TOF);
    tree229->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);
    
    tree223->SetBranchAddress("PSP2X", &d_PSP2X);
    tree223->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree223->SetBranchAddress("PSP3X", &d_PSP3X);
    tree223->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree223->SetBranchAddress("Fib10X", &d_Fib10X);
    tree223->SetBranchAddress("Fib11X", &d_Fib11X);
    tree223->SetBranchAddress("Fib12X", &d_Fib12X);
    tree223->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree223->SetBranchAddress("TOF", &d_TOF);
    tree223->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);    

    double par1 = par[0]; double par2 = par[1];
    double par3 = par[2]; double par4 = par[3];
    double par5 = par[4]; double par6 = par[5];
    double par7 = par[6]; double par8 = par[7];
    double par9 = par[8]; double par10 = par[9];


    TH1D* dist;
    dist = new TH1D("dist","dist",1500,2.,3.5);
    TH1D* angle227;
    angle227 = new TH1D("angle227", "angle227", 200,-20.,-15.);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        Double_t fi10_to_TP   = 255. + 201.3 + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2166.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18.05 ) * TMath::Pi() / 180.;
        psp2_lab.SetXYZ(d_PSP2X-0.5+par8, d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X-0.5+par10, d_PSP3Y , 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9+ par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , par7 , 0);
        psp3_ang_offset.SetXYZ(0 , par9 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
        
        double pq_227 = CalculatePoQ1012_11(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist->Fill(pq_227);
        double ang = atan(((fi12_lab.X()-fi10_lab.X())/(fi12_lab.Z()-fi10_lab.Z())))*(180/3.14159);
        angle227->Fill(ang);
        }
    
    TH1D* dist229;
    dist229 = new TH1D("dist229","dist229",1800,2.,3.8);
    TH1D* angle229;
    angle229 = new TH1D("angle229", "angle229", 200,-20.,-15.);
    Long64_t nentries229 = tree229->GetEntries();
    for (Long64_t i = 0; i < nentries229; ++i) {
        tree229->GetEntry(i);

        Double_t fi10_to_TP   = 255. + 201.3  + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2641.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18.05 ) * TMath::Pi() / 180.;
        psp2_lab.SetXYZ(d_PSP2X-0.5+par8, d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X-0.5+par10, d_PSP3Y , 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9+ par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , par7 , 0);
        psp3_ang_offset.SetXYZ(0 , par9 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3);      
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
        
        double pq_227 = CalculatePoQ1012_11(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist229->Fill(pq_227);
        double ang = atan(((fi12_lab.X()-fi10_lab.X())/(fi12_lab.Z()-fi10_lab.Z())))*(180/3.14159);
        angle229->Fill(ang);
        }
    
    TH1D* dist223;
    dist223 = new TH1D("dist223","dist223",1800,2.,3.8);
    TH1D* angle223;
    angle223 = new TH1D("angle223", "angle223", 200,-20.,-15.);
    Long64_t nentries223 = tree223->GetEntries();
    for (Long64_t i = 0; i < nentries223; ++i) {
        tree223->GetEntry(i);

        Double_t fi10_to_TP   = 255. + 201.3 +  par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3  ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2527.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18.05 ) * TMath::Pi() / 180.;
        psp2_lab.SetXYZ(d_PSP2X-0.5+par8, d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X-0.5+par10, d_PSP3Y , 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9+ par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , par7 , 0);
        psp3_ang_offset.SetXYZ(0 , par9 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3);     
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
        
        double pq_227 = CalculatePoQ1012_11(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist223->Fill(pq_227);
        double ang = atan(((fi12_lab.X()-fi10_lab.X())/(fi12_lab.Z()-fi10_lab.Z())))*(180/3.14159);
        angle223->Fill(ang);
        }
    
    std::vector<double> expected_peaks = {2.886};
    double total_penalty = 0.0;

    double peak_pos = 2.891;
    double peak_pos229 = 3.4560406;
    double peak_pos223 = 3.3233172;
    //for (double peak_pos : expected_peaks) {
        TF1 gauss("gauss", "gaus", peak_pos - 0.03, peak_pos + 0.03);
        TF1 gauss229("gauss229", "gaus", peak_pos229 - 0.03, peak_pos229 + 0.03);
        TF1 gauss223("gauss223", "gaus", peak_pos223 - 0.03, peak_pos223 + 0.03);
        int fit_status = dist->Fit(&gauss, "QNR");
        int fit_status229 = dist229->Fit(&gauss229, "QNR");
        int fit_status223 = dist223->Fit(&gauss223, "QNR");
        TF1 gauss_ang227("gauss_ang227", "gaus", -18.5, -17.0);
        TF1 gauss_ang229("gauss_ang229", "gaus", -18.5, -17.0);
        TF1 gauss_ang223("gauss_ang223", "gaus", -18.5, -17.0);
        int fit_status_ang227 = angle227->Fit(&gauss_ang227, "QNR");
        int fit_status_ang229 = angle229->Fit(&gauss_ang229, "QNR");
        int fit_status_ang223 = angle223->Fit(&gauss_ang223, "QNR");

        /*if (fit_status != 0) {
            std::cerr << "Fit failed at position: " << peak_pos << std::endl;
            total_penalty += 1000.0;
            continue;
        }*/
        double chi2_ndf = gauss.GetChisquare() / std::max(1.0, (double)gauss.GetNDF());
        double mean_error = gauss.GetParError(1);
        double amp = gauss.GetParameter(0);
        double amp229 = gauss229.GetParameter(0);
        double amp223 = gauss223.GetParameter(0);
        double amp_err = gauss.GetParError(0);
        double sigma = gauss.GetParameter(2);
        double sigma229 = gauss229.GetParameter(2);
        double sigma223 = gauss223.GetParameter(2);
        double mean = gauss.GetParameter(1);
        double mean229 = gauss229.GetParameter(1);
        double mean223 = gauss223.GetParameter(1);
        
        double mean_ang227 = gauss_ang227.GetParameter(1);
        double mean_ang223 = gauss_ang223.GetParameter(1);
        double mean_ang229 = gauss_ang229.GetParameter(1);
        
        //double max_error = std::max({sigma,sigma229,sigma223}); //,fabs(mean-peak_pos),fabs(mean229-peak_pos229),fabs(mean223-peak_pos223)
        double max_error = std::max({fabs(mean_ang223+17.871),fabs(mean_ang229+18.068),fabs(mean_ang227+17.657)});
        if (fit_status != 0) {
            std::cerr << "Fit failed at position for 227: " << peak_pos << std::endl;
            max_error += 1000.0;
            
        }
        if (fit_status229 != 0) {
            std::cerr << "Fit failed at position for 229: " << peak_pos229 << std::endl;
            max_error += 1000.0;
            
        }
        if (fit_status223 != 0) {
            std::cerr << "Fit failed at position for 223: " << peak_pos223 << std::endl;
            max_error += 1000.0;
            
        }
        if(amp<200){
            std::cout << "Small amplitude 227" << std::endl;
            max_error += 2000.0;
        } 
        if(amp229<200){
            std::cout << "Small amplitude 229" << std::endl;
            max_error += 2000.0;
        }
        if(amp223<200){
            std::cout << "Small amplitude 223" << std::endl;
            max_error += 2000.0;
        }
        std::cout << "Error = " << max_error << " sigma227 = " << sigma << " - sigma223 = " << sigma223 << " - sigma229 = " << sigma229 << std::endl;
        std::cout << "Angle 227 = " << mean_ang227 << " - Angle 223 = " << mean_ang223 << " - Angle 229 = " << mean_ang229 << std::endl;
        std::cout << "PoQ 227 = " << mean << " - PoQ 223 = " << mean223 << " - PoQ 229 = " << mean229 << std::endl;
        delete dist;
        delete dist229;
        delete dist223;
        delete angle227;
        delete angle223;
        delete angle229;
        return max_error;
        
        /*if (gauss.GetChisquare() <= 0 || gauss.GetParameter(1) == 0) {
            total_penalty += 1000.0; // Große Strafe für fehlgeschlagenen Fit
            continue;
        }*/

        // Gewichtung kannst du anpassen je nach Stabilität
        //total_penalty += std::abs(chi2_ndf - 1.0);
        //total_penalty += mean_error * 10.0;
        //if (amp_err > 0) total_penalty += 1.0 / (amp / amp_err);
        total_penalty += std::abs(sigma);
    //}

    
    //return total_penalty;

}

double MinimizeGaus09(const double* par) {

    double d_PSP2X, d_PSP2Y, d_PSP3X, d_PSP3Y, d_Fib10X, d_Fib11X, d_Fib12X, d_Fib13Y, d_TOF;
    int i_TOFD_Bar;

    // Branches verbinden

    tree->SetBranchAddress("PSP2X", &d_PSP2X);
    tree->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree->SetBranchAddress("PSP3X", &d_PSP3X);
    tree->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree->SetBranchAddress("Fib10X", &d_Fib10X);
    tree->SetBranchAddress("Fib11X", &d_Fib11X);
    tree->SetBranchAddress("Fib12X", &d_Fib12X);
    tree->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree->SetBranchAddress("TOF", &d_TOF);
    tree->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);
    
    tree229->SetBranchAddress("PSP2X", &d_PSP2X);
    tree229->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree229->SetBranchAddress("PSP3X", &d_PSP3X);
    tree229->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree229->SetBranchAddress("Fib10X", &d_Fib10X);
    tree229->SetBranchAddress("Fib11X", &d_Fib11X);
    tree229->SetBranchAddress("Fib12X", &d_Fib12X);
    tree229->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree229->SetBranchAddress("TOF", &d_TOF);
    tree229->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);
    
    tree223->SetBranchAddress("PSP2X", &d_PSP2X);
    tree223->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree223->SetBranchAddress("PSP3X", &d_PSP3X);
    tree223->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree223->SetBranchAddress("Fib10X", &d_Fib10X);
    tree223->SetBranchAddress("Fib11X", &d_Fib11X);
    tree223->SetBranchAddress("Fib12X", &d_Fib12X);
    tree223->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree223->SetBranchAddress("TOF", &d_TOF);
    tree223->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);    

    double par1 = par[0]; double par2 = par[1];
    double par3 = par[2]; double par4 = par[3];
    double par5 = par[4]; double par6 = par[5];
    double par7 = par[6]; double par8 = par[7];
    double par9 = par[8];


    TH1D* dist;
    dist = new TH1D("dist","dist",1500,2.,3.5);
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);

        Double_t fi10_to_TP   = 255. + 201.3 + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2166.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18.05 ) * TMath::Pi() / 180.;
        psp2_lab.SetXYZ(d_PSP2X-0.5, d_PSP2Y , 0);
        psp3_lab.SetXYZ((d_PSP3X-0.5+par8)*par9, d_PSP3Y , 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9+ par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3+par7);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
        
        double pq_227 = CalculatePoQ1012_11(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist->Fill(pq_227);
        }
    
    TH1D* dist229;
    dist229 = new TH1D("dist229","dist229",1800,2.,3.8);
    Long64_t nentries229 = tree229->GetEntries();
    for (Long64_t i = 0; i < nentries229; ++i) {
        tree229->GetEntry(i);

        Double_t fi10_to_TP   = 255. + 201.3  + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2641.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18.05 ) * TMath::Pi() / 180.;
        psp2_lab.SetXYZ(d_PSP2X-0.5, d_PSP2Y , 0);
        psp3_lab.SetXYZ((d_PSP3X-0.5+par8)*par9, d_PSP3Y , 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9+ par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3+par7);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
        
        double pq_227 = CalculatePoQ1012_11(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist229->Fill(pq_227);
        }
    
    TH1D* dist223;
    dist223 = new TH1D("dist223","dist223",1800,2.,3.8);
    Long64_t nentries223 = tree223->GetEntries();
    for (Long64_t i = 0; i < nentries223; ++i) {
        tree223->GetEntry(i);

        Double_t fi10_to_TP   = 255. + 201.3 +  par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3  ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2527.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18.05 ) * TMath::Pi() / 180.;
        psp2_lab.SetXYZ(d_PSP2X-0.5, d_PSP2Y , 0);
        psp3_lab.SetXYZ((d_PSP3X-0.5+par8)*par9, d_PSP3Y , 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9+ par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3+par7);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
        
        double pq_227 = CalculatePoQ1012_11(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist223->Fill(pq_227);
        }
    
    std::vector<double> expected_peaks = {2.886};
    double total_penalty = 0.0;

    double peak_pos = 2.891;
    double peak_pos229 = 3.4560406;
    double peak_pos223 = 3.3233172;
    //for (double peak_pos : expected_peaks) {
        TF1 gauss("gauss", "gaus", peak_pos - 0.03, peak_pos + 0.03);
        TF1 gauss229("gauss229", "gaus", peak_pos229 - 0.03, peak_pos229 + 0.03);
        TF1 gauss223("gauss223", "gaus", peak_pos223 - 0.03, peak_pos223 + 0.03);
        int fit_status = dist->Fit(&gauss, "QNR");
        int fit_status229 = dist229->Fit(&gauss229, "QNR");
        int fit_status223 = dist223->Fit(&gauss223, "QNR");

        /*if (fit_status != 0) {
            std::cerr << "Fit failed at position: " << peak_pos << std::endl;
            total_penalty += 1000.0;
            continue;
        }*/
        double chi2_ndf = gauss.GetChisquare() / std::max(1.0, (double)gauss.GetNDF());
        double mean_error = gauss.GetParError(1);
        double amp = gauss.GetParameter(0);
        double amp229 = gauss229.GetParameter(0);
        double amp223 = gauss223.GetParameter(0);
        double amp_err = gauss.GetParError(0);
        double sigma = gauss.GetParameter(2);
        double sigma229 = gauss229.GetParameter(2);
        double sigma223 = gauss223.GetParameter(2);
        
        double max_error = std::max({sigma,sigma229,sigma223}); //, errors[5], errors[6] , errors[7], errors[8]
        if (fit_status != 0) {
            std::cerr << "Fit failed at position for 227: " << peak_pos << std::endl;
            max_error += 1000.0;
            
        }
        if (fit_status229 != 0) {
            std::cerr << "Fit failed at position for 229: " << peak_pos229 << std::endl;
            max_error += 1000.0;
            
        }
        if (fit_status223 != 0) {
            std::cerr << "Fit failed at position for 223: " << peak_pos223 << std::endl;
            max_error += 1000.0;
            
        }
        if(amp<200){
            std::cout << "Small amplitude 227" << std::endl;
            max_error += 2000.0;
        } 
        if(amp229<200){
            std::cout << "Small amplitude 229" << std::endl;
            max_error += 2000.0;
        }
        if(amp223<200){
            std::cout << "Small amplitude 223" << std::endl;
            max_error += 2000.0;
        }
        std::cout << "Error = " << max_error << std::endl;
        delete dist;
        delete dist229;
        delete dist223;
        return max_error;
        
        /*if (gauss.GetChisquare() <= 0 || gauss.GetParameter(1) == 0) {
            total_penalty += 1000.0; // Große Strafe für fehlgeschlagenen Fit
            continue;
        }*/

        // Gewichtung kannst du anpassen je nach Stabilität
        //total_penalty += std::abs(chi2_ndf - 1.0);
        //total_penalty += mean_error * 10.0;
        //if (amp_err > 0) total_penalty += 1.0 / (amp / amp_err);
        total_penalty += std::abs(sigma);
    //}

    
    //return total_penalty;

}

double MinimizeAngles(const double* par) {

    double d_PSP2X, d_PSP2Y, d_PSP3X, d_PSP3Y, d_Fib10X, d_Fib11X, d_Fib12X, d_Fib13Y, d_TOF;
    int i_TOFD_Bar;

    // Branches verbinden

    tree227->SetBranchAddress("PSP2X", &d_PSP2X);
    tree227->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree227->SetBranchAddress("PSP3X", &d_PSP3X);
    tree227->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree227->SetBranchAddress("Fib10X", &d_Fib10X);
    tree227->SetBranchAddress("Fib11X", &d_Fib11X);
    tree227->SetBranchAddress("Fib12X", &d_Fib12X);
    tree227->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree227->SetBranchAddress("TOF", &d_TOF);
    tree227->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);

    double par1 = par[0]; double par2 = par[1];
    double par3 = par[2]; double par4 = par[3];
    double par5 = par[4]; double par6 = par[5];
    double par7 = par[6]; double par8 = par[7];
    double par9 = par[8];


    TH1D* dist10;
    TH1D* dist12;
    TH1D* distDiff;
    dist10 = new TH1D("dist10","dist10",100,-1.2,1.2);
    dist12 = new TH1D("dist12","dist12",300,-1.2,1.2);
    distDiff = new TH1D("distDiff","distDiff",300,-1.2,1.2);
    Long64_t nentries227 = tree227->GetEntries();
    for (Long64_t i = 0; i < nentries227; ++i) {
        tree227->GetEntry(i);
    
        Double_t fi10_to_TP   = 255. + 201.3 + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2166.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18.0 ) * TMath::Pi() / 180.;        
        psp2_lab.SetXYZ(d_PSP2X-0.5 , d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X-0.5 , d_PSP3Y, 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X + par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
    
        
            
        double angle12 = atan(((fi12_lab.X()-fi10_lab.X())/(fi12_lab.Z()-fi10_lab.Z()))*(180/3.14159));
        double angle10 = atan(((fi10_lab.X()-fi11_lab.X())/(fi10_lab.Z()-fi11_lab.Z()))*(180/3.14159));
        
        dist10->Fill(angle10);
        dist12->Fill(angle12);
        //distDiff->Fill(angle10-angle12);

    
    }
    
    double total_penalty = 0.0;
    
    double peak_pos10, peak_pos12;

    
    TSpectrum s(1);  // Suche nach 1 Peak
    s.Search(dist10, 2, "", 0.05);  // σ = 2 bins, Schwelle = 5%

    Double_t* peaks = s.GetPositionX();
    if (s.GetNPeaks() > 0) {
    peak_pos10 = peaks[0];
    }
    s.Search(dist12, 2, "goff", 0.05);  // σ = 2 bins, Schwelle = 5%

    peaks = s.GetPositionX();
    if (s.GetNPeaks() > 0) {
    peak_pos12 = peaks[0];
    }
    
    TF1 gauss("gauss", "gaus", -0.04, 0.04);
    int fit_status = distDiff->Fit(&gauss, "QNR");
    double sigma = gauss.GetParameter(2);
    
    total_penalty += fabs(sigma-0.02);    
    total_penalty += fabs(peak_pos10-0.32);
    total_penalty += fabs(peak_pos12-0.32);
    
    delete dist12;
    delete dist10;
    delete distDiff;
    return total_penalty;    
        
        
}        

double MinimizeArea(const double* par) {

    double d_PSP2X, d_PSP2Y, d_PSP3X, d_PSP3Y, d_Fib10X, d_Fib11X, d_Fib12X, d_Fib13Y, d_TOF;
    int i_TOFD_Bar;

    // Branches verbinden

    tree->SetBranchAddress("PSP2X", &d_PSP2X);
    tree->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree->SetBranchAddress("PSP3X", &d_PSP3X);
    tree->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree->SetBranchAddress("Fib10X", &d_Fib10X);
    tree->SetBranchAddress("Fib11X", &d_Fib11X);
    tree->SetBranchAddress("Fib12X", &d_Fib12X);
    tree->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree->SetBranchAddress("TOF", &d_TOF);
    tree->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);
    
    tree229->SetBranchAddress("PSP2X", &d_PSP2X);
    tree229->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree229->SetBranchAddress("PSP3X", &d_PSP3X);
    tree229->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree229->SetBranchAddress("Fib10X", &d_Fib10X);
    tree229->SetBranchAddress("Fib11X", &d_Fib11X);
    tree229->SetBranchAddress("Fib12X", &d_Fib12X);
    tree229->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree229->SetBranchAddress("TOF", &d_TOF);
    tree229->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);

    double par1 = par[0]; double par2 = par[1];
    double par3 = par[2]; double par4 = par[3];
    double par5 = par[4]; double par6 = par[5];
    double par7 = par[6]; double par8 = par[7];
    double par9 = par[8];

    double sum_area227 = 0., sum_area229;
    double total_penalty = 0.0;
    double counter227 = 0., counter229 = 0.;

    TH1D* dist10;
    dist10 = new TH1D("dist10","dist10",1000,2.,3.0);
    Long64_t nentries227 = tree->GetEntries();
    for (Long64_t i = 0; i < nentries227; ++i) {
        tree->GetEntry(i);
    
        Double_t fi10_to_TP   = 255. + 201.3 + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 + par7;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2166.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18. ) * TMath::Pi() / 180.;        
        psp2_lab.SetXYZ(d_PSP2X , d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X , d_PSP3Y, 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X + par8, 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X + par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, par9, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
    
        if(d_PSP2X>-0.05 && d_PSP3X>-0.05 && d_PSP2X<0.05 && d_PSP3X<0.05 && d_PSP2Y>-1.4 && d_PSP3Y>-1.4 && d_PSP2Y<-1.0 && d_PSP3Y<-1.0){
            
        double angle12 = atan(((fi12_lab.X()-fi10_lab.X())/(fi12_lab.Z()-fi10_lab.Z()))*(180/3.14159));
        double angle10 = atan(((fi10_lab.X()-fi11_lab.X())/(fi10_lab.Z()-fi11_lab.Z()))*(180/3.14159));
        
        double area = 0.5 * fabs( fi11_lab.X() * (fi10_lab.Z() - fi12_lab.Z()) + fi10_lab.X() * (fi12_lab.Z() - fi11_lab.Z()) + fi12_lab.X() * (fi11_lab.Z() - fi10_lab.Z()) );
        counter227++;
        sum_area227 += area;
        
        double pq_227 = CalculatePoQ1012(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist10->Fill(pq_227);
        
        }
    
    }

    
    TH1D* dist229;
    dist229 = new TH1D("dist229","dist229",1000,2.7,3.7);
    Long64_t nentries229 = tree229->GetEntries();
    for (Long64_t i = 0; i < nentries229; ++i) {
        tree229->GetEntry(i);
    
        Double_t fi10_to_TP   = 255. + 201.3 + par1;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 65.3 + par7;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18. + par4;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2641.0;
        
        Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-18. ) * TMath::Pi() / 180.;        
        psp2_lab.SetXYZ(d_PSP2X , d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X , d_PSP3Y, 0);
        fi10_lab.SetXYZ(d_Fib10X + par2 , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X + par8, 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X + par5, 0 , 0 );
        fi13_lab.SetXYZ(0 , d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, par9, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, par6, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, par3, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2);
        psp3_position.SetXYZ(0,0,targ_to_psp3);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
    
        if(d_PSP2X>-0.05 && d_PSP3X>-0.05 && d_PSP2X<0.05 && d_PSP3X<0.05 && d_PSP2Y>-1.4 && d_PSP3Y>-1.4 && d_PSP2Y<-1.0 && d_PSP3Y<-1.0){
            
        double angle12 = atan(((fi12_lab.X()-fi10_lab.X())/(fi12_lab.Z()-fi10_lab.Z()))*(180/3.14159));
        double angle10 = atan(((fi10_lab.X()-fi11_lab.X())/(fi10_lab.Z()-fi11_lab.Z()))*(180/3.14159));
        
        double area = 0.5 * fabs( fi11_lab.X() * (fi10_lab.Z() - fi12_lab.Z()) + fi10_lab.X() * (fi12_lab.Z() - fi11_lab.Z()) + fi12_lab.X() * (fi11_lab.Z() - fi10_lab.Z()) );
        counter229++;
        sum_area229 += area;
        
        double pq_229 = CalculatePoQ1012(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist229->Fill(pq_229);
        
        }
    
    }
    
        double peak_pos10, peak_pos229;

    /*
    TSpectrum s(1);  // Suche nach 1 Peak
    
    s.Search(dist10, 2, "", 0.05);  // σ = 2 bins, Schwelle = 5%

    Double_t* peaks = s.GetPositionX();
    if (s.GetNPeaks() > 0) {
    peak_pos10 = peaks[0];
    }*/
    TF1 gauss("gauss", "gaus", 2.8 - 0.04, 2.9 + 0.04);
    int fit_status = dist10->Fit(&gauss, "QNR");
    double mean = gauss.GetParameter(1);
    double sigma = gauss.GetParameter(2);
    TF1 gauss229("gauss", "gaus", 3.3 - 0.04, 3.6 + 0.04);
    int fit_status229 = dist229->Fit(&gauss229, "QNR");
    double mean229 = gauss229.GetParameter(1);
    double sigma229 = gauss229.GetParameter(2);
    //std::cout << "mean227 = " << mean << " - Error = " << ((fabs(mean-Reference_PoQ_227))*10500.) << " mean229 = " << mean229 << " - Error = " << ((fabs(mean229-Reference_PoQ_229))*10500.) << " Error area = " << total_penalty/counter << std::endl;
    //std::cout << "mean227 = " << mean << " mean229 = " << mean229 << " - Error = " << ((fabs(mean229-Reference_PoQ_229)-fabs(mean-Reference_PoQ_227))*10500.) << " Error area 227 = " << sum_area227/counter227 << " Error area 229 = " << sum_area229/counter229 << std::endl;
    std::cout << "mean227 = " << mean << " sigma227 = " << sigma << " mean229 = " << mean229 << " sigma229 = " << sigma229 << " - Error = " << ((fabs(mean229-Reference_PoQ_229)-fabs(mean-Reference_PoQ_227))*10500.) << std::endl;
    
    //total_penalty=((sum_area229/counter229)+(sum_area227/counter227))/2.;
    //total_penalty=total_penalty+((fabs(mean-Reference_PoQ_227))*10500.);
    //total_penalty=total_penalty+((fabs(mean229-Reference_PoQ_229))*10500.);
    total_penalty=total_penalty+fabs((fabs(mean229-Reference_PoQ_229)-fabs(mean-Reference_PoQ_227))*10500.);
    //total_penalty=total_penalty/2.;
    delete dist10;
    delete dist229;
    return total_penalty;    
        
        
}        

double MinimizePSP(const double* par) {

    double d_PSP2X, d_PSP2Y, d_PSP3X, d_PSP3Y, d_Fib10X, d_Fib11X, d_Fib12X, d_Fib13Y, d_TOF;
    int i_TOFD_Bar;

    // Branches verbinden

    tree->SetBranchAddress("PSP2X", &d_PSP2X);
    tree->SetBranchAddress("PSP2Y", &d_PSP2Y);
    tree->SetBranchAddress("PSP3X", &d_PSP3X);
    tree->SetBranchAddress("PSP3Y", &d_PSP3Y);
    tree->SetBranchAddress("Fib10X", &d_Fib10X);
    tree->SetBranchAddress("Fib11X", &d_Fib11X);
    tree->SetBranchAddress("Fib12X", &d_Fib12X);
    tree->SetBranchAddress("Fib13Y", &d_Fib13Y);
    tree->SetBranchAddress("TOF", &d_TOF);
    tree->SetBranchAddress("TOFD_Bar", &i_TOFD_Bar);

    double par11 = par[0]; double par12 = par[1];

/*
double par1 = 0.977852;
double par2 = 0.914766;
double par3 = -0.0084096;
double par4 = 1.98389;
double par5 = 0.988075;
double par6 = 0.024969;
double par7 = 1.50106;
double par8 = 0.0864224;
double par9 = -0.0246014;
double par10 = 0.9801;*/

double par1 = 1.43553;
double par2 = 2.00796;
double par3 = -0.00981122;
double par4 = -1.87933;
double par5 = 1.34969;
double par6 = -0.000663831;
double par7 = 1.65553;
double par8 = 1.35315;
double par9 = 0.00180196;
double par10 = 0.9801;

    TH1D* dist;
    dist = new TH1D("dist","dist",1000,2.,3.0);

    Long64_t nentries227 = tree->GetEntries();
    for (Long64_t i = 0; i < nentries227; ++i) {
        tree->GetEntry(i);
    
        Double_t fi10_to_TP   = 255. + 191.5 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi11_to_TP   = 255. + 62.1 ;//(TP to vac.cham. internal wall) + (wall to det)
        Double_t fi12_to_TP   = 1301.17 + 18.;//(TP to exit flange of vac. pipe) + (flange to det)
        double Actual_Current = 2166.0;
        
        Double_t AngleAl = (-1.*(18.0) ) * TMath::Pi() / 180.;
        Double_t AngleAl12 = (-1.*(18.0) ) * TMath::Pi() / 180.;        
        psp2_lab.SetXYZ(d_PSP2X -0.5 , d_PSP2Y , 0);
        psp3_lab.SetXYZ(d_PSP3X -0.5, d_PSP3Y, 0);
        fi10_lab.SetXYZ(d_Fib10X , 0 , 0 ); 
        fi11_lab.SetXYZ(d_Fib11X , 0 , 0 );
        fi12_lab.SetXYZ(d_Fib12X -0.9, 0 , 0 );
        fi13_lab.SetXYZ(0 , -1.*d_Fib13Y , 0);
        
        fib11_ang_offset.SetXYZ(0, 0, 0);
        fib11_pos_offset.SetXYZ(0, 0 ,0);  

        fib12_ang_offset.SetXYZ(0, 0, 0);
        fib12_pos_offset.SetXYZ(0, 0, 0);  

        fib10_ang_offset.SetXYZ(0, 0, 0);
        fib10_pos_offset.SetXYZ(0, 0, 0);          
        
        psp2_ang_offset.SetXYZ(0 , 0 , 0);
        psp3_ang_offset.SetXYZ(0 , 0 , 0);
        fib13_ang_offset.SetXYZ(0 , 0 , 0);
        psp2_pos_offset.SetXYZ(0 , 0 , 0);
        psp3_pos_offset.SetXYZ(0 , 0 , 0);
        fib13_pos_offset.SetXYZ(0 , 0 , 0);
        
        psp2_position.SetXYZ(0,0,targ_to_psp2+par11);
        psp3_position.SetXYZ(0,0,targ_to_psp3+par12);    
        fib10_position.SetXYZ(fi10_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi10_to_TP * TMath::Cos(AngleAl));
        fib11_position.SetXYZ(fi11_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi11_to_TP * TMath::Cos(AngleAl));
        fib12_position.SetXYZ(fi12_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi12_to_TP * TMath::Cos(AngleAl12));
        fib13_position.SetXYZ(fi13_to_TP * TMath::Sin(AngleAl12), 0, targ_to_TP + fi13_to_TP * TMath::Cos(AngleAl12));      
            
        fib10_angles.SetXYZ(0, AngleAl, 0);
        fib11_angles.SetXYZ(0, AngleAl, 0);
        fib12_angles.SetXYZ(0, AngleAl12, 0);
        fib13_angles.SetXYZ(0, AngleAl12, 0);
        
        TransformPoint(psp2_lab, psp2_angles + psp2_ang_offset , psp2_position + psp2_pos_offset );
        TransformPoint(psp3_lab, psp3_angles + psp3_ang_offset , psp3_position + psp3_pos_offset );
        TransformPoint(fi10_lab, fib10_angles + fib10_ang_offset , fib10_position + fib10_pos_offset );
        TransformPoint(fi11_lab, fib11_angles + fib11_ang_offset , fib11_position + fib11_pos_offset );
        TransformPoint(fi12_lab, fib12_angles + fib12_ang_offset , fib12_position + fib12_pos_offset );
        TransformPoint(fi13_lab, fib13_angles + fib13_ang_offset , fib13_position + fib13_pos_offset );
    
        double angle12 = atan(((fi12_lab.X()-fi10_lab.X())/(fi12_lab.Z()-fi10_lab.Z()))*(180/3.14159));
        double angle10 = atan(((fi10_lab.X()-fi11_lab.X())/(fi10_lab.Z()-fi11_lab.Z()))*(180/3.14159));

        double pq_227 = CalculatePoQ1012(fi10_lab, fi12_lab, psp2_lab, psp3_lab, fi13_lab, Actual_Current, Reference_Current, MDF_PoQ1012);
        dist->Fill(pq_227);
    }
    
    double total_penalty = 0.0;
    
    double peak_pos;

    
    TF1 gauss("gauss", "gaus", 2.885, 2.895);
    int fit_status = dist->Fit(&gauss, "QNR");
    double sigma = gauss.GetParameter(2);
    
    total_penalty = sigma;    
    std::cout << "Sigma = " << sigma << " - Error = " << total_penalty << std::endl;
    
    delete dist;
    return total_penalty;    
        
        
}        
