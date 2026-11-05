#include "libs.hh"
#include "definitions.hh"
#include "R3BGladFieldMap.h"
#include "R3BTPropagator.h"
#include "R3BTrackingParticle.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include <iostream>
#include <limits>
#include <vector>

using namespace std;

R3BGladFieldMap* field = nullptr;
R3BTPropagator* propagator = nullptr;
const double TOL = 0.4;

TH1D* h_par1;
TH1D* h_par2;
TH1D* h_par3;

void TransformPoint(TVector3& point, TVector3 rot, TVector3 trans)
{
    TRotation r;
    r.RotateY(rot.Y());
    TVector3 localX(1, 0, 0);
    localX.SetTheta(r.ThetaX());
    localX.SetPhi(r.PhiX());
    r.Rotate(rot.X(), localX);
    TVector3 localZ(0, 0, 1);
    localZ.SetTheta(r.ThetaZ());
    localZ.SetPhi(r.PhiZ());
    r.Rotate(rot.Z(), localZ);
    point.Transform(r);
    point += trans;
}

R3BTrackingParticle Generate_Primary(double ekin_u, double par)
{
    const Double_t Ekin_Total = ekin_u * BEAM_MASS_GEV_C2 / AMU;
    Double_t P = sqrt(Ekin_Total * (Ekin_Total + 2 * BEAM_MASS_GEV_C2));
    TVector3 mom(0, 0, P);
    TVector3 pos(0. + par, 0, -20.);
    Double_t beta = 1. / TMath::Sqrt(1 + TMath::Power(BEAM_MASS_GEV_C2 / mom.Mag(), 2));
    return R3BTrackingParticle(BEAM_Z, pos.X(), pos.Y(), pos.Z(), mom.X(), mom.Y(), mom.Z(), beta, BEAM_MASS_GEV_C2);
}

bool Propagate(R3BTrackingParticle& particle, double z, double x, double z0, double x0, TVector3& result)
{
    TVector3 v1(x, 0, z), v2(x, 10., z), v3(x0, 0, z0);
    double tof = 0;
    if (!propagator->PropagateToPlaneRK_eloss(&particle, v1, v2, v3, tof, 0.1, false)) return false;
    //if (!propagator->PropagateToPlaneRK(&particle, v1, v2, v3)) return false;
    result = particle.GetPosition();
    return true;
}

double ObjectiveFunction(const double* p)
{
    double a11 = (-18. ) * TMath::DegToRad();
    double z11 = fi11_to_TP ;
    double a10 = (-18. ) * TMath::DegToRad();
    double z10 = fi10_to_TP ;
    double a12 = (-18. + p[0]) * TMath::DegToRad();
    double z12 = fi12_to_TP + p[1];

    TVector3 pos11(z11 * sin(a11), 0, Target_to_TP + z11 * cos(a11));
    TVector3 pos10(z10 * sin(a10), 0, Target_to_TP + z10 * cos(a10));
    TVector3 pos12(z12 * sin(a12), 0, Target_to_TP + z12 * cos(a12));

    TVector3 rot11(0, a11, 0), rot10(0, a10, 0), rot12(0, a12, 0);
    
    TVector3 fi11_lab_223(-1.*(-3.20910 ) , 0, 0);
    TVector3 fi10_lab_223(-1.*(-3.20786)  , 0, 0);
    TVector3 fi12_lab_223(-1.*(-11.2588)  +p[3], 0, 0);
    TVector3 fi11_lab_227(-1.*(-4.60759) , 0, 0);
    TVector3 fi10_lab_227(-1.*(-5.20989) , 0, 0);
    TVector3 fi12_lab_227(-1.*(-16.74) +p[3], 0, 0);
    TVector3 fi11_lab_229(-1.*(-2.73996) , 0, 0);
    TVector3 fi10_lab_229(-1.*(-2.53834) , 0, 0);
    TVector3 fi12_lab_229(-1.*(-9.38697) +p[3], 0, 0);
    TVector3 fi11_lab_202(-1.*(-2.80716) , 0, 0);
    TVector3 fi10_lab_202(-1.*(-2.63749) , 0, 0);
    TVector3 fi12_lab_202(-1.*(-9.61590) +p[3], 0, 0);
    TVector3 fi11_lab_201(-1.*(-3.13696) , 0, 0);
    TVector3 fi10_lab_201(-1.*(-3.09914) , 0, 0);
    TVector3 fi12_lab_201(-1.*(-10.3187) +p[3], 0, 0);
    
    TVector3 fi10_lab_zero(0., 0., 0.);
    TVector3 fi11_lab_zero(0., 0., 0.);
    TVector3 fi12_lab_zero(0., 0., 0.);    

    TVector3 t11_229 = fi11_lab_229, t10_229 = fi10_lab_229, t12_229 = fi12_lab_229;
    TVector3 t11_227 = fi11_lab_227, t10_227 = fi10_lab_227, t12_227 = fi12_lab_227;
    TVector3 t11_223 = fi11_lab_223, t10_223 = fi10_lab_223, t12_223 = fi12_lab_223;
    TVector3 t11_202 = fi11_lab_202, t10_202 = fi10_lab_202, t12_202 = fi12_lab_202;
    TVector3 t11_201 = fi11_lab_201, t10_201 = fi10_lab_201, t12_201 = fi12_lab_201;
    TVector3 t11_zero = fi11_lab_zero, t10_zero = fi10_lab_zero, t12_zero = fi12_lab_zero;

    TransformPoint(t11_229, rot11, pos11);
    TransformPoint(t10_229, rot10, pos10);
    TransformPoint(t12_229, rot12, pos12);
    TransformPoint(t11_227, rot11, pos11);
    TransformPoint(t10_227, rot10, pos10);
    TransformPoint(t12_227, rot12, pos12);
    TransformPoint(t11_223, rot11, pos11);
    TransformPoint(t10_223, rot10, pos10);
    TransformPoint(t12_223, rot12, pos12);
    TransformPoint(t11_202, rot11, pos11);
    TransformPoint(t10_202, rot10, pos10);
    TransformPoint(t12_202, rot12, pos12);
    TransformPoint(t11_201, rot11, pos11);
    TransformPoint(t10_201, rot10, pos10);
    TransformPoint(t12_201, rot12, pos12);
    
    TransformPoint(t11_zero, rot11, pos11);
    TransformPoint(t10_zero, rot10, pos10);
    TransformPoint(t12_zero, rot12, pos12);    

    R3BTrackingParticle ion229 = Generate_Primary(0.784853, p[2]);
    R3BTrackingParticle ion227 = Generate_Primary(0.592354, p[2]);
    R3BTrackingParticle ion223 = Generate_Primary(0.73866, p[2]);
    R3BTrackingParticle ion202 = Generate_Primary(0.380468, p[2]);
    R3BTrackingParticle ion201 = Generate_Primary(0.317634, p[2]);

    TVector3 hit11_229, hit10_229, hit12_229;
    TVector3 hit11_227, hit10_227, hit12_227;
    TVector3 hit11_223, hit10_223, hit12_223;
    TVector3 hit11_202, hit10_202, hit12_202;
    TVector3 hit11_201, hit10_201, hit12_201;

    field->SetScale(-2641.0 / 3583.81);
    //if (!Propagate(ion229, t11_229.Z(), t11_229.X(), t11_zero.Z(), t11_zero.X(), hit11_229)) return 1e6;
    //if (!Propagate(ion229, t10_229.Z(), t10_229.X(), t10_zero.Z(), t10_zero.X(), hit10_229)) return 1e6;
    if (!Propagate(ion229, t12_229.Z(), t12_229.X(), t12_zero.Z(), t12_zero.X(), hit12_229)) return 1e6;

    field->SetScale(-2166.0 / 3583.81);
    //if (!Propagate(ion227, t11_227.Z(), t11_227.X(), t11_zero.Z(), t11_zero.X(), hit11_227)) return 1e6;
    //if (!Propagate(ion227, t10_227.Z(), t10_227.X(), t10_zero.Z(), t10_zero.X(), hit10_227)) return 1e6;
    if (!Propagate(ion227, t12_227.Z(), t12_227.X(), t12_zero.Z(), t12_zero.X(), hit12_227)) return 1e6;

    field->SetScale(-2527.0 / 3583.81);
    //if (!Propagate(ion223, t11_223.Z(), t11_223.X(), t11_zero.Z(), t11_zero.X(), hit11_223)) return 1e6;
    //if (!Propagate(ion223, t10_223.Z(), t10_223.X(), t10_zero.Z(), t10_zero.X(), hit10_223)) return 1e6;
    if (!Propagate(ion223, t12_223.Z(), t12_223.X(), t12_zero.Z(), t12_zero.X(), hit12_223)) return 1e6;
        
    field->SetScale(-1691.0 / 3583.81);
    //if (!Propagate(ion202, t11_202.Z(), t11_202.X(), t11_zero.Z(), t11_zero.X(), hit11_202)) return 1e6;
    //if (!Propagate(ion202, t10_202.Z(), t10_202.X(), t10_zero.Z(), t10_zero.X(), hit10_202)) return 1e6;
    if (!Propagate(ion202, t12_202.Z(), t12_202.X(), t12_zero.Z(), t12_zero.X(), hit12_202)) return 1e6;
    
    field->SetScale(-1520.0 / 3583.81);
    //if (!Propagate(ion201, t11_201.Z(), t11_201.X(), t11_zero.Z(), t11_zero.X(), hit11_201)) return 1e6;
    //if (!Propagate(ion201, t10_201.Z(), t10_201.X(), t10_zero.Z(), t10_zero.X(), hit10_201)) return 1e6;
    if (!Propagate(ion201, t12_201.Z(), t12_201.X(), t12_zero.Z(), t12_zero.X(), hit12_201)) return 1e6;


    vector<pair<TVector3, TVector3>> checks = {
        //{hit10_229, t10_229}, {hit10_223, t10_223},
        //{hit10_227, t10_227}, 
        //{hit10_202, t10_202},
        //{hit10_201, t10_201},
        {hit12_229, t12_229}, {hit12_223, t12_223},
        {hit12_227, t12_227}, 
        {hit12_202, t12_202},        
        {hit12_201, t12_201},
        //{hit10_227, t10_227}, {hit12_227, t12_227}
        //{hit11_229, t11_229}, {hit10_229, t10_229}, {hit12_229, t12_229},
        //{hit11_227, t11_227}, {hit10_227, t10_227}, {hit12_227, t12_227}        
    };

    double total_error = 0;
    for (auto& [hit, target] : checks)
    {
        double dx = hit.X() - target.X();
        double dz = hit.Z() - target.Z();
        if (fabs(dx) > TOL || fabs(dz) > TOL)
            return 1e6;
        total_error += dx * dx + dz * dz;
    }

    return total_error;
}

int main()
{
    h_par1 = new TH1D("par1","par1",1000,-0.2,0.2);
    h_par2 = new TH1D("par2","par2",1000,-100,20);
    h_par3 = new TH1D("par3","par3",1000,-10,10);
    gRandom = new TRandom3(0);
    field = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-4to13_step10mm", "R");
    //field = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-3to5_step5mm", "R");
    //field = new R3BGladFieldMap("R3BGladMap", "A");
    field->SetPosition(0.7871, 1.75 - 1.526, Target_to_GLAD_flange + 54.05 - 0.55580);
    field->SetXAngle(-0.113);
    field->SetYAngle(-14.08);
    field->SetZAngle(-0.83);
    field->Init();

    propagator = new R3BTPropagator(field, kFALSE);

    // Parametergrenzen und Schrittweiten
    /*const int NPAR = 6;
    vector<string> names = {"angle11", "pos11", "angle10", "pos10", "angle12", "pos12", "xpos10", "xpos12"};
    vector<double> par_min = {-0.01, -20., -0.25, -1., -0.25, -70.};
    vector<double> par_max = {0.01, 0., 0.25, 1., 0.25, -30.};
    vector<double> step_size = {0.001, 0.1, 0.001, 0.1, 0.001, 0.1};*/

    const int NPAR = 4;
    vector<string> names = { "angle12", "zpos12", "offPSP", "xpos12"};
    vector<double> par_min = { -0.05 , -1.1, -0.8, -6.};
    vector<double> par_max = { 0.05 , -1., 0., 6.};
    vector<double> step_size = { 0.01, 0.01, 0.01, 0.01};

    std::cout << " min = " << par_min[0] << " , " << par_min[1] << " , " << par_min[2] << std::endl;
    std::cout << " max = " << par_max[0] << " , " << par_max[1] << " , " << par_max[2] << std::endl; 
    std::cout << " step_size = " << step_size[0] << " , " << step_size[1] << " , " << step_size[2] << std::endl;     

    int num = 0;
    while(num<3){
    bool success = false;
    int attempts = 0;

    while (!success)
    {
        Double_t precision = 1e-10; // 0 - default precision will be automaticalle determined
        Double_t tolerance = 0.001;
        ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
        min->SetMaxFunctionCalls(1000000000);
        min->SetMaxIterations(1000000);
        min->SetTolerance(tolerance);
        min->SetPrecision(precision);
        min->SetStrategy(2);
        
        ROOT::Math::Functor f(&ObjectiveFunction, NPAR);
        min->SetFunction(f);

        double current_start[NPAR];

        for (int i = 0; i < NPAR; ++i)
        {
            current_start[i] = gRandom->Uniform(par_min[i] * 0.9, par_max[i] * 0.9);
            min->SetLimitedVariable(i, names[i].c_str(), current_start[i], step_size[i], par_min[i], par_max[i]);
        }

        bool minim_ok = min->Minimize();
        double errval = min->MinValue();
        const double* best = min->X();

        // Startwert-Ausgabe
        cout << "Versuch #" << attempts << ":\n";
        cout << "  Startwerte:\n";
        for (int i = 0; i < NPAR; ++i)
            //cout << "    " << names[i] << " = " << current_start[i] << (i % 2 == 0 ? " deg" : " cm") << endl;
            cout << "    " << names[i] << " = " << best[i] << (i % 2 == 0 ? " deg" : " cm") << endl;

        // Ergebnis
        cout << "  Minimierer sagt: " << (minim_ok ? "erfolgreich" : "fehlgeschlagen") << endl;
        cout << "  Gefundener Fehlerwert: " << errval << " cm²\n";

        // Zielvergleich (Detailanalyse)
        double a11 = (-18.0 ) * TMath::DegToRad();
        double z11 = fi11_to_TP ;
        double a10 = (-18.0 ) * TMath::DegToRad();
        double z10 = fi10_to_TP ;
        double a12 = (-18. + best[0]) * TMath::DegToRad();
        double z12 = fi12_to_TP + best[1];

        TVector3 pos11(z11 * sin(a11), 0, Target_to_TP + z11 * cos(a11));
        TVector3 pos10(z10 * sin(a10), 0, Target_to_TP + z10 * cos(a10));
        TVector3 pos12(z12 * sin(a12), 0, Target_to_TP + z12 * cos(a12));

        TVector3 rot11(0, a11, 0), rot10(0, a10, 0), rot12(0, a12, 0);

        /*TVector3 fi11_lab_223(-1.*(-3.20910 ) , 0, 0);
        TVector3 fi10_lab_223(-1.*(-3.20786) + best[6], 0, 0);
        TVector3 fi12_lab_223(-1.*(-11.2588) + best[7], 0, 0);
        TVector3 fi11_lab_227(-1.*(-4.60759), 0, 0);
        TVector3 fi10_lab_227(-1.*(-5.20989) + best[6], 0, 0);
        TVector3 fi12_lab_227(-1.*(-16.74) + best[7], 0, 0);
        TVector3 fi11_lab_229(-1.*(-2.73996), 0, 0);
        TVector3 fi10_lab_229(-1.*(-2.53834) + best[6], 0, 0);
        TVector3 fi12_lab_229(-1.*(-9.38697) + best[7], 0, 0);
        TVector3 fi11_lab_202(-1.*(-2.80716) , 0, 0);
        TVector3 fi10_lab_202(-1.*(-2.63749) + best[6], 0, 0);
        TVector3 fi12_lab_202(-1.*(-9.61590) + best[7], 0, 0);
        TVector3 fi11_lab_201(-1.*(-3.13696) , 0, 0);
        TVector3 fi10_lab_201(-1.*(-3.09914) + best[6], 0, 0);
        TVector3 fi12_lab_201(-1.*(-10.3187) + best[7], 0, 0);*/
        
        TVector3 fi11_lab_223(-1.*(-3.20910) , 0, 0);
        TVector3 fi10_lab_223(-1.*(-3.20786) , 0, 0);
        TVector3 fi12_lab_223(-1.*(-11.2588) +best[3], 0, 0);
        TVector3 fi11_lab_227(-1.*(-4.60759) , 0, 0);
        TVector3 fi10_lab_227(-1.*(-5.20989) , 0, 0);
        TVector3 fi12_lab_227(-1.*(-16.74) +best[3], 0, 0);
        TVector3 fi11_lab_229(-1.*(-2.73996), 0, 0);
        TVector3 fi10_lab_229(-1.*(-2.53834) , 0, 0);
        TVector3 fi12_lab_229(-1.*(-9.38697) +best[3], 0, 0);
        TVector3 fi11_lab_202(-1.*(-2.80716) , 0, 0);
        TVector3 fi10_lab_202(-1.*(-2.63749) , 0, 0);
        TVector3 fi12_lab_202(-1.*(-9.61590) +best[3], 0, 0);
        TVector3 fi11_lab_201(-1.*(-3.13696) , 0, 0);
        TVector3 fi10_lab_201(-1.*(-3.09914) , 0, 0);
        TVector3 fi12_lab_201(-1.*(-10.3187) +best[3], 0, 0); 
        
        TVector3 fi10_lab_zero(0., 0., 0.);
        TVector3 fi11_lab_zero(0., 0., 0.);
        TVector3 fi12_lab_zero(0., 0., 0.);       

        TVector3 t11_229 = fi11_lab_229, t10_229 = fi10_lab_229, t12_229 = fi12_lab_229;
        TVector3 t11_227 = fi11_lab_227, t10_227 = fi10_lab_227, t12_227 = fi12_lab_227;
        TVector3 t11_223 = fi11_lab_223, t10_223 = fi10_lab_223, t12_223 = fi12_lab_223;
        TVector3 t11_202 = fi11_lab_202, t10_202 = fi10_lab_202, t12_202 = fi12_lab_202;
        TVector3 t11_201 = fi11_lab_201, t10_201 = fi10_lab_201, t12_201 = fi12_lab_201;
        TVector3 t11_zero = fi11_lab_zero, t10_zero = fi10_lab_zero, t12_zero = fi12_lab_zero;

        TransformPoint(t11_229, rot11, pos11);
        TransformPoint(t10_229, rot10, pos10);
        TransformPoint(t12_229, rot12, pos12);
        TransformPoint(t11_227, rot11, pos11);
        TransformPoint(t10_227, rot10, pos10);
        TransformPoint(t12_227, rot12, pos12);
        TransformPoint(t11_223, rot11, pos11);
        TransformPoint(t10_223, rot10, pos10);
        TransformPoint(t12_223, rot12, pos12);
        TransformPoint(t11_202, rot11, pos11);
        TransformPoint(t10_202, rot10, pos10);
        TransformPoint(t12_202, rot12, pos12);
        TransformPoint(t11_201, rot11, pos11);
        TransformPoint(t10_201, rot10, pos10);
        TransformPoint(t12_201, rot12, pos12);
        
        TransformPoint(t11_zero, rot11, pos11);
        TransformPoint(t10_zero, rot10, pos10);
        TransformPoint(t12_zero, rot12, pos12);

        R3BTrackingParticle ion229 = Generate_Primary(0.784798, best[2]);
        R3BTrackingParticle ion227 = Generate_Primary(0.592354, best[2]);
        R3BTrackingParticle ion223 = Generate_Primary(0.73866, best[2]);
        R3BTrackingParticle ion202 = Generate_Primary(0.380468, best[2]);
        R3BTrackingParticle ion201 = Generate_Primary(0.317634, best[2]);

        TVector3 hit11_229, hit10_229, hit12_229;
        TVector3 hit11_227, hit10_227, hit12_227;
        TVector3 hit11_223, hit10_223, hit12_223;
        TVector3 hit11_202, hit10_202, hit12_202;
        TVector3 hit11_201, hit10_201, hit12_201;

        bool valid = true;

        field->SetScale(-2641.0 / 3583.81);
        valid &= Propagate(ion229, t11_229.Z(), t11_229.X(), t11_zero.Z(), t11_zero.X(), hit11_229);
        valid &= Propagate(ion229, t10_229.Z(), t10_229.X(), t10_zero.Z(), t10_zero.X(), hit10_229);
        valid &= Propagate(ion229, t12_229.Z(), t12_229.X(), t12_zero.Z(), t12_zero.X(), hit12_229);

        field->SetScale(-2166.0 / 3583.81);
        valid &= Propagate(ion227, t11_227.Z(), t11_227.X(), t11_zero.Z(), t11_zero.X(), hit11_227);
        valid &= Propagate(ion227, t10_227.Z(), t10_227.X(), t10_zero.Z(), t10_zero.X(), hit10_227);
        valid &= Propagate(ion227, t12_227.Z(), t12_227.X(), t12_zero.Z(), t12_zero.X(), hit12_227);

        field->SetScale(-2527.0 / 3583.81);
        valid &= Propagate(ion223, t11_223.Z(), t10_223.X(), t11_zero.Z(), t11_zero.X(), hit11_223);
        valid &= Propagate(ion223, t10_223.Z(), t10_223.X(), t10_zero.Z(), t10_zero.X(), hit10_223);
        valid &= Propagate(ion223, t12_223.Z(), t12_223.X(), t12_zero.Z(), t12_zero.X(), hit12_223);        
                
        field->SetScale(-1691.0 / 3583.81);
        valid &= Propagate(ion202, t11_202.Z(), t11_202.X(), t11_zero.Z(), t11_zero.X(), hit11_202);
        valid &= Propagate(ion202, t10_202.Z(), t10_202.X(), t10_zero.Z(), t10_zero.X(), hit10_202);
        valid &= Propagate(ion202, t12_202.Z(), t12_202.X(), t12_zero.Z(), t12_zero.X(), hit12_202);
        
        field->SetScale(-1520.0 / 3583.81);
        valid &= Propagate(ion201, t11_201.Z(), t11_201.X(), t11_zero.Z(), t11_zero.X(), hit11_201);
        valid &= Propagate(ion201, t10_201.Z(), t10_201.X(), t10_zero.Z(), t10_zero.X(), hit10_201);
        valid &= Propagate(ion201, t12_201.Z(), t12_201.X(), t12_zero.Z(), t12_zero.X(), hit12_201);

        if (valid)
        {
            vector<pair<string, pair<TVector3, TVector3>>> checks = {
                //{"229 @ F11", {hit11_229, t11_229}},
                {"229 @ F10", {hit10_229, t10_229}},
                {"229 @ F12", {hit12_229, t12_229}},
                //{"227 @ F11", {hit11_227, t11_227}},
                {"227 @ F10", {hit10_227, t10_227}},
                {"227 @ F12", {hit12_227, t12_227}},
                //{"223 @ F11", {hit11_223, t11_223}},
                {"223 @ F10", {hit10_223, t10_223}},
                {"223 @ F12", {hit12_223, t12_223}},
                {"202 @ F10", {hit10_202, t10_202}},
                {"202 @ F12", {hit12_202, t12_202}},
                {"201 @ F10", {hit10_201, t10_201}},
                {"201 @ F12", {hit12_201, t12_201}},                
            };
            if(minim_ok==true){
            h_par1->Fill(best[0]);
            h_par2->Fill(best[1]);
            //h_par3->Fill(best[2]);
            }
            if(minim_ok==true){num++;}
            for (auto& [label, pair] : checks)
            {
                auto hit = pair.first;
                auto target = pair.second;
                double dx = hit.X() - target.X();
                double dz = hit.Z() - target.Z();

                cout << "  " << label << ":\n";
                cout << "    Hit:    X=" << hit.X() << ", Z=" << hit.Z() << "\n";
                cout << "    Target: X=" << target.X() << ", Z=" << target.Z() << "\n";
                cout << "    ΔX=" << dx << ", ΔZ=" << dz;
                if (fabs(dx) <= TOL && fabs(dz) <= TOL)
                    cout << " ✅ innerhalb Toleranz\n";
                else
                    cout << " ❌ außerhalb Toleranz\n";
            }
        }
        else
        {
            cout << "  ❌ Propagation fehlgeschlagen\n";
        }

        cout << "----------------------------------------------------------\n";

        // Erfolgskriterium: minim_ok und Fehlerwert realistisch
        success = errval < 1e4;

        delete min;
    }
    }
    
    TCanvas* c_par = new TCanvas();
    c_par->Divide(3,1);
    c_par->cd(1);
    h_par1->Draw("hist");
    c_par->cd(2);
    h_par2->Draw("hist"); 
    c_par->cd(3);
    h_par3->Draw("hist");    
    
    c_par->SaveAs("parameter_scan.pdf");
    
    TApplication* app = new TApplication("App", 0, nullptr);
    app->Run();

    return 0;
    
}
