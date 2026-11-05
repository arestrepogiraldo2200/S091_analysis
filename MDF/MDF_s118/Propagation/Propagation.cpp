#include"libs.hh"
#include"definitions.hh"
#include"R3BTrackingParticle.h"
#include"R3BGladFieldMap.h"
#include"R3BTPropagator.h"
#include"R3BMCTrack.h"
using namespace std;

R3BTrackingParticle Generate_Primary()
{
    const Double_t  Ekin_Total = BEAM_EKIN_GEV_U * BEAM_MASS_GEV_C2/AMU;
    Double_t P, Theta, RCos, Phi, vx, vy, vz, Rand_R, beta;
    TVector3 PVector;
    P     = sqrt(Ekin_Total*(Ekin_Total+2*BEAM_MASS_GEV_C2)); //initial momentum
    auto Ptot  = P;
    P       = (Ptot+gRandom->Uniform((-1)*Ptot*SCALE_MIN, Ptot*SCALE_MAX));//min and max momentum spread

    while(1)
    {
        RCos    = gRandom->Uniform(TMath::Cos(BEAM_DTHETA_DEG*TMath::Pi()/180.),1.); 
        Theta   = TMath::ACos(RCos); 
        Phi     = gRandom->Uniform(0., 2.*TMath::Pi());
        Rand_R  = gRandom->Uniform(0, BEAM_RADIUS_CM);
        vx	    = Rand_R * cos(Phi);
        vy	    = Rand_R * sin(Phi);
        vz = -20;//somehow R3BTPropagator fails at z=0
        PVector.SetMagThetaPhi(P, Theta, Phi);
        cout << "\nX: " << vx + (Target_to_RP-vz) * PVector.X()/PVector.Z();
        cout << "\nY: " << vy + (Target_to_RP-vz) * PVector.Y()/PVector.Z();
        //Do not use events which may potentially hit the wall of GLAD chamber
        if(fabs(vx + (Target_to_RP-vz) * PVector.X()/PVector.Z() )<13. && 
                fabs(vy + (Target_to_RP-vz) * PVector.Y()/PVector.Z() )<13.)
            break;
    }
    beta = 1. / TMath::Sqrt(1 + TMath::Power(BEAM_MASS_GEV_C2 / PVector.Mag(), 2));
    R3BTrackingParticle ion(BEAM_Z,
            vx,
            vy,
            vz,
            PVector.X(),
            PVector.Y(),
            PVector.Z(),
            beta,
            BEAM_MASS_GEV_C2);
    return ion;
}

void Propagation()
{
    TApplication* theApp = new TApplication("App", 0, 0);
    cout << "\n\t*************************************************" << endl;
    cout << "\t*                                               *" << endl;
    cout << "\t*  Propagating particles through GLAD field     *" << endl;
    cout << "\t*                                               *" << endl;
    cout << "\t*************************************************" << endl;

    cout << "\n-- Initializing GLAD field" << endl;
    //R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap","A");
    R3BGladFieldMap* magField = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-4to13_step10mm","R");
    magField->SetPosition(0.7871, 1.75-1.526, Target_to_GLAD_flange + 54.05 - 0.55580);//x,y,z in cm 
    magField->SetXAngle(-0.113); //deg
    magField->SetYAngle(-14.08); //deg
    magField->SetZAngle(-0.83); //deg
    magField->SetScale(GLAD_CURRENT/3583.81);
    magField->Init();
    magField->Print();

    cout << "\n-- Initializing field propagator";
    R3BTPropagator * fPropagator = new R3BTPropagator(magField, kFALSE);

    cout << "\n-- Plotting field map";
    double x,y,z,B;
    TH2D* h2 = new TH2D("h2","h2", NbinsZ, Zmin, Zmax, NbinsX, Xmin, Xmax);
    for (auto i = 0; i < NbinsZ; ++i)
    {
        z = Zmin + (i + 0.5) * (Zmax-Zmin)/NbinsZ;
        y=0;
        for (auto j = 0; j < NbinsX; ++j)
        {
            x = Xmin + (j + 0.5) * (Xmax-Xmin)/NbinsX;
            B = 0.1 * sqrt(pow(magField->GetBx(x, y, z),2) + pow(magField->GetBy(x, y, z),2)+ pow(magField->GetBz(x, y, z),2));//0.1 for kG->Tesla
            h2->Fill(z, x, B);
        }
    }

    //=========== Drawing setup configuration ======================
    TCanvas * c1 = new TCanvas("c1","c1",800,800);
    TPaletteAxis *palette;
    c1->cd(1);
    h2->Draw("colz Z");
    h2->GetZaxis()->SetRangeUser(1e-6, 3);
    c1->Update();
    c1->SetLogz();
    palette=(TPaletteAxis*)h2->FindObject("palette");
    palette->Draw();

    //Drawing entrance flange
    TLine* line_entrance_flange = new TLine(Target_to_GLAD_flange, -50., Target_to_GLAD_flange, 50.);
    line_entrance_flange->Draw();
    line_entrance_flange->SetLineWidth(1);
    c1->Update();

    //Drawing Rotation point
    TLine* line_rp = new TLine(Target_to_RP, 0., Target_to_RP, 20.);
    line_rp->SetLineWidth(2);
    line_rp->SetLineColor(kBlue);
    line_rp->Draw();
    c1->Update();

    //Drawing Turning point 
    TLine* line_tp = new TLine(Target_to_TP, 0., Target_to_TP, 20.);
    line_tp->SetLineWidth(2);
    line_tp->SetLineColor(kBlue);
    line_tp->Draw();
    c1->Update();

    //Drawing F11 (first, X coord)
    TVector3 p0(fi11_to_TP * TMath::Sin(Angle), 0, Target_to_TP + fi11_to_TP * TMath::Cos(Angle));
    TVector3 p1(-25., 0., 0.);
    TVector3 p2(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_f11 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_f11->SetLineWidth(2);
    line_f11->Draw();
    c1->Update();

    //Drawing F10 (second, X coord)
    p0.SetXYZ(fi10_to_TP * TMath::Sin(Angle), 0, Target_to_TP + fi10_to_TP * TMath::Cos(Angle));
    p1.SetXYZ(-25., 0., 0.);
    p2.SetXYZ(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_f10 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_f10->SetLineWidth(2);
    line_f10->Draw();
    c1->Update();

    //Drawing F12 (3rd, X coord)
    p0.SetXYZ(fi12_to_TP * TMath::Sin(Angle), 0, Target_to_TP + fi12_to_TP * TMath::Cos(Angle));
    p1.SetXYZ(-25., 0., 0.);
    p2.SetXYZ(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_f12 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_f12->SetLineWidth(2);
    line_f12->Draw();
    c1->Update();

    //Drawing F13 (4th, Y coord)
    p0.SetXYZ(fi13_to_TP * TMath::Sin(Angle), 0, Target_to_TP + fi13_to_TP * TMath::Cos(Angle));
    p1.SetXYZ(-25., 0., 0.);
    p2.SetXYZ(25., 0., 0.);
    p1.RotateY(Angle);
    p2.RotateY(Angle);
    p1 += p0;
    p2 += p0;
    TLine* line_f13 = new TLine(p1.Z(),p1.X(),p2.Z(),p2.X());
    line_f13->SetLineWidth(2);
    line_f13->Draw();
    c1->Update();    

    cout << "\n-- Creating particle tracks" << endl;
    TVector3 momentum, position;
    TVector3 track_pos;
    TVector3 v1, v2, v3;//points on the the plane for propagation

    //Trajectory object
    double Xtraj[NbinsZ];
    double Ztraj[NbinsZ];
    TGraph * trajectory;
    R3BTrackingParticle * particle;
    R3BTrackingParticle primary;
    
    //particle = new R3BTrackingParticle(BEAM_Z,
    //            position.X(),
    //            position.Y(),
    //            position.Z(),
    //            momentum.X(),
    //            momentum.Y(),
    //            momentum.Z(),
    //            primary.GetBeta(),
    //            primary.GetMass());
    particle = new R3BTrackingParticle(0,0,0,0,0,0,0,0,0);

    for(int ev=0; ev<NEVENTS; ev++)
    {
        for(auto k=0; k<NbinsZ; k++)
        {
            Xtraj[k]=0;
            Ztraj[k]=0;
        }

        //cout << "\r-- Working on entry " << ev << flush;
        cout << "\n-- Working on entry " << ev;
        primary = Generate_Primary();
        position = primary.GetStartPosition();
        momentum = primary.GetStartMomentum();
        Double_t position_array[3];
        position_array[0] = position.X();
        position_array[1] = position.Y();
        position_array[2] = position.Z();

        for (auto i = 0; i < NbinsZ; i++)
        {
            z = (i + 1) * Zmax / NbinsZ;
            if(z<position.Z()) continue;

            particle->SetCharge(BEAM_Z);
            particle->SetPosition(position_array);
            particle->SetMomentum(momentum);
            particle->SetStartBeta(primary.GetBeta());
            particle->SetBeta(primary.GetBeta());
            particle->SetMass(primary.GetMass());

            v1.SetXYZ(0,0,z);
            v2.SetXYZ(1,0,z);
            v3.SetXYZ(0,1,z);

            double tof=0;
            //if(!fPropagator->PropagateToPlaneRK(particle,v1,v2,v3))
            if(!fPropagator->PropagateToPlaneRK_eloss(particle,v1,v2,v3,tof,10,false)){
                cout << "\nFailed propagation...";
                return;
            }
            track_pos = particle->GetPosition();
            if(isnan(track_pos.X()) || isnan(track_pos.Z()))
            {
                cout << "\n\nERROR! Propagated NaN values! @ Z=" << z << "cm\n";
                cout << "\nStart position X, Y, Z: " << position.X() << "\t" << position.Y() << "\t" << position.Z()<< " cm";
                cout << "\nStart angles TX, TY: " << momentum.X()/momentum.Z() << "\t" << momentum.Y()/momentum.Z() << " \n";

                return;
            }
            Xtraj[i]=track_pos.X();
            Ztraj[i]=track_pos.Z();

        }
        //delete particle;

        trajectory  = new TGraph(NbinsZ,Ztraj,Xtraj);
        trajectory->SetMarkerStyle(kFullDotSmall);
        trajectory->Draw("same");
        c1->Update();
    }

    //Beam axis (Z-axis)
    /*TLine* line_beamline = new TLine(Zmin, 0., Zmax, 0.);
    line_beamline->SetLineWidth(2);
    line_beamline->SetLineColor(kBlue);
    line_beamline->Draw();*/

    //Drawing 14 deg line
    /*TLine* line_18deg = new TLine(Target_to_TP, 0., Zmax, (Zmax - Target_to_TP) * TMath::Tan(-18.*TMath::DegToRad()));
    line_18deg->SetLineWidth(2);
    line_18deg->SetLineColor(kBlue);
    line_18deg->Draw();
    c1->Update();*/

    cout << "\n-- FINISH --" << endl;
    c1->Update();
    theApp->Run();
    return;
}
int main(Int_t argc, Char_t* argv[])
{
    //Char_t *in_filename=0;
    //Bool_t is_PCA = kFALSE;
    //Bool_t NeedHelp = kTRUE;

    //if (argc > 1){
    //    for (Int_t i = 0; i < argc; i++){
    //        if (strncmp(argv[i],"--file=",7) == 0){
    //            in_filename = argv[i]+7;
    //            NeedHelp = kFALSE;
    //        }

    //        else if (strncmp(argv[i],"--pca",5) == 0){
    //            is_PCA = kTRUE;
    //        }
    //    }
    //}
    //if (NeedHelp){
    //    cout << "\nUse the following flags: ";
    //    cout << "\n--file=/path/to/your/file/filename.root   : (mandatory) simulation input file";
    //    return 0;
    //}

    gRandom = new TRandom3();
    gRandom->SetSeed(0);
    //gROOT->Macro("~/rootlogon.C");
    gStyle->SetPalette(kRainBow);

    Propagation();
    return 0;
}
