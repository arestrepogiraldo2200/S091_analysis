#include "libs.hh"
#include "definitions.hh"
#include "R3BGladFieldMap.h"
#include "R3BTPropagator.h"
#include "R3BTrackingParticle.h"
#include <iostream>
#include <limits>

using namespace std;

R3BGladFieldMap* field = nullptr;
R3BTPropagator* propagator = nullptr;
const double TOL = 0.3;

double ObjectiveFunction(const double* p);

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

R3BTrackingParticle Generate_Primary(double ekin_u)
{
    const Double_t Ekin_Total = ekin_u * BEAM_MASS_GEV_C2 / AMU;
    Double_t P = sqrt(Ekin_Total * (Ekin_Total + 2 * BEAM_MASS_GEV_C2));
    TVector3 mom(0, 0, P);
    TVector3 pos(0, 0, 0);
    Double_t beta = 1. / TMath::Sqrt(1 + TMath::Power(BEAM_MASS_GEV_C2 / mom.Mag(), 2));
    return R3BTrackingParticle(BEAM_Z, pos.X(), pos.Y(), pos.Z(), mom.X(), mom.Y(), mom.Z(), beta, BEAM_MASS_GEV_C2);
}

bool Propagate(R3BTrackingParticle& particle, double z, TVector3& result)
{
    TVector3 v1(0, 0, z), v2(1, 0, z), v3(0, 1, z);
    double tof = 0;
    if (!propagator->PropagateToPlaneRK_eloss(&particle, v1, v2, v3, tof, 0.1, false)) return false;
    result = particle.GetPosition();
    return true;
}

int main()
{
    gRandom = new TRandom3(0);
    field = new R3BGladFieldMap("R3BGladMap_Bxyz_X-3to3_Y-1to1_Z-4to13_step10mm", "R");
    field->SetPosition(0.7871, 1.75 - 1.526, Target_to_GLAD_flange + 54.05 - 0.55580);
    field->SetXAngle(-0.113);
    field->SetYAngle(-14.08);
    field->SetZAngle(-0.83);
    field->Init();

    propagator = new R3BTPropagator(field, kFALSE);

    double best_error = 1e10;
    double best_angle = 0, best_zoffset = 0, best_xoffset = 0;

    for (double angle12 = -0.25; angle12 <= 0.25; angle12 += 0.01)
    {
        for (double z12 = -10; z12 <= 10; z12 += 0.5)
        {
            for (double x12 = 2.0; x12 <= 6.0; x12 += 0.1)
            {
                double p[8] = {0, 0, 0, 0, angle12, z12, 0, x12};
                double error = ObjectiveFunction(p);

                if (error < best_error)
                {
                    best_error = error;
                    best_angle = angle12;
                    best_zoffset = z12;
                    best_xoffset = x12;
                }
            }
        }
    }

    cout << "\nBeste Parameterkombination (nur Fi12 gescannt):" << endl;
    cout << "  Winkel [deg]:  " << best_angle << endl;
    cout << "  Z-Offset [cm]: " << best_zoffset << endl;
    cout << "  X-Offset [cm]: " << best_xoffset << endl;
    cout << "  Fehlerwert:    " << best_error << " cm^2" << endl;

    // Trefferpositionen für beste Konfiguration anzeigen
    double a12 = (-18. + best_angle) * TMath::DegToRad();
    double z12 = fi12_to_TP + best_zoffset;
    TVector3 pos12(z12 * sin(a12), 0, Target_to_TP + z12 * cos(a12));
    TVector3 rot12(0, a12, 0);

    TVector3 fi12_lab_227(-1.*(-16.882), 0, 0);
    TVector3 fi12_lab_229(-1.*(-9.38697), 0, 0);
    fi12_lab_227.SetX(fi12_lab_227.X() + best_xoffset);
    fi12_lab_229.SetX(fi12_lab_229.X() + best_xoffset);

    TransformPoint(fi12_lab_227, rot12, pos12);
    TransformPoint(fi12_lab_229, rot12, pos12);

    R3BTrackingParticle ion229 = Generate_Primary(0.784798);
    R3BTrackingParticle ion227 = Generate_Primary(0.594);

    TVector3 hit12_229, hit12_227;
    field->SetScale(-2641.0 / 3583.81);
    bool valid_229 = Propagate(ion229, fi12_lab_229.Z(), hit12_229);
    field->SetScale(-2166.0 / 3583.81);
    bool valid_227 = Propagate(ion227, fi12_lab_227.Z(), hit12_227);

    cout << "\nTrefferpositionen Fi12:" << endl;
    if (valid_229)
    {
        cout << "  229Sn Hit:    X=" << hit12_229.X() << " cm, Z=" << hit12_229.Z() << " cm" << endl;
        cout << "  229Sn Target: X=" << fi12_lab_229.X() << " cm, Z=" << fi12_lab_229.Z() << " cm" << endl;
    }
    if (valid_227)
    {
        cout << "  227Sn Hit:    X=" << hit12_227.X() << " cm, Z=" << hit12_227.Z() << " cm" << endl;
        cout << "  227Sn Target: X=" << fi12_lab_227.X() << " cm, Z=" << fi12_lab_227.Z() << " cm" << endl;
    }

    return 0;
}

double ObjectiveFunction(const double* p)
{
    double a12 = (-18. + p[4]) * TMath::DegToRad();
    double z12 = fi12_to_TP + p[5];

    TVector3 pos12(z12 * sin(a12), 0, Target_to_TP + z12 * cos(a12));
    TVector3 rot12(0, a12, 0);

    TVector3 fi12_lab_227(-1.*(-16.882), 0, 0);
    TVector3 fi12_lab_229(-1.*(-9.38697), 0, 0);
    fi12_lab_227.SetX(fi12_lab_227.X() + p[7]);
    fi12_lab_229.SetX(fi12_lab_229.X() + p[7]);

    TransformPoint(fi12_lab_227, rot12, pos12);
    TransformPoint(fi12_lab_229, rot12, pos12);

    R3BTrackingParticle ion229 = Generate_Primary(0.784798);
    R3BTrackingParticle ion227 = Generate_Primary(0.594);

    TVector3 hit12_229, hit12_227;
    field->SetScale(-2641.0 / 3583.81);
    if (!Propagate(ion229, fi12_lab_229.Z(), hit12_229)) return 1e6;
    field->SetScale(-2166.0 / 3583.81);
    if (!Propagate(ion227, fi12_lab_227.Z(), hit12_227)) return 1e6;

    vector<pair<TVector3, TVector3>> checks = {
        {hit12_229, fi12_lab_229},
        {hit12_227, fi12_lab_227}
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
