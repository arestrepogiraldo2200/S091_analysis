#include <algorithm>
#include <array>
#include <map>
#include <vector>
#include <fstream>

R3BMDFWrapper* MDF_PoQ;
 
TVector3 fi32_lab, fi31_lab, fi33_lab, mwpc1_lab;
TVector3 fi32_position, fi31_position, fi33_position, mwpc1_position;
TVector3 fi32_angles  , fi31_angles  , fi33_angles  , mwpc1_angles  ;
TVector3 fi32_pos_offset, fi31_pos_offset, fi33_pos_offset, mwpc1_pos_offset;
TVector3 fi32_ang_offset, fi31_ang_offset, fi33_ang_offset, mwpc1_ang_offset;

// Distances in cm ======================================================================================
const Double_t Target_to_GLAD_flange  = 127.9;                                         // <-------------------
const Double_t targ_to_TP   = 293.1;                                                   // <-------------------
const Double_t TP_to_GLAD_exit_flange = 385.0 - 125.0;                                 // <-------------------
const Double_t Target_to_RP   = Target_to_GLAD_flange +  54.05 - 0.55580;              // <-------------------

const Double_t targ_to_mwpc1 = 95.7;                              // target to MWPC1                 // <-------------------
const Double_t fi_32_to_TP   = 385.0;                            // turning point to fib32           // <-------------------
const Double_t fi_31_to_TP   = TP_to_GLAD_exit_flange + 242.5;   // turning point to fib31           // <-------------------
const Double_t fi_33_to_TP   = TP_to_GLAD_exit_flange + 225.0;   // turning point to fib33           // <-------------------
const Double_t tofd_to_TP    = TP_to_GLAD_exit_flange + 661.8;   // turning point to fib33           // <-------------------

// PoQ reference values ======================================================================================
double Reference_PoQ_10 = 1.527872;
double Reference_PoQ_11 = 1.7092484;
double Reference_PoQ_12 = 1.9902324;
double Reference_PoQ_13 = 2.0717421;
double Reference_PoQ_14 = 2.2532898;

// Transform point ======================================================================================
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

// Function to calculate PoQ ======================================================================================
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

// Minimization function ======================================================================================
// mwpc1x -> Must be aligned to 0
// fib32x -> Must be aligned to 0
// fib33x -> Must be aligned to 0 of 32
// fib31x -> Must be aligned to 0 of 32. Perpendicular projection.
double MinimizeDelta(const double* hit_x_positions, const double* par, double Reference_PoQ, TString MDF_PoQ_file, double Actual_Current, double Reference_Current) {

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

    double mwpc1x = hit_x_positions[0];
    double fib32x = hit_x_positions[1];
    double fib33x = hit_x_positions[2];
    double fib31x = hit_x_positions[3];

    bool hit_33 = true;
    bool hit_31 = true;

    // Flags to tell if there was no hit in a certain fiber
    if (fib33x == 0.00000){
        hit_33 = false;
    }

    if (fib31x == 0.00000){
        hit_31 = false;
    }

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

    Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;  
    
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
    fi32_pos_offset.SetXYZ( par3, 0,  par4);
    fi33_pos_offset.SetXYZ( par7, 0,  par8);
    fi31_pos_offset.SetXYZ(par11, 0, par12);

    mwpc1_ang_offset.SetXYZ(0, 0, 0);
    fi32_ang_offset.SetXYZ(0,  par5,  par6 );
    fi33_ang_offset.SetXYZ(0,  par9,  par10);
    fi31_ang_offset.SetXYZ(0, par13,  par14);

    // Set the positions of hits 
    mwpc1_lab.SetXYZ(mwpc1x, 0, 0);
    fi32_lab.SetXYZ(fib32x, 0, 0);
    fi33_lab.SetXYZ(fib33x, 0, 0);
    fi31_lab.SetXYZ(fib31x, 0, 0);

    if ( (hit_31 == false && hit_33 == true) || (hit_31 == true && hit_33 == true) ){

        // Up to this point the x and z positions of the fibers and mwpc are being placed in the lab frame 
        TransformPoint(mwpc1_lab, mwpc1_angles + mwpc1_ang_offset, mwpc1_position + mwpc1_pos_offset);
        TransformPoint( fi32_lab,   fi32_angles + fi32_ang_offset,   fi32_position + fi32_pos_offset);
        TransformPoint( fi33_lab,   fi33_angles + fi33_ang_offset,   fi33_position + fi33_pos_offset);

        double delta = 0.0;
        double PoQ_calculated = CalculatePoQ(mwpc1_lab, fi33_lab, fi32_lab, Actual_Current, Reference_Current, MDF_PoQ);
        // cout << "Calculated: " << PoQ_calculated << "   Reference: " << Reference_PoQ << endl;
        delta = fabs(Reference_PoQ - PoQ_calculated);
        return delta;
    }

    if (hit_31 == true && hit_33 == false){

        // Up to this point the x and z positions of the fibers and mwpc are being placed in the lab frame 
        TransformPoint(mwpc1_lab, mwpc1_angles + mwpc1_ang_offset, mwpc1_position + mwpc1_pos_offset);
        TransformPoint( fi32_lab,   fi32_angles + fi32_ang_offset,   fi32_position + fi32_pos_offset);
        TransformPoint( fi31_lab,   fi31_angles + fi31_ang_offset,   fi31_position + fi31_pos_offset);

        double delta = 0.0;
        double PoQ_calculated = CalculatePoQ(mwpc1_lab, fi31_lab, fi32_lab, Actual_Current, Reference_Current, MDF_PoQ);
        // cout << "Calculated: " << PoQ_calculated << "   Reference: " << Reference_PoQ << endl;
        delta = fabs(Reference_PoQ - PoQ_calculated);
        return delta;
    }
}

// Several runs ======================================================================================
double MinimizeDelta_several_runs(const std::vector<std::array<double, 4>>& hit_x_positions, const double* par, const std::vector<double>& Reference_PoQ, TString MDF_PoQ_file, const std::vector<double>& Actual_Current, const std::vector<double>& Reference_Current) {

    const std::size_t num_entries = hit_x_positions.size();
    std::vector<double> errors(num_entries);

    for (std::size_t i = 0; i < num_entries; ++i) {

        errors[i] = MinimizeDelta(hit_x_positions[i].data(), par, Reference_PoQ[i], MDF_PoQ_file, Actual_Current[i], Reference_Current[i]);
        // cout << "Error " << i << " is: " << errors[i] << endl;
    }

    return *std::max_element(errors.begin(), errors.end());
}

// Loads the detector positions ======================================================================================
std::vector<std::array<double, 4>> LoadDetectorPositions(const TString& fileName)
{
    std::ifstream input(fileName.Data());
    if (!input.is_open()) {
        throw std::runtime_error(Form("Failed to open %s", fileName.Data()));
    }

    std::vector<std::array<double, 4>> rows;
    std::array<double, 4> row{};

    while (input >> row[0] >> row[1] >> row[2] >> row[3]) {
        rows.push_back(row);
    }

    return rows;
}

// Return the average of the columns ======================================================================================
std::array<double, 4> AverageColumns(const std::vector<std::array<double, 4>>& rows)
{
    std::array<double, 4> averages{0.0, 0.0, 0.0, 0.0};

    if (rows.empty()) {
        return averages;
    }

    for (const auto& row : rows) {
        for (std::size_t i = 0; i < 4; ++i) {
            averages[i] += row[i];
        }
    }

    for (double& value : averages) {
        value /= static_cast<double>(rows.size());
    }

    return averages;
}

// Return the mode of the columns ======================================================================================
std::array<double, 4> ModeColumns(const std::vector<std::array<double, 4>>& rows)
{
    std::array<double, 4> modes{0.0, 0.0, 0.0, 0.0};

    if (rows.empty()) {
        return modes;
    }

    std::array<std::map<double, std::size_t>, 4> counts;

    for (const auto& row : rows) {
        for (std::size_t i = 0; i < 4; ++i) {
            ++counts[i][row[i]];
        }
    }

    for (std::size_t i = 0; i < 4; ++i) {
        std::size_t best_count = 0;
        double best_value = 0.0;

        for (const auto& entry : counts[i]) {
            if (entry.second > best_count) {
                best_count = entry.second;
                best_value = entry.first;
            }
        }

        modes[i] = best_value;
    }

    return modes;
}
