#include <array>
#include <vector>

#include "function.h"

void test_minimizer(){

    // Single hit trajectory ======================================================================================
    double Actual_Current = 1490.0;
    double Reference_Current = 1490.0;
    TString MDF_PoQ_file="./MDF_params.txt";
    double Reference_PoQ = 1.527872;
    // double par[14] = {0.0, 0.0,    0.0,  0.0,  0.0,  0.0,     0.0,  0.0,  0.0,  0.0,    0.0,  0.0,  0.0,   0.0};
    // double par[14] = {0.05, 0.05,    0.05,  0.05,  0.01,  0.01,     0.05,  0.05,  0.01,  0.01,    0.05,  0.05,  0.01,   0.01}
    // double par[14] = {-3.38931, 4.81549, -0.907715, -4.98417, -0.488445, -0.498736, -3.21726, 0.369934, -0.49921, 0.370774, -3.05182, 1.6021, -0.213565, -0.48408};
    double par[14] = {-2.56985, 0.815912, 0.991082, 2.46896, 0.406464, 0.259915, 0.78244, 3.42896, 0.433834, 0.124068, 1.15385, 2.11223, 0.337565, 0.282816};

    double hit_x_positions[4] = {-0., -26, -30, -1111};

    double Minim = MinimizeDelta(hit_x_positions, par, Reference_PoQ, MDF_PoQ_file, Actual_Current, Reference_Current);
    cout << "Test of output: " << Minim << endl;

    // Several hits isotopes ======================================================================================
    std::vector<std::array<double, 4>> hit_x_positions_several = {
        std::array<double, 4>{-0.163635, -23.337, 0.0, -30.6629},
        std::array<double, 4>{-0.000375727, -12.161, 0.0, -15.6697},
        std::array<double, 4>{0.0212206, 5.84064, 7.62374, 0.0},
        std::array<double, 4>{0.352435, 10.2366, 13.1563, 0.0},
        std::array<double, 4>{0.186022, 17.6546, 22.6421, 0.0}
    };

    std::vector<double> Actual_Current_several = {1490.0, 1490.0, 1490.0, 1490.0, 1490.0};
    std::vector<double> Reference_Current_several = {1490.0, 1490.0, 1490.0, 1490.0, 1490.0};
    std::vector<double> Reference_PoQ_several = {1.527872, 1.7092484, 1.9902324, 2.0717421, 2.2532898};

    double Minim_several = MinimizeDelta_several_runs(hit_x_positions_several, par, Reference_PoQ_several, MDF_PoQ_file, Actual_Current_several, Reference_Current_several);

    cout << "Test of output several: " << Minim_several << endl;

    // Average hits ======================================================================================
    std::vector<std::array<double, 4>> C10 = LoadDetectorPositions("./10C_detector_positions.txt");
    std::vector<std::array<double, 4>> C11 = LoadDetectorPositions("./11C_detector_positions.txt");
    std::vector<std::array<double, 4>> C12 = LoadDetectorPositions("./12C_detector_positions.txt");
    std::vector<std::array<double, 4>> C13 = LoadDetectorPositions("./13C_detector_positions.txt");
    std::vector<std::array<double, 4>> C14 = LoadDetectorPositions("./14C_detector_positions.txt");

    std::array<double, 4> C10_av = AverageColumns(C10);
    std::array<double, 4> C11_av = AverageColumns(C11);
    std::array<double, 4> C12_av = AverageColumns(C12);
    std::array<double, 4> C13_av = AverageColumns(C13);
    std::array<double, 4> C14_av = AverageColumns(C14);

    std::array<double, 4> C10_mode = ModeColumns(C10);
    std::array<double, 4> C11_mode = ModeColumns(C11);
    std::array<double, 4> C12_mode = ModeColumns(C12);
    std::array<double, 4> C13_mode = ModeColumns(C13);
    std::array<double, 4> C14_mode = ModeColumns(C14);

    // cout << "C10_av: {";
    // for (double value : C10_av) {
    //     cout << value << ", ";
    // }
    // cout << "}" << endl;

    // cout << "C11_av: {";
    // for (double value : C11_av) {
    //     cout << value << ", ";
    // }
    // cout << "}" << endl;

    // cout << "C12_av: {";
    // for (double value : C12_av) {
    //     cout << value << ", ";
    // }
    // cout << "}" << endl;

    // cout << "C13_av: {";
    // for (double value : C13_av) {
    //     cout << value << ", ";
    // }
    // cout << "}" << endl;

    // cout << "C14_av: {";
    // for (double value : C14_av) {
    //     cout << value << ", ";
    // }
    // cout << "}" << endl;

    // Several hits averages ======================================================================================
    std::vector<std::array<double, 4>> hit_x_positions_several_avr = {C10_av,
                                                                      C11_av,
                                                                      C12_av,
                                                                      C13_av,
                                                                      C14_av};
    
    //double Minim_several_avr = MinimizeDelta_several_runs(hit_x_positions_several_avr, par, Reference_PoQ_several, MDF_PoQ_file, Actual_Current_several, Reference_Current_several);
    // cout << "Test of output several: " << Minim_several_avr << endl;

    // Several hits modes ======================================================================================
    std::vector<std::array<double, 4>> hit_x_positions_several_mode = {C10_mode,
                                                                       C11_mode,
                                                                       C12_mode,
                                                                       C13_mode,
                                                                       C14_mode};
    
    //double Minim_several_mode = MinimizeDelta_several_runs(hit_x_positions_several_mode, par, Reference_PoQ_several, MDF_PoQ_file, Actual_Current_several, Reference_Current_several);
    // cout << "Test of output several: " << Minim_several_mode << endl;

}