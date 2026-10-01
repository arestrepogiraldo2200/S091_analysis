#include "Math/Functor.h"
#include "function.h"

char histname[200];
TH1D* h_pars[14];

TVector3 fi32_lab_opt, fi31_lab_opt, fi33_lab_opt, mwpc1_lab_opt;
TVector3 fi32_position_opt, fi31_position_opt, fi33_position_opt, mwpc1_position_opt;
TVector3 fi32_angles_opt  , fi31_angles_opt  , fi33_angles_opt  , mwpc1_angles_opt  ;
TVector3 fi32_pos_offset_opt, fi31_pos_offset_opt, fi33_pos_offset_opt, mwpc1_pos_offset_opt;
TVector3 fi32_ang_offset_opt, fi31_ang_offset_opt, fi33_ang_offset_opt, mwpc1_ang_offset_opt;

void minimizePars() {
    
    MDF_PoQ = new R3BMDFWrapper("./MDF_params.txt");

    const UInt_t num_pars = 14;
    const double target_delta = 0.5;

    // Limits for the parameters       MWPC                 FIB32                         FIB33                       FIB31
    double par_min[num_pars]   = {  -5,   -5,      -5,    -5,  -0.5,  -0.5,       -5,    -5,  -0.5,  -0.5,      -5,    -5,  -0.5,   -0.5};
    double par_max[num_pars]   = {   5,    5,       5,     5,   0.5,   0.5,        5,     5,   0.5,   0.5,       5,     5,   0.5,    0.5};    
    double step_size[num_pars] = {0.05, 0.05,    0.05,  0.05,  0.01,  0.01,     0.05,  0.05,  0.01,  0.01,    0.05,  0.05,  0.01,   0.01};

    // Histogram of parameters is created
    for(int i = 0; i < num_pars; i++){
        sprintf(histname, "hist par %i",i);
        h_pars[i] = new TH1D(histname, histname, 100, par_min[i], par_max[i]);   
    }

    int num = 0;
    while(num<20){

        // Initializes randomness
        std::srand(std::time(nullptr));

        bool success = false;
        int attempts = 0;

        num++;
        while (!success) {
            
            attempts++;
            std::cout << "\n--- Attempt " << attempts << " ---" << std::endl;

            // Minimizer setup
            Double_t precision = 1e-10; 
            Double_t tolerance = 0.001;
            ROOT::Math::Minimizer* minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");
            minimizer->SetMaxFunctionCalls(10000000);
            minimizer->SetMaxIterations(100000);
            minimizer->SetTolerance(tolerance);
            minimizer->SetPrecision(precision);
            minimizer->SetStrategy(2);

            // Parameters to minimize
            TString MDF_PoQ_file="./MDF_params.txt";
            std::vector<std::array<double, 4>> hit_x_position_several = {
                std::array<double, 4>{-0.163635, -23.337, 0.0, -30.6629},
                std::array<double, 4>{-0.000375727, -12.161, 0.0, -15.6697},
                std::array<double, 4>{0.0212206, 5.84064, 7.62374, 0.0},
                std::array<double, 4>{0.352435, 10.2366, 13.1563, 0.0},
                std::array<double, 4>{0.186022, 17.6546, 22.6421, 0.0}};
            std::vector<double> Actual_Current_several = {1490.0, 1490.0, 1490.0, 1490.0, 1490.0};
            std::vector<double> Reference_Current_several = {1490.0, 1490.0, 1490.0, 1490.0, 1490.0};
            std::vector<double> Reference_PoQ_several = {1.527872, 1.7092484, 1.9902324, 2.0717421, 2.2532898};

            auto objective = [&](const double* par) {
                return MinimizeDelta_several_runs(
                    hit_x_position_several,
                    par,
                    Reference_PoQ_several,
                    MDF_PoQ_file,
                    Actual_Current_several,
                    Reference_Current_several
                );
            };
            ROOT::Math::Functor f(objective, num_pars);
            minimizer->SetFunction(f);

            // Set initial values of the minimizer as random numbers in the interval of the object
            for (unsigned int i = 0; i < num_pars; ++i) {

                double start = gRandom->Uniform( par_min[i]*0.9,  par_max[i]*0.9);
                minimizer->SetLimitedVariable(i, ("par" + std::to_string(i+1)).c_str(), start, step_size[i], par_min[i], par_max[i]);
            }

            // Minimize
            minimizer->Minimize();
            double min_val = minimizer->MinValue();
            std::cout << "Found minimum: " << min_val << std::endl;

            if (min_val < target_delta) {

                const double* result = minimizer->X();

                for (int i = 0; i < num_pars; i++){
                    if ( fabs((result[i] - par_min[i]) / par_min[i]) < 0.1 || fabs((par_max[i] - result[i]) / par_max[i]) < 0.1){
                        break;
                    }
                }
                
                success = true;

                double par1 = result[0];      // mwpc x
                double par2 = result[1];      // mwpc z

                double par3 = result[2];      // fib32 x
                double par4 = result[3];      // fib32 z
                double par5 = result[4];      // fib32 y'
                double par6 = result[5];      // fib32 z'

                double par7 = result[6];      // fib33 x 
                double par8 = result[7];      // fib33 z 
                double par9 = result[8];      // fib33 y'
                double par10 = result[9];     // fib33 z'

                double par11 = result[10];    // fib31 x 
                double par12 = result[11];    // fib31 z 
                double par13 = result[12];    // fib31 y'
                double par14 = result[13];    // fib31 z'


                for (int hit=0; hit<=4; hit++){

                    double PoQ=0;

                    double mwpc1x = hit_x_position_several[hit][0];
                    double fib32x = hit_x_position_several[hit][1];
                    double fib33x = hit_x_position_several[hit][2];
                    double fib31x = hit_x_position_several[hit][3];

                    bool hit_33 = true;
                    bool hit_31 = true;

                    // Flags to tell if there was no hit in a certain fiber
                    if (fib33x == 0.00000){
                        hit_33 = false;
                    }

                    if (fib31x == 0.00000){
                        hit_31 = false;
                    }

                    for(int i = 0; i < 10; i++){
                        h_pars[i]->Fill(result[i]);
                    }

                    // Distance of fibers to TP
                    const Double_t fi_32_to_TP   = 385.0;                              // turning point to fib32           // <-------------------
                    const Double_t fi_33_to_TP   = TP_to_GLAD_exit_flange + 225.0;     // turning point to fib33           // <-------------------
                    const Double_t fi_31_to_TP   = TP_to_GLAD_exit_flange + 242.5;     // turning point to fib31           // <-------------------

                    // Current of the run
                    double Actual_Current = 1490.0;
                    double Reference_Current = 1490.0;
                    Double_t AngleAl = (-18.0 ) * TMath::Pi() / 180.;

                    // Set the positions of central values of the planes of the detectors with the optimized parameters
                    mwpc1_position_opt.SetXYZ(0,0,targ_to_mwpc1);
                    fi32_position_opt.SetXYZ(fi_32_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_32_to_TP * TMath::Cos(AngleAl));
                    fi33_position_opt.SetXYZ(fi_33_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_33_to_TP * TMath::Cos(AngleAl));
                    fi31_position_opt.SetXYZ(fi_31_to_TP * TMath::Sin(AngleAl), 0, targ_to_TP + fi_31_to_TP * TMath::Cos(AngleAl)); 
                    
                    mwpc1_angles_opt.SetXYZ(0, 0, 0);
                    fi32_angles_opt.SetXYZ(0, AngleAl, 0);
                    fi33_angles_opt.SetXYZ(0, AngleAl, 0);
                    fi31_angles_opt.SetXYZ(0, AngleAl, 0);

                    // Set the variations of positions and angles as the parameters
                    mwpc1_pos_offset_opt.SetXYZ(par1, 0,  par2);
                    fi32_pos_offset_opt.SetXYZ( par3,  0,  par4);
                    fi33_pos_offset_opt.SetXYZ( par7, 0,  par8);
                    fi31_pos_offset_opt.SetXYZ(par11, 0, par12);

                    mwpc1_ang_offset_opt.SetXYZ(0, 0, 0);
                    fi32_ang_offset_opt.SetXYZ(0,  par5,  par6 );
                    fi33_ang_offset_opt.SetXYZ(0,  par9,  par10);
                    fi31_ang_offset_opt.SetXYZ(0, par13,  par14);

                    // Set the positions of hits 
                    mwpc1_lab_opt.SetXYZ(mwpc1x, 0, 0);
                    fi32_lab_opt.SetXYZ(fib32x, 0, 0);
                    fi33_lab_opt.SetXYZ(fib33x, 0, 0);
                    fi31_lab_opt.SetXYZ(fib31x, 0, 0);

                    // Up to this point the x and z positions of the fibers and mwpc are being placed in the lab frame 
                    TransformPoint(mwpc1_lab_opt, mwpc1_angles_opt + mwpc1_ang_offset_opt, mwpc1_position_opt + mwpc1_pos_offset_opt);
                    TransformPoint( fi32_lab_opt,   fi32_angles_opt + fi32_ang_offset_opt,   fi32_position_opt + fi32_pos_offset_opt);
                    TransformPoint( fi33_lab_opt,   fi33_angles_opt + fi33_ang_offset_opt,   fi33_position_opt + fi33_pos_offset_opt);
                    TransformPoint( fi31_lab_opt,   fi31_angles_opt + fi31_ang_offset_opt,   fi31_position_opt + fi31_pos_offset_opt);

                    if ( (hit_31 == false && hit_33 == true) || (hit_31 == true && hit_33 == true) ){

                        // Up to this point the x and z positions of the fibers and mwpc are being placed in the lab frame 
                        TransformPoint(mwpc1_lab_opt, mwpc1_angles_opt + mwpc1_ang_offset_opt, mwpc1_position_opt + mwpc1_pos_offset_opt);
                        TransformPoint( fi32_lab_opt,   fi32_angles_opt + fi32_ang_offset_opt,   fi32_position_opt + fi32_pos_offset_opt);
                        TransformPoint( fi33_lab_opt,   fi33_angles_opt + fi33_ang_offset_opt,   fi33_position_opt + fi33_pos_offset_opt);

                        PoQ =  CalculatePoQ(mwpc1_lab_opt, fi33_lab_opt, fi32_lab_opt, Actual_Current, Reference_Current, MDF_PoQ);
                        cout << "WEIRD VALUE: " << CalculatePoQ(mwpc1_lab_opt, fi31_lab_opt, fi32_lab_opt, Actual_Current, Reference_Current, MDF_PoQ) << endl; 
                    }

                    if (hit_31 == true && hit_33 == false){

                        // Up to this point the x and z positions of the fibers and mwpc are being placed in the lab frame 
                        TransformPoint(mwpc1_lab_opt, mwpc1_angles_opt + mwpc1_ang_offset_opt, mwpc1_position_opt + mwpc1_pos_offset_opt);
                        TransformPoint( fi32_lab_opt,   fi32_angles_opt + fi32_ang_offset_opt,   fi32_position_opt + fi32_pos_offset_opt);
                        TransformPoint( fi31_lab_opt,   fi31_angles_opt + fi31_ang_offset_opt,   fi31_position_opt + fi31_pos_offset_opt);

                        PoQ =  CalculatePoQ(mwpc1_lab_opt, fi31_lab_opt, fi32_lab_opt, Actual_Current, Reference_Current, MDF_PoQ);
                        cout << "WEIRD VALUE: " << CalculatePoQ(mwpc1_lab_opt, fi31_lab_opt, fi32_lab_opt, Actual_Current, Reference_Current, MDF_PoQ) << endl; 

                    }

                    std::cout << "\nFinal PoQ values with optimized parameters:\n";
                    std::cout << "PoQ"<< hit << " = " << PoQ << " | Ref = " << Reference_PoQ_several[hit] << " | Diff = " << PoQ - Reference_PoQ_several[hit] << "\n";
                }

                std::cout << "Best parameters: \n";
                for (int i = 0; i < 14; ++i) {
                    std::cout << "double par" << (i+1) << " = " << result[i] << ";" << "\n";
                }
                std::cout << "Best delta: " << minimizer->MinValue() << std::endl;
                
            }
        }
    }

    TCanvas* c_histos = new TCanvas();
    c_histos->Divide(4,4);
    for(int i = 0; i < 10; i++){
        c_histos->cd(i+1);
        h_pars[i]->Draw("hist");
    }
}