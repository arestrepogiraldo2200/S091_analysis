#include "Math/Functor.h"
// #include "Minuit2/Minuit2Minimizer.h"
#include "function.h"

char histname[200];
TH1D* h_pars[14];
// TH2D* h_229vs202;

void minimizePars() {
    
    // MDF_PoQ10 = new R3BMDFWrapper("/home/lukas/software/r3broot/R3BRoot/macros/r3b/Analyse/S473/MDF/output/MDF_params_poq_9pars200Terms180Fib10.txt");
    // MDF_PoQ = new R3BMDFWrapper("/home/lukas/software/r3broot/R3BRoot/macros/r3b/Analyse/S473/MDF/output/MDF_params_poq_9pars200TermsUncert180.txt");
    // MDF_PoQ1012 = new R3BMDFWrapper("/home/lukas/software/r3broot/R3BRoot/macros/r3b/Analyse/S473/MDF/output/MDF_params_poq_10parUncert2.txt");

    MDF_PoQ = new R3BMDFWrapper("./MDF_params.txt");

    const UInt_t num_pars = 14;
    const double target_delta = 0.5;

    // // Limits for the parameters 
    // double par_min[num_pars]   = { 0.  ,  0.  ,  -0.015,  -2.9,  -2.0,     0.0,   0.985,   0.985,  0.985,  0.985 };
    // double par_max[num_pars]   = { 0.01,  0.01,   0.015,  -2.8,  2.01,  0.0001,   1.015,   1.015,  1.015,  1.015 };        
    // double step_size[num_pars] = { 0.01,  0.01,  0.0001,  0.01,  0.01,  0.0001,  0.0001,  0.0001, 0.0001,  0.0001};


    // Limits for the parameters       MWPC                 FIB32                         FIB33                       FIB31
    double par_min[num_pars]   = {  -5,   -5,     -5,    -5,  -0.5,  -0.5,      -5,    -5,  -0.5,  -0.5,      -5,    -5,  -0.5,   -0.5};
    double par_max[num_pars]   = {   5,    5,      5,     5,   0.5,   0.5,       5,     5,   0.5,   0.5,       5,     5,   0.5,    0.5};    
    double step_size[num_pars] = {0.05, 0.05,   0.05,  0.05,  0.01,  0.01,    0.05,  0.05,  0.01,  0.01,    0.05,  0.05,  0.01,   0.01};


    // Histogram of parameters is created
    for(int i = 0; i < num_pars; i++){
        sprintf(histname, "hist par%i",i);
        h_pars[i] = new TH1D(histname, histname, 100, par_min[i], par_max[i]);   
    }
    
    // h_229vs202 = new TH2D("229vs202","229 vs 202",1000,3.43,3.47,1000,2.19,2.230);


    int num = 0;
    while(num<1000){

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

        ROOT::Math::Functor f(&MinimizeDelta4, num_pars, "./MDF_params.txt");
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

        // 
        if (min_val < target_delta) {

            const double* result = minimizer->X();

            for (int i = 0; i < num_pars; i++)
            {
                if ( fabs((result[i] - par_min[i]) / par_min[i]) < 0.1 || fabs((par_max[i] - result[i]) / par_max[i]) < 0.1){
                    break;
                }
            }

            success = true;

            double par1 = result[0];
            double par2 = result[1];
            double par3 = result[2];
            double par4 = result[3];
            double par5 = result[4];
            double par6 = result[5];
            double par7 = result[6];
            double par8 = result[7];
            double par9 = result[8];
            double par10 = result[9];
            double par11 = result[10];
            double par12 = result[11];
            double par13 = result[12];
            double par14 = result[13];

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


            double pq_227 = CalculatePoQ1012(fi10_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ1012);


            std::cout << "\nFinal PoQ values with optimized parameters:\n";
            std::cout << "PoQ 227 = " << pq_227 << " | Ref = " << Reference_PoQ_227 << " | Diff = " << pq_227 - Reference_PoQ_227 << "\n";


            // fi11_lab_227.SetXYZ(-1.*(-4.60759 ) , 0 , 0 );
            // fi12_lab_227.SetXYZ(-1.*(-16.8725 + par5) , 0 , 0 );
            // fi10_lab_227.SetXYZ(-1.*(-5.20989 + par2) , 0 , 0 );  

            // fi11_lab_223.SetXYZ(-1.*(-3.20910 ) , 0 , 0 );
            // fi12_lab_223.SetXYZ(-1.*(-11.2588 + par5) , 0 , 0 );         
            // fi10_lab_223.SetXYZ(-1.*(-3.20786 + par2) , 0 , 0 ); 

            // fi11_lab_225.SetXYZ(-1.*(-3.64903 ) , 0 , 0 );
            // fi12_lab_225.SetXYZ(-1.*(-13.0462 + par5) , 0 , 0 );
            // fi10_lab_225.SetXYZ(-1.*(-3.84497 + par2) , 0 , 0 ); 
            
            // fi11_lab_229.SetXYZ(-1.*(-2.73996 ) , 0 , 0 );
            // fi12_lab_229.SetXYZ(-1.*(-9.26699 + par5) , 0 , 0 );
            // fi10_lab_229.SetXYZ(-1.*(-2.53834 + par2) , 0 , 0 ); 
            
            // fi11_lab_212.SetXYZ(-1.*(-2.83415 ) , 0 , 0 );
            // fi12_lab_212.SetXYZ(-1.*(-9.74823 + par5) , 0 , 0 ); 
            // fi10_lab_212.SetXYZ(-1.*(-2.68217 + par2) , 0 , 0 );
                            
            // fi11_lab_220.SetXYZ(-1.*(-2.83545 ) , 0 , 0 );
            // fi12_lab_220.SetXYZ(-1.*(-10.0403 + par5) , 0 , 0 );
            // fi10_lab_220.SetXYZ(-1.*(-2.71734 + par2) , 0 , 0 );
            
            // fi11_lab_201.SetXYZ(-1.*(-2.93167 ) , 0 , 0 );
            // fi12_lab_201.SetXYZ(-1.*(-10.3565 + par5) , 0 , 0 );
            // fi10_lab_201.SetXYZ(-1.*(-2.85664 + par2) , 0 , 0 );
            
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
            
            // double pq_227 = CalculatePoQ1012(fi10_lab_227, fi12_lab_227, psp2_lab, psp3_lab, fi13_lab, Actual_Current_227, Reference_Current, MDF_PoQ1012);
            // double pq_223 = CalculatePoQ1012(fi10_lab_223, fi12_lab_223, psp2_lab, psp3_lab, fi13_lab, Actual_Current_223, Reference_Current, MDF_PoQ1012);
            // double pq_225 = CalculatePoQ1012(fi10_lab_225, fi12_lab_225, psp2_lab, psp3_lab, fi13_lab, Actual_Current_225, Reference_Current, MDF_PoQ1012);
            // double pq_229 = CalculatePoQ1012(fi10_lab_229, fi12_lab_229, psp2_lab, psp3_lab, fi13_lab, Actual_Current_229, Reference_Current, MDF_PoQ1012);
            // double pq_212 = CalculatePoQ1012(fi10_lab_212, fi12_lab_212, psp2_lab, psp3_lab, fi13_lab, Actual_Current_212, Reference_Current, MDF_PoQ1012);
            // double pq_220 = CalculatePoQ1012(fi10_lab_220, fi12_lab_220, psp2_lab, psp3_lab, fi13_lab, Actual_Current_220, Reference_Current, MDF_PoQ1012);
            // double pq_201 = CalculatePoQ1012(fi10_lab_201, fi12_lab_201, psp2_lab, psp3_lab, fi13_lab, Actual_Current_201, Reference_Current, MDF_PoQ1012);
            // double pq_202 = CalculatePoQ1012(fi10_lab_202, fi12_lab_202, psp2_lab, psp3_lab, fi13_lab, Actual_Current_202, Reference_Current, MDF_PoQ1012);

            // // 🖨️ Ausgabe:
            // std::cout << "\nFinal PoQ values with optimized parameters:\n";
            // std::cout << "PoQ 227 = " << pq_227 << " | Ref = " << Reference_PoQ_227 << " | Diff = " << pq_227 - Reference_PoQ_227 << "\n";
            // std::cout << "PoQ 223 = " << pq_223 << " | Ref = " << Reference_PoQ_223 << " | Diff = " << pq_223 - Reference_PoQ_223 << "\n";
            // std::cout << "PoQ 225 = " << pq_225 << " | Ref = " << Reference_PoQ_225 << " | Diff = " << pq_225 - Reference_PoQ_225 << "\n";
            // std::cout << "PoQ 229 = " << pq_229 << " | Ref = " << Reference_PoQ_229 << " | Diff = " << pq_229 - Reference_PoQ_229 << "\n";
            // std::cout << "PoQ 212 = " << pq_212 << " | Ref = " << Reference_PoQ_212 << " | Diff = " << pq_212 - Reference_PoQ_212 << "\n";
            // std::cout << "PoQ 220 = " << pq_220 << " | Ref = " << Reference_PoQ_220 << " | Diff = " << pq_220 - Reference_PoQ_220 << "\n";
            // std::cout << "PoQ 201 = " << pq_201 << " | Ref = " << Reference_PoQ_201 << " | Diff = " << pq_201 - Reference_PoQ_201 << "\n";
            // std::cout << "PoQ 202 = " << pq_202 << " | Ref = " << Reference_PoQ_202 << " | Diff = " << pq_202 - Reference_PoQ_202 << "\n";
                

            // h_229vs202->Fill(pq_229,pq_202);
            std::cout << "Best parameters: \n";

            for (int i = 0; i < 14; ++i) {
                std::cout << "double par" << (i+1) << " = " << result[i] << ";" << "\n";
            }
            std::cout << "Best delta: " << minimizer->MinValue() << std::endl;
            
    }
     // delete minimizer; // WICHTIG: Aufräumen!
    }
    }

    //std::cout << "\nOptimization completed after " << attempts << " attempt(s)." << std::endl;
    TCanvas* c_histos = new TCanvas();
    c_histos->Divide(4,4);
    for(int i = 0; i < 10; i++){
        c_histos->cd(i+1);
        h_pars[i]->Draw("hist");
    }
    TCanvas* c_compare = new TCanvas();
    h_229vs202->Draw("COLZ");
    
}