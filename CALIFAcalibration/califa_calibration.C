
void califa_calibration(){

    // Load file
    auto calibration_file = new TFile("califa_calibration.root");

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // 60Co 

    /////////////////////////////////////////////////////////////////////////////////////////
    // Get crystal key and set its range
    auto crystal4700 = (TH1F*) calibration_file->Get("crystal_4700_60Co");
    
    crystal4700->GetXaxis()->SetRangeUser(0,3500);

    TCanvas* c_crystal = new TCanvas("","Spectrum 60Co");
    c_crystal->cd(1);
    crystal4700->Draw();

    /////////////////////////////////////////////////////////////////////////////////////////
    // Use spectrum to estimate the background and search for peaks
    auto spectrum = new TSpectrum();
    // auto crystal4700_back = (TH1F*) spectrum->Background(crystal4700);

    spectrum->Search(crystal4700,4);     // Second argument is the width of the expected peak
    for (int i = 0; i < spectrum->GetNPeaks(); i++) cout << "Position of peak " << i+1 << ": " << spectrum->GetPositionX()[i] << "\n";

    /////////////////////////////////////////////////////////////////////////////////////////
    auto fitfunc = new TF1("fitfunc","[0]*TMath::Gaus(x,[1],[2],false) + [3]*TMath::Gaus(x,[4],[5],false) + expo(6)",65,140);
    fitfunc->SetParameters(1000,90.5,3, 1000,102.5,3, 15,-5);
    fitfunc->SetParLimits(5,0,3.5);                              // Limits the values of the parameters
    fitfunc->SetParLimits(2,0,3.5);
    TFitResultPtr fit_values = crystal4700->Fit(fitfunc,"RS");
    fitfunc->Draw("same");

    /////////////////////////////////////////////////////////////////////////////////////////
    auto gauss1 = new TF1("fitfunc","[0]*TMath::Gaus(x,[1],[2], false)", 0,140);
    gauss1->SetParameters( fit_values->Parameter(0), fit_values->Parameter(1), fit_values->Parameter(2));
   
    auto gauss2 = new TF1("fitfunc2","[0]*TMath::Gaus(x,[1],[2], false)", 0,140);
    gauss2->SetParameters( fit_values->Parameter(3), fit_values->Parameter(4), fit_values->Parameter(5));
    
    auto lin    = new TF1("fitfunc3","expo",  0,140);
    lin->SetParameters( fit_values->Parameter(6), fit_values->Parameter(7));

    gauss1->SetLineColor(kBlack);
    gauss2->SetLineColor(kBlack);
    lin->SetLineColor(kBlack);

    gauss1->SetLineStyle(2);
    gauss2->SetLineStyle(2);
    lin->SetLineStyle(2);

    gauss1->Draw("same");
    gauss2->Draw("same");
    lin->Draw("same");

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // 22Na

    auto crystal4700Na = (TH1F*) calibration_file->Get("crystal_4700_22Na");
    crystal4700Na->GetXaxis()->SetRangeUser(0,3500);

    TCanvas* c_crystalNa = new TCanvas("","Spectrum 22Na");
    c_crystalNa->cd(1);
    crystal4700Na->Draw();

    auto spectrumNa = new TSpectrum();
    spectrumNa->Search(crystal4700Na, 4);

    auto fitfuncNa = new TF1("fitfuncNa","[0]*TMath::Gaus(x,[1],[2],false)", 90, 110);
    fitfuncNa->SetParameters(1200, 90.5, 3);
    TFitResultPtr fit_valuesNa = crystal4700Na->Fit(fitfuncNa,"RS");

    auto gauss3 = new TF1("fitfunc","[0]*TMath::Gaus(x,[1],[2],false)", 0, 140);
    gauss3->SetParameters(fit_valuesNa->Parameter(0), fit_valuesNa->Parameter(1), fit_valuesNa->Parameter(2));
    gauss3->SetLineColor(kBlack);
    gauss3->SetLineStyle(2);
    gauss3->Draw("same");

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // Calibration 

    double channels[4] = {0,fit_values->Parameter(1), fit_values->Parameter(4), fit_valuesNa->Parameter(1)};
    double energies[4] = {0,1173.228, 1332.492, 1274.537};

    cout << "\nChannel  " << "Energies" << "\n"; 
    for (int i = 0; i < 3; i++){
        cout << channels[i] << "  " << energies[i] << "\n";
    }

    auto calibration = new TGraph(3,channels,energies);
    calibration->Fit("pol1");

    TCanvas* c_calibration = new TCanvas("Channel vs. Energy","Channel vs. Energy");
    c_calibration->cd(1);
    calibration->SetMarkerSize(2);
    calibration->GetXaxis()->SetLimits(80,110);
    calibration->SetMarkerStyle(21);
    calibration->Draw("AP"); 


    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // Rescaling of the original signal


    // auto linear_calib = new TF1("myfunc", "pol1",0,150);
    // linear_calib->SetParameters(29.5924, 12.7102);

    // cout<<crystal4700->GetNbinsX()<<"\n";

    // auto calib_hist = new TH1D("Calibrated","Calibrated",10,0,1500);
    
    // for (int i = 1; i <= crystal4700->GetNbinsX(); ++i) {

    //     double x = crystal4700->GetBinCenter(i);
    //     cout << i << "  " << x << "\n"; 
    //     double scaled = linear_calib->Eval(x);
    //     calib_hist->Fill(scaled);

    // }

    // TCanvas* c_calibrated = new TCanvas("","Spectrum 60Co energy calibrated");
    // c_calibrated->cd(1);
    // calib_hist->Draw();

    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // ABORTED, ONLY THE TOTAL ENTRIES PER BIN ARE STORED IN THE TH1F OBJECT, THE ORIGINAL ENTRIES ARE LOST
    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    
    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // 22Na 511 keV peak

    auto calibration_file2 = new TFile("califa_calibration.root");

    auto crystal4700Na2 = (TH1F*) calibration_file2->Get("crystal_4700_22Na");
    crystal4700Na2->GetXaxis()->SetRangeUser(0,140);

    TCanvas* c_crystalNa2 = new TCanvas("","Spectrum 22Na 511 keV peak");
    c_crystalNa2->cd(1);
    crystal4700Na2->Draw();

    auto spectrumNa2 = new TSpectrum();
    spectrumNa2->Search(crystal4700Na2, 4);

    auto fitfuncNa2 = new TF1("fitfuncNa","[0]*TMath::Gaus(x,[1],[2],false) + [3]*TMath::Gaus(x,[4],[5],false) + [6]*TMath::Gaus(x,[7],[8],false)", 0, 110);
    fitfuncNa2->SetParameters(6000, 40, 2, 800, 70, 6 , 1200, 90.5, 3);
    TFitResultPtr fit_valuesNa2 = crystal4700Na2->Fit(fitfuncNa2,"RS");

    auto gauss4 = new TF1("fitfunc","[0]*TMath::Gaus(x,[1],[2],false)", 0, 140);
    gauss4->SetParameters(fit_valuesNa2->Parameter(0), fit_valuesNa2->Parameter(1), fit_valuesNa2->Parameter(2));
    gauss4->SetLineColor(kBlack);
    gauss4->SetLineStyle(2);
    gauss4->Draw("same");

    auto gauss5 = new TF1("fitfunc","[0]*TMath::Gaus(x,[1],[2],false)", 0, 140);
    gauss5->SetParameters(fit_valuesNa2->Parameter(3), fit_valuesNa2->Parameter(4), fit_valuesNa2->Parameter(5));
    gauss5->SetLineColor(kBlack);
    gauss5->SetLineStyle(2);
    gauss5->Draw("same");

    auto gauss6 = new TF1("fitfunc","[0]*TMath::Gaus(x,[1],[2],false)", 0, 140);
    gauss6->SetParameters(fit_valuesNa2->Parameter(6), fit_valuesNa2->Parameter(7), fit_valuesNa2->Parameter(8));
    gauss6->SetLineColor(kBlack);
    gauss6->SetLineStyle(2);
    gauss6->Draw("same");


    /////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////
    // Calibration 4 points

    double channels4[5] = {0,fit_values->Parameter(1), fit_values->Parameter(4), fit_valuesNa2->Parameter(1), fit_valuesNa2->Parameter(7)};
    double energies4[5] = {0,1173.228, 1332.492, 510.99895069, 1274.537};

    cout << "\nChannel  " << "Energies" << "\n"; 
    for (int i = 0; i < 3; i++){
        cout << channels4[i] << "  " << energies4[i] << "\n";
    }

    auto calibration4 = new TGraph(4,channels4,energies4);
    calibration4->Fit("pol1");

    TCanvas* c_calibration4 = new TCanvas("Channel vs. Energy 4 points","Channel vs. Energy 4 points");
    c_calibration4->cd(1);
    calibration4->SetMarkerSize(2);
    calibration4->GetXaxis()->SetLimits(20,110);
    calibration4->SetMarkerStyle(21);
    calibration4->Draw("AP"); 



}