void Gaussians() {


    auto ideal_spec = new TF1("gauss6", "[2]*exp(-0.5*((x-3)/0.15)^2)+(x-3)**2*([0]*exp(-0.5*((x-1)/0.15)^2) +  [1]*exp(-0.5*((x-2)/0.15)^2) +  [2]*exp(-0.5*((x-3)/0.15)^2) +  [3]*exp(-0.5*((x-4)/0.15)^2) +  [4]*exp(-0.5*((x-5)/0.15)^2) +  [5]*exp(-0.5*((x-6)/0.15)^2))", 0, 7);
    ideal_spec->SetParameters(1000, 1000, 1000, 1000, 1000, 1000); 
    TCanvas* c_ideal = new TCanvas();
    c_ideal->cd(1);
    ideal_spec->Draw();


    auto gauss1 = new TF1("gauss1", "exp(-0.5*((x-3)/0.15)^2) + exp(-0.5*((x-3.4)/0.15)^2)", 0, 7);
    TCanvas* c_gauss1 = new TCanvas();
    gauss1->SetLineColor(kGreen);
    gauss1->Draw();
    
    auto gauss2 = new TF1("gauss2", "exp(-0.5*((x-2.8)/0.15)^2)+exp(-0.5*((x-4)/0.15)^2)", 0, 7);
    gauss2->Draw("same");

    auto gaussprod = new TF1("gaussprod", "gauss1*gauss2", 0, 7);
    gaussprod->SetLineColor(kBlue);
    gaussprod->Draw("same");
    

    auto gaussprod2 = new TF1("gaussprod2", "gaussprod*gauss2", 0, 7);

    double max_val = gaussprod2->GetMaximum(0, 7);
    gaussprod2->SetParameter(0, gaussprod2->GetParameter(0) / max_val);

    gaussprod2->SetLineColor(kYellow);
    gaussprod2->Draw("same");
    


}
