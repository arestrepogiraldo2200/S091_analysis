#include <array>
#include <iostream>

void plot_MWPC1_by_isotope() {
    gStyle->SetOptStat(0);

    const std::array<TString, 5> fileNames = {
        "/shared/storage/physnp/gmw537/s118/MDF/main0326_0001_unpacked_MDF_v3.root", // 10C
        "/shared/storage/physnp/gmw537/s118/MDF/main0321_0001_unpacked_MDF_v3.root", // 11C
        "/shared/storage/physnp/gmw537/s118/MDF/main0327_0001_unpacked_MDF_v3.root", // 12C
        "/shared/storage/physnp/gmw537/s118/MDF/main0324_0001_unpacked_MDF_v3.root", // 13C
        "/shared/storage/physnp/gmw537/s118/MDF/main0325_0001_unpacked_MDF_v3.root"  // 14C
    };

    const std::array<TString, 5> isotopeLabels = {
        "^{10}C",
        "^{11}C",
        "^{12}C",
        "^{13}C",
        "^{14}C"
    };

    const std::array<Int_t, 5> colors = {kRed + 1, kAzure + 1, kGreen + 2, kOrange + 7, kMagenta + 1};

    std::array<TH1D*, 5> hMwpc1X = {nullptr, nullptr, nullptr, nullptr, nullptr};

    for (size_t i = 0; i < fileNames.size(); ++i) {
        hMwpc1X[i] = new TH1D(
            Form("hMwpc1X_%zu", i),
            Form("MWPC1 x position for %s;MWPC1 x [mm];Counts", isotopeLabels[i].Data()),
            200,
            -150.0,
            0.0);
        hMwpc1X[i]->SetLineColor(colors[i]);
        hMwpc1X[i]->SetLineWidth(2);
    }

    for (size_t i = 0; i < fileNames.size(); ++i) {
        TChain chain("evt");
        chain.Add(fileNames[i]);

        TClonesArray* hitMWPC1 = new TClonesArray("R3BMwpcHitData", 5);
        chain.SetBranchAddress("Mwpc1HitData", &hitMWPC1);

        const Long64_t nEntries = chain.GetEntries();
        Long64_t nFilled = 0;

        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            hitMWPC1->Clear();
            chain.GetEntry(entry);

            const Long64_t nHits = hitMWPC1->GetEntriesFast();
            for (Long64_t hitIdx = 0; hitIdx < nHits; ++hitIdx) {
                auto* hit = static_cast<R3BMwpcHitData*>(hitMWPC1->At(hitIdx));
                if (!hit) {
                    continue;
                }

                const double x = hit->GetX();
                if (!std::isfinite(x)) {
                    continue;
                }

                hMwpc1X[i]->Fill(x);
                ++nFilled;
            }
        }

        std::cout << "[INFO] " << isotopeLabels[i] << " MWPC1 mean x = "
                  << hMwpc1X[i]->GetMean()
                  << " from " << nFilled << " hit(s)" << std::endl;
    }

    TCanvas* cOverlay = new TCanvas("cOverlay", "MWPC1 x by isotope", 1200, 800);
    hMwpc1X[0]->SetTitle("MWPC1 x position by isotope;MWPC1 x [mm];Counts");
    hMwpc1X[0]->SetMaximum(1.15 * std::max_element(hMwpc1X.begin(), hMwpc1X.end(),
        [](TH1D* a, TH1D* b) { return a->GetMaximum() < b->GetMaximum(); })->GetMaximum());
    hMwpc1X[0]->Draw("hist");
    for (size_t i = 1; i < hMwpc1X.size(); ++i) {
        hMwpc1X[i]->Draw("hist same");
    }

    TLegend* legend = new TLegend(0.70, 0.68, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    for (size_t i = 0; i < hMwpc1X.size(); ++i) {
        legend->AddEntry(hMwpc1X[i], isotopeLabels[i], "l");
    }
    legend->Draw();

    TCanvas* cPads = new TCanvas("cPads", "MWPC1 x by isotope", 1400, 900);
    cPads->Divide(3, 2);
    for (size_t i = 0; i < hMwpc1X.size(); ++i) {
        cPads->cd(i + 1);
        hMwpc1X[i]->Draw("hist");
        TLatex label;
        label.SetNDC(true);
        label.SetTextSize(0.06);
        label.DrawLatex(0.18, 0.85, isotopeLabels[i]);
    }

    TFile* out = TFile::Open("MWPC1_by_isotope.root", "RECREATE");
    if (out && !out->IsZombie()) {
        for (auto* h : hMwpc1X) {
            h->Write();
        }
        cOverlay->Write();
        cPads->Write();
        out->Close();
    }
}