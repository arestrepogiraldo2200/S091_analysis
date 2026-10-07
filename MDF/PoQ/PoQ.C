
#include <iostream>
#include <memory>
#include <vector>
#include <algorithm>
#include <cmath>

namespace {

void PrintPoQResolution(TH1D* hist, const TString& label)
{
    if (!hist || hist->GetEntries() <= 0) {
        std::cout << "[WARN] Cannot estimate resolution for " << label << " because the histogram is empty." << std::endl;
        return;
    }

    const Int_t peakBin = hist->GetMaximumBin();
    const Double_t peakX = hist->GetXaxis()->GetBinCenter(peakBin);
    const Double_t fitHalfWidth = 0.15;
    const Double_t fitMin = std::max(hist->GetXaxis()->GetXmin(), peakX - fitHalfWidth);
    const Double_t fitMax = std::min(hist->GetXaxis()->GetXmax(), peakX + fitHalfWidth);

    TF1 gauss(Form("gauss_%s", hist->GetName()), "gaus", fitMin, fitMax);
    const int fitStatus = hist->Fit(&gauss, "QNR");
    if (fitStatus != 0) {
        std::cout << "[WARN] Gaussian fit failed for " << label << " (status " << fitStatus << ")" << std::endl;
        return;
    }

    const Double_t mean = gauss.GetParameter(1);
    const Double_t sigma = std::abs(gauss.GetParameter(2));
    if (mean == 0.0) {
        std::cout << "[WARN] Gaussian mean is zero for " << label << ", cannot compute DeltaP/P." << std::endl;
        return;
    }

    const Double_t fwhm = 2.0 * std::sqrt(2.0 * std::log(2.0)) * sigma;
    const Double_t deltaPOverP = fwhm / std::abs(mean);
    std::cout << "[INFO] PoQ resolution for " << label
              << ": DeltaP/P (FWHM/mean) = " << deltaPOverP
              << " (" << deltaPOverP * 100.0 << "%)"
              << ", mean = " << mean
              << ", sigma = " << sigma
              << ", FWHM = " << fwhm
              << std::endl;
}

void MergeHistogramFromFiles(
    const std::vector<TString>& inputFiles,
    const TString& histName,
    const TString& objectPrefix,
    const TString& outFileName,
    const TString& canvasName,
    const TString& canvasTitle,
    const TString& overlayTitle,
    const TString& xAxisTitle,
    const Int_t targetNBins,
    const Double_t targetXMin,
    const Double_t targetXMax,
    const Int_t* colors,
    const Int_t nColors,
    const std::vector<Double_t>* verticalMarkers = nullptr,
    const Double_t xScale = 1.0)
{
    TH1D* hMerged = nullptr;
    std::vector<TH1D*> hPerFile;
    std::vector<TString> labels;
    Int_t nUsed = 0;
    Double_t yMax = 0.0;
    const Double_t markerLineAlpha = 0.35;

    for (size_t i = 0; i < inputFiles.size(); ++i) {
        const auto& fileName = inputFiles[i];
        TFile* f = TFile::Open(fileName, "READ");
        if (!f || f->IsZombie()) {
            std::cout << "[WARN] Cannot open file: " << fileName << std::endl;
            if (f) {
                f->Close();
            }
            continue;
        }

        TH1D* h = dynamic_cast<TH1D*>(f->Get(histName));
        if (!h) {
            std::cout << "[WARN] Histogram '" << histName << "' not found in " << fileName << std::endl;
            f->Close();
            continue;
        }

        TString fileLabel = fileName;
        if (fileLabel.Contains("10C")) {
            fileLabel = "^{10}C";
        } else if (fileLabel.Contains("11C")) {
            fileLabel = "^{11}C";
        } else if (fileLabel.Contains("12C")) {
            fileLabel = "^{12}C";
        } else if (fileLabel.Contains("13C")) {
            fileLabel = "^{13}C";
        } else if (fileLabel.Contains("14C")) {
            fileLabel = "^{14}C";
        }

        TH1D* hClone = new TH1D(Form("%s_%zu", objectPrefix.Data(), i), h->GetTitle(), targetNBins, targetXMin, targetXMax);
        hClone->SetDirectory(nullptr);
        hClone->Sumw2();

        for (Int_t srcBin = 1; srcBin <= h->GetNbinsX(); ++srcBin) {
            const Double_t srcLow = h->GetXaxis()->GetBinLowEdge(srcBin) * xScale;
            const Double_t srcWidth = h->GetXaxis()->GetBinWidth(srcBin) * xScale;
            const Double_t srcUp = srcLow + srcWidth;
            const Double_t srcContent = h->GetBinContent(srcBin);
            const Double_t srcError = h->GetBinError(srcBin);

            if (srcWidth <= 0.0 || (srcContent == 0.0 && srcError == 0.0)) {
                continue;
            }

            const Double_t clipLow = std::max(srcLow, targetXMin);
            const Double_t clipUp = std::min(srcUp, targetXMax);
            if (clipUp <= clipLow) {
                continue;
            }

            const Int_t firstTargetBin = hClone->GetXaxis()->FindBin(clipLow + 1e-12);
            const Int_t lastTargetBin = hClone->GetXaxis()->FindBin(clipUp - 1e-12);

            for (Int_t dstBin = firstTargetBin; dstBin <= lastTargetBin; ++dstBin) {
                const Double_t dstLow = hClone->GetXaxis()->GetBinLowEdge(dstBin);
                const Double_t dstUp = dstLow + hClone->GetXaxis()->GetBinWidth(dstBin);
                const Double_t overlap = std::min(clipUp, dstUp) - std::max(clipLow, dstLow);
                if (overlap <= 0.0) {
                    continue;
                }

                const Double_t fraction = overlap / srcWidth;
                const Double_t addContent = srcContent * fraction;
                const Double_t addError = srcError * fraction;
                const Double_t oldContent = hClone->GetBinContent(dstBin);
                const Double_t oldError = hClone->GetBinError(dstBin);

                hClone->SetBinContent(dstBin, oldContent + addContent);
                hClone->SetBinError(dstBin, std::sqrt(oldError * oldError + addError * addError));
            }
        }

        hClone->SetEntries(h->GetEntries());
        hClone->SetLineColor(colors[i % nColors]);
        hClone->SetLineWidth(2);
        hClone->SetFillStyle(0);

        std::cout << "[INFO] Mean of " << fileLabel << " " << objectPrefix << " distribution: "
              << hClone->GetMean() << std::endl;
          PrintPoQResolution(hClone, fileLabel + " " + objectPrefix);

        hPerFile.push_back(hClone);
        labels.push_back(fileLabel);
        if (hClone->GetMaximum() > yMax) {
            yMax = hClone->GetMaximum();
        }

        if (!hMerged) {
            hMerged = dynamic_cast<TH1D*>(hClone->Clone(Form("%s_merged", objectPrefix.Data())));
            hMerged->SetDirectory(nullptr);
            hMerged->SetTitle(overlayTitle);
        } else {
            hMerged->Add(hClone);
        }

        ++nUsed;
        f->Close();
    }

    if (!hMerged) {
        std::cout << "[ERROR] No valid histogram was found to merge for " << histName << "." << std::endl;
        return;
    }

    std::cout << "[INFO] Mean of merged " << objectPrefix << " distribution: "
              << hMerged->GetMean() << std::endl;

    TCanvas* cOverlay = nullptr;
    if (!hPerFile.empty()) {
        cOverlay = new TCanvas(canvasName, canvasTitle, 1200, 800);
        hPerFile[0]->SetTitle(overlayTitle + ";" + xAxisTitle + ";Counts");
        hPerFile[0]->GetXaxis()->SetRangeUser(targetXMin, targetXMax);
        hPerFile[0]->SetMaximum(yMax * 1.15);
        hPerFile[0]->Draw("hist");
        for (size_t i = 1; i < hPerFile.size(); ++i) {
            hPerFile[i]->Draw("hist same");
        }

        std::vector<TLine*> markerLines;
        if (verticalMarkers) {
            const Double_t markerYMax = yMax * 1.12;
            const size_t nMarkers = std::min(verticalMarkers->size(), hPerFile.size());
            markerLines.reserve(nMarkers);
            for (size_t i = 0; i < nMarkers; ++i) {
                const Double_t xValue = verticalMarkers->at(i) * xScale;
                TLine* line = new TLine(xValue, 0.0, xValue, markerYMax);
                line->SetLineColorAlpha(colors[i % nColors], markerLineAlpha);
                line->SetLineWidth(3);
                line->SetLineStyle(7);
                line->Draw("same");
                markerLines.push_back(line);
            }
        }

        TLegend* legend = new TLegend(0.12, 0.12, 0.30, 0.32);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextSize(0.032);

        for (size_t i = 0; i < hPerFile.size(); ++i) {
            legend->AddEntry(hPerFile[i], labels[i], "l");
        }
        legend->Draw();
    }

    TFile* outFile = TFile::Open(outFileName, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cout << "[ERROR] Cannot create output file: " << outFileName << std::endl;
        if (outFile) {
            outFile->Close();
        }
        return;
    }

    hMerged->Write();
    for (auto* h : hPerFile) {
        h->Write();
    }
    if (cOverlay) {
        cOverlay->Write();
    }
    outFile->Close();

    std::cout << "Merged histogram from " << nUsed << " file(s)." << std::endl;
    std::cout << "Output written to: " << outFileName << std::endl;
}

} // namespace

void PoQ() {

    // // Input files to merge 
    // const std::vector<TString> inputFiles = {
    //     "./10C_n_mdf.root",
    //     "./11C_n_mdf.root",
    //     "./12C_n_mdf.root",
    //     "./13C_n_mdf.root",
    //     "./14C_n_mdf.root"
    // };

    // Input files to merge 
    const std::vector<TString> inputFiles = {
        "./10C_opt_mdf.root",
        "./11C_opt_mdf.root",
        "./12C_opt_mdf.root",
        "./13C_opt_mdf.root",
        "./14C_opt_mdf.root"
    };

    const Int_t targetNBins = 1000;
    const Double_t targetXMin = 0.0;
    const Double_t targetXMax = 4;
    const Int_t colors[] = {kRed, kAzure, kGreen, kOrange, kMagenta};
    const Int_t nColors = sizeof(colors) / sizeof(colors[0]);
    constexpr Double_t kRadToDeg = 180.0 / TMath::Pi();

    // Edit these values to place one vertical marker per isotope in the same order as inputFiles.
    // Leave the vector empty if you do not want markers.
    const std::vector<Double_t> poqExperimentalValues = {
        1.527872,
        1.7092484,
        1.9902324,
        2.0717421,
        2.2532898
    };
    
    MergeHistogramFromFiles(
        inputFiles,
        "PoQ from MDF",
        "PoQ",
        "./PoQ_merged.root",
        "cPoQ_overlay",
        "PoQ by isotope",
        "PoQ from MDF by isotope",
        "PoQ [GeV/c]",
        targetNBins,
        targetXMin,
        targetXMax,
        colors,
        nColors,
        &poqExperimentalValues);

    const Int_t thetaNBins = 200;
    MergeHistogramFromFiles(
        inputFiles,
        "Incidence angle between plane of FIB32 and FIB31-FIB33",
        "IncidenceAngle",
        "./IncidenceAngle_merged.root",
        "cIncidenceAngle_overlay",
        "Incidence angle by isotope",
        "Incidence angle between plane of FIB32 and FIB31-FIB33 by isotope",
        "Theta [deg]",
        thetaNBins,
        -0.1 * kRadToDeg,
        0.1 * kRadToDeg,
        colors,
        nColors,
        nullptr,
        kRadToDeg);
}
