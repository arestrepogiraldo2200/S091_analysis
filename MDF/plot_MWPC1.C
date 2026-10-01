#include <algorithm>
#include <array>
#include <cctype>
#include <iostream>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TFile.h"
#include "TIterator.h"
#include "TList.h"
#include "TKey.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TLatex.h"
#include "TSystem.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"

namespace {

const TString kDirectorySuffix = "C_mdf";
const int kMinDirectoryNumber = 10;
const int kMaxDirectoryNumber = 14;

bool IsSelectedDirectory(const TString& name, const TString& suffix) {
    if (!name.EndsWith(suffix)) {
        return false;
    }

    const int prefixLength = name.Length() - suffix.Length();
    if (prefixLength <= 0) {
        return false;
    }

    TString prefix(name(0, prefixLength));
    if (!std::all_of(prefix.Data(), prefix.Data() + prefix.Length(), [](unsigned char c) {
        return std::isdigit(c) != 0;
    })) {
        return false;
    }

    const int directoryNumber = std::stoi(prefix.Data());
    return directoryNumber >= kMinDirectoryNumber && directoryNumber <= kMaxDirectoryNumber;
}

struct InputSample {
    TString label;
    TString filePath;
    Int_t color;
    Int_t order;
};

bool IsDirectoryCandidate(const TSystemFile* file) {
    if (!file || !file->IsDirectory()) {
        return false;
    }

    const TString name = file->GetName();
    return IsSelectedDirectory(name, kDirectorySuffix) && name != "." && name != "..";
}

InputSample MakeSample(const TString& baseDir, const TString& dirName) {
    InputSample sample;
    sample.filePath = Form("%s/%s/%s.root", baseDir.Data(), dirName.Data(), dirName.Data());
    sample.label = dirName;
    sample.color = kGray + 2;
    sample.order = 999;

    if (dirName.Contains("10C")) {
        sample.label = "^{10}C";
        sample.color = kRed;
        sample.order = 10;
    } else if (dirName.Contains("11C")) {
        sample.label = "^{11}C";
        sample.color = kAzure;
        sample.order = 11;
    } else if (dirName.Contains("12C")) {
        sample.label = "^{12}C";
        sample.color = kGreen;
        sample.order = 12;
    } else if (dirName.Contains("13C")) {
        sample.label = "^{13}C";
        sample.color = kOrange;
        sample.order = 13;
    } else if (dirName.Contains("14C")) {
        sample.label = "^{14}C";
        sample.color = kMagenta;
        sample.order = 14;
    }

    return sample;
}

std::vector<InputSample> CollectSamples(const TString& baseDir) {
    std::vector<InputSample> samples;

    TSystemDirectory directory("mdf_scan", baseDir);
    TList* entries = directory.GetListOfFiles();
    if (!entries) {
        return samples;
    }

    TIter next(entries);
    while (auto* obj = next()) {
        auto* file = dynamic_cast<TSystemFile*>(obj);
        if (!IsDirectoryCandidate(file)) {
            continue;
        }

        const TString dirName = file->GetName();
        const InputSample sample = MakeSample(baseDir, dirName);
        if (gSystem->AccessPathName(sample.filePath)) {
            std::cout << "[WARN] Missing file: " << sample.filePath << std::endl;
            continue;
        }

        samples.push_back(sample);
    }

    std::sort(samples.begin(), samples.end(), [](const InputSample& a, const InputSample& b) {
        if (a.order != b.order) {
            return a.order < b.order;
        }
        return a.label < b.label;
    });

    return samples;
}

TH1D* CloneHistogram(TH1D* source, const TString& name) {
    if (!source) {
        return nullptr;
    }

    auto* clone = dynamic_cast<TH1D*>(source->Clone(name));
    if (clone) {
        clone->SetDirectory(nullptr);
    }
    return clone;
}

} // namespace

void plot_MWPC1(const char* baseDir = ".") {
    gStyle->SetOptStat(0);

    const TString histName = "MWPC1 x position";
    const TString outFileName = "MWPC1_from_mdf_dirs.root";

    const std::vector<InputSample> samples = CollectSamples(baseDir);
    if (samples.empty()) {
        std::cout << "[ERROR] No 10C_mdf to 14C_mdf directories with matching ROOT files were found in "
                  << baseDir << std::endl;
        return;
    }

    std::vector<TH1D*> perSampleHistograms;
    std::vector<TString> labels;
    TH1D* mergedHistogram = nullptr;
    Double_t maxY = 0.0;

    for (size_t i = 0; i < samples.size(); ++i) {
        const auto& sample = samples[i];
        TFile* inputFile = TFile::Open(sample.filePath, "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cout << "[WARN] Cannot open file: " << sample.filePath << std::endl;
            if (inputFile) {
                inputFile->Close();
            }
            continue;
        }

        auto* source = dynamic_cast<TH1D*>(inputFile->Get(histName));
        if (!source) {
            std::cout << "[WARN] Histogram '" << histName << "' not found in " << sample.filePath << std::endl;
            inputFile->Close();
            continue;
        }

        TH1D* hist = CloneHistogram(source, Form("hMWPC1_%zu", i));
        if (!hist) {
            std::cout << "[WARN] Could not clone histogram from " << sample.filePath << std::endl;
            inputFile->Close();
            continue;
        }

        hist->SetTitle(Form("MWPC1 x position for %s;MWPC1 x [mm];Counts", sample.label.Data()));
        hist->SetLineColor(sample.color);
        hist->SetLineWidth(2);
        hist->SetFillStyle(0);

        perSampleHistograms.push_back(hist);
        labels.push_back(sample.label);

        if (hist->GetMaximum() > maxY) {
            maxY = hist->GetMaximum();
        }

        if (!mergedHistogram) {
            mergedHistogram = CloneHistogram(hist, "hMWPC1_merged");
            if (mergedHistogram) {
                mergedHistogram->SetTitle("Merged MWPC1 x position;MWPC1 x [mm];Counts");
            }
        } else {
            mergedHistogram->Add(hist);
        }

        std::cout << "[INFO] " << sample.label << " mean MWPC1 x = " << hist->GetMean()
                  << " from " << hist->GetEntries() << " entries" << std::endl;

        inputFile->Close();
    }

    if (perSampleHistograms.empty() || !mergedHistogram) {
        std::cout << "[ERROR] No valid MWPC1 histograms were found to merge." << std::endl;
        return;
    }

    TCanvas* cOverlay = new TCanvas("cMWPC1_overlay", "MWPC1 x by isotope", 1200, 800);
    perSampleHistograms.front()->SetTitle("MWPC1 x position by isotope;MWPC1 x [mm];Counts");
    perSampleHistograms.front()->GetXaxis()->SetRangeUser(-150.0, 0.0);
    perSampleHistograms.front()->SetMaximum(maxY * 1.15);
    perSampleHistograms.front()->Draw("hist");

    for (size_t i = 1; i < perSampleHistograms.size(); ++i) {
        perSampleHistograms[i]->Draw("hist same");
    }

    TLegend* legend = new TLegend(0.70, 0.68, 0.88, 0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);
    for (size_t i = 0; i < perSampleHistograms.size(); ++i) {
        legend->AddEntry(perSampleHistograms[i], labels[i], "l");
    }
    legend->Draw();

    TCanvas* cMerged = new TCanvas("cMWPC1_merged", "Merged MWPC1 x position", 1200, 800);
    mergedHistogram->Draw("hist");

    TFile* outFile = TFile::Open(outFileName, "RECREATE");
    if (!outFile || outFile->IsZombie()) {
        std::cout << "[ERROR] Cannot create output file: " << outFileName << std::endl;
        if (outFile) {
            outFile->Close();
        }
        return;
    }

    mergedHistogram->Write();
    for (auto* hist : perSampleHistograms) {
        hist->Write();
    }
    cOverlay->Write();
    cMerged->Write();
    outFile->Close();

    std::cout << "[INFO] Wrote merged histogram and overlay canvas to " << outFileName << std::endl;
}