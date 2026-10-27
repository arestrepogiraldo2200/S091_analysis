
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <vector>

char fname1[200], fname2[200], fname3[200], fname4[200], fname5[200];

// Sentinel: set to false to draw the MWPC1->Fib32 leg as a straight line instead of an arc
const bool kUseArcMwpc1ToFib32 = true;
const Double_t kArcSagitta = 20.0; // bulge of the arc away from the straight chord
const Int_t kArcPoints = 30;      // number of points used to sample the arc

// Keeps an angle within +-pi of a reference angle, so a monotonic sweep can be built
Double_t UnwrapAngle(Double_t angle, Double_t ref) {
    while (angle - ref > TMath::Pi()) angle -= TMath::TwoPi();
    while (angle - ref < -TMath::Pi()) angle += TMath::TwoPi();
    return angle;
}

// Samples a circular arc from (z0,x0) to (z1,x1) that bulges by `sagitta` at its midpoint
void ComputeArcPoints(Double_t z0, Double_t x0, Double_t z1, Double_t x1, Double_t sagitta, Int_t nPoints,
                       std::vector<Double_t>& zOut, std::vector<Double_t>& xOut) {
    const Double_t dz = z1 - z0;
    const Double_t dx = x1 - x0;
    const Double_t chordLen = std::sqrt(dz * dz + dx * dx);

    if (chordLen <= 0 || std::fabs(sagitta) < 1e-9) {
        zOut = {z0, z1};
        xOut = {x0, x1};
        return;
    }

    const Double_t halfChord = 0.5 * chordLen;
    const Double_t radius = (halfChord * halfChord + sagitta * sagitta) / (2.0 * sagitta);

    const Double_t midZ = 0.5 * (z0 + z1);
    const Double_t midX = 0.5 * (x0 + x1);
    const Double_t perpZ = -dx / chordLen;
    const Double_t perpX = dz / chordLen;

    const Double_t centerZ = midZ - (radius - sagitta) * perpZ;
    const Double_t centerX = midX - (radius - sagitta) * perpX;
    const Double_t apexZ = midZ + sagitta * perpZ;
    const Double_t apexX = midX + sagitta * perpX;

    const Double_t angle0 = std::atan2(x0 - centerX, z0 - centerZ);
    const Double_t angleApex = UnwrapAngle(std::atan2(apexX - centerX, apexZ - centerZ), angle0);
    const Double_t angle1 = UnwrapAngle(std::atan2(x1 - centerX, z1 - centerZ), angleApex);

    zOut.resize(nPoints);
    xOut.resize(nPoints);
    for (Int_t k = 0; k < nPoints; ++k) {
        const Double_t t = static_cast<Double_t>(k) / (nPoints - 1);
        const Double_t theta = angle0 + t * (angle1 - angle0);
        zOut[k] = centerZ + radius * std::cos(theta);
        xOut[k] = centerX + radius * std::sin(theta);
    }
}

void Propagation_MDF_training(){

    // Load the root files /////////////////////////////////////////////////
    TChain* chain = new TChain("detectors");

    const std::array<TString, 5> isotopeLabels = {"^{10}C", "^{11}C", "^{12}C", "^{13}C", "^{14}C"};
    const std::array<Int_t, 5> isotopeColors = {kRed, kAzure, kGreen, kOrange, kMagenta};
    std::array<TH1D*, 5> hMwpc1PerIsotope = {nullptr, nullptr, nullptr, nullptr, nullptr};
    std::array<TH1D*, 5> hFib32PerIsotope = {nullptr, nullptr, nullptr, nullptr, nullptr};
    std::array<TH1D*, 5> hFib3133PerIsotope = {nullptr, nullptr, nullptr, nullptr, nullptr};

    // sprintf(fname1,"./10C_n_lab.root");
    // sprintf(fname2,"./11C_n_lab.root");
    // sprintf(fname3,"./12C_n_lab.root");
    // sprintf(fname4,"./13C_n_lab.root");
    // sprintf(fname5,"./14C_n_lab.root");

    sprintf(fname1,"./10C_lab.root");
    sprintf(fname2,"./11C_lab.root");
    sprintf(fname3,"./12C_lab.root");
    sprintf(fname4,"./13C_lab.root");
    sprintf(fname5,"./14C_lab.root");

    chain->Reset();
    chain->Add(fname1);
    chain->Add(fname2);
    chain->Add(fname3);
    chain->Add(fname4);
    chain->Add(fname5);

    ////////////////////////////////////////////////////////////////////////

    // Struct to load the branches data
    struct DetectorPosition {
        Double_t x;
        Double_t z;
    };

    DetectorPosition mwpc1{};
    DetectorPosition fib32{};
    DetectorPosition fib3133{};

    chain->SetBranchAddress("MWPC1", &mwpc1);
    chain->SetBranchAddress("Fib32", &fib32);
    chain->SetBranchAddress("Fib31_Fib33", &fib3133);

    for (size_t i = 0; i < hMwpc1PerIsotope.size(); ++i) {
        hMwpc1PerIsotope[i] = new TH1D(
            Form("hMwpc1_%zu", i),
            Form("MWPC1 x positions;MWPC1 x [cm];Counts"),
            150,
            -10,
            +10);
        hMwpc1PerIsotope[i]->SetLineColor(isotopeColors[i]);
        hMwpc1PerIsotope[i]->SetLineWidth(2);

        hFib32PerIsotope[i] = new TH1D(
            Form("hFib32_%zu", i),
            Form("Fib32 x positions;Fib32 x [cm];Counts"),
            400,
            -200,
            -50);
        hFib32PerIsotope[i]->SetLineColor(isotopeColors[i]);
        hFib32PerIsotope[i]->SetLineWidth(2);

        hFib3133PerIsotope[i] = new TH1D(
            Form("hFib3133_%zu", i),
            Form("Fib31_Fib33 x positions;Fib31_Fib33 x [cm];Counts"),
            400,
            -200,
            -50);
        hFib3133PerIsotope[i]->SetLineColor(isotopeColors[i]);
        hFib3133PerIsotope[i]->SetLineWidth(2);
    }
    
    ////////////////////////////////////////////////////////////////////////

    // Arrays to store data
    Long64_t nEntries = chain->GetEntries();

    // Maximum number of points in the trajectories
    const Int_t maxPointsPerFile = 2000;                             // <----------------------------

    std::vector<std::array<Double_t, 3>> allZ;
    std::vector<std::array<Double_t, 3>> allX;
    std::vector<Int_t> allFileIdx;
    std::map<Int_t, Int_t> pointsPerFile;

    allZ.reserve(nEntries);
    allX.reserve(nEntries);
    allFileIdx.reserve(nEntries);

    Double_t zMin = std::numeric_limits<Double_t>::max();
    Double_t zMax = std::numeric_limits<Double_t>::lowest();
    Double_t xMin = std::numeric_limits<Double_t>::max();
    Double_t xMax = std::numeric_limits<Double_t>::lowest();
    ////////////////////////////////////////////////////////////////////////

    for (Long64_t entry = 0; entry < nEntries; ++entry) {

        chain->GetEntry(entry);

        const Int_t fileIdx = chain->GetTreeNumber();
        if (pointsPerFile[fileIdx] >= maxPointsPerFile) {
            continue;
        }
        pointsPerFile[fileIdx]++;

        const std::array<Double_t, 3> zTrack = {mwpc1.z, fib32.z, fib3133.z};
        const std::array<Double_t, 3> xTrack = {mwpc1.x, fib32.x, fib3133.x};

        bool validTrack = true;
        for (Int_t i = 0; i < 3; ++i) {
            if (!std::isfinite(zTrack[i]) || !std::isfinite(xTrack[i])) {
                validTrack = false;
                break;
            }
        }
        if (!validTrack) {
            continue;
        }

        if (fileIdx >= 0 && fileIdx < static_cast<Int_t>(hMwpc1PerIsotope.size())) {
            hMwpc1PerIsotope[fileIdx]->Fill(mwpc1.x);
            hFib32PerIsotope[fileIdx]->Fill(fib32.x);
            hFib3133PerIsotope[fileIdx]->Fill(fib3133.x);
        }

        allZ.push_back(zTrack);
        allX.push_back(xTrack);

        // Gets the id of a file in the chain
        allFileIdx.push_back(fileIdx);

        // Set limits for plot ////////////////////////////////////////////////
        for (Int_t i = 0; i < 3; ++i) {
            if (zTrack[i] < zMin) zMin = zTrack[i];
            if (zTrack[i] > zMax) zMax = zTrack[i];
            if (xTrack[i] < xMin) xMin = xTrack[i];
            if (xTrack[i] > xMax) xMax = xTrack[i];
        }
    }

    TCanvas* c2 = new TCanvas("c2", "Trajectories", 1200, 900);

    // To adjust the frame limits ////////////////////////////////////////////
    const Double_t zRange = (zMax - zMin > 0 ? zMax - zMin : 1.0);
    const Double_t xRange = (xMax - xMin > 0 ? xMax - xMin : 1.0);

    // Add a wider margin so thick markers/lines stay visibly inside the frame.
    const Double_t zPad = 0.12 * zRange + 1.0;
    const Double_t xPad = 0.12 * xRange + 1.0;

    TH2F* hFrame = new TH2F(
        "hFrame",
        "Fragments trajectories in the X-Z plane;Z;X",
        100,
        zMin - zPad,
        zMax + zPad,
        100,
        xMin - xPad,
        xMax + xPad
    );
    hFrame->SetStats(0);
    hFrame->Draw();
    ////////////////////////////////////////////////////////////////////////

    const Int_t fileColors[] = {kRed , kAzure , kGreen , kOrange , kMagenta };
    const Int_t nFileColors = sizeof(fileColors) / sizeof(fileColors[0]);
    std::map<Int_t, TPolyLine*> firstLinePerFile;

    for (size_t i = 0; i < allZ.size(); ++i) {

        // Sets the colors according to the file 
        const Int_t fileIdx = allFileIdx[i];
        const Int_t color = fileColors[fileIdx % nFileColors];

        TPolyLine* legendLine = nullptr;

        if (kUseArcMwpc1ToFib32) {
            // MWPC1 -> Fib32 (arc)
            std::vector<Double_t> zArc, xArc;
            ComputeArcPoints(allZ[i][0], allX[i][0], allZ[i][1], allX[i][1], kArcSagitta, kArcPoints, zArc, xArc);
            TPolyLine* arc = new TPolyLine(zArc.size(), zArc.data(), xArc.data());
            arc->SetLineColorAlpha(color, 0.05);
            arc->SetLineWidth(1);
            arc->Draw("same");
            legendLine = arc;

            // Fib32 -> Fib31_Fib33 (straight)
            std::array<Double_t, 2> zSeg2 = {allZ[i][1], allZ[i][2]};
            std::array<Double_t, 2> xSeg2 = {allX[i][1], allX[i][2]};
            TPolyLine* seg2 = new TPolyLine(2, zSeg2.data(), xSeg2.data());
            seg2->SetLineColorAlpha(color, 0.05);
            seg2->SetLineWidth(1);
            seg2->Draw("same");
        } else {
            // Draws a line between the x-z coordinates of the detectors
            TPolyLine* line = new TPolyLine(3, allZ[i].data(), allX[i].data());
            line->SetLineColorAlpha(color, 0.05);
            line->SetLineWidth(1);
            line->Draw("same");
            legendLine = line;
        }

        if (firstLinePerFile.find(fileIdx) == firstLinePerFile.end()) {
            firstLinePerFile[fileIdx] = legendLine;
        }

        // Draws the coordinates
        TPolyMarker* points = new TPolyMarker(3, allZ[i].data(), allX[i].data());
        points->SetMarkerStyle(20);
        points->SetMarkerColor(kBlack);
        points->SetMarkerSize(1.3);
        points->Draw("same");
    }

    // Plot the legend
    TLegend* legend = new TLegend(0.12, 0.12, 0.38, 0.28);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextSize(0.03);

    TObjArray* files = chain->GetListOfFiles();
    for (const auto& entry : firstLinePerFile) {
        const Int_t idx = entry.first;
        const TChainElement* fileElement = dynamic_cast<TChainElement*>(files->At(idx));
        TString fileLabel = fileElement ? fileElement->GetTitle() : Form("file %d", idx);

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

        legend->AddEntry(entry.second, fileLabel, "l");
    }
    legend->Draw();

    TCanvas* cMwpc1 = new TCanvas("cMwpc1", "MWPC1 x positions per isotope", 1200, 900);
    cMwpc1->SetGrid();

    Double_t maxY = 0.0;
    for (const auto* hist : hMwpc1PerIsotope) {
        if (hist && hist->GetMaximum() > maxY) {
            maxY = hist->GetMaximum();
        }
    }

    bool firstHist = true;
    for (size_t i = 0; i < hMwpc1PerIsotope.size(); ++i) {
        if (!hMwpc1PerIsotope[i]) {
            continue;
        }

        hMwpc1PerIsotope[i]->SetMaximum(maxY > 0.0 ? 1.15 * maxY : 1.0);
        if (firstHist) {
            hMwpc1PerIsotope[i]->Draw("hist");
            firstHist = false;
        } else {
            hMwpc1PerIsotope[i]->Draw("hist same");
        }
    }

    TLegend* mwpc1Legend = new TLegend(0.62, 0.68, 0.88, 0.88);
    mwpc1Legend->SetBorderSize(0);
    mwpc1Legend->SetFillStyle(0);
    for (size_t i = 0; i < hMwpc1PerIsotope.size(); ++i) {
        if (hMwpc1PerIsotope[i]) {
            mwpc1Legend->AddEntry(hMwpc1PerIsotope[i], isotopeLabels[i], "l");
        }
    }
    mwpc1Legend->Draw();

    TCanvas* cFib32 = new TCanvas("cFib32", "Fib32 x positions per isotope", 1200, 900);
    cFib32->SetGrid();
    firstHist = true;
    maxY = 0.0;
    for (const auto* hist : hFib32PerIsotope) {
        if (hist && hist->GetMaximum() > maxY) {
            maxY = hist->GetMaximum();
        }
    }
    for (size_t i = 0; i < hFib32PerIsotope.size(); ++i) {
        if (!hFib32PerIsotope[i]) {
            continue;
        }

        hFib32PerIsotope[i]->SetMaximum(maxY > 0.0 ? 1.15 * maxY : 1.0);
        if (firstHist) {
            hFib32PerIsotope[i]->Draw("hist");
            firstHist = false;
        } else {
            hFib32PerIsotope[i]->Draw("hist same");
        }
    }

    TLegend* fib32Legend = new TLegend(0.62, 0.68, 0.88, 0.88);
    fib32Legend->SetBorderSize(0);
    fib32Legend->SetFillStyle(0);
    for (size_t i = 0; i < hFib32PerIsotope.size(); ++i) {
        if (hFib32PerIsotope[i]) {
            fib32Legend->AddEntry(hFib32PerIsotope[i], isotopeLabels[i], "l");
        }
    }
    fib32Legend->Draw();

    TCanvas* cFib3133 = new TCanvas("cFib3133", "Fib31_Fib33 x positions per isotope", 1200, 900);
    cFib3133->SetGrid();
    firstHist = true;
    maxY = 0.0;
    for (const auto* hist : hFib3133PerIsotope) {
        if (hist && hist->GetMaximum() > maxY) {
            maxY = hist->GetMaximum();
        }
    }
    for (size_t i = 0; i < hFib3133PerIsotope.size(); ++i) {
        if (!hFib3133PerIsotope[i]) {
            continue;
        }

        hFib3133PerIsotope[i]->SetMaximum(maxY > 0.0 ? 1.15 * maxY : 1.0);
        if (firstHist) {
            hFib3133PerIsotope[i]->Draw("hist");
            firstHist = false;
        } else {
            hFib3133PerIsotope[i]->Draw("hist same");
        }
    }

    TLegend* fib3133Legend = new TLegend(0.62, 0.68, 0.88, 0.88);
    fib3133Legend->SetBorderSize(0);
    fib3133Legend->SetFillStyle(0);
    for (size_t i = 0; i < hFib3133PerIsotope.size(); ++i) {
        if (hFib3133PerIsotope[i]) {
            fib3133Legend->AddEntry(hFib3133PerIsotope[i], isotopeLabels[i], "l");
        }
    }
    fib3133Legend->Draw();

    c2->Update();
    c2->Draw();

    cMwpc1->Update();
    cMwpc1->Draw();

    cFib32->Update();
    cFib32->Draw();

    cFib3133->Update();
    cFib3133->Draw();
}
