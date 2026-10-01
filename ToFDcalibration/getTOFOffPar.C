
#include <iostream>
#include <vector>
#include "TFile.h"
#include "TDirectory.h"
#include "TH2D.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TSpectrum.h"
#include "TGraph.h"
#include "TF1.h"

// Funktion zum Laden eines Histograms
TH1D* loadHistogram(TFile* file, int i, int j) {

    const char* dirName = Form("hist_TOFD-Plane%d-ToF-Bar",i+1);
    const char* histName = Form("TOFD Plane%i - Bar%i - ToF Bar", i + 1, j + 1);
    //const char* histName = Form("TOFD Plane1 - Bar1 - ToT vs ToF Corr");
    TDirectory* dir = file->GetDirectory(dirName);
    if (!dir) {
        std::cerr << "Verzeichnis nicht gefunden: " << dirName << std::endl;
        file->Close();
        return nullptr;
    }
    dir->cd();    
    std::cout << "Ordner erfolgreich geöffnet: " << dirName << std::endl;
    //dir->ls();
    TH1D* hist = dynamic_cast<TH1D*>(dir->Get(histName));
    if (!hist) {
        std::cerr << "TH2D Histogramm nicht gefunden: " << histName << std::endl;
        return nullptr;
    }

    std::cout << "TH2D Histogramm erfolgreich geladen: " << histName << std::endl;
    return hist;
}

// Hauptfunktion
void getTOFOffPar() {
    double TOFPar[4][44];
    const char* fileName = "./Hists/oldHists/hists_Analysis_0197_Sep24_Zcorr_hist.root";

    TFile* file = TFile::Open(fileName);
    if (!file || file->IsZombie()) {
        std::cerr << "Fehler beim Öffnen der Datei: " << fileName << std::endl;
        return;
    }
    std::cout << "Datei erfolgreich geöffnet: " << fileName << std::endl;

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 44; j++) {
            // Lade das Histogramm
            TH1D* hist = loadHistogram(file, i, j);
            if (!hist) continue;

            //TString name = Form("ProjY Plane%d - Bar%d", i, j);
            //TH1D* hist = hist2D->ProjectionY(name, 0, 1000);


            // TSpectrum zur Peak-Suche
            TSpectrum* s = new TSpectrum(1, 0.5);
            int nfound = s->Search(hist, 1.0, "goff", 0.7);
            if (nfound > 0) {
                double* xpeak = s->GetPositionX();
                std::cout << "Gefundener Peak bei: " << xpeak[0] << std::endl;
                TF1* fgauss = new TF1("fgaus", "gaus(0)", xpeak[0] - 0.05, xpeak[0] + 0.05);
                hist->Fit("fgaus", "QR");
                auto tof_offset = fgauss->GetParameter(1);
                //TOFPar[i][j] = tof_offset;
                TOFPar[i][j] = xpeak[0];
            }
            delete s; // Speicher freigeben
        }
    }
    file->Close();
    
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 44; j++){
            char OutName[100];
            sprintf(OutName, "./parameter_tof_offset_197");
            ofstream savefile(OutName, ios_base::app);
            
            savefile << i+1 << " " << j+1 << " 1 10 " <<  TOFPar[i][j] << endl;
        }
    }            
    
}
getTOFOffPar.C
Displaying getTOFOffPar.C.lukas.ponnath@york.ac.uk