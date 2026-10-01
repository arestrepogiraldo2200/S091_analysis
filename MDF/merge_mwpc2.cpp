// merge_mwpc2.C
// root -l -q merge_mwpc2.C

#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TString.h>
#include <iostream>

void merge_mwpc2()
{
    for (int run = 320; run <= 328; ++run)
    {
        TString base = Form("main%04d_0001_unpacked_MDF", run);

        TString origName = base + ".root";
        TString v2Name   = base + "_v2.root";
        TString outName  = base + "_v3.root";

        std::cout << "\nProcessing " << base << std::endl;

        TFile* fOrig = TFile::Open(origName);
        TFile* fV2   = TFile::Open(v2Name);

        if (!fOrig || !fV2) {
            std::cout << "Cannot open input files, skipping.\n";
            continue;
        }

        TTree* tOrig = (TTree*)fOrig->Get("evt");
        TTree* tV2   = (TTree*)fV2->Get("evt");

        if (!tOrig || !tV2) {
            std::cout << "Cannot find tree evt, skipping.\n";
            continue;
        }

        Long64_t nentries = tOrig->GetEntries();

        // Open output file
        TFile* fOut = new TFile(outName, "RECREATE");

        // Copy all existing branches from the original tree
        TTree* tOut = tOrig->CloneTree(-1, "fast");

        // Objects that will receive data from v2
        TClonesArray* mapped = nullptr;
        TClonesArray* cal    = nullptr;
        TClonesArray* hit    = nullptr;

        tV2->SetBranchAddress("Mwpc2MappedData", &mapped);
        tV2->SetBranchAddress("Mwpc2CalData",    &cal);
        tV2->SetBranchAddress("Mwpc2HitData",    &hit);

        // Create new branches in output tree
        TBranch* bMap =
            tOut->Branch("Mwpc2MappedData", &mapped, 32000, 99);

        TBranch* bCal =
            tOut->Branch("Mwpc2CalData", &cal, 32000, 99);

        TBranch* bHit =
            tOut->Branch("Mwpc2HitData", &hit, 32000, 99);

        // Fill the new branches entry by entry
        for (Long64_t i = 0; i < nentries; ++i) {
            tV2->GetEntry(i);

            bMap->Fill();
            bCal->Fill();
            bHit->Fill();
        }

        tOut->Write();
        fOut->Close();

        fOrig->Close();
        fV2->Close();

        std::cout << "Created " << outName << std::endl;
    }
}