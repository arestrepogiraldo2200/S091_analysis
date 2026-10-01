// merge_root_files.C
// ROOT macro to merge multiple ROOT files containing the same TTree

#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <iostream>

using namespace std;

void merger()
{
    // Input file names
    const char* fname1 = "./main0320_0001_unpacked_MDF_v3.root";
    const char* fname2 = "./main0321_0001_unpacked_MDF_v3.root";
    const char* fname3 = "./main0323_0001_unpacked_MDF_v3.root";
    const char* fname4 = "./main0324_0001_unpacked_MDF_v3.root";
    const char* fname5 = "./main0325_0001_unpacked_MDF_v3.root";
    const char* fname6 = "./main0326_0001_unpacked_MDF_v3.root";
    const char* fname7 = "./main0327_0001_unpacked_MDF_v3.root";

    // Name of the tree inside the ROOT files
    TChain* chain = new TChain("evt");

    // Add files to chain
    // chain->Add(fname1);
    chain->Add(fname2);
    // chain->Add(fname3);
    chain->Add(fname4);
    chain->Add(fname5);
    chain->Add(fname6);
    chain->Add(fname7);

    cout << "Total entries in chain: "
         << chain->GetEntries() << endl;

    // Create output file
    TFile* outfile = new TFile("merged_MDF_runs_v3.root", "RECREATE");

    // Clone all entries into a new tree
    TTree* newtree = chain->CloneTree(-1, "fast");

    // Write tree
    newtree->Write();

    // Close file
    outfile->Close();

    cout << "Merged tree saved to merged_MDF_runs_v3.root" << endl;
}