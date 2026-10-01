// File: merge_first1000_safe.C
//
// Usage:
// root -l -b -q 'merge_first1000_safe.C("merged.root")'

#include <TFile.h>
#include <TTree.h>
#include <TList.h>
#include <TKey.h>
#include <TClass.h>

#include <vector>
#include <string>
#include <iostream>

TTree* FindFirstTree(TFile* f)
{
    TIter next(f->GetListOfKeys());
    TKey* key;

    while ((key = (TKey*)next())) {
        TClass* cl = gROOT->GetClass(key->GetClassName());

        if (cl && cl->InheritsFrom(TTree::Class())) {
            return (TTree*)key->ReadObj();
        }
    }

    return nullptr;
}

void merger_partial(const char* outputFile = "merged_partial.root")
{
    std::vector<std::string> files = {
        "./main0321_0001_unpacked_MDF.root",
        "./main0323_0001_unpacked_MDF.root",
        "./main0324_0001_unpacked_MDF.root",
        "./main0325_0001_unpacked_MDF.root",
        "./main0326_0001_unpacked_MDF.root",
        "./main0327_0001_unpacked_MDF.root"
    };

    TList trees;

    for (const auto& fname : files) {

        std::cout << "Opening: " << fname << std::endl;

        TFile* f = TFile::Open(fname.c_str());

        if (!f || f->IsZombie()) {
            std::cerr << "ERROR opening " << fname << std::endl;
            continue;
        }

        TTree* tree = FindFirstTree(f);

        if (!tree) {
            std::cerr << "ERROR: no tree found in "
                      << fname << std::endl;
            continue;
        }

        std::cout << "Tree: " << tree->GetName()
                  << " | entries = "
                  << tree->GetEntries()
                  << std::endl;

        // Copy only first 1000 entries
        TTree* smallTree = tree->CopyTree("", "", 5000, 0);

        trees.Add(smallTree);
    }

    TFile out(outputFile, "RECREATE");

    // Merge all reduced trees
    TTree* merged = TTree::MergeTrees(&trees);

    if (merged) {
        merged->Write();
        std::cout << "Merged tree written to "
                  << outputFile << std::endl;
    } else {
        std::cerr << "ERROR: merge failed"
                  << std::endl;
    }

    out.Close();
}