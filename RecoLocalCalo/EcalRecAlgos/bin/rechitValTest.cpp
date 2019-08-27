#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TPaveStats.h>

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"

int main(int argc, char *argv[]) {

// *** File setup to test EE; It can be used to test EB just by changing all EE to EB. ***
  
// Set the pointers to GPU and CPU pointers for both EB and EE. 
    edm::Wrapper<ecal::RecHit<ecal::Tag::soa>> *wgpuEE=nullptr;
    edm::Wrapper<EERecHitCollection> *wcpuEE = nullptr;

// The input file containing the data to be validated (i.e. result.root)
    std::string fileName = argv[1]; 

// output file setup
    std::ofstream test_ee("test_ee.txt", std::ios_base::app | std::ios_base::out);

// input file setup for tree
    TFile rf{fileName.c_str()};
    TTree *rt = (TTree*)rf.Get("Events");
// Allocating the appropriate data to their respective pointers
    rt->SetBranchAddress("ecalTagsoaecalRecHit_ecalCPURecHitProducer_EcalRecHitsEE_RECO.", &wgpuEE);
    rt->SetBranchAddress("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_RECO.", &wcpuEE);

// accumulate sizes for events and sizes of each event on both GPU and CPU
    auto const nentries = rt->GetEntries();
    for (int ie=0; ie<nentries; ++ie) {
        rt->GetEntry(ie);
	
    // Check size of EE in both CPU and GPU:
        auto cpu_ee_size = wcpuEE->bareProduct().size();
        auto gpu_ee_size = wgpuEE->bareProduct().did.size();
	test_ee<<"cpu size:\t"<<cpu_ee_size<<"\t gpu size:\t"<< gpu_ee_size<<std::endl;
	
    // Check EE GPU components:
	for (uint32_t i=0; i<gpu_ee_size; ++i) {
 
	auto const did_gpu = wgpuEE->bareProduct().did[i];  
        auto const ener_gpu = wgpuEE->bareProduct().energy[i];
	auto const chi2_gpu = wgpuEE->bareProduct().chi2[i]; 
	auto const flag_gpu = wgpuEE->bareProduct().flagBits[i]; 
	auto const extra_gpu = wgpuEE->bareProduct().extra[i]; 
	
	// Choose investigated component(s) to print: 
        if (did_gpu != 0) {
	test_ee<<"Event: \t"<<ie<<"\t rechit: \t"<<i <<std::endl; 
	test_ee<<"******** EE GPU DetID:*********** \t"<< did_gpu<<std::endl;
	}
	else {
	test_ee<<"DetID is 0"<<std::endl;
       }
      }
     }
    

    // Close all open files
    rf.Close();

    return 0;
}
