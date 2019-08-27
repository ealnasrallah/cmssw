#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include <TCanvas.h>
#include <TStyle.h>
#include <TPad.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TPaveStats.h>

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"

int main(int argc, char *argv[]) {
    if (argc<3) {
        std::cout << "run with: ./validateGPU <path to input file> <output file>\n";
        exit(0);
    }
    
  // Set the pointers to GPU and CPU pointers for both EB and EE. 
    edm::Wrapper<ecal::RecHit<ecal::Tag::soa>> *wgpuEB=nullptr;
    edm::Wrapper<ecal::RecHit<ecal::Tag::soa>> *wgpuEE=nullptr;
    edm::Wrapper<EBRecHitCollection> *wcpuEB = nullptr;
    edm::Wrapper<EERecHitCollection> *wcpuEE = nullptr;

  // Set the input and output files
    std::string fileName = argv[1]; // The input file containing the data to be validated (i.e. result.root)
    std::string outFileName = argv[2]; //The output file in which the validation results will be saved (i.e. output.root)

  // output file setup
    TFile rfout{outFileName.c_str(), "recreate"};

    int nbins = 300;
    float last = 2.;
    
  // Histograms for the reconstructed energies for EB and EE for each event on both GPU and CPU
    auto hSOIEnergiesEBGPU = new TH1D("hSOIEnergiesEBGPU", "hSOIEnergiesEBGPU; Energy [GeV]", nbins, 0, last);
    auto hSOIEnergiesEEGPU = new TH1D("hSOIEnergiesEEGPU", "hSOIEnergiesEEGPU; Energy [GeV]", nbins, 0, last);
    auto hSOIEnergiesEBCPU = new TH1D("hSOIEnergiesEBCPU", "hSOIEnergiesEBCPU; Energy [GeV]", nbins, 0, last);
    auto hSOIEnergiesEECPU = new TH1D("hSOIEnergiesEECPU", "hSOIEnergiesEECPU; Energy [GeV]", nbins, 0, last);
    auto hSOIEnergiesEBGPUvsCPU = new TH2D("hSOIEnergiesEBGPUvsCPU", "hSOIEnergiesEBGPUvsCPU; CPU; GPU", nbins, 0, last, nbins, 0, last);
    auto hSOIEnergiesEEGPUvsCPU = new TH2D("hSOIEnergiesEEGPUvsCPU", "hSOIEnergiesEEGPUvsCPU; CPU; GPU", nbins, 0, last, nbins, 0, last);
    auto hSOIEnergiesEBGPUCPUratio = new TH1D("EnergiesEBGPU/CPUratio", "EnergiesEBGPU/CPUraio; GPU/CPU", 100, 0, 4);
    auto hSOIEnergiesEEGPUCPUratio = new TH1D("EnergiesEEGPU/CPUratio", "EnergiesEEGPU/CPUratio; GPU/CPU", 100, 0, 1);
    auto hSOIEnergiesEEGPUCPUratio_eta1 = new TH1D("EnergiesEEGPU/CPUratio_eta1", "(0-20 Eta)EnergiesEEGPU/CPUratio; GPU/CPU", 100, 0, 1.5);
    auto hSOIEnergiesEEGPUCPUratio_eta2 = new TH1D("EnergiesEEGPU/CPUratio_eta2", "(20-35 Eta)EnergiesEEGPU/CPUratio; GPU/CPU", 100, 0, 1.5);
    auto hSOIEnergiesEEGPUCPUratio_eta3 = new TH1D("EnergiesEEGPU/CPUratio_eta3", "(35-49 Eta)EnergiesEEGPU/CPUratio; GPU/CPU", 100, 0, 1.5);


  // input file setup for tree
    std::cout << "validating file " << fileName << std::endl;
    TFile rf{fileName.c_str()};
    TTree *rt = (TTree*)rf.Get("Events");
    
  // Allocating the appropriate data to their respective pointers
    rt->SetBranchAddress("ecalTagsoaecalRecHit_ecalCPURecHitProducer_EcalRecHitsEB_RECO.", &wgpuEB);
    rt->SetBranchAddress("ecalTagsoaecalRecHit_ecalCPURecHitProducer_EcalRecHitsEE_RECO.", &wgpuEE);
    rt->SetBranchAddress("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.", &wcpuEB);
    rt->SetBranchAddress("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_RECO.", &wcpuEE);

  // accumulate sizes for events and sizes of each event on both GPU and CPU
    auto const nentries = rt->GetEntries();
    std::cout << "#events to validate over: " << nentries << std::endl;
    for (int ie=0; ie<nentries; ++ie) {
        rt->GetEntry(ie);
        
        auto cpu_eb_size = wcpuEB->bareProduct().size();
        auto cpu_ee_size = wcpuEE->bareProduct().size();
        auto gpu_eb_size = wgpuEB->bareProduct().energy.size();
        auto gpu_ee_size = wgpuEE->bareProduct().energy.size();
	
    /*    
     const char* ordinal[] = { "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th" };
     // condition that sizes on GPU and CPU should be the same for EB or EE
	if (cpu_eb_size != gpu_eb_size or cpu_ee_size != gpu_ee_size) {
          std::cerr << ie << ordinal[ie % 10] << " entry:\n"
                    << "  EB size: " << std::setw(4) << cpu_eb_size << " (cpu) vs " << std::setw(4) << gpu_eb_size << " (gpu)\n"
                    << "  EE size: " << std::setw(4) << cpu_ee_size << " (cpu) vs " << std::setw(4) << gpu_ee_size << " (gpu)" << std::endl;
		   
          continue;
        }
        assert(wgpuEB->bareProduct().energy.size() == wcpuEB->bareProduct().size());
        assert(wgpuEE->bareProduct().energy.size() == wcpuEE->bareProduct().size()); 
        auto const neb = wcpuEB->bareProduct().size(); //like cpu_eb_size but set to constant
        auto const nee = wcpuEE->bareProduct().size(); //like cpu_ee_size but set to constant
        
     */
	
    // EB:
        for (uint32_t i=0; i<gpu_eb_size; ++i) {
            auto const did_gpu_eb = wgpuEB->bareProduct().did[i];
            auto const soi_enr_gpu_eb = wgpuEB->bareProduct().energy[i];
	    auto const cpu_iter_eb = wcpuEB->bareProduct().find(DetId{did_gpu_eb});
	    
	  // In case of error failing to find the corresponding CPU input, print message and move on.
            if (cpu_iter_eb == wcpuEB->bareProduct().end()) {
             //   std::cerr << ie << ordinal[ie % 10] << " entry\n"
             //             << "  Did not find a DetId " << did_gpu
             //             << " in a CPU collection\n";
                continue;
            }
          // Set the CPU's energy
            auto const soi_enr_cpu_eb = cpu_iter_eb->energy();

	  // Fill the energy histograms for GPU and CPU
            hSOIEnergiesEBGPU->Fill(soi_enr_gpu_eb);
            hSOIEnergiesEBCPU->Fill(soi_enr_cpu_eb);
            hSOIEnergiesEBGPUvsCPU->Fill(soi_enr_cpu_eb, soi_enr_gpu_eb);
	    hSOIEnergiesEBGPUCPUratio->Fill(soi_enr_gpu_eb/soi_enr_cpu_eb);
        }
        
        
    // EE:
	      for (uint32_t i=0; i<gpu_ee_size; ++i) {
            auto const did_gpu_ee = wgpuEE->bareProduct().did[i]; 
            auto const soi_enr_gpu_ee = wgpuEE->bareProduct().energy[i];   
	    auto const cpu_iter_ee = wcpuEE->bareProduct().find(DetId{did_gpu_ee}); 
	    
	  // In case of error failing to find the corresponding CPU input, print message and move on.
            if (cpu_iter_ee == wcpuEE->bareProduct().end()) {
             //   std::cerr << ie << ordinal[ie % 10] << " entry\n"
             //             << "  Did not find a DetId " << did_gpu
             //             << " in a CPU collection\n";
                continue;
            }
          // Set the CPU's energy
            auto const soi_enr_cpu_ee = cpu_iter_ee->energy();

	  // Fill the energy and Chi2 histograms for GPU and CPU
            hSOIEnergiesEEGPU->Fill(soi_enr_gpu_ee);
            hSOIEnergiesEECPU->Fill(soi_enr_cpu_ee);
            hSOIEnergiesEEGPUvsCPU->Fill(soi_enr_cpu_ee, soi_enr_gpu_ee);
	    hSOIEnergiesEEGPUCPUratio->Fill(soi_enr_gpu_ee/soi_enr_cpu_ee);
	    
	  // ******* eta test for the endcap *******
	    
	    // Get x and y values:
	    auto const soi_x_cpu_ee = EEDetId{did_gpu_ee}.ix();
	    auto const soi_y_cpu_ee = EEDetId{did_gpu_ee}.iy();
	    
	    // Equation for calculating eta:
	    auto const soi_eta_cpu_ee = sqrt((soi_x_cpu_ee - 50.5f) * (soi_x_cpu_ee - 50.5f) + (soi_y_cpu_ee - 50.5f) * (soi_y_cpu_ee - 50.5f));
	    std::cout<<"eta is = \t"<<soi_eta_cpu_ee<<std::endl;
	    
	    // Plot 3 regions of eta:
	    if (soi_eta_cpu_ee >= 0 && soi_eta_cpu_ee < 20) {
	      hSOIEnergiesEEGPUCPUratio_eta1->Fill(soi_enr_gpu_ee/soi_enr_cpu_ee);
	    } else if (soi_eta_cpu_ee >= 20 && soi_eta_cpu_ee < 35) {
	      hSOIEnergiesEEGPUCPUratio_eta2->Fill(soi_enr_gpu_ee/soi_enr_cpu_ee);
	    } else if (soi_eta_cpu_ee >= 35 && soi_eta_cpu_ee <= 49) {
	      hSOIEnergiesEEGPUCPUratio_eta3->Fill(soi_enr_gpu_ee/soi_enr_cpu_ee);
	    }
        }
    }

    // Plotting the results:
    {
      // Canvas for values of energy
      TCanvas cEnergies("Energies", "Energies", 1750, 860);
      cEnergies.Divide(3, 2);
      //to monitor EE GPU/CPU ratio over eta
      TCanvas cEta("Eta", "Energies", 1750, 600);
      cEta.Divide(3, 1);
      
  // Energy
      cEnergies.cd(1);
      {
          gPad->SetLogy();
          hSOIEnergiesEBCPU->SetLineColor(kBlack);
          hSOIEnergiesEBCPU->SetLineWidth(1.);
          hSOIEnergiesEBCPU->Draw("");
          hSOIEnergiesEBGPU->SetLineColor(kBlue);
          hSOIEnergiesEBGPU->SetLineWidth(1.);
          hSOIEnergiesEBGPU->Draw("sames");
          gPad->Update();
          auto stats = (TPaveStats*)hSOIEnergiesEBGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
      }
      cEnergies.cd(4);
      {
          gPad->SetLogy();
          hSOIEnergiesEECPU->SetLineColor(kBlack);
          hSOIEnergiesEECPU->SetLineWidth(1.);
          hSOIEnergiesEECPU->Draw("");
          hSOIEnergiesEEGPU->SetLineColor(kBlue);
          hSOIEnergiesEEGPU->SetLineWidth(1.);
          hSOIEnergiesEEGPU->Draw("sames");
          gPad->Update();
          auto stats = (TPaveStats*)hSOIEnergiesEEGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
      }
      cEnergies.cd(2); {
      hSOIEnergiesEBGPUvsCPU->Draw("COLZ");
      }
      cEnergies.cd(5); {
      hSOIEnergiesEEGPUvsCPU->Draw("COLZ");
      }
      cEnergies.cd(3); {
      hSOIEnergiesEBGPUCPUratio->Draw("");
      gPad->SetLogy();
      }
      cEnergies.cd(6); {
      hSOIEnergiesEEGPUCPUratio->Draw("");
      gPad->SetLogy();
      }
      cEnergies.SaveAs("ecal-energies.root");

  // Energy for different eta's
      cEta.cd(1); {
      hSOIEnergiesEEGPUCPUratio_eta1->Draw("");
      gPad->SetLogy();
      }
      cEta.cd(2); {
      hSOIEnergiesEEGPUCPUratio_eta2->Draw("");
      gPad->SetLogy();
      }
      cEta.cd(3); {
      hSOIEnergiesEEGPUCPUratio_eta3->Draw("");
      gPad->SetLogy();
      }
      cEta.SaveAs("ecal-EE-eta.root");
      }
    
    
    // Close all open files
    rf.Close();
    rfout.Write();
    rfout.Close();

    return 0;
}
