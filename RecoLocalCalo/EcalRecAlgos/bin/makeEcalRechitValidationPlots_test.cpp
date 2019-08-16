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
    if (argc<3) {
        std::cout << "run with: ./validateGPU <path to input file> <output file>\n";
        exit(0);
    }
// Set the pointers to GPU and CPU pointers for both EB and EE. 
    edm::Wrapper<ecal::RecHit<ecal::Tag::soa>> *wgpuEB=nullptr;
    edm::Wrapper<ecal::RecHit<ecal::Tag::soa>> *wgpuEE=nullptr;
    edm::Wrapper<EBRecHitCollection> *wcpuEB = nullptr;
    edm::Wrapper<EERecHitCollection> *wcpuEE = nullptr;

    std::string fileName = argv[1]; // The input file containing the data to be validated (i.e. result.root)
    std::string outFileName = argv[2]; //The output file in which the validation results will be saved (i.e. output.root)

    // output file setup
    // Eissa: to test output and debugging
    std::ofstream test_r("test_r.txt", std::ios_base::app | std::ios_base::out);
    
    TFile rfout{outFileName.c_str(), "recreate"};

    int nbins = 300;
    float last = 50.;

    int nbins_chi2 = 1000;
    float last_chi2 = 1000.;

    int nbins_delta = 201;  // use an odd number to center around 0
    float delta = 0.2;
    
    // Histograsm to show the sizes for EB and EE for each event on both GPU and CPU
    auto hSOIRechitsEBGPU = new TH1D("RechitsEBGPU", "RechitsEBGPU", 100, 0, 3000);
    auto hSOIRechitsEBCPU = new TH1D("RechitsEBCPU", "RechitsEBCPU", 100, 0, 3000);
    auto hSOIRechitsEEGPU = new TH1D("RechitsEEGPU", "RechitsEEGPU", 100, 0, 3000);
    auto hSOIRechitsEECPU = new TH1D("RechitsEECPU", "RechitsEECPU", 100, 0, 3000);
    
    // Histograms for the reconstructed energies for EB and EE for each event on both GPU and CPU
    auto hSOIEnergiesEBGPU = new TH1D("hSOIEnergiesEBGPU", "hSOIEnergiesEBGPU", nbins, 0, last);
    auto hSOIEnergiesEEGPU = new TH1D("hSOIEnergiesEEGPU", "hSOIEnergiesEEGPU", nbins, 0, last);
    auto hSOIEnergiesEBCPU = new TH1D("hSOIEnergiesEBCPU", "hSOIEnergiesEBCPU", nbins, 0, last);
    auto hSOIEnergiesEECPU = new TH1D("hSOIEnergiesEECPU", "hSOIEnergiesEECPU", nbins, 0, last);

    // Chi2 plot for EB and EE on both GPU and CPU
    auto hChi2EBGPU = new TH1D("hChi2EBGPU", "hChi2EBGPU", nbins_chi2, 0, last_chi2);
    auto hChi2EEGPU = new TH1D("hChi2EEGPU", "hChi2EEGPU", nbins_chi2, 0, last_chi2);
    auto hChi2EBCPU = new TH1D("hChi2EBCPU", "hChi2EBCPU", nbins_chi2, 0, last_chi2);
    auto hChi2EECPU = new TH1D("hChi2EECPU", "hChi2EECPU", nbins_chi2, 0, last_chi2);

    // 2D histograms to compare GPU and CPU, and delta and CPU energies for both EB and EE 
    auto hSOIEnergiesEBGPUvsCPU = new TH2D("hSOIEnergiesEBGPUvsCPU", "hSOIEnergiesEBGPUvsCPU", nbins, 0, last, nbins, 0, last);
    auto hSOIEnergiesEEGPUvsCPU = new TH2D("hSOIEnergiesEEGPUvsCPU", "hSOIEnergiesEEGPUvsCPU", nbins, 0, last, nbins, 0, last);
    auto hSOIEnergiesEBdeltavsCPU = new TH2D("hSOIEnergiesEBdeltavsCPU", "hSOIEnergiesEBdeltavsCPU", nbins, 0, last, nbins_delta, -delta, delta);
    auto hSOIEnergiesEEdeltavsCPU = new TH2D("hSOIEnergiesEEdeltavsCPU", "hSOIEnergiesEEdeltavsCPU", nbins, 0, last, nbins_delta, -delta, delta);

    // 2D histograms to compare GPU and CPU and delta and CPU Chi2 for both EB and EE
    auto hChi2EBGPUvsCPU = new TH2D("hChi2EBGPUvsCPU", "hChi2EBGPUvsCPU", nbins_chi2, 0, last_chi2, nbins_chi2, 0, last_chi2);
    auto hChi2EEGPUvsCPU = new TH2D("hChi2EEGPUvsCPU", "hChi2EEGPUvsCPU", nbins_chi2, 0, last_chi2, nbins_chi2, 0, last_chi2);
    auto hChi2EBdeltavsCPU = new TH2D("hChi2EBdeltavsCPU", "hChi2EBdeltavsCPU", nbins_chi2, 0, last_chi2, nbins_delta, -delta, delta);
    auto hChi2EEdeltavsCPU = new TH2D("hChi2EEdeltavsCPU", "hChi2EEdeltavsCPU", nbins_chi2, 0, last_chi2, nbins_delta, -delta, delta);

    // input file setup for tree
    std::cout << "validating file " << fileName << std::endl;
    TFile rf{fileName.c_str()};
    TTree *rt = (TTree*)rf.Get("Events");
    // Allocating the appropriate data to their respective pointers
    rt->SetBranchAddress("ecalTagsoaecalRecHit_ecalCPURecHitProducer_EcalRecHitsEB_RECO.", &wgpuEB);
    rt->SetBranchAddress("ecalTagsoaecalRecHit_ecalCPURecHitProducer_EcalRecHitsEE_RECO.", &wgpuEE);
    rt->SetBranchAddress("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEB_RECO.", &wcpuEB);
    rt->SetBranchAddress("EcalRecHitsSorted_ecalRecHit_EcalRecHitsEE_RECO.", &wcpuEE);

    constexpr float eps_diff = 1e-3;

    // accumulate sizes for events and sizes of each event on both GPU and CPU
    auto const nentries = rt->GetEntries();
    std::cout << "#events to validate over: " << nentries << std::endl;
    for (int ie=0; ie<nentries; ++ie) {
        rt->GetEntry(ie);

        const char* ordinal[] = { "th", "st", "nd", "rd", "th", "th", "th", "th", "th", "th" };
        auto cpu_eb_size = wcpuEB->bareProduct().size();
        auto cpu_ee_size = wcpuEE->bareProduct().size();
        auto gpu_eb_size = wgpuEB->bareProduct().energy.size();
        auto gpu_ee_size = wgpuEE->bareProduct().energy.size();
	
	// Filling up the histograms on events sizes for EB and EE on both GPU and CPU
	std::cerr<<"EBGPU size ="<< gpu_eb_size<<std::endl;
	hSOIRechitsEBGPU->Fill(gpu_eb_size);
        hSOIRechitsEBCPU->Fill(cpu_eb_size);
        hSOIRechitsEEGPU->Fill(gpu_ee_size);
        hSOIRechitsEECPU->Fill(cpu_ee_size);
    /*    
     // condition that sizes on GPU and CPU should be the same for EB or EE
	if (cpu_eb_size != gpu_eb_size or cpu_ee_size != gpu_ee_size) {
          std::cerr << ie << ordinal[ie % 10] << " entry:\n"
                    << "  EB size: " << std::setw(4) << cpu_eb_size << " (cpu) vs " << std::setw(4) << gpu_eb_size << " (gpu)\n"
                    << "  EE size: " << std::setw(4) << cpu_ee_size << " (cpu) vs " << std::setw(4) << gpu_ee_size << " (gpu)" << std::endl;
		   
          continue;
        }
        assert(wgpuEB->bareProduct().energy.size() == wcpuEB->bareProduct().size());
        assert(wgpuEE->bareProduct().energy.size() == wcpuEE->bareProduct().size()); */
        auto const neb = wcpuEB->bareProduct().size(); //like cpu_eb_size but set to constant
        auto const nee = wcpuEE->bareProduct().size(); //like cpu_ee_size but set to constant
	
	
	// EB:
	// Find the energies on GPU and CPU reflecting the same did ID
        for (uint32_t i=0; i<neb; ++i) {
            auto const did_gpu = wgpuEB->bareProduct().did[i]; // set the did for the current RecHit
        //    test_r<<"\n EB did \t"<< did_gpu<<"\n";
	//    std::cerr<<"The iteration is "<<i<<" and the did_gpu "<< did_gpu<<std::endl; 
            auto const soi_enr_gpu = wgpuEB->bareProduct().energy[i]; //set the energy corresponding to the did value.
            auto const cpu_iter = wcpuEB->bareProduct().find(DetId{did_gpu}); // find the Rechit on CPU reflecting the same did
	    // In case of error failing to find the corresponding CPU input, print message and move on.
            if (cpu_iter == wcpuEB->bareProduct().end()) {
             //   std::cerr << ie << ordinal[ie % 10] << " entry\n"
             //             << "  Did not find a DetId " << did_gpu
             //             << " in a CPU collection\n";
                continue;
            }
            // Set the variables for cpu's energy on CPU and Chi2 for both CPU and GPU
            auto const soi_enr_cpu = cpu_iter->energy();
            auto const chi2_gpu = wgpuEB->bareProduct().chi2[i];
            auto const chi2_cpu = cpu_iter->chi2();

	    // Fill the energy and Chi2 histograms for GPU and CPU and their comparisons with delta
            hSOIEnergiesEBGPU->Fill(soi_enr_gpu);
            hSOIEnergiesEBCPU->Fill(soi_enr_cpu);
            hSOIEnergiesEBGPUvsCPU->Fill(soi_enr_cpu, soi_enr_gpu);
            hSOIEnergiesEBdeltavsCPU->Fill(soi_enr_cpu, soi_enr_gpu-soi_enr_cpu);
            hChi2EBGPU->Fill(chi2_gpu);
            hChi2EBCPU->Fill(chi2_cpu);
            hChi2EBGPUvsCPU->Fill(chi2_cpu, chi2_gpu);
            hChi2EBdeltavsCPU->Fill(chi2_cpu, chi2_gpu-chi2_cpu);

	    // Check if abs difference between GPU and CPU values for energy and Chi2 are smaller than eps, if not print message
      /*      if ((std::abs(soi_enr_gpu - soi_enr_cpu) >= eps_diff) or
                (std::abs(chi2_gpu - chi2_cpu) >= eps_diff) or std::isnan(chi2_gpu))
            {
                printf("EB eventid = %d chid = %d energy_gpu = %f energy_cpu %f chi2_gpu = %f chi2_cpu = %f\n",
                    ie, i, soi_enr_gpu, soi_enr_cpu, chi2_gpu, chi2_cpu);
                if (std::isnan(chi2_gpu))
                  printf("*** nan ***\n");
            } */
        }

        // EE:
	// Find the energies on GPU and CPU reflecting the same did ID
	std::cerr<<"nee is = "<< nee;
        for (uint32_t i=0; i<nee; ++i) {
	  test_r<<"EE i = \t"<< i<<"\n";
            auto const did_gpu = wgpuEE->bareProduct().did[i]; // set the did for the current RecHit
            test_r<<"EE did \t"<< did_gpu<<"\n";
            auto const soi_enr_gpu = wgpuEE->bareProduct().energy[i]; //set the energy corresponding to the did value.
            test_r<<"EE energy \t"<< soi_enr_gpu<<"\n";
            auto const cpu_iter = wcpuEE->bareProduct().find(DetId{did_gpu});  // find the Rechit on CPU reflecting the same did
	    // In case of error failing to find the corresponding CPU input, print message and move on.
            if (cpu_iter == wcpuEE->bareProduct().end()) {
              //  std::cerr << ie << ordinal[ie % 10] << " entry\n"
              //            << "  did not find a DetId " << did_gpu
              //            << " in a CPU collection\n";
                continue;
            }
            // Set the variables for cpu's energy on CPU and Chi2 for both CPU and GPU
            auto const soi_enr_cpu = cpu_iter->energy();
            auto const chi2_gpu = wgpuEE->bareProduct().chi2[i];
            auto const chi2_cpu = cpu_iter->chi2();

	    // Fill the energy and Chi2 histograms for GPU and CPU and their comparisons with delta
            hSOIEnergiesEEGPU->Fill(soi_enr_gpu);
            hSOIEnergiesEECPU->Fill(soi_enr_cpu);
            hSOIEnergiesEEGPUvsCPU->Fill(soi_enr_cpu, soi_enr_gpu);
            hSOIEnergiesEEdeltavsCPU->Fill(soi_enr_cpu, soi_enr_gpu-soi_enr_cpu);
            hChi2EEGPU->Fill(chi2_gpu);
            hChi2EECPU->Fill(chi2_cpu);
            hChi2EEGPUvsCPU->Fill(chi2_cpu, chi2_gpu);
            hChi2EEdeltavsCPU->Fill(chi2_cpu, chi2_gpu-chi2_cpu);

	     // Check if abs difference between GPU and CPU values for energy and Chi2 are smaller than eps, if not print message
            if ((std::abs(soi_enr_gpu - soi_enr_cpu) >= eps_diff) or
                (std::abs(chi2_gpu - chi2_cpu) >= eps_diff) or std::isnan(chi2_gpu))
            {
                printf("EE eventid = %d chid = %d enr_gpu = %f enr_cpu %f chi2_gpu = %f chi2_cpu = %f\n",
                    ie, static_cast<int>(neb+i), soi_enr_gpu, soi_enr_cpu, chi2_gpu, chi2_cpu);
                if (std::isnan(chi2_gpu))
                  printf("*** nan ***\n");
            }
        }
    }

    // Plotting the results:
    {
      // Canvas for sizes of EB and EE
      TCanvas cRechits("plots", "plots", 8200, 7200);
      cRechits.Divide(1, 2);
      
      // Canvas for values of energy and Chi2
      TCanvas c1("plots", "plots", 4200, 6200);
      c1.Divide(2, 3);
      TCanvas c2("plots", "plots", 4200, 6200);
      c2.Divide(2, 3);
      

      
      // Plotting the sizes of GPU vs CPU for each event of EB 
           cRechits.cd(1);
      {
          gPad->SetLogy();
	  gPad->SetTitle("Rechits for EB");
          hSOIRechitsEBCPU->SetLineColor(kRed);
          hSOIRechitsEBCPU->SetLineWidth(4.);
          hSOIRechitsEBCPU->Draw("");
          hSOIRechitsEBGPU->SetLineColor(kBlue);
          hSOIRechitsEBGPU->SetLineWidth(4.);
          hSOIRechitsEBGPU->Draw("sames");
          cRechits.Update();
          auto stats = (TPaveStats*)hSOIRechitsEBGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
	  
	  
      }
      // Plotting the sizes of GPU vs CPU for each event of EE 
	   cRechits.cd(2);
      {
          gPad->SetLogy();
	  gPad->SetTitle("Rechits for EE");
          hSOIRechitsEECPU->SetLineColor(kRed);
          hSOIRechitsEECPU->SetLineWidth(4.);
          hSOIRechitsEECPU->Draw("");
          hSOIRechitsEEGPU->SetLineColor(kBlue);
          hSOIRechitsEEGPU->SetLineWidth(4.);
          hSOIRechitsEEGPU->Draw("sames");
          cRechits.Update();
          auto stats = (TPaveStats*)hSOIRechitsEEGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
	  
	 
      }
      // Saving the canvas as .png for plots
      cRechits.SaveAs("ecal-rechit_cpu_gpu.png");

      // Pad for EB energy values on CPU and GPU
      c1.cd(1);
      {
          gPad->SetLogy();
          hSOIEnergiesEBCPU->SetLineColor(kBlack);
          hSOIEnergiesEBCPU->SetLineWidth(2.);
          hSOIEnergiesEBCPU->Draw("");
          hSOIEnergiesEBGPU->SetLineColor(kBlue);
          hSOIEnergiesEBGPU->SetLineWidth(2.);
          hSOIEnergiesEBGPU->Draw("sames");
          gPad->Update();
          auto stats = (TPaveStats*)hSOIEnergiesEBGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
      }
      
      // Pad for EE energy values on CPU and GPU
      c1.cd(2);
      {
          gPad->SetLogy();
          hSOIEnergiesEECPU->SetLineColor(kBlack);
          hSOIEnergiesEECPU->SetLineWidth(2.);
          hSOIEnergiesEECPU->Draw("");
          hSOIEnergiesEEGPU->SetLineColor(kBlue);
          hSOIEnergiesEEGPU->SetLineWidth(2.);
          hSOIEnergiesEEGPU->Draw("sames");
          gPad->Update();
          auto stats = (TPaveStats*)hSOIEnergiesEEGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
      }
      
      // Comparing GPU and delta with CPU for EB and EE energies
      c1.cd(3);
      hSOIEnergiesEBGPUvsCPU->Draw("COLZ");
      c1.cd(4);
      hSOIEnergiesEEGPUvsCPU->Draw("COLZ");
      c1.cd(5);
      hSOIEnergiesEBdeltavsCPU->Draw("COLZ");
      c1.cd(6);
      hSOIEnergiesEEdeltavsCPU->Draw("COLZ");

      c1.SaveAs("ecal-energies.png");

      // Pad for EB Chi2 values on CPU and GPU
      c2.cd(1);
      {
          gPad->SetLogy();
          hChi2EBCPU->SetLineColor(kBlack);
          hChi2EBCPU->SetLineWidth(1.);
          hChi2EBCPU->Draw("");
          hChi2EBGPU->SetLineColor(kBlue);
          hChi2EBGPU->SetLineWidth(1.);
          hChi2EBGPU->Draw("sames");
          gPad->Update();
          auto stats = (TPaveStats*)hChi2EBGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
      }
      
      // Pad for EE Chi2 values on CPU and GPU
      c2.cd(2);
      {
          gPad->SetLogy();
          hChi2EECPU->SetLineColor(kBlack);
          hChi2EECPU->SetLineWidth(1.);
          hChi2EECPU->Draw("");
          hChi2EEGPU->SetLineColor(kBlue);
          hChi2EEGPU->SetLineWidth(1.);
          hChi2EEGPU->Draw("sames");
          gPad->Update();
          auto stats = (TPaveStats*)hChi2EEGPU->FindObject("stats");
          auto y2 = stats->GetY2NDC();
          auto y1 = stats->GetY1NDC();
          stats->SetY2NDC(y1);
          stats->SetY1NDC(y1 - (y2-y1));
      }
      
      // Comparing GPU and delta with CPU for EB and EE Chi
      c2.cd(3);
      hChi2EBGPUvsCPU->Draw("COLZ");
      c2.cd(4);
      hChi2EEGPUvsCPU->Draw("COLZ");
      c2.cd(5);
      hChi2EBdeltavsCPU->Draw("COLZ");
      c2.cd(6);
      hChi2EEdeltavsCPU->Draw("COLZ");

      c2.SaveAs("ecal-chi2.png");
    }

    // Close all open files
    rf.Close();
    rfout.Write();
    rfout.Close();

    return 0;
}
