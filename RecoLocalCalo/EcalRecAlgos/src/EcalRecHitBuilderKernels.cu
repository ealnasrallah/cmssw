#include "cuda.h"

#include "KernelHelpers.h"

#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalRecHit_soa.h"

//
//
#include "EcalRecHitBuilderKernels.h"


namespace ecal {
  namespace rechit {
    
    __global__
    void kernel_create_ecal_rehit(
                    uint32_t const* did_eb,
                    uint32_t const* did_ee,
                    ::ecal::reco::StorageScalarType const* amplitude_eb,   // in adc counts  
                    ::ecal::reco::StorageScalarType const* amplitude_ee,   // in adc counts  
                    ::ecal::reco::StorageScalarType* energy,   // in energy [GeV]  
                    int const nchannels
         ) {
      
      
//       
//    NB: energy   "type_wrapper<reco::StorageScalarType, L>::type" most likely std::vector<float>
//       
      
      int idx = threadIdx.x + blockDim.x*blockIdx.x;
      
      if (idx < nchannels) {
        
        // simple copy
        energy[idx] = amplitude_eb[idx];
        
      }
      
    }
    
    
    
    // host version, to be called by the plugin
    void create_ecal_rehit(
                  EventInputDataGPU const& eventInputGPU,
                  EventOutputDataGPU&      eventOutputGPU,
                  //     eventDataForScratchGPU_,
                  //     conditions,
                  //     configParameters_,
                  cuda::stream_t<>& cudaStream
             ){
    
      int nchannels = 10;
      
      unsigned int totalChannels = 10; //eventInputGPU.ebUncalibRecHits.nchannels +
//       eventInputGPU.eeUncalibRecHits.nchannels;
      
      unsigned int nchannels_per_block = 32;
      unsigned int threads_1d = 10 * nchannels_per_block;
      //   unsigned int blocks_1d = threads_1d > 10*totalChannels  ? 1 : (totalChannels*10 + threads_1d - 1) / threads_1d;
      unsigned int blocks_1d = 2;
      

//       kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
//         eventInputGPU.ebUncalibRecHits.did,
//         eventInputGPU.eeUncalibRecHits.did,
//         eventInputGPU.ebUncalibRecHits.amplitude, 
//         eventInputGPU.eeUncalibRecHits.amplitude, 
//         eventOutputGPU.energy,
//         nchannels
//       );
      
        
      
      
      

      // 
      // kernel
      //
        
//       kernel_create_ecal_rehit <<< blocks_1d, threads_1d >>> (
//                                amplitude,
//                                energy,
//                                nchannels
//       );
      
    }
    
    
//     
// //     error: cannot convert 'ecal::type_wrapper<float, ecal::Tag::soa>::type {aka std::vector<float, std::allocator<float> >}' to 'float*' for argument '2' to 'void ecal::rechit::create_ecal_rehit(const float*, float*, int)'
//     
// //     error: cannot convert 'ecal::type_wrapper<float, ecal::Tag::soa>::type {aka std::vector<float, CUDAHostAllocator<float, 0> >}' to 'float*' for argument '2' to 'void ecal::rechit::create_ecal_rehit(const float*, float*, int)'
//     
// //     error: cannot convert 'ecal::type_wrapper<float, ecal::Tag::soa>::type {aka std::vector<float, CUDAHostAllocator<float, 0> >}' to 'const float*' for argument '2' to 'void ecal::rechit::create_ecal_rehit(const float*, const float*, int)'
//     
// //     error: cannot convert 'ecal::type_wrapper<float, ecal::Tag::soa>::type {aka std::vector<float, CUDAHostAllocator<float, 0> >}' to 'ecal::type_wrapper<float, ecal::Tag::soa>::type* {aka std::vector<float, CUDAHostAllocator<float, 0> >*}'    
// 
//  error: no operator "=" matches these operands
//  operand types are: std::vector<ecal::reco::StorageScalarType, CUDAHostAllocator<ecal::reco::StorageScalarType, 0U>> = const float
// 
//    


    
  }
  
}

