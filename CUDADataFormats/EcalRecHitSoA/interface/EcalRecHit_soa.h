#ifndef CUDADataFormats_EcalRecHitSoA_interface_EcalRecHit_soa_h
#define CUDADataFormats_EcalRecHitSoA_interface_EcalRecHit_soa_h

#include <vector>
#include <array>

// #include "DataFormats/EcalDigi/interface/EcalDataFrame.h"   -> why?

#include "CUDADataFormats/EcalRecHitSoA/interface/RecoTypes.h"
#include "HeterogeneousCore/CUDAUtilities/interface/CUDAHostAllocator.h"

// needed for "soa" definition
#include "CUDADataFormats/EcalRecHitSoA/interface/EcalUncalibratedRecHit_soa.h"

namespace ecal {
  
  template<typename L = Tag::soa>
  struct RecHit  : public Detail::Base<L> {
    
    RecHit() = default;
    RecHit(const RecHit&) = default;
    RecHit& operator=(const RecHit&) = default;
    
    RecHit(RecHit&&) = default;
    RecHit& operator=(RecHit&&) = default;
    
    typename type_wrapper<reco::StorageScalarType, L>::type energy;
    typename type_wrapper<reco::StorageScalarType, L>::type time;
    typename type_wrapper<reco::StorageScalarType, L>::type chi2;
    typename type_wrapper<uint32_t, L>::type flagBits; // store rechit condition (see Flags enum) in a bit-wise way
    typename type_wrapper<uint32_t, L>::type extra;    // packed uint32_t for timeError, chi2, energyError
    
    // refrence https://github.com/cms-sw/cmssw/blob/master/DataFormats/EcalRecHit/interface/EcalRecHit.h

    typename type_wrapper<uint32_t, L>::type did;
    
    
    template<typename U = L>
    typename std::enable_if<std::is_same<U, Tag::soa>::value, void>::type 
    resize(size_t size) {
      energy.resize(size);
      time.resize(size);
      chi2.resize(size);
      flagBits.resize(size);
      extra.resize(size);
      did.resize(size);
    }
  };
  
  using SoARecHitCollection = RecHit<Tag::soa>;
  
}

#endif 
// RecoLocalCalo_EcalRecAlgos_interface_EcalRecHit_soa_h

