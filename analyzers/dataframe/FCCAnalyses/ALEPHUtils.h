#ifndef ALEPHUTILS_H
#define ALEPHUTILS_H

/*
  ALEPHUtils — ALEPH analysis utilities for FCCAnalyses.

  Namespace hierarchy:
    FCCAnalyses::ALEPH::Preprocess  — event/run filtering, class bits, event type from MC
    FCCAnalyses::ALEPH::Particle    — charged selection, jet constituent PID masks,
                                      constituent-type helpers
    FCCAnalyses::ALEPH::Track       — track quality selection, impact-parameter cuts,
                                      primary vertex, D0 flip
    FCCAnalyses::ALEPH::dEdx        — Bethe-Bloch, PID hypothesis p-values, dEdx getters

  Shared type aliases (FCCAnalyses::ALEPH namespace):
    FCCAnalysesJetConstituents      — RVec<ReconstructedParticleData>  (one jet)
    FCCAnalysesJetConstituentsData  — RVec<float>                      (per-constituent scalar)

  Example RDataFrame usage:
    namespace AL = FCCAnalyses::ALEPH;

    df = df.Filter(AL::Preprocess::sel_class_filter(16), {"ClassBitset"})
           .Define("charged", AL::Particle::sel_charged(1), {"ReconstructedParticles"})
           .Define("isMu",    AL::Particle::get_isMu,       {"JetConstituents"});
*/

#include "edm4hep/ReconstructedParticleCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/MCParticleData.h"
#include "edm4hep/TrackData.h"
#include "edm4hep/TrackState.h"
#include "edm4hep/ReconstructedParticleData.h"
#include "edm4hep/ParticleIDData.h"
#include "edm4hep/RecDqdx.h"
#include "edm4hep/EDM4hepVersion.h"

#include "FCCAnalyses/ReconstructedParticle.h"

#include <ROOT/RVec.hxx>
#include <TLorentzVector.h>

#include <array>
#include <bitset>
#include <set>
#include <string>
#include <vector>
#include <cstdint>


// =============================================================================
// Shared type aliases — FCCAnalyses::ALEPH
// =============================================================================
namespace FCCAnalyses { namespace ALEPH {

  namespace rv = ROOT::VecOps;

  using FCCAnalysesJetConstituents     = rv::RVec<edm4hep::ReconstructedParticleData>;
  using FCCAnalysesJetConstituentsData = rv::RVec<float>;

}} // namespace FCCAnalyses::ALEPH


// =============================================================================
// Preprocess — event/run filtering and classification
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace Preprocess {

/// Returns |PDG| of the first light quark (1–5) in the MC list if class bit 15
/// is set in ClassBit[0]. Returns -1 otherwise.
float getJetPID(const rv::RVec<uint32_t>&                ClassBit,
                const rv::RVec<edm4hep::MCParticleData>& particles);

/// Functor. Returns true if bit (m_class - 1) is set in ClassBit[0].
struct sel_class_filter {
  const int m_class;
  explicit sel_class_filter(int arg_class);
  bool operator()(const rv::RVec<uint32_t>& bitset_coll) const;
};

/// Returns a vector of all 1-based class indices whose bit is set in bitset_coll[0].
std::vector<int> bitsetToIndices(const rv::RVec<uint32_t>& bitset_coll);

/// Functor. Returns true if the event's run number is in the allowed set.
struct sel_runs_filter {
  const std::set<int>& m_runs_set;
  explicit sel_runs_filter(const std::set<int>& arg_runs_set);
  bool operator()(const rv::RVec<edm4hep::EventHeader>& event_header) const;
};

/// Returns |PDG| of the first particle with 0 < |PDG| < 6 in the MC list.
/// Returns -1 if none found.
float get_EventType(const rv::RVec<edm4hep::MCParticleData>& in);

}}} // namespace FCCAnalyses::ALEPH::Preprocess


// =============================================================================
// Particle — reconstructed-particle selection and jet-constituent utilities
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace Particle {

using FCCAnalysesJetConstituents     = ALEPH::FCCAnalysesJetConstituents;
using FCCAnalysesJetConstituentsData = ALEPH::FCCAnalysesJetConstituentsData;

/// Functor. Returns a subset collection of particles with |charge| == m_charge.
struct sel_charged {
  const int m_charge;
  explicit sel_charged(int arg_charge);
  edm4hep::ReconstructedParticleCollection
  operator()(const edm4hep::ReconstructedParticleCollection& in_coll) const;
};

/// Electrons: |charge| > 0 and mass ≈ 0.511 MeV
rv::RVec<FCCAnalysesJetConstituentsData>
get_isEl(const rv::RVec<FCCAnalysesJetConstituents>& jcs);

/// Muons: |charge| > 0 and mass ≈ 105.658 MeV
rv::RVec<FCCAnalysesJetConstituentsData>
get_isMu(const rv::RVec<FCCAnalysesJetConstituents>& jcs);

/// Charged hadrons: |charge| > 0 and mass ≈ 139.57 MeV (pion mass)
rv::RVec<FCCAnalysesJetConstituentsData>
get_isChargedHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs);

/// Neutral hadrons: PDG == 130 (KL0)
rv::RVec<FCCAnalysesJetConstituentsData>
get_isNeutralHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs);

/// Photons: PDG == 22
rv::RVec<FCCAnalysesJetConstituentsData>
get_isGamma(const rv::RVec<FCCAnalysesJetConstituents>& jcs);

/// Generic mask: returns 1 for constituents whose type value equals `type`.
rv::RVec<FCCAnalysesJetConstituentsData>
get_isType(const rv::RVec<FCCAnalysesJetConstituentsData>& jcs, float type);

/// Functor. Builds per-jet, per-constituent type vectors from ParticleIDData
/// and a jet-constituent index table.
struct build_constituents_Types {
  rv::RVec<FCCAnalysesJetConstituentsData>
  operator()(const rv::RVec<edm4hep::ParticleIDData>& rpid,
             const std::vector<std::vector<int>>&      indices) const;
};

}}} // namespace FCCAnalyses::ALEPH::Particle


// =============================================================================
// Track — track quality selection and primary vertex utilities
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace Track {

/// Bundles surviving TrackData and TrackStates after any selection step.
/// tracks[i] and trackStates[i] always correspond to the same track.
struct SelectedTracks {
  rv::RVec<edm4hep::TrackData>  tracks;
  rv::RVec<edm4hep::TrackState> trackStates;
};

rv::RVec<float> get_track_chi2      (const rv::RVec<edm4hep::TrackData>& tracks_in);
rv::RVec<float> get_track_ndf       (const rv::RVec<edm4hep::TrackData>& tracks_in);
rv::RVec<float> get_track_chi2_o_ndf(const rv::RVec<edm4hep::TrackData>& tracks_in);

/// Quality cuts: ndf > 0, chi2/ndf <= 10, covariance elements [0][2][9] positive
/// and finite. Throws std::runtime_error if a track has != 1 TrackState.
SelectedTracks
select_tracks_baseline(const rv::RVec<edm4hep::TrackData>&  tracks_in,
                       const rv::RVec<edm4hep::TrackState>& trackstates_in);

/// Retains only tracks with |D0| <= d0_upper_bound and |Z0| <= z0_upper_bound.
SelectedTracks
select_tracks_impactparameters(const SelectedTracks& input,
                               float d0_upper_bound,
                               float z0_upper_bound);

/// Functor. Returns the generator-level primary vertex as a TLorentzVector
/// (x, y, z [mm], t [mm/c]) from MCParticle data.
struct get_EventPrimaryVertexP4 {
  int m_genstatus = 21;
  get_EventPrimaryVertexP4();
  TLorentzVector operator()(const rv::RVec<edm4hep::MCParticleData>& in) const;
};

/// Returns a copy of the track-state collection with D0 and omega sign-flipped,
/// along with correlated covariance elements [1,4,6,8,10,12].
/// Index convention: https://bib-pubdb1.desy.de/record/81214, sec 5
rv::RVec<edm4hep::TrackState>
flipSign(const rv::RVec<edm4hep::TrackState>& tracks);

}}} // namespace FCCAnalyses::ALEPH::Track


// =============================================================================
// dEdx — Bethe-Bloch parametrisation and PID hypothesis testing
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace dEdx {

using FCCAnalysesJetConstituentsData = ALEPH::FCCAnalysesJetConstituentsData;

/// Evaluates a*(b + c*ln(p) + d*p^e_pow) for momentum p [GeV].
/// par = {a, b, c, d, e_pow}. Returns 0 for p <= 0.
double bethe_bloch(double p, const std::vector<double>& par);

/// Signed p-value for one dE/dx measurement against one particle hypothesis.
/// hypothesis: "e","mu","pi","K","p". is_wires: true for wire sensor, false for pads.
/// Returns -9 for unknown hypothesis, err <= 0, or non-finite residual.
double signed_p_value(double p, double dedx, double err,
                      const std::string& hypothesis, bool is_wires);

/// Returns {p_e, p_mu, p_pi, p_K, p_p} signed p-values in one call.
std::array<double, 5>
all_hypotheses_pvalues(double p, double dedx, double err, bool is_wires);

/// Functor. Matches jet constituents to their dEdx measurements via track links,
/// fills dEdx objects and 5-element PID p-value arrays per constituent.
/// Neutral particles and invalid measurements (dQdx.type != 0) get dummy values of -9.
struct build_constituents_dEdx_PIDhypo {

  struct dEdx_and_PID_result {
    rv::RVec<rv::RVec<edm4hep::RecDqdxData>>  dedx_constituents;
    rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_constituents;
  };

  dEdx_and_PID_result
  operator()(const rv::RVec<edm4hep::ReconstructedParticleData>& recoParticles,
             const rv::RVec<int>&                                  _recoParticlesIndices,
             const rv::RVec<edm4hep::RecDqdxData>&                dEdxCollection,
             const rv::RVec<int>&                                  _dEdxIndicesCollection,
             const std::vector<std::vector<int>>&                  jet_indices,
             bool                                                  is_wires) const;
};

rv::RVec<rv::RVec<float>>
get_dEdx_type(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>>& dedx_vec);

rv::RVec<rv::RVec<float>>
get_dEdx_value(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>>& dedx_vec);

rv::RVec<rv::RVec<float>>
get_dEdx_error(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>>& dedx_vec);

/// Extracts p-values for one hypothesis index: 0=e, 1=mu, 2=pi, 3=K, 4=p.
rv::RVec<rv::RVec<float>>
get_PID_pvalue(const rv::RVec<rv::RVec<std::array<double, 5>>>& pid_array_vec,
               int particle_index);

}}} // namespace FCCAnalyses::ALEPH::dEdx

#endif // ALEPHUTILS_H
