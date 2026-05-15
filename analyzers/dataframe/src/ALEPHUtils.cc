#include "ALEPHUtils.h"

#include "edm4hep/Track.h"
#include "edm4hep/TrackData.h"

#include "FCCAnalyses/ReconstructedParticle2Track.h"
#include "FCCAnalyses/ReconstructedParticle2MC.h"
#include "FCCAnalyses/JetClusteringUtils.h"
#include "FCCAnalyses/MCParticle.h"

#include <cmath>
#include <stdexcept>
#include <unordered_map>


// =============================================================================
// Preprocess
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace Preprocess {

float getJetPID(const rv::RVec<uint32_t>&                ClassBit,
                const rv::RVec<edm4hep::MCParticleData>& particles)
{
    if (ClassBit.empty() || !std::bitset<32>(ClassBit[0])[15])
        return -1.0f;
    for (const auto& p : particles)
        if (std::abs(p.PDG) > 0 && std::abs(p.PDG) < 6)
            return static_cast<float>(std::abs(p.PDG));
    return -1.0f;
}

sel_class_filter::sel_class_filter(int arg_class) : m_class(arg_class) {}

bool sel_class_filter::operator()(const rv::RVec<uint32_t>& bitset_coll) const
{
    if (bitset_coll.empty()) return false;
    return std::bitset<32>(bitset_coll[0])[m_class - 1];
}

std::vector<int> bitsetToIndices(const rv::RVec<uint32_t>& bitset_coll)
{
    std::vector<int> indices;
    std::bitset<32> bits(bitset_coll[0]);
    for (size_t i = 0; i < bits.size(); ++i)
        if (bits.test(i))
            indices.push_back(static_cast<int>(i) + 1);
    return indices;
}

sel_runs_filter::sel_runs_filter(const std::set<int>& arg_runs_set)
    : m_runs_set(arg_runs_set) {}

bool sel_runs_filter::operator()(const rv::RVec<edm4hep::EventHeader>& event_header) const
{
    if (event_header.empty()) return false;
    return m_runs_set.count(event_header[0].getRunNumber()) > 0;
}

float get_EventType(const rv::RVec<edm4hep::MCParticleData>& in)
{
    for (const auto& p : in)
        if (std::abs(p.PDG) > 0 && std::abs(p.PDG) < 6)
            return static_cast<float>(std::abs(p.PDG));
    return -1.0f;
}

}}} // namespace FCCAnalyses::ALEPH::Preprocess


// =============================================================================
// Particle
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace Particle {

sel_charged::sel_charged(int arg_charge) : m_charge(arg_charge) {}

edm4hep::ReconstructedParticleCollection
sel_charged::operator()(const edm4hep::ReconstructedParticleCollection& in_coll) const
{
    edm4hep::ReconstructedParticleCollection result;
    result.setSubsetCollection();
    for (const auto& p : in_coll)
        if (std::abs(p.getCharge()) == m_charge)
            result.push_back(p);
    return result;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isEl(const rv::RVec<FCCAnalysesJetConstituents>& jcs)
{
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());
    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());
        for (const auto& c : jet)
            mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.000511) < 1e-5) ? 1.f : 0.f);
        out.push_back(std::move(mask));
    }
    return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isMu(const rv::RVec<FCCAnalysesJetConstituents>& jcs)
{
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());
    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());
        for (const auto& c : jet)
            mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.105658) < 1e-3) ? 1.f : 0.f);
        out.push_back(std::move(mask));
    }
    return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isChargedHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs)
{
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());
    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());
        for (const auto& c : jet)
            mask.push_back((std::abs(c.charge) > 0 && std::abs(c.mass - 0.13957) < 1e-3) ? 1.f : 0.f);
        out.push_back(std::move(mask));
    }
    return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isNeutralHad(const rv::RVec<FCCAnalysesJetConstituents>& jcs)
{
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());
    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());
        for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
            mask.push_back((c.PDG == 130) ? 1.f : 0.f);
#else
            mask.push_back((c.type == 130) ? 1.f : 0.f);
#endif
        out.push_back(std::move(mask));
    }
    return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isGamma(const rv::RVec<FCCAnalysesJetConstituents>& jcs)
{
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());
    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());
        for (const auto& c : jet)
#if edm4hep_VERSION > EDM4HEP_VERSION(0, 10, 5)
            mask.push_back((c.PDG == 22) ? 1.f : 0.f);
#else
            mask.push_back((c.type == 22) ? 1.f : 0.f);
#endif
        out.push_back(std::move(mask));
    }
    return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
get_isType(const rv::RVec<FCCAnalysesJetConstituentsData>& jcs, float type)
{
    rv::RVec<FCCAnalysesJetConstituentsData> out;
    out.reserve(jcs.size());
    for (const auto& jet : jcs) {
        FCCAnalysesJetConstituentsData mask;
        mask.reserve(jet.size());
        for (const auto& c : jet)
            mask.push_back(c == type ? 1.f : 0.f);
        out.push_back(std::move(mask));
    }
    return out;
}

rv::RVec<FCCAnalysesJetConstituentsData>
build_constituents_Types::operator()(const rv::RVec<edm4hep::ParticleIDData>& rpid,
                                     const std::vector<std::vector<int>>&      indices) const
{
    rv::RVec<FCCAnalysesJetConstituentsData> jcs;
    for (const auto& jet_index : indices) {
        FCCAnalysesJetConstituentsData jc;
        for (const auto& idx : jet_index)
            jc.push_back(rpid.at(idx).type);
        jcs.push_back(jc);
    }
    return jcs;
}

}}} // namespace FCCAnalyses::ALEPH::Particle


// =============================================================================
// Track
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace Track {

rv::RVec<float> get_track_chi2(const rv::RVec<edm4hep::TrackData>& tracks_in)
{
    rv::RVec<float> out;
    for (const auto& t : tracks_in) out.push_back(t.chi2);
    return out;
}

rv::RVec<float> get_track_ndf(const rv::RVec<edm4hep::TrackData>& tracks_in)
{
    rv::RVec<float> out;
    for (const auto& t : tracks_in) out.push_back(t.ndf);
    return out;
}

rv::RVec<float> get_track_chi2_o_ndf(const rv::RVec<edm4hep::TrackData>& tracks_in)
{
    rv::RVec<float> out;
    for (const auto& t : tracks_in) out.push_back(static_cast<float>(t.chi2 / t.ndf));
    return out;
}

SelectedTracks
select_tracks_baseline(const rv::RVec<edm4hep::TrackData>&  tracks_in,
                       const rv::RVec<edm4hep::TrackState>& trackstates_in)
{
    SelectedTracks selected;
    for (const auto& track : tracks_in) {
        if (track.ndf == 0)               continue;
        if (track.chi2 / track.ndf > 10.) continue;

        const auto n_states = track.trackStates_end - track.trackStates_begin;
        if (n_states != 1)
            throw std::runtime_error(
                "ALEPHUtils::Track::select_tracks_baseline: "
                "expected exactly one TrackState per Track");

        for (unsigned idx = track.trackStates_begin; idx < track.trackStates_end; ++idx) {
            const auto& state = trackstates_in[idx];
            const auto& cov   = state.covMatrix;
            if (cov[0] <= 1e-12 || cov[2] <= 1e-12 || cov[9] <= 1e-12)          continue;
            if (!std::isfinite(cov[0]) || !std::isfinite(cov[2]) || !std::isfinite(cov[9])) continue;
            selected.trackStates.push_back(state);
        }
        selected.tracks.push_back(track);
    }
    return selected;
}

SelectedTracks
select_tracks_impactparameters(const SelectedTracks& input,
                               float d0_upper_bound,
                               float z0_upper_bound)
{
    SelectedTracks selected;
    for (size_t i = 0; i < input.tracks.size(); ++i) {
        const auto& state = input.trackStates[i];
        if (std::abs(state.D0) > d0_upper_bound) continue;
        if (std::abs(state.Z0) > z0_upper_bound) continue;
        selected.tracks.push_back(input.tracks[i]);
        selected.trackStates.push_back(state);
    }
    return selected;
}

get_EventPrimaryVertexP4::get_EventPrimaryVertexP4() {}

TLorentzVector
get_EventPrimaryVertexP4::operator()(const rv::RVec<edm4hep::MCParticleData>& in) const
{
    TLorentzVector result(-1e12, -1e12, -1e12, -1e12);

    // Primary: take the first particle (Pythia8 hard-process incoming, status 21)
    for (const auto& p : in) {
        result = TLorentzVector(p.vertex.x, p.vertex.y, p.vertex.z,
                                p.time * 1.0e3 * 2.99792458e+8);
        return result;
    }

    // Fallback: first particle with generatorStatus == 2 and non-zero z vertex
    for (const auto& p : in) {
        if (p.generatorStatus == 2 && std::abs(p.vertex.z) > 1.e-12) {
            result = TLorentzVector(p.vertex.x, p.vertex.y, p.vertex.z,
                                    p.time * 1.0e3 * 2.99792458e+8);
            break;
        }
    }
    return result;
}

rv::RVec<edm4hep::TrackState>
flipSign(const rv::RVec<edm4hep::TrackState>& tracks)
{
    rv::RVec<edm4hep::TrackState> out;
    out.reserve(tracks.size());
    for (const auto& t : tracks) {
        edm4hep::TrackState tt = t;
        tt.D0    = -tt.D0;
        tt.omega = -tt.omega;
        tt.covMatrix[1]  = -tt.covMatrix[1];   // cov(D0, phi)
        tt.covMatrix[4]  = -tt.covMatrix[4];   // cov(phi, omega)
        tt.covMatrix[6]  = -tt.covMatrix[6];   // cov(D0, z0)
        tt.covMatrix[8]  = -tt.covMatrix[8];   // cov(omega, z0)
        tt.covMatrix[10] = -tt.covMatrix[10];  // cov(D0, tanLambda)
        tt.covMatrix[12] = -tt.covMatrix[12];  // cov(omega, tanLambda)
        out.push_back(std::move(tt));
    }
    return out;
}

}}} // namespace FCCAnalyses::ALEPH::Track


// =============================================================================
// dEdx
// =============================================================================
namespace FCCAnalyses { namespace ALEPH { namespace dEdx {

// Bethe-Bloch parameter table: hypothesis -> sensor -> {a, b, c, d, e_pow}
// Fits by Matteo Reale for pads and wires.
static const std::unordered_map<std::string,
             std::unordered_map<std::string, std::vector<double>>> k_bb_params = {
  {"e",  {{"pads",  {0.9932659976792287,  1.7094419426691188,  0.07384141776786941, 0.0,           -2.0        }},
          {"wires", {0.7930482536047483,  2.1221617291314856,  0.048686287341931526, 0.0,           -2.0        }}}},
  {"mu", {{"pads",  {0.8590747,           1.32919203,          0.18728265,           0.01248365,    -1.96898002 }},
          {"wires", {0.55374855,          1.97190215,          0.26416447,           0.01549576,    -1.972309   }}}},
  {"pi", {{"pads",  {0.536889582,         2.10964282,          0.269949484,          0.00310166058, -3.49991524 }},
          {"wires", {0.7922792,           1.30952175,          0.19162276,           0.01523584,    -2.4295921  }}}},
  {"K",  {{"pads",  {0.25619823,          3.84049172,          0.53610855,           0.65090594,    -2.44656219 }},
          {"wires", {0.45920637,          1.82346531,          0.34181439,           0.52306056,    -2.21132773 }}}},
  {"p",  {{"pads",  {0.73189858,          1.05891917,          0.256201,             1.34293618,    -2.02610177 }},
          {"wires", {0.68811606,          1.03354162,          0.24500224,           1.44362792,    -2.08656854 }}}}
};

static const edm4hep::RecDqdxData k_dEdx_dummy = []() {
    edm4hep::RecDqdxData d{};
    d.dQdx.value = -9.0f;
    d.dQdx.error = -9.0f;
    d.dQdx.type  = -9.0f;
    return d;
}();

static const std::array<double, 5> k_pid_dummy{{-9.0, -9.0, -9.0, -9.0, -9.0}};

double bethe_bloch(double p, const std::vector<double>& par)
{
    if (p <= 0.0) return 0.0;
    return par[0] * (par[1] + par[2] * std::log(p) + par[3] * std::pow(p, par[4]));
}

double signed_p_value(double p, double dedx, double err,
                      const std::string& hypothesis, bool is_wires)
{
    if (err <= 0.0) return -9.0;

    auto part_it = k_bb_params.find(hypothesis);
    if (part_it == k_bb_params.end()) return -9.0;

    const std::string sensor = is_wires ? "wires" : "pads";
    auto sens_it = part_it->second.find(sensor);
    if (sens_it == part_it->second.end()) return -9.0;

    const double residual = (dedx - bethe_bloch(p, sens_it->second)) / err;
    if (!std::isfinite(residual)) return -9.0;

    const double sf = 0.5 * std::erfc(std::fabs(residual) / std::sqrt(2.0));
    return 2.0 * sf * (residual > 0 ? 1.0 : -1.0);
}

std::array<double, 5>
all_hypotheses_pvalues(double p, double dedx, double err, bool is_wires)
{
    static const std::string hypos[5] = {"e", "mu", "pi", "K", "p"};
    std::array<double, 5> results;
    for (int i = 0; i < 5; ++i)
        results[i] = signed_p_value(p, dedx, err, hypos[i], is_wires);
    return results;
}

build_constituents_dEdx_PIDhypo::dEdx_and_PID_result
build_constituents_dEdx_PIDhypo::operator()(
    const rv::RVec<edm4hep::ReconstructedParticleData>& recoParticles,
    const rv::RVec<int>&                                 _recoParticlesIndices,
    const rv::RVec<edm4hep::RecDqdxData>&               dEdxCollection,
    const rv::RVec<int>&                                 _dEdxIndicesCollection,
    const std::vector<std::vector<int>>&                 jet_indices,
    bool                                                 is_wires) const
{
    rv::RVec<rv::RVec<edm4hep::RecDqdxData>>  dedx_constituents;
    rv::RVec<rv::RVec<std::array<double, 5>>> pid_array_constituents;

    // Build track_index -> dEdx lookup to avoid O(N*M) inner loops
    std::unordered_map<int, edm4hep::RecDqdxData> track_to_dEdx;
    track_to_dEdx.reserve(_dEdxIndicesCollection.size());
    for (size_t i = 0; i < _dEdxIndicesCollection.size(); ++i)
        track_to_dEdx[_dEdxIndicesCollection[i]] = dEdxCollection[i];

    for (const auto& jet_const_indices : jet_indices) {
        rv::RVec<edm4hep::RecDqdxData>     jet_dEdx;
        rv::RVec<std::array<double, 5>>    jet_pid;

        for (int constituent_index : jet_const_indices) {
            const auto& reco = recoParticles[constituent_index];
            TLorentzVector tlv;
            tlv.SetXYZM(reco.momentum.x, reco.momentum.y, reco.momentum.z, reco.mass);

            bool found = false;
            for (int t = reco.tracks_begin; t < reco.tracks_end; ++t) {
                auto it = track_to_dEdx.find(_recoParticlesIndices[t]);
                if (it == track_to_dEdx.end()) continue;

                const auto& d = it->second;
                if (d.dQdx.type == 0) {
                    jet_dEdx.push_back(d);
                    jet_pid.push_back(all_hypotheses_pvalues(tlv.P(), d.dQdx.value,
                                                             d.dQdx.error, is_wires));
                } else {
                    jet_dEdx.push_back(k_dEdx_dummy);
                    jet_pid.push_back(k_pid_dummy);
                }
                found = true;
                break;
            }
            if (!found) {
                jet_dEdx.push_back(k_dEdx_dummy);
                jet_pid.push_back(k_pid_dummy);
            }
        }
        dedx_constituents.push_back(std::move(jet_dEdx));
        pid_array_constituents.push_back(std::move(jet_pid));
    }
    return {dedx_constituents, pid_array_constituents};
}

rv::RVec<rv::RVec<float>>
get_dEdx_type(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>>& dedx_vec)
{
    rv::RVec<rv::RVec<float>> out;
    for (const auto& inner : dedx_vec) {
        rv::RVec<float> v;
        for (const auto& d : inner) v.push_back(d.dQdx.type);
        out.push_back(std::move(v));
    }
    return out;
}

rv::RVec<rv::RVec<float>>
get_dEdx_value(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>>& dedx_vec)
{
    rv::RVec<rv::RVec<float>> out;
    for (const auto& inner : dedx_vec) {
        rv::RVec<float> v;
        for (const auto& d : inner) v.push_back(d.dQdx.value);
        out.push_back(std::move(v));
    }
    return out;
}

rv::RVec<rv::RVec<float>>
get_dEdx_error(const rv::RVec<rv::RVec<edm4hep::RecDqdxData>>& dedx_vec)
{
    rv::RVec<rv::RVec<float>> out;
    for (const auto& inner : dedx_vec) {
        rv::RVec<float> v;
        for (const auto& d : inner) v.push_back(d.dQdx.error);
        out.push_back(std::move(v));
    }
    return out;
}

rv::RVec<rv::RVec<float>>
get_PID_pvalue(const rv::RVec<rv::RVec<std::array<double, 5>>>& pid_array_vec,
               int particle_index)
{
    rv::RVec<rv::RVec<float>> out;
    for (const auto& inner : pid_array_vec) {
        rv::RVec<float> v;
        for (const auto& d : inner) v.push_back(static_cast<float>(d[particle_index]));
        out.push_back(std::move(v));
    }
    return out;
}

}}} // namespace FCCAnalyses::ALEPH::dEdx
