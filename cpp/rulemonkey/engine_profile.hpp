#pragma once

// =============================================================================
// Dev profiler / invariant infrastructure for engine.cpp.
//
// This header carries the gate constants, sampling cadences, and file-scope
// profile structs/instances that the engine uses to instrument hot paths
// during optimization sprints.  None of it is part of the runtime contract
// or the public API — pure dev-time scaffolding.
//
// In default builds the master macro `RM_DEV_PROFILES` is undefined, every
// gate is `false`, and every `if constexpr (k*Profile)` block in engine.cpp
// is dead-stripped by the compiler.  Opt in with:
//
//   cmake --preset release -DRULEMONKEY_ENABLE_DEV_PROFILES=ON
//
// Then flip individual gates here back to `false` to disable just that
// profiler — common during a sprint when you want only one phase
// instrumented.  Invariant gates (the second block below) stay `false`
// even when profilers are on; they run an O(N) reference-vs-fast-path
// equality check on every call and are intended for correctness
// verification during refactors of the matching/sampling fast paths,
// not for production builds.
//
// The reporting blocks at the end of run_ssa (still in engine.cpp because
// they touch private Impl state) emit a `[<category> profile]` block to
// stderr when their gate is on.  See engine.cpp for the report formats.
// =============================================================================

#include <array>
#include <cstdint>

namespace rulemonkey {

// ---------------------------------------------------------------------------
// Master gate.  Driven by the RULEMONKEY_ENABLE_DEV_PROFILES cmake option.
// ---------------------------------------------------------------------------
#ifdef RM_DEV_PROFILES
inline constexpr bool kDevProfilesEnabled = true;
#else
inline constexpr bool kDevProfilesEnabled = false;
#endif

// ---------------------------------------------------------------------------
// Profile gates.  When the master gate is on, these are individually
// editable: change `kDevProfilesEnabled` to `false` here to mute one
// profiler without touching engine.cpp.
// ---------------------------------------------------------------------------

// AgentPool::split_complex_if_needed: per-call counters (singleton /
// cycle-bond / tree-bond) and complex-size histograms always-on, plus
// sampled per-phase wall-time every Kth split call.
inline constexpr bool kRemoveBondProfile = kDevProfilesEnabled;
inline constexpr int kRemoveBondProfileSampleEvery = 8;

// count_multi_mol_fast dispatcher and generic body: always-on call counters
// (fm_hits vs generic, zero-seed vs singleton-pattern shares, seed-emb /
// bfs-visited sums + maxes) plus sampled per-sub-phase wall-time every
// Kth generic call.
inline constexpr bool kCountMultiProfile = kDevProfilesEnabled;
inline constexpr int kCountMultiProfileSampleEvery = 8;

// count_2mol_1bond_fc (the engaged-set hot spot specialized matcher,
// 215-269 ns/call): always-on reject taxonomy + size-stat counters,
// plus per-call-rotated sampled chrono timing for 5 sub-phases (rotor %
// 5 keeps chrono overhead <=10% when SampleEvery is large enough — the
// 215 ns-per-call baseline puts a hard ceiling on sampling density).
inline constexpr bool kCmmFcProfile = kDevProfilesEnabled;
inline constexpr int kCmmFcProfileSampleEvery = 32;

// select_reactants: classifies each call into one of five inner paths
// (zero / uni_single / uni_multi_fm / uni_multi_gen / bimol), times the
// whole call K-sampled per path, and accumulates always-on counters for
// the inner sampler (Fenwick vs linear scan), per-path work widths, and
// null outcomes.  Default K=16.
inline constexpr bool kSelectReactantsProfile = kDevProfilesEnabled;
inline constexpr int kSelectReactantsProfileSampleEvery = 16;

// fire_rule: per-op-type call counts (always, cheap) plus sampled per-
// sub-phase wall-time every Kth fire_rule call.  K chosen to keep chrono
// overhead below ~5% of fire_rule wall.
inline constexpr bool kFireRuleProfile = kDevProfilesEnabled;
inline constexpr int kFireRuleProfileSampleEvery = 8;

// incremental_update: two-level sampling.  Every Kth call is fully sampled
// (outer phase brackets); within a sampled call, every Mth per-mid
// inner-loop entry is inner-phase sampled (8 sub-phases inside the
// recompute body).
inline constexpr bool kIncrUpdateProfile = kDevProfilesEnabled;
inline constexpr int kIncrUpdateProfileSampleEvery = 8;
inline constexpr int kIncrUpdateProfileInnerSample = 32;

// record_at: splits wall into compute_observables vs result.record_time_point
// halves; evaluate_observable accumulates per-observable wall, branch-taken
// counts (Species / Molecules), complex-walk counts, and embedding call
// counts.  Sampling period K=1 — record_at fires only n_steps+1 times so
// full-sample chrono overhead is < 0.01% of wall.
inline constexpr bool kRecordAtProfile = kDevProfilesEnabled;

// incremental_update_observables / flush_species_incr_observables:
// sampled wall-time per sub-phase and per-obs evaluate counts/ns.  At
// init, a one-shot inventory dump flags obs that share a pattern
// signature so dedup opportunities can be audited.
inline constexpr bool kObsIncrProfile = kDevProfilesEnabled;
inline constexpr int kObsIncrProfileSampleEvery = 16;

// ---------------------------------------------------------------------------
// Invariant gates.  Slow correctness checks that run a reference path
// alongside a fast path and abort on divergence.  OFF by default even when
// the master profiler gate is on — the workflow is "flip on while
// refactoring the fast path, run the suite, flip off."
// ---------------------------------------------------------------------------

// sample_molecule_weighted: verifies Fenwick.find(r) returns a mid whose
// mol-id-ordered cumulative weight range contains r.  O(N) per sampler
// call; ship off.
inline constexpr bool kFenwickInvariant = false;

// FastMatchSlot specialization: runs BOTH the specialized and generic
// embedding paths and aborts on mismatch.  With the gate on, count_multi
// does strictly MORE work, so wall time will NOT improve — A/B benchmarks
// must be taken with the gate off.
inline constexpr bool kFastMatchInvariant = false;

// Cycle-bond-counter product-molecularity short-circuit: whenever
// cycle_bond_count==0 says "skip the BFS," ALSO run the BFS and abort
// on inconsistency (the BFS should always come back found_b==false on
// tree complexes).
inline constexpr bool kProductMolInvariant = false;

// FastMatchSlot-powered select_reactants fast path: every eligible
// unimolecular multi-mol select runs BOTH select_2mol_1bond_fc_match
// and the generic select_multi_mol_unimolecular, compares mol_ids /
// comp_ids and the post-call RNG state, aborts on divergence.
inline constexpr bool kFastSelectInvariant = false;

// Species-observable incremental tracker: runs both the incremental and
// full-walk paths at each sample point and asserts equality.  See
// kSpeciesIncrObs in engine.cpp for the feature flag itself.
inline constexpr bool kSpeciesIncrObsInvariant = false;

// 2-mol-1-bond-fc specialization for the obs-tracking path: compares both
// paths at every call and aborts on mismatch.  See kObsFastMatch in
// engine.cpp for the feature flag itself.
inline constexpr bool kObsFastMatchInvariant = false;

// ---------------------------------------------------------------------------
// File-scope profile structs.
//
// The two file-scope instances (cm_profile_, cmm_fc_profile_) are valid
// for the single-threaded driver in use today — switching to multi-engine
// concurrency would require per-Engine instances.
//
// SrProfile is defined here for symmetry but its only instance lives as
// a member of Engine::Impl in engine.cpp (the select_reactants profiler
// runs entirely inside the engine and accumulates per-engine state).
// ---------------------------------------------------------------------------

struct CountMultiProfile {
  // -- always on when gate is on --
  uint64_t calls = 0;                   // dispatcher entries
  uint64_t fm_hits = 0;                 // FastMatchSlot specialization taken
  uint64_t generic_calls = 0;           // fell through to generic body
  uint64_t sampled_calls = 0;           // generic calls bracketed with chrono
  uint64_t singleton_pattern_calls = 0; // pat_end - pat_start <= 1
  uint64_t zero_seed_calls = 0;         // seed_embs.empty() early return
  uint64_t disjoint_calls = 0;          // entered the unassigned (disjoint) path
  uint64_t seed_emb_sum = 0;            // sum seed_embs.size() over generic calls
  uint64_t seed_emb_max = 0;
  std::array<uint64_t, 7> seed_emb_hist = {0, 0, 0, 0, 0, 0, 0};
  uint64_t bfs_visited_sum = 0; // sum nodes popped from bfs_queue (all seed iters)
  uint64_t bfs_visited_max = 0;
  std::array<uint64_t, 7> bfs_visited_hist = {0, 0, 0, 0, 0, 0, 0};
  uint64_t n_pat_mols_sum = 0; // sum pat_end-pat_start (multi-mol patterns)
  std::array<uint64_t, 7> n_pat_mols_hist = {0, 0, 0, 0, 0, 0, 0};
  // -- sampled (xK at report) --
  uint64_t total_ns = 0;    // whole generic-body span
  uint64_t seed_emb_ns = 0; // initial count_embeddings_single(seed)
  uint64_t bfs_ns = 0;      // sum bfs while-loop walls across seed_embs
  uint64_t disjoint_ns = 0; // unassigned-recursion path
};

struct CmmFcProfile {
  // -- always on (when gate is on) --
  uint64_t fc_calls = 0;                 // dispatcher entries into count_2mol_1bond_fc
  uint64_t fc_inactive_seed = 0;         // !seed.active -> return 0
  uint64_t fc_seed_type_mismatches = 0;  // seed.type_index != fm.seed_type -> return 0
  uint64_t fc_seed_non_bond_rejects = 0; // any reject inside fm.seed_non_bond_checks loop
  uint64_t fc_candidate_iters = 0;       // sum seed_bond_candidates iterations entered
  // Reject taxonomy (sum + fc_total_matches == fc_candidate_iters)
  uint64_t fc_candidate_oob = 0;            // seed_ci >= n_seed
  uint64_t fc_candidate_state_rejects = 0;  // seed_bond_state_req mismatch
  uint64_t fc_candidate_bond_rejects = 0;   // seed_bond_bond_req mismatch OR partner_cid<0
  uint64_t fc_candidate_self_bond = 0;      // partner_mol_id == seed_mol_id
  uint64_t fc_partner_inactive = 0;         // !partner.active
  uint64_t fc_partner_type_mismatches = 0;  // partner.type_index != fm.partner_type
  uint64_t fc_plocal_not_found = 0;         // plocal < 0 after partner comp_ids scan
  uint64_t fc_plocal_not_ok = 0;            // plocal not in partner_bond_candidates
  uint64_t fc_partner_state_rejects = 0;    // partner_bond_state_req mismatch
  uint64_t fc_partner_non_bond_rejects = 0; // partner_non_bond_checks loop reject
  uint64_t fc_total_matches = 0;            // ++total fired

  // Size sums (/ fc_calls -> mean per call)
  uint64_t fc_seed_bond_candidates_sum = 0;
  uint64_t fc_partner_bond_candidates_sum = 0;
  uint64_t fc_partner_non_bond_checks_sum = 0;
  uint64_t fc_seed_non_bond_checks_sum = 0;
  uint64_t fc_seed_bond_candidates_max = 0;
  uint64_t fc_partner_bond_candidates_max = 0;
  uint64_t fc_partner_non_bond_checks_max = 0;
  std::array<uint64_t, 7> fc_partner_non_bond_checks_hist = {0, 0, 0, 0, 0, 0, 0}; // 0,1,2,3,4,5,6+

  // -- sampled per-phase (xKxN_rotor at report) --
  // Phases: 0=seed_checks, 1=partner_trace, 2=partner_local_scan,
  //         3=plocal_ok_scan, 4=partner_non_bond_checks.
  uint64_t sampled_calls = 0;
  std::array<uint64_t, 5> phase_ns = {0, 0, 0, 0, 0};
  std::array<uint64_t, 5> phase_hits = {0, 0, 0, 0, 0};
};

struct SrProfile {
  // Path identity (classification is exclusive; sum == total calls).
  static constexpr int kPathZero = 0;
  static constexpr int kPathUniSingle = 1;
  static constexpr int kPathUniMultiFm = 2;
  static constexpr int kPathUniMultiGen = 3;
  static constexpr int kPathBimol = 4;
  static constexpr int kNPaths = 5;

  // -- always on (when gate is on) --
  uint64_t calls = 0;
  std::array<uint64_t, kNPaths> path_calls = {0, 0, 0, 0, 0};
  // Outcomes — summed over all paths.
  std::array<uint64_t, kNPaths> path_null_no_seed = {0, 0, 0, 0, 0}; // sampler returned -1
  std::array<uint64_t, kNPaths> path_null_post = {0, 0, 0, 0, 0};    // sampled but match empty
  std::array<uint64_t, kNPaths> path_success = {0, 0, 0, 0, 0};

  // Sampler (sample_molecule_weighted / sample_molecule_by_local_propensity)
  // sub-decisions — global, not per-path.  Counted inside the sampler body.
  uint64_t sampler_calls = 0;          // total sample_molecule_weighted entries
  uint64_t sampler_fenwick_uses = 0;   // Fenwick find returned valid mol
  uint64_t sampler_fenwick_drifts = 0; // Fenwick returned invalid, linear fallback
  // Tri-state drift classifier (sum == sampler_fenwick_drifts).
  uint64_t sampler_drift_invalid_mid = 0;   // mid<0 or mid>=pool.molecule_count()
  uint64_t sampler_drift_inactive_mol = 0;  // mid valid but pool.molecule(mid).active==false
  uint64_t sampler_drift_type_mismatch = 0; // mid active but type_index != requested
  // Sub-classifier on invalid_mid: target vs tree.sum().
  uint64_t sampler_drift_target_eq_sum = 0; // r >= tree_sum (weight-loss in tree)
  uint64_t sampler_drift_target_lt_sum = 0; // r < tree_sum but find() still returned n
  // Magnitude of weight loss when target>=sum (sum of (r - tree_sum) divided by total).
  double sampler_drift_excess_sum = 0.0;
  double sampler_drift_total_sum = 0.0;  // sum of `total` over invalid_mid events
  uint64_t sampler_linear_calls = 0;     // linear scan (no Fenwick path)
  uint64_t sampler_empty_pool = 0;       // mols empty or total <= 0
  uint64_t sampler_local_prop_calls = 0; // sample_molecule_by_local_propensity

  // Per-path work counters (describe the loop widths inside each path).
  uint64_t uni_single_embs_sum = 0;        // embs.size() inside uni_single body
  uint64_t uni_single_embs_empty = 0;      // embs empty (bailed as null)
  uint64_t uni_mm_fm_success = 0;          // fm path returned non-empty
  uint64_t uni_mm_fm_null = 0;             // fm path returned empty
  uint64_t uni_mm_gen_seed_embs_sum = 0;   // seed_embs.size() sum
  uint64_t uni_mm_gen_seed_embs_empty = 0; // generic seed_embs empty
  uint64_t uni_mm_gen_success = 0;
  uint64_t uni_mm_gen_null = 0;
  uint64_t bimol_embs_a_sum = 0;       // bimol embs_a.size() sum
  uint64_t bimol_embs_b_sum = 0;       // bimol embs_b.size() sum
  uint64_t bimol_same_mol_rejects = 0; // mol_a == mol_b
  uint64_t bimol_same_cx_rejects = 0;  // block_same_complex_binding
  uint64_t bimol_embs_empty = 0;       // either embs_a or embs_b empty
  uint64_t bimol_resolve_calls = 0;    // resolve_pattern_context invocations
  uint64_t bimol_resolve_failures = 0; // resolve_pattern_context returned false

  // K-sampled whole-call chrono, bucketed per path.
  uint64_t sampled_calls = 0;
  std::array<uint64_t, kNPaths> path_ns = {0, 0, 0, 0, 0};
  std::array<uint64_t, kNPaths> path_hits = {0, 0, 0, 0, 0};
};

// Single-threaded-driver-only file-scope instances.  `inline` so the
// header can be included by multiple TUs without ODR violation; today
// only engine.cpp consumes them.
inline CountMultiProfile cm_profile_;
inline CmmFcProfile cmm_fc_profile_;

} // namespace rulemonkey
