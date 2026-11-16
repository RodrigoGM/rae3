// rae -- standalone, non-interactive RAE for .seg files.
// One required argument (the .seg); every parameter has a preset; the output set is fixed.
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include "rae/report.hpp"

static void usage() {
  std::cout <<
    "rae " << rae::RAE_VERSION << " -- recurrent copy-number alterations (reimplementation of RAE, Taylor et al. 2008)\n\n"
    "usage: rae <cohort.seg> [options]\n\n"
    "input\n"
    "  <cohort.seg>            ID chrom loc.start loc.end num.mark seg.mean [seg.sd]   (log2 ratios)\n"
    "  --bins FILE             marker/bin table: chrom start end. Default: grid reconstructed from the segments\n"
    "  --sample-stats FILE     per-sample bin-level statistics: ID dn q01 q025 (dn = MAPD), as written by\n"
    "                          rae3::as_seg(x, bin_noise = TRUE). Default: estimated from the segmentation\n"
    "  --purity-ploidy FILE    per-sample purity and ploidy: ID purity ploidy (e.g. from ACE, Rascal, ichorCNA).\n"
    "                          Used by --preset sequencing to place the detectors; otherwise estimated (G-13)\n"
    "  --normals FILE          .seg of non-tumour samples for the germline flag\n"
    "  --normal-stats FILE     sample statistics for the normals\n"
    "  --cnv default|none|FILE population CNV catalogue for the overlap columns of roi.tsv. default: bundled\n"
    "                          (hg19 DGV Gold Standard; hg38 gnomAD-SV v4.1; chm13 HPRC deletions + lifted gnomAD gains).\n"
    "                          FILE: chrom start end type(gain|loss) [frequency]. Overlap is reported, never filtered\n"
    "  --build NAME            hg38 (default), hg19, chm13, auto (best-fitting bundled build), custom\n"
    "output\n"
    "  --out PREFIX            default: the input path without its extension, plus .rae\n"
    "                          per tumour, PREFIX.samples.tsv has rae_fga (fraction of the genome altered by\n"
    "                          the RAE calls; the Hieronymus 2014 CNA burden), fga (by a log2 threshold, as in\n"
    "                          cBioPortal) and n_altered; see the fga options below\n"
    "method\n"
    "  --preset NAME           sequencing (default): detectors placed on the copy-number ladder, one-copy\n"
    "                          step from the data or from --purity-ploidy (G-13; docs/spec/08_britroc.md).\n"
    "                          For array CGH, detectors from the noise model: cgh, the released code\n"
    "                          (RAE beta 0.91; docs/spec/03_original_code_diff.md); taylor2008, the method\n"
    "                          as the text of Taylor et al. 2008 describes it. The preset is applied first,\n"
    "                          wherever it appears; the other options then override it\n"
    "  --set name=value        set any parameter by name (all names are listed in run.json)\n"
    "sensitivity (defaults of the sequencing and taylor2008 presets)\n"
    "  --fdr X                 significance cutoff on q            [0.01]\n"
    "  --peak-q X              eligibility for peak detection      [0.001]\n"
    "  --shoulder X            peak sensitivity, x median error    [2]\n"
    "  --n-perm N              permutations                        [10000]\n"
    "  --seed N                random seed                         [1]\n"
    "inferred choices (see docs/spec/02_gap_register.md)\n"
    "  --null-positions N      positions scored per permutation, 0 = all   [100]   G-08f\n"
    "  --recleave-every N      permutations per hotspot randomisation      [1000]  G-08g\n"
    "  --ubp-outlier-k X       segment-count outlier fence, 0 = off        [1.5]   G-07a\n"
    "  --beta-fallback-loss X  fallback slope on the loss side             [20]    G-03c\n"
    "  --beta-fallback-gain X  fallback slope on the gain side             [20]\n"
    "  --marker-sd X           per-marker SD when the seg has none         [estimated] G-01a\n"
    "  --sd-c X                constant of the SD estimate                 [11.5]  G-01a\n"
    "  --grid-spacing X        marker spacing in bp without --bins         [estimated] G-01d\n"
    "  --error-norm rms|sem    normalisation of the analytical error       [rms]   G-10b\n"
    "  --split-at-centromere   force a region boundary at centromeres      [off]   G-07e\n"
    "  --density-noise sem|sd  noise added back to segment means           [sem]   G-02g\n"
    "  --peak-window-lo X      gain/loss peak search window, x DN          [0.5]   G-03a\n"
    "  --peak-window-hi X                                                  [2]     G-03a\n"
    "later-paper thresholds (not RAE 2008)\n"
    "  --call-threshold X [0.9]   --amp-threshold X [0.5; 0.25 for array presets]   --burden-threshold X [0.75]\n"
    "fraction of the genome altered (.samples.tsv: rae_fga by RAE scores, fga by log2 ratio, n_altered)\n"
    "  --fga-threshold X       fga counts |log2 ratio| > X             [0.2]   G-12e\n"
    "  --fga-threshold-loss X  losses use log2 ratio < -|X| instead    [= --fga-threshold]  G-12e\n"
    "  --fga-assembly          divide rae_fga and fga by the whole autosomal assembly\n"
    "                          instead of each tumour's segment-covered span [off]   G-12a\n"
    "opt-in analyses (off by default; default output is unchanged)\n"
    "  --ploidy-prior          choose the ladder's base copy number by a soft prior on ploidy instead of the\n"
    "                          lowest-level rule (estimated ladders only); adds ploidy_hard, ploidy_alt and\n"
    "                          ploidy_prior_ratio to .samples.tsv, flag ploidy_ambiguous  [off]   G-14a, G-14b\n"
    "  --ploidy-prior-wgd X    weight of the doubled-genome component      [0.3]   G-14a\n"
    "  --ploidy-prior-mu1 X    mean ploidy of the diploid-like component   [1.99]  G-14a\n"
    "  --ploidy-prior-mu2 X    mean ploidy of the doubled component        [3.31]  G-14a\n"
    "  --ploidy-prior-sd X     SD of both components                       [0.4]   G-14a\n"
    "  --ploidy-prior-ambig X  flag ploidy_ambiguous when best/runner-up < X [3]   G-14b\n"
    "  --arm-baseline M        how an arm's baseline level is estimated, median or twolevel [median]  G-14c\n"
    "                          (opt-in: takes effect only with an analysis that uses arm baselines)\n"
    "  --arm-min-markers N     an arm needs at least N markers to be eligible           [50]    G-14c\n"
    "  --focal-only            score each segment on its deviation from its arm's baseline, so significance\n"
    "                          and the null see only focal events; adds focal_excluded_frac to .samples.tsv [off]  G-14d\n"
    "  --focal-max-mb X        longest focal deviation in Mb; longer ones score 0     [20]    G-14d\n"
    "  --arm-test              arm-level cohort test: each tumour's detectors at its arm baselines, tested\n"
    "                          against a within-tumour permutation of the arms; writes .arms.tsv  [off]  G-14e\n"
    "  --arm-perm N            permutations of the arm-level null                     [10000] G-14e\n"
    "  --peak-independence     test each later peak of an event for independence from the higher ones: its score\n"
    "                          is recomputed with the tumours' runs through their apexes masked; adds resid_score,\n"
    "                          resid_p and independent to .roi.tsv  [off]  G-14f\n"
    "  --peak-run-min X        a run is a stretch of segments with detector score >= X  [0.5]   G-14f\n";
}

int main(int argc, char** argv) {
  if (argc < 2 || !std::strcmp(argv[1], "-h") || !std::strcmp(argv[1], "--help")) { usage(); return argc < 2 ? 1 : 0; }
  if (!std::strcmp(argv[1], "--version")) { std::cout << rae::RAE_VERSION << "\n"; return 0; }
  rae::Params P;
  std::string seg, bins, stats, normals, nstats, pp, build = "hg38", out, cnv = "default";
  try {
    for (int i = 1; i + 1 < argc; ++i)   // the preset resets every parameter, so it goes first
      if (!std::strcmp(argv[i], "--preset") && !rae::set_preset(P, argv[i + 1]))
        throw std::runtime_error("--preset must be sequencing, cgh or taylor2008");
    for (int i = 1; i < argc; ++i) {
      std::string a = argv[i];
      auto val = [&]() -> std::string {
        if (i + 1 >= argc) throw std::runtime_error("missing value for " + a);
        return argv[++i];
      };
      auto num = [&]() { std::string v = val(); char* e; double d = std::strtod(v.c_str(), &e); if (*e) throw std::runtime_error("not a number for " + a + ": " + v); return d; };
      if (a == "--preset") val();          // applied above
      else if (a == "--set") {
        std::string kv = val(); size_t eq = kv.find('=');
        char* e = nullptr;
        double d = eq == std::string::npos ? 0 : std::strtod(kv.c_str() + eq + 1, &e);
        if (eq == std::string::npos || (e && *e) || !rae::set_param(P, kv.substr(0, eq), d))
          throw std::runtime_error("--set expects name=number with a known parameter name: " + kv);
      }
      else if (a == "--bins") bins = val();
      else if (a == "--cnv") cnv = val();
      else if (a == "--sample-stats") stats = val();
      else if (a == "--purity-ploidy") pp = val();
      else if (a == "--normals") normals = val();
      else if (a == "--normal-stats") nstats = val();
      else if (a == "--build") build = val();
      else if (a == "--out") out = val();
      else if (a == "--fdr") P.fdr = num();
      else if (a == "--peak-q") P.peak_q = num();
      else if (a == "--shoulder") P.shoulder_factor = num();
      else if (a == "--n-perm") P.n_perm = (int)num();
      else if (a == "--seed") P.seed = (uint64_t)num();
      else if (a == "--null-positions") P.null_positions = (int)num();
      else if (a == "--recleave-every") P.recleave_every = (int)num();
      else if (a == "--ubp-outlier-k") P.ubp_outlier_k = num();
      else if (a == "--beta-fallback-loss") P.beta_fallback_loss = num();
      else if (a == "--beta-fallback-gain") P.beta_fallback_gain = num();
      else if (a == "--marker-sd") P.marker_sd = num();
      else if (a == "--sd-c") P.sd_from_breakpoints_c = num();
      else if (a == "--grid-spacing") P.grid_spacing = num();
      else if (a == "--error-norm") { std::string v = val(); if (v != "rms" && v != "sem") throw std::runtime_error("--error-norm must be rms or sem"); P.error_sem = v == "sem"; }
      else if (a == "--split-at-centromere") P.split_at_centromere = true;
      else if (a == "--density-noise") { std::string v = val(); if (v != "sem" && v != "sd") throw std::runtime_error("--density-noise must be sem or sd"); P.density_noise_sem = v == "sem"; }
      else if (a == "--peak-window-lo") P.peak_window_lo = num();
      else if (a == "--peak-window-hi") P.peak_window_hi = num();
      else if (a == "--call-threshold") P.call_threshold = num();
      else if (a == "--amp-threshold") P.amp_threshold = num();
      else if (a == "--burden-threshold") P.burden_threshold = num();
      else if (a == "--fga-threshold") P.fga_threshold = num();
      else if (a == "--fga-threshold-loss") P.fga_threshold_loss = num();
      else if (a == "--fga-assembly") P.fga_assembly = true;
      else if (a == "--ploidy-prior") P.ploidy_prior = true;
      else if (a == "--ploidy-prior-wgd") P.ploidy_prior_wgd = num();
      else if (a == "--ploidy-prior-mu1") P.ploidy_prior_mu1 = num();
      else if (a == "--ploidy-prior-mu2") P.ploidy_prior_mu2 = num();
      else if (a == "--ploidy-prior-sd") P.ploidy_prior_sd = num();
      else if (a == "--ploidy-prior-ambig") P.ploidy_prior_ambig = num();
      else if (a == "--arm-baseline") { std::string v = val(); if (v != "median" && v != "twolevel") throw std::runtime_error("--arm-baseline must be median or twolevel"); P.arm_baseline = v == "twolevel"; }
      else if (a == "--arm-min-markers") P.arm_min_markers = (int)num();
      else if (a == "--focal-only") P.focal_only = true;
      else if (a == "--focal-max-mb") P.focal_max_bp = num() * 1e6;
      else if (a == "--arm-test") P.arm_test = true;
      else if (a == "--arm-perm") P.arm_n_perm = (int)num();
      else if (a == "--peak-independence") P.peak_independence = true;
      else if (a == "--peak-run-min") P.peak_run_min = num();
      else if (a.size() > 1 && a[0] == '-' && a[1] == '-') throw std::runtime_error("unknown option " + a + " (see --help)");
      else if (seg.empty()) seg = a;
      else throw std::runtime_error("unexpected argument " + a);
    }
    if (seg.empty()) throw std::runtime_error("no .seg file given");
    if (out.empty()) {
      out = seg;
      size_t dot = out.find_last_of('.'), slash = out.find_last_of('/');
      if (dot != std::string::npos && (slash == std::string::npos || dot > slash)) out = out.substr(0, dot);
      out += ".rae";
    }
    auto t0 = std::chrono::steady_clock::now();
    rae::InputLog log;
    std::vector<rae::SegRow> rows = rae::read_seg(seg, log);
    std::vector<rae::BinRow> b;
    if (!bins.empty()) b = rae::read_bins(bins);
    std::map<std::string, rae::StatRow> st, nst;
    if (!stats.empty()) st = rae::read_sample_stats(stats);
    if (!pp.empty())  // purity and ploidy per sample, merged into the sample statistics  [G-13a]
      for (const auto& kv : rae::read_sample_stats(pp)) { st[kv.first].purity = kv.second.purity; st[kv.first].ploidy = kv.second.ploidy; }
    std::vector<rae::SegRow> nrows;
    rae::InputLog nlog;
    if (!normals.empty()) nrows = rae::read_seg(normals, nlog);
    if (!nstats.empty()) nst = rae::read_sample_stats(nstats);
    rae::Result R = rae::run(rows, b, st, build, P, log, nrows, nst, cnv);
    std::vector<std::string> files = rae::write_outputs(R, out, seg);
    double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    int nsig = 0;
    for (int r = 0; r < R.ubp.size(); ++r) nsig += (R.qA[r] <= P.fdr) || (R.qD[r] <= P.fdr);
    std::cerr << "rae: " << R.S.size() << " samples, " << R.ubp.size() << " regions, " << nsig << " significant, "
              << R.roi.size() << " ROI rows, build " << R.G.build << ", " << secs << " s\n";
    for (const std::string& w : R.log.warnings) std::cerr << "rae: warning: " << w << "\n";
    for (const std::string& f : files) std::cout << f << "\n";
  } catch (const std::exception& e) {
    std::cerr << "rae: error: " << e.what() << "\n";
    return 2;
  }
  return 0;
}
