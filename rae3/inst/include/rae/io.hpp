// RAE 3.0 core -- input: .seg rows, bins, sample statistics; genome and marker-grid construction.
#pragma once
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "annotation_data.hpp"
#include "cnv_data.hpp"
#include "density.hpp"
#include "types.hpp"

namespace rae {

struct SegRow { std::string id; int chrom; int64_t start, end; int nmark; double mean, sd; };
struct BinRow { int chrom; int64_t start, end; };
struct StatRow { double dn = NaN, q01 = NaN, q025 = NaN, purity = NaN, ploidy = NaN; };
// Limits on values taken from input files, so a malformed file cannot drive allocation.
const double MAX_COORD = 1e10;        // larger than any chromosome
const double MAX_NUM_MARK = 5e7;      // markers in one segment
const size_t MAX_GRID_MARKERS = 20000000;  // markers in a reconstructed grid
inline bool row_in_range(double st, double en, double nm) {
  return st >= 0 && en >= st && en <= MAX_COORD && (is_nan(nm) || (nm >= 0 && nm <= MAX_NUM_MARK));
}
struct InputLog {
  int rows_read = 0, rows_non_autosomal = 0, rows_bad = 0, segs_dropped_unmappable = 0, segs_dropped_duplicate = 0;
  bool has_sd = false;
  std::vector<std::string> warnings;
};

// "chr7" / "7" -> 6 ; anything that is not an autosome -> -1.  [INFERRED G-11a]
inline int parse_chrom(std::string s) {
  if (s.size() > 3 && (s.compare(0, 3, "chr") == 0 || s.compare(0, 3, "Chr") == 0)) s = s.substr(3);
  if (s.empty() || s.size() > 2) return -1;
  for (char c : s) if (c < '0' || c > '9') return -1;
  int v = std::atoi(s.c_str());
  return (v >= 1 && v <= 22) ? v - 1 : -1;
}

inline std::vector<std::string> split_fields(const std::string& line) {
  std::vector<std::string> f;
  std::string cur;
  bool tab = line.find('\t') != std::string::npos;
  for (char c : line) {
    if (c == '\r') continue;
    if (tab ? c == '\t' : (c == ' ' || c == ',')) { f.push_back(cur); cur.clear(); }
    else cur += c;
  }
  f.push_back(cur);
  return f;
}
inline bool to_double(const std::string& s, double& v) {
  if (s.empty() || s == "NA" || s == "NaN" || s == "nan") { v = NaN; return true; }
  char* e = nullptr;
  v = std::strtod(s.c_str(), &e);
  return e && *e == '\0';
}

// Standard six-column seg (ID, chrom, loc.start, loc.end, num.mark, seg.mean), optional seventh
// column seg.sd. A header row is detected by a non-numeric start column.
inline std::vector<SegRow> read_seg(const std::string& path, InputLog& log) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open seg file: " + path);
  std::vector<SegRow> rows;
  std::string line;
  bool first = true;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> f = split_fields(line);
    if (f.size() < 6) { if (!first) log.rows_bad++; first = false; continue; }
    double st, en, nm, mean, sd = NaN;
    bool ok = to_double(f[2], st) && to_double(f[3], en) && to_double(f[4], nm) && to_double(f[5], mean);
    if (first && !ok) { first = false; continue; }  // header
    first = false;
    if (!ok || is_nan(st) || is_nan(en) || is_nan(mean) || !std::isfinite(mean) || !row_in_range(st, en, nm)) { log.rows_bad++; continue; }
    if (!is_nan(sd) && !(sd >= 0 && std::isfinite(sd))) sd = NaN;
    if (f.size() >= 7 && to_double(f[6], sd) && !is_nan(sd)) log.has_sd = true;
    log.rows_read++;
    int c = parse_chrom(f[1]);
    if (c < 0) { log.rows_non_autosomal++; continue; }
    rows.push_back({f[0], c, (int64_t)st, (int64_t)en, is_nan(nm) ? 1 : (int)nm, mean, sd});
  }
  return rows;
}

inline std::vector<BinRow> read_bins(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open bins file: " + path);
  std::vector<BinRow> rows;
  std::string line;
  while (std::getline(in, line)) {
    if (line.empty() || line[0] == '#') continue;
    std::vector<std::string> f = split_fields(line);
    double st, en;
    if (f.size() < 3 || !to_double(f[1], st) || !to_double(f[2], en) || is_nan(st) || is_nan(en)) continue;  // header / bad
    if (!row_in_range(st, en, NaN)) continue;
    int c = parse_chrom(f[0]);
    if (c >= 0) rows.push_back({c, (int64_t)st, (int64_t)en});
  }
  return rows;
}

// ID, dn, q01, q025, and optionally purity, ploidy (header required; columns located by name; any may be absent).
inline std::map<std::string, StatRow> read_sample_stats(const std::string& path) {
  std::ifstream in(path);
  if (!in) throw std::runtime_error("cannot open sample-stats file: " + path);
  std::map<std::string, StatRow> m;
  std::string line;
  int ci = 0, cd = -1, c1 = -1, c25 = -1, cp = -1, cy = -1;
  bool header = true;
  while (std::getline(in, line)) {
    if (line.empty()) continue;
    std::vector<std::string> f = split_fields(line);
    if (header) {
      for (size_t i = 0; i < f.size(); ++i) {
        if (f[i] == "ID" || f[i] == "id" || f[i] == "sample") ci = (int)i;
        else if (f[i] == "dn") cd = (int)i;
        else if (f[i] == "q01") c1 = (int)i;
        else if (f[i] == "q025") c25 = (int)i;
        else if (f[i] == "purity") cp = (int)i;
        else if (f[i] == "ploidy") cy = (int)i;
      }
      header = false;
      continue;
    }
    if (ci >= (int)f.size()) continue;  // short row: no ID field
    StatRow s;
    if (cd >= 0 && cd < (int)f.size()) to_double(f[cd], s.dn);
    if (c1 >= 0 && c1 < (int)f.size()) to_double(f[c1], s.q01);
    if (c25 >= 0 && c25 < (int)f.size()) to_double(f[c25], s.q025);
    if (cp >= 0 && cp < (int)f.size()) to_double(f[cp], s.purity);
    if (cy >= 0 && cy < (int)f.size()) to_double(f[cy], s.ploidy);
    m[f[ci]] = s;
  }
  return m;
}

inline const annot::BuildData* find_build(const std::string& name) {
  std::string n = name;
  if (n == "GRCh37" || n == "b37") n = "hg19";
  if (n == "GRCh38") n = "hg38";
  if (n == "hs1" || n == "t2t" || n == "T2T-CHM13" || n == "chm13v2") n = "chm13";
  for (int i = 0; i < annot::n_builds; ++i) if (n == annot::builds[i].name) return &annot::builds[i];
  return nullptr;
}

// Pick the bundled build whose chromosome lengths fit the data best ("pulled from the data").
// A build is consistent if no coordinate exceeds its chromosome length; among consistent builds the
// one with the least total slack between chromosome ends and the largest coordinate wins.
inline std::string detect_build(const std::vector<SegRow>& rows) {
  int64_t mx[22] = {0};
  for (const SegRow& r : rows) mx[r.chrom] = std::max(mx[r.chrom], r.end);
  std::string best = "";
  double best_slack = 1e300;
  for (int b = 0; b < annot::n_builds; ++b) {
    const annot::BuildData& B = annot::builds[b];
    if (std::string(B.name) == "hg17") continue;  // kept for reproducing the paper, never auto-selected
    bool ok = true;
    double slack = 0;
    for (int c = 0; c < 22; ++c) {
      if (mx[c] > B.chrom_len[c]) { ok = false; break; }
      if (mx[c] > 0) slack += (double)(B.chrom_len[c] - mx[c]);
    }
    if (ok && slack < best_slack) { best_slack = slack; best = B.name; }
  }
  return best;
}

// Population CNV catalogue. `spec`: "default" = the catalogue bundled for the build (DGV Gold Standard
// for hg19, gnomAD-SV v4.1 for hg38, gnomAD-SV v4.1 lifted for CHM13); "none"; or a file with columns
// chrom, start, end, type (gain / loss / anything else = unknown) and optionally allele frequency.
// [STATED: M&M "Human structural variation"]: overlap is reported; nothing is excluded.
inline void load_cnv(Genome& G, const std::string& spec, InputLog& log) {
  G.cnv.clear();
  if (spec == "none") { G.cnv_source = "none"; return; }
  if (spec == "default") {
    for (int i = 0; i < annot::n_cnv_sets; ++i) {
      const annot::CnvSet& S = annot::cnv_sets[i];
      if (G.build != S.build) continue;
      for (int k = 0; k < S.n; ++k) G.cnv.push_back({S.chrom[k] - 1, S.start[k], S.end[k], S.type[k], S.af[k]});
      G.cnv_source = S.source;
    }
    if (G.cnv.empty()) G.cnv_source = "none (no catalogue bundled for build " + G.build + ")";
  } else {
    std::ifstream in(spec);
    if (!in) throw std::runtime_error("cannot open CNV catalogue: " + spec);
    std::string line;
    while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#') continue;
      std::vector<std::string> f = split_fields(line);
      double st, en, af = NaN;
      if (f.size() < 3 || !to_double(f[1], st) || !to_double(f[2], en) || is_nan(st) || is_nan(en) || !row_in_range(st, en, NaN)) continue;
      int c = parse_chrom(f[0]);
      if (c < 0) continue;
      int8_t t = 0;
      if (f.size() > 3) { if (f[3] == "gain" || f[3] == "DUP" || f[3] == "Gain") t = 1; else if (f[3] == "loss" || f[3] == "DEL" || f[3] == "Loss") t = -1; }
      if (f.size() > 4) to_double(f[4], af);
      G.cnv.push_back({c, (int64_t)st, (int64_t)en, t, (float)af});
    }
    G.cnv_source = spec;
    if (G.cnv.empty()) log.warnings.push_back("CNV catalogue " + spec + " has no usable autosomal rows");
  }
  std::sort(G.cnv.begin(), G.cnv.end(), [](const CnvRow& a, const CnvRow& b) { return a.chrom != b.chrom ? a.chrom < b.chrom : a.start < b.start; });
}

// Overlap of one interval with the catalogue: variants overlapping, and the fraction of the interval
// covered by gain-type, loss-type, unknown-type and any variants (unions, so overlaps are not double counted).
struct CnvOverlap { int n = 0; double gain = 0, loss = 0, unknown = 0, any = 0; };
inline CnvOverlap cnv_overlap(const Genome& G, int chrom, int64_t s, int64_t e) {
  CnvOverlap o;
  std::vector<std::pair<int64_t, int64_t>> iv[4];
  for (const CnvRow& r : G.cnv) {  // catalogues are small; a linear scan per region of interest is enough
    if (r.chrom < chrom) continue;
    if (r.chrom > chrom) break;
    if (r.end < s || r.start > e) continue;
    std::pair<int64_t, int64_t> p(std::max(r.start, s), std::min(r.end, e));
    ++o.n;
    iv[r.type > 0 ? 0 : (r.type < 0 ? 1 : 2)].push_back(p);
    iv[3].push_back(p);
  }
  auto cover = [&](std::vector<std::pair<int64_t, int64_t>>& v) {
    std::sort(v.begin(), v.end());
    int64_t tot = 0, cs = 0, ce = -1;
    for (auto& p : v) {
      if (p.first > ce) { if (ce >= cs) tot += ce - cs + 1; cs = p.first; ce = p.second; }
      else ce = std::max(ce, p.second);
    }
    if (ce >= cs) tot += ce - cs + 1;
    return (double)tot / (double)(e - s + 1);
  };
  o.gain = cover(iv[0]); o.loss = cover(iv[1]); o.unknown = cover(iv[2]); o.any = cover(iv[3]);
  return o;
}

// Genome with marker grid. `bins` empty -> virtual uniform grid from the data.  [INFERRED G-01d]
inline Genome make_genome(const std::string& build, const std::vector<SegRow>& rows,
                          const std::vector<BinRow>& bins, const Params& P, InputLog& log) {
  Genome G;
  G.chrom_len.assign(22, 0); G.cen_start.assign(22, 0); G.cen_end.assign(22, 0);
  G.hotspots.assign(22, {});
  std::string b = build;
  if (b == "auto") {
    b = detect_build(rows);
    if (b.empty()) { b = "custom"; log.warnings.push_back("no bundled build is consistent with the coordinates; running as custom genome"); }
  }
  const annot::BuildData* B = find_build(b);
  if (B) {
    G.build = B->name;
    for (int c = 0; c < 22; ++c) { G.chrom_len[c] = B->chrom_len[c]; G.cen_start[c] = B->cen_start[c]; G.cen_end[c] = B->cen_end[c]; }
    for (int i = 0; i < B->n_hot; ++i) G.hotspots[B->hot_chrom[i] - 1].push_back(B->hot_pos[i]);
    G.has_hotspots = true;
    if (G.build == "hg17")
      log.warnings.push_back("hg17 is kept only for reproducing paper-era data during development; do not use it for current analyses");
    int over = 0;
    for (const SegRow& r : rows) if (r.end > G.chrom_len[r.chrom]) ++over;
    if (over > 0) log.warnings.push_back(std::to_string(over) + " segments end beyond the chromosome length of build " + G.build + "; check --build");
  } else {
    G.build = "custom";
    for (const SegRow& r : rows) G.chrom_len[r.chrom] = std::max(G.chrom_len[r.chrom], r.end);
    log.warnings.push_back("custom genome: no centromere or hotspot table; uniform cleavage sites used in the null model [G-08k]");
  }
  G.chrom_first.assign(23, 0);
  if (!bins.empty()) {
    std::vector<BinRow> s(bins);
    std::sort(s.begin(), s.end(), [](const BinRow& a, const BinRow& c) { return a.chrom != c.chrom ? a.chrom < c.chrom : a.start < c.start; });
    for (const BinRow& r : s) { G.m_chrom.push_back(r.chrom); G.m_start.push_back(r.start); G.m_end.push_back(r.end); }
  } else {
    G.virtual_grid = true;
    double d = P.grid_spacing;
    if (is_nan(d)) {
      std::vector<double> v, all;
      for (const SegRow& r : rows) {
        if (r.nmark < 1) continue;
        double x = (double)(r.end - r.start + 1) / r.nmark;
        all.push_back(x);
        if (r.nmark >= 10) v.push_back(x);
      }
      d = median(v.empty() ? all : v);
    }
    int64_t di = std::max<int64_t>(1, (int64_t)std::llround(d));
    G.grid_spacing = (double)di;
    int64_t lo[22], hi[22];
    for (int c = 0; c < 22; ++c) { lo[c] = INT64_MAX; hi[c] = -1; }
    for (const SegRow& r : rows) {
      lo[r.chrom] = std::min(lo[r.chrom], r.start / di);
      hi[r.chrom] = std::max(hi[r.chrom], (r.end + di - 1) / di - 1);
    }
    size_t total = 0;
    for (int c = 0; c < 22; ++c) if (hi[c] >= 0) total += (size_t)(hi[c] - lo[c] + 1);
    if (total > MAX_GRID_MARKERS)
      throw std::runtime_error("reconstructed marker grid would have " + std::to_string(total) +
                               " markers (spacing " + std::to_string(di) + " bp); supply --bins or --grid-spacing");
    for (int c = 0; c < 22; ++c)
      for (int64_t k = lo[c]; hi[c] >= 0 && k <= hi[c]; ++k) {
        G.m_chrom.push_back(c); G.m_start.push_back(k * di + 1); G.m_end.push_back((k + 1) * di);
      }
  }
  for (int c = 0; c < 22; ++c) {
    G.chrom_first[c] = (int32_t)(std::lower_bound(G.m_chrom.begin(), G.m_chrom.end(), c) - G.m_chrom.begin());
  }
  G.chrom_first[22] = (int32_t)G.m_chrom.size();
  return G;
}

// Map rows onto the marker index and assemble samples (order of first appearance).
inline std::vector<Sample> make_samples(const std::vector<SegRow>& rows, const Genome& G,
                                        const std::map<std::string, StatRow>& stats, InputLog& log) {
  std::vector<Sample> S;
  std::map<std::string, int> idx;
  for (const SegRow& r : rows) {
    auto it = idx.find(r.id);
    int j;
    if (it == idx.end()) { j = (int)S.size(); idx[r.id] = j; S.emplace_back(); S.back().id = r.id; }
    else j = it->second;
    int lo = G.chrom_first[r.chrom], hi = G.chrom_first[r.chrom + 1];
    // first marker ending at or after loc.start; last marker starting at or before loc.end
    int a = (int)(std::lower_bound(G.m_end.begin() + lo, G.m_end.begin() + hi, r.start) - G.m_end.begin());
    int b = (int)(std::upper_bound(G.m_start.begin() + lo, G.m_start.begin() + hi, r.end) - G.m_start.begin()) - 1;
    if (a >= hi || b < lo || a > b) { log.segs_dropped_unmappable++; continue; }
    Segment g;
    g.start = a; g.end = b; g.chrom = r.chrom; g.bp_start = r.start; g.bp_end = r.end;
    g.nmark = std::max(1, r.nmark); g.mean = r.mean; g.sd = r.sd;
    S[j].segs.push_back(g);
  }
  for (Sample& s : S) {
    std::sort(s.segs.begin(), s.segs.end(), [](const Segment& a, const Segment& b) { return a.start < b.start; });
    std::vector<Segment> keep;
    for (const Segment& g : s.segs) {
      if (!keep.empty() && keep.back().start == g.start) { log.segs_dropped_duplicate++; continue; }
      keep.push_back(g);
    }
    // Each segment extends to the marker before the next segment's start; the first one on a
    // chromosome extends back to the chromosome start, the last one to its end.  [INFERRED G-07d]
    for (size_t i = 0; i < keep.size(); ++i) {
      bool first_on_chrom = i == 0 || keep[i - 1].chrom != keep[i].chrom;
      bool last_on_chrom = i + 1 == keep.size() || keep[i + 1].chrom != keep[i].chrom;
      if (first_on_chrom) keep[i].start = G.chrom_first[keep[i].chrom];
      keep[i].end = last_on_chrom ? G.chrom_first[keep[i].chrom + 1] - 1 : keep[i + 1].start - 1;
    }
    s.segs.swap(keep);
    auto it = stats.find(s.id);
    if (it != stats.end()) {
      s.dn = it->second.dn; s.q01 = it->second.q01; s.q025 = it->second.q025;
      s.in_purity = it->second.purity; s.in_ploidy = it->second.ploidy;
    }
  }
  return S;
}

}  // namespace rae
