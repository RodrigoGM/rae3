// RAE 3.0 core -- figures: dependency-free inline SVG collected in one self-contained HTML report.
// The figure set is fixed, whatever the input, and mirrors the paper's figures:
//   genome-wide significance (Fig. 5), alteration frequency, summary scores, ROI zooms with error
//   bars (Fig. 4C), per-sample density with detectors (Fig. 2C / S5), region-size histogram (Fig. 3C),
//   cohort E_k distribution (Fig. S10).
// Included from report.hpp.
#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include "pipeline.hpp"

namespace rae {

namespace fig {
// colours: gain / loss are a two-pole pair; ink and grid are recessive
static const char* const GAIN = "#eb6834";
static const char* const LOSS = "#2a78d6";
static const char* const INK = "#0b0b0b";
static const char* const INK2 = "#52514e";
static const char* const MUTED = "#898781";
static const char* const GRID = "#e1e0d9";
static const char* const SURF = "#fcfcfb";

// Everything that originates in the input (sample IDs, file names, warnings) is escaped before it
// is written into HTML or SVG.
inline std::string esc(const std::string& in) {
  std::string o;
  o.reserve(in.size());
  for (char c : in) {
    switch (c) {
      case '&': o += "&amp;"; break;
      case '<': o += "&lt;"; break;
      case '>': o += "&gt;"; break;
      case '"': o += "&quot;"; break;
      case '\'': o += "&#39;"; break;
      default: o += c;
    }
  }
  return o;
}
inline std::string n(double v) { char b[32]; std::snprintf(b, sizeof b, "%.1f", v); return b; }
inline std::string g(double v) { char b[32]; std::snprintf(b, sizeof b, "%.3g", v); return b; }

struct Svg {
  std::ostringstream o;
  double W, H;
  Svg(double w, double h, const std::string& label) : W(w), H(h) {
    o << "<svg xmlns='http://www.w3.org/2000/svg' viewBox='0 0 " << w << ' ' << h << "' width='100%' role='img' aria-label='"
      << label << "' style='max-width:" << w << "px;background:" << SURF << ";font-family:system-ui,sans-serif'>\n";
  }
  void line(double x1, double y1, double x2, double y2, const char* col, double sw = 1, const char* dash = "") {
    o << "<line x1='" << n(x1) << "' y1='" << n(y1) << "' x2='" << n(x2) << "' y2='" << n(y2) << "' stroke='" << col
      << "' stroke-width='" << sw << "'" << (*dash ? std::string(" stroke-dasharray='") + dash + "'" : std::string()) << "/>\n";
  }
  void text(double x, double y, const std::string& s, const char* col = INK2, int size = 11, const char* anchor = "start", int rot = 0) {
    o << "<text x='" << n(x) << "' y='" << n(y) << "' fill='" << col << "' font-size='" << size << "' text-anchor='" << anchor << "'";
    if (rot) o << " transform='rotate(" << rot << ' ' << n(x) << ' ' << n(y) << ")'";
    o << ">" << s << "</text>\n";
  }
  void rect(double x, double y, double w, double h, const char* fill, const std::string& title = "", double op = 1) {
    o << "<rect x='" << n(x) << "' y='" << n(y) << "' width='" << n(std::max(0.0, w)) << "' height='" << n(std::max(0.0, h))
      << "' fill='" << fill << "' fill-opacity='" << op << "'>";
    if (!title.empty()) o << "<title>" << title << "</title>";
    o << "</rect>\n";
  }
  void path(const std::string& d, const char* stroke, const char* fill = "none", double sw = 1.5, double fop = 1) {
    o << "<path d='" << d << "' stroke='" << stroke << "' fill='" << fill << "' stroke-width='" << sw << "' fill-opacity='" << fop
      << "' stroke-linejoin='round'/>\n";
  }
  std::string str() { return o.str() + "</svg>\n"; }
};

// Genome x-axis: chromosomes laid end to end.
struct GenomeAxis {
  std::vector<double> off; double total = 0;
  explicit GenomeAxis(const Genome& G) {
    off.assign(G.n_chrom + 1, 0);
    for (int c = 0; c < G.n_chrom; ++c) {
      double len = (double)G.chrom_len[c];
      if (len <= 0 && G.chrom_first[c] < G.chrom_first[c + 1]) len = (double)G.m_end[G.chrom_first[c + 1] - 1];
      off[c + 1] = off[c] + len;
    }
    total = off[G.n_chrom];
  }
};

// Mirrored step plot along the genome: `up` drawn above the baseline, `down` below.
inline std::string genome_mirror(const Result& R, const std::vector<double>& up, const std::vector<double>& down,
                                 double ymax, const std::string& ylab, const std::string& label,
                                 double ref = NaN, const std::string& reflab = "") {
  const double W = 1100, H = 300, L = 62, Rm = 12, T = 22, B = 30;
  Svg s(W, H, label);
  GenomeAxis ax(R.G);
  const double pw = W - L - Rm, mid = T + (H - T - B) / 2, half = (H - T - B) / 2;
  auto X = [&](int c, double bp) { return L + pw * (ax.off[c] + bp) / ax.total; };
  for (int c = 0; c < R.G.n_chrom; ++c) {
    if (c % 2 == 1) s.rect(X(c, 0), T, X(c + 1, 0) - X(c, 0), H - T - B, "#f1f0ea");
    s.text((X(c, 0) + X(c + 1, 0)) / 2, H - B + 14, std::to_string(c + 1), MUTED, 10, "middle");
  }
  for (int t = 1; t <= 2; ++t) {
    double v = ymax * t / 2;
    s.line(L, mid - half * t / 2, W - Rm, mid - half * t / 2, GRID);
    s.line(L, mid + half * t / 2, W - Rm, mid + half * t / 2, GRID);
    s.text(L - 6, mid - half * t / 2 + 4, g(v), MUTED, 10, "end");
    s.text(L - 6, mid + half * t / 2 + 4, g(v), MUTED, 10, "end");
  }
  for (int dir = 0; dir < 2; ++dir) {
    const std::vector<double>& v = dir == 0 ? up : down;
    std::ostringstream d;
    int cur = -1;
    for (int r = 0; r < R.ubp.size(); ++r) {
      int c = R.ubp.chrom[r];
      double x0 = X(c, (double)reg_bp_start(R, r)), x1 = X(c, (double)reg_bp_end(R, r));
      double h = half * std::min(1.0, std::max(0.0, v[r] / ymax));
      double y = dir == 0 ? mid - h : mid + h;
      if (c != cur) { if (cur >= 0) d << "V" << n(mid) << "Z "; d << "M" << n(x0) << ' ' << n(mid) << ' '; cur = c; }
      d << "V" << n(y) << "H" << n(x1);
    }
    if (cur >= 0) d << "V" << n(mid) << "Z";
    s.path(d.str(), "none", dir == 0 ? GAIN : LOSS, 0, 0.9);
  }
  s.line(L, mid, W - Rm, mid, INK2);
  if (!is_nan(ref)) {
    double h = half * ref / ymax;
    s.line(L, mid - h, W - Rm, mid - h, INK2, 1, "4 3");
    s.line(L, mid + h, W - Rm, mid + h, INK2, 1, "4 3");
    s.text(W - Rm - 2, mid - h - 3, reflab, INK2, 10, "end");
  }
  s.text(14, mid, ylab, INK2, 11, "middle", -90);
  s.rect(L + 4, 6, 10, 10, GAIN); s.text(L + 18, 15, "gain / amplification", INK2, 11);
  s.rect(L + 150, 6, 10, 10, LOSS); s.text(L + 164, 15, "loss / deletion", INK2, 11);
  return s.str();
}

inline std::string histogram_pair(const std::vector<double>& a, const std::vector<double>& b, const std::string& la,
                                  const std::string& lb) {
  const double W = 1100, H = 230, L = 50, T = 26, B = 34, gap = 40;
  Svg s(W, H, "Region size histograms");
  const int nb = 40;
  const double lo = 4.0, hi = 9.0;  // log10 bp
  const double pw = (W - L - 12 - gap) / 2;
  for (int k = 0; k < 2; ++k) {
    const std::vector<double>& v = k == 0 ? a : b;
    std::vector<int> h(nb, 0);
    for (double x : v) { int i = (int)((std::log10(std::max(1.0, x)) - lo) / (hi - lo) * nb); h[std::min(nb - 1, std::max(0, i))]++; }
    int mx = std::max(1, *std::max_element(h.begin(), h.end()));
    double x0 = L + k * (pw + gap);
    for (int i = 0; i < nb; ++i) {
      double bh = (H - T - B) * h[i] / (double)mx;
      s.rect(x0 + pw * i / nb + 1, H - B - bh, pw / nb - 2, bh, k == 0 ? INK2 : MUTED,
             g(std::pow(10, lo + (hi - lo) * i / nb)) + " bp: " + std::to_string(h[i]) + " regions");
    }
    s.line(x0, H - B, x0 + pw, H - B, INK2);
    for (int e = 4; e <= 9; ++e) {
      double x = x0 + pw * (e - lo) / (hi - lo);
      static const char* lab[] = {"10 kb", "100 kb", "1 Mb", "10 Mb", "100 Mb", "1 Gb"};
      s.text(x, H - B + 14, lab[e - 4], MUTED, 10, "middle");
    }
    s.text(x0, 16, (k == 0 ? la : lb) + " (n = " + std::to_string(v.size()) + ", tallest bin " + std::to_string(mx) + ")", INK2, 11);
  }
  return s.str();
}

inline std::string ek_distribution(const Result& R) {
  const double W = 1100, H = 260, L = 62, T = 22, B = 22;
  Svg s(W, H, "Per-sample E for single-copy gain and hemizygous loss");
  std::vector<double> a, d;
  for (const Sample& x : R.S) { a.push_back(x.E[A0]); d.push_back(-x.E[D0]); }
  std::sort(a.begin(), a.end()); std::sort(d.begin(), d.end());
  double ymax = std::max(0.2, std::max(a.back(), d.back()));
  const double mid = T + (H - T - B) / 2, half = (H - T - B) / 2, pw = W - L - 12;
  const size_t N = a.size();
  for (int t = 1; t <= 2; ++t) {
    s.line(L, mid - half * t / 2, W - 12, mid - half * t / 2, GRID);
    s.line(L, mid + half * t / 2, W - 12, mid + half * t / 2, GRID);
    s.text(L - 6, mid - half * t / 2 + 4, g(ymax * t / 2), MUTED, 10, "end");
    s.text(L - 6, mid + half * t / 2 + 4, g(-ymax * t / 2), MUTED, 10, "end");
  }
  for (int k = 0; k < 2; ++k) {
    const std::vector<double>& v = k == 0 ? a : d;
    std::ostringstream p;
    p << "M" << n(L) << ' ' << n(mid);
    for (size_t i = 0; i < N; ++i) {
      double y = k == 0 ? mid - half * v[i] / ymax : mid + half * v[i] / ymax;
      p << "V" << n(y) << "H" << n(L + pw * (i + 1) / N);
    }
    p << "V" << n(mid) << "Z";
    s.path(p.str(), "none", k == 0 ? GAIN : LOSS, 0, 0.9);
  }
  s.line(L, mid, W - 12, mid, INK2);
  s.text(14, mid, "E (log2 ratio)", INK2, 11, "middle", -90);
  s.text(L + 4, 15, "E for single-copy gain (A0), samples sorted; median " + g(quantile7_sorted(a, 0.5)), INK2, 11);
  s.text(L + 4, H - 6, "E for hemizygous loss (D0), samples sorted; median " + g(-quantile7_sorted(d, 0.5)), INK2, 11);
  return s.str();
}

inline std::string sample_density(const Sample& x) {
  const double W = 264, H = 170, L = 8, T = 30, B = 22;
  Svg s(W, H, "Density and detectors for " + esc(x.id));
  const double pw = W - L - 8, ph = H - T - B;
  double ymax = *std::max_element(x.dens_y.begin(), x.dens_y.end());
  auto X = [&](double v) { return L + pw * (v + 1) / 2; };
  auto Y = [&](double f) { return T + ph * (1 - f); };
  s.line(X(0), T, X(0), H - B, GRID);
  s.line(X(x.Ha), T, X(x.Ha), H - B, MUTED, 1, "2 2");
  s.line(X(x.Hd), T, X(x.Hd), H - B, MUTED, 1, "2 2");
  std::ostringstream p;
  for (size_t i = 0; i < x.dens_x.size(); i += 4) p << (i ? "L" : "M") << n(X(x.dens_x[i])) << ' ' << n(Y(x.dens_y[i] / ymax));
  s.path(p.str(), INK, "none", 1.5);
  for (int k : {A0, D0, A1, D1}) {
    std::ostringstream q;
    for (int i = 0; i <= 200; ++i) {
      double v = -1 + 2.0 * i / 200;
      q << (i ? "L" : "M") << n(X(v)) << ' ' << n(Y(score_component(k, v, x.E, x.B, x.a1_logistic)));
    }
    s.o << "<path d='" << q.str() << "' stroke='" << ((k == A0 || k == A1) ? GAIN : LOSS) << "' fill='none' stroke-width='1.5'"
        << ((k == A1 || k == D1) ? " stroke-dasharray='5 3'" : "") << "/>\n";
  }
  s.line(L, H - B, W - 8, H - B, INK2);
  for (double t : {-1.0, -0.5, 0.0, 0.5, 1.0}) s.text(X(t), H - B + 13, g(t), MUTED, 9, "middle");
  s.text(L, 12, esc(x.id), INK, 11);
  s.text(L, 24, "E: D0 " + g(x.E[D0]) + "  A0 " + g(x.E[A0]), INK2, 10);
  return s.str();
}

inline std::string roi_zoom(const Result& R, const Roi& ev, const std::vector<const Roi*>& peaks) {
  const double W = 540, H = 230, L = 46, T = 30, B = 30;
  Svg s(W, H, "Region of interest");
  const std::vector<double>& sc = ev.direction > 0 ? R.st.Ap : R.st.Dp;
  const std::vector<double>& er = ev.direction > 0 ? R.st.errA : R.st.errD;
  const std::vector<double>& q = ev.direction > 0 ? R.qA : R.qD;
  const char* col = ev.direction > 0 ? GAIN : LOSS;
  int c = R.ubp.chrom[ev.r_first];
  int a = ev.r_first, b = ev.r_last;
  // show up to 3 flanking regions for context
  for (int k = 0; k < 3 && a > 0 && R.ubp.chrom[a - 1] == c; ++k) --a;
  for (int k = 0; k < 3 && b + 1 < R.ubp.size() && R.ubp.chrom[b + 1] == c; ++k) ++b;
  double x0 = (double)reg_bp_start(R, a), x1 = (double)reg_bp_end(R, b);
  double ymax = 0;
  for (int r = a; r <= b; ++r) ymax = std::max(ymax, sc[r] + er[r]);
  ymax = std::max(ymax * 1.08, 1e-6);
  const double pw = W - L - 10, ph = H - T - B;
  auto X = [&](double bp) { return L + pw * (bp - x0) / std::max(1.0, x1 - x0); };
  auto Y = [&](double v) { return T + ph * (1 - v / ymax); };
  for (int t = 0; t <= 2; ++t) { s.line(L, Y(ymax * t / 2), W - 10, Y(ymax * t / 2), GRID); s.text(L - 5, Y(ymax * t / 2) + 4, g(ymax * t / 2), MUTED, 10, "end"); }
  for (const Roi* p : peaks) s.rect(X((double)reg_bp_start(R, p->r_first)), T, X((double)reg_bp_end(R, p->r_last)) - X((double)reg_bp_start(R, p->r_first)), ph, col, "peak", 0.15);
  std::ostringstream d;
  for (int r = a; r <= b; ++r) {
    double xm = X(0.5 * ((double)reg_bp_start(R, r) + (double)reg_bp_end(R, r)));
    d << (r == a ? "M" : "L") << n(xm) << ' ' << n(Y(sc[r]));
  }
  s.path(d.str(), MUTED, "none", 1);
  for (int r = a; r <= b; ++r) {
    double xm = X(0.5 * ((double)reg_bp_start(R, r) + (double)reg_bp_end(R, r)));
    bool sig = q[r] <= R.P.fdr;
    s.line(xm, Y(std::min(ymax, sc[r] + er[r])), xm, Y(std::max(0.0, sc[r] - er[r])), sig ? col : MUTED, 1.5);
    s.o << "<circle cx='" << n(xm) << "' cy='" << n(Y(sc[r])) << "' r='4' fill='" << (sig ? col : MUTED) << "' stroke='" << SURF
        << "' stroke-width='1.5'><title>chr" << c + 1 << ":" << reg_bp_start(R, r) << "-" << reg_bp_end(R, r) << " score " << g(sc[r])
        << " +/- " << g(er[r]) << " q " << g(q[r]) << "</title></circle>\n";
  }
  s.line(L, H - B, W - 10, H - B, INK2);
  s.text(L, H - B + 14, n(x0 / 1e6) + " Mb", MUTED, 10);
  s.text(W - 10, H - B + 14, n(x1 / 1e6) + " Mb", MUTED, 10, "end");
  s.text(L, 14, std::string(ev.direction > 0 ? "Gain" : "Loss") + " chr" + std::to_string(c + 1) + ":" +
                    n(reg_bp_start(R, ev.r_first) / 1e6) + "-" + n(reg_bp_end(R, ev.r_last) / 1e6) + " Mb, " +
                    std::to_string(ev.n_regions) + " regions, " + std::to_string(peaks.size()) + " peak(s), min q " + g(ev.q), INK, 11);
  s.text(L, 26, std::string("summary score ") + (ev.direction > 0 ? "A'" : "D'") + " with analytical error; shaded = peak; grey = not significant", INK2, 10);
  return s.str();
}
}  // namespace fig

inline void write_report(const Result& R, const std::string& path, const std::string& input) {
  std::ofstream o(path);
  const Params& P = R.P;
  const int NR = R.ubp.size();
  std::vector<double> lqA(NR), lqD(NR), fg(NR, 0), fl(NR, 0), usz, csz;
  double qfloor = 1e-12;
  for (int r = 0; r < NR; ++r) qfloor = std::max(qfloor, 0.0);
  double lmax = 2.5;
  for (int r = 0; r < NR; ++r) {
    lqA[r] = -std::log10(std::max(R.qA[r], 1e-300)); lqD[r] = -std::log10(std::max(R.qD[r], 1e-300));
    lmax = std::max(lmax, std::max(lqA[r], lqD[r]));
    int n = 0, g = 0, l = 0;
    for (int j = 0; j < R.M.N; ++j) {
      float a0 = R.M.at(A0, r, j);
      if (a0 != a0) continue;
      ++n;
      int st = discrete_state(a0, R.M.at(A1, r, j), R.M.at(D0, r, j), R.M.at(D1, r, j), P);
      g += st > 0; l += st < 0;
    }
    if (n) { fg[r] = (double)g / n; fl[r] = (double)l / n; }
    usz.push_back((double)(reg_bp_end(R, r) - reg_bp_start(R, r) + 1));
  }
  for (int r = 0; r < R.cleaved_example.size(); ++r)
    csz.push_back((double)(R.G.m_end[R.cleaved_example.end[r]] - R.G.m_start[R.cleaved_example.start[r]] + 1));
  double fmax = 0.1, smax = 0.05;
  for (int r = 0; r < NR; ++r) { fmax = std::max(fmax, std::max(fg[r], fl[r])); smax = std::max(smax, std::max(R.st.Ap[r], R.st.Dp[r])); }
  int nev = 0, npk = 0;
  for (const Roi& x : R.roi) { if (x.is_peak) ++npk; else ++nev; }

  o << "<!DOCTYPE html>\n<html lang='en'><head><meta charset='utf-8'><title>RAE report</title>\n"
       "<meta name='viewport' content='width=device-width, initial-scale=1'>\n<style>\n"
       "body{font-family:system-ui,sans-serif;color:#0b0b0b;background:#f9f9f7;margin:0;padding:24px;line-height:1.45}\n"
       "main{max-width:1120px;margin:0 auto}h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 6px}\n"
       "p,li{font-size:13px;color:#52514e;margin:4px 0}table{border-collapse:collapse;font-size:12px;width:100%}\n"
       "th,td{text-align:right;padding:3px 8px;border-bottom:1px solid #e1e0d9}th{color:#52514e;font-weight:600}\n"
       "td:nth-child(-n+3),th:nth-child(-n+3){text-align:left}.tiles{display:flex;gap:12px;flex-wrap:wrap;margin:12px 0}\n"
       ".tile{background:#fcfcfb;border:1px solid #e1e0d9;border-radius:6px;padding:8px 14px}.tile b{display:block;font-size:20px}\n"
       ".tile span{font-size:12px;color:#52514e}.grid{display:flex;flex-wrap:wrap;gap:8px}.grid>svg{flex:0 0 auto}\n"
       ".scroll{overflow-x:auto}.note{background:#fcfcfb;border-left:3px solid #898781;padding:6px 10px}\n"
       "</style></head><body><main>\n";
  o << "<h1>RAE report</h1><p>" << fig::esc(input) << " &middot; genome " << R.G.build << " &middot; rae " << RAE_VERSION
    << " &middot; seed " << P.seed << "</p>\n";
  o << "<div class='tiles'>"
    << "<div class='tile'><b>" << R.S.size() << "</b><span>tumours</span></div>"
    << "<div class='tile'><b>" << NR << "</b><span>regions in the unified breakpoint profile</span></div>"
    << "<div class='tile'><b>" << nev << "</b><span>significant events (q &le; " << fig::g(P.fdr) << ")</span></div>"
    << "<div class='tile'><b>" << npk << "</b><span>peaks</span></div>"
    << "<div class='tile'><b>" << R.null.A.size() << "</b><span>null scores</span></div></div>\n";
  o << "<p class='note'>Reimplementation of RAE (Taylor et al. 2008, PLoS ONE 3:e3179) from the publication. Choices the paper does not "
       "determine are listed in the gap register; the ones that fired in this run are in the run.json manifest";
  if (R.G.virtual_grid) o << ". Marker positions were reconstructed from the segments (no bins file)";
  if (!R.G.has_hotspots) o << ". No hotspot table for this genome: uniform cleavage sites were used in the null model";
  o << ".</p>\n";
  for (const std::string& w : R.log.warnings) o << "<p class='note'>Warning: " << fig::esc(w) << "</p>\n";

  o << "<h2>Genome-wide significance</h2><p>&minus;log10 of the false discovery rate (q) per region of the unified breakpoint profile; "
       "dashed line: q = " << fig::g(P.fdr) << ". Values are capped by the smallest p-value the permutation can resolve.</p>\n"
    << fig::genome_mirror(R, lqA, lqD, lmax, "-log10 q", "Genome-wide significance", -std::log10(P.fdr), "q = " + fig::g(P.fdr));
  o << "<h2>Summary scores</h2><p>A&prime; (above) and D&prime; (below), the L2 norm of the cohort-mean component scores.</p>\n"
    << fig::genome_mirror(R, R.st.Ap, R.st.Dp, smax, "summary score", "Summary scores");
  o << "<h2>Alteration frequency</h2><p>Fraction of tumours called gained or lost in each region (component score &ge; "
    << fig::g(P.call_threshold) << "; thresholds from Taylor et al. 2010, not part of RAE 2008).</p>\n"
    << fig::genome_mirror(R, fg, fl, fmax, "fraction of tumours", "Alteration frequency");

  o << "<h2>Regions of interest</h2>\n";
  std::vector<const Roi*> events;
  for (const Roi& x : R.roi) if (!x.is_peak) events.push_back(&x);
  std::sort(events.begin(), events.end(), [](const Roi* a, const Roi* b) { return a->score > b->score; });
  o << "<p>The " << std::min<size_t>(12, events.size()) << " highest-scoring events of " << events.size() << ". All are listed in the roi.tsv table.</p><div class='grid'>\n";
  for (size_t i = 0; i < events.size() && i < 12; ++i) {
    std::vector<const Roi*> pk;
    for (const Roi& x : R.roi) if (x.is_peak && x.direction == events[i]->direction && x.event == events[i]->event) pk.push_back(&x);
    o << fig::roi_zoom(R, *events[i], pk);
  }
  o << "</div>\n<div class='scroll'><table><tr><th>Type</th><th>Direction</th><th>Location</th><th>Size (Mb)</th><th>Regions</th><th>Score</th><th>Error</th><th>q</th></tr>\n";
  int shown = 0;
  for (const Roi& x : R.roi) {
    if (++shown > 200) break;
    int c = R.ubp.chrom[x.r_first];
    int64_t s = reg_bp_start(R, x.r_first), e = reg_bp_end(R, x.r_last);
    o << "<tr><td>" << (x.is_peak ? "&nbsp;&nbsp;peak" : "event") << "</td><td>" << (x.direction > 0 ? "gain" : "loss") << "</td><td>chr"
      << c + 1 << ":" << s << "-" << e << "</td><td>" << fig::g((e - s + 1) / 1e6) << "</td><td>" << x.n_regions << "</td><td>"
      << fig::g(x.score) << "</td><td>" << fig::g(x.error) << "</td><td>" << fig::g(x.q) << "</td></tr>\n";
  }
  o << "</table></div>\n";
  if (R.roi.size() > 200) o << "<p>Table truncated at 200 rows; see roi.tsv.</p>\n";

  o << "<h2>Per-tumour detectors</h2><p>Twelve tumours spaced evenly through the input order. Black: density of noise-added segment means "
       "after mode centring; dotted: half-maximum points; orange: single-copy gain (solid) and amplification (dashed) detectors; "
       "blue: hemizygous loss (solid) and homozygous deletion (dashed). Parameters for every tumour are in samples.tsv.</p><div class='grid'>\n";
  size_t ns = R.S.size(), show = std::min<size_t>(12, ns);
  for (size_t i = 0; i < show; ++i) o << fig::sample_density(R.S[show > 1 ? i * (ns - 1) / (show - 1) : 0]);
  o << "</div>\n<h2>Detector locations across the cohort</h2>\n" << fig::ek_distribution(R);
  o << "<h2>Region sizes</h2><p>Unified breakpoint profile, and the profile after hotspot cleavage in the first null-model round.</p>\n"
    << fig::histogram_pair(usz, csz, "Unified breakpoint profile", "Hotspot-cleaved");
  o << "</main></body></html>\n";
}

}  // namespace rae
