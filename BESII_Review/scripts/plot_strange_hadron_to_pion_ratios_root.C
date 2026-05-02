#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

struct ParticleStyle {
  const char *title;
  const char *file_tag;
  int marker;
  Color_t color;
  double marker_size;
};

struct Series {
  std::vector<double> data_x;
  std::vector<double> data_y;
  std::vector<double> data_ey;
  std::vector<double> gce_x;
  std::vector<double> gce_y;
  std::vector<double> sce_x;
  std::vector<double> sce_y;
};

const std::vector<std::string> kParticles = {
    "K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi",
};

const std::map<std::string, ParticleStyle> kStyle = {
    {"K+", {"K^{+}", "Kplus", 20, kBlue + 1, 1.25}},
    {"K-", {"K^{-}", "Kminus", 21, kAzure + 2, 1.25}},
    {"Ks0", {"K^{0}_{S}", "Ks0", 22, kGreen + 2, 1.35}},
    {"Lambda", {"#Lambda", "Lambda", 29, kRed + 1, 1.55}},
    {"Lambda_bar", {"#bar{#Lambda}", "Lambda_bar", 30, kOrange + 1, 1.45}},
    {"Xi", {"#Xi^{-}", "Xi", 33, kMagenta + 1, 1.50}},
    {"Xi_bar", {"#bar{#Xi}^{+}", "Xi_bar", 34, kViolet + 5, 1.45}},
    {"phi", {"#phi", "phi", 47, kGray + 2, 1.35}},
};

std::vector<std::string> splitCsvLine(const std::string &line)
{
  std::vector<std::string> fields;
  std::string field;
  bool in_quotes = false;
  for (std::size_t i = 0; i < line.size(); ++i) {
    const char c = line[i];
    if (c == '"') {
      if (in_quotes && i + 1 < line.size() && line[i + 1] == '"') {
        field.push_back('"');
        ++i;
      } else {
        in_quotes = !in_quotes;
      }
    } else if (c == ',' && !in_quotes) {
      fields.push_back(field);
      field.clear();
    } else {
      field.push_back(c);
    }
  }
  fields.push_back(field);
  return fields;
}

bool parseDouble(const std::string &text, double &value)
{
  if (text.empty()) return false;
  char *end = nullptr;
  value = std::strtod(text.c_str(), &end);
  return end != text.c_str();
}

int findColumn(const std::vector<std::string> &header, const char *name)
{
  for (std::size_t i = 0; i < header.size(); ++i) {
    if (header[i] == name) return static_cast<int>(i);
  }
  return -1;
}

std::string fieldAt(const std::vector<std::string> &fields, int index)
{
  if (index < 0 || static_cast<std::size_t>(index) >= fields.size()) return "";
  return fields[index];
}

bool knownParticle(const std::string &particle)
{
  return std::find(kParticles.begin(), kParticles.end(), particle) != kParticles.end();
}

void sortPairs(std::vector<double> &x, std::vector<double> &y)
{
  std::vector<std::pair<double, double>> points;
  for (std::size_t i = 0; i < x.size(); ++i) points.emplace_back(x[i], y[i]);
  std::sort(points.begin(), points.end());
  for (std::size_t i = 0; i < points.size(); ++i) {
    x[i] = points[i].first;
    y[i] = points[i].second;
  }
}

void sortTriples(std::vector<double> &x, std::vector<double> &y, std::vector<double> &ey)
{
  struct Point {
    double x;
    double y;
    double ey;
  };
  std::vector<Point> points;
  for (std::size_t i = 0; i < x.size(); ++i) points.push_back({x[i], y[i], ey[i]});
  std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
  for (std::size_t i = 0; i < points.size(); ++i) {
    x[i] = points[i].x;
    y[i] = points[i].y;
    ey[i] = points[i].ey;
  }
}

void applyFreezeOutStyle()
{
  gStyle->SetOptFit(0);
  gStyle->SetOptStat(0);
  gStyle->SetEndErrorSize(0.01);
  gStyle->SetTickLength(0.04, "X");
  gStyle->SetTickLength(0.04, "Y");
  gStyle->SetGridWidth(1);
  gStyle->SetGridStyle(2);
  gStyle->SetGridColor(kGray + 1);

  gStyle->SetFillColor(10);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadBorderSize(2);
  gStyle->SetFrameFillColor(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetPadLeftMargin(0.18);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadTopMargin(0.03);
  gStyle->SetPadBottomMargin(0.18);

  gStyle->SetNdivisions(108, "X");
  gStyle->SetNdivisions(108, "Y");
  gStyle->SetLabelOffset(0.005, "X");
  gStyle->SetLabelOffset(0.01, "Y");
  gStyle->SetLabelSize(0.05, "X");
  gStyle->SetLabelSize(0.05, "Y");
  gStyle->SetLabelFont(42, "X");
  gStyle->SetLabelFont(42, "Y");
  gStyle->SetTitleOffset(1.0, "X");
  gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetTitleSize(0.07, "X");
  gStyle->SetTitleSize(0.07, "Y");
  gStyle->SetTitleFont(42, "X");
  gStyle->SetTitleFont(42, "Y");

  gStyle->SetMarkerSize(1.8);
  gStyle->SetMarkerStyle(20);
  gStyle->SetLegendFillColor(10);
  gROOT->ForceStyle();
}

void drawFrameBox(double x1, double x2, double y1, double y2, int width = 2)
{
  TLine *bottom = new TLine(x1, y1, x2, y1);
  TLine *top = new TLine(x1, y2, x2, y2);
  TLine *left = new TLine(x1, y1, x1, y2);
  TLine *right = new TLine(x2, y1, x2, y2);
  for (TLine *line : {bottom, top, left, right}) {
    line->SetLineColor(kBlack);
    line->SetLineWidth(width);
    line->Draw("same");
  }
}

void drawText(double x, double y, const char *text, double size = 0.05)
{
  TLatex *latex = new TLatex(x, y, text);
  latex->SetTextFont(42);
  latex->SetTextSize(size);
  latex->Draw("same");
}

double maxValue(const Series &series)
{
  double y_max = 0.0;
  for (double y : series.data_y) y_max = std::max(y_max, y);
  for (std::size_t i = 0; i < series.data_y.size(); ++i) {
    if (i < series.data_ey.size()) y_max = std::max(y_max, series.data_y[i] + series.data_ey[i]);
  }
  for (double y : series.gce_y) y_max = std::max(y_max, y);
  for (double y : series.sce_y) y_max = std::max(y_max, y);
  return y_max;
}

TGraphErrors *makeDataGraph(const Series &series, const ParticleStyle &style)
{
  if (series.data_x.empty()) return nullptr;
  std::vector<double> ex(series.data_x.size(), 0.0);
  TGraphErrors *graph = new TGraphErrors(
      static_cast<int>(series.data_x.size()),
      const_cast<double *>(series.data_x.data()),
      const_cast<double *>(series.data_y.data()),
      ex.data(),
      const_cast<double *>(series.data_ey.data()));
  graph->SetMarkerStyle(style.marker);
  graph->SetMarkerSize(style.marker_size);
  graph->SetMarkerColor(style.color);
  graph->SetLineColor(style.color);
  graph->SetLineWidth(2);
  return graph;
}

TGraph *makeGceGraph(const Series &series, const ParticleStyle &style)
{
  if (series.gce_x.empty()) return nullptr;
  TGraph *graph = new TGraph(
      static_cast<int>(series.gce_x.size()),
      const_cast<double *>(series.gce_x.data()),
      const_cast<double *>(series.gce_y.data()));
  graph->SetLineColor(style.color);
  graph->SetLineWidth(2);
  graph->SetLineStyle(1);
  return graph;
}

void loadRows(const TString &input_path, std::map<std::string, Series> &rows)
{
  std::ifstream input(input_path.Data());
  if (!input) {
    std::cerr << "Cannot open " << input_path << std::endl;
    return;
  }

  std::string line;
  if (!std::getline(input, line)) return;
  const std::vector<std::string> header = splitCsvLine(line);
  const int i_energy = findColumn(header, "energy_GeV");
  const int i_particle = findColumn(header, "particle");
  const int i_data = findColumn(header, "data_ratio");
  const int i_err = findColumn(header, "data_ratio_err");
  const int i_gce = findColumn(header, "gc_ratio");
  const int i_sce = findColumn(header, "sce_ratio");

  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> fields = splitCsvLine(line);
    double energy = 0.0;
    if (!parseDouble(fieldAt(fields, i_energy), energy)) continue;

    const std::string particle = fieldAt(fields, i_particle);
    if (!knownParticle(particle)) continue;

    Series &series = rows[particle];
    double value = 0.0;
    if (parseDouble(fieldAt(fields, i_data), value) && value > 0.0) {
      double err = 0.0;
      parseDouble(fieldAt(fields, i_err), err);
      series.data_x.push_back(energy);
      series.data_y.push_back(value);
      series.data_ey.push_back(err);
    }
    if (parseDouble(fieldAt(fields, i_gce), value) && value > 0.0) {
      series.gce_x.push_back(energy);
      series.gce_y.push_back(value);
    }
    if (parseDouble(fieldAt(fields, i_sce), value) && value > 0.0) {
      series.sce_x.push_back(energy);
      series.sce_y.push_back(value);
    }
  }

  for (auto &entry : rows) {
    sortTriples(entry.second.data_x, entry.second.data_y, entry.second.data_ey);
    sortPairs(entry.second.gce_x, entry.second.gce_y);
    sortPairs(entry.second.sce_x, entry.second.sce_y);
  }
}

void drawRatioPanel(const std::string &particle, const Series &series, bool show_legend, bool large_labels)
{
  const ParticleStyle style = kStyle.at(particle);
  const double xmin = 2.5;
  const double xmax = 250.0;
  const double ymin = 0.0;
  double ymax = maxValue(series) * 1.35;
  if (ymax <= 0.0) ymax = 0.05;
  if (ymax < 0.012) ymax = 0.012;

  TH1D *frame = new TH1D(
      Form("frame_ratio_%s_%p", style.file_tag, gPad),
      "",
      1,
      xmin,
      xmax);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  frame->GetYaxis()->SetTitle(large_labels ? "Ratio to #frac{1}{2}(#pi^{+}+#pi^{-})" : "Ratio");

  if (large_labels) {
    frame->GetXaxis()->SetLabelSize(0.045);
    frame->GetYaxis()->SetLabelSize(0.045);
    frame->GetXaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetXaxis()->SetTitleOffset(1.1);
    frame->GetYaxis()->SetTitleOffset(1.25);
  } else {
    frame->GetXaxis()->SetLabelSize(0.065);
    frame->GetYaxis()->SetLabelSize(0.065);
    frame->GetXaxis()->SetTitleSize(0.075);
    frame->GetYaxis()->SetTitleSize(0.075);
    frame->GetXaxis()->SetTitleOffset(0.95);
    frame->GetYaxis()->SetTitleOffset(0.90);
  }
  frame->Draw();

  TGraph *gce = makeGceGraph(series, style);
  if (gce) gce->Draw("C SAME");

  for (std::size_t i = 0; i < series.sce_x.size(); ++i) {
    TLine *line = new TLine(series.sce_x[i] * 0.93, series.sce_y[i], series.sce_x[i] * 1.07, series.sce_y[i]);
    line->SetLineColor(style.color);
    line->SetLineWidth(3);
    line->Draw("same");
  }

  TGraphErrors *data = makeDataGraph(series, style);
  if (data) data->Draw("P SAME");

  drawText(3.1, ymax * 0.84, style.title, large_labels ? 0.050 : 0.085);

  if (show_legend) {
    TGraphErrors *dummy_data = new TGraphErrors();
    dummy_data->SetMarkerStyle(style.marker);
    dummy_data->SetMarkerSize(style.marker_size);
    dummy_data->SetMarkerColor(style.color);
    dummy_data->SetLineColor(style.color);

    TGraph *dummy_gce = new TGraph();
    dummy_gce->SetLineColor(style.color);
    dummy_gce->SetLineWidth(2);

    TLine *dummy_sce = new TLine();
    dummy_sce->SetLineColor(style.color);
    dummy_sce->SetLineWidth(3);

    TLegend *legend = new TLegend(large_labels ? 0.47 : 0.42, large_labels ? 0.67 : 0.14, 0.92, large_labels ? 0.90 : 0.38);
    legend->SetFillStyle(4000);
    legend->SetBorderSize(0);
    legend->SetTextSize(large_labels ? 0.034 : 0.060);
    legend->AddEntry(dummy_data, "Measured ratio", "p");
    legend->AddEntry(dummy_gce, "THERMUS GCE", "l");
    if (!series.sce_x.empty()) legend->AddEntry(dummy_sce, "THERMUS SCE", "l");
    legend->Draw();
  }

  drawFrameBox(xmin, xmax, ymin, ymax, large_labels ? 3 : 2);
}

void drawSingleFigure(const TString &out_dir, const std::string &particle, const Series &series)
{
  const ParticleStyle style = kStyle.at(particle);
  TCanvas *canvas = new TCanvas(Form("cRatio_%s", style.file_tag), "", 900, 700);
  canvas->SetLeftMargin(0.13);
  canvas->SetBottomMargin(0.13);
  canvas->SetRightMargin(0.03);
  canvas->SetTopMargin(0.03);
  canvas->SetLogx(1);
  canvas->Draw();
  drawRatioPanel(particle, series, true, true);
  canvas->Update();
  canvas->SaveAs(out_dir + "/" + style.file_tag + "_to_pion_ratio_root_style.pdf");
  canvas->SaveAs(out_dir + "/" + style.file_tag + "_to_pion_ratio_root_style.png");
  delete canvas;
}

}  // namespace

void plot_strange_hadron_to_pion_ratios_root()
{
  applyFreezeOutStyle();

  TString macro_path = __FILE__;
  if (!macro_path.BeginsWith("/")) {
    macro_path = TString(gSystem->WorkingDirectory()) + "/" + macro_path;
  }
  const TString script_dir = gSystem->DirName(macro_path);
  const TString repo_dir = gSystem->DirName(script_dir);
  const TString input_path = repo_dir + "/data/particle_to_pi_ratios_with_thermus.csv";
  const TString out_pdf = repo_dir + "/data/strange_hadron_to_pion_ratios_root_panels.pdf";
  const TString out_png = repo_dir + "/data/strange_hadron_to_pion_ratios_root_panels.png";
  const TString out_dir = repo_dir + "/data/strange_hadron_to_pion_ratio_root_figures";
  gSystem->mkdir(out_dir, kTRUE);

  std::map<std::string, Series> rows;
  loadRows(input_path, rows);
  if (rows.empty()) {
    std::cerr << "No rows loaded from " << input_path << std::endl;
    return;
  }

  TCanvas *canvas = new TCanvas("cStrangeHadronToPionRatios", "", 1100, 1400);
  canvas->Divide(2, 4, 0.001, 0.001);
  for (std::size_t i = 0; i < kParticles.size(); ++i) {
    canvas->cd(static_cast<int>(i + 1));
    gPad->SetLeftMargin(0.17);
    gPad->SetRightMargin(0.04);
    gPad->SetTopMargin(0.06);
    gPad->SetBottomMargin(0.17);
    gPad->SetLogx(1);

    const std::string &particle = kParticles[i];
    const auto it = rows.find(particle);
    if (it == rows.end()) continue;
    drawRatioPanel(particle, it->second, i == 0, false);
  }
  canvas->Update();
  canvas->SaveAs(out_pdf);
  canvas->SaveAs(out_png);

  for (const std::string &particle : kParticles) {
    const auto it = rows.find(particle);
    if (it == rows.end()) continue;
    drawSingleFigure(out_dir, particle, it->second);
  }

  std::cout << "wrote " << out_pdf << std::endl;
  std::cout << "wrote " << out_png << std::endl;
  std::cout << "wrote individual figures in " << out_dir << std::endl;
}
