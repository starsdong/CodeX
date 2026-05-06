#include "TAxis.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TBox.h"
#include "TObject.h"
#include "TPad.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

struct ParticleStyle {
  int marker;
  Color_t color;
  double size;
};

struct Series {
  std::vector<double> data_x;
  std::vector<double> data_y;
  std::vector<double> data_ey;
  std::vector<double> gce_x;
  std::vector<double> gce_y;
  std::vector<double> gce_ey_low;
  std::vector<double> gce_ey_high;
  std::vector<double> sce_x;
  std::vector<double> sce_y;
  std::vector<double> sce_ey_low;
  std::vector<double> sce_ey_high;
};

const std::vector<std::string> kParticles = {
    "K+", "K-", "Ks0", "Lambda", "Lambda_bar", "Xi", "Xi_bar", "phi",
};

const std::map<std::string, ParticleStyle> kStyle = {
    {"K+", {20, kBlue + 1, 1.45}},
    {"K-", {21, kAzure + 2, 1.45}},
    {"Ks0", {22, kGreen + 2, 1.55}},
    {"Lambda", {29, kRed + 1, 1.85}},
    {"Lambda_bar", {30, kOrange + 1, 1.65}},
    {"Xi", {33, kMagenta + 1, 1.75}},
    {"Xi_bar", {34, kViolet + 5, 1.65}},
    {"phi", {47, kGray + 2, 1.55}},
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

int findColumn(const std::vector<std::string> &header, const std::vector<std::string> &names)
{
  for (const auto &name : names) {
    for (std::size_t i = 0; i < header.size(); ++i) {
      if (header[i] == name) return static_cast<int>(i);
    }
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

std::string displayLabel(const std::string &particle)
{
  if (particle == "K+") return "K^{+}";
  if (particle == "K-") return "K^{-}";
  if (particle == "Ks0") return "K^{0}_{S}";
  if (particle == "Lambda") return "#Lambda";
  if (particle == "Lambda_bar") return "#bar{#Lambda}";
  if (particle == "Xi") return "#Xi^{-}";
  if (particle == "Xi_bar") return "#bar{#Xi}^{+}";
  if (particle == "phi") return "#phi";
  return particle;
}

double displayEnergy(double energy, const std::string &particle)
{
  if (std::fabs(energy - 3.0) > 1.0e-6) return energy;
  static const std::map<std::string, double> offsets = {
      {"K+", 0.88},
      {"K-", 0.92},
      {"Ks0", 0.96},
      {"Lambda", 1.00},
      {"Lambda_bar", 1.04},
      {"Xi", 1.01},
      {"Xi_bar", 1.12},
      {"phi", 1.04},
  };
  const auto it = offsets.find(particle);
  return energy * (it == offsets.end() ? 1.0 : it->second);
}

void sortPairs(std::vector<double> &x, std::vector<double> &y)
{
  std::vector<std::pair<double, double>> points;
  for (std::size_t i = 0; i < x.size(); ++i) {
    points.emplace_back(x[i], y[i]);
  }
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
  for (std::size_t i = 0; i < x.size(); ++i) {
    points.push_back({x[i], y[i], ey[i]});
  }
  std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
  for (std::size_t i = 0; i < points.size(); ++i) {
    x[i] = points[i].x;
    y[i] = points[i].y;
    ey[i] = points[i].ey;
  }
}

void sortModelBand(std::vector<double> &x, std::vector<double> &y, std::vector<double> &ey_low,
                   std::vector<double> &ey_high)
{
  struct Point {
    double x;
    double y;
    double ey_low;
    double ey_high;
  };
  std::vector<Point> points;
  for (std::size_t i = 0; i < x.size(); ++i) {
    points.push_back({x[i], y[i], ey_low[i], ey_high[i]});
  }
  std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
  for (std::size_t i = 0; i < points.size(); ++i) {
    x[i] = points[i].x;
    y[i] = points[i].y;
    ey_low[i] = points[i].ey_low;
    ey_high[i] = points[i].ey_high;
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

void drawFrameBox(double x1, double x2, double y1, double y2)
{
  TLine *bottom = new TLine(x1, y1, x2, y1);
  TLine *top = new TLine(x1, y2, x2, y2);
  TLine *left = new TLine(x1, y1, x1, y2);
  TLine *right = new TLine(x2, y1, x2, y2);
  for (TLine *line : {bottom, top, left, right}) {
    line->SetLineColor(kBlack);
    line->SetLineWidth(3);
    line->Draw("same");
  }
}

void drawText(double x, double y, const char *text, double size = 0.034)
{
  TLatex *latex = new TLatex(x, y, text);
  latex->SetTextFont(42);
  latex->SetTextSize(size);
  latex->Draw("same");
}

void drawNdcText(double x, double y, const char *text, double size = 0.034)
{
  TLatex *latex = new TLatex(x, y, text);
  latex->SetNDC();
  latex->SetTextFont(42);
  latex->SetTextSize(size);
  latex->Draw("same");
}

void drawLogXLabels(double xmin, double xmax)
{
  const double ticks[] = {3.0, 5.0, 10.0, 30.0, 100.0, 200.0};
  const char *labels[] = {"3", "5", "10", "30", "100", "200"};
  const double left = gPad->GetLeftMargin();
  const double right = 1.0 - gPad->GetRightMargin();
  const double log_min = std::log10(xmin);
  const double log_max = std::log10(xmax);
  const double label_y = std::max(0.035, gPad->GetBottomMargin() - 0.028);

  TLatex text;
  text.SetNDC(kTRUE);
  text.SetTextFont(42);
  text.SetTextSize(0.045);
  text.SetTextAlign(22);
  for (int i = 0; i < 6; ++i) {
    const double frac = (std::log10(ticks[i]) - log_min) / (log_max - log_min);
    const double xndc = left + frac * (right - left);
    text.DrawLatex(xndc, label_y, labels[i]);
  }
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
  graph->SetMarkerSize(style.size * 1.2);
  graph->SetMarkerColor(style.color);
  graph->SetLineColor(style.color);
  graph->SetLineWidth(2);
  return graph;
}

TGraph *makeLineGraph(const Series &series, const ParticleStyle &style)
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

TGraphAsymmErrors *makeBandGraph(const Series &series, const ParticleStyle &style)
{
  if (series.gce_x.empty()) return nullptr;
  double max_err = 0.0;
  for (std::size_t i = 0; i < series.gce_x.size(); ++i) {
    max_err = std::max(max_err, std::max(series.gce_ey_low[i], series.gce_ey_high[i]));
  }
  if (max_err <= 0.0) return nullptr;

  std::vector<double> ex_low(series.gce_x.size(), 0.0);
  std::vector<double> ex_high(series.gce_x.size(), 0.0);
  TGraphAsymmErrors *graph = new TGraphAsymmErrors(
      static_cast<int>(series.gce_x.size()),
      const_cast<double *>(series.gce_x.data()),
      const_cast<double *>(series.gce_y.data()),
      ex_low.data(),
      ex_high.data(),
      const_cast<double *>(series.gce_ey_low.data()),
      const_cast<double *>(series.gce_ey_high.data()));
  graph->SetFillColorAlpha(style.color, 0.18);
  graph->SetLineColor(style.color);
  graph->SetLineWidth(0);
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
  const int i_energy = findColumn(header, {"energy_GeV"});
  const int i_particle = findColumn(header, {"particle"});
  const int i_data = findColumn(header, {"data_yield"});
  const int i_err = findColumn(header, {"data_err"});
  const int i_gce = findColumn(header, {"gc_effective_model_yield", "gc_model_yield", "gce_param_model_yield"});
  const int i_gce_err_low = findColumn(header, {"gc_model_err_low", "gce_param_model_err_low"});
  const int i_gce_err_high = findColumn(header, {"gc_model_err_high", "gce_param_model_err_high"});
  const int i_sce = findColumn(header, {"sce_blended_model_yield", "sce_blended_yield"});
  const int i_sce_err_low = findColumn(header, {"sce_blended_model_err_low", "sce_model_err_low"});
  const int i_sce_err_high = findColumn(header, {"sce_blended_model_err_high", "sce_model_err_high"});

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
      series.data_x.push_back(displayEnergy(energy, particle));
      series.data_y.push_back(value);
      series.data_ey.push_back(err);
    }
    if (energy >= 7.7 && parseDouble(fieldAt(fields, i_gce), value) && value > 0.0) {
      double err_low = 0.0;
      double err_high = 0.0;
      parseDouble(fieldAt(fields, i_gce_err_low), err_low);
      parseDouble(fieldAt(fields, i_gce_err_high), err_high);
      series.gce_x.push_back(energy);
      series.gce_y.push_back(value);
      series.gce_ey_low.push_back(std::max(0.0, err_low));
      series.gce_ey_high.push_back(std::max(0.0, err_high));
    }
    if (std::fabs(energy - 3.0) < 1.0e-6 && parseDouble(fieldAt(fields, i_sce), value) && value > 0.0) {
      double err_low = 0.0;
      double err_high = 0.0;
      parseDouble(fieldAt(fields, i_sce_err_low), err_low);
      parseDouble(fieldAt(fields, i_sce_err_high), err_high);
      series.sce_x.push_back(displayEnergy(energy, particle));
      series.sce_y.push_back(value);
      series.sce_ey_low.push_back(std::max(0.0, err_low));
      series.sce_ey_high.push_back(std::max(0.0, err_high));
    }
  }

  for (auto &entry : rows) {
    sortTriples(entry.second.data_x, entry.second.data_y, entry.second.data_ey);
    sortModelBand(entry.second.gce_x, entry.second.gce_y, entry.second.gce_ey_low, entry.second.gce_ey_high);
    sortModelBand(entry.second.sce_x, entry.second.sce_y, entry.second.sce_ey_low, entry.second.sce_ey_high);
  }
}

}  // namespace

void plot_strange_hadron_yields_root(const char *variant = "exact")
{
  applyFreezeOutStyle();

  TString macro_path = __FILE__;
  if (!macro_path.BeginsWith("/")) {
    macro_path = TString(gSystem->WorkingDirectory()) + "/" + macro_path;
  }
  const TString script_dir = gSystem->DirName(macro_path);
  const TString repo_dir = gSystem->DirName(script_dir);
  const TString variant_text = variant ? variant : "exact";
  const bool use_parametrized = variant_text.Contains("param");
  const TString input_path = use_parametrized
      ? repo_dir + "/data/strange_hadron_yields_parametrized_thermus.csv"
      : repo_dir + "/data/strange_hadron_yields_with_thermus.csv";
  const TString out_pdf = repo_dir + "/data/strange_hadron_yields_root_style.pdf";
  const TString out_png = repo_dir + "/data/strange_hadron_yields_root_style.png";
  const TString param_out_pdf = repo_dir + "/data/strange_hadron_yields_parametrized_thermus_root_style.pdf";
  const TString param_out_png = repo_dir + "/data/strange_hadron_yields_parametrized_thermus_root_style.png";

  std::map<std::string, Series> rows;
  loadRows(input_path, rows);
  if (rows.empty()) {
    std::cerr << "No rows loaded from " << input_path << std::endl;
    return;
  }

  const double xmin = 2.5;
  const double xmax = 300.0;
  const double ymin = 3.0e-3;
  const double ymax = 1.0e3;

  TCanvas *canvas = new TCanvas("cStrangeHadronYields", "", 904, 928);
  canvas->SetLeftMargin(0.13);
  canvas->SetRightMargin(0.04);
  canvas->SetTopMargin(0.03);
  canvas->SetBottomMargin(0.13);
  canvas->SetLogx(1);
  canvas->SetLogy(1);
  canvas->Draw();

  TH1D *frame = new TH1D("frameStrangeHadronYields", "", 1, xmin, xmax);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetXaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  frame->GetXaxis()->SetMoreLogLabels(kTRUE);
  frame->GetXaxis()->SetNoExponent(kTRUE);
  frame->GetXaxis()->SetLabelSize(0.0);
  frame->GetXaxis()->SetTickLength(0.03);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetXaxis()->SetTitleSize(0.055);
  frame->GetYaxis()->CenterTitle();
  frame->GetYaxis()->SetTitle("dN/dy");
  frame->GetYaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleSize(0.057);
  frame->GetYaxis()->SetLabelSize(0.043);
  frame->Draw();

  std::vector<TGraph *> model_graphs;
  std::vector<TGraphAsymmErrors *> model_bands;
  std::vector<TGraphErrors *> data_graphs;

  for (const std::string &particle : kParticles) {
    const auto row_it = rows.find(particle);
    if (row_it == rows.end()) continue;
    const ParticleStyle style = kStyle.at(particle);
    TGraphAsymmErrors *band = makeBandGraph(row_it->second, style);
    if (band) {
      model_bands.push_back(band);
      band->Draw("3 SAME");
    }
  }

  for (const std::string &particle : kParticles) {
    const auto row_it = rows.find(particle);
    if (row_it == rows.end()) continue;
    const ParticleStyle style = kStyle.at(particle);
    TGraph *model = makeLineGraph(row_it->second, style);
    if (model) {
      model_graphs.push_back(model);
      model->Draw("C SAME");
    }
  }

  for (const std::string &particle : kParticles) {
    const auto row_it = rows.find(particle);
    if (row_it == rows.end()) continue;
    const ParticleStyle style = kStyle.at(particle);
    for (std::size_t i = 0; i < row_it->second.sce_x.size(); ++i) {
      const double x = row_it->second.sce_x[i];
      const double y = row_it->second.sce_y[i];
      const double ey_low = row_it->second.sce_ey_low[i];
      const double ey_high = row_it->second.sce_ey_high[i];
      const double y_low = std::max(ymin * 1.03, y - ey_low);
      const double y_high = std::max(y_low * 1.01, y + ey_high);
      const double half_width = 0.12;
      TBox *band = new TBox(x - half_width, y_low, x + half_width, y_high);
      band->SetFillColorAlpha(style.color, 0.18);
      band->SetLineColor(style.color);
      band->SetLineWidth(1);
      band->Draw("same");
      TLine *line = new TLine(x - half_width, y, x + half_width, y);
      line->SetLineColor(style.color);
      line->SetLineWidth(4);
      line->Draw("same");
    }
  }

  for (const std::string &particle : kParticles) {
    const auto row_it = rows.find(particle);
    if (row_it == rows.end()) continue;
    const ParticleStyle style = kStyle.at(particle);
    TGraphErrors *data = makeDataGraph(row_it->second, style);
    if (data) {
      data_graphs.push_back(data);
      data->Draw("P SAME");
    }
  }

  drawNdcText(0.18, 0.91, "STAR Au+Au 0-5%", 0.043);

  TLegend *leg_model = new TLegend(0.62, 0.83, 0.96, 0.95);
  leg_model->SetFillStyle(4000);
  leg_model->SetBorderSize(0);
  leg_model->SetTextSize(0.025);

  TGraph *dummy_gce = new TGraph();
  dummy_gce->SetLineColor(kBlack);
  dummy_gce->SetLineWidth(2);
  leg_model->AddEntry(dummy_gce, "THERMUS", "l");

  leg_model->AddEntry((TObject *)0, "GCE #pm1#sigma, 7.7-200 GeV", "");
  leg_model->AddEntry((TObject *)0, "SCE #pm1#sigma, 3 GeV", "");
  leg_model->Draw();

//  drawNdcText(0.65, 0.36, "Data", 0.022);
  TLegend *leg_data_left = new TLegend(0.68, 0.17, 0.79, 0.37);
  TLegend *leg_data_right = new TLegend(0.80, 0.17, 0.94, 0.37);
  for (TLegend *legend : {leg_data_left, leg_data_right}) {
    legend->SetFillStyle(4000);
    legend->SetBorderSize(0);
    legend->SetTextSize(0.030);
  }
  const std::vector<std::string> left_data_order = {"K+", "Ks0", "K-", "phi"};
  const std::vector<std::string> right_data_order = {"Lambda", "Xi", "Lambda_bar", "Xi_bar"};
  for (const std::string &particle : left_data_order) {
    const ParticleStyle style = kStyle.at(particle);
    TGraphErrors *dummy = new TGraphErrors();
    dummy->SetMarkerStyle(style.marker);
    dummy->SetMarkerSize(style.size * 0.75);
    dummy->SetMarkerColor(style.color);
    dummy->SetLineColor(style.color);
    leg_data_left->AddEntry(dummy, displayLabel(particle).c_str(), "p");
  }
  for (const std::string &particle : right_data_order) {
    const ParticleStyle style = kStyle.at(particle);
    TGraphErrors *dummy = new TGraphErrors();
    dummy->SetMarkerStyle(style.marker);
    dummy->SetMarkerSize(style.size * 0.75);
    dummy->SetMarkerColor(style.color);
    dummy->SetLineColor(style.color);
    leg_data_right->AddEntry(dummy, displayLabel(particle).c_str(), "p");
  }
  leg_data_left->Draw();
  leg_data_right->Draw();

  drawFrameBox(xmin, xmax, ymin, ymax);
  drawLogXLabels(xmin, xmax);

  canvas->Update();
  canvas->SaveAs(out_pdf);
  canvas->SaveAs(out_png);
  if (use_parametrized) {
    canvas->SaveAs(param_out_pdf);
    canvas->SaveAs(param_out_png);
  }

  std::cout << "wrote " << out_pdf << std::endl;
  std::cout << "wrote " << out_png << std::endl;
  if (use_parametrized) {
    std::cout << "wrote " << param_out_pdf << std::endl;
    std::cout << "wrote " << param_out_png << std::endl;
  }
}
