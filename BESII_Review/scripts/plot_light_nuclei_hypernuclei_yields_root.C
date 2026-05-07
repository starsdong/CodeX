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
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr double kH3LambdaBR = 0.23;
constexpr double kAliceH3LambdaBR = 0.25;

constexpr double kRatioXMax = 300.0;
constexpr double kRatioYMin = 8.0e-6;
constexpr double kTpVisibilityScale = 5.0;

struct CsvRow {
  std::map<std::string, std::string> fields;
};

struct Series {
  std::string label;
  int marker = 20;
  Color_t color = kBlack;
  double marker_size = 1.4;
  int line_width = 2;
  bool draw_line = false;
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> ey;
};

struct ModelCurve {
  std::string label;
  int color = kBlack;
  int style = 1;
  int width = 2;
  bool show_in_legend = true;
  std::vector<double> x;
  std::vector<double> y;
};

struct YieldPoint {
  double x = 0.0;
  double y = 0.0;
  double ey = 0.0;
  bool ok = false;
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

std::string field(const CsvRow &row, const std::string &name)
{
  const auto it = row.fields.find(name);
  if (it == row.fields.end()) return "";
  return it->second;
}

bool hasObservable(const CsvRow &row, const std::set<std::string> &observables)
{
  return observables.count(field(row, "observable")) > 0;
}

double rowValue(const CsvRow &row, const std::string &name, double fallback = 0.0)
{
  double value = fallback;
  parseDouble(field(row, name), value);
  return value;
}

double rowError(const CsvRow &row)
{
  double total = 0.0;
  if (parseDouble(field(row, "total_err"), total)) return total;
  const double stat = rowValue(row, "stat_err", 0.0);
  const double sys = rowValue(row, "sys_err", 0.0);
  return std::sqrt(stat * stat + sys * sys);
}

std::string normalizedObservable(const CsvRow &row)
{
  const std::string particle = field(row, "particle");
  const std::string obs = field(row, "observable");
  if (particle == "H3Lambda" && obs == "BR*dN/dy" && field(row, "source_id") == "ins1946124") return "dN/dy";
  if (particle == "H3Lambda" && obs == "BR*dN/dy" && field(row, "source_id") == "ALICE-PLB754-2016") return "dN/dy";
  if (particle == "H3Lambda" && obs == "dN/dy (BR=0.23 assumed)") return "dN/dy";
  return obs;
}

bool isXinH3LambdaDndy(const CsvRow &row)
{
  return field(row, "particle") == "H3Lambda" &&
         field(row, "observable") == "dN/dy (BR=0.23 assumed)" &&
         field(row, "source_table") == "data/data_for_Xin/BES_II_H3L_dNdy.txt";
}

bool isStar2110H3LambdaDndy(const CsvRow &row)
{
  return field(row, "particle") == "H3Lambda" &&
         field(row, "observable") == "BR*dN/dy" &&
         field(row, "source_id") == "ins1946124" &&
         field(row, "source_table") == "Table 6";
}

bool isAlice2760H3LambdaDndy(const CsvRow &row)
{
  return field(row, "particle") == "H3Lambda" &&
         field(row, "observable") == "BR*dN/dy" &&
         field(row, "source_id") == "ALICE-PLB754-2016";
}

double rapidityYieldScale(const CsvRow &row)
{
  if (isStar2110H3LambdaDndy(row)) return 1.0 / kH3LambdaBR;
  if (isAlice2760H3LambdaDndy(row)) return 1.0 / kAliceH3LambdaBR;
  return 1.0;
}

bool isDisplayedRapidityParticle(const CsvRow &row)
{
  const std::string particle = field(row, "particle");
  return particle == "p" || particle == "Lambda" ||
         particle == "d" || particle == "t" || particle == "He3" || particle == "H3Lambda";
}

double multiplicityHigh(const CsvRow &row)
{
  const std::string note = field(row, "kinematic_note");
  const std::size_t open = note.find('[');
  const std::size_t comma = note.find(',', open == std::string::npos ? 0 : open);
  const std::size_t close = note.find(']', comma == std::string::npos ? 0 : comma);
  if (open == std::string::npos || comma == std::string::npos || close == std::string::npos) {
    return -1.0;
  }
  double value = -1.0;
  parseDouble(note.substr(comma + 1, close - comma - 1), value);
  return value;
}

bool parseCentralityRange(const std::string &cent, double &low, double &high)
{
  const std::size_t dash = cent.find('-');
  if (dash == std::string::npos) return false;
  return parseDouble(cent.substr(0, dash), low) && parseDouble(cent.substr(dash + 1), high);
}

bool isCentralSummaryRow(const CsvRow &row)
{
  const std::string cent = field(row, "centrality");
  if (cent == "max Npart row (most central)") return true;

  double low = -1.0;
  double high = -1.0;
  if (parseCentralityRange(cent, low, high)) {
    return std::fabs(low) < 1.0e-6 && high <= 20.0;
  }

  return false;
}

std::string multiplicityGroupKey(const CsvRow &row)
{
  return field(row, "experiment") + "|" + field(row, "system") + "|" + field(row, "sqrt_sNN_GeV") + "|" +
         field(row, "particle") + "|" + field(row, "matter") + "|" + field(row, "observable");
}

std::vector<CsvRow> readCsv(const TString &path)
{
  std::vector<CsvRow> rows;
  std::ifstream input(path.Data());
  if (!input) {
    std::cerr << "Cannot open " << path << std::endl;
    return rows;
  }

  std::string line;
  if (!std::getline(input, line)) return rows;
  const std::vector<std::string> header = splitCsvLine(line);
  while (std::getline(input, line)) {
    if (line.empty()) continue;
    const std::vector<std::string> values = splitCsvLine(line);
    CsvRow row;
    for (std::size_t i = 0; i < header.size(); ++i) {
      row.fields[header[i]] = i < values.size() ? values[i] : "";
    }
    rows.push_back(row);
  }
  return rows;
}

std::vector<CsvRow> selectSummaryRows(const std::vector<CsvRow> &rows)
{
  const std::set<std::string> rapidity_observables = {
      "dN/dy", "BR*dN/dy", "dN/dy (BR=0.23 assumed)", "dN/dy (BR=0.25)",
  };

  std::vector<CsvRow> selected;
  std::map<std::string, std::vector<CsvRow>> multiplicity_groups;

  for (const CsvRow &row : rows) {
    const std::string obs = field(row, "observable");
    if (obs == "4pi_yield") {
      selected.push_back(row);
      continue;
    }
    if (!rapidity_observables.count(obs)) continue;

    if (field(row, "centrality").empty() && multiplicityHigh(row) >= 0.0) {
      multiplicity_groups[multiplicityGroupKey(row)].push_back(row);
      continue;
    }
    if (isCentralSummaryRow(row)) selected.push_back(row);
  }

  for (const auto &entry : multiplicity_groups) {
    const auto best = std::max_element(
        entry.second.begin(),
        entry.second.end(),
        [](const CsvRow &a, const CsvRow &b) { return multiplicityHigh(a) < multiplicityHigh(b); });
    if (best != entry.second.end()) selected.push_back(*best);
  }

  return selected;
}

std::string seriesKey(const CsvRow &row)
{
  return field(row, "particle") + "|" + field(row, "matter") + "|" + normalizedObservable(row);
}

std::string displayLabel(const CsvRow &row)
{
  const std::string particle = field(row, "particle");
  const std::string obs = normalizedObservable(row);

  if (particle == "p") return "p";
  if (particle == "Lambda") return "#Lambda";
  if (particle == "d") return "d";
  if (particle == "anti-d") return "anti-d";
  if (particle == "t" && obs == "4pi_yield") return "t 4#pi";
  if (particle == "t") return "t";
  if (particle == "anti-t") return "anti-t";
  if (particle == "He3") return "{}^{3}He";
  if (particle == "anti-He3") return "anti-{}^{3}He";
  if (particle == "He4") return "{}^{4}He";
  if (particle == "anti-He4") return "anti-{}^{4}He";
  if (particle == "H4Lambda") return "{}_{#Lambda}^{4}H BR#timesdN/dy";
  if (particle == "H3Lambda" && obs == "dN/dy") return "{}_{#Lambda}^{3}H";
  if (particle == "(H3Lambda+anti-H3Lambda)/2") return "({}_{#Lambda}^{3}H+anti)/2, BR=0.25";
  return particle;
}

void applyStyleForRow(const CsvRow &row, Series &series)
{
  const std::string particle = field(row, "particle");
  const std::string obs = normalizedObservable(row);
  series.label = displayLabel(row);
  series.marker_size = 1.35;
  series.draw_line = false;

  if (particle == "p") {
    series.color = kRed + 1;
    series.marker = 34;
    series.marker_size = 1.45;
  } else if (particle == "Lambda") {
    series.color = kAzure + 2;
    series.marker = 33;
    series.marker_size = 1.45;
  } else if (particle == "d") {
    series.color = kBlue + 1;
    series.marker = 20;
  } else if (particle == "anti-d") {
    series.color = kBlue + 1;
    series.marker = 24;
  } else if (particle == "t") {
    series.color = kOrange + 7;
    series.marker = 21;
    series.marker_size = 1.55;
    series.line_width = 3;
  } else if (particle == "anti-t") {
    series.color = kOrange + 7;
    series.marker = 25;
  } else if (particle == "He3") {
    series.color = kGreen + 2;
    series.marker = 22;
  } else if (particle == "anti-He3") {
    series.color = kGreen + 2;
    series.marker = 26;
  } else if (particle == "He4") {
    series.color = kRed + 1;
    series.marker = 33;
    series.marker_size = 1.55;
  } else if (particle == "anti-He4") {
    series.color = kRed + 1;
    series.marker = 27;
    series.marker_size = 1.55;
  } else if (particle == "H4Lambda") {
    series.color = kOrange + 3;
    series.marker = 34;
    series.marker_size = 1.45;
  } else if (particle == "H3Lambda") {
    series.color = kViolet + 5;
    series.marker = 23;
    series.marker_size = 1.45;
  } else if (particle == "(H3Lambda+anti-H3Lambda)/2") {
    series.color = kViolet + 5;
    series.marker = 30;
    series.marker_size = 1.45;
  }
}

void sortSeries(Series &series)
{
  struct Point {
    double x;
    double y;
    double ey;
  };
  std::vector<Point> points;
  for (std::size_t i = 0; i < series.x.size(); ++i) {
    points.push_back({series.x[i], series.y[i], series.ey[i]});
  }
  std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
  for (std::size_t i = 0; i < points.size(); ++i) {
    series.x[i] = points[i].x;
    series.y[i] = points[i].y;
    series.ey[i] = points[i].ey;
  }
}

std::map<std::string, Series> buildSeries(const std::vector<CsvRow> &rows, bool four_pi)
{
  const std::set<std::string> rapidity_observables = {
      "dN/dy", "BR*dN/dy", "dN/dy (BR=0.23 assumed)", "dN/dy (BR=0.25)",
  };

  std::map<std::string, Series> series_map;
  for (const CsvRow &row : rows) {
    if (!four_pi && !isDisplayedRapidityParticle(row)) continue;

    const std::string obs = field(row, "observable");
    if (four_pi) {
      if (obs != "4pi_yield") continue;
    } else if (!rapidity_observables.count(obs)) {
      continue;
    }
    if (!four_pi && field(row, "particle") == "H3Lambda" &&
        !isXinH3LambdaDndy(row) && !isStar2110H3LambdaDndy(row) && !isAlice2760H3LambdaDndy(row)) {
      continue;
    }

    double x = 0.0;
    double y = 0.0;
    if (!parseDouble(field(row, "sqrt_sNN_GeV"), x) || !parseDouble(field(row, "value"), y) || y <= 0.0) {
      continue;
    }
    if (x > 6000.0) continue;

    const std::string key = seriesKey(row);
    Series &series = series_map[key];
    if (series.label.empty()) applyStyleForRow(row, series);
    const double scale = rapidityYieldScale(row);
    series.x.push_back(x);
    series.y.push_back(scale * y);
    series.ey.push_back(scale * rowError(row));
  }

  for (auto &entry : series_map) sortSeries(entry.second);
  return series_map;
}

bool isOverlayYieldParticle(const std::string &particle)
{
  return particle == "Lambda";
}

void appendOverlayYieldData(const TString &path, std::map<std::string, Series> &series_map)
{
  const std::vector<CsvRow> rows = readCsv(path);
  for (const CsvRow &row : rows) {
    const std::string particle = field(row, "particle");
    if (!isOverlayYieldParticle(particle)) continue;

    double x = 0.0;
    double y = 0.0;
    double ey = 0.0;
    if (!parseDouble(field(row, "energy_GeV"), x) || !parseDouble(field(row, "data_yield"), y)) continue;
    parseDouble(field(row, "data_err"), ey);
    if (x <= 0.0 || y <= 0.0 || x > 6000.0) continue;

    CsvRow style_row;
    style_row.fields["particle"] = particle;
    style_row.fields["matter"] = "matter";
    style_row.fields["observable"] = "dN/dy";

    const std::string key = seriesKey(style_row);
    Series &series = series_map[key];
    if (series.label.empty()) applyStyleForRow(style_row, series);
    series.x.push_back(x);
    series.y.push_back(y);
    series.ey.push_back(ey);
  }

  for (auto &entry : series_map) sortSeries(entry.second);
}

std::string energyKey(double energy)
{
  std::ostringstream out;
  out << std::fixed << std::setprecision(3) << energy;
  return out.str();
}

void ratioStyle(const std::string &ratio, Series &series)
{
  series.label = ratio;
  series.marker_size = 1.35;
  series.draw_line = false;

  if (ratio == "d/p") {
    series.label = "d/p (prompt p)";
    series.color = kBlue + 1;
    series.marker = 20;
  } else if (ratio == "t/p") {
    series.label = "5#times t/p";
    series.color = kOrange + 7;
    series.marker = 21;
    series.marker_size = 1.55;
  } else if (ratio == "{}^{3}He/p") {
    series.color = kGreen + 2;
    series.marker = 22;
  } else if (ratio == "{}_{#Lambda}^{3}H/#Lambda") {
    series.color = kViolet + 5;
    series.marker = 23;
    series.marker_size = 1.45;
  }
}

int colorForRatio(const std::string &ratio)
{
  if (ratio == "d/p") return kBlue + 1;
  if (ratio == "t/p") return kOrange + 7;
  if (ratio == "{}^{3}He/p") return kGreen + 2;
  if (ratio == "{}_{#Lambda}^{3}H/#Lambda") return kViolet + 5;
  return kGray + 2;
}

std::string ratioLabelForNumerator(const std::string &particle)
{
  if (particle == "d") return "d/p";
  if (particle == "t") return "t/p";
  if (particle == "He3") return "{}^{3}He/p";
  if (particle == "H3Lambda") return "{}_{#Lambda}^{3}H/#Lambda";
  return particle;
}

double ratioDisplayScale(const std::string &ratio)
{
  return ratio == "t/p" ? kTpVisibilityScale : 1.0;
}

void addRatioPoint(
    std::map<std::string, Series> &series_map,
    const std::string &ratio,
    const YieldPoint &num,
    const YieldPoint &den)
{
  if (!num.ok || !den.ok || num.y <= 0.0 || den.y <= 0.0) return;
  const double scale = ratioDisplayScale(ratio);
  const double value = scale * num.y / den.y;
  const double rel_num = num.ey > 0.0 ? num.ey / num.y : 0.0;
  const double rel_den = den.ey > 0.0 ? den.ey / den.y : 0.0;
  const double err = value * std::sqrt(rel_num * rel_num + rel_den * rel_den);

  Series &series = series_map[ratio];
  if (series.label.empty()) ratioStyle(ratio, series);
  series.x.push_back(num.x);
  series.y.push_back(value);
  series.ey.push_back(err);
}

std::map<std::string, std::map<std::string, YieldPoint>> readLightYieldMaps(const std::vector<CsvRow> &summary)
{
  std::map<std::string, std::map<std::string, YieldPoint>> values;
  for (const CsvRow &row : summary) {
    const std::string particle = field(row, "particle");
    if (particle != "p" && particle != "d" && particle != "t" && particle != "He3") continue;
    if (normalizedObservable(row) != "dN/dy") continue;

    double x = 0.0;
    double y = 0.0;
    if (!parseDouble(field(row, "sqrt_sNN_GeV"), x) || !parseDouble(field(row, "value"), y)) continue;
    if (x <= 0.0 || y <= 0.0) continue;
    values[particle][energyKey(x)] = {x, y, rowError(row), true};
  }
  return values;
}

std::map<std::string, Series> buildRatioSeries(
    const std::vector<CsvRow> &summary,
    const std::vector<CsvRow> &all_rows)
{
  std::map<std::string, Series> ratios;
  const auto light = readLightYieldMaps(summary);
  const auto proton_it = light.find("p");
  if (proton_it == light.end()) return ratios;
  const auto &proton = proton_it->second;

  for (const std::string &particle : {"d", "t", "He3"}) {
    const auto lit = light.find(particle);
    if (lit == light.end()) continue;
    for (const auto &entry : lit->second) {
      const auto pit = proton.find(entry.first);
      if (pit == proton.end()) continue;
      if (particle == "He3" && std::fabs(entry.second.x - 3.0) < 1.0e-6) continue;
      addRatioPoint(ratios, ratioLabelForNumerator(particle), entry.second, pit->second);
    }
  }

  const std::string h3l_ratio_label = ratioLabelForNumerator("H3Lambda");
  for (const CsvRow &row : all_rows) {
    if (field(row, "particle") != "H3Lambda" || field(row, "observable") != "H3Lambda/Lambda") continue;
    double x = 0.0;
    double y = 0.0;
    if (!parseDouble(field(row, "sqrt_sNN_GeV"), x) || !parseDouble(field(row, "value"), y)) continue;
    if (x <= 0.0 || y <= 0.0) continue;
    Series &series = ratios[h3l_ratio_label];
    if (series.label.empty()) ratioStyle(h3l_ratio_label, series);
    series.x.push_back(x);
    series.y.push_back(y);
    series.ey.push_back(rowError(row));
  }

  for (auto &entry : ratios) sortSeries(entry.second);
  return ratios;
}

std::vector<std::string> sortedKeys(const std::map<std::string, Series> &series_map)
{
  std::vector<std::string> keys;
  for (const auto &entry : series_map) keys.push_back(entry.first);
  std::sort(keys.begin(), keys.end(), [&](const std::string &a, const std::string &b) {
    const std::string &la = series_map.at(a).label;
    const std::string &lb = series_map.at(b).label;
    return la < lb;
  });
  return keys;
}

TGraphErrors *makeGraph(const Series &series)
{
  if (series.x.empty()) return nullptr;
  std::vector<double> ex(series.x.size(), 0.0);
  TGraphErrors *graph = new TGraphErrors(
      static_cast<int>(series.x.size()),
      const_cast<double *>(series.x.data()),
      const_cast<double *>(series.y.data()),
      ex.data(),
      const_cast<double *>(series.ey.data()));
  graph->SetMarkerStyle(series.marker);
  graph->SetMarkerSize(series.marker_size);
  graph->SetMarkerColor(series.color);
  graph->SetLineColor(series.color);
  graph->SetLineWidth(series.line_width);
  return graph;
}

bool isThermusDisplayedParticle(const std::string &particle)
{
  return particle == "p" || particle == "Lambda" ||
         particle == "d" || particle == "t" || particle == "He3" || particle == "H3Lambda";
}

int thermusColorForParticle(const std::string &particle)
{
  if (particle == "p") return kRed + 1;
  if (particle == "Lambda") return kAzure + 2;
  if (particle == "d") return kBlue + 1;
  if (particle == "t") return kOrange + 7;
  if (particle == "He3") return kGreen + 2;
  if (particle == "H3Lambda") return kViolet + 5;
  return kGray + 2;
}

bool isThermusSceModel(const std::string &model_set)
{
  return model_set.find("THERMUS SCE") != std::string::npos;
}

bool isThermusGceModel(const std::string &model_set)
{
  return model_set.find("THERMUS GCE") != std::string::npos;
}

std::string thermusLegendLabel(const std::string &model_set)
{
  if (isThermusSceModel(model_set)) return "THERMUS SCE";
  if (isThermusGceModel(model_set)) return "THERMUS GCE";
  return model_set;
}

int thermusLineStyle(const std::string &model_set)
{
  if (isThermusSceModel(model_set)) return 2;
  if (isThermusGceModel(model_set)) return 9;
  return 1;
}

void sortModelCurve(ModelCurve &curve)
{
  struct Point {
    double x;
    double y;
  };
  std::vector<Point> points;
  for (std::size_t i = 0; i < curve.x.size(); ++i) {
    points.push_back({curve.x[i], curve.y[i]});
  }
  std::sort(points.begin(), points.end(), [](const Point &a, const Point &b) { return a.x < b.x; });
  for (std::size_t i = 0; i < points.size(); ++i) {
    curve.x[i] = points[i].x;
    curve.y[i] = points[i].y;
  }
}

std::vector<ModelCurve> readThermusCurves(const TString &path)
{
  std::map<std::string, ModelCurve> curve_map;
  const std::vector<CsvRow> rows = readCsv(path);

  for (const CsvRow &row : rows) {
    const std::string model_set = field(row, "model_set");
    const std::string particle = field(row, "particle");
    if (!isThermusDisplayedParticle(particle)) continue;
    if (!isThermusGceModel(model_set)) continue;

    double x = 0.0;
    double y = 0.0;
    if (!parseDouble(field(row, "energy_GeV"), x) || !parseDouble(field(row, "model_yield"), y)) continue;
    if (x <= 0.0 || y <= 0.0 || x > 6000.0) continue;

    const std::string key = model_set + "|" + particle;
    ModelCurve &curve = curve_map[key];
    if (curve.label.empty()) {
      curve.label = thermusLegendLabel(model_set);
      curve.color = thermusColorForParticle(particle);
      curve.style = thermusLineStyle(model_set);
      curve.width = 3;
    }
    curve.x.push_back(x);
    curve.y.push_back(y);
  }

  std::vector<ModelCurve> curves;
  std::set<std::string> seen_legend_labels;
  for (auto &entry : curve_map) {
    sortModelCurve(entry.second);
    entry.second.show_in_legend = seen_legend_labels.insert(entry.second.label).second;
    curves.push_back(entry.second);
  }
  return curves;
}

void addModelRatioPoint(
    std::map<std::string, ModelCurve> &curve_map,
    const std::string &ratio,
    double x,
    double y)
{
  if (x <= 0.0 || y <= 0.0) return;
  ModelCurve &curve = curve_map[ratio];
  if (curve.label.empty()) {
    curve.label = "THERMUS GCE d/p";
    curve.color = colorForRatio(ratio);
    curve.style = 9;
    curve.width = 3;
  }
  curve.x.push_back(x);
  curve.y.push_back(y);
}

std::vector<ModelCurve> readThermusRatioCurves(const TString &path)
{
  std::map<std::string, std::map<std::string, YieldPoint>> by_energy;
  for (const CsvRow &row : readCsv(path)) {
    if (!isThermusGceModel(field(row, "model_set"))) continue;
    const std::string particle = field(row, "particle");
    if (particle != "p" && particle != "d") {
      continue;
    }

    double x = 0.0;
    double y = 0.0;
    if (!parseDouble(field(row, "energy_GeV"), x) || !parseDouble(field(row, "model_yield"), y)) continue;
    if (x <= 0.0 || y <= 0.0 || x > kRatioXMax) continue;
    by_energy[energyKey(x)][particle] = {x, y, 0.0, true};
  }

  std::map<std::string, ModelCurve> curve_map;
  for (const auto &entry : by_energy) {
    const auto &values = entry.second;
    const auto pit = values.find("p");
    if (pit != values.end()) {
      for (const std::string &particle : {"d"}) {
        const auto nit = values.find(particle);
        if (nit == values.end()) continue;
        addModelRatioPoint(curve_map, ratioLabelForNumerator(particle), nit->second.x, nit->second.y / pit->second.y);
      }
    }
  }

  std::vector<ModelCurve> curves;
  bool show_thermus_legend = true;
  for (auto &entry : curve_map) {
    sortModelCurve(entry.second);
    entry.second.show_in_legend = show_thermus_legend;
    show_thermus_legend = false;
    curves.push_back(entry.second);
  }
  return curves;
}

ModelCurve readModelCurve(const TString &path, const std::string &label, int color, int style)
{
  ModelCurve curve;
  curve.label = label;
  curve.color = color;
  curve.style = style;
  curve.width = 3;

  std::ifstream input(path.Data());
  if (!input) return curve;

  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;

    std::istringstream row(line);
    std::string a;
    std::string b;
    row >> a >> b;

    double x = 0.0;
    double y = 0.0;
    if (parseDouble(a, x) && parseDouble(b, y) && y > 0.0) {
      curve.x.push_back(x);
      curve.y.push_back(y);
    }
  }
  return curve;
}

std::map<std::string, YieldPoint> readModelPointMap(const TString &path)
{
  std::map<std::string, YieldPoint> values;
  std::ifstream input(path.Data());
  if (!input) return values;

  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t\r\n");
    if (first == std::string::npos || line[first] == '#') continue;

    std::istringstream row(line);
    std::string a;
    std::string b;
    row >> a >> b;

    double x = 0.0;
    double y = 0.0;
    if (parseDouble(a, x) && parseDouble(b, y) && x > 0.0 && y > 0.0) {
      values[energyKey(x)] = {x, y, 0.0, true};
    }
  }
  return values;
}

ModelCurve readModelProductCurve(
    const TString &left_path,
    const TString &right_path,
    const std::string &label,
    int color,
    int style)
{
  ModelCurve curve;
  curve.label = label;
  curve.color = color;
  curve.style = style;
  curve.width = 3;

  const auto left = readModelPointMap(left_path);
  const auto right = readModelPointMap(right_path);
  for (const auto &entry : left) {
    const auto rit = right.find(entry.first);
    if (rit == right.end()) continue;
    curve.x.push_back(entry.second.x);
    curve.y.push_back(entry.second.y * rit->second.y);
  }
  sortModelCurve(curve);
  return curve;
}

ModelCurve scaleModelCurve(ModelCurve curve, double scale)
{
  for (double &y : curve.y) y *= scale;
  return curve;
}

TGraph *makeModelGraph(const ModelCurve &curve)
{
  if (curve.x.empty()) return nullptr;
  TGraph *graph = new TGraph(
      static_cast<int>(curve.x.size()),
      const_cast<double *>(curve.x.data()),
      const_cast<double *>(curve.y.data()));
  graph->SetLineColor(curve.color);
  graph->SetLineStyle(curve.style);
  graph->SetLineWidth(curve.width);
  return graph;
}

int legendColorForCurve(const ModelCurve &curve)
{
  if (curve.label.find("THERMUS") == 0) return kBlack;
  return curve.color;
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

void drawText(double x, double y, const char *text, double size = 0.045)
{
  TLatex *latex = new TLatex(x, y, text);
  latex->SetTextFont(42);
  latex->SetTextSize(size);
  latex->Draw("same");
}

void drawPanel(
    const char *frame_name,
    const std::map<std::string, Series> &series_map,
    const std::vector<ModelCurve> &models,
    double xmin,
    double xmax,
    double ymin,
    double ymax,
    const char *ylabel,
    const char *title,
    bool draw_x_title,
    double title_x,
    double title_y)
{
  TH1D *frame = new TH1D(frame_name, "", 1, xmin, xmax);
  frame->SetMinimum(ymin);
  frame->SetMaximum(ymax);
  frame->GetXaxis()->CenterTitle();
  frame->GetYaxis()->CenterTitle();
  frame->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  frame->GetYaxis()->SetTitle(ylabel);
  frame->GetXaxis()->SetTitleSize(draw_x_title ? 0.055 : 0.0);
  frame->GetXaxis()->SetLabelSize(draw_x_title ? 0.045 : 0.0);
  frame->GetXaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetTitleSize(0.057);
  frame->GetYaxis()->SetTitleOffset(1.05);
  frame->GetYaxis()->SetLabelSize(0.043);
  frame->Draw();

  drawText(title_x, title_y, title, 0.038);

  std::vector<TGraph *> model_graphs;
  for (const ModelCurve &curve : models) {
    TGraph *graph = makeModelGraph(curve);
    if (!graph) continue;
    model_graphs.push_back(graph);
    graph->Draw("C SAME");
  }

  std::vector<TGraphErrors *> data_graphs;
  for (const std::string &key : sortedKeys(series_map)) {
    const Series &series = series_map.at(key);
    TGraphErrors *graph = makeGraph(series);
    if (!graph) continue;
    data_graphs.push_back(graph);
    graph->Draw(series.draw_line && series.x.size() > 1 ? "LP SAME" : "P SAME");
  }

  drawFrameBox(xmin, xmax, ymin, ymax);
}

void drawLegends(
    const std::map<std::string, Series> &rapidity,
    const std::vector<ModelCurve> &models)
{
  TLegend *model_leg = new TLegend(0.58, 0.76, 0.94, 0.94);
  model_leg->SetFillStyle(4000);
  model_leg->SetBorderSize(0);
  model_leg->SetTextSize(0.027);
  model_leg->SetHeader("Model curves", "C");
  for (const ModelCurve &curve : models) {
    if (!curve.show_in_legend) continue;
    TGraph *dummy = new TGraph();
    dummy->SetLineColor(legendColorForCurve(curve));
    dummy->SetLineStyle(curve.style);
    dummy->SetLineWidth(curve.width);
    model_leg->AddEntry(dummy, curve.label.c_str(), "l");
  }
  model_leg->Draw();

  TLegend *data_leg = new TLegend(0.58, 0.49, 0.94, 0.73);
  data_leg->SetFillStyle(4000);
  data_leg->SetBorderSize(0);
  data_leg->SetTextSize(0.024);
  data_leg->SetNColumns(2);
  data_leg->SetHeader("Data", "C");

  std::set<std::string> seen;
  for (const std::string &key : sortedKeys(rapidity)) {
    const Series &series = rapidity.at(key);
    if (seen.count(series.label)) continue;
    seen.insert(series.label);
    TGraphErrors *dummy = new TGraphErrors();
    dummy->SetMarkerStyle(series.marker);
    dummy->SetMarkerSize(series.marker_size * 0.85);
    dummy->SetMarkerColor(series.color);
    dummy->SetLineColor(series.color);
    dummy->SetLineWidth(series.line_width);
    data_leg->AddEntry(dummy, series.label.c_str(), "p");
  }
  data_leg->Draw();
}

}  // namespace

void plot_light_nuclei_hypernuclei_yields_root()
{
  applyFreezeOutStyle();

  TString macro_path = __FILE__;
  if (!macro_path.BeginsWith("/")) {
    macro_path = TString(gSystem->WorkingDirectory()) + "/" + macro_path;
  }
  const TString script_dir = gSystem->DirName(macro_path);
  const TString repo_dir = gSystem->DirName(script_dir);
  const TString input_path = repo_dir + "/data/light_nuclei_hypernuclei_yields_vs_energy.csv";
  const TString strange_input_path = repo_dir + "/data/strange_hadron_yields_gc_sce_overlay.csv";
  const TString model_dir = repo_dir + "/data/data_for_Xin/model";
  const TString thermus_path = repo_dir + "/data/thermus_light_hypernuclei_predictions/light_hypernuclei_thermus_predictions.csv";
  const TString out_pdf = repo_dir + "/data/light_nuclei_hypernuclei_yields_root_style.pdf";
  const TString out_png = repo_dir + "/data/light_nuclei_hypernuclei_yields_root_style.png";
  const TString ratio_out_pdf = repo_dir + "/data/light_nuclei_to_hadron_ratios_root_style.pdf";
  const TString ratio_out_png = repo_dir + "/data/light_nuclei_to_hadron_ratios_root_style.png";
  const TString ratio_sqrts_out_pdf = repo_dir + "/data/light_nuclei_to_hadron_ratios_vovchenko_sqrts.pdf";
  const TString ratio_sqrts_out_png = repo_dir + "/data/light_nuclei_to_hadron_ratios_vovchenko_sqrts.png";

  const std::vector<CsvRow> rows = readCsv(input_path);
  const std::vector<CsvRow> summary = selectSummaryRows(rows);
  if (summary.empty()) {
    std::cerr << "No summary rows loaded from " << input_path << std::endl;
    return;
  }

  std::map<std::string, Series> rapidity = buildSeries(summary, false);
  appendOverlayYieldData(strange_input_path, rapidity);
  std::vector<ModelCurve> models = readThermusCurves(thermus_path);
  models.push_back(readModelCurve(model_dir + "/Thermal_fist_H3L_dNdy.txt", "Thermal-FIST H3L", kBlack, 1));
  models.push_back(readModelCurve(model_dir + "/UrQMD_coal_H3L_dNdy.txt", "UrQMD+coal. H3L", kTeal + 3, 3));

  TCanvas *canvas = new TCanvas("cLightNucleiHypernucleiYields", "", 1100, 700);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.13);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  drawPanel(
      "frameLightNucleiRapidity",
      rapidity,
      models,
      2.4,
      6000.0,
      1.0e-4,
      3.0e2,
      "dN/dy",
      "Hadron and light-nuclei rapidity-density yields",
      true,
      4.2,
      135.0);
  drawLegends(rapidity, models);

  canvas->Update();
  canvas->SaveAs(out_pdf);
  canvas->SaveAs(out_png);

  std::cout << "wrote " << out_pdf << std::endl;
  std::cout << "wrote " << out_png << std::endl;

  std::map<std::string, Series> ratios = buildRatioSeries(summary, rows);
  std::vector<ModelCurve> ratio_models;
  ratio_models.push_back(readModelCurve(
      model_dir + "/Thermal_model_deuteron2proton_star2311_fig12.txt",
      "Thermal model d/p",
      kBlue + 1,
      1));
  ratio_models.push_back(scaleModelCurve(
      readModelCurve(
          model_dir + "/Thermal_model_triton2proton_star2311_fig12.txt",
          "Thermal model 5#times t/p",
          kOrange + 7,
          1),
      kTpVisibilityScale));
  ratio_models.push_back(scaleModelCurve(
      readModelCurve(
          model_dir + "/UrQMD_coal_triton2proton.txt",
          "UrQMD+coal. 5#times t/p",
          kOrange + 7,
          3),
      kTpVisibilityScale));
  ratio_models.push_back(readModelProductCurve(
      model_dir + "/Thermal_fist_S3.txt",
      model_dir + "/Thermal_fist_triton2proton.txt",
      "Thermal-FIST {}_{#Lambda}^{3}H/#Lambda",
      kViolet + 5,
      1));
  ratio_models.push_back(readModelCurve(
      model_dir + "/UrQMD_coal_H3L2Lambda.txt",
      "UrQMD+coal. {}_{#Lambda}^{3}H/#Lambda",
      kViolet + 5,
      3));

  TCanvas *ratio_canvas = new TCanvas("cLightNucleiToHadronRatios", "", 904, 928);
  gPad->SetLeftMargin(0.13);
  gPad->SetRightMargin(0.04);
  gPad->SetTopMargin(0.03);
  gPad->SetBottomMargin(0.13);
  gPad->SetLogx(1);
  gPad->SetLogy(1);
  drawPanel(
      "frameLightNucleiRatios",
      ratios,
      ratio_models,
      2.4,
      kRatioXMax,
      kRatioYMin,
      4.0e-1,
      "yield ratio",
      "Light-nuclei-to-hadron ratios",
      true,
      4.2,
      1.7e-1);
  drawLegends(ratios, ratio_models);

  ratio_canvas->Update();
  ratio_canvas->SaveAs(ratio_out_pdf);
  ratio_canvas->SaveAs(ratio_out_png);
  ratio_canvas->SaveAs(ratio_sqrts_out_pdf);
  ratio_canvas->SaveAs(ratio_sqrts_out_png);

  std::cout << "wrote " << ratio_out_pdf << std::endl;
  std::cout << "wrote " << ratio_out_png << std::endl;
  std::cout << "wrote " << ratio_sqrts_out_pdf << std::endl;
  std::cout << "wrote " << ratio_sqrts_out_png << std::endl;
}
