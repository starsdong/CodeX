#include <TCanvas.h>
#include <TClass.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TKey.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct SliceFitResult {
  bool ok = false;
  int fitStatus = -1;
  int iterations = 0;
  double entries = 0.0;
  double mean = std::numeric_limits<double>::quiet_NaN();
  double meanErr = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double sigmaErr = std::numeric_limits<double>::quiet_NaN();
  double chi2Ndf = std::numeric_limits<double>::quiet_NaN();
};

struct FwhmResult {
  bool ok = false;
  double fwhm = std::numeric_limits<double>::quiet_NaN();
  double sigma = std::numeric_limits<double>::quiet_NaN();
  double left = std::numeric_limits<double>::quiet_NaN();
  double right = std::numeric_limits<double>::quiet_NaN();
};

std::string CleanName(const std::string &name)
{
  std::string out;
  out.reserve(name.size());
  bool lastWasUnderscore = false;
  for (unsigned char ch : name) {
    if (ch == '+') {
      out += "plus";
      lastWasUnderscore = false;
      continue;
    }
    if (ch == '-') {
      out += "minus";
      lastWasUnderscore = false;
      continue;
    }

    const bool keep = std::isalnum(ch) || ch == '_';
    if (keep) {
      out.push_back(static_cast<char>(ch));
      lastWasUnderscore = false;
    } else if (!lastWasUnderscore) {
      out.push_back('_');
      lastWasUnderscore = true;
    }
  }
  while (!out.empty() && out.front() == '_') out.erase(out.begin());
  while (!out.empty() && out.back() == '_') out.pop_back();
  return out.empty() ? "hist" : out;
}

std::string StripRootExtension(const char *path)
{
  TString base = gSystem->BaseName(path);
  if (base.EndsWith(".root")) base.Resize(base.Length() - 5);
  return CleanName(base.Data());
}

std::vector<std::string> SplitRequestedHistograms(const std::string &requested)
{
  std::vector<std::string> names;
  std::stringstream stream(requested);
  std::string item;

  while (std::getline(stream, item, ',')) {
    const auto begin = item.find_first_not_of(" \t\n\r");
    const auto end = item.find_last_not_of(" \t\n\r");
    if (begin == std::string::npos || end == std::string::npos) continue;
    names.push_back(item.substr(begin, end - begin + 1));
  }

  return names;
}

std::string CapitalizeFirst(std::string text)
{
  if (!text.empty()) {
    text[0] = static_cast<char>(std::toupper(static_cast<unsigned char>(text[0])));
  }
  return text;
}

std::vector<std::string> GaussianWidthCandidates(const std::string &histName)
{
  std::vector<std::string> candidates;
  candidates.push_back("Gaussian_Width_of_" + histName);

  if (histName.rfind("h_", 0) == 0 && histName.size() > 2) {
    candidates.push_back("Gaussian_Width_of_" + CapitalizeFirst(histName.substr(2)));
  }

  if (histName.rfind("h", 0) == 0 && histName.size() > 1) {
    candidates.push_back("Gaussian_Width_of_" + histName.substr(1));
  }

  return candidates;
}

TH1 *FindInputGaussianWidth(TDirectory *dir, const std::string &histName)
{
  if (!dir) return nullptr;

  for (const auto &candidate : GaussianWidthCandidates(histName)) {
    TH1 *hist = dynamic_cast<TH1 *>(dir->Get(candidate.c_str()));
    if (hist) return hist;
  }

  return nullptr;
}

double RobustInitialSigma(TH1 *h)
{
  const double rms = h->GetRMS();
  if (std::isfinite(rms) && rms > 0.0) return rms;

  const int maxBin = h->GetMaximumBin();
  const double binWidth = h->GetXaxis()->GetBinWidth(maxBin);
  return (binWidth > 0.0) ? binWidth : 1.0;
}

bool HasEnoughFitContent(TH1 *h, double fitMin, double fitMax, int minPopulatedBins = 3)
{
  if (!h || !(fitMax > fitMin)) return false;

  const int firstBin = std::max(1, h->GetXaxis()->FindFixBin(fitMin));
  const int lastBin = std::min(h->GetNbinsX(), h->GetXaxis()->FindFixBin(fitMax));
  double sum = 0.0;
  int populatedBins = 0;

  for (int bin = firstBin; bin <= lastBin; ++bin) {
    const double content = h->GetBinContent(bin);
    if (content > 0.0) {
      sum += content;
      ++populatedBins;
    }
  }

  return sum > 0.0 && populatedBins >= minPopulatedBins;
}

double InterpolateX(double x1, double y1, double x2, double y2, double y)
{
  if (std::fabs(y2 - y1) <= 0.0) return 0.5 * (x1 + x2);
  return x1 + (y - y1) * (x2 - x1) / (y2 - y1);
}

FwhmResult CalculateFwhmSigma(TH1D *slice, double minEntries = 50.0)
{
  FwhmResult result;
  if (!slice || slice->GetEntries() < minEntries || slice->Integral() <= 0.0) return result;

  const int maxBin = slice->GetMaximumBin();
  const double maximum = slice->GetBinContent(maxBin);
  if (!(maximum > 0.0)) return result;

  const double halfMaximum = 0.5 * maximum;
  int leftAbove = maxBin;
  while (leftAbove > 1 && slice->GetBinContent(leftAbove - 1) >= halfMaximum) {
    --leftAbove;
  }

  int rightAbove = maxBin;
  const int nBins = slice->GetNbinsX();
  while (rightAbove < nBins && slice->GetBinContent(rightAbove + 1) >= halfMaximum) {
    ++rightAbove;
  }

  if (leftAbove <= 1 || rightAbove >= nBins) return result;

  result.left = InterpolateX(slice->GetBinCenter(leftAbove - 1),
                             slice->GetBinContent(leftAbove - 1),
                             slice->GetBinCenter(leftAbove),
                             slice->GetBinContent(leftAbove),
                             halfMaximum);
  result.right = InterpolateX(slice->GetBinCenter(rightAbove),
                              slice->GetBinContent(rightAbove),
                              slice->GetBinCenter(rightAbove + 1),
                              slice->GetBinContent(rightAbove + 1),
                              halfMaximum);
  result.fwhm = result.right - result.left;
  if (!(result.fwhm > 0.0)) return result;

  result.sigma = result.fwhm / 2.354820045;
  result.ok = std::isfinite(result.sigma) && result.sigma > 0.0;
  return result;
}

SliceFitResult FitSliceIterative(TH1D *slice,
                                 double nSigma = 2.0,
                                 double relSigmaTolerance = 0.005,
                                 int maxIterations = 10,
                                 double minEntries = 50.0)
{
  SliceFitResult result;
  result.entries = slice ? slice->GetEntries() : 0.0;

  if (!slice || result.entries < minEntries || slice->Integral() <= 0.0) {
    return result;
  }

  const double axisMin = slice->GetXaxis()->GetXmin();
  const double axisMax = slice->GetXaxis()->GetXmax();
  const double axisWidth = axisMax - axisMin;
  if (!(axisWidth > 0.0)) return result;

  double fitMin = axisMin;
  double fitMax = axisMax;
  double previousSigma = -1.0;
  const double absSigmaTolerance = axisWidth * 1.0e-6;

  TF1 fitFunc(Form("fit_%s", slice->GetName()), "gaus", fitMin, fitMax);
  fitFunc.SetNpx(500);
  fitFunc.SetParameters(slice->GetMaximum(), slice->GetMean(), RobustInitialSigma(slice));
  fitFunc.SetParLimits(2, axisWidth * 1.0e-9, axisWidth * 10.0);

  for (int iter = 0; iter < maxIterations; ++iter) {
    if (!HasEnoughFitContent(slice, fitMin, fitMax)) break;

    fitFunc.SetRange(fitMin, fitMax);
    const TFitResultPtr fitResult = slice->Fit(&fitFunc, "Q0RS");
    const int status = static_cast<int>(fitResult);

    const double mean = fitFunc.GetParameter(1);
    const double sigma = std::fabs(fitFunc.GetParameter(2));
    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.0) {
      result.fitStatus = status;
      result.iterations = iter + 1;
      return result;
    }

    result.ok = (status == 0);
    result.fitStatus = status;
    result.iterations = iter + 1;
    result.mean = mean;
    result.meanErr = fitFunc.GetParError(1);
    result.sigma = sigma;
    result.sigmaErr = fitFunc.GetParError(2);
    result.chi2Ndf = (fitFunc.GetNDF() > 0) ? fitFunc.GetChisquare() / fitFunc.GetNDF()
                                            : std::numeric_limits<double>::quiet_NaN();

    if (iter > 0) {
      const double allowedChange = std::max(absSigmaTolerance, relSigmaTolerance * previousSigma);
      if (std::fabs(sigma - previousSigma) <= allowedChange) break;
    }

    previousSigma = sigma;
    fitMin = std::max(axisMin, mean - nSigma * sigma);
    fitMax = std::min(axisMax, mean + nSigma * sigma);

    if (!(fitMax > fitMin)) break;
    if (slice->GetXaxis()->FindBin(fitMax) - slice->GetXaxis()->FindBin(fitMin) < 2) break;
    if (!HasEnoughFitContent(slice, fitMin, fitMax)) break;

    fitFunc.SetParameters(fitFunc.GetParameter(0), mean, sigma);
  }

  return result;
}

TH2D *BuildGroupedXHistogram(TH2 *h2,
                             const std::string &cleanName,
                             const std::vector<double> &edges,
                             int xBinGroup)
{
  if (!h2 || edges.size() < 2) return nullptr;

  const int nOutBins = static_cast<int>(edges.size()) - 1;
  const int nY = h2->GetYaxis()->GetNbins();
  TH2D *grouped = nullptr;
  if (h2->GetYaxis()->GetXbins()->GetSize() > 0) {
    grouped = new TH2D(("h2_grouped_" + cleanName).c_str(), h2->GetTitle(),
                       nOutBins, edges.data(), nY, h2->GetYaxis()->GetXbins()->GetArray());
  } else {
    grouped = new TH2D(("h2_grouped_" + cleanName).c_str(), h2->GetTitle(),
                       nOutBins, edges.data(),
                       nY, h2->GetYaxis()->GetXmin(), h2->GetYaxis()->GetXmax());
  }
  grouped->SetDirectory(nullptr);
  grouped->Sumw2();
  grouped->GetXaxis()->SetTitle(h2->GetXaxis()->GetTitle());
  grouped->GetYaxis()->SetTitle(h2->GetYaxis()->GetTitle());

  const int nX = h2->GetXaxis()->GetNbins();
  for (int outXBin = 1; outXBin <= nOutBins; ++outXBin) {
    const int firstXBin = (outXBin - 1) * xBinGroup + 1;
    const int lastXBin = std::min(nX, outXBin * xBinGroup);
    for (int yBin = 1; yBin <= nY; ++yBin) {
      double content = 0.0;
      double error2 = 0.0;
      for (int xBin = firstXBin; xBin <= lastXBin; ++xBin) {
        content += h2->GetBinContent(xBin, yBin);
        const double error = h2->GetBinError(xBin, yBin);
        error2 += error * error;
      }
      grouped->SetBinContent(outXBin, yBin, content);
      grouped->SetBinError(outXBin, yBin, std::sqrt(error2));
    }
  }
  grouped->SetEntries(h2->GetEntries());

  return grouped;
}

TH1D *BuildRootFitSlicesYSigma(TH2 *h2,
                               const std::string &cleanName,
                               const std::vector<double> &edges,
                               int xBinGroup,
                               double minEntries)
{
  TH1D *hFitSlicesSigma = new TH1D(("h_fitslicesy_sigma_" + cleanName).c_str(),
                                   Form("%s;%.120s;ROOT FitSlicesY Gaussian #sigma of %.120s",
                                        h2->GetTitle(),
                                        h2->GetXaxis()->GetTitle(),
                                        h2->GetYaxis()->GetTitle()),
                                   static_cast<int>(edges.size()) - 1, edges.data());
  hFitSlicesSigma->SetDirectory(nullptr);

  TH2D *grouped = BuildGroupedXHistogram(h2, cleanName, edges, xBinGroup);
  if (!grouped) return hFitSlicesSigma;

  TObjArray fitSlices;
  fitSlices.SetOwner(kTRUE);
  grouped->FitSlicesY(nullptr, 1, grouped->GetXaxis()->GetNbins(),
                      static_cast<int>(std::max(0.0, minEntries)), "QNR", &fitSlices);

  TH1 *sigma = dynamic_cast<TH1 *>(fitSlices.At(2));
  if (sigma) {
    for (int bin = 1; bin <= hFitSlicesSigma->GetNbinsX(); ++bin) {
      const double value = std::fabs(sigma->GetBinContent(bin));
      if (std::isfinite(value) && value > 0.0) {
        hFitSlicesSigma->SetBinContent(bin, value);
        hFitSlicesSigma->SetBinError(bin, sigma->GetBinError(bin));
      }
    }
  }

  delete grouped;
  return hFitSlicesSigma;
}

TH1D *BuildInputGaussianWidthHist(TH1 *source,
                                  const std::string &cleanName,
                                  const std::vector<double> &edges,
                                  const char *xTitle,
                                  const char *yTitle)
{
  TH1D *hist = new TH1D(("h_input_gaussian_width_" + cleanName).c_str(),
                        Form("%s;%.120s;Input Gaussian width of %.120s",
                             source ? source->GetTitle() : "Input Gaussian width",
                             xTitle ? xTitle : "",
                             yTitle ? yTitle : ""),
                        static_cast<int>(edges.size()) - 1, edges.data());
  hist->SetDirectory(nullptr);

  if (!source) return hist;

  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const double x = hist->GetBinCenter(bin);
    const int sourceBin = source->GetXaxis()->FindFixBin(x);
    if (sourceBin < 1 || sourceBin > source->GetNbinsX()) continue;

    const double value = source->GetBinContent(sourceBin);
    if (!(value > 0.0) || !std::isfinite(value)) continue;

    hist->SetBinContent(bin, value);
    hist->SetBinError(bin, source->GetBinError(sourceBin));
  }

  return hist;
}

void CollectTH2(TDirectory *dir, std::vector<std::string> &names, const std::string &prefix = "")
{
  if (!dir) return;

  TIter next(dir->GetListOfKeys());
  while (TKey *key = static_cast<TKey *>(next())) {
    TClass *cl = gROOT->GetClass(key->GetClassName());
    if (!cl) continue;

    const std::string keyName = key->GetName();
    const std::string fullName = prefix.empty() ? keyName : prefix + "/" + keyName;

    if (cl->InheritsFrom(TDirectory::Class())) {
      TDirectory *subdir = dynamic_cast<TDirectory *>(key->ReadObj());
      CollectTH2(subdir, names, fullName);
      continue;
    }

    if (cl->InheritsFrom(TH2::Class())) {
      names.push_back(fullName);
    }
  }
}

std::vector<double> MakeGroupedEdges(const TAxis *axis, int xBinGroup)
{
  std::vector<double> edges;
  const int nX = axis->GetNbins();
  for (int first = 1; first <= nX; first += xBinGroup) {
    edges.push_back(axis->GetBinLowEdge(first));
  }
  edges.push_back(axis->GetBinUpEdge(nX));
  return edges;
}

double HistogramMaxWithErrors(const TH1 *hist)
{
  double maximum = 0.0;
  if (!hist) return maximum;

  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const double value = hist->GetBinContent(bin);
    if (value <= 0.0 || !std::isfinite(value)) continue;
    maximum = std::max(maximum, value + hist->GetBinError(bin));
  }
  return maximum;
}

TGraphErrors *MakeGraphFromHistogram(const TH1D *hist, const std::string &name,
                                     int markerStyle, int color)
{
  auto *graph = new TGraphErrors();
  graph->SetName(name.c_str());
  graph->SetMarkerStyle(markerStyle);
  graph->SetMarkerSize(0.9);
  graph->SetMarkerColor(color);
  graph->SetLineColor(color);
  graph->SetLineWidth(2);

  if (!hist) return graph;
  for (int bin = 1; bin <= hist->GetNbinsX(); ++bin) {
    const double value = hist->GetBinContent(bin);
    if (!(value > 0.0) || !std::isfinite(value)) continue;

    const int point = graph->GetN();
    graph->SetPoint(point, hist->GetBinCenter(bin), value);
    graph->SetPointError(point, 0.0, hist->GetBinError(bin));
  }

  return graph;
}

void DrawResolutionCanvas(TH1D *hIterativeSigma,
                          TH1D *hFitSlicesSigma,
                          TH1D *hFwhmSigma,
                          TH1D *hInputGaussianWidth,
                          const std::string &canvasName,
                          const std::string &outDir,
                          double yAxisMax = -1.0)
{
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);

  TCanvas canvas(("c_" + canvasName).c_str(), "", 900, 700);
  canvas.SetLeftMargin(0.13);
  canvas.SetRightMargin(0.04);
  canvas.SetTopMargin(0.07);
  canvas.SetBottomMargin(0.12);
  canvas.SetGridy();

  TH1D *frame = dynamic_cast<TH1D *>(hIterativeSigma->Clone(("frame_" + canvasName).c_str()));
  frame->Reset("ICES");
  frame->SetDirectory(nullptr);
  frame->SetMinimum(0.0);
  const double autoMaximum = 1.2 * std::max({HistogramMaxWithErrors(hIterativeSigma),
                                             HistogramMaxWithErrors(hFitSlicesSigma),
                                             HistogramMaxWithErrors(hFwhmSigma),
                                             HistogramMaxWithErrors(hInputGaussianWidth),
                                             1.0e-12});
  frame->SetMaximum((yAxisMax > 0.0) ? yAxisMax : autoMaximum);
  frame->SetLineColor(kWhite);
  frame->SetMarkerColor(kWhite);
  frame->Draw("HIST");

  TGraphErrors *gIterative = MakeGraphFromHistogram(hIterativeSigma, "g_iterative_" + canvasName,
                                                    20, kBlack);
  TGraphErrors *gFitSlices = MakeGraphFromHistogram(hFitSlicesSigma, "g_fitslicesy_" + canvasName,
                                                    24, kBlue + 1);
  TGraphErrors *gFwhm = MakeGraphFromHistogram(hFwhmSigma, "g_fwhm_" + canvasName,
                                               25, kRed + 1);
  TGraphErrors *gInput = MakeGraphFromHistogram(hInputGaussianWidth,
                                                "g_input_gaussian_width_" + canvasName,
                                                26, kGreen + 2);

  gIterative->Draw("P SAME");
  gFitSlices->Draw("P SAME");
  gFwhm->Draw("P SAME");
  if (gInput->GetN() > 0) gInput->Draw("P SAME");

  TLegend legend(0.51, 0.66, 0.90, 0.88);
  legend.SetBorderSize(0);
  legend.SetFillStyle(0);
  legend.SetTextSize(0.030);
  legend.AddEntry(gIterative, "Iterative Gaussian", "lep");
  legend.AddEntry(gFitSlices, "ROOT FitSlicesY()", "lep");
  legend.AddEntry(gFwhm, "FWHM / 2.355", "lep");
  if (gInput->GetN() > 0) {
    legend.AddEntry(gInput, "Input Gaussian width", "lep");
  }
  legend.Draw();
  frame->Draw("AXIS SAME");

  canvas.Print(Form("%s/%s_resolution_comparison_vs_x.pdf", outDir.c_str(), canvasName.c_str()));
  canvas.Print(Form("%s/%s_resolution_comparison_vs_x.png", outDir.c_str(), canvasName.c_str()));
  canvas.Write(("c_" + canvasName + "_resolution_comparison_vs_x").c_str());

  delete gIterative;
  delete gFitSlices;
  delete gFwhm;
  delete gInput;
  delete frame;
}

void ProcessOneTH2(TH2 *h2,
                   TDirectory *inputDir,
                   TFile *outFile,
                   const std::string &plotDir,
                   int xBinGroup,
                   double nSigma,
                   double relSigmaTolerance,
                   int maxIterations,
                   double minEntries,
                   double yAxisMax)
{
  if (!h2 || !outFile) return;

  const std::string cleanName = CleanName(h2->GetName());
  const std::vector<double> edges = MakeGroupedEdges(h2->GetXaxis(), xBinGroup);
  const int nOutBins = static_cast<int>(edges.size()) - 1;

  TH1D *hIterativeSigma = new TH1D(("h_iterative_sigma_" + cleanName).c_str(),
                                   Form("%s;%.120s;#sigma of %.120s",
                                        h2->GetTitle(),
                                        h2->GetXaxis()->GetTitle(),
                                        h2->GetYaxis()->GetTitle()),
                                   nOutBins, edges.data());
  TH1D *hFitSlicesSigma = BuildRootFitSlicesYSigma(h2, cleanName, edges, xBinGroup, minEntries);
  TH1D *hFwhmSigma = new TH1D(("h_fwhm_sigma_" + cleanName).c_str(),
                              Form("%s;%.120s;FWHM / 2.355 of %.120s",
                                   h2->GetTitle(),
                                   h2->GetXaxis()->GetTitle(),
                                   h2->GetYaxis()->GetTitle()),
                              nOutBins, edges.data());
  TH1D *hFwhm = new TH1D(("h_fwhm_" + cleanName).c_str(),
                         Form("%s;%.120s;FWHM of %.120s",
                              h2->GetTitle(),
                              h2->GetXaxis()->GetTitle(),
                              h2->GetYaxis()->GetTitle()),
                         nOutBins, edges.data());
  TH1D *hInputGaussianWidth = BuildInputGaussianWidthHist(
      FindInputGaussianWidth(inputDir, h2->GetName()),
      cleanName, edges, h2->GetXaxis()->GetTitle(), h2->GetYaxis()->GetTitle());
  TH1D *hMean = new TH1D(("h_mean_" + cleanName).c_str(),
                         Form("%s;%.120s;Core Gaussian mean of %.120s",
                              h2->GetTitle(),
                              h2->GetXaxis()->GetTitle(),
                              h2->GetYaxis()->GetTitle()),
                         nOutBins, edges.data());
  TH1D *hChi2Ndf = new TH1D(("h_chi2ndf_" + cleanName).c_str(),
                            Form("%s;%.120s;#chi^{2}/NDF",
                                 h2->GetTitle(), h2->GetXaxis()->GetTitle()),
                            nOutBins, edges.data());
  TH1D *hEntries = new TH1D(("h_entries_" + cleanName).c_str(),
                            Form("%s;%.120s;Slice entries",
                                 h2->GetTitle(), h2->GetXaxis()->GetTitle()),
                            nOutBins, edges.data());
  TH1D *hIterations = new TH1D(("h_iterations_" + cleanName).c_str(),
                               Form("%s;%.120s;Fit iterations",
                                    h2->GetTitle(), h2->GetXaxis()->GetTitle()),
                               nOutBins, edges.data());
  TH1D *hFitStatus = new TH1D(("h_fit_status_" + cleanName).c_str(),
                              Form("%s;%.120s;Fit status",
                                   h2->GetTitle(), h2->GetXaxis()->GetTitle()),
                              nOutBins, edges.data());

  std::ofstream csv(Form("%s/%s_resolution_comparison_vs_x.csv", plotDir.c_str(), cleanName.c_str()));
  csv << "x_low,x_high,x_center,entries,"
      << "iterative_mean,iterative_mean_error,iterative_sigma,iterative_sigma_error,"
      << "iterative_chi2_ndf,iterative_iterations,iterative_fit_status,"
      << "fitslicesy_sigma,fitslicesy_sigma_error,fwhm,fwhm_sigma,"
      << "input_gaussian_width,input_gaussian_width_error\n";
  csv << std::setprecision(12);

  int outBin = 1;
  const int nX = h2->GetXaxis()->GetNbins();
  for (int firstXBin = 1; firstXBin <= nX; firstXBin += xBinGroup, ++outBin) {
    const int lastXBin = std::min(nX, firstXBin + xBinGroup - 1);
    TH1D *slice = h2->ProjectionY(Form("slice_%s_xbin_%d_%d",
                                       cleanName.c_str(), firstXBin, lastXBin),
                                  firstXBin, lastXBin, "e");
    slice->SetDirectory(nullptr);

    const SliceFitResult fit = FitSliceIterative(slice, nSigma, relSigmaTolerance,
                                                 maxIterations, minEntries);
    const FwhmResult fwhm = CalculateFwhmSigma(slice, minEntries);
    const double xLow = h2->GetXaxis()->GetBinLowEdge(firstXBin);
    const double xHigh = h2->GetXaxis()->GetBinUpEdge(lastXBin);
    const double xCenter = 0.5 * (xLow + xHigh);

    hEntries->SetBinContent(outBin, fit.entries);
    hIterations->SetBinContent(outBin, fit.iterations);
    hFitStatus->SetBinContent(outBin, fit.fitStatus);

    if (fit.ok) {
      hIterativeSigma->SetBinContent(outBin, fit.sigma);
      hIterativeSigma->SetBinError(outBin, fit.sigmaErr);
      hMean->SetBinContent(outBin, fit.mean);
      hMean->SetBinError(outBin, fit.meanErr);
      hChi2Ndf->SetBinContent(outBin, fit.chi2Ndf);
    }

    if (fwhm.ok) {
      hFwhm->SetBinContent(outBin, fwhm.fwhm);
      hFwhm->SetBinError(outBin, 0.0);
      hFwhmSigma->SetBinContent(outBin, fwhm.sigma);
      hFwhmSigma->SetBinError(outBin, 0.0);
    }

    csv << xLow << ',' << xHigh << ',' << xCenter << ','
        << fit.entries << ',' << fit.mean << ',' << fit.meanErr << ','
        << fit.sigma << ',' << fit.sigmaErr << ',' << fit.chi2Ndf << ','
        << fit.iterations << ',' << fit.fitStatus << ','
        << hFitSlicesSigma->GetBinContent(outBin) << ','
        << hFitSlicesSigma->GetBinError(outBin) << ','
        << fwhm.fwhm << ',' << fwhm.sigma << ','
        << hInputGaussianWidth->GetBinContent(outBin) << ','
        << hInputGaussianWidth->GetBinError(outBin) << '\n';

    delete slice;
  }

  outFile->cd();
  TDirectory *dir = outFile->mkdir(cleanName.c_str());
  if (!dir) dir = outFile;
  dir->cd();
  hIterativeSigma->Write("resolution_sigma");
  hIterativeSigma->Write("iterative_sigma");
  hFitSlicesSigma->Write("fitslicesy_sigma");
  hFwhmSigma->Write("fwhm_equivalent_sigma");
  hFwhm->Write("fwhm");
  hInputGaussianWidth->Write("input_gaussian_width");
  hMean->Write("gaussian_mean");
  hChi2Ndf->Write("chi2_ndf");
  hEntries->Write("slice_entries");
  hIterations->Write("fit_iterations");
  hFitStatus->Write("fit_status");
  DrawResolutionCanvas(hIterativeSigma, hFitSlicesSigma, hFwhmSigma, hInputGaussianWidth,
                       cleanName, plotDir, yAxisMax);

  delete hIterativeSigma;
  delete hFitSlicesSigma;
  delete hFwhmSigma;
  delete hFwhm;
  delete hInputGaussianWidth;
  delete hMean;
  delete hChi2Ndf;
  delete hEntries;
  delete hIterations;
  delete hFitStatus;
}

}  // namespace

void ListTH2InFile(const char *inputFile)
{
  TFile inFile(inputFile, "READ");
  if (inFile.IsZombie()) {
    std::cerr << "Could not open " << inputFile << std::endl;
    return;
  }

  std::vector<std::string> names;
  CollectTH2(&inFile, names);
  std::sort(names.begin(), names.end());

  std::cout << "TH2 histograms in " << inputFile << ":\n";
  for (const auto &name : names) {
    std::cout << "  " << name << '\n';
  }
}

void ExtractResolutionVsX(const char *inputFile,
                          const char *histName = "",
                          const char *outDir = "resolution_vs_x",
                          int xBinGroup = 1,
                          double nSigma = 2.0,
                          double relSigmaTolerance = 0.005,
                          int maxIterations = 10,
                          double minEntries = 50.0,
                          double yAxisMax = -1.0)
{
  TH1::AddDirectory(kFALSE);
  gSystem->mkdir(outDir, kTRUE);

  TFile inFile(inputFile, "READ");
  if (inFile.IsZombie()) {
    std::cerr << "Could not open " << inputFile << std::endl;
    return;
  }

  const std::string requested = histName ? histName : "";
  if (requested == "LIST" || requested == "list") {
    ListTH2InFile(inputFile);
    return;
  }

  std::vector<std::string> histNames;
  if (requested.empty()) {
    CollectTH2(&inFile, histNames);
    std::sort(histNames.begin(), histNames.end());
  } else {
    histNames = SplitRequestedHistograms(requested);
  }

  if (histNames.empty()) {
    std::cerr << "No TH2 histograms found in " << inputFile << std::endl;
    return;
  }

  const std::string outRootName = std::string(outDir) + "/" + StripRootExtension(inputFile) +
                                  "_resolution_methods.root";
  const std::string plotDir = std::string(outDir) + "/" + StripRootExtension(inputFile);
  gSystem->mkdir(plotDir.c_str(), kTRUE);

  TFile outFile(outRootName.c_str(), "RECREATE");
  if (outFile.IsZombie()) {
    std::cerr << "Could not create " << outRootName << std::endl;
    return;
  }

  std::cout << "Writing resolution-method output to " << outRootName << std::endl;
  std::cout << "Fit settings: xBinGroup=" << xBinGroup
            << ", nSigma=" << nSigma
            << ", relSigmaTolerance=" << relSigmaTolerance
            << ", maxIterations=" << maxIterations
            << ", minEntries=" << minEntries
            << ", yAxisMax=" << yAxisMax << std::endl;

  for (const auto &name : histNames) {
    TH2 *h2 = dynamic_cast<TH2 *>(inFile.Get(name.c_str()));
    if (!h2) {
      std::cerr << "Skipping " << name << ": not found or not a TH2" << std::endl;
      continue;
    }

    std::cout << "Processing " << name << std::endl;
    ProcessOneTH2(h2, &inFile, &outFile, plotDir, std::max(1, xBinGroup), nSigma,
                  relSigmaTolerance, maxIterations, minEntries, yAxisMax);
  }

  outFile.Close();
  std::cout << "Done." << std::endl;
}
