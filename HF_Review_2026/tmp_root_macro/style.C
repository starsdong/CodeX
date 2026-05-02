//
//  Style file
// version 0.1
//

#include <iostream>
#include "TROOT.h"

//===========================
TStyle* RHIC_style() 
{
  TStyle *style = new TStyle("RHIC_style","RHIC_style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  style->SetFrameBorderMode(icol);
  style->SetFrameFillColor(icol);
  style->SetCanvasBorderMode(icol);
  style->SetCanvasColor(icol);
  style->SetPadBorderMode(icol);
  style->SetPadColor(icol);
  style->SetStatColor(icol);

  // set margin sizes
  style->SetPadBottomMargin(0.12);
  style->SetPadTopMargin(0.05);
  style->SetPadLeftMargin(0.12);
  style->SetPadRightMargin(0.05);

  // set title offsets (for axis label)
  style->SetTitleXOffset(1.1);
  style->SetTitleYOffset(1.1);
  style->SetTitleOffset(1.1,"z"); // Set the offset for Z axis titles expliticly to avoid it being cut off

  // use large fonts
  //Int_t font=72; // Helvetica italics
  const Int_t g_font = 42; // Helvetica
  const Double_t g_tsize = 0.045;
  const Double_t g_tsize_label = 0.042;
  style->SetTextFont(g_font);
  style->SetTextSize(g_tsize);
  
  style->SetLabelFont(g_font,"x");
  style->SetTitleFont(g_font,"x");
  style->SetLabelFont(g_font,"y");
  style->SetTitleFont(g_font,"y");
  style->SetLabelFont(g_font,"z");
  style->SetTitleFont(g_font,"z");
  
  style->SetLabelSize(g_tsize_label,"x");
  style->SetTitleSize(g_tsize,"x");
  style->SetLabelSize(g_tsize_label,"y");
  style->SetTitleSize(g_tsize,"y");
  style->SetLabelSize(g_tsize_label,"z");
  style->SetTitleSize(g_tsize,"z");

  // use bold lines and markers
  style->SetMarkerStyle(20);
  style->SetMarkerSize(1.2);
  style->SetHistLineWidth(2.);
  style->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of error bar caps
  style->SetEndErrorSize(0.);

  // statistics
  style->SetOptStat(0);
  style->SetOptFit(1);
  style->SetOptDate(0);
  style->SetStatFontSize(g_tsize);

  // get rid of grid
  style->SetPadGridX(0);
  style->SetPadGridY(0);

  // legend modification
  style->SetLegendBorderSize(0);
  style->SetLegendFillColor(0);
  style->SetLegendFont(g_font);

#if ROOT_VERSION_CODE >= ROOT_VERSION(6,00,0)
  std::cout << "style: ROOT6 mode" << std::endl;
  style->SetLegendTextSize(g_tsize);
  style->SetPalette(kBird);
#else
  std::cout << "style: ROOT5 mode" << std::endl;
  // color palette - manually define 'kBird' palette only available in ROOT 6
  Int_t alpha = 0;
  Double_t stops[9] = { 0.0000, 0.1250, 0.2500, 0.3750, 0.5000, 0.6250, 0.7500, 0.8750, 1.0000};
  Double_t red[9]   = { 0.2082, 0.0592, 0.0780, 0.0232, 0.1802, 0.5301, 0.8186, 0.9956, 0.9764};
  Double_t green[9] = { 0.1664, 0.3599, 0.5041, 0.6419, 0.7178, 0.7492, 0.7328, 0.7862, 0.9832};
  Double_t blue[9]  = { 0.5293, 0.8684, 0.8385, 0.7914, 0.6425, 0.4662, 0.3499, 0.1968, 0.0539};
  TColor::CreateGradientColorTable(9, stops, red, green, blue, 255, alpha);
#endif

  style->SetNumberContours(80);

  return style;

}

//===========================
void set_style()
{
  static TStyle* style = 0;
  std::cout << "Style: Applying nominal settings." << std::endl ;
  if ( style==0 ) style = RHIC_style();
  gROOT->SetStyle("RHIC_style");
  gROOT->ForceStyle();
}
