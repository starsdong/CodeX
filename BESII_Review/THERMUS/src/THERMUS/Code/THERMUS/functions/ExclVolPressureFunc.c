// Author: Spencer Wheaton 14 July 2004 //

#include <TThermalModelBSQ.h>

extern TTMThermalModelBSQ *gModelBSQExcVol;

Float_t ExclVolPressureFunc(Float_t x)
{

  Double_t shift = (Double_t)x;

  Double_t pid = gModelBSQExcVol->ExclVolShiftedPressure(shift);

  return (Float_t)(pid - shift)/pid;

}
