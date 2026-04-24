// Author: Spencer Wheaton 14 July 2004 //

#include <TROOT.h>
#include <TObject.h>
#include <TString.h>

TString Int_2_String(Int_t x)
{
char name[20];
sprintf(name,"%010d",x);
TString string(name);
return string;
}
