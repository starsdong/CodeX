void SampleFit(){

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<"Output to compare to that of Section 3.9.5"<<endl;
cout<<endl;

TTMParticleSet set("/home/wheaton/THERMUSPackage/THERMUS/particles/PartList_PPB2002.txt");
set.InputDecays("/home/wheaton/THERMUSPackage/THERMUS/particles/",true);

TTMParameterSetBSQ par(0.160,0.05,0.,0.,1.);
par.FitT(0.160);
par.FitMuB(0.05);
par.FitMuS(0.);

TTMThermalFitBSQ fit(&set,&par,"ExpData.txt");
fit.SetQStats(kFALSE);
fit.SetWidth(kFALSE);
fit.GenerateYields();
fit.ListYields();

fit.FitData(0);
par.List();

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<endl;

}
