void SampleMacro(){

// Section 3.3.2 //

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<"Output to compare to that of Section 3.3.2"<<endl;
cout<<endl;

TTMParticleSet set("/home/wheaton/THERMUSPackage/THERMUS/particles/PartList_PPB2002.txt");
set.InputDecays("/home/wheaton/THERMUSPackage/THERMUS/particles/",true);

TTMParticle *part = set.GetParticle(32114);
part->List();

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<endl;

// Section 3.5.4 //

cout<<"#####################################################################################"<<endl;
cout<<"Output to compare to that of Section 3.5.4"<<endl;
cout<<endl;

TTMParameterSetBQ parBQ(0.160,0.2,-0.01,0.8,6.,6.);
parBQ.FitT(0.160);
parBQ.FitMuB(0.2);
parBQ.ConstrainMuQ(1.2683);
parBQ.List();

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<endl;

// Section 3.6.4 //

cout<<"#####################################################################################"<<endl;
cout<<"Output to compare to that of Section 3.6.4"<<endl;
cout<<endl;

TTMThermalParticleBQ therm_delta(part,&parBQ,1.);
cout<<"n_Boltz (No Width) = "<<therm_delta.DensityBoltzmannNoWidth()<<endl;
cout<<"e_Boltzmann (Width) = "<<therm_delta.EnergyBoltzmannWidth()<<endl;

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<endl;

// Section 3.7.9 //

cout<<"#####################################################################################"<<endl;
cout<<"Output to compare to that of Section 3.7.9"<<endl;
cout<<endl;

TTMThermalModelBQ modBQ(&set,&parBQ);
modBQ.GenerateParticleDens();
parBQ.List();

cout<<endl;

modBQ.GenerateEnergyDens();
modBQ.GenerateEntropyDens();
modBQ.GeneratePressure();

TTMDensObj *delta_dens = modBQ.GetDensities(32114);
delta_dens->List();

cout<<endl;

TTMDensObj *piplus_dens = modBQ.GetDensities(211);
piplus_dens->List();

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<endl;

// Section 3.9.5 //

cout<<"#####################################################################################"<<endl;
cout<<"Output to compare to that of Section 3.9.5"<<endl;
cout<<endl;

TTMParameterSetBSQ par(0.160,0.05,0.,0.,1.);
par.FitT(0.160);
par.FitMuB(0.05);
par.FitMuS(0.);

TTMThermalFitBSQ fit(&set,&par,"ExpData.txt");
fit.SetQStats(kFALSE);
fit.SetWidth(kFALSE);
fit.GenerateYields();
fit.ListYields();

cout<<endl;
cout<<"#####################################################################################"<<endl;
cout<<endl;

}
