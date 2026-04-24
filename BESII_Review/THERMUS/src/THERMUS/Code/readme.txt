
50 files in /main (.h and .cxx are C++ header and class implementation files respectively):
----------------------------------------------------------------------------------------

 makefile                                             :makefile to compile classes in \main and build shared object file                  
 THERMUSClassesLinkDef.h                              :required by ROOT to generate dictionary of classes

 TTMIDObj.h & TTMIDObj.cxx                            :basic class for storing particle ID                 

 TTMDecay.h & TTMDecay.cxx                            :class storing parent ID, daughter ID and branching ratio                  
 TTMDecayChannel.h & TTMDecayChannel.cxx              :class storing list of daughters and branching ratio of 
                                                       a decay channel         

 TTMParticle.h & TTMParticle.cxx                      :class storing particle properties relevant to thermal model

 TTMParticleSet.h & TTMParticleSet.cxx                :class grouping a collection of particles

 TTMParameter.h & TTMParameter.cxx                    :class storing details of a thermal parameter

 TTMParameterSet.h & TTMParameterSet.cxx              :base class for parameter set classes
 
 TTMParameterSetBSQ.h & TTMParameterSetBSQ.cxx        :derived parameter set class applicable to grand-canonical 
                                                       treatment of B, S and Q.
 
 TTMParameterSetBQ.h & TTMParameterSetBQ.cxx          :derived parameter set class applicable to the grand-canonical 
                                                       treatment of B and Q and the canonical treatment of S
  
 TTMParameterSetCanBSQ.h & TTMParameterSetCanBSQ.cxx  :derived parameter set class applicable to the canonical 
                                                       treatment of B, S and Q

 TTMDensObj.h & TTMDensObj.cxx                        :object containing densities for storage in a container class              

 TThermalParticle.h & TThermalParticle.cxx             :base class for thermal particle classes      
 
 TThermalParticleBSQ.h & TThermalParticleBSQ.cxx       :derived thermal particle class applicable to grand-canonical 
                                                        treatment of B, S and Q.

 TThermalParticleBQ.h & TThermalParticleBQ.cxx         :derived thermal particle class applicable to the grand-canonical 
                                                        treatment of B and Q and the canonical treatment of S

 TThermalParticleCanBSQ.h & TThermalParticleCanBSQ.cxx :derived thermal particle class applicable to the canonical 
                                                        treatment of B, S and Q

 TThermalModel.h & TThermalModel.cxx                   :base class for thermal model classes           
 
 TThermalModelBSQ.h & TThermalModelBSQ.cxx             :derived thermal model class applicable to grand-canonical 
                                                        treatment of B, S and Q.
 
 TThermalModelBQ.h & TThermalModelBQ.cxx               :derived thermal model class applicable to the grand-canonical 
                                                        treatment of B and Q and the canonical treatment of S

 TThermalModelCanBSQ.h & TThermalModelCanBSQ.cxx       :derived thermal model class applicable to the canonical 
                                                        treatment of B, S and Q

 TTMYield.h & TTMYield.cxx                             :object storing experimental yield or ratio of yields and relevant 
                                                        decay chains

 TThermalFit.h & TThermalFit.cxx                       :base class for thermal fit classes 

 TThermalFitBSQ.h & TThermalFitBSQ.cxx                 :derived thermal fit class applicable to grand-canonical 
                                                        treatment of B, S and Q.

 TThermalFitBQ.h & TThermalFitBQ.cxx                   :derived thermal fit class applicable to the grand-canonical 
                                                        treatment of B and Q and the canonical treatment of S

 TThermalFitCanBSQ.h & TThermalFitCanBSQ.cxx           :derived thermal model class applicable to the canonical 
                                                        treatment of B, S and Q


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

1 file in /nrc:
---------------
 makefile               :makefile to compile nrc functions added by user to /nrc and build shared object file

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

2 files in /include:
--------------------

 FncsConstrain.h        :contains declarations of all of the constraining functions used in THERMUS
 FncsThermalModel.h     :contains declarations of all the additional non-member functions used in THERMUS 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

52 files in /functions:
-----------------------

 makefile                                  :makefile to compile functions in \functions and build shared object file
 Minuit_function.c                         :function calculating the function to be minimised in a TMinuit fit
 fit_function.c                            :function initiating and running a fit of a thermal fit object using TMinuit
 Functions.c                               :non-member functions used by THERMUS
 Int_2_String.c                            :function converting integer to string
 IntegrationRoutines2D.c                   :2D integration routines used by THERMUS
 IntegrationRoutines.c                     :1D integration routines used by THERMUS
 ExclVolPressureFunc.c                     :function used in finding the excluded volume pressure
 FindExclVolPressure.c                     :function called to find excluded volume pressure (calls ExclVolPressureFunc)
 BSQCanonicalFncs.c                        :functions used in the BSQ canonical ensemble
 

 Constraining functions for the model treating B and Q grand-canonically and S canonically:
 ------------------------------------------------------------------------------------------

 BQfuncQBDens.c & BQConstrainQBDens.c      :constrains muQ and muB given nb + nb-bar                      
             
 BQfuncQEN.c & BQConstrainQEN.c            :constrains muQ and muB given E/N

 BQfuncQNetBDens.c & BQConstrainQNetBDens.c:constrains muQ and muB given nb - nb-bar

 BQfuncQPercolation.c & BQConstrainQPercolation.c :constrains muQ and muB given percolation model

 BQfuncST3.c & BQConstrainST3.c                   :constrains muB given S/T^3 

 BQfuncQST3.c & BQConstrainQST3.c                 :constrains muQ and muB given S/T^3    

 BQfuncBDens.c & BQConstrainBDens.c               :constrains muB given nb + nb-bar

 BQfuncEN.c & BQConstrainEN.c                     :constrains muB given E/N
           
 BQfuncQ.c & BQConstrainQ.c                       :constrains muQ     

 Constraining functions for the model treating B, S and Q grand-canonically:
 ---------------------------------------------------------------------------

 BSQfuncSQEN.c & BSQConstrainSQEN.c     :constrains muS, muQ and muB given E/N

 BSQfuncSQ.c & BSQConstrainSQ.c         :constrains muS and muQ

 BSQfuncSQBDens.c & BSQConstrainSQBDens.c :constrains muS, muQ and muB given nb + nb-bar

 BSQfuncSQNetBDens.c & BSQConstrainSQNetBDens.c :constrains muS, muQ and muB given nb-nb-bar

 BSQfuncSQST3.c & BSQConstrainSQST3.c           :constrains muS, muQ and muB given S/T^3
 
 BSQfuncSST3.c & BSQConstrainSST3.c             :constrains muS and muB given S/T^3

 BSQfuncSQPercolation.c & BSQConstrainSQPercolation.c :constrains muS, muQ and muB given percolation model 

 BSQfuncQ.c & BSQConstrainQ.c                         :constrains muQ 

 BSQfuncS.c & BSQConstrainS.c                         :constrains muS

 BSQfuncSPercolation.c & BSQConstrainSPercolation.c   :constrains muS and muB given percolation model

 BSQfuncSEN.c & BSQConstrainSEN.c                     :constrains muS and muB given E/N

 BSQfuncSBDens.c & BSQConstrainSBDens.c               :constrains muS and muB given nb + nb-bar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 File containing particle properties in \particles:
 --------------------------------------------------

 PartList_PPB2002.txt


 Particle decay files in \particles:
 -----------------------------------

 a0(1450)0_decay.txt
 a0(1450)+_decay.txt
 a0(980)0_decay.txt
 a0(980)+_decay.txt
 a1(1260)0_decay.txt
 a1(1260)+_decay.txt
 a2(1320)0_decay.txt
 a2(1320)+_decay.txt
 a4(2040)0_decay.txt
 a4(2040)+_decay.txt
 b1(1235)0_decay.txt
 b1(1235)+_decay.txt
 Delta(1232)0_decay.txt
 Delta(1232)-_decay.txt
 Delta(1232)+_decay.txt
 Delta(1232)++_decay.txt
 Delta(1600)0_decay.txt
 Delta(1600)-_decay.txt
 Delta(1600)+_decay.txt
 Delta(1600)++_decay.txt
 Delta(1620)0_decay.txt
 Delta(1620)-_decay.txt
 Delta(1620)+_decay.txt
 Delta(1620)++_decay.txt
 Delta(1700)0_decay.txt
 Delta(1700)-_decay.txt
 Delta(1700)+_decay.txt
 Delta(1700)++_decay.txt
 Delta(1905)0_decay.txt
 Delta(1905)-_decay.txt
 Delta(1905)+_decay.txt
 Delta(1905)++_decay.txt
 Delta(1910)0_decay.txt
 Delta(1910)-_decay.txt
 Delta(1910)+_decay.txt
 Delta(1910)++_decay.txt
 Delta(1920)0_decay.txt
 Delta(1920)-_decay.txt
 Delta(1920)+_decay.txt
 Delta(1920)++_decay.txt
 Delta(1930)0_decay.txt
 Delta(1930)-_decay.txt
 Delta(1930)+_decay.txt
 Delta(1930)++_decay.txt
 Delta(1950)0_decay.txt
 Delta(1950)-_decay.txt
 Delta(1950)+_decay.txt
 Delta(1950)++_decay.txt
 Delta(2420)0_decay.txt
 Delta(2420)-_decay.txt
 Delta(2420)+_decay.txt
 Delta(2420)++_decay.txt
 eta(1295)_decay.txt
 eta(1440)_decay.txt
 eta'(958)_decay.txt
 eta_decay.txt
 f0(1370)_decay.txt
 f0(1500)_decay.txt
 f0(1710)_decay.txt
 f0(600)_decay.txt
 f0(980)_decay.txt
 f1(1285)_decay.txt
 f1(1420)_decay.txt
 f2(1270)_decay.txt
 f2'(1525)_decay.txt
 f2(2010)_decay.txt
 f2(2300)_decay.txt
 f2(2340)_decay.txt
 f4(2050)_decay.txt
 h1(1170)_decay.txt
 K0*(1430)0_decay.txt
 K0*(1430)+_decay.txt
 K0_decay.txt
 K0S_decay.txt
 K1(1270)0_decay.txt
 K1(1270)+_decay.txt
 K1(1400)0_decay.txt
 K1(1400)+_decay.txt
 K*(1410)0_decay.txt
 K*(1410)+_decay.txt
 K*(1680)0_decay.txt
 K*(1680)+_decay.txt
 K2*(1430)0_decay.txt
 K2*(1430)+_decay.txt
 K2(1770)0_decay.txt
 K2(1770)+_decay.txt
 K2(1820)0_decay.txt
 K2(1820)+_decay.txt
 K3*(1780)0_decay.txt
 K3*(1780)+_decay.txt
 K4*(2045)0_decay.txt
 K4*(2045)+_decay.txt
 K*(892)0_decay.txt
 K*(892)+_decay.txt
 Ksi0_decay.txt
 Ksi(1530)0_decay.txt
 Ksi(1530)-_decay.txt
 Ksi(1690)0_decay.txt
 Ksi(1690)-_decay.txt
 Ksi(1820)0_decay.txt
 Ksi(1820)-_decay.txt
 Ksi(1950)0_decay.txt
 Ksi(1950)-_decay.txt
 Ksi(2030)0_decay.txt
 Ksi(2030)-_decay.txt
 Ksi-_decay.txt
 Lambda(1405)_decay.txt
 Lambda(1520)_decay.txt
 Lambda(1600)_decay.txt
 Lambda(1670)_decay.txt
 Lambda(1690)_decay.txt
 Lambda(1800)_decay.txt
 Lambda(1810)_decay.txt
 Lambda(1820)_decay.txt
 Lambda(1830)_decay.txt
 Lambda(1890)_decay.txt
 Lambda(2100)_decay.txt
 Lambda(2110)_decay.txt
 Lambda(2350)_decay.txt
 Lambda_decay.txt
 N(1440)0_decay.txt
 N(1440)+_decay.txt
 N(1520)0_decay.txt
 N(1520)+_decay.txt
 N(1535)0_decay.txt
 N(1535)+_decay.txt
 N(1650)0_decay.txt
 N(1650)+_decay.txt
 N(1675)0_decay.txt
 N(1675)+_decay.txt
 N(1680)0_decay.txt
 N(1680)+_decay.txt
 N(1700)0_decay.txt
 N(1700)+_decay.txt
 N(1710)0_decay.txt
 N(1710)+_decay.txt
 N(1720)0_decay.txt
 N(1720)+_decay.txt
 N(2190)0_decay.txt
 N(2190)+_decay.txt
 N(2220)0_decay.txt
 N(2220)+_decay.txt
 N(2250)0_decay.txt
 N(2250)+_decay.txt
 N(2600)0_decay.txt
 N(2600)+_decay.txt
 omega(1420)_decay.txt
 omega(1650)_decay.txt
 omega3(1670)_decay.txt
 omega(782)_decay.txt
 Omega_decay.txt
 phi(1020)_decay.txt
 phi(1680)_decay.txt
 phi3(1850)_decay.txt
 pi(1300)0_decay.txt
 pi(1300)+_decay.txt
 pi(1800)0_decay.txt
 pi(1800)+_decay.txt
 pi2(1670)0_decay.txt
 pi2(1670)+_decay.txt
 rho0_decay.txt
 rho(1450)0_decay.txt
 rho(1450)+_decay.txt
 rho(1700)0_decay.txt
 rho(1700)+_decay.txt
 rho3(1690)0_decay.txt
 rho3(1690)+_decay.txt
 rho+_decay.txt
 Sigma0_decay.txt
 Sigma(1385)0_decay.txt
 Sigma(1385)-_decay.txt
 Sigma(1385)+_decay.txt
 Sigma(1660)0_decay.txt
 Sigma(1660)-_decay.txt
 Sigma(1660)+_decay.txt
 Sigma(1670)0_decay.txt
 Sigma(1670)-_decay.txt
 Sigma(1670)+_decay.txt
 Sigma(1750)0_decay.txt
 Sigma(1750)-_decay.txt
 Sigma(1750)+_decay.txt
 Sigma(1775)0_decay.txt
 Sigma(1775)-_decay.txt
 Sigma(1775)+_decay.txt
 Sigma(1915)0_decay.txt
 Sigma(1915)-_decay.txt
 Sigma(1915)+_decay.txt
 Sigma(1940)0_decay.txt
 Sigma(1940)-_decay.txt
 Sigma(1940)+_decay.txt
 Sigma(2030)0_decay.txt
 Sigma(2030)-_decay.txt
 Sigma(2030)+_decay.txt
 Sigma(2250)0_decay.txt
 Sigma(2250)-_decay.txt  
 Sigma(2250)+_decay.txt  
 Sigma-_decay.txt
 Sigma+_decay.txt

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

0 files in /lib (directory is populated only once makefiles are run)
---------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
