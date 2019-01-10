$PARAM @annotated
  CL   :  2 : Clearance (volume/time)
  Vc    : 20 : Central volume (volume)
  KA1   : 1 : Absorption rate constant 1 (1/time)
  KA2   : 1 : Absorption rate constant 2 (1/time)
  F1:    0.8 : bioavailability 1
  F2:    0.2 : bioavailability 2
  ALAG1  : 0 : Lag time depot 1 (time)
  ALAG2  : 0 : Lag time depot 2 (time)
  Q    :  2 : Inter-compartmental clearance (volume/time)
  Vp   : 10 : Peripheral volume of distribution (volume)
  
  $CMT  @annotated
  EV1    : First extravascular compartment
  EV2    : Second extravascular compartment
  CENT   : Central compartment
  PER   : Peripheral compartment
  
  
  $MAIN
  double CLi = exp(log(CL) + ETA_CL);
  double V3i = exp(log(Vc)  + ETA_Vc);
  double KA1i = exp(log(KA1)  + ETA_KA1);
  double KA2i = exp(log(KA2)  + ETA_KA2);
  double Qi = exp(log(Q)  + ETA_Q);
  double V4i = exp(log(Vp)  + ETA_Vp);

ALAG_EV1 = ALAG1;
F_EV1 = F1;
ALAG_EV2 = ALAG2;
F_EV2 = F2;

$OMEGA @labels ETA_CL ETA_Vc ETA_KA1 ETA_KA2 ETA_Q ETA_Vp
  0 0 0 0 0 0

$GLOBAL
#define CP (CENT/V3i)
  
  $ODE
  dxdt_EV1  = -KA1i*EV1;
  dxdt_EV2  = -KA2i*EV2;
  dxdt_CENT = KA1i*EV1 + KA2i*EV2 - (CLi/V3i)*CENT + (Qi/V4i)*PER - (Qi/V3i)*CENT;
  dxdt_PER  = - (Qi/V4i)*PER + (Qi/V3i)*CENT;

$CAPTURE @annotated
  CP : Plasma concentration (mass/volume)