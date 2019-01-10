$PARAM @annotated
  TVCL   :  2 : Clearance (volume/time)
  TVV    : 20 : Central volume (volume)
  KA1   : 1 : Absorption rate constant 1 (1/time)
  KA2   : 1 : Absorption rate constant 2 (1/time)
  F1: 0.8 : bioavailability 1
  F2: 0.2 : bioavailability 2
  ALAG1  : 0 : Lag time depot 1 (time)
  ALAG2  : 0 : Lag time depot 2 (time)
  
  $CMT  @annotated
  EV1    : First extravascular compartment
  EV2    : Second extravascular compartment
  CENT   : Central compartment

  
  $MAIN
  double CL = exp(log(TVCL) + ETA_CL);
  double V = exp(log(TVV)  + ETA_V);
  double KA1i = exp(log(KA1)  + ETA_KA1);
  double KA2i = exp(log(KA2)  + ETA_KA2);

ALAG_EV1 = ALAG1;
F_EV1 = F1;
ALAG_EV2 = ALAG2;
F_EV2 = F2;

$OMEGA @labels ETA_CL ETA_V ETA_KA1 ETA_KA2
  0 0 0 0

$GLOBAL
#define CP (CENT/V)

  $ODE
  dxdt_EV1 = -KA1i*EV1;
dxdt_EV2 = -KA2i*EV2;
dxdt_CENT = KA1i*EV1 + KA2i*EV2 - (CL)*CP;
  
  $CAPTURE @annotated
  CP : Plasma concentration (mass/volume)