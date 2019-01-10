$PARAM @annotated
  TVCL   :  2 : Clearance (volume/time)
  TVV    : 20 : Central volume (volume)
  TVKA   : 1 : Absorption rate constant (1/time)
  F: 1 : bioavailability
  ALAG  : 0 : Lag time (time)
  
  $CMT  @annotated
  EV   : Extravascular compartment
  CENT : Central compartment
  
  $MAIN
  double CL = exp(log(TVCL) + ETA_CL);
  double V = exp(log(TVV)  + ETA_V);
  double KA = exp(log(TVKA)  + ETA_KA);

ALAG_EV = ALAG;
F_EV = F;

$OMEGA @labels ETA_CL ETA_V ETA_KA
  0 0 0

$GLOBAL
#define CP (CENT/V)
  
  $PKMODEL ncmt = 1, depot = TRUE
  
  $CAPTURE @annotated
  CP : Plasma concentration (mass/volume)