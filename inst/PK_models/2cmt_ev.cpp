$PARAM @annotated
  TVCL   :  2 : Clearance (volume/time)
  TVV    : 20 : Central volume (volume)
  TVKA   : 1 : Absorption rate constant (1/time)
  Q    :  2 : Inter-compartmental clearance (volume/time)
  Vp   : 10 : Peripheral volume of distribution (volume)
  F: 1 : bioavailability
  ALAG  : 0 : Lag time (time)
  
  $CMT  @annotated
  EV   : Extravascular compartment
  CENT : Central compartment
  PERIPH : Peripheral compartment
  
  $MAIN
  double CL = exp(log(TVCL) + ETA_CL);
  double V2 = exp(log(TVV)  + ETA_V);
  double KA = exp(log(TVKA)  + ETA_KA);
  double V3 = exp(log(Vp)  + ETA_Vp);

ALAG_EV = ALAG;
F_EV = F;

$OMEGA @labels ETA_CL ETA_V ETA_KA ETA_Vp
  0 0 0 0

$GLOBAL
#define CP (CENT/V2)
  
  $PKMODEL ncmt = 2, depot = TRUE
  
  $CAPTURE @annotated
  CP : Plasma concentration (mass/volume)