$PARAM @annotated
  TVCL   :  1 : Clearance (volume/time)
  TVV    : 20 : Central volume (volume)
  Q    :  2 : Inter-compartmental clearance (volume/time)
  V2   : 10 : Peripheral volume of distribution (volume)
  
  $CMT  @annotated
  CENT : Central compartment
  PERIPH : Peripheral compartment
  
  $MAIN
  double CL = exp(log(TVCL) + ETA_CL);
  double V1 = exp(log(TVV)  + ETA_V);


$OMEGA @labels ETA_CL ETA_V
  0 0


$GLOBAL
#define CP (CENT/V1)
  
  $PKMODEL ncmt = 2, depot = FALSE
  
  $CAPTURE @annotated
  CP : Plasma concentration (mass/volume)