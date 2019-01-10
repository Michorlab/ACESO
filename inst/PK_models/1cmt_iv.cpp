$PARAM @annotated
  TVCL   :  1 : Clearance (volume/time)
  TVV    : 20 : Central volume (volume)
  
  $CMT  @annotated
  CENT : Central compartment
  
  $MAIN
  double CL = exp(log(TVCL) + ETA_CL);
  double V = exp(log(TVV)  + ETA_V);


$OMEGA @labels ETA_CL ETA_V
  0 0


$GLOBAL
#define CP (CENT/V)
  
  $PKMODEL ncmt = 1, depot = FALSE
  
  $CAPTURE @annotated
  CP : Plasma concentration (mass/volume)