$PARAM @annotated
  CL   :  1 : Clearance (volume/time)
  V    : 20 : Central volume (volume)
  KA   : 1 : Absorption rate constant (1/time)
  Q1    :  2 : First inter-compartmental clearance (volume/time)
  VP1   : 10 : First peripheral volume of distribution (volume)
  Q2   :   2 : Second inter-compartmental clearance (volume/time)
  VP2  : 100 : Second peripheral volume (volume) 
  
  $CMT  @annotated
  EV   : Extravascular compartment
  CENT : Central compartment
  PERIPH : First Peripheral compartment
  PERIPH2 : Second peripheral compartment
  
  $MAIN
  double CL1 = exp(log(CL) + ETA_CL);
double V1 = exp(log(V)  + ETA_V);
double KA1 = exp(log(KA)  + ETA_KA);


$OMEGA @labels ETA_CL ETA_V ETA_KA
  0 0 0


$GLOBAL
#define CP (CENT/V1)
#define CT (PERIPH/VP1)
#define CT2 (PERIPH2/VP2)
  
  $ODE
  dxdt_EV= -KA1*EV;
  dxdt_CENT = KA1*EV - (CL1+Q1+Q2)*CP  + Q1*CT + Q2*CT2;
dxdt_PERIPH = Q1*CP - Q1*CT;
dxdt_PERIPH2 = Q2*CP - Q2*CT2;

$CAPTURE @annotated
  CP : Plasma concentration (mass/volume)
