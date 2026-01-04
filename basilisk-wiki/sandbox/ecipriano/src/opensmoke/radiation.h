double opensmoke_optically_thin (void * p) {
  ThermoState * ts = (ThermoState *)p;
  double T = ts->T, P = ts->P, * x = ts->x;

  OpenSMOKE_GasProp_SetTemperature (T);
  OpenSMOKE_GasProp_SetPressure (P);
  double kPlanckMix = OpenSMOKE_GasProp_kPlanckMix (x);

  return -4.*5.669e-8*kPlanckMix*(pow (T, 4.) - pow (300, 4.));
}
