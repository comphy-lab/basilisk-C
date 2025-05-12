
static double small_volume_removed = 0.0;

void remove_mdot_for_small_vols() {
  scalar f_pre[];
  foreach () { f_pre[] = f[]; }
#ifdef BUBBLE_VOL
  bool bubbles = true;
#else
  bool bubbles = false;
#endif

  remove_droplets(f, 3.0, 1.0e-3, bubbles);
  foreach (reduction(+:small_volume_removed)) {
    small_volume_removed += (f_pre[] - f[]) * dv();
  }

  foreach () {
    if (!(interfacial(point, f))) {
      m_dot[] = 0.0;
    }
  }
  return;
}
