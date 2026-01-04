PhaseMassBalance * liq_balance = NULL, * gas_balance = NULL;
FILE * fbalances = NULL;

#define FORMAT_S "%-15s"
#define FORMAT_G "%-15g"

event init (i = 0) {
  if (nboundary != 4)
    fprintf (ferr,
        "src/balances/two-phase.h:%d: warning: balances considers just 4 boundaries\n",
        LINENO), fflush (ferr);
  liq_balance = new_phase_mass_balance (liq);
  gas_balance = new_phase_mass_balance (gas);

  char name[80];
  sprintf (name, "balances-%d", grid->maxdepth);
  fbalances = fopen (name, "w");

  fprintf (fbalances, "time[s](1) totmass1(2) totmass2(3) evapmass(4) ");
  int counter = 5;
  foreach_species_in (liq)
    fprintf (fbalances, "mL_%s(%d) ", liq->species[i], counter++);
  foreach_species_in (gas)
    fprintf (fbalances, "mG_%s(%d) ", gas->species[i], counter++);
  foreach_species_in (gas)
    fprintf (fbalances, "mE_%s(%d) ", gas->species[i], counter++);
  fprintf (fbalances, "\n");
}

event cleanup (t = end) {
  delete_phase_mass_balances (liq_balance), liq_balance = NULL;
  delete_phase_mass_balances (gas_balance), gas_balance = NULL;
  fclose (fbalances), fbalances = NULL;
}

static void write_balances (void) {
  fprintf (fbalances, FORMAT_G FORMAT_G FORMAT_G FORMAT_G, t,
      liq_balance->mtot0 - (liq_balance->mtot + liq_balance->mftot),
      (gas_balance->mtot + gas_balance->mftot) - gas_balance->mtot0,
      gas_balance->mevaptot);
  foreach_species_in (liq)
    fprintf (fbalances, FORMAT_G,
        liq_balance->m0[i] - (liq_balance->m[i] + liq_balance->mf[i]));
  foreach_species_in (gas)
    fprintf (fbalances, FORMAT_G,
        (gas_balance->m[i] + gas_balance->mf[i]) - gas_balance->m0[i]);
  foreach_species_in (gas)
    fprintf (fbalances, FORMAT_G, gas_balance->mevap[i]);
  fprintf (fbalances, "\n");
}

event balances (i++) {
  face vector ufl = (nv > 1) ? uflist[1] : uflist[0];
  face vector ufg = uflist[0];
  phase_mass_balance (liq_balance, liq, NULL, dt, ufl, f, boundaries = true);
  phase_mass_balance (gas_balance, gas, mEvapList, dt, ufg, f, boundaries = true);

  if (pid() == 0)
    write_balances();
}

