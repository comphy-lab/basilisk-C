int NCCLASS;

void discrete_blues (double cmap[NCMAP][3])
{
  double startColor[3] = {0.,0.,0.5};
  double endColor[3] = {0.95,0.95,1.};
  double classWidth = (-1.+(double)NCMAP)/((double)NCCLASS);

  for (int i = 0; i < NCMAP; i++) {
    double f = i<(NCMAP-1) ? floor( ((double) i) / classWidth ) / (-1+(double)NCCLASS) : 1.0;
    cmap[i][0] = f * endColor[0] + (1.-f) * startColor[0];
    cmap[i][1] = f * endColor[1] + (1.-f) * startColor[1];
    cmap[i][2] = f * endColor[2] + (1.-f) * startColor[2];
  }
}

void discrete_reds (double cmap[NCMAP][3])
{
  double startColor[3] = {0.45,0.,0.};
  double endColor[3] = {1.,0.9,0.9};
  double classWidth = (-1.+(double)NCMAP)/((double)NCCLASS);

  for (int i = 0; i < NCMAP; i++) {
    double f = i<(NCMAP-1) ? floor( ((double) i) / classWidth ) / (-1+(double)NCCLASS) : 1.0;
    cmap[i][0] = f * endColor[0] + (1.-f) * startColor[0];
    cmap[i][1] = f * endColor[1] + (1.-f) * startColor[1];
    cmap[i][2] = f * endColor[2] + (1.-f) * startColor[2];
  }
}

void discrete_grays (double cmap[NCMAP][3])
{
  double startColor[3] = {0.2,0.2,0.2};
  double endColor[3] = {0.9,0.9,0.9};
  double classWidth = (-1.+(double)NCMAP)/((double)NCCLASS);

  for (int i = 0; i < NCMAP; i++) {
    double f = i<(NCMAP-1) ? floor( ((double) i) / classWidth ) / (-1+(double)NCCLASS) : 1.0;
    cmap[i][0] = f * endColor[0] + (1.-f) * startColor[0];
    cmap[i][1] = f * endColor[1] + (1.-f) * startColor[1];
    cmap[i][2] = f * endColor[2] + (1.-f) * startColor[2];
  }
}

void continuous_blues (double cmap[NCMAP][3])
{
  double startColor[3] = {0.,0.,0.5};
  double endColor[3] = {0.95,0.95,1.};

  for (int i = 0; i < NCMAP; i++) {
    double f = ((double) i) / (-1.+(double)NCMAP);
    cmap[i][0] = f * endColor[0] + (1.-f) * startColor[0];
    cmap[i][1] = f * endColor[1] + (1.-f) * startColor[1];
    cmap[i][2] = f * endColor[2] + (1.-f) * startColor[2];
  }
}

void continuous_reds (double cmap[NCMAP][3])
{
  double startColor[3] = {0.45,0.,0.};
  double endColor[3] = {1.,0.9,0.9};

  for (int i = 0; i < NCMAP; i++) {
    double f = ((double) i) / (-1.+(double)NCMAP);
    cmap[i][0] = f * endColor[0] + (1.-f) * startColor[0];
    cmap[i][1] = f * endColor[1] + (1.-f) * startColor[1];
    cmap[i][2] = f * endColor[2] + (1.-f) * startColor[2];
  }
}

void continuous_grays (double cmap[NCMAP][3])
{
  double startColor[3] = {0.2,0.2,0.2};
  double endColor[3] = {0.9,0.9,0.9};

  for (int i = 0; i < NCMAP; i++) {
    double f = ((double) i) / (-1.+(double)NCMAP);
    cmap[i][0] = f * endColor[0] + (1.-f) * startColor[0];
    cmap[i][1] = f * endColor[1] + (1.-f) * startColor[1];
    cmap[i][2] = f * endColor[2] + (1.-f) * startColor[2];
  }
}

void topocolormap (double cmap[NCMAP][3])
{
   double classColors[5][3] = {{0.,0.,0.6},
			     {0.6,0.6,1.},
                             {1.,1.,1.},
                             {1.,0.6,0.6},
                             {0.6,0.,0.}};

  double classWidth = (-1.+(double)NCMAP)/5.;

  for (int i = 0; i < NCMAP; i++) {
    int nc = i<(NCMAP-1) ? (int) floor( ((double) i) / classWidth ) : 5;
    cmap[i][0] = classColors[nc][0];
    cmap[i][1] = classColors[nc][1];
    cmap[i][2] = classColors[nc][2];
  }
   
}

void uniform_gray(double cmap[NCMAP][3])
{
  for (int i = 0; i < NCMAP; i++) {
    cmap[i][0] = 0.7;
    cmap[i][1] = 0.7;
    cmap[i][2] = 0.7;
  }
}
