/**
# Computing the sun position
 
 
Focasting the positions of the Sun, Moon and planets has a long story.
This required long nights of observations recorded by the Babylonians (circa -1000).
More recently, during the Hellenistic period of Greece, the ancient greeks were able to create the first analog computer
 (the   [Antikythera](http://www.breves-de-maths.fr/trois-nombres-et-un-siecle-pour-decrypter-anticythere/)   [mechanism](http://dlib.nyu.edu/awdl/isaw/isaw-papers/4/)  circa -100).
 <br> <br>
 
 The next mile stone was [Newton](https://en.wikipedia.org/wiki/Isaac_Newton)'s law, fundamental in Mechanics, which is maybe a result of confinement:
 "Soon after Newton had obtained his BA degree in August 1665, the university temporarily closed as a precaution against the Great Plague. Although he had been undistinguished as a Cambridge student, Newton's private studies at his home in Woolsthorpe over the subsequent two years saw the development of his theories on calculus, optics, and the law of gravitation."
 Proof of Kepler's equation came later in 1679.<br><br>
 
 First computers were of course human computers (e.g. [bavarian gunners](https://link.springer.com/book/10.1007%2F978-3-642-02992-9) of [Petzal](https://www.transfert-films-dvd.com/lobjectif-photographique/) lens 1840; organized according to  a program, a "calculation plan"), or Marine officers for computation of "Droite de hauteur" with Sextant.
 But, the famous one were the "female computers" in the context of WWII and
  space conquest.
 On of these women is
 [Katherine Johnson](https://en.wikipedia.org/wiki/Katherine_Johnson) (deceased recently in feb 2020
 and popularized by the film "Hidden Figures").
 John Glenn  refused to fly unless Johnson
 verified the calculation by hand of his 1962 trajectory (computed by a
 FORTRAN 54 computer).<br><br>
 
Computation of the positions of planets was the computational challlenge in the early days of hand held computing. From hp RPN (starting from [hp 65](https://www.hpmuseum.org/hp65.htm) 1974),
 to BASIC (Meeus used  a [BASIC HP-85](http://articles.adsabs.harvard.edu//full/1986LAstr.100..571M/0000573.000.html),
he did a computation which took 470 hours and dissipated 120kWH, many books and programs have been written
around the 80' (Sérane prosed programs for popular BASIC computers in the late 80' etc).
  I guess that the formulas used by Meeus, Bouiges and others come from Danjon "Astronomie générale", but I have to check.
 The final tables of Meeus, have been coded in web pages  by many astromers
 ([here](http://xjubier.free.fr/site_pages/astronomy/ephemerides.html) or [there](http://f1rzv.free.fr/calculs/AstroCalc.php))
 googling the numerical values [33.231+%2B+13.17639653*N](https://www.google.fr/search?q=33.231+%2B+13.17639653*N)
 leeds to some page of [code](https://www.hpmuseum.org/forum/archive/index.php?thread-4262.html)...
 <br>
 <br>
 Of course the "Bureau des Longitudes" provides every year its "Ephémérides Nautiques"
 and  "[Ephémérides Astronomiques](https://www.imcce.fr/content/medias/publications/publications-institutionnelles/CDT_2020_ebook.pdf)"
 with rigorous position of all the planets. This is the reference.

 
 
 
 <br>
 <hr>
  <br> <br>
This program (in standard `C`) computes the position (cartesian and angular) of the Sun at any time (given UTC time, day, month and year) and any where in the word (for given latitude and longitude) using Bouiges 1981 and Meeus 1988,2014 simplified algorithms.
   An approximation for the Moon is given as well (to be continuated, only a first approximation).
 
 * The Kepler equation is solved. Some changes (rotations) of frameworks (ecliptic, equatorial...) are done to obtain final azimuth and hight of Sun and Moon.
 

 * This program can be used to find the true North using the azimuth of the sun for use with a [Bagnold sun-compas](https://www.iwm.org.uk/collections/item/object/30028094) (Bagnold and King 1931, Gross 2011).
 
 * This program can be used  to find the height of the sun for nautical celestial navigation with a sextant (Sacaze 1980, Caradec 1977).
 
 * As a result we plot a sundial (position of projection of
 the top of the gnomon on the soil, with analemma at noon)
 or "Cadran Solaire" (Camus & Gotteland 1993), we plot it with "brute force" without any simplification (as found in books on sundials, Opizzo  1990...).
 
 * We compute as well by "brute force" the equation of time (difference between time on a sundial and
 real time).
 
<br>
<br>

##  Tools
*/
 #include <stdio.h>
 #include <stdlib.h>
 #include <math.h>
 #include <string.h>
 #include <time.h>
/** Caution, all angles in degree (historical use), pay attention to Kepler Equation where there is a mix
*/
double sindeg(double theta){
    return sin(pi*theta/180);
}
double cosdeg(double theta){
    return cos(pi*theta/180);
}
double tandeg(double theta){
    return tan(pi*theta/180);
}
/**
 conversion from decimal to sexagesimal
 `HMS` function of HP calculator.
*/
char* HMS(double angle,char *symb){
    static char chaine[15];
    double s=1;
    if(angle<0) s=-1;
    angle = s*angle;
    double dd = (int)(angle);
    double mm = (int)(60*(angle - dd));
    double ss = ((angle - dd) - mm/60)*60*60;
    if (ss==60) {ss=0;mm=mm+1;}
    if (mm==60) {dd=dd+1;}
    sprintf (chaine,"%.0f%s%2.fmn%2.1fs", s*dd,symb,mm,ss);
    return chaine;
}
/** computation itself
 */
double hauteur,azimut;
int check;

void calc(double J, double M, double A, double heureTU, double lat, double longit){
/** number of days since  00/01/1901 ,
 from "SHOM Table des Marées des grands ports du Monde 1984"
see as well
 Bouiges p 32 and  "Date Algorithms" By Peter Baum
 )
     hours in UTC, latitude (positive North) and longitude (positive W) in decimal degrees
 */
    double Nj =
    floor(30.6001*(1 + M + 12*floor(0.7 + 1/(1 + M))))+
    floor(365.25*(A - floor(0.7 + 1/(1 + M))))
    -694403 + heureTU/24. + J ;
    ;
/**
## Kepler Equation for Sun position
 position of Sun, see Bouiges p 37:<br>
$L$ longitude moyenne hypothétique,
 note that $360/365.2425=0.985647342$, one year for a complete turn<br>
$\bar \omega$ longitude du périhélie.<br>
   with a given excentricity  $e=0.01675$<br>
 Values of Meeus are nearly the same, but he proposes quadratic corrections
 
 Cartesian position $x$ and $y$ of the Sun,
   $u$  anomalie
 excentrique et $v$ l'anomalie vraie  :
 ~~~gnuplot
 reset
 set size square
 set parametric
 set trange [0:2*pi]
 
 e=.8
 fx(t) = cos(t)
 fy(t) = sin(t)
 
 v(t)=arcos(cos(t)-e)/(1-e*cos(u))
 fxe(t) = cos(t)
 fye(t) = sqrt(1-e*e)*sin(t)
 
 a=1.
 set arrow 1 nohead from 0,0 to cos(a),sin(a)
 set arrow 2 nohead from cos(a),0 to cos(a),sin(a)
 set arrow 3 nohead from e,0 to cos(a),sqrt(1-e*e)*sin(a)
 set arrow 4 nohead from 0,0 to e,0
 set label 1 "u" at .1,.1
 set label 2 "v" at .05+e,.1
 set label 3 "O" at .05,-.1
 set label 4 "F" at .05+e,-.1
 set xrange[-1.5:1.5]
 set yrange[-1.5:1.5]
 plot fx(t),fy(t) not ,fxe(t),fye(t) t'Earth trajectory'
 
 unset arrow 1
 unset arrow 2
 unset arrow 3
 unset label 1
 unset label 3
 unset label 4
 
 reset
 
 ~~~
 
 
*/
    double L = 278.965 + 0.985647342*Nj;
    double omb = 281.235 + 0.0000469*Nj;
    double e = 0.01675;
/**
Kepler Equation (see details in Capderou 2012):
 $$u-e \sin(u) = M_t$$
 $M_t$ anomalie moyenne, c'est $M_t=L-\bar \omega$<br>
 $u$ est l'anomalie
 excentrique <br>
Note that as we are in degree for $u$ and $M_t$ there is an extra factor $\frac{180 }{\pi}$ :
 $$u-\frac{180 e}{\pi} \sin(u) = M_t$$
 Newton iteration to solve $K(u)=0$ with $K(u)=u-e sin(u) -M_t$ and the derivative (note that there is no $\frac{180 }{\pi}$ in it):
 $$K'(u)=\frac{dK}{du}= 1 - e \cos(u)$$
 Iteration, starting from 
  $u^0 = M_t$, 
  one or two is enough for the Sun as $e$ is small:
 $$u^{n+1}=u^{n} - \frac{K(u^n)}{K'(u^n)}$$
 */
    double Mt = fmod((L - omb),360);
    double u;
    u = Mt;
    u = u  - (u - e *180/pi*sindeg(u) - Mt)/(1 - e * cosdeg(u) );
    u = u  - (u - e *180/pi*sindeg(u) - Mt)/(1 - e * cosdeg(u) );
/**
## Coordonnées écliptiques géocentriques du SOleil

 Cartesian position $x$ and $y$ of the Sun (Coordonnées écliptiques géocentriques: $l,r$).<br>
 Relation entre $u$ l'anomalie
 excentrique et $v$ l'anomalie vraie et $r$ :<br>
 
 
 $r= a(1-e \cos(u))$ et   $x = a (\cos u -e)$ et $y= a \sqrt{1-e^2}\sin(u)$ <br>
 puis  $r=a (1-e^2)/(1+e \cos(v))$ et 
 $$\tan(\frac{v}{2})=\sqrt{\frac{1+e}{1-e}} \tan (\frac{u}{2})$$
 */
    double x =  (cosdeg(u) - e);
    double y = sqrt(1 - e*e)* sindeg(u);
    double v = atan2(y,x)*180/pi;
    double r = 0.999721/(1 + e*cosdeg(v));
    double l = fmod(v + omb, 360);
/**
## Coordonnées Astronomiques (coordonnées Equatoriales horaires)
 We now turn on the earth, wich axis of rotation is inclined by almost
 23.437463405244976 ° (Meeus)<br>
  the coordinates are "ascension droite"  $Asd$ and "déclinaison" $\delta$  */
    double atrop = 23.437463405244976;
    double delta = asin(sindeg(l)*sindeg(atrop))/pi*180;
    double Asd = fmod(24./360.*atan2(cosdeg(atrop)*sindeg(l),cosdeg(l))*180/pi, 24);
    double Asddeg = 360./24*Asd;
/**
 temps Sidéral  B p 84<br>
Local sideral time due to rotation of the Earth

every 23h56min04.1s the earth is aligned with same stars
  24 365.2425/(365.2425  + 1)  = "23h56m4.09073s"
 indeed <br>
$0.9856473415 1 + (23 + 56/60 + 4/3600)/24*360 \simeq 360$
 */
    double TSG360 = fmod(98.95958334 + 0.9856473415*(Nj) + heureTU/24*360, 360);
    double TSG24 = TSG360/360.*24 ;
    double TS360 = fmod(TSG360 - longit,360);
    double TS24 = TS360/360*24 ;
    /**  Angle Horaire à 00:00 TU */
    double AHv0 = fmod(TS360 - Asddeg, 360);
    

/**
## Azimutal Coordinates of Sun
 Changement de coordonnées Equatoriales en Azimutales page 87 (et 91) Bouiges:<br>
 $H$ angle horaire $TS-\alpha$ <br>
 $\delta$ déclinaison<br>
 $Az$ Azimut<br>
 $h$ hauteur<br>
 $\phi$ latitude du lieu d'observation
 $$\cos(\pi/2-h)=\sin(\phi)\sin(\delta) + \cos(\phi)  \cos(\delta) \cos(H)$$
    $$\sin(Az)=(\cos(\delta)\sin(H))/\sin(\pi/2-h)$$
    $$ \cos(Az)= ( -\cos(\phi)\sin(\delta) + \sin(\phi)  \cos(\delta) \cos(H))/(\sin(\pi/2-h))$$
 Calcul de la hauteur et Azimut (faire une procédure)
 */
    double haut = asin(
                  sindeg(lat)*sindeg(delta) +
                  cosdeg(lat)*cosdeg(delta)*cosdeg(TS360 - Asddeg))*180/pi;
    double yy1 = (-cosdeg(delta)* sindeg(TS360 - Asddeg));
    double xx1 = (
           cosdeg(lat)*sindeg(delta) - sindeg(lat)*cosdeg(delta)*cosdeg(TS360 - Asddeg) );
    double Az = fmod(360+atan2(yy1,xx1)*180/pi, 360);
    
    hauteur = haut;
    azimut = Az;
/**
## Check Sun
control and comparisons to check */

    if(check==1){
    fprintf(stdout, "#* * * * * * * * * * * *\n");
    fprintf(stdout, "#-----------------------\n");
    fprintf(stdout, "#The %2.f/%2.f/%4.f at UTC %s",J,M,A,HMS(heureTU,"h"));
    fprintf(stdout, " is Nj=%.2f Days since 00/01/1901 \n",Nj );
    fprintf(stdout, "#Temps Sidéral G = %s ",HMS(TSG360,"°"));
    fprintf(stdout, " Temps Sidéral Heures = %s",HMS(TSG24,"h"));
    fprintf(stdout, " Temps Sidéral local = %s\n",HMS(TS24,"h"));
    fprintf(stdout, "#-----------------------\n");
    fprintf(stdout, "# HP [EQ]  : r=%lf l=%lf H=%s",r,l,HMS(AHv0,"°"));
    fprintf(stdout, "  delta=%s \n",HMS(delta,"°"));
    fprintf(stdout, "# HP [E->A]: haut h=%s",HMS(haut,"°"));
    fprintf(stdout, " Azimut Z= %s\n",HMS(Az,"°"));
    fprintf(stdout, "#-----------------------\n");
    fprintf(stdout, "# Dec. val.: ");
    fprintf(stdout, "# Long=%.3f° Lat=0.000° haut=%.3f° Az=%.3f°  Asc drt=%.3fh  dec=%.3f° dist=%lf UA\n",
            l,haut,Az,Asd,delta,r);
    }
/**
## Moon
 
 Moon is rising soon!<br>
Bouiges page 46<br>
Meeus page 111 ...<br>
 
Many terms in the perturbation, others to be added soon!
 
To determine moonLongitude,  Bouiges uses 13 terms, Meeus 60
 */
    double Lp= 33.231 + 13.17639653*Nj;
    double Op=239.882 - 0.052953922*Nj;
    double Mp=18.294  + 13.06499245*Nj;
 
    double moonLong = Lp + 6.28875*sindeg(Mp)
                         + 0.2136*sindeg(2*Mp)
                         + 0.6583*sindeg(2*(Lp-L))
                         - 0.1856*sindeg(2*Mt)
                         + 1.2740*sindeg(2*(Lp-L)-Mp)
                         - 0.1143*sindeg(2*(Lp-Op))
                         + 0.00001;
    
    double moonLat = 5.1280*sindeg(Lp-Op)
                   + 0.2806*sindeg(Mp+Lp-Op)
                   + 0.2777*sindeg(Mp-Lp+Op)
                   + 0.1732*sindeg( 2*(Lp-L)-(Lp-Op))
                   + 0.00001;
/**
 
## Azimutal Coordinates of Moon
 
  Equatorial to Azimutal page 87 (and 91) Bouiges, or Meeus page 38 :<br>
 $H_{Moon}$   horar angle of Moon  <br>
 $\delta_{Moon}$ Moon's declinaison<br>
 $Az_{Moon}$ Azimuth of Moon<br>
 $h_{Moon}$ height of Moon<br>
 $\phi$ latitude
 $$\sin(h_{Moon})=\sin(\phi)\sin(\delta_{Moon}) + \cos(\phi)  \cos(\delta_{Moon}) \cos(H_{Moon})$$
 $$\sin(Az_{Moon})=  (-\cos(\delta_{Moon})\sin(H_{Moon}))/\sin(h_{Moon}-\pi/2)$$
 $$ \cos(Az_{Moon})= ( \cos(\phi)\sin(\delta_{Moon}) - \sin(\phi)  \cos(\delta_{Moon}) \cos(H_{Moon}))/(\sin(h_{Moon}-\pi/2))$$
  Compute Moon's height and azimuth (should do a procedure)
    */
    double moondelta=asin(sindeg(moonLat)*cosdeg(atrop) + cosdeg(moonLat)*sindeg(atrop)*sindeg(moonLong))/pi*180;
    double moonalpha = fmod(24./360.*atan2(cosdeg(atrop)*sindeg(moonLong)-tandeg(moonLat)*sindeg(atrop) ,cosdeg(moonLong))*180/pi, 24);
    double moonalphadeg=360./24*moonalpha;
    haut = asin(
                       sindeg(lat)*sindeg(moondelta) +
                       cosdeg(lat)*cosdeg(moondelta)*cosdeg(TS360 - moonalphadeg))*180/pi;
    yy1 = (-cosdeg(moondelta)* sindeg(TS360 - moonalphadeg));
    xx1 = (
                  cosdeg(lat)*sindeg(moondelta) - sindeg(lat)*cosdeg(moondelta)*cosdeg(TS360 - moonalphadeg) );
    Az = fmod(360+atan2(yy1,xx1)*180/pi, 360);
/**
## Check Moon
    control and comparisons to check */
     if(check==1){
         moonLong=fmod(moonLong,360);
         fprintf(stdout, "# Moon: Long=%.3f° lat=%.3f° --- alpha=%.3fh delta=%.3f° --- Z_moon=%.3f° h_moon =%.3f°  \n",
                 moonLong,moonLat,fmod(24+moonalpha,24),moondelta,Az,haut);
       
     }
    return;
}
/**
## Main
*/
int main() {
/**
check the values (comparison with Ephémérides etc) in file `log`
 */
    check=1;
    fprintf(stdout, " # check 26/04/2020 at UTC 16:00\n");
    calc(26,04,2020,16,48+48./60+45./3600,-(2+20./60+33./3600));
    fprintf(stdout, "\n\n # check 10/05/2020 at UTC 13:45\n");
    calc(10,05,2020,13.75,48+48./60+45./3600,-(2+20./60+33./3600));
    //calc(12,04,2011,0,0,0);
    //calc(10,05,1989,0,0,0);
    //  exit(0);
    check=0;
/**
 to be done: compute at present time
*/
    time_t timestamp = time(NULL);
    struct tm * now = localtime( & timestamp );
    int AAAA,MM,JJ,HH,mim;
    AAAA=1900+now->tm_year;
    MM=now->tm_mon+1;
    JJ=now->tm_mday;
    HH=now->tm_hour-2.; //summertime
    mim=now->tm_min;
    fprintf(stdout, "#-----------------------\n");
    fprintf(stdout, "# today is %2d/%2d/%4d at UTC=%2d:%2d\n",JJ,MM,AAAA,HH,mim);
    calc(JJ,MM,AAAA,HH+mim/60.,48+48./60+45./3600,-(2+20./60+33./3600));
   // fprintf(stdout, "# hauteur=%s ", HMS(hauteur,"°"));
   // fprintf(stdout, ", azimut=%s \n\n",HMS(azimut,"°"));
    fprintf(stdout, "# hauteur= %lf azimut=%lf\n", hauteur,azimut);

/**
 position of the shadow every month in 2020, rotation of the plate 43° (axis of Cité Verte), if the plate is vertical, other rotations are required.
 Mind the minus sign (clockwise rotation!)
*/
   
    double rot=43;
    for(double m=1;m<=12;m++)
    {
    for(double h=8;h<=16;h+=.25)
    {
     calc(21,m,2020,h,48+48./60+45./3600,-(2+20./60+33./3600));
     fprintf(stdout, " %lf, %lf %lf\n",fabs(1/tandeg(hauteur))*cosdeg(-azimut  + rot),
                   fabs(1/tandeg(hauteur))*sindeg(-azimut  + rot),m);
    }
         fprintf(stdout, "\n");
    }
/**
save hours for 3 representative months
*/
    FILE * fp = fopen ("hours.txt", "w");
    for(double h=8;h<=16;h+=1)
    {
        calc(21,12,2020,h,48+48./60+45./3600,-(2+20./60+33./3600));
        fprintf(fp, " %lf, %lf %2.2f \n",fabs(1/tandeg(hauteur))*cosdeg(-azimut  + rot),
                fabs(1/tandeg(hauteur))*sindeg(-azimut  + rot),h);
        calc(21,6,2020,h,48+48./60+45./3600,-(2+20./60+33./3600));
        fprintf(fp, " %lf, %lf %2.2f \n",fabs(1/tandeg(hauteur))*cosdeg(-azimut  + rot),
                fabs(1/tandeg(hauteur))*sindeg(-azimut  + rot),h);
        calc(20,3,2020,h,48+48./60+45./3600,-(2+20./60+33./3600));
        fprintf(fp, " %lf, %lf %2.2f \n",fabs(1/tandeg(hauteur))*cosdeg(-azimut  + rot),
                fabs(1/tandeg(hauteur))*sindeg(-azimut  + rot),h);
    }
    fclose (fp);
    fp = fopen ("hoursTD.txt", "w");
    for(double h=8;h<=16;h+=.5)
    {
        calc(JJ,MM,AAAA,h,48+48./60+45./3600,-(2+20./60+33./3600));
        fprintf(fp, " %lf, %lf %2.2f \n",fabs(1/tandeg(hauteur))*cosdeg(-azimut  + rot),
                fabs(1/tandeg(hauteur))*sindeg(-azimut  + rot),h);
    }
    fclose (fp);
 /**
The  position of the shadow at UTC 12:00 draws the analemma
*/ 
    for(double m=0;m<=12;m+=.25){
        calc(1,m,2020,12,48+48./60+45./3600,-(2+20./60+33./3600));
        fprintf(stdout, " %lf, %lf \n",fabs(1/tandeg(hauteur))*cosdeg(-azimut  + rot),
                fabs(1/tandeg(hauteur))*sindeg(-azimut  + rot));
    }
    
/**
 The  equation of time at Greewich:  one Newton step to find the time $T$ for which  azimut is South: $Z=180$), difference  $(T-12)$ is the "equation of time".
*/
    fp = fopen ("edt.txt", "w");
    for(double j=0;j<=365;j++){
        calc(j,1,2020,12,48+48./60+45./3600,0);
        double azimut12=azimut;
        calc(j,1,2020,12.1,48+48./60+45./3600,0);
        double dadt=(azimut-azimut12)/.1;
        double heure180 = 12 - (azimut12-180)/dadt;
        calc(j,1,2020,heure180,48+48./60+45./3600,0);
        heure180 = heure180 - (azimut-180)/dadt;
        fprintf(fp, " %lf, %lf \n",j,heure180);
    }
    fclose (fp);
    fprintf(stderr," no error, just to compute today the %2d/%2d/%4d at UTC=%2d:%2d\n",JJ,MM,AAAA,HH,mim);
    //exit (1); // always fails (forces rerun)
}
/**
# Results

## Exemple of Values for a given position and date
 
 Exemple of computation from Meeus 1986 implementation [http://xjubier.free.fr/site_pages/astronomy/ephemerides.html]()<br>
 26 April 2020 16:00 UTC,
 position 48°48'45'' N, 2°20'33'' E, so negative:
 -(2+20./60+33./3600)
 at the "Cité Verte" (close to Parisian meridian 2°20'14')<br>

 
~~~code
 #* * * * * * * * * * * *
 #-----------------------
 #The 26/ 4/2020 at UTC 16h 0mn0.0s is Nj=43581.67 Days since 00/01/1901
 #Temps Sidéral G = 95° 6mn48.5s  Temps Sidéral Heures = 6h20mn27.2s Temps Sidéral local = 6h29mn49.4s
 #-----------------------
 # HP [EQ]  : r=1.006472 l=36.886502 H=62°54mn25.7s  delta=13°48mn44.1s
 # HP [E->A]: haut h=28° 5mn35.9s Azimut Z= 258°31mn3.4s
 #-----------------------
 # Dec. val.: # Long=36.886502 Lat=0.00 haut=28.093306 Az=258.517622  Asc drt=2.303256  dec=13.812241 dist=1.006472
 
~~~
 
 
 
 <table>
 <tr>
 <td>  from  Meeus implementation </td>
 <td>----></td>
 <td>  this program</td>
 <td> Error</td>
 </tr>
 <tr>
 <td>  Long 36°52'27,03" = 36.8742  </td>
 <td>----></td>
 <td>  l = 36.886502 </td>
 <td>0.01° </td>
 </tr>
 <tr>
 <td>Lat 00°00'00,00"</td>
 <td>----></td>
 <td> 0°</td>
 <td> 0°</td>
 </tr>
 <tr>
 <td>hauteur 28°04m48.s = 28,08° </td>
 <td>----></td>
 <td> haut = 28°05mn35.9s =28.093306°</td>
 <td> 0.01°  (1111m) </td>
 </tr>
 <tr>
 <td> Azimut 258°31m12.s = 258,52°  </td>
 <td>----></td>
 <td> Azimut Z=  258°31mn3.4s = 258.517622° </td>
  <td> 0.002°</td>
 </tr>
 <tr>
 <td> Asc Drt 2h18m08,89s 2.303256  </td>
 <td>----></td>
 <td> Asc drt=2.303256° </td>
 <td> 0.00°</td>
 </tr>
 <tr>
 <td> déclinaison 13° 48' 27,65" </td>
 <td>----></td>
 <td> delta=13°48mn44.1s  </td>
 <td> 017''</td>
 </tr>
 <tr>
 <td>  dist 1,0064877 U </td>
 <td>----></td>
 <td> r=1.006472  </td>
 <td> 0.000157</td>
 </tr>
 </table>
## Example of anual sundial plot
 

Plot of path of Sun, the 21 of each month 1,2,3...12 with the UTC hours the
21/12, left bottom,
 20/03
 and 21/06 right top (2 and 10, 1 and 11, 3 and 9 (march and september), 4 and 8, 5 and 7 are almost superposed),
 with an analemma at noon.
 
The motto "transibunt et augebitur Scientia" is on the sundial of Cuvier's house,
 Jardin des Plantes  (see Camus Gotteland 1993, [http://michel.lalos.free.fr/cadrans_solaires/autres_depts/paris/cs_paris_05.html]())

~~~gnuplot
 reset
 set size square
 set title "transibunt et augebitur Scientia "
 set label 1 at 0,0 "  Gnomon"
 set label 12 at -2.25,-4.75 "December"
 set label 6 at 1.25,-1.35 "June"
 set label 3 at 1.25,-3 "March"
 set object circle at first 0,0 radius char 0.5 \
 fillstyle empty border lc rgb '#aa1100' lw 2
 p[-5:2][-5:2]'out' u 1:2 w l not,\
 'hours.txt' lc rgb '#aa1100' pt 7  ps 1  not,\
 '' with labels center offset 3.4 not,\
 'hoursTD.txt' w lp lc 3 t'today from 08:00 every 30min"
 unset label 12
 unset label 1
 unset label 3
 unset label 6
~~~

## Example of  sundial for one day plot

 The sundial at the day of the computation, for Cité Verte orientation.
~~~gnuplot
 reset
 set size square
 set title "transibunt et augebitur Scientia "
 set label 1 at 0.15,0.25 "Gnomon length,\nand gnomon position"
 set arrow 1  nohead from 0,0 to 0,1
 set arrow 2  nohead from -.5,0 to 0.5,0
 set object circle at 0,0 radius  .05
 p[:1.5][:1.5]'hoursTD.txt' using 1:2:($3-int($3)<0.5? sprintf("%d:00", $3): " ") with labels point  pt 7  offset 2.5,0 notitle,'' w l linec -1 not
 unset arrow 1
 unset arrow 2
 unset label 1

 ~~~
 
 
## Equation of time.


 This corresponds to the Analemma: azimut is not exactly 180° at UTC=12:00, there is a variation in time
 called "equation of time" visible in practice through the plot of Analemma.
 This is the difference in time between "true" an "apparent" sun.
 Here is the procedure we use to compute equation of time in Greewich:  find the hour at which
 $Z=180$ substract 12h and put in minutes.
 It depens on the obliquity of earth axis $\varepsilon$ and on the excentricty $e$ of the path of the earth.
The Meeus 1988,2014 formula is presented as a good approximation
 $$\Delta=\frac{5}{4} e^2 \sin \left(\frac{\pi  M}{90}\right)-4 e y_o \cos \left(\frac{\pi
 L}{90}\right) \sin \left(\frac{\pi  M}{180}\right)+2 e \sin \left(\frac{\pi
 M}{180}\right)+\frac{1}{2} y_o^2 \sin \left(\frac{\pi  L}{45}\right)-y_o
 \sin \left(\frac{\pi L}{90}\right)\text{ with }
 y_o= \tan ^2\left(\frac{\pi  \varepsilon}{360}\right)
 $$
 where $$L = 280.46646 + 36000.76983 T + 0.0003032 T^2$$
$$M=357.52911 + 35999.05029 T - 0.0001537 T^2,$$
$$\varepsilon = 23 + (26*60 + (21.448 -
 T (46.815 + T (0.00059 - T 0.001813))))/60/60,$$
 $$
 e =(0.016708634 - 0.000042037 T + 0.0000001267 T^2);$$
 where $T$ the number of centuries since 01/01/2000, so number of days divided bay
 36525. It is not plotted here, but it has been checked that the difference between Meeus formula
 and this programm is 15 seconds (of hour) at most.
 
 Meeus formula is simplified  in a simple approximation tacking into account obliquity and excentricty
 as a function of $j$ the day of the year, the first of january $j=1$, and $j=81$ at spring equinox.
 the approximate "equation of time" ([https://fr.wikipedia.org/wiki/Équation_du_temps]()) is just:
 $$\Delta=7.678 \sin(B(j) + 1.374) - 9.87 \sin( 2 B(j)) \text {  with  } B(j)=2 \pi (j - 81)/365$$
 It represents well the actual result:
~~~gnuplot
 set title " equation of time (french sign)"
 set xlabel "day in year"
 set ylabel "difference in min"
 B(x)=2*pi*(x - 81)/365.
 edt(x)=7.678*sin(B(x) + 1.374) - 9.87*sin( 2*B(x))
 p[0:365] edt(x) t'simple approximation','edt.txt'u 1:(($2-12)*60) t'computed'
~~~

## Conclusion
 
 Position of sun and moon (not fully finshed yet) usefull for tides (see maree_bretagne.c), navigation.
 
 Position of sun : Bagnold sun compass, invented by Bagnold.
 
 
 
# Biblio
 
 * Serge Bouiges "Calcul Astronomique pour Amateur, adapté à l'emploi d'un calculateur ou d'un micro ordinateur" Masson 1981
 
 * Jean Meeus "Astronomical algorithms" 2nd Edition, Willmann-Bell inc., 1998
 
 * Jean Meeus "Calculs Astronomiques" Nouvelle Edition 2014 (thanks Jérôme H.)
 
 * Guy Sérane "Astronomie & ordinateur" Dunod 1987
 
 * [http://xjubier.free.fr/site_pages/astronomy/ephemerides.html]() calcul en ligne
 
 * [http://f1rzv.free.fr/calculs/AstroCalc.php]() calcul en ligne
 
 * Bureau des Longitudes,  "Ephémérides Nautiques 2020"
 
 * Bureau des Longitudes,  "Ephémérides Astronomiques"
 [https://www.imcce.fr/content/medias/publications/publications-institutionnelles/CDT_2020_ebook.pdf]()
 
 *  Michel Capderou "Satellites - De Kepler au GPS" Springer 2011, 844 p
 
 * PYL copie perso de "table des marées des grands ports du monde" du SHOM  (1984) [tabledesmareesdesgrandsportsdumonde.pdf](https://mycore.core-cloud.net/index.php/s/5Y5mKjr1pWmj7aD)
 
 * [https://fr.wikipedia.org/wiki/Équation_du_temps]()
 
 * [http://freveille.free.fr/Equation_du_temps.html]()
 
 * Maillart & Millet "Cosmographie, classe de Mathématiques", Hachette 1955
 
 * Y. Opizzo "Cadrans solaires de précision" Masson 1990

 * Georges Camus, Andrée Gotteland "Cadrans solaires de Paris" 1993
 [https://www.cnrseditions.fr/catalogue/histoire/cadrans-solaires-paris-georges-camus/]()
 
 * [http://michel.lalos.free.fr/cadrans_solaires/autres_depts/paris/cs_paris_05.html]() images et positions de cadrans parisiens

 * Kuno Gross "The Bagnold sun-compass, history and utilization" Books on demand 2011
 
 * R. A. Bagnold and W. J. Harding King "Journeys in the Libyan Desert 1929 and 1930",
 The Geographical Journal, Vol. 78, No. 6 (Dec., 1931), pp. 524+526-535
 [https://www.jstor.org/stable/1784967]()
 
 * Loic Caradec, "une calculatrice à bord, pourquoi et comment?"
 Editions Maritimes et d'Outre Mer 1977
 
 * Amiral Sacaze, "Navigation Astronomique simplifiée",  Editions Maritimes et d'Outre Mer 1980
 
 * P Posth navigation astronomique Navigation astronomique et calculatrices programmables [https://books.google.fr/books?id=J-l2QgGTgSQC&printsec=frontcover&dq=posth+navigation+astronomique&hl=fr&sa=X&ved=0ahUKEwjd2Om7vojpAhVLxYUKHZzWCmMQ6AEIKDAA#v=onepage&q=posth%20navigation%20astronomique&f=false]()
 
 * [http://www.breves-de-maths.fr/trois-nombres-et-un-siecle-pour-decrypter-anticythere/]() machine d'Anticythère
 
 * [http://dlib.nyu.edu/awdl/isaw/isaw-papers/4/]() machine d'Anticythère
 
 * [https://www.hpmuseum.org/forum/archive/index.php?thread-4262.html]() 59 of the 60 terms of Meuus coded in HP-BASIC
 
 * [http://www.geoastro.de/elevaz/basics/meeus.htm]() Meuus Formulas for sun
 
 * [http://www.stargazing.net/kepler/sun.html]() sun
 
 * [http://www.stargazing.net/kepler/moon.html]() 6 terms for Moon
 

 
 <br><br>
conf Avril 2020
*/
