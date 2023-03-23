#include <stdio.h> 
#include "seyfert.h" 

/* program to get column densities for various molecular tranistions (wc95,wc96b) */

/* THE REAL WAY */
main()
{char mol[MAXSTRING],ans[MAXSTRING],line[MAXSTRING];
  float frest, T, T1, z, N, tau, kT, f, Q, h, blah, A, E, v, c, E1, x, mu, J,  g, gu,E_1612,E_1665,E_1667,E_1720,f_1612,f_1665,f_1667,f_1720;
 int  i,j; /* J,  g, gu where originally integers (N-int.c) but want to put in OH */


printf("\nEnter molecule (e.g. co): ");
scanf("%s", mol);
 if ( strcmp(mol, "co") == 0  ||  strcmp(mol, "cs") == 0 ||  strcmp(mol, "hco+") == 0 ||  strcmp(mol, "hcn") == 0 ||  strcmp(mol, "oh") == 0){
  if (strcmp(mol, "oh") != 0){

printf("\nEnter lower level (e.g. 1 for J=1-2): ");
scanf("%f", &J);

 int top;
 if (strcmp(mol, "co") == 0) top = 9;
 else top = 7;   // top (lower) levels I've got to adding

 if (J > top){ 
    printf("\nThis isn't in yet!\n");
 exit(1);
 }
 else{

printf("\nEnter excitation temperature: ");
scanf("%f", &T1);

printf("Constant - enter C\nRedshift dependent - enter R:\n");
scanf("%s", ans);

 if ( strcmp(ans, "C") == 0 ){
   T = T1;

 }
 T = T1;

if ( strcmp(ans, "R") == 0 ){
printf("\nEnter the redshift: ");
scanf("%f", &z);
    T = ((1+z)*2.7)+(T1-3);   /* does the trick for z up to 1000 anyway */
  }

/* printf("\nT = %1.0f K", T); */
printf("\nEnter the line width (km/s):");
scanf("%f", &v);
printf("\n...and the optical depth: ");
scanf("%f", &tau);
 kT = T*1.381e-23; /* multipling by Boltzmann constant */
 
                      /*    printf("\nkT = %1.3e J\n", kT); */
 h = 6.626e-34; /* Planck constant */

 /* to generalise */ E = (0.5*J*(J+1))*E1;

if ( strcmp(mol, "co") == 0 ){

 E1 = 3.8450*1.9878e-23; /* energy of J=1 - convering JPL cm^-1 to Joules */

 E = (0.5*J*(J+1))*E1;
  if (J==0){
 
    frest = 115271.2018e6; // from http://spec.jpl.nasa.gov/ftp/pub/catalog/c028001.cat
    /*   A = 7.4e-8;*/ /* from Rohlfs p 390 */
    A = 7.7e-8; /* http://cdsarc.u-strasbg.fr/viz-bin/Cat?J/A%2bAS/117/557  - cms96 */
  }               // SEARCH VIA ROTATION QUANTUM NUMBERS AND SETTING Vibrational quantum number TO 0
 if (J==1){
 
    frest = 230.538e9; 
   
  /*  A = 7.1e-7; *//* from Rohlfs p 390 */
   A = 7.4e-7; /* cms96 */
}
 if (J==2){
 
   frest = 345.7959899e9;
   A = 2.7e-06;
 }
 if (J==3){
 
   frest = 461.0407682e9;
   A = 6.5e-06;
 }
  if (J==4){
 
    frest = 576267.9305e6;
   A = 1.3e-05;
 } 
  if (J==5){
 
    frest = 691473.0763e6;
   A = 2.3e-05;
 }
if (J==6){
 
  frest = 806651.8060e6;
   A = 3.7e-05;
 }
if (J==7){
 
  frest = 921799.7000e6; 
   A = 5.5e-05;
 }

if (J==8){
 
  frest = 1036912.3930e6;
   A = 7.8e-05;
 }

if (J==9){
 
  frest = 1151985.4520e6;
   A = 1.1e-04;
 }


}

if ( strcmp(mol, "cs") == 0 ){

  E1 = 1.6342*1.9878e-23;
  E = (0.5*J*(J+1))*E1;
 if (J==0){
 
  frest = 48.9909780e9;
  A = 1.74e-06;
 }
 if (J==1){
 
    frest = 97.9809500e9;
  A = 1.68e-05;
 }
 if (J==2){
 
   frest = 146.9690330e9;
  A = 6.05e-05;
 }
  if (J==3){
 
   frest = 195.9542260e9;
  A = 1.49e-04;
 }
  if (J==4){
 
   frest = 244.93564350e9;
  A = 2.97e-04;
 }
  if (J==5){
 
   frest = 293.9122440e9;
  A = 5.21e-04;
 }
  if (J==6){
 
   frest = 342.8830000e9;
  A = 8.37e-04;
 }
  if (J==7){
  
   frest = 391.8470300e9;
   A = 1.26e-03;          /* all of these also from CDS cklh95 */
 }

}


if ( strcmp(mol, "hco+") == 0 ){
  
  E1 = 2.975008479*1.9878e-23; // http://home.strw.leidenuniv.nl/~moldata/molformat.html
  E = (0.5*J*(J+1))*E1;      // also from http://spec.jpl.nasa.gov/ftp/pub/catalog/c029002.cat
                            //  178375.0100  0.0500 -1.3699 2    2.9750  5 -290021303 2 0 0       1 0 0 
  // For DCO+ would this be 144077.2144  0.0072 -2.7245 2    2.4030  3 -300031304 2 0 0 1     1 0 0 1 
  // i.e. 2.403 


   /* have to work out A as not given in literature */
   //  mu = 4.48; /* Debyes - see MGB's sheet and Rohlfs p379 */  
   mu = 3.888; //http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d029002.cat
   /// can check with RADEX online  http://www.sron.rug.nl/~vdtak/radex/radex.php


 if (J==0){
 
    frest = 89.188518e9; 
        
     A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
   
    //this gives A = 5.530e-05, cf. 3.0e-5 from p390 of Rohlfs // printf("\nA = %1.3e \n", A); 

 }
 if (J==1){
 
   frest = 178.3750650e9; 
    A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
       // printf("\nA = %1.3e \n", A);

 }
  if (J==2){
  frest = 267.5576190e9; 
  
   A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
   // A =  1.0e-3;
  printf("\nA = %1.3e \n", A);
 
/*  A = 1.0e-3; *//* from Rohlfs p 390 */
  }
  if (J==3){
  
   frest = 356.7342880e9; 
  A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
  }
  if (J==4){
  
   frest = 445.9029960e9; 
  A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
  }
  if (J==5){
  
   frest = 535.0617755e9; 
  A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
  }
  if (J==6){
  
   frest = 624.2086733e9; 
  A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
  }
  if (J==7){
  
   frest = 713.3420900e9; 
  A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
  }


}
if ( strcmp(mol, "hcn") == 0 ){

  E1 = 2.9564*1.9878e-23;
   E = (0.5*J*(J+1))*E1;

   //mu = 2.98; /* Debyes - see MGB's sheet and Rohlfs p379 */  
   mu = 2.942; //http://spec.jpl.nasa.gov/ftp/pub/catalog/doc/d027003.cat

 if (J==0){
  
    frest = 88.6318470e9;
    A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
 }
 if (J==1){
  
   frest = 177.2612230e9;
   A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
 }
 if (J==2){
  
  frest = 265.8861800e9;
  A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
 }
 if (J==3){
  
   frest = 354.5054759e9;
   A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
 }
 if (J==4){
  
   frest = 443.1161554e9;
   A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
 }
 if (J==5){
  
   frest = 531.7163875e9;
 A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
 }
 if (J==6){
  
   frest = 620.30409525e9;
 A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
 }
 if (J==7){
  
   frest = 708.8772081e9;
   A = (1.165e-11)*(pow((frest*1e-9),3))*(pow(mu,2))*(J+1)/(2*J+3);
}
  
 
  /* printf("\nA = %1.3e \n", A); */
 
}
/* printf("\nfrest = %1.3f \n", frest);*/

 blah = exp(E/kT)/(1-(exp((-h*frest)/(kT))));
 
  gu = 2*(J+1)+1;  /* upper level so e.g. should be 5 for 1-2 */
  g = (2*J)+1;

  //////// PARTITION FUNCTION /////////////////////////////// 

  x = -E1/kT;
  Q = 1; 
    j = 1;
    
     while (((2*j)+1)*exp(0.5*j*(j+1)*x) > 0.001){  /* summing to this tolerance */
       //   printf("%d %e %e \n", j, x, ((2*j)+1)*exp(0.5*j*(j+1)*x));
   Q = Q + (((2*j)+1)*exp(0.5*j*(j+1)*x));
  j = j+1; /*this is J, but use j to save further use of J */

  //  printf("%f  \n", ((2*j)+1)*exp(0.5*j*(j+1)*x));
}
  f = Q*blah;


 c = 3e5; /* km/s */

 N = (8*3.1416*(pow(frest,3))*f*(1.06*tau*v))/((pow(c,3))*A*gu*1e10); /* last number converts from km to cm */


printf("\nFrom %s %1.0f-%1.0f at %1.0f K, tau = %1.3f and a FWHM of %1.0f km/s\nthe total column density of %s is N = %1.1e cm^-2\n", mol, J, J+1, T, tau, v, mol, N);
 printf("(Used Einstein coefficient of A%1.0f,%1.0f = %1.3e)",J, J+1, A);

printf("\nArbitrarily high T_x now also applied to papers/22/limits.c\n");
 }
  }
  
  
  if ( strcmp(mol, "oh") == 0 ){ 

    printf("All 18-cm just now, but will add 6-cm soon\n");
    
    // E_1612 = 0.0018*1.9878e-23; //1612 MHz
    // E_1665 = 0.0000; //1665
    // E_1667 = 0.0018*1.9878e-23;  /* 1667 MHz from JPL */ //convering JPL cm^-1 to Joules
    // E_1720 = 0.0000; //1720 MHz all from JPL

    //  f_1612 = 1612.2310e6;
    //  f_1665 = 1665.4018e6;
    // f_1667 = 1667.3590e6;
    // f_1720 = 1720.5300e6;

printf("\nEnter excitation temperature: ");
scanf("%f", &T1);

printf("Constant - enter C\nRedshift dependent - enter R:\n");
scanf("%s", ans);

 if ( strcmp(ans, "C") == 0 ){
   T = T1;

 }
 T = T1;

if ( strcmp(ans, "R") == 0 ){
printf("\nEnter the redshift: ");
scanf("%f", &z);
    T = ((1+z)*2.7)+(T1-3);   /* does the trick for z up to 1000 anyway */
  }

printf("\nEnter the line width (km/s):");
scanf("%f", &v);
printf("\n...and the optical depth: ");
scanf("%f", &tau);
 kT = T*1.381e-23; /* multipling by Boltzmann constant */
  h = 6.626e-34; /* Planck constant */


  printf("\nEnter the line (1612, 1665, 1667 or 1720 MHz) : ");
      scanf("%s", line);
      if ( strcmp(line, "1612") == 0 ){
	E = 0.0018*1.9878e-23; //1612 MHz from JPL
	A = 1.29e-11; // from cm67
	frest = 1612.2310e6;
	// in emission is 1 -> 2 (lid67) so lower level is 2
	J = 2; 
      }
      if ( strcmp(line, "1665") == 0 ){
	E = 0.0000; //1665 from JPL
	A = 7.11e-11; // from cm67
        frest = 1665.4018e6;
	J = 1; // in emission is 1 -> 1 (lid67) so lower level is 1
      }	
      if ( strcmp(line, "1667") == 0 ){
        E = 0.0018*1.9878e-23; //1667 MHz from JPL
	A = 7.71e-11; // from cm67
       frest = 1667.3590e6;
	J = 2; // in emission is 2 -> 2 (lid67) so lower level is 2
      }	
      if ( strcmp(line, "1720") == 0 ){
	E = 0.0000; //1720 from JPL
	A = 0.94e-11; // from cm67
       frest = 1720.5300e6;
	J = 1; // in emission is 2 -> 1 (lid67) so lower level is 1
      }
      c = 3e5; /* km/s */
      
      // don't worry about partitition function for OH, lower level given by (pkn+04)
      gu = 1; //set
      
       blah = exp(E/kT)/(1-(exp((-h*frest)/(kT))));
       //      printf("blah = %1.3g\n",blah); 
      
       // FOR THE 4 TRANSITIONS LOWER LEVELS GIVE g = 3+3+5+5=16; and 
       g = (2*J)+1;
      f = 16/g; // multiplying by this weighting gives total column density
      //  printf("f = %1.3f",f);

      N = (8*3.1416*(pow(frest,3))*f*blah*(1.06*tau*v))/((pow(c,3))*A*gu*1e10); //last number converts from km to cm

      printf("\nAt T = %1.0f K, tau = %1.4f and a FWHM of %1.0f km/s\nthe total column density from %s MHz OH is N = %1.2e cm^-2\n",  T, tau, v, line, N);

      // printf("\nCF. T = %1.0f K, tau = %1.4f and a FWHM of %1.0f km/s\nthe total column density from 1667 MHz OH is N = %1.2e cm^-2\n",  T, tau, v, 2.38*1e14*T*tau*v);
      // printf("\nand from 1665 MHz OH is N = %1.2e cm^-2\n",4.30*1e14*T*tau*v);
// printf("\n[See Plume et al., 2004 (astro-ph/0401482)]\n");


  }
 }
 
 else{
    printf("\nThis isn't in yet!\n");
 exit(1);
 }
 
}
 
