/** 
   \file utils.c
   Utility definitions & functions (miscellaneous).
	
   Contains some utility functions and definitions.
	
*/
#include "MultiCoreUtils.h"

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>


int debug;
int verbose;
char parseFileDelims[] = " \t\n";



/***********************************************************************************
 *	logSumOfExponents
 * 	- takes an array of (log) values and returns the log of the sums of the exponents of the log-values
 *	- does this by exponentiating only ratios, to prevent under/over flow.
 *	- indsToSum is a boolean array indicating for each index whether to sum it or not. (if NULL, then sum all)
 ***********************************************************************************/
double logSumOfExponents(double* logArray, int arrayLength, unsigned short* indsToSum) {

  double MAX_RATIO = 10.0;
  int ind, maxInd;
	
  double maxLog;			// records maximum value of array
  double sumExpRatios;	// sum of the ratios  e^{A[ind]-maxLog}
  double diff;			// difference maxLog - A[ind] 
	
  if(arrayLength <= 0) {
    return logArray[0];
  }
	
  // perform first loop to figure out maxLog (and potentially avoid taking too many exponents
  maxLog = 0.0;
  maxInd = -1;
  for(ind=0; ind<arrayLength; ind++) {
    if(indsToSum != NULL && !indsToSum[ind]) {
      continue;
    }
    if(maxInd< 0 || logArray[ind] > maxLog) {
      maxLog = logArray[ind];
      maxInd = ind;
    }
  }
	
  if(maxInd < 0) {
    //		printf("\n\n----- no valid values to sum in log array -----\n\n");
    return -1e34;
  }
	
  sumExpRatios = 1;	// contribution of maxInd
  for(ind=0; ind<arrayLength; ind++) {
    if((ind == maxInd) || (indsToSum != NULL && !indsToSum[ind])) {
      continue;
    }
		
    diff = maxLog-logArray[ind];
    if(diff < MAX_RATIO) {
      sumExpRatios += exp(-diff);
    } // otherwise logArray[ind] is negligible
  }// end of for(ind)
	
  // return sum of maxlog and log of sumExpRatios
  // COULD POTENTIALLY AVOID LOG, IF SUFFICIENTLY CLOSE TO 1 !!!
  return maxLog + log(sumExpRatios);
}
/** end of logSumOfExponents **/


/**  Safely read in next string from file
 *    @param stream File descriptor to read from
 *	 @param MaxLength Length of buffer to read into
 *    @param destination Char buffer of size MaxLength, 
 *	        which will be filled with string read from file
 *		returns number of characters read (even if not into buffer
 */
int readStringFromFile(FILE *stream, int bufferLength, char *destination) {
	int length;
	char c;
  
	//Clear destination buffer  
	for(length=0; length<bufferLength; length++) {
		destination[length] = '\0';
	}

	c = ' ';
	length = 0;
	//Read by character either 'MaxLength characters' or 
	//'until white space', whichever comes first
	while(length<bufferLength-1) {
		c = fgetc(stream);
		if(feof(stream)) {
			return length;
		}
		
		if (!isspace(c)) {
			destination[length++] = c;
		} else if(length > 0) {
			return length;
		}
	}
	
	// finish reading the srting, and update length, but not buffer 
	//We have to read until a whitespace, just don't add to the buffer
	while(!isspace(fgetc(stream)) && !feof(stream)) {
		length++;
	}

	return length;
}
/** end of readStringFromFile **/


/***********************************************************************************
 *	mergeArrays
 * 	- merges two sorted portions of a single array	(first part should typically be shorter than second)
 *	- borderPoint indicates the start point of the second array
 *	- tmpSpace array is provided as a temporary workspace (of size borderPoint)
 ***********************************************************************************/
void mergeArrays(double array[], int numEntries, int borderPoint, double tmpArray[]) {
	
  int i,j;
	
  // copy first half onto temp space
  for(i=0; i<borderPoint; i++) {
    tmpArray[i] = array[i];
  }
  // note that the size of available space in array[] is always
  // no shorter than the size of the unvisited portion of tmpSpace[]
  for(i=0, j=borderPoint; i<borderPoint; /*either i or j gets advanced within loop*/) {
		
		if(j == numEntries || tmpArray[i] < array[j]) //previous formatting was causing warnings
		{  
		  array[i+j-borderPoint] = tmpArray[i];
		  i++;
		} else {
		  array[i+j-borderPoint] = array[j];
		  j++;
		}
    //array[i+j-borderPoint] = (j == numEntries || tmpArray[i] < array[j]) ?  (tmpArray[i++]) : (array[j++]);
  }
  return;
}



/***********************************************************************************
 *	mergeSrot
 * 	- implementation of merge sort
 *	- receives an array of doubles and its length and sorts the array in non-descending order
 *	- a recursive procedure
 *	- tmpSpace array is provided as a temporary workspace (of size numEntries)
 *	- returns the sum of all values in array
 ***********************************************************************************/
double mergeSort(double array[], int numEntries, double tmpArray[]) {
  int i, j;
  double sum;

  if(numEntries < 6) {
    sum = array[0];
    // perform bubble sort
    for(i=1; i<numEntries;i++) {
      sum += array[i];
      for(j=0;j<i; j++) {
        if(array[j]>array[i]) {
          swap2(array[j],array[i]);
        }
      }
    }
    return sum;
  }
	
  i = numEntries/2;
  sum  = mergeSort(array  ,      i       , tmpArray);
  sum += mergeSort(array+i, numEntries-i , tmpArray);
	
  mergeArrays(array, numEntries, i, tmpArray);
	
  return sum;
}

void turnOnBooleanArray(unsigned short* booleanArray, int arrayLength) {
  for( arrayLength--; arrayLength>=0; arrayLength--) {
    booleanArray[arrayLength] = 1;
  }
  return;
}



// time functions

static time_t time_start;

void starttime (void)
{
	time_start=time(NULL);
}


static TIMERS runTimes;

void setStartTimeMethod(enum METHOD_NAME method)
{
	time_t tic = time(NULL);
	switch(method)
	{
		case T_UpdateGB_InternalNode: runTimes.UpdateGB_InternalNode.start = tic; break;
		case T_UpdateGB_MigrationNode: runTimes.UpdateGB_MigrationNode.start = tic;break;
		case T_UpdateGB_MigSPR: runTimes.UpdateGB_MigSPR.start = tic;break;
		case T_UpdateTau: runTimes.UpdateTau.start = tic;break;
		case T_UpdateMigRates: runTimes.UpdateMigRates.start = tic;break;
		case T_mixing: runTimes.mixing.start = tic;break;
		case T_UpdateTheta: runTimes.UpdateTheta.start = tic;break;
		case T_UpdateSampleAge: runTimes.UpdateSampleAge.start = tic;break;
		case T_UpdateLocusRate: runTimes.UpdateLocusRate.start = tic;break;
		case T_MCMCIterations: runTimes.MCMCIterations.start = tic;break;
	}
}

void setEndTimeMethod(enum METHOD_NAME method)
{
	time_t toc = time(NULL);
	switch(method)
	{
		case T_UpdateGB_InternalNode: runTimes.UpdateGB_InternalNode.accumulated += (toc - runTimes.UpdateGB_InternalNode.start); break;
		case T_UpdateGB_MigrationNode: runTimes.UpdateGB_MigrationNode.accumulated += (toc - runTimes.UpdateGB_MigrationNode.start); break;
		case T_UpdateGB_MigSPR: runTimes.UpdateGB_MigSPR.accumulated += (toc - runTimes.UpdateGB_MigSPR.start); break;
		case T_UpdateTau: runTimes.UpdateTau.accumulated += (toc - runTimes.UpdateTau.start); break;
		case T_UpdateMigRates: runTimes.UpdateMigRates.accumulated += (toc - runTimes.UpdateMigRates.start); break;
		case T_mixing: runTimes.mixing.accumulated += (toc - runTimes.mixing.start); break;
		case T_UpdateTheta: runTimes.UpdateTheta.accumulated += (toc - runTimes.UpdateTheta.start); break;
		case T_UpdateSampleAge: runTimes.UpdateSampleAge.accumulated += (toc - runTimes.UpdateSampleAge.start); break;
		case T_UpdateLocusRate: runTimes.UpdateLocusRate.accumulated += (toc - runTimes.UpdateLocusRate.start); break;
		case T_MCMCIterations:  runTimes.MCMCIterations.accumulated += (toc - runTimes.MCMCIterations.start); break;

	}
}
void printMethodTimes()
{
#ifdef RECORD_METHOD_TIMES

	char timeString[STRING_LENGTH];
	printf("===== METHOD RUN TIME ======\n");
	printf("UpdateGB_InternalNode (sec):            %s\n", printtime_i(runTimes.UpdateGB_InternalNode.accumulated , timeString));
	printf("UpdateGB_MigrationNode (sec):           %s\n", printtime_i(runTimes.UpdateGB_MigrationNode.accumulated, timeString));
	printf("UpdateGB_MigSPR (sec):                  %s\n", printtime_i(runTimes.UpdateGB_MigSPR.accumulated, timeString));
	printf("UpdateTau (sec):                        %s\n", printtime_i(runTimes.UpdateTau.accumulated, timeString));
	printf("UpdateMigRates (sec):                   %s\n", printtime_i(runTimes.UpdateMigRates.accumulated, timeString));
	printf("mixing (sec):                           %s\n", printtime_i(runTimes.mixing.accumulated, timeString));
	printf("UpdateTheta (sec):                      %s\n", printtime_i(runTimes.UpdateTheta.accumulated, timeString));
	printf("UpdateSampleAge (sec):                  %s\n", printtime_i(runTimes.UpdateSampleAge.accumulated, timeString));
	printf("UpdateLocusRate NO MT(sec):             %s\n", printtime_i(runTimes.UpdateLocusRate.accumulated, timeString));
	time_t total = 0;
	total += runTimes.UpdateGB_InternalNode.accumulated;
	total += runTimes.UpdateGB_MigrationNode.accumulated;
	total += runTimes.UpdateGB_MigSPR.accumulated;
	total += runTimes.UpdateTau.accumulated;
	total += runTimes.UpdateMigRates.accumulated;
	total += runTimes.mixing.accumulated;
	total += runTimes.UpdateTheta.accumulated;
	total += runTimes.UpdateSampleAge.accumulated;
	printf("===== Total in MT Methods:              %s\n" , printtime_i(total, timeString));
	total += runTimes.UpdateLocusRate.accumulated;
	printf("===== Total in All MCMCM Methods        %s\n" , printtime_i(total, timeString));
	printf("===== Total in MCMCM Iterations:        %s\n" , printtime_i(runTimes.MCMCIterations.accumulated, timeString));
#endif
}
char *printtime_i(int t , char timestr[])
{
	  int h, m, s;

	  h=t/3600; m=(t%3600)/60; s=t-(t/60)*60;
	  if(h)  sprintf(timestr,"%d:%02d:%02d", h,m,s);
	  else   sprintf(timestr,"00:%02d:%02d", m,s);
	  return(timestr);
}
char* printtime (char timestr[])
{
  /* print time elapsed since last call to starttime()
   */
  time_t t;
  int h, m, s;

  t=time(NULL)-time_start;
  h=t/3600; m=(t%3600)/60; s=t-(t/60)*60;
  if(h)  sprintf(timestr,"%d:%02d:%02d", h,m,s);
  else   sprintf(timestr,"%2d:%02d", m,s);
  return(timestr);
}


/*-----------------------------------------------------------------------------
 * This returns a variable in the range (a,b) by reflecting
 * x back into the range
 */
double reflect(double x, double a, double b )
{

  // safety margins for upper and lower bounds
  static double slack = 0.000000001;
  double xnew, double_interval;

  a += slack;
  b -= slack;

  // if interval is empty (due to slackness),
  // return middle of interval
  if( b <= a )
  {
    //fprintf(stderr, "very small interval in reflect(%g,%g,%g)"
    //        " [slackness = %g].\n",x,a-slack,b+slack,slack);
    return (a+b)/2.;
  }

  if( x < b && x > a )
    return x;

  // reflect upwards, if necessary
  xnew = x;
  if( xnew <= a )
    xnew = 2. * a - xnew;

  // "fold twice" as many time as needed
  double_interval = 2. * ( b - a );
  xnew = xnew - double_interval*floor( (xnew-a) / double_interval );

  // reflect downwards one last time, if necessary
  if( xnew >= b )
    xnew = 2. * b - xnew;

  // value should be within interval at this stage, but numerical
  // precision might put it slightly outside
  while( xnew <= a || xnew >= b )
  {
    if( xnew >= b )
    {
      if( debug )
      {
        fprintf(stderr, "reflect percision in reflect(%g,%g,%g): "
            "obtaining reflection at %g (%g greater than upper "
            "bound)\n", x,
                        a - slack,
                        b + slack,
                        xnew, xnew-b);
      }
      xnew = 2. * b - xnew;
    }
    else
    {
      if( debug )
      {
        fprintf(stderr, "reflect percision in reflect(%g,%g,%g): obtaining "
                        "reflection at %g (%g smaller than lower bound)\n",
                        x, a - slack,
                        b + slack, xnew, a-xnew);
      }
      xnew = 2*a - xnew;
    }
  }
  return xnew;
}

//================= Random Generator related functions ========================
RandGeneratorContext RndCtx;

//-----------------------------------------------------------------------------
#define MALLOC_AND_ASSIGN(ptr, t, sz, val) \
  ptr = (t*) malloc(sz); \
  if (NULL == ptr ) \
    printf("Error on Random context allocation"); \
  for(i=0; i < RndCtx.nOfSlots; ++i) \
    ptr[i] = val;

void initRandomGenerator( int nNumLoci, unsigned int unSeed )
{
  if( 4 != sizeof(int) )
    puts("oh-oh, we are in trouble. int is not 32-bit?");

  //The last extra slot is for general purpose computations
  RndCtx.nOfSlots = nNumLoci + 1;

  int bs = sizeof(int) * RndCtx.nOfSlots;
  int i = 0;
  unsigned int v = 170 * (unSeed % 178) + 137;
  MALLOC_AND_ASSIGN( RndCtx.rndu_z, unsigned int, bs, v ) //137
  v = unSeed*127773;
  MALLOC_AND_ASSIGN( RndCtx.rndu_w, unsigned int, bs, v ) //123456757
  MALLOC_AND_ASSIGN( RndCtx.rndu_x, unsigned int, bs, 11 )
  MALLOC_AND_ASSIGN( RndCtx.rndu_y, unsigned int, bs, 23 )
  bs = sizeof(double) * RndCtx.nOfSlots;
  MALLOC_AND_ASSIGN( RndCtx.m2s2_kernel, double, bs, 8. )
  double dv = sqrt(RndCtx.m2s2_kernel[0]/(RndCtx.m2s2_kernel[0] + 1.));
  MALLOC_AND_ASSIGN( RndCtx.m2N_kernel, double, bs, dv )
  dv = sqrt(1./(RndCtx.m2s2_kernel[0] + 1.));
  MALLOC_AND_ASSIGN( RndCtx.s2N_kernel, double, bs, dv )

  MALLOC_AND_ASSIGN( RndCtx.rndgamma2_b, double, bs,  0.0 )
  MALLOC_AND_ASSIGN( RndCtx.rndgamma2_h, double, bs,  0.0 )
  MALLOC_AND_ASSIGN( RndCtx.rndgamma2_ss, double, bs, 0.0 ) //0.0

  MALLOC_AND_ASSIGN( RndCtx.rndgamma1_a, double, bs,  0.0 )
  MALLOC_AND_ASSIGN( RndCtx.rndgamma1_p, double, bs,  0.0 )
  MALLOC_AND_ASSIGN( RndCtx.rndgamma1_uf, double, bs, 0.0 )
  MALLOC_AND_ASSIGN( RndCtx.rndgamma1_ss, double, bs, 10.0 ) //10.0
  MALLOC_AND_ASSIGN( RndCtx.rndgamma1_d, double, bs,  0.0 )
/*
  rndu_z=170*(seed%178)+137;
  rndu_w = seed*127773;

  m2N_kernel = sqrt(m2s2_kernel/(m2s2_kernel+1));  
  s2N_kernel = sqrt(1/(m2s2_kernel+1));
*/
}

/*-----------------------------------------------------------------------------
  standard normal variate, using the Box-Muller method (1958), improved by
  Marsaglia and Bray (1964).  The method generates a pair of random
  variates, and only one used.
  See N. L. Johnson et al. (1994), Continuous univariate distributions,
  vol 1. p.153.
*/
double rndnormal( int nLocusIdx )
{
  double u, v, s;
  while( 1 )
  {
    u = 2*rndu( nLocusIdx ) - 1;
    v = 2*rndu( nLocusIdx ) - 1;
    s = u * u + v * v;
    if( s > 0 && s < 1 )
      break;
  }
  s = sqrt( -2. * log( s ) / s );
  return u * s;
}

/*-----------------------------------------------------------------------------
   This returns a variate from the mixture of two normals
   N(-m, s2) and N(m, s2), with mean 0 and variance m^2 + s2 = 1
   and m^2/s^2 = 8.

   Let this standard variate be z.  Then mean + z * sigma will be a variate
   with mean mean and SD sigma.  This is useful for generating MCMC proposals
*/
double rnd2normal8( int nLocusIdx )
{
  double z =   RndCtx.m2N_kernel[nLocusIdx]
             + rndnormal( nLocusIdx ) * RndCtx.s2N_kernel[nLocusIdx];
  z = rndu( nLocusIdx ) < 0.5 ? z : -z;
  return z;
}

/*-----------------------------------------------------------------------------
   U(0,1): AS 183: Appl. Stat. 31:188-190
   Wichmann BA & Hill ID.  1982.  An efficient and portable
   pseudo-random number generator.  Appl. Stat. 31:188-190

   x, y, z are any numbers in the range 1-30000.  Integer operation up
   to 30323 required.
*/
double rndu( int nLocusIdx )
{
  double r;

  RndCtx.rndu_x[nLocusIdx] =    171 * ( RndCtx.rndu_x[nLocusIdx] % 177 )
                             -  2 * ( RndCtx.rndu_x[nLocusIdx] / 177 );
  RndCtx.rndu_y[nLocusIdx] =    172 * ( RndCtx.rndu_y[nLocusIdx] % 176 )
                             -  35 * ( RndCtx.rndu_y[nLocusIdx] / 176 );
  RndCtx.rndu_z[nLocusIdx] =    170 * ( RndCtx.rndu_z[nLocusIdx] % 178 )
                             -  63 * ( RndCtx.rndu_z[nLocusIdx] / 178 );
  r =   RndCtx.rndu_x[nLocusIdx] / 30269.0
      + RndCtx.rndu_y[nLocusIdx] / 30307.0
      + RndCtx.rndu_z[nLocusIdx] / 30323.0;
  r = ( r - (int)r );
  return r;
}

/*-----------------------------------------------------------------------------
   random standard gamma (Mean=Var=s,  with shape parameter=s, scale para=1)
   r^(s-1)*exp(-r)
   J. Dagpunar (1988) Principles of random variate generation,
   Clarendon Press, Oxford
   calling rndgamma1() if s<1 or
   rndgamma2() if s>1 or
   exponential if s=1
   @@TODO: (still actual?)
   This is unsafe, and is found to return 0 when s is very small.
*/
double rndgamma( int nLocusIdx, double s )
{
  double r=0;
  if ( s <= 0 )
    puts ("jgl gamma..");
  else if( s < 1)
    r = rndgamma1( nLocusIdx, s );
  else if( s > 1 )
    r = rndgamma2( nLocusIdx, s );
  else
    r = -log( rndu( nLocusIdx ) );
  return r;
}

//-----------------------------------------------------------------------------
double rndgamma1( int nLocusIdx, double s )
{
  /* random standard gamma for s<1
     switching method
  */
  double r, x=0,small=1e-37,w;
  if( s != RndCtx.rndgamma1_ss[nLocusIdx] )
  {
    RndCtx.rndgamma1_a[nLocusIdx] = 1 - s;
    RndCtx.rndgamma1_p[nLocusIdx] =
                   RndCtx.rndgamma1_a[nLocusIdx]/(RndCtx.rndgamma1_a[nLocusIdx]
                   + s * exp(-RndCtx.rndgamma1_a[nLocusIdx]));
    RndCtx.rndgamma1_uf[nLocusIdx] =   RndCtx.rndgamma1_p[nLocusIdx]
                                     * pow( small/RndCtx.rndgamma1_a[nLocusIdx],
                                            s );
    RndCtx.rndgamma1_d[nLocusIdx] =   RndCtx.rndgamma1_a[nLocusIdx]
                                    * log(RndCtx.rndgamma1_a[nLocusIdx]);
    RndCtx.rndgamma1_ss[nLocusIdx] = s;
  }
  while( 1 )
  {
    r = rndu( nLocusIdx );
    if( r > RndCtx.rndgamma1_p[nLocusIdx] )
    {
      x =   RndCtx.rndgamma1_a[nLocusIdx]
          - log((1 - r) / (1 - RndCtx.rndgamma1_p[nLocusIdx]));
      w =   RndCtx.rndgamma1_a[nLocusIdx] * log(x)
          - RndCtx.rndgamma1_d[nLocusIdx];
    }
    else if( r > RndCtx.rndgamma1_uf[nLocusIdx] )
    {
      x =   RndCtx.rndgamma1_a[nLocusIdx]
          * pow(r / RndCtx.rndgamma1_p[nLocusIdx], 1 / s);
      w = x;
    }
    else
      return (0);

    r = rndu( nLocusIdx );
    if( (1. - r) <= w && r > 0. )
      if( r * ( w + 1 ) >= 1 || -log( r ) <= w )
        continue;
    break;
  }
  return x;
}

//-----------------------------------------------------------------------------
// random standard gamma for s>1
//   Best's (1978) t distribution method

double rndgamma2( int nLocusIdx, double s )
{
  double r,d,f,g,x;
  if( s != RndCtx.rndgamma2_ss[nLocusIdx])
  {
    RndCtx.rndgamma2_b[nLocusIdx]  = s-1;
    RndCtx.rndgamma2_h[nLocusIdx]  = sqrt(3*s-0.75);
    RndCtx.rndgamma2_ss[nLocusIdx] = s;
  }
  while( 1 )
  {
    r = rndu( nLocusIdx );
    g = r - r * r;
    f = (r -0.5) * RndCtx.rndgamma2_h[nLocusIdx]/sqrt(g);
    x = RndCtx.rndgamma2_b[nLocusIdx] + f;
    if (x <= 0)
      continue;
    r=rndu( nLocusIdx );
    d = 64 * r * r * g * g * g;
    if(    d * x < x - 2 * f * f
        || log(d) < 2 * ( RndCtx.rndgamma2_b[nLocusIdx]
                          * log ( x / RndCtx.rndgamma2_b[nLocusIdx]) - f ) )
      break;
  }
  return x;
}
//-----------------------------------------------------------------------------

/*
  double PointChi2 (double prob, double v)
  {
.  returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
  returns -1 if in error.   0.000002<prob<0.999998
  RATNEST FORTRAN by
  Best DJ & Roberts DE (1975) The percentage points of the 
  Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
  Converted into C by Ziheng Yang, Oct. 1993.

  double e=.5e-6, aa=.6931471805, p=prob, g, small=1e-6;
  double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

  if (p<small)   return(0);
  if (p>1-small) return(9999);
  if (v<=0)      return (-1);

  g = LnGamma (v/2);
  xx=v/2;   c=xx-1;
  if (v >= -1.24*log(p)) goto l1;

  ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
  if (ch-e<0) return (ch);
  goto l4;
  l1:
  if (v>.32) goto l3;
  ch=0.4;   a=log(1-p);
  l2:
  q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
  t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
  ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
  if (fabs(q/ch-1)-.01 <= 0) goto l4;
  else                       goto l2;
  
  l3: 
  x=InverseCDFNormal(p);
  p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
  if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
  l4:
  q=ch;   p1=.5*ch;
  if ((t=IncompleteGamma (p1, xx, g))<0)
  printf ("\nIncompleteGamma!!!\n\n");
  p2=p-t;
  t=p2*exp(xx*aa+g+p1-c*log(ch));   
  b=t/ch;  a=0.5*t-b*c;

  s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
  s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
  s3=(210+a*(462+a*(707+932*a)))/2520;
  s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
  s5=(84+264*a+c*(175+606*a))/2520;
  s6=(120+c*(346+127*c))/5040;
  ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
  if (fabs(q/ch-1) > e) goto l4;

  return (ch);
  }
*/


//-----------------------------------------------------------------------------
void flushLine(FILE* readFile)
{
  static char restOfLine[16000] = {'\0'};
  char* p_res = fgets(restOfLine, 16000, readFile);
  if(NULL == p_res)
	// Just to please the compiler. We might do some error logging here.
    return;
}

/*-----------------------------------------------------------------------------
 * Comment Aware version of strtok
 * If string retrieved starts with '#' it is a comment
 * and this function returns NULL
 ----------------------------------------------------------------------------*/
char *strtokCS( char * str, const char * delimiters)
{
  char * result;
  result = strtok(str, delimiters);
  if (result == NULL || result[0] == '#')
    return NULL;
  else
    return result;
}

//============================= END OF FILE ===================================
