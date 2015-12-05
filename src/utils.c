/** 
   \file utils.c
   Utility definitions & functions (miscellaneous).
	
   Contains some utility functions and definitions.
	
*/

#include "utils.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include <limits.h>
#include <float.h>
#include <time.h>


double rndgamma1 (double s);
double rndgamma2 (double s);

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


void test(const char* message){
  printf("test %s\n",message);
  return;
}


void resetBooleanArray1(unsigned short* booleanArray, int arrayLength){
  int i;
  for(i=0; i<arrayLength; i++) {
    booleanArray[i] = 0;
  }
  return;
}

void resetBooleanArray(unsigned short* booleanArray, int arrayLength){
  static unsigned short zeros[200] = {0};

  if(arrayLength > 200) {
    resetBooleanArray1(booleanArray,arrayLength);
  } else {
    memcpy(booleanArray, zeros, arrayLength*sizeof(unsigned short));
  }
  return;
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

double reflect(double x, double a, double b)
{
  /* This returns a variable in the range (a,b) by reflecting x back into the range
   */

  static double slack = 0.000000001;		// safety margins for upper and lower bounds
  double xnew, double_interval;
	
  a += slack;
  b -= slack;
	
  // if interval is empty (due to slackness), return middle of interval
  if(b<=a) {
//		fprintf(stderr, "very small interval in reflect(%g,%g,%g) [slackness = %g].\n",x,a-slack,b+slack,slack);
    return (a+b)/2;
  }
	
  if(x<b && x>a)
    return x;
	
  // reflect upwards, if necessary
  xnew = x;
  if(xnew<=a)
    xnew = 2*a - xnew;
		
  // "fold twice" as many time as needed
  double_interval = 2*(b-a);
  xnew = xnew - double_interval*floor( (xnew-a) / double_interval );
	
  // reflect downwards one last time, if necessary 
  if(xnew>=b)
    xnew = 2*b - xnew;
	
  // value should be within interval at this stage, but numerical percision might put it slightly outside
  while(xnew <=a || xnew >= b) {
	if(xnew>=b) {
		if(debug) {
			fprintf(stderr, "reflect percision in reflect(%g,%g,%g): obtaining reflection at %g (%g greater than upper bound)\n",x,a-slack,b+slack,xnew, xnew-b);
		}
    	xnew = 2*b - xnew;
	} else {
		if(debug) {
			fprintf(stderr, "reflect percision in reflect(%g,%g,%g): obtaining reflection at %g (%g smaller than lower bound)\n",x,a-slack,b+slack,xnew, a-xnew);
		}
    	xnew = 2*a - xnew;
	}
  }
		
	  
  return xnew;
	
}


// random sampling functions

static unsigned int z_rndu=137;
static unsigned int w_rndu=123456757;
static double m2s2_kernel=8, m2N_kernel, s2N_kernel;;

void setSeed (unsigned int seed) {
  if(sizeof(int) != 4) 
    puts("oh-oh, we are in trouble.  int not 32-bit?");
  z_rndu=170*(seed%178)+137;
  w_rndu = seed*127773;

  m2N_kernel = sqrt(m2s2_kernel/(m2s2_kernel+1));  
  s2N_kernel = sqrt(1/(m2s2_kernel+1));
}

double rndnormal (void)
{
  /* standard normal variate, using the Box-Muller method (1958), improved by 
     Marsaglia and Bray (1964).  The method generates a pair of random
     variates, and only one used.
     See N. L. Johnson et al. (1994), Continuous univariate distributions, 
     vol 1. p.153.
  */
  double u, v, s;

  for (; ;) {
    u = 2*rndu() - 1;
    v = 2*rndu() - 1;
    s = u*u + v*v;
    if (s>0 && s<1) break;
  }
  s = sqrt(-2*log(s)/s);
  return (u*s);
}


double rnd2normal8 (void)
{
  /* This returns a variate from the mixture of two normals
     N(-m, s2) and N(m, s2), with mean 0 and variance m^2 + s2 = 1 and m^2/s^2 = 8.

     Let this standard variate be z.  Then mean + z * sigma will be a variate 
     with mean mean and SD sigma.  This is useful for generating MCMC proposals
  */
  double z = m2N_kernel + rndnormal()*s2N_kernel;
  ;
  z = (rndu()<0.5 ? z : -z);
  return (z);
}


double rndu (void)
{
  /* U(0,1): AS 183: Appl. Stat. 31:188-190 
     Wichmann BA & Hill ID.  1982.  An efficient and portable
     pseudo-random number generator.  Appl. Stat. 31:188-190

     x, y, z are any numbers in the range 1-30000.  Integer operation up
     to 30323 required.
  */
  static unsigned int x_rndu=11, y_rndu=23;
  double r;

  x_rndu = 171*(x_rndu%177) -  2*(x_rndu/177);
  y_rndu = 172*(y_rndu%176) - 35*(y_rndu/176);
  z_rndu = 170*(z_rndu%178) - 63*(z_rndu/178);
  /*
    if (x_rndu<0) x_rndu+=30269;
    if (y_rndu<0) y_rndu+=30307;
    if (z_rndu<0) z_rndu+=30323;
  */
  r = x_rndu/30269.0 + y_rndu/30307.0 + z_rndu/30323.0;
  return (r-(int)r);
}



double rndgamma (double s)
{
  /* random standard gamma (Mean=Var=s,  with shape parameter=s, scale para=1)
     r^(s-1)*exp(-r)
     J. Dagpunar (1988) Principles of random variate generation,
     Clarendon Press, Oxford
     calling rndgamma1() if s<1 or
     rndgamma2() if s>1 or
     exponential if s=1

     This is unsafe, and is found to return 0 when s is very small.
  */
  double r=0;

  if (s<=0)      puts ("jgl gamma..");
  else if (s<1)  r=rndgamma1 (s);
  else if (s>1)  r=rndgamma2 (s);
  else           r=-log(rndu());
  return (r);
}


double rndgamma1 (double s)
{
  /* random standard gamma for s<1
     switching method
  */
  double r, x=0,small=1e-37,w;
  static double a,p,uf,ss=10,d;

  if (s!=ss) {
    a=1-s;
    p=a/(a+s*exp(-a));
    uf=p*pow(small/a,s);
    d=a*log(a);
    ss=s;
  }
  for (;;) {
    r=rndu();
    if (r>p)        x=a-log((1-r)/(1-p)), w=a*log(x)-d;
    else if (r>uf)  x=a*pow(r/p,1/s), w=x;
    else            return (0);
    r=rndu ();
    if (1-r<=w && r>0)
      if (r*(w+1)>=1 || -log(r)<=w)  continue;
    break;
  }
  return (x);
}

double rndgamma2 (double s)
{
  /* random standard gamma for s>1
     Best's (1978) t distribution method
  */
  double r,d,f,g,x;
  static double b,h,ss=0;
  if (s!=ss) {
    b=s-1;
    h=sqrt(3*s-0.75);
    ss=s;
  }
  for (;;) {
    r=rndu ();
    g=r-r*r;
    f=(r-0.5)*h/sqrt(g);
    x=b+f;
    if (x <= 0) continue;
    r=rndu();
    d=64*r*r*g*g*g;
    if (d*x < x-2*f*f || log(d) < 2*(b*log(x/b)-f))  break;
  }
  return (x);
}



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



void flushLine(FILE* readFile) {
  static char restOfLine[16000]={'\0'};
  // char* res;
  fgets(restOfLine,16000,readFile);
  return;
}


  /***********************************************************************************
   *    Comment Aware version of strtok
   *    -If string retrieved starts with '#' it is a comment and this function returns NULL
   ***********************************************************************************/
  char *strtokCS( char * str, const char * delimiters) {
    char * result;
    result = strtok(str, delimiters);
    if (result == NULL || result[0] == '#')
      return NULL;
    else
      return result;
  }
