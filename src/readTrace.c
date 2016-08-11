/** 
    \file readTrace.c
    Post-run analysis of trace output file
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <unistd.h>
#include <limits.h>
#include <getopt.h>


static struct option long_options[] =
  {
    /* These options don't set a flag.
       We distinguish them by their indices. */
    {"block-size",     required_argument,       0, 'b'},
    {"discard",  required_argument,  0, 'd'},
    {"sub-sampling",  required_argument, 0, 's'},
    {"help",  no_argument, 0, 'h'},
    {0, 0, 0, 0}
  };

void printHelp() {
  printf("-b, --block-size  SIZE     Blocksize\n");
  printf("-d, --discard  NUMBER      Number of samples from to discard from beginning of file\n");
  printf("-h, --help                 This help page\n");
//  printf("-s, --subsample  NUMBER    Subsampling, sample every x lines\n");
//  printf("\nReport bugs to <GPhoCS-help-L@cornell.edu>\n");
}

void printUsage(char *filename) {
  printf("Usage: %s <trace-file-name> [options]\n", filename);
  printHelp();
}

int main (int argc, char*argv[]) {
  FILE *traceFile = NULL;
  const int bufferLen = 4096;
  char **colNames;
  char line[bufferLen], *token,  tempStr[bufferLen];
  long double *sums;
  float temp;
  char delims[] = " \t";
  int numCols, ancestorNum, i, numLines, count, blockCount;
  int discardXFromBeginning=0;
  int blockSize=-1; 
  char valueLine[bufferLen];
  char titleLine[bufferLen];
  char formatStr[bufferLen];
  double **data;

  int *width;
  
  int option_index;
  int c;

  opterr = 0;
     
  while (1)
    {
      option_index  = 0;
      c = getopt_long(argc, argv, "b:d:s:ht:", long_options, &option_index);

      if (c == -1)
        break;
           
      switch (c)
        {
        case 'b': //Block size
          blockSize = atoi(optarg);
          break;
        case 'd': //Discard # samples from beginning
          discardXFromBeginning = atoi(optarg);
          break;
//        case 's': //Subsampling
//          subSample = atoi(optarg);
//          break;
        case 'h':
          printUsage(argv[0]);
          printHelp();
          return 0;
          break;
        case '?':
          if ((optopt == 'b') || (optopt == 'd') || (optopt == 's'))
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          else if (isprint (optopt))
            fprintf (stderr, "Unknown option `-%c'.\n", optopt);
          else
            fprintf (stderr,
                     "Unknown option character `\\x%x'.\n",
                     optopt);
          return 1;
        default:
          abort ();
        }
    }
     
  
  if(argv[optind] == NULL) {
    fprintf(stderr, "Missing trace filename.\n");
    printUsage(argv[0]);
    return 1;
  }

  //Open trace file
  traceFile = fopen(argv[optind], "r");
  //Verify trace file was opened successfully
  if(traceFile == NULL) {
    fprintf(stderr, "Could not find trace file '%s' specified.\n", argv[optind]);
    printUsage(argv[0]);
    return 1;
  }

  numLines = 0; //Determine number of lines in the file
  fgets(line, bufferLen, traceFile);	// discard header
  while(!feof(traceFile)) {
    fgets(line, bufferLen, traceFile);
    numLines++;
  } 
  //If user didn't specify a block size, then set the block size to be the number of lines in the file
  if (blockSize < 0) {
    blockSize = numLines;
  }
  if(discardXFromBeginning >= numLines) {
    fprintf(stderr, "%d lines specified to discard, but trace file contains only %d lines.\n", discardXFromBeginning , numLines);
    return 1;
  }
  fseek(traceFile, 0, SEEK_SET);
  //Get first line of trace file
  fgets(line, bufferLen, traceFile);
  strcpy(tempStr,  line);  

  //Read each entry in the first line of the trace file
  //  to determine the number of columns
  token = strtok(line, delims);
  numCols = 0;
  while (token != NULL) {
//	  printf("column %s.\n",token);
    numCols++;
	token = strtok(NULL, delims);
  }
  // remove first column from considerations
  numCols--;   

//	printf("Reading trace for file %s with %d columns.\n",argv[optind], numCols);


  width = (int*)malloc(sizeof(int) * numCols);
  //Allocate space to hold each float from the file
  data = (double**)malloc(sizeof(double*) * (((numLines-1)/blockSize)+1));
  for(i=0;i<=((numLines-1)/blockSize);i++) {
    data[i] = (double*)malloc(sizeof(double) * numCols);
  }
  //Allocate space for column names
  colNames = (char**)malloc(sizeof(char*) * numCols);  
  //Allocate memory for sums of numbers read in from file
  sums = (long double*)malloc(sizeof(long double) * numCols);

  if(tempStr[strlen(tempStr)-1] == '\n') //Strip newline off title line
    tempStr[strlen(tempStr)-1] = ' ';

  //Parse each column name (header) from the line we read in
  token = strtok(tempStr, delims);
  token = strtok(NULL, delims);
  for(i=0;i<numCols;i++) {
    colNames[i] = malloc(sizeof(char) * (strlen(token)+1));
    strcpy(colNames[i],token);
	token = strtok(NULL, delims);
  }   

  
  //Discard specified number of samples
  for(i=0;i<discardXFromBeginning;i++)
    fgets(line, bufferLen, traceFile);
  
  //For each line in trace file
  count = 0;
  blockCount = 0;
  while (1) {
    //Get the next line from the file
    fgets(line, bufferLen, traceFile);
    
    //If we are at the end of the file stop
    if (feof(traceFile))
      break;

    token = strtok(line, delims);
    token = strtok(NULL, delims);
    count++;
    for(ancestorNum=0; ancestorNum < numCols; ancestorNum++) {
      //Convert next entry to float
      sscanf(token, "%f", &temp);
      //Add newly read value to column sum
      sums[ancestorNum] = sums[ancestorNum] + temp;  
      if(count == blockSize) {
        //Format the float as a string
        sprintf(tempStr, "%.6Lf    ", sums[ancestorNum] / blockSize); 
        //If this new float (in characters) is longer than any previous, update the width this column will display at
        if(strlen(tempStr) > width[ancestorNum])
          width[ancestorNum] = strlen(tempStr);
        //Save the average we just computed to memory so we can display it later once we know all column widths
        data[blockCount][ancestorNum] = sums[ancestorNum] / blockSize;
        //Reset the sums
        sums[ancestorNum] = 0;
      }
      
      //read in next entry
      token = strtok(NULL, delims);
    }
    if(count == blockSize) {
      count = 0;
      blockCount++;
    }
    
    //Update line number
  }

  //All trace file in memory
  fclose(traceFile); 

  if(count > 0) {
    for(ancestorNum=0; ancestorNum < numCols; ancestorNum++) {
      //Format the float as a string
      sprintf(tempStr, "%.6Lf    ", sums[ancestorNum] / count); 
      //If this new float (in characters) is longer than any previous, update the width this column will display at
      if(strlen(tempStr) > width[ancestorNum]) {
        width[ancestorNum] = strlen(tempStr);
        //Save the average we just computed to memory so we can display it later once we know all column widths
        data[blockCount][ancestorNum] = sums[ancestorNum] / count;
      }
    }
    blockCount++;
  }

  sprintf(valueLine, "%s", "");
  sprintf(titleLine, "%s", "");
  //Print results to screen
    for(i=0; i < blockCount; i++) {
      for(ancestorNum=0; ancestorNum < numCols; ancestorNum++) {   
        //Just as before, format the float as a string w/ precision of 6
        sprintf(tempStr, "%.6f    ", data[i][ancestorNum]); 
        //Create the format string that takes the column width into account
        sprintf(formatStr, "%%s%%-%ds", width[ancestorNum]);
        //Add the float to the string to be printed formatted for the right width
        sprintf(valueLine, formatStr, valueLine, tempStr);
        //Add the title to the string to be printed for the right width
        if(i==0) //Only create the title once
          sprintf(titleLine, formatStr, titleLine, colNames[ancestorNum]);
        if((((ancestorNum+1) % 90) == 0) || ancestorNum == numCols-1) {
          if(i==0)//Only print the title on the first line
            printf("%s\n", titleLine);
          printf("%s\n", valueLine); //Print the line of formatted values
          sprintf(titleLine, "%s", ""); //Clear the titleLine & valueLine strings
          sprintf(valueLine, "%s", "");
        }
      }
    }


  printf("\n");


  //Completed successfully
  return 0;
}
