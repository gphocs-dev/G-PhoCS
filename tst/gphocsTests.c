#include <stdio.h>
#include <stdlib.h>


#define TRUE 1
#define FALSE 0

int areFilesEqual(char fname1[40], char fname2[40]) {
	FILE *fp1, *fp2;
	int ch1, ch2, result;

	fp1 = fopen(fname1, "r");
	fp2 = fopen(fname2, "r");
	if (fp1 == NULL) {
		fprintf(stderr , "Cannot open %s for reading\n", fname1);
		result = FALSE;
	} else if (fp2 == NULL) {
		fprintf(stderr, "Cannot open %s for reading\n", fname2);
		result = FALSE;
	} else {
		ch1 = getc(fp1);
		ch2 = getc(fp2);
		while ((ch1 != EOF) && (ch2 != EOF) && (ch1 == ch2)) {
			ch1 = getc(fp1);
			ch2 = getc(fp2);
		}
		if (ch1 == ch2){
			result = TRUE;
		}
		else if (ch1 != ch2){
			result = FALSE;
		}
		fclose(fp1);
		fclose(fp2);
	}
	return result;
}


int main(int argc, char* argv[]) {

	//TODO - delete old test files!!

	system("bin/G-PhoCS-1-2-3.exe tst/test-control-file.ctl");


	printf("\n");
	printf("===============================================\n");
	printf("=======             TESTING             =======\n");
	printf("===============================================\n\n");

	printf("1) Compare Trace files:\t\t");
	if (areFilesEqual("./tst/logs/expected/sample-data.trace.tsv", "./tst/logs/test-data.trace.tsv")){
		printf("SUCCESS\n");
	} else {
		printf("FAILURE\n");
	}

	printf("2) Compare FlatStats files:\t");
	if (areFilesEqual("./tst/logs/expected/sample-data.flatstats.tsv", "./tst/logs/test-data.flatstats.tsv")){
		printf("SUCCESS\n");
	} else {
		printf("FAILURE\n");
	}

	return 0;
}
