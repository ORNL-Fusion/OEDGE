#include <stdio.h>
#include <stdlib.h>
#include <string.h>


char input[50],target[50],output[50];
char marker[50],chop[50],replace[50];

int code;

FILE *fpinput;

void getinputdata(void) {

  fgets(marker ,50,fpinput);
  fgets(chop   ,50,fpinput);
  if (fgets(replace,50,fpinput)==NULL) code = 1;

  marker [strlen(marker )-1] = NULL;
  chop   [strlen(chop   )-1] = NULL;
  replace[strlen(replace)-1] = NULL;

  if (code==0) {
    printf("Marker : \"%s\"\n",marker);
    printf("Chop   : \"%s\"\n",chop);
    printf("Replace: \"%s\"\n",replace);
  }

}




int idline(char *line) {

  int i;

  if (code==1) return NULL;

  for (i=0;i<strlen(marker);i++) {
    if (marker[i]!=line[i]) return NULL;
  }

  return 1;
}

int findchop(char *line) {

  int i,j;
  
  for (i=0;i<strlen(line)-strlen(chop);i++) {
    for (j=0;j<strlen(chop);j++) if (line[i+j]!=chop[j]) break;
    if (j==strlen(chop)) return i;
  }

  return -1;
}


void splicefile(void) {

  int i,ind;
  char line[200],lineout[200];
  FILE *fpout,*fp;

  if ((fp=fopen(target,"r"))==NULL) {
    printf("Can't find \"%s\".\n",target);
    exit(1);
  }

  fpout = fopen(output,"w");

  getinputdata();

  while (fgets(line,200,fp)!=NULL) {
    if (idline(line)==1) {
      printf("Line   : %s",line);

      ind = findchop(line);

      if (ind==-1) {
	printf("Can't find chop string.\n");
	exit(1);
      }

      for (i=0;i<ind;i++) lineout[i] = line[i];

      for (i=ind;i<ind+strlen(replace);i++) 
	lineout[i] = replace[i-ind];

      for (i=ind+strlen(chop);i<strlen(line);i++)
	lineout[i+strlen(replace)-strlen(chop)] = line[i];

      lineout[strlen(line)-strlen(chop)+strlen(replace)] = NULL;

      printf("Lineout: %s",lineout);

      fprintf(fpout,"%s",lineout);

      getinputdata();
    }
    else {
      fprintf(fpout,"%s",line);
    }
  }

  fclose(fp);
  fclose(fpout);

}


void main(int argc, char **argv) {

  if (argc!=4) {
    printf("Sorry, need [input][target][dest].\n");
    exit(1);
  }

  code = 0;

  printf("\n");

  sprintf(input ,"%s",argv[1]);
  sprintf(target,"%s",argv[2]);
  sprintf(output,"%s",argv[3]);

  printf("Input  : \"%s\"\n",input);
  printf("Target : \"%s\"\n",target);
  printf("Output : \"%s\"\n",output);

  if ((fpinput=fopen(input,"r"))==NULL) {
    printf("Can't find \"%s\".\n",input);
    exit(1);
  }

  splicefile();

  exit (1);
}
