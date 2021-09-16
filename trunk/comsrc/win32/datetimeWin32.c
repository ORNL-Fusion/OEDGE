/* C Subroutines to act as interface for FORTRAN to the 
   C run-time library */

#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>

int DATETEST(char *string)
{
	time_t tloc;
	struct tm *tres;
	int  day,month,year;

	time(&tloc);

	tres = localtime(&tloc);

	day = (*tres).tm_mday ;
	month = (*tres).tm_mon+1;
	year = (*tres).tm_year;
	if (year>99) { year = year -100; };

	sprintf(string,"%2.2d/%2.2d/%2.2d",month,day,year);

	return(1);
}

int CLOCKTEST(char *string)
{
	time_t tloc;
	time(&tloc);
	strncpy(string,ctime(&tloc)+11,8);

	return 1;
}
/*
int setenv(char *EnvFile)
{
	int Res = 0;
	int len = 0;
	FILE* fptr = fopen(EnvFile, "rt");
	if(fptr != NULL)
	{
		char CurLine[256];
		while(fgets(CurLine, 256, fptr) != NULL)
		{
			printf("SetEnv(): %\n", CurLine);
			len = (int)strlen(CurLine);
			CurLine[len - 1] = '\0';
			Res = _putenv(CurLine);
			if(Res < 0)
			{
				printf("SetEnv(): _putenv() failed\n");
				break;
			}
		}
		fclose(fptr);
	}
	
	return Res;
}
*/


