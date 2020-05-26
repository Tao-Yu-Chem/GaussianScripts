#include "stdio.h"
#include "string.h"
#include "math.h"
#include "vector"
#include "iostream"
using namespace std;

#define PRINT_LOWEST_ENERGY if (!LowestEnergyPrinted) LowestEnergyPrinted=1,printf("\033[33mLowest Energy: AbsE = %.6f Hartree, RelaE = %.2f kJ/mol\033[0m\n",E_lowest,2625.5*(E_lowest-Energy_Level))

void ReadLine(FILE *f,int &EOLN, char *Line)
{
	char ch;
	int p=0;
	
	EOLN=0;
	while(1)
	{
		ch=fgetc(f);
		if (ch=='\r') continue;
		if (ch=='\n' || ch==EOF || feof(f)) {EOLN=1; Line[p]=0; break;}
		Line[p++]=ch;
	}
	Line[p]=0;
}

void DualPrint(double x,double tol,const char *ifYes,const char *ifNo)
{
	char Ctrl[50]; const char *p;
	p=(x<tol)?ifYes:ifNo;
	sprintf(Ctrl,"\t\033[%sm%%.6f\t%%.6f\033[0m",p);
	printf(Ctrl,x,tol);
}	//https://blog.csdn.net/virtual_func/article/details/47707601

void PrintScanList(const vector <double> ScanList, const vector < const char * > ScanEWList)
{
	int i,ScanL;
	if ((ScanL=ScanList.size())>1)
	{
		puts("\nScan summary:");
		for (i=0;i<ScanL;i++)
			printf("%d\t AbsE = %.6f Hartree, RelaE = %.2f kJ/mol",i,ScanList[i],2625.5*(ScanList[i]-ScanList[0])),
			printf(" \t\033[31m%s\033[0m\n",ScanEWList[i]);
		puts("");
	}
}

int main(int argc, char *argv[])
{
	FILE *f;
	f=fopen(argv[1],"r");
	if (f==NULL)
	{
		printf("Error reading file with error code %d\n",ferror(f));
		return 0;
	}
	char Line[9999];
	int LowestEnergyPrinted=0;
	int EOLN,i,N=0,flag=1;	// flag=1 indicates the new starting of an opt step or what.
	char *p;
	double MAX_F,TOL_MAX_F;
	double RMS_F,TOL_RMS_F;
	double MAX_D,TOL_MAX_D;
	double RMS_D,TOL_RMS_D;
	double Energy,Energy_Level,E_lowest=1000;
	vector <double> ScanList;
	vector < const char * > ScanEWList;

	i=1;
	printf("\n\t\t\t\t\t\t\t\t\t\t\033[33mAbsE(Hartree)\tRelaE(kJ/mol)\033[0m\nn\tMAX_Force\tThreshold\tRMS_Force\tThreshold\tMAX_Dsplsmnt\tThreshold\tRMS_Dsplsmnt\tThreshold\n");
	while (!feof(f))
	{
		ReadLine(f,EOLN,Line);
		p=strstr(Line,"Converged?");
		if (p)
		{
			ReadLine(f,EOLN,Line);
			sscanf(strstr(Line,".")-1,"%lf%lf",&MAX_F,&TOL_MAX_F);
			ReadLine(f,EOLN,Line);
			sscanf(strstr(Line,".")-1,"%lf%lf",&RMS_F,&TOL_RMS_F);
			ReadLine(f,EOLN,Line);
			sscanf(strstr(Line,".")-1,"%lf%lf",&MAX_D,&TOL_MAX_D);
			ReadLine(f,EOLN,Line);
			sscanf(strstr(Line,".")-1,"%lf%lf",&RMS_D,&TOL_RMS_D);
			printf("%d",i);
			DualPrint(MAX_F,TOL_MAX_F,"4","0");
			DualPrint(RMS_F,TOL_RMS_F,"4","0");
			DualPrint(MAX_D,TOL_MAX_D,"4","0");
			DualPrint(RMS_D,TOL_RMS_D,"4","0");
			puts("");
			i++;
		}
		p=strstr(Line,"SCF Done:");
		if (p)
		{
			LowestEnergyPrinted=0;
			sscanf(strstr(Line,"=")+1,"%lf",&Energy);
			if (flag)
			{
				flag=0;
				E_lowest=Energy_Level=Energy;
			}
			printf("\t\t\t\t\t\t\t\t\t\t\033[33m%.6f\t%.2f\033[0m\n",Energy,2625.5*(Energy-Energy_Level));
			if (E_lowest>Energy) E_lowest=Energy;
		}
		p=strstr(Line,"Optimization completed");
		if (p)
		{
			if (flag) puts("Optimization Completed.\n"); else
			{
				printf("\033[32mOptimization Completed with ");
				PRINT_LOWEST_ENERGY;
				ScanList.push_back(E_lowest);
				ScanEWList.push_back("");
			}
			flag=1;
			i=1;
		}
		p=strstr(Line,"Optimization stopped");
		if (p)
		{
			if (flag)
				printf("\033[31m%s\033[0m\n",p);
			else
			{
				printf("\033[31mOptimization stopped with ");
				PRINT_LOWEST_ENERGY;
				ScanList.push_back(E_lowest);
				ScanEWList.push_back("Optimization stopped");
			}
			flag=1,i=1;
			ReadLine(f,EOLN,Line);
			puts(Line);
			ReadLine(f,EOLN,Line);
			puts(Line);
		}
		p=strstr(Line,"Error termination");
		if (p)
		{
			PRINT_LOWEST_ENERGY;
			printf("\033[31m%s\033[0m\n",p);
			PrintScanList(ScanList,ScanEWList);
			ScanList.clear();
			ScanEWList.clear();
		}
		p=strstr(Line,"Normal termination");
		if (p)
		{
			printf("\033[32m%s\033[0m\n",p);
			flag=1;
			PRINT_LOWEST_ENERGY;
			PrintScanList(ScanList,ScanEWList);
			ScanList.clear();
			i=1;
		}

	}
	puts("");
	if (flag==0)
	{
		PRINT_LOWEST_ENERGY;
		ScanList.push_back(E_lowest);
		ScanEWList.push_back("Optimization still in progress");
	}
	PrintScanList(ScanList,ScanEWList);
	ScanList.clear();
	fclose(f);
	puts("");
	return 0;
}

