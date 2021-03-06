#include "stdio.h"
#include "string.h"
#include "math.h"
#include "vector"
#include "iostream"
using namespace std;

#define PRINT_LOWEST_ENERGY printf("\033[33mLowest Energy: AbsE = %.6f Hartree, RelaE = %.2f kcal/mol in step %d\033[0m\n",E_lowest,627.5095*(E_lowest-Energy_Level),i_min_E)
#define PRINT_LOWEST_MAX_F printf("Minimum MAX_F = %.6f in step %d\t\t",min_MAX_F,i_min_MAX_F)
#define PRINT_LOWEST_RMS_F printf("Minimum RMS_F = %.6f in step %d\n",min_RMS_F,i_min_RMS_F)
#define PRINT_LOWEST_MAX_D printf("Minimum MAX_D = %.6f in step %d\t\t",min_MAX_D,i_min_MAX_D)
#define PRINT_LOWEST_RMS_D printf("Minimum RMS_D = %.6f in step %d\n",min_RMS_D,i_min_RMS_D)
#define PRINT_LOWEST if (!LowestEnergyPrinted) {LowestEnergyPrinted=1; PRINT_LOWEST_ENERGY; PRINT_LOWEST_MAX_F; PRINT_LOWEST_RMS_F; PRINT_LOWEST_MAX_D; PRINT_LOWEST_RMS_D; }

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
	sprintf(Ctrl,"\t\033[%sm%%.6f  %%.6f\033[0m",p);
	printf(Ctrl,x,tol);
}	//https://blog.csdn.net/virtual_func/article/details/47707601

void PrintScanList(const vector <double> ScanList, const vector < const char * > ScanEWList)
{
	int i,ScanL;
	if ((ScanL=ScanList.size())>1)
	{
		puts("\nScan summary:");
		for (i=0;i<ScanL;i++)
			printf("%d\t AbsE = %.6f Hartree, RelaE = %.2f kcal/mol",i,ScanList[i],627.5095*(ScanList[i]-ScanList[0])),
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
	int i_min_E,i_min_MAX_F,i_min_RMS_F,i_min_MAX_D,i_min_RMS_D;
	int EOLN,i,N=0,flag=1;	// flag=1 indicates the new starting of an opt step or what.
	int MoreThanSCF=0;
	char *p,*p2;
	double MAX_F,TOL_MAX_F,min_MAX_F;
	double RMS_F,TOL_RMS_F,min_RMS_F;
	double MAX_D,TOL_MAX_D,min_MAX_D;
	double RMS_D,TOL_RMS_D,min_RMS_D;
	double Energy,Energy_Level,E_lowest=1000;
	vector <double> ScanList;
	vector < const char * > ScanEWList;

	i=1;
	printf("\n  \033[33mAbsE(Hartree)\tRelaE(kcal/mol)\033[0m\nn\tMAX_Force  Threshold\tRMS_Force  Threshold\tMAX_Dsplsmnt  Threshold\tRMS_Dsplsmnt  Threshold\n");
	while (!feof(f))
	{
		ReadLine(f,EOLN,Line);
		p=strstr(Line,"Converged?");
		if (p)
		{
			ReadLine(f,EOLN,Line);
			if (sscanf(strstr(Line,".")-1,"%lf%lf",&MAX_F,&TOL_MAX_F)!=2) continue;
			ReadLine(f,EOLN,Line);
			if (sscanf(strstr(Line,".")-1,"%lf%lf",&RMS_F,&TOL_RMS_F)!=2) continue;
			ReadLine(f,EOLN,Line);
			if (sscanf(strstr(Line,".")-1,"%lf%lf",&MAX_D,&TOL_MAX_D)!=2) continue;
			ReadLine(f,EOLN,Line);
			if (sscanf(strstr(Line,".")-1,"%lf%lf",&RMS_D,&TOL_RMS_D)!=2) continue;
			printf("%d",i);
			DualPrint(MAX_F,TOL_MAX_F,"4","0");
			DualPrint(RMS_F,TOL_RMS_F,"4","0");
			DualPrint(MAX_D,TOL_MAX_D,"4","0");
			DualPrint(RMS_D,TOL_RMS_D,"4","0");
			puts("");

			if (i==1)
			{
				min_MAX_F=MAX_F, min_RMS_F=RMS_F, min_MAX_D=MAX_D, min_RMS_D=RMS_D;
				i_min_MAX_F = i_min_RMS_F = i_min_MAX_D = i_min_RMS_D = i;
			}
			if (min_MAX_F>MAX_F) min_MAX_F=MAX_F, i_min_MAX_F=i;
			if (min_RMS_F>RMS_F) min_RMS_F=RMS_F, i_min_RMS_F=i;
			if (min_MAX_D>MAX_D) min_MAX_D=MAX_D, i_min_MAX_D=i;
			if (min_RMS_D>RMS_D) min_RMS_D=RMS_D, i_min_RMS_D=i;

			i++;
		}
		p=strstr(Line,"SCF Done:");
		if (p && (!MoreThanSCF))
		{
			LowestEnergyPrinted=0;
			if (sscanf(strstr(Line,"=")+1,"%lf",&Energy)!=1) continue;
			if (flag)
			{
				flag=0;
				i_min_E = 1;
				E_lowest=Energy_Level=Energy;	//	This would not interfere a MP2 or double hybrid functional job, since MP2 always gives a negative correction in energy.
			}
			printf("  \033[33m%.6f(SCF)\t%.2f\033[0m\n",Energy,627.5095*(Energy-Energy_Level));
			if (E_lowest>Energy) E_lowest=Energy, i_min_E=i;
		}
		p2=strstr(Line,"E2(");
		if (p2)
		{
			p=strstr(p2,"E(");
			if (!p) continue;
			int Walker=0;
			while (p2[Walker])
			{
				if (p2[Walker]=='(' || p2[Walker]==')') 
					p2[Walker]=' ';	//	to enable the following sscanf
				Walker++;
			}
			char Fktional1[300]="",Fktional2[300]="",EnergyStr[300]="",MP2Corr[300]="";
			if (sscanf(p2,"E2 %s  = %s E %s  = %s",Fktional1,MP2Corr,Fktional2,EnergyStr)!=4) continue;
			if (strcmp(Fktional1,Fktional2)) continue;
			int DtoE=0;
			while(EnergyStr[DtoE])
			{
				if (EnergyStr[DtoE]=='d' || EnergyStr[DtoE]=='D') EnergyStr[DtoE]='E';
				DtoE++;
			}
			if (sscanf(EnergyStr,"%lf",&Energy)!=1) continue;
			if (!MoreThanSCF)
			{
				flag=0;
				i_min_E=1;
				MoreThanSCF=1;
				E_lowest=Energy_Level=Energy;
			}
			printf("  \033[33m%.6f(%s)\t%.2f\033[0m\n",Energy,Fktional1,627.5095*(Energy-Energy_Level));
			if (E_lowest>Energy) E_lowest=Energy, i_min_E=i;
		}
		p=strstr(Line,"Optimization completed");
		if (p)
		{
			if (flag) puts("Optimization Completed.\n"); else
			{
				printf("\033[32mOptimization Completed with ");
				PRINT_LOWEST;
				ScanList.push_back(E_lowest);
				ScanEWList.push_back("");
			}
			flag=1;
			MoreThanSCF=0;
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
				PRINT_LOWEST;
				ScanList.push_back(E_lowest);
				ScanEWList.push_back("Optimization stopped");
			}
			flag=1,i=1;
			MoreThanSCF=0;
			ReadLine(f,EOLN,Line);
			puts(Line);
			ReadLine(f,EOLN,Line);
			puts(Line);
		}
		p=strstr(Line,"Error termination");
		if (p)
		{
			PRINT_LOWEST;
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
			PRINT_LOWEST;
			PrintScanList(ScanList,ScanEWList);
			ScanList.clear();
			i=1;
			MoreThanSCF=0;
		}

	}
	puts("");
	if (flag==0)
	{
		PRINT_LOWEST;
		ScanList.push_back(E_lowest);
		ScanEWList.push_back("Optimization still in progress");
	}
	PrintScanList(ScanList,ScanEWList);
	ScanList.clear();
	fclose(f);
	puts("");
	return 0;
}

