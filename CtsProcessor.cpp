
/*
 *	The usage of the program is: 
 *	./CtsProcessor.out A.gjf B.gjf
 *	The program processes (rotate and translate) the cartesian coordinates found in A.gjf 
 *	in such way that the processed coordinates are closest to those in B.gjf.
 *	Then the program gives the processed coordinates and standard deviation.
 *	
 *	FORMAT in *.gjf:
 *	ANY line that scanf("%s%lf%lf%lf") can read 4 vars (1 str, 3 double) are processed. Other lines are DISCARDED.
 *	Coordinates in A.gjf and B.gjf MUST be in the same sequence.
 */

#include "stdio.h"
#include "math.h"
#include "string.h"
#include "Optimization.h"

const double pi=3.1415926536;
const int NAtomMax=300;	// A lower number of actual atoms is tolerated
int nAtom;

typedef double TransFormMatrix[3][3];

class AtomCartesian
{
public:
	double x[3];
	char Element[5];

	AtomCartesian() {x[0]=x[1]=x[2]=0, Element[0]=0;}

	int ReadFromLine(const char *Line)
	{
		return sscanf(Line,"%s%lf%lf%lf",Element,x,x+1,x+2)==4;
	}

	AtomCartesian operator * (TransFormMatrix T)
	{
		AtomCartesian z;
		z=*this;
		z.x[0]=z.x[1]=z.x[2]=0;
		int i,j;
		for(i=0;i<3;i++) for(j=0;j<3;j++) z.x[i]+=T[i][j]*x[j];
		return z;
	}

	AtomCartesian operator -= (AtomCartesian b)
	{
		for(int i=0;i<3;i++) x[i]-=b.x[i];
		return *this;
	}

	AtomCartesian operator - (AtomCartesian b)
	{
		AtomCartesian z=*this;
		z-=b;
		return z;
	}

	double SqrMod()
	{
		double z=0;
		int i,j;
		for(j=0;j<3;j++) z+=x[j]*x[j];
		return z;
	}

	void print()
	{
		printf("%s  %.6lf  %.6lf  %.6lf\n",Element,x[0],x[1],x[2]);
	}

} AtomInitList[NAtomMax],AtomRotatedList[NAtomMax],AtomTargetList[NAtomMax];

void CalcUnitary(double fai, double theda, double chai, TransFormMatrix Unitary)
{	//	the matrix of rotation
	double cf=cos(fai),sf=sin(fai),
		ct=cos(theda),st=sin(theda),
		cc=cos(chai),sc=sin(chai);
	Unitary[0][0]=cf*ct*cc-sf*sc;
	Unitary[0][1]=sf*ct*cc+cf*sc;
	Unitary[0][2]=-st*cc;
	Unitary[1][0]=-cf*ct*sc-sf*cc;
	Unitary[1][1]=cf*cc-sf*ct*sc;
	Unitary[1][2]=st*sc;
	Unitary[2][0]=cf*st;
	Unitary[2][1]=sf*st;
	Unitary[2][2]=ct;
}

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

int ReadAtomListFromFile(AtomCartesian *AtomList, const char *fname)
{
	FILE *f=fopen(fname,"r");
	if (f==NULL)
	{
		printf("Cannot open %s\n",fname);
		return 0;
	}
	char Line[2000]; int EOLN,i=0;
	while (!feof(f))
	{
		ReadLine(f,EOLN,Line);
		if (AtomList[i].ReadFromLine(Line)) i++;
	}

	printf("Found %d atoms in %s\n",i,fname);
	fclose(f);
	return i;
}

/*
void PrintUnitary()
{
	int i,j,k;
	double sum;
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++) printf("%.4lf\t",Unitary[i][j]);
		puts("");
	}
	puts("");
	for (i=0;i<3;i++)
	{
		for (j=0;j<3;j++)
		{
			sum=0;
			for(k=0;k<3;k++) sum+=Unitary[i][k]*Unitary[j][k];	// = check if U is unitary
			printf("%.4lf\t",sum);
		}
		puts("");
	}
	puts("");
}
*/

double ErrorFunc(double *v, int n)
{	// args: double fai,double theda,double chai, double xOffset, double yOffset, double zOffset
	TransFormMatrix T;
	CalcUnitary(v[0],v[1],v[2],T);
	double z=0;
	int i,j,k;
	AtomCartesian temp;

	for (i=0;i<nAtom;i++)
	{
		temp=AtomInitList[i]*T;
		temp.x[0]+=v[3];
		temp.x[1]+=v[4];
		temp.x[2]+=v[5];
		AtomRotatedList[i]=temp;
		temp-=AtomTargetList[i];
		z+=temp.SqrMod();
		//	Weight
		if (strcmp(AtomInitList[i].Element,"W")==0 || strcmp(AtomInitList[i].Element,"Mo")==0)
			z+=40000*temp.SqrMod();
		if (strcmp(AtomInitList[i].Element,"C")==0 || strcmp(AtomInitList[i].Element,"O")==0)
			z+=100*temp.SqrMod();
	}
	return z;
}

double LooseErrorFunc(double *v, int n)
{       // args: double fai,double theda,double chai, double xOffset, double yOffset, double zOffset
        TransFormMatrix T;
        CalcUnitary(v[0],v[1],v[2],T);
        double z=0;
        int i,j,k,nAtomUsed;
        AtomCartesian temp;
	nAtomUsed=6;
	if (nAtomUsed>nAtom) nAtomUsed=nAtom;

        for (i=0;i<nAtomUsed;i++)
        {
                temp=AtomInitList[i]*T;
                temp.x[0]+=v[3];
                temp.x[1]+=v[4];
                temp.x[2]+=v[5];
                AtomRotatedList[i]=temp;
                temp-=AtomTargetList[i];
                z+=temp.SqrMod();
		//	Weight
		if (strcmp(AtomInitList[i].Element,"W")==0 || strcmp(AtomInitList[i].Element,"Mo")==0) z+=100*temp.SqrMod();
        }
        return z;
}

int main(int argc, char *argv[])
{
	int i,j;
	if((nAtom=ReadAtomListFromFile(AtomInitList,argv[1]))!=ReadAtomListFromFile(AtomTargetList,argv[2]))
	{
		puts("Error: Different number of atoms found");
		return 0;
	}
	for (i=0;i<nAtom;i++) if(strcmp(AtomInitList[i].Element,AtomTargetList[i].Element))
	{
		printf("Atom %d found different: %s in %s but %s in %s\n",i+1,AtomInitList[i].Element,argv[1],AtomTargetList[i].Element,argv[2]);
		return 0;
	}

	double OptArg[]={0,0,0,0,0,0};
	double OptStep[]={0.5,0.5,0.5,1,1,1};
	double stdev=Opt(LooseErrorFunc,6,OptArg,OptStep);	//	Pre-optimization
	for(j=0;j<6;j++) OptStep[j]*=200;
	stdev=Opt(ErrorFunc,6,OptArg,OptStep);
	for (i=0;i<nAtom;i++) AtomRotatedList[i].print();
	stdev=sqrt(stdev/i);
	printf("stdev=%.6lf\n",stdev);
	return 0;
}


