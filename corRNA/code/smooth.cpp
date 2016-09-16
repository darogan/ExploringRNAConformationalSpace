/*
c++ smooth.cpp -o smooth sims/util.o sims/geom.o -lm

code/smooth model.pdb pairs.dat N (cycles)

*/
#include "sims/util.hpp"
#include "sims/geom.hpp"

typedef struct {
	char	*res;	// residue type
	float	*acc;	// surface accesssibility
	int	*sec;	// secondary structure type (0=c, 1=A, 2=B, 3=3ten)
	int	*rid;	// residue number
	Vec	*cas;	// CA position (runs 1...N with dummy 0 and N+1)
	Vec	*cbs;	// dummy CB/centroid (extended <bond>A from CA)
	int	len;	// chain length
	int	gaps;	// number of chain breaks
	float	pcta, pctb;	// percent alpha and beta structure
} Prot;

void pdbout ( Prot*, char* ); 
void smooth ( Vec*, int, int );

char	aa3[2222][4];

void addres ( Vec* cas, int i, int j, int k, int add ) {
// Add a new residue position <add> off the end of the chain away from mid <j> and <k>
Vec	m = cas[j] & cas[k];
Vec	v = m - cas[i];
	cas[add] = m+v;
}

void extend ( Prot* chain ) {
// Add dummy CA atoms to the end of the chain
int	n = chain->len, m = n+1; 
        addres(chain->cas,  3,  2,1, 0);	// add res 0 to Nterm
        addres(chain->cas,n-2,n-1,n, m);	// add res N+1 to Cterm
	chain->res[0] = 'n';
	chain->res[m] = 'c';
}

void add_cb ( Prot *chain, float bond ) {
// Add dummy CB atoms to the chain (but not to 0 or N+1)
        FIR(i,chain->len) { Vec n,c,b; float d;
		n = (chain->cas[i] - chain->cas[i-1]).getNorm();
		c = (chain->cas[i] - chain->cas[i+1]).getNorm();
		b = (n+c).getNorm()*bond;
		chain->cbs[i] = chain->cas[i]+b;
        }
}

void getpdb ( Prot* chain, FILE *pdb ) { 
// Read a PDB format file from stream <pdb> into structure Prot 1..N
int	io, n = 1;
char	line[222], junk[30];
Vec	xyz[2222];
float	acc[2222];
        while(1) { float x,y,z, s,t; int a,b; char c;
		io = read_line(pdb,line);
		if ( io < 1 ) break;
		if (!strncmp(line,"END",3)) break;
		if (!strncmp(line,"TER",3)) continue; // LINKs may follow
		if (strncmp(line,"ATOM   ",7) && strncmp(line,"HETATM ",7)) continue; // MSE is a HETATM
		if (strncmp(line+13,"P  ",3)) continue;
		sscanf(line,"%30c %f%f%f%f%f", junk, &x, &y, &z, &s, &t);
		xyz[n].x = x; xyz[n].y = y; xyz[n].z = z; acc[n] = t;
		strncpy(aa3[n],junk+17,3);
		n++;
	}
	chain->len = n-1;
	n++; // 2 extra positions for dummy chain extensions
	chain->cas = new Vec[n];
	chain->cbs = new Vec[n];
	chain->sec = new int[n];
	chain->rid = new int[n];
	chain->res = new char[n];
	chain->acc = new float[n];
	FIR(i,chain->len) {
		chain->cas[i].x = xyz[i].x;
		chain->cas[i].y = xyz[i].y;
		chain->cas[i].z = xyz[i].z;
		chain->acc[i] = acc[i];
	}
	extend(chain);		// add 0 and N+1 positions
	add_cb(chain,2.0);	// add pseudo CB/centroid
}

float** getmat ( Prot *chain )
{
int	len = chain->len, lenn = len+2; // allow for dummy atoms 0 and N+1
float	**mat = new float*[lenn];
	FOR(i,lenn) mat[i] = new float[lenn];
	FOR(i,lenn) FOR(j,lenn) {
		mat[i][j] = 0.0;
		if (i<j) mat[i][j] = chain->cas[i] | chain->cas[j];
		if (i>j) mat[i][j] = chain->cbs[i] | chain->cbs[j];
	}
	return mat;
}

// 1 = pdb file, 2 = constraints file, 3 = cycles

int main (int argc, char** argv) {
// XYZ from simprot realigned with tube axis along Z (1st pdb line = prot.dat top GROUP line)
Vec	rms, rog;
float	sum0, sum1, sum2;
int	n, len, in=500;
char	line[111];
Pairs	cons[1111];
float	**dij;
Prot	*chain = new Prot;
Prot	*known = new Prot;
FILE	*pdb = fopen(argv[1],"r");
FILE	*con = fopen(argv[2],"r");
int	cycles;
	getpdb(chain,pdb);
	fclose(pdb);
	len = chain->len;
	Pi(len) NL
	n = 0;
	FOR(i,in) { int a,b,c,d; float s;
		if (read_line(con,line) <= 0) break;
		sscanf(line,"%d %d %d %d %f", &a, &b, &c, &d, &s);
		if (a>b) { int e=a; a=b; b=e; } // keep a<b
		if (d>9000) d -= 9000; // added to stop early breaks in sim
		cons[n].a = a;
		cons[n].b = b;
		cons[n].s = s;
		n++;
	}
	sscanf(argv[3],"%d", &cycles);
	Pi(cycles) NL
	smooth(chain->cas,len,cycles);
	pdbout(chain,"smooth.pdb");
}

void pdbout ( Prot *a, char *file ) 
{
FILE	*pdb;
	pdb = fopen(file,"w");
	FIR(i,a->len) { float w = 0.0;
		fprintf(pdb,"ATOM%7d  P     G A%4d     %7.3f %7.3f %7.3f  0.00 %5.2f\n", i,i,
			a->cas[i].x, a->cas[i].y, a->cas[i].z, w);
	}
	fprintf(pdb,"TER\n");
	fclose(pdb);
}

void smooth ( Vec *a, int len, int cycles)
{
	FOR(i,cycles)
	{ Vec	p1 = a[0], p2 = a[1], p3 = a[2];
		FIR(j,len) { Vec q;
			q = (p1+p2+p3)/3;
			p1 = a[j]; p2 = a[j+1]; p3 = a[j+2];
			a[j] = q;
		}
	}
}
