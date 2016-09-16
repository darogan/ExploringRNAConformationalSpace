/*
c++ scoreWrms.cpp -o scoreWrms sims/util.o sims/geom.o -lm

./scoreWrms model.cas psicov.top (needs true.pdb)

*/
#include "sims/util.hpp"
#include "sims/geom.hpp"

#include "util/aa/incl/matrix.h"
#include "util/aa/incl/bestrot.h"
#include "util/aa/incl/siva.h"
#include "util/aa/incl/ql.h"

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

Vec	rms121 ( Prot*, Prot* );
void	smooth ( Vec*, int, int );
double	supermac ( double*, double**, double**, int, double*, double*, double** );
void pdbout ( Prot*, Prot*, double**, double**, int, double*, double*, double**, char* ); 

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

Vec rogs ( Prot *A, float *w ) {
// get the weighted RoG around a rough bundle axis
int	m, n, in, len = A->len;
Vec	*c = new Vec[len];
Pairs	*p = new Pairs[999];	
int	q[999];
Seg	axis;
Vec	cog, move;
Vec	sum, rog;
	m = n = 0;
	cog.zero();
	FIR(i,len) { // caps have weight = 1
		if ( w[i]>0.5) { cog += A->cas[i]; m++; }
		if ( w[i]>0.5 && w[i]<1.5) c[n++] = A->cas[i];
	}
	cog /= (float)m;
	m = 0;
	FOR(i,n-8) { // for pairs of cap residues...
		for (int j=i-1; j<i+8; j++) // between adjacent caps
		{ float	d = c[i] | c[j];
			if (d < 30.0) continue;
			p[m].a = i; p[m].b = j; p[m].s = d; m++;
			if (m==999) break;
		}
	}
	sort(p,q,m); // reverse sort with rank in <q[]>
	axis.A = c[p[q[0]].a]; axis.B = c[p[q[0]].b]; // first axis guess is widest cap pair
	in = 1;
	FOR(j,m)
	{ int	i = q[j], pa = p[i].a, pb = p[i].b;
	  Vec	a = c[pa], b = c[pb];
	  Seg	now = (axis.A/(float)j, axis.B/(float)j);
	  float dA = a|now.A, dB = a|now.B;
		if (dA < dB) {
			axis.A += a; axis.B += b; 
		} else {
			axis.A += b; axis.B += a;
		} 
	}
	axis.A /= (float)m;
	axis.B /= (float)m;
	move = cog-(axis.A & axis.B);
	axis.A += move; axis.B += move;
	// the bundle axis is now centred on the CoG of the TM segments
	rog.zero(); sum.zero();
	FIR(i,len)
	{ float d = axis.vec_to_line(A->cas[i]),
		dd = d*d, wi = w[i];
		rog.x += dd;    sum.x += 1.0;
		rog.y += dd*wi; sum.y += wi;
		wi += 0.5;
		rog.z += dd*wi; sum.z += wi;
	}
	rog.x /= sum.x; rog.x = sqrt(rog.x);
	rog.y /= sum.y; rog.y = sqrt(rog.y);
	rog.z /= sum.z; rog.z = sqrt(rog.z);
	return rog;
}

// 1 = pdb file, 2 = constraints file, 3 = known str., 4 = weights

int main (int argc, char** argv) {
// XYZ from simprot realigned with tube axis along Z (1st pdb line = prot.dat top GROUP line)
Vec	rms, rog;
float	bad, bad0, bad1, bad2;
int	n, len, in=500;
char	line[111];
Pairs	cons[1111];
float	**dij;
Prot	*chain = new Prot;
Prot	*known = new Prot;
FILE	*con = fopen(argv[1],"r");
FILE	*pdb = fopen(argv[2],"r");
float	want, give, bump = 4.0;
int	nbump = 0;
	getpdb(chain,pdb);
	fclose(pdb);
	len = chain->len;
	dij = getmat(chain); // i<j=CA i>j=CB
	n = 0;
	FOR(i,in) { int a,b,c,d; float s;
		if (read_line(con,line) <= 0) break;
		sscanf(line,"%d %d %d %d", &a, &b, &c, &d);
		if (a>b) { int e=a; a=b; b=e; } // keep a<b
		cons[n].a = a;
		cons[n].b = b;
		cons[n].s = d;
		n++;
	}
	Pi(n) NL
	bad = bad0 = bad1 = bad2 = 0;
	FOR(i,n) { int a = cons[i].a, b = cons[i].b; float s = cons[i].s, d, p;
		d = dij[a][b];  // CB=[b][a], CA=[a][b] = P
		if (d < 21.0) continue;
		Pi(i) Pi(a) Pi(b) Pr(d) NL
		bad++;
		if (bad==1) bad0 = i;	// 1st bad
		p = 100.0*(float)bad/(float)len;
		if (bad1==0 && p>2.0) bad1 = i;
		if (bad2==0 && p>5.0) bad2 = i;
	}
	Pi(bad0) Pi(bad1) Pi(bad2) NL
}
