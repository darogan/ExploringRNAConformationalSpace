/*
c++ -ggdb stems.cpp -o stems sims/util.o sims/geom.o -lm


*/
#include "sims/util.hpp"
#include "sims/geom.hpp"

typedef struct {
	Pairs	ends[6]; // start,end for tail-Watson-loop-loop-Crick-tail
	int	n;
} Group;

typedef struct {
	char	*res;	// residue type
	float	*acc;	// surface accesssibility
	int	*sec;	// secondary structure type (0=c, 1=A, 2=B, 3=3ten)
	int	*rid;	// residue number
	Vec	*cas;	// P position (runs 1...N with dummy 0 and N+1)
	Vec	*cbs;	// dummy CB/centroid (extended <bond>A from CA)
	int	len;	// chain length
	int	gaps;	// number of chain breaks
	float	pcta, pctb;	// percent alpha and beta structure
} Prot;

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
// Add dummy atom to the chain (but not to 0 or N+1)
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

// 1 = pdb file, 2 = constraints file

int main (int argc, char** argv) {
// XYZ from simprot realigned with tube axis along Z (1st pdb line = prot.dat top GROUP line)
Vec	rms, rog;
float	sum0, sum1, sum2;
int	n, len, in=500;
char	line[111];
Pairs	cons[1111];
float	**dij;
Prot	*chain = new Prot;
FILE	*con = fopen(argv[1],"r");
FILE	*pdb = fopen(argv[2],"r");
FILE	*seg = fopen("line.plot","r");
FILE	*dat;
float	want, give, bump = 4.0;
int	ncon, nseg, nbump = 0;
int	segW[22][2], segC[22][2], five[22], rank[22];
int	*done, *list;
Group	*segs;
	getpdb(chain,pdb);
	fclose(pdb);
	dij = getmat(chain); // i<j=CA i>j=CB
	len = chain->len;
	FOR(i,len) FOR(j,len) { float d;
		if (i<j) continue;
		if (i-j < 5) continue;
		d = chain->cas[i] | chain->cas[j];
		if (d < bump) nbump++;
	}
	Pi(nbump) NL
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
	ncon = n;
	nseg = 1;
	n = 0;
	DO { int io, i, j;
		io = read_line(seg,line);
		if (io < 0) break;
		if (io < 9) {	// next seg
			n = 0;
			continue;
		} else {	// get ends
			sscanf(line,"%d %d", &i, &j);
			segW[nseg][n] = i; segC[nseg][n] = j;
			n++;
		}
		if (n>1) nseg++;
	}
	Pi(nseg) NL
	five[0] = 0;
	FIR(i,nseg-1)
	{ int	w0 = segW[i][0], w1 = segW[i][1],
		c0 = segC[i][0], c1 = segC[i][1];
		// make all 5'--> 3'
		if (w1 < w0) { n = w0; w0 = w1; w1 = n; }
		if (c1 < c0) { n = c0; c0 = c1; c1 = n; }
		// lowest start first
		if (w0 < c0) {
			segW[i][0] = w0; segW[i][1] = w1;
			segC[i][0] = c0; segC[i][1] = c1;
		} else { // swap
			segW[i][0] = c0; segW[i][1] = c1;
			segC[i][0] = w0; segC[i][1] = w1;
		}
		five[i] = segW[i][0];
	}
	five[nseg] = len;
	sort(five,rank,-(nseg+1));
	// add dummy terminal markers
	segC[0][0] = segW[0][0] = segC[0][1] = segW[0][1] = 0;
	segC[nseg][0] = segW[nseg][0] = segC[nseg][1] = segW[nseg][1] = len;
	FIR(j,nseg-1)
	{ int	i = rank[j],
		w0 = segW[i][0], w1 = segW[i][1],
		c0 = segC[i][0], c1 = segC[i][1];
		Pi(w0) Pi(w1) Pt(---) Pi(c0) Pi(c1) NL
	} NL
	list = new int[len+1];
	done = new int[len+1];
	FOR(i,len) list[i] = done[i] = 0;
	segs = new Group[nseg];
	n = 0;
	FIR(j,nseg-1)
	{ int	i = rank[j], last, next, in, m=j-1,
		w0 = segW[i][0], w1 = segW[i][1],
		c0 = segC[i][0], c1 = segC[i][1];
		FOR(k,6) segs[m].ends[k].a = segs[m].ends[k].b = -1;
		segs[m].n = 0;
		if (w0>1) {  int gap, mingap = 999;	// something before w0
			FOR(k,nseg+2) {	// look for closest end to w0
				if (k==i) continue;
				//      k seg finish - start of next strand i
				gap = w0 - segW[k][1];
				if (gap>0 && gap<mingap) { last = segW[k][1]; mingap = gap; }
				gap = w0 - segC[k][1];
				if (gap>0 && gap<mingap) { last = segC[k][1]; mingap = gap; }
			} 
			if (mingap > 1) {
				last += mingap/2;
				printf("%d end1 = %3d -%3d\n",j,last,w0-1);
				for (int k=last; k<w0; k++) {
					if (done[k]) continue;
					if (segs[m].ends[0].a < 0) segs[m].ends[0].a = k;
					segs[m].ends[0].b = k;
					list[n++] = k;
					done[k] = 1;
				}
				segs[m].n++;
				Pt(:) NL
			}
		}
		printf("%d segW = %3d -%3d\n",j,w0,w1);
		segs[m].ends[1].a = w0;
		segs[m].ends[1].b = w1;
		segs[m].n++;
		for (int k=w0; k<=w1; k++) {
			chain->sec[k] = 2;
			list[n++] = k;
			done[k] = 1;
		}
		in = c0-w1;
		if (in>1) {  int gap, mingap = in;	// something between
			next = c0;
			FOR(k,nseg+2) {	// look for closest end to w1
				if (k==i) continue;
				//      k seg starts - end of last strand i
				gap = segW[k][0] - w1;
				if (gap>0 && gap<mingap) { next = segW[k][0]; mingap = gap; }
				gap = segC[k][0] - w1;
				if (gap>0 && gap<mingap) { next = segC[k][0]; mingap = gap; }
			} 
			next -= mingap/2;
			printf("%d gap1 = %3d -%3d\n",j,w1+1,next);
			for (int k=w1+1; k<=next; k++) {
				if (done[k]) continue;
				if (segs[m].ends[2].a < 0) segs[m].ends[2].a = k;
				segs[m].ends[2].b = k;
				list[n++] = k;
				done[k] = 1;
			}
			segs[m].n++;
		}
		Pt(:) NL
		if (in>2) {  int gap, mingap = in;	// write the other half
			last = w1;
			FOR(k,nseg+2) {	// look for closest end to c0
				if (k==i) continue;
				//      k seg finish - start of next strand i
				gap = c0 - segW[k][1];
				if (gap>0 && gap<mingap) { last = segW[k][1]; mingap = gap; }
				gap = c0 - segC[k][1];
				if (gap>0 && gap<mingap) { last = segC[k][1]; mingap = gap; }
			} 
			last += mingap/2;
			if (last==next) last++;
			printf("%d gap2 = %3d -%3d\n",j,last,c0-1);
			for (int k=last; k<c0; k++) {
				if (done[k]) continue;
				if (segs[m].ends[3].a < 0) segs[m].ends[3].a = k;
				segs[m].ends[3].b = k;
				list[n++] = k;
				done[k] = 1;
			}
			if (segs[m].ends[3].b+1 != segs[m].ends[3].a) segs[m].n++;
		}
		printf("%d segC = %3d -%3d\n",j,c0,c1);
		segs[m].ends[4].a = c0;
		segs[m].ends[4].b = c1;
		segs[m].n++;
		for (int k=c0; k<=c1; k++) {
			list[n++] = k;
			done[k] = 1;
		}
		if (c1 < len) {  int gap, mingap = 999;	// something following
			FOR(k,nseg+2) {	// look for closest end after c1
				if (k==i) continue;
				//      k seg starts - end of last strand i
				gap = segW[k][0] - c1;
				if (gap>0 && gap<mingap) { next = segW[k][0]; mingap = gap; }
				gap = segC[k][0] - c1;
				if (gap>0 && gap<mingap) { next = segC[k][0]; mingap = gap; }
			} 
			if (mingap > 2) {
				Pt(:) NL
				next -= mingap/2;
				printf("%d end2 = %3d -%3d\n",j,c1+1,next);
				for (int k=c1+1; k<=next; k++) {
					if (done[k]) continue;
					if (segs[m].ends[5].a < 0) segs[m].ends[5].a = k;
					segs[m].ends[5].b = k;
					list[n++] = k;
					done[k] = 1;
				}
			}
		}
		Pt(----------------) NL
	}
	if ( n != len ) { Pt(Length mismatch) Pi(n) Pi(len) NL }
	FIR(i,len) {
		if (done[i]) continue;
		Pt(Missing) Pi(i) NL
	}
	pdb = fopen("stems.pdb","w");
	dat = fopen("stems.dat","w");
	nseg--;
	fprintf(dat,"GROUP 0 %d\n", nseg);
	FOR(i,nseg) { int a,b; Vec end5,end3;
		a = segs[i].ends[1].a; b = segs[i].ends[4].b;
		end5 = chain->cas[a] & chain->cas[b];
		a = segs[i].ends[1].b; b = segs[i].ends[4].a;
		end3 = chain->cas[a] & chain->cas[b];
		Pi(i) Pi(segs[i].n) NL
		fprintf(dat,"  GROUP -1 %d", segs[i].n-1);  // double strand counts as 1
		fprintf(dat,"   %5.1f %5.1f %5.1f    %5.1f %5.1f %5.1f\n", end5.x,end5.y,end5.z, end3.x,end3.y,end3.z);
		FOR(j,6) {
			if (segs[i].ends[j].a < 0) continue;
			Pi(j) Pi(segs[i].ends[j].a) Pi(segs[i].ends[j].b) NL
			n = segs[i].ends[j].b - segs[i].ends[j].a +1;
			if (j!=1 && j!=4) { 
				fprintf(dat,"    GROUP 0 %d\n", n);
			} else {
				if (j==1) fprintf(dat,"    DOUBL %d %d\n", i+1,n);
				if (j==4) fprintf(dat,"    DOUBL %d %d\n", -(i+1),n);
			}
			for (int k = segs[i].ends[j].a; k <= segs[i].ends[j].b; k++)
			{ float sec = (float)(j+1);
				fprintf(dat,"      ATOM%7d  P     G A%4d     %7.3f %7.3f %7.3f  1.00  0.00\n", k,k,
					chain->cas[k].x, chain->cas[k].y, chain->cas[k].z);
				fprintf(pdb,"ATOM%7d  P     G A%4d     %7.3f %7.3f %7.3f  1.00 %5.2f\n", k,k,
					chain->cas[k].x, chain->cas[k].y, chain->cas[k].z, sec);
			}
		}
		fprintf(pdb,"TER\n");
	}
	fclose(dat);
	fclose(pdb);
}
