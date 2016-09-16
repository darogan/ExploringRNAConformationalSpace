/*
c++ hbonds.cpp -o hbonds sims/util.o sims/geom.o -lm

./hbonds 94 100
         |   |
        len cycles

*/
#include "sims/util.hpp"
#include "sims/geom.hpp"

typedef struct { Vec xyz; int id; } Atom;

int median ( int, int*, int* );
int avrage ( int, int*, int* );
int smooth ( int, int*, int* );

int main (int argc, char** argv) {
Atom	*hyd, *hev;
char	line[111];
FILE	*in,*out;
int	n, cyc, len, hydn, hevn, **bonds, *secsi, *secsj, *secsn;
float	x,y,z;
	sscanf(argv[1],"%d", &len);
	sscanf(argv[2],"%d", &cyc);
	Pi(len) Pi(cyc) NL
	len++;
	hyd = new Atom[len*10];
	hev = new Atom[len*10];
	secsi = new int[len];
	secsj = new int[len];
	secsn = new int[len];
	bonds = new int*[len]; FOR(i,len) bonds[i] = new int[len];
	FOR(i,len) {
		secsi[i] = secsj[i] = secsn[i] = 0;
		FOR(j,len) bonds[i][j] = 0;
	}
	in = fopen("hydro.dat","r");
	hydn = 0;
	DO { int io = read_line(in,line); 
		if (io <= 0) break;
		sscanf(line,"%d %f%f%f", &n, &x,&y,&z);
		hyd[hydn].xyz.x = x; hyd[hydn].xyz.y = y; hyd[hydn].xyz.z = z;
		hyd[hydn].id = n;
		hydn++;
	}
	Pt(Hydrogens) Pi(hydn) NL
	fclose(in);
	in = fopen("heavy.dat","r");
	hevn = 0;
	DO { int io = read_line(in,line); 
		if (io <= 0) break;
		sscanf(line,"%d %f%f%f", &n, &x,&y,&z);
		hev[hevn].xyz.x = x; hev[hevn].xyz.y = y; hev[hevn].xyz.z = z;
		hev[hevn].id = n;
		hevn++;
	}
	Pt(Heavy atm) Pi(hevn) NL
	fclose(in);
	// count Hbonds
	n = 0;
	FOR(i,hevn) { float a,b,c;
		FOR(j,hydn) {
			if (hev[i].id != hyd[j].id) continue;
			a = hev[i].xyz|hyd[j].xyz;		// a = A--H
			if (a > 1.1) continue;
			FOR(k,hevn) { int id = hev[i].id, kd = hev[k].id;
				if (abs(id-kd) < 4 ) continue;
				b = hev[k].xyz|hyd[j].xyz;	// b = H...B
				if (b > 3.5) continue;
				c = hev[i].xyz|hev[k].xyz;	// c = A--H...B
				if (c < 0.8*(a+b)) continue;
				bonds[id][kd]++;
				bonds[kd][id]++;
				n++;
			}
		}
	}
	Pi(n) NL
	// count bonds/res
	FOR(i,len) {
		FOR(j,len) {
			if (bonds[i][j] <  2) {
				bonds[i][j] = 0;
				continue;
			}
			secsi[i] += j; secsn[i]++;
			secsi[j] += i; secsn[j]++;
		}
	}
	FOR(i,len) {
		if (secsn[i]>1) secsi[i] /= secsn[i];
	}
	/* smooth?
	FOR(j,1) {
		FIR(i,len) secsj[i] = smooth(i,secsi,secsn);
		FIR(i,len) secsi[i] = smooth(i,secsj,secsn);
	}
	*/
	FOR(j,cyc) {
		FIR(i,len) secsj[i] = median(i,secsi,secsn);
		FIR(i,len) secsi[i] = avrage(i,secsj,secsn);
	}
	out = fopen("secs.out","w");
	FOR(i,len) {
		if (secsi[i]==0) continue;
		fprintf(out,"%d %d\n", i, secsi[i]);
	}
	fclose(out);
	out = fopen("secs.dat","w");
	FOR(i,len) {
		fprintf(out,"%d %d\n", i, secsi[i]);
	}
	fclose(out);
	out = fopen("pairs.out","w");
	FOR(i,len) FOR(j,len) {
		if (i > j) continue;
		if (bonds[i][j]==0) continue;
		fprintf(out,"%d %d %d\n", i,j,bonds[i][j]);
	}
	fclose(out);
}

int median ( int i, int *s, int *w ) {
// return weighted median
int	m, n, mid, med, num[99];
	if (w[i-1]==0 && w[i]==0) return 0;
	if (w[i+1]==0 && w[i]==0) return 0;
	if (w[i-1]==0 && w[i+1]==0) return s[i];
	n = 0;
	FOR(j,w[i-1]) num[n++] = s[i-1];
	FOR(j,w[i]  ) num[n++] = s[i];
	FOR(j,w[i+1]) num[n++] = s[i+1];
	sort(num,n);
	mid = n/2;
	med = num[mid];
	if (w[i]>1) w[i]--;
	if (w[i]<1) w[i] = 1;
	return med;
}

int avrage ( int i, int *s, int *w ) {
// return weighted average
int	a, m, n;
	if (w[i-1]==0 && w[i]==0) return 0;
	if (w[i+1]==0 && w[i]==0) return 0;
	if (w[i-1]==0 && w[i+1]==0) return s[i];
	a = s[i-1]*w[i-1] + s[i]*w[i] + s[i+1]*w[i+1];
	n = w[i-1] + w[i] + w[i+1];
	m = a/n;
	return m;
}

int smooth ( int i, int *s, int *w ) {
// return weighted smoothed value
int	a, m, n;
	n = w[i-1] + w[i] + w[i+1];
	if (n==0) return 0;
	a = s[i-1]*w[i-1] + s[i]*w[i] + s[i+1]*w[i+1];
	m = a/n;
	w[i]++;
//Pi(i) Pi(s[i-1]) Pi(s[i]) Pi(s[i+1]) Pi(m) NL
	if (s[i]<m) return s[i]+1;
	if (s[i]>m) return s[i]-1;
}
