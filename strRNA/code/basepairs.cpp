/*
c++ basepairs.cpp -o basepairs sims/util.o sims/geom.o -lm

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
int	n, len, lenn, **bonds, *secsi, *secsj, *secsn;
float	x,y,z;
	in = fopen("pairs.ct","r");
	fscanf(in,"%d", &len);
	Pi(len) NL
	len++;
	hyd = new Atom[len*10];
	hev = new Atom[len*10];
	lenn = len+2;
	secsi = new int[lenn];
	secsj = new int[lenn];
	secsn = new int[lenn];
	bonds = new int*[lenn]; FOR(i,lenn) bonds[i] = new int[lenn];
	FOR(i,lenn) {
		secsi[i] = secsj[i] = secsn[i] = 0;
		FOR(j,len) bonds[i][j] = 0;
	}
	// read basepairs
	n = 0;
	next_line(in);
	DO { int io = read_line(in,line); char b; int a,c,d,e,f;
		if (io <= 10) break;
		sscanf(line,"%d %c %d %d %d %d", &a,&b,&c,&d,&e,&f);
		if (e==0) continue;
		bonds[a][e] = bonds[e][a] = 1;
		n++;
	}
	Pi(n) NL
	// count bonds/res
	FOR(i,len) {
		FOR(j,len) {
			if (bonds[i][j] == 0 ) continue;
			secsi[i] += j; secsn[i]++;
			secsi[j] += i; secsn[j]++;
		}
	}
	FOR(i,len) {
		if (secsn[i]>1) secsi[i] /= secsn[i];
	}
	FIR(i,len) {
		if (secsi[i-1]==0 || secsi[i]==0) continue;
		if (abs(secsi[i-1]-secsi[i])<10) continue;
		Pt(break) Pi(secsi[i-1]) Pi(secsi[i]) NL
		secsi[i-1] = secsi[i] = 0;
		secsn[i-1] = secsn[i] = 0;
	}
	FOR(j,100) {
		//for (int i=50; i<len; i++) printf("%3d", secsi[i]); NL
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
