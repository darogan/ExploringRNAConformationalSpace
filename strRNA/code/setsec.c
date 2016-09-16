/*
cc setsec.c -o setsec util/wt/util.o -lm

*/
#include "util/wt/incl/util.h"
#include "util/wt/incl/geom.h"

float parse ( float, int );
float score ( FILE*, int, int, int, int, int );
float stripe ( float**, int, int, float[99][99], int, int );
float moment ( float**, int, int, int, int );
float momat ( int, int, float**, int, int );

int	inbox = 0;
int	len, segs;
float	**con, **mat, **sum, *emax;
int	**top[2], *edge[2], *ends[2], *base[2];
int	sec[999], set[999];

main ( int argc, char** argv )
{
float	cut = 0.0, gap = 1.0;
int	i, j, n, tm, last, lenn, cycles = 1;//3;
FILE	*pred, *grid, *diag, *line, *boxA, *boxP, *pack, *term;
FILE	*pin = fopen("pairs.dat","r");
FILE	*sin = fopen("secs.dat","r");
	for (i=0; i<999; i++) sec[i] = 0;
	// read TM/sec segment consensus
	n = 0;
	len = 0;
	last = 0;
	segs = 0;
	pred = fopen("pred.dat","w");
	// count TM predictions (M) for each res in topcons (m = added helix)
	while (1) { int io; char secin[22];
		io = read_line(sin,secin);
		if (io <=0) break;
		sscanf(secin,"%d %d", &i, &j);
		sec[len] = j;
		if (len) last = sec[len-1];
		if (sec[len]) { // in a segment
			if (last==0) { // in a new segment
				n++; last = sec[n];
			}			
		}
		set[len] = n;
		if (sec[len]) fprintf(pred,"%d %d %d\n", len+1,len+1,sec[len]);
		len++;
	}
	segs = n;
	for (i=0; i<len; i++) if (sec[i]==0) set[i] = 0;
	if (argc>1) sscanf(argv[1],"%f", &cut);
	if (argc>2) sscanf(argv[2],"%f", &gap);
	if (cut < 0.0) { inbox = 1; cut = -cut; cycles = 1; }
	Pi(len) Pr(cut) Pr(gap) Pi(inbox) NL
	// allocate arrays
	lenn = len+1;
	con = (float**)malloc(sizeof(float*)*lenn); TEST(con)
	mat = (float**)malloc(sizeof(float*)*lenn); TEST(mat)
	sum = (float**)malloc(sizeof(float*)*lenn); TEST(sum)
	top[0] = (int**)malloc(sizeof(int*)*lenn);  TEST(top[0])
	top[1] = (int**)malloc(sizeof(int*)*lenn);  TEST(top[1])
	for (i=0; i<lenn; i++) {
		con[i] = (float*)malloc(sizeof(float)*lenn); TEST(con[i])
		mat[i] = (float*)malloc(sizeof(float)*lenn); TEST(mat[i])
		sum[i] = (float*)malloc(sizeof(float)*lenn); TEST(sum[i])
		top[0][i] = (int*)malloc(sizeof(int)*lenn);  TEST(top[0][i])
		top[1][i] = (int*)malloc(sizeof(int)*lenn);  TEST(top[1][i])
		for (j=0; j<lenn; j++) con[i][j] = 0.0;
	}
	emax = (float*)malloc(sizeof(float)*lenn);
	edge[0] = (int*)malloc(sizeof(int)*lenn);
	edge[1] = (int*)malloc(sizeof(int)*lenn);
	ends[0] = (int*)malloc(sizeof(int)*lenn);
	ends[1] = (int*)malloc(sizeof(int)*lenn);
	base[0] = (int*)malloc(sizeof(int)*lenn);
	base[1] = (int*)malloc(sizeof(int)*lenn);
        last = 0;
        n = -1;
        for (i=1; i<len; i++) {
                if (set[i] && last==0) {
                        n++;
                        ends[0][n] = i;
                }
                if (set[i]==0 && last) {
                        ends[1][n] = i;
                }
                last = set[i];
        }
	if (last) ends[1][n] = len;
        n++;
        //for (i=0; i<len; i++) printf("%2d",set[i]); NL
        for (i=0; i<n; i++) { Pi(i) Pi(ends[0][i]) Pi(ends[1][i]) NL }
	// read predicted contacts
	n = 0;
	while (1) { int io, a,b,c,d; float x,e; char conin[22];
		io = read_line(pin,conin);
		if (io <=0) break;
		sscanf(conin,"%d %d %d %d %f", &a,&b,&c,&d,&e);
		if (a>len || b>len) continue;
		a--; b--; // pairs read in range 1..N
		// NB values accumulate and 1/2 weight for reverse indices
		// x = (float)n;
		// e = exp(-x*x*0.01/(float)len); // so all souces score same
		//if (a>b) e *= 0.5;
		con[a][b] += e;
		con[b][a] += e; 
		n++;
	}
	Pi(n) NL
/* print contact matrix
for(i=0; i<len; i++) { for (j=0; j<len; j++) {
	if (i==j) printf("%3d",set[i]); else printf("%3d",(int)(100.0*con[i][j]));
}  NL } NL
*/
	for (i=0; i<cycles; i++) { float s;
		s = parse(gap,i);
		Pi(i) Pr(s) NL
	}
	NLL
	for (i=0; i<lenn; i++) for (j=0; j<lenn; j++) mat[i][j] = sum[i][j] = 0.0;
	term = fopen("ends.dat","w");
	pack = fopen("pack.dat","w");
	grid = fopen("grid.plot","w");
	boxA = fopen("boxA.plot","w");
	boxP = fopen("boxP.plot","w");
	for (i=0; i<segs; i++)
	{	int m, w, ai,bi, aj,bj, zero=0;
		ai = ends[0][i]; bi = ends[1][i];
		w = (bi-ai)/2;
		m = ai+w-1;
		fprintf(term,"%5d %5d\n", ai,bi);
		printf("\t%d-%d\t",ai,bi); Pi(m) Pi(w) NL
		fprintf(grid,"%4d %4d\n",  zero,ai);
		fprintf(grid,"%4d %4d\n\n", len,ai);
		fprintf(grid,"%4d %4d\n",  zero,bi);
		fprintf(grid,"%4d %4d\n\n", len,bi);
		fprintf(grid,"%4d %4d\n",  ai,zero);
		fprintf(grid,"%4d %4d\n\n", ai,len);
		fprintf(grid,"%4d %4d\n",  bi,zero);
		fprintf(grid,"%4d %4d\n\n", bi,len);
		for (j=0; j<segs; j++) // score each box
		{ float s; int a,b; FILE *box;
			if (i <= j) continue;
			aj = ends[0][j]; bj = ends[1][j];
			s = score(pack,m,w,i,j,-1);
			if (s < 0.0) { s = -s; box = boxP; } else { box = boxA; }
			if (s < cut) continue;
			fprintf(box,"%4d %4d %f\n",  ai,aj,s);
			fprintf(box,"%4d %4d\n",  ai,bj);
			fprintf(box,"%4d %4d\n",  bi,bj);
			fprintf(box,"%4d %4d\n",  bi,aj);
			fprintf(box,"%4d %4d\n\n",ai,aj);
			fprintf(box,"%4d %4d\n",  aj,ai);
			fprintf(box,"%4d %4d\n",  aj,bi);
			fprintf(box,"%4d %4d\n",  bj,bi);
			fprintf(box,"%4d %4d\n",  bj,ai);
			fprintf(box,"%4d %4d\n\n",aj,ai);
			for (a=ai; a<=bi; a++) for (b=aj; b<=bj; b++) sum[a][b] = mat[a][b];
		}
	}
	diag = fopen("diag.plot","w");
	line = fopen("line.dat","w");
	for (i=0; i<lenn; i++) { for (j=0; j<lenn; j++) {
		if (sum[i][j] < 0.0) { int sup; float s = -sum[i][j];
			// separate score and mark 
			sup = (int)(s/1000.0);
			s -= 1000.0*(float)sup;
			fprintf(line,"%5d %5d %7.1f\n", i,j,s);
			sum[i][j] = (s-sum[i][j])/1000.0-0.2;
		}
		if (sum[i][j] > 1.1) {
			if (i>j) fprintf(diag,"%5d %5d\n", i,j);
			if (j>i) fprintf(diag,"%5d %5d\n", i,j);
		}
	} }
}

float parse ( float gap, int cycle )
// cycle = 0...N (cycle = -1 on final rescoring pass)
{
int	maxseg = 50;
int	minseg = 1; // ie 1+1+1 = 3
int	mingap = 4;
int	i, j, n;
int	m, w, topm, topw;
float	tops;
	for (m=0; m<len; m++) for (w=0; w<len; w++) mat[m][w] = sum[m][w] = 0.0;
	// fill mat with the raw scores
	for (m=2; m<len-2; m++) { int skip; float d, e; // m = 5..len-5 allows w=5 at each end (min seg = 11)
		for (w=2; w<maxseg/2; w++) { float s;
			if (m-w < 0) break;
			if (m+w > len-1) break;
			//mat[m][w] = -gap*(float)(w+1);
			if (w>mingap) mat[m][w] = -gap*sqrt((float)(w-mingap));
			mat[m][w] += score(0,m,w,0,0,cycle);
		}
	}
	for (m=0; m<len; m++) {
		if (set[m]) printf("%3d | %2d %2d", m+1,sec[m],set[m]); else printf("%3d |      ", m+1);
		sum[m][0] = mat[m][0] = 0.0;
		sum[m][1] = mat[m][1] = 0.0;
		for (w=1; w<20; w++) printf("%3d", (int)(mat[m][w])); NL
	}
	// fill sum with the accummulated scores
	for (w=2; w<maxseg/2; w++) {
		//WAS  for (m=2; m<len-2; m++) {
		for (m=2; m<len-2; m++) {
			if (m-w < 0) continue;
			if (m+w > len-1) continue;
			if (w==2 && sum[m][1] > 0.0) {
					sum[m][w] = mat[m][w]+sum[m][1]+sum[m-1][1]+sum[m+1][1];
			} else {
					sum[m][w] = sum[m-1][w-1]+sum[m+1][w-1]-sum[m][w-2]+mat[m][w]+mat[m][w-1];
			}
		}
	}
	NLL
	for (m=0; m<len; m++) {
		if (set[m]) printf("%3d | %2d %2d", m+1,sec[m],set[m]); else printf("%3d |      ", m+1);
		for (w=1; w<20; w++) {
			if (inbox) { // wipe values outside box (for printing)
				if (set[m-w] != set[m+w]) sum[m][w] = -1.0;
				if (set[m-w]==0 || set[m+w]==0) sum[m][w] = -1.0;
			}
			if (sum[m][w] > 0.0) printf("%4d", (int)sum[m][w]);
					else printf("    ");
		} NL
	}
	for (m=0; m<=len; m++) {
		emax[m] = 0.0;
		edge[0][m] = edge[1][m] = 0;
		for (w=0; w<len; w++) {
			top[0][m][w] = top[1][m][w] = -1;
			if (inbox==0) continue;
			// wipe all scores not inside a segment (set[]>0)
			if (set[m-w] != set[m+w]) sum[m][w] = -1.0;
			if (set[m-w]==0 || set[m+w]==0) sum[m][w] = -1.0;
		}
	}
	// seq runs 0..N-1
	tops = 0.0;
	for (n=3; n<len; n++) {
		m = n+1;
		for (w=1; w<maxseg/2; w++)
		{ float max;
		  int	p, q, r;
			m--;	// move back centre as window gets bigger
			if (m-w < 0) break;
			if (m+w > len) break;
			max = 0.0;
			// r = i = m-w-3; /* -3 forces 3 gaps  between segments */
			// r = i = m-w-2; /* -2 forces 2 gap between segments */
			r = i = m-w-1; /* -1 forces 1 gap  between segments */
			if (r < 0) continue;
			p = q = j = 0;
			while (i>0 && j<len/2) {
				/* scan leading edge for max */
				if (sum[i][j] > max) {
					if (j>minseg) {
						max = sum[i][j];
						p = i; q = j;
					}
				}
				i--; j++;
			}
			// for each diagonal (r) hold max and its position
			emax[r] = max;
			edge[0][r] = p;
			edge[1][r] = q;
			max += sum[m][w];	// score = current (m,w) + last diag max
			for (i=1; i<=r; i++)	// then scan earlier edges for better max
			{ float s = sum[m][w]+emax[i];
				if (s > max) {
					max = s;
					p = edge[0][i];
					q = edge[1][i];
				}
			}
			top[0][m][w] = p;
			top[1][m][w] = q;
                        sum[m][w] = max;	// overwrite sum to hold max
                        if (sum[m][w] > tops) {
                                tops = sum[m][w];
                                topm = m; topw = w;
                        }
		}
	}
	n = 0;
	NLL // Pr(tops) NLL
	i = topm; j = topw;
	while (i>0 && j>0) { int ii, jj;
		Pi(i) Pi(j) Pr(sum[i][j]) NL
		edge[0][n] = i;
		edge[1][n] = j;
		ii = top[0][i][j]; jj = top[1][i][j];
		i = ii; j = jj;
		n++;
	}
	NL
	if (cycle==0) for (i=0; i<len; i++) set[i] = 0; // reset <set> to starting segments
	segs = 0;
	for (i=n-1; i>-1; i--)
	{ int	nn = edge[0][i]-edge[1][i]+1,
		cc = edge[0][i]+edge[1][i]+1;
		j = n-1-i;
		segs++;
		ends[0][j] = nn;
		ends[1][j] = cc;
		if (cycle>0) continue;
		base[0][j] = nn;
		base[1][j] = cc;
		for (j=nn; j<=cc; j++) set[j-1] = segs;
	}
	for (i=0; i<segs; i++) {
		printf("   %d-%d",ends[0][i],ends[1][i]);
	}
	Pr(tops) NL
	return tops;
}

float score ( FILE *out, int m, int w, int si, int sj, int cycle ) {
float	dia[99][99], dip[99][99]; // holds best diagonals in stripe()
float	score, tmpred = 0.0;
float	condif, conmax = 0.0, consum = 0.0;
float	boxsum, pack, pars, ants;
float	**box, **xob; 
int	maxa,maxb;
int	i, j, k, n = 0;
int	a,b, seg1,seg2, sec1,sec2, mixed;
int	am = m-w, pm = m+w;
int	print = 0;
int	one = 0;
int	id = 0;
// print = 1;
if (am<40 && pm>40) print = 1;
	box = (float**)alloca(sizeof(float*)*99);
	xob = (float**)alloca(sizeof(float*)*99);
	for (i=am; i<=pm; i++) { // sum TM prediction for the segment
		if (set[i]) {	   // in a segment
			if (id==0) { // hold the first seg ID
				id = set[i];
			} else {
				if (set[i] != id) return -1.0; // two segs in one window
			}
		}
		tmpred += sqrt((float)abs(sec[i]));
		n++;
	}
	tmpred /= (float)n;
//	if (cycle==0) return tmpred*2.0; // first pass is just TM prediction score
	a = w+w+1; // current seg length
//	if (a < 15) return tmpred;
//	if (a > 55) return tmpred;
	if (cycle < 0) {	// -ve = flag for one box run (segIDs = si,sj)
		one = 100*si+sj; 	// unique id for later sort
		seg1 = sj; seg2 = seg1+1;
	} else {	// default = run over all segs (bar self) 
		seg1 = 0; seg2 = segs;
	}
	if (print) { Pi(id) Pi(m) Pi(w) Pi(am) Pi(pm) Pi(seg1) Pi(seg2) Pi(segs) NL }
	maxa = maxb = -1;
	for (k=seg1; k<seg2; k++) // for each segment (NB ends = 1..N)
	{ int	ak = ends[0][k]-1, pk = ends[1][k]-1, ss = 0;
	  float bias;
		if ((ak>=am && ak<=pm) || (pk<=pm && pk>=am)) continue;
		b = pk-ak+1;
		sec1 = sec2 = mixed = 0;
		for (i=ak; i<=pk; i++) sec1 += sec[i];
		for (i=am; i<=pm; i++) sec2 += sec[i];
		if (sec1<0 && sec2>0) mixed =  1;  // 1=beta, 2=alpha
		if (sec2<0 && sec1>0) mixed = -1;  // 2=beta, 1=alpha
		boxsum = 0.0;
		for (j=ak; j<=pk; j++) { // sum pairs in box (not self)
			for (i=am; i<=pm; i++) {
				boxsum += con[j][i];
				if (sec[i] && sec[j]) ss++;
			}
		}
		if (print) { SPP SPP Pi(k) Pi(ak) Pi(pk) Pi(a) Pi(b) Pr(boxsum) Pi(ss) NL }
		if (boxsum < 1.0) continue; // too weak to score
		if (ss==0) continue;	// no pverlapping helix prediction
		boxsum=1.0;	// don't use in score (optional)
		if (a >= b) {
			for (j=0; j<a; j++) box[j] = xob[a-j-1] = con[j+am]+ak;
			if (print) { for (i=0;i<a;i++){ for(j=0;j<b;j++) printf("%3d",(int)(box[i][j]*100.0)); NL} NL }
			bias = moment(box,a,b,one,mixed);	// bias<0 = anti
			if (one) { // copy the diagonal scores into the full plot
				if (bias>0.0) {
					pack = stripe(box,a,b,dip,one,mixed);
					for (i=0; i<a; i++) for(j=0; j<b; j++) mat[i+am+1][j+ak+1] = dip[i][j];
				} else {
					pack = stripe(xob,a,b,dia,one,mixed);
					for (i=0; i<a; i++) for(j=0; j<b; j++) mat[i+am+1][j+ak+1] = dia[a-i-1][j];
				}
			}
		} else {
			mixed = -mixed;
			for (j=0; j<b; j++) box[j] = xob[b-j-1] = con[j+ak]+am;
			if (print) { for (i=0;i<b;i++){ for(j=0;j<a;j++) printf("%3d",(int)(box[i][j]*100.0)); NL} NL }
			bias = moment(box,b,a,one,mixed);
			if (one) { // copy the diagonal scores into the full plot
				if (bias>0.0) {
					pack = stripe(box,b,a,dip,one,mixed);
					for (i=0; i<b; i++) for(j=0; j<a; j++) mat[j+am+1][i+ak+1] = dip[i][j];
				} else {
					pack = stripe(xob,b,a,dia,one,mixed);
					for (i=0; i<b; i++) for(j=0; j<a; j++) mat[j+am+1][i+ak+1] = dia[b-i-1][j];
				}
			}
		}
		if (print) { Pr(bias) Pr(pack) NLL }
		if (one) { char pa; // score for single box (-/+ = ant/par)
			if (bias<0.0) pa = 'A'; else pa = 'P';
			fprintf(out,"%5d %5d   %f %f %f  %c %f\n", si+1,sj+1, boxsum,ants,pars,pa,bias,pack);
			if (bias<0.0) pack = -pack;
			return pack;
		}
		score = bias;
		if (score < 0.0) score = -score;
		if (score > conmax) { // hold best box score and its size
			conmax = score;
			maxa = a; maxb = b;
		}
		consum += score;
	}
	if (one) return 0.0;
	//condif = conmax - (consum-conmax);	// subtract scores in other boxes
	condif = conmax/(1.0+10.0*(consum-conmax));	// divide by scores in other boxes
	condif /= (1.0 + 1.0*(float)abs(maxa-maxb));
	if (print) { Pi(am) Pi(pm) Pi(maxa) Pi(maxb) Pr(consum) Pr(conmax) Pr(condif) NLL }
	return condif;
}

float moment ( float **boxin, int ain, int bin, int seg, int sec )
// calculate the inertial moment ratio for the points in the box
// NB boxin is an array of pointers to the raw data (do not change)
// sec<0: a=beta, sec>0 b=beta
{
int	i,j, a,b, m,n;
float	score, sum;
float	midX, midY, cenX,cenY, momh,momg; // centroid and moments about each diagonal
float	**box;
	box = (float**)alloca(sizeof(float*)*99);
	for (i=0; i<99; i++) {
		box[i] = (float*)alloca(sizeof(float)*99);
		for (j=0; j<99; j++) box[i][j] = 0.0;
	}
	a = ain; b = bin;
	m = 0;
	for (i=0; i<a; i++) {
		n = 0;
		for (j=0; j<b; j++) {
			box[m][n] = boxin[i][j];
			if (sec<0) n++; n++;
		}
		if (sec>0) m++; m++;
	}
	a = m; b = n;
	cenX = cenY = sum = 0.0;
	midX = 0.5*(float)a; midY = 0.5*(float)b;
	for (i=0; i<a; i++) {
		for (j=0; j<b; j++)
		{ float x,y, dd,g, w = box[i][j];
			if (w < NOISE) continue;
			x = (float)i-midX; y = (float)j-midY;
			dd = x*x+y*y;	// sq.d to crntre of box
			g = exp(-dd*0.01);
			cenX += g*w*(float)i;
			cenY += g*w*(float)j;
			sum += g*w;	// centre weighted sum
		}
	}
	cenX /= sum; cenY /= sum;	// weighted centroid of box
	sum = momh = momg = NOISE;	// avoid -0.0000      
	for (i=0; i<a; i++) {
		for (j=0; j<b; j++)
		{ float x,y, dd,hh,gg,xx,yy, w = box[i][j];
			if (w < NOISE) continue;
			x = (float)i-cenX; xx = x*x;
			y = (float)j-cenY; yy = y*y;
			dd = xx+yy;
			hh = dd*0.5 - x*y;
			gg = dd*0.5 + x*y;
			momh += hh*w;
			momg += gg*w;
			sum += w;
		}
	}
	momh = sqrt(momh/sum)+1.0;
	momg = sqrt(momg/sum)+1.0;
	if (momh > momg) {
		score = -momh/momg +1.0;
	} else {
		score =  momg/momh -1.0;
	}
	return score;
}

float stripe ( float **boxin, int ain, int bin, float dig[99][99], int seg, int sec )
// score the diagonal stripe for each full diagonal (NB: a>b)
// NB: <inbox> points into <con> so a change to <box> -> change to <con>
// sec<0: a=beta, sec>0 b=beta
{
float	mint = 999.9;
int	i,j, k, m,n, mini,minj;
int	iup,jup, idn,jdn;
int	e1i,e1j, e2i,e2j;
int	a,b, mida,midb;
float   **box, mark;
	box = (float**)alloca(sizeof(float*)*99);
	for (i=0; i<99; i++) {
		box[i] = (float*)alloca(sizeof(float)*99);
		for (j=0; j<99; j++) box[i][j] = 0.0;
	}
	// copy box, adding padding if single beta
	a = ain; b = bin;
	m = 0;
	for (i=0; i<ain; i++) {
		n = 0;
		for (j=0; j<bin; j++) {
			box[m][n] = boxin[i][j];
			if (sec>0) n++; n++;
		}
		if (sec<0) m++; m++;
	}
	a = m; b = n;
	// find minimum (weighted) inertial axis
	mida=a/2; midb=b/2;
	iup = idn = mida; jup = jdn = midb;
	while (1) { float t;	// for each full diagonal
		t = momat(iup,jup,box,a,b);
		if (t<mint) { mint=t; mini=iup; minj=jup; }
		iup++;
		t = momat(iup,jup,box,a,b);
		if (t<mint) { mint=t; mini=iup; minj=jup; }
		jup--;
		if (jdn-idn >= b/2) break; 
		idn--;
		t = momat(idn,jdn,box,a,b);
		if (t<mint) { mint=t; mini=idn; minj=jdn; }
		jdn++;
		t = momat(idn,jdn,box,a,b);
		if (t<mint) { mint=t; mini=idn; minj=jdn; }
	}
	if (seg==0) return mint;
	for (i=0; i<a; i++) for (j=0; j<b; j++) dig[i][j] = 0.0;
	for (i=0; i<a; i++) {
		for (j=0; j<b; j++) // add marks on the diagonal
		{ float x,y, dd,hh,gg,g,xx,yy, w = box[i][j];
			if (w < NOISE) continue;
			x = (float)(i-mini); xx = x*x;
			y = (float)(j-minj); yy = y*y;
			dd = xx+yy;
			hh = dd*0.5 - x*y;
			gg = dd-hh;
			g = sqrt(gg*0.5);
			k = (int)(g+0.5);
			xx = (float)(i*i+j*j);
			yy = (float)(mini*mini+minj*minj);
			if (xx < yy) k = -k; 
			dig[mini+k][minj+k] = 2.0;
		}
	}
	// run centre point out to each edge
	a--; b--;
	e1i = mini; e1j = minj;
	while (1) { e1i--; e1j--; if (e1i==0 || e1j==0) break; } 
	e2i = mini; e2j = minj;
	while (1) { e2i++; e2j++; if (e2i==a || e2j==b) break; } 
	// protect value (if over 1) then add mark at line endpoints
	mark = (float)(-seg); // seg = sort key for lines
	dig[e1i][e1j] = (float)((int)(dig[e1i][e1j]+NOISE));
	dig[e2i][e2j] = (float)((int)(dig[e2i][e2j]+NOISE));
	if (dig[e1i][e1j] > 0.9) dig[e1i][e1j] = -1000.0*dig[e1i][e1j];
	if (dig[e2i][e2j] > 0.9) dig[e2i][e2j] = -1000.0*dig[e2i][e2j];
	dig[e1i][e1j] += mark-0.1; dig[e2i][e2j] += mark-0.2; // 0.1-0.2 used to sort ends
	if (sec==0) return mint;
	// strip beta padding from dig
	m = 0;
	for (i=0; i<ain; i++) {
		n = 0;
		for (j=0; j<bin; j++) {
			dig[i][j] = dig[m][n]; 
			if (sec<0) dig[i][j] += dig[m+1][n]; // add in lost row
			n++;
			if (sec>0) dig[i][j] += dig[m][n++]; // add in lost col
		}
		if (sec<0) m++; m++;
	}
	return mint;
}

float momat ( int p, int q, float **box, int a, int b )
// calculate the inertial moment points in the box about p,q
{
int	i,j;
float	mom = 0.0, sum = 0.0; // centroid and moments about each diagonal
	for (i=0; i<a; i++) {
		for (j=0; j<b; j++)
		{ float x,y, dd,hh,g,xx,yy, w = box[i][j];
			if (w < NOISE) continue;
			x = (float)(i-p); xx = x*x;
			y = (float)(j-q); yy = y*y;
			dd = xx+yy;
			hh = dd*0.5 - x*y;
			g = exp(-hh*0.01);
			mom += g*w*hh;
			sum += g*w;
		}
	}
	mom = sqrt(mom/sum)+1.0;
	return mom;
}
