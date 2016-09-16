/*
c++ -ggdb score.cpp -o score sims/util.o sims/geom.o -lm

*/
#include "sims/util.hpp"
#include "sims/geom.hpp"

int main (int argc, char** argv) {
char	line[111];
float	avd[555], nad[555];
Vec	ave[555], nat[555];
FILE	*pdb = fopen(argv[1],"r");
int	len, n, m;
	// read native (lies at the end as chain B)
	n = 0;
	DO
	{ float x,y,z;
	  int	io = read_line(pdb,line);
		if (io <= 0) break;
		if (line[21] != 'B') continue;
		sscanf(line+26,"%f %f %f", &x, &y, &z);
		nat[n].x = x; nat[n].y = y; nat[n].z = z; 
		n++;
	}
	len = n;
	Pi(len) NL
	fclose(pdb);
	// read models
	n = 0; m = 0;
	FOR(i,555) ave[i].zero();
	pdb = fopen(argv[1],"r");
	DO
	{ float x,y,z;
	  int	io = read_line(pdb,line);
		if (io <= 0) break;
		if (line[21] == 'B') break;
		if (line[0] == 'T') { m++; n = 0; continue; }
		if (line[0] != 'A') continue;
		sscanf(line+26,"%f %f %f", &x, &y, &z);
		ave[n].x += x; ave[n].y += y; ave[n].z += z; 
		n++;
	}
	Pi(m) NL
	FOR(i,len) ave[i] /= (float)m;
	fclose(pdb);
	// reread models
	n = 0;
	FOR(i,555) avd[i] = nad[i] = 0.0;
	pdb = fopen(argv[1],"r");
	DO
	{ float x,y,z, dx,dy,dz;
	  int	io = read_line(pdb,line);
		if (io <= 0) break;
		if (line[21] == 'B') break;
		if (line[0] == 'T') { n = 0; continue; }
		if (line[0] != 'A') continue;
		sscanf(line+26,"%f %f %f", &x, &y, &z);
		dx = x - nat[n].x; dx *= dx;
		dy = y - nat[n].y; dy *= dy;
		dz = z - nat[n].z; dz *= dz;
		nad[n] += (dx+dy+dz);
		dx = x - ave[n].x; dx *= dx;
		dy = y - ave[n].y; dy *= dy;
		dz = z - ave[n].z; dz *= dz;
		avd[n] += (dx+dy+dz);
		n++;
	}
	FOR(i,len) {
		avd[i] /= (float)m; avd[i] = sqrt(avd[i]);
		nad[i] /= (float)m; nad[i] = sqrt(nad[i]);
		Pi(i) Pr(avd[i]) Pr(nad[i]) NL
	}
}
