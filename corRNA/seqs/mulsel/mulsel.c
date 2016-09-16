/*
                               MULTAS
                Multiple sequence alignment program,
                copyright William Ramsay Taylor 1987

cc -O mulsel.c -o mulsel seqsort.o util/wt/util.o util/wt/sort.o -lm

*/
#include "util/wt/incl/util.h"
#include "util/multal/multal.h"
#include "util/multal/gcgseq.h"

#define SEQ 100000
#define TITLE 200
#define TREE 100
#define CODE 100

#define MAX 999999999
#define MIN -MAX


/*	STRUCTURES FOR ALIGNMENTS	*/

typedef	struct	{
		char	*residue,
			title[TITLE+1],
			code[CODE+1],
			tree[TREE+1];
		int	id, odd,
			length;
		}
	Seqs;
Seqs	**sequence;
int  	nseqs;

typedef struct	{
		unsigned char	key;
		int	members, length, grouped;
		int	score_to_last;
		Seqs	**seq;
		}
	Blocks;
Blocks	**block;
int	nblocks, nkey;

typedef	struct	{
		int	last, member[MAXLISLEN];
		}
	Halfs;
typedef struct {
		Halfs	right[1], left[1];
		int	closed, dead;
		}
	Lists;
Lists	*list;
int	nlists, list_alloc;

typedef struct	{
		Blocks	*block_a, *block_b;
		int	position_a[MAXALNLEN], position_b[MAXALNLEN],
			last_link[MAXALNLEN], next_link[MAXALNLEN],
			length;
		int	score;
		}
	Alns;	
Alns	*aligned_pairs;

typedef struct	{
		int	members;
		Blocks	*blk[MAXLISLEN*2];
		}
	Groups;
Groups	*group;
int	ngroups, group_alloc;

typedef struct	{
		int	block_a, block_b;
		}
	Pairs;
Pairs	*pair;
int	npairs, pair_alloc;

struct sequence seq;

char	*alignment;
int	pat[MAXWINLEN][MAXALNLEN];
int	maxcol, mincol, colmax[MAXNALN],
	matrix_a[NACID][NACID],
	matrix_b[NACID][NACID],
	matrix[NACID][NACID],
	nres = 0, mat_wt = 0,
	table, gap_pen, window, mul_wt, damp,
	span, minscore, list_limit, old_list_limit = 0,
	pep_span, last_span, circles, output, outline;
int	aln_len, pep_len, pepscore, soft,
	region = 0, delayed = 0;
FILE	*run_data;

int	seedseq, seedlen, lowscore;

main (argc, argv) int argc; char *argv[];
{
	seedseq = seedlen = lowscore = 0;
	if (argc > 1) sscanf(argv[1],"%d", &seedlen);
	if (argc > 2) sscanf(argv[2],"%d", &lowscore);
	Pi(seedlen) Pi(lowscore) NL
	matin("data/md.mat",matrix_b);
	matin("data/id.mat",matrix_a);
	seqin();
	last_span = nblocks;
	if (pep_len)
	{ int	*pep_order;
	  int	i;
		pep_order = (int*)malloc(sizeof(int)*(nblocks+1));
		TEST(pep_order)
		for (i=0; i<=nblocks; i++) pep_order[i] = i;
		pept_out(nblocks);
		for (i=0; i<abs(pep_span); i++)
		{ int	aspan, bspan, minp, hold, diag;
			fscanf(run_data,"%d%d%d%d%d", &aspan,&bspan,&minp,&hold,&diag);
			if (aspan < 0) aspan = nblocks/(i+2);
			if (bspan < 0) bspan = last_span;
			npairs = count(nblocks,aspan,bspan,pep_order,minp,hold,diag);
			last_span = aspan;
			pre_order(pep_order);
			pep_order[0]++;
		}
		if (pep_span < 0) exit(1);
	}
	while(align());
	write_one_seq();
	exit(1);
}

write_one_seq()
{	int	i, j, k;
	int	total, kept, seedseq, seedblk;
	char	*c;
	FILE	*aln, *lis, *tre;
	aln = fopen("final.seq","w");
	lis = fopen("final.lis","w");
	tre = fopen("final.tre","w");
	for ( i = 0; i < nblocks; i++ ) {
                for ( j=0; j<block[i]->members; j++ ) { Seqs  *s = block[i]->seq[j];
			if (strstr(s->code,"SEED:")) {
				seedblk = i;
				seedseq = j;
				break;
			}
		}
	}
	Pi(nblocks) Pi(seedblk) Pi(seedseq) NL
	kept = 0;
	for ( i = 0; i < nblocks; i++ ) 
	{ 	Blocks	*b = block[i], *b0 = block[0];
		int	b_memb = b->members, n = 1;
		float	mean, best = 999.9;
		int	nseed, seq = -1;
		mean = 0.0;
		if (!b_memb) continue;
		if (lowscore) { int s;
			s = score( block[seedblk], b, 0 )*100/block[seedblk]->length;
			if ( s<lowscore ) {
				printf("Block %d score %d below lowscore %d\n", i,s,lowscore);
				continue;
			}
		}
		for ( j=0; j<b_memb; j++ ) mean += (float)b->seq[j]->length;
		mean /= (float)b_memb;
		if (seedlen) {
			printf("Filtering to match seed length (%d)\n", seedlen);
			mean = (float)seedlen;
		}
		for ( j=0; j<b_memb; j++ )
		{ Seqs	*s = b->seq[j];
                  char  ends[10], resn[10], bval[10];
                  int   nend, cend, nres;
		  float resol, bvalue;
			fprintf(lis,"%s =%d= %s\n", s->code,s->odd,s->title);
			fprintf(tre,"%s %s =%d= %s\n",
				s->tree,s->code,s->odd,s->title);
			resol = mean - (float)(s->length);
			resol = log(resol*resol+1.0);
                        resol += (float)s->odd;
                        if (strstr(s->title,"PROBABLE")) resol += 1.0;
                        if (strstr(s->title,"probable")) resol += 1.0;
                        if (strstr(s->title,"PRECURSOR")) resol += 2.0;
                        if (strstr(s->title,"precursor")) resol += 2.0;
                        if (strstr(s->title,"HYPOTHETICAL")) resol += 5.0;
                        if (strstr(s->title,"hypothetical")) resol += 5.0;
                        if (strstr(s->title,"FRAGMENT")) resol += 50.0;
                        if (strstr(s->title,"fragment")) resol += 50.0;
                        if (strstr(s->title,"MUTANT")) resol += 40.0;
                        if (strstr(s->title,"mutant")) resol += 40.0;
                        if (strstr(s->code,"SEED:")) resol -= 100.0;
                        if (strstr(s->code,"seed:")) resol -= 40.0;
                        if (strstr(s->code,"DOM:")) resol -= 50.0;
                        if (strstr(s->code,"CAS:")) resol -= 70.0;
                        if (strstr(s->code,"PDB:")) resol -= 60.0;
                        if (strstr(s->code,"PDB;")) resol -= 60.0;
                        if (strstr(s->code,"pdb;")) resol -= 60.0;
                        if (strstr(s->code,"pdb:")) resol -= 60.0;
                        if (strstr(s->code,"|pdb|")) resol -= 60.0;
                        if (strstr(s->code,"SPROT:")) resol -= 20.0;
                        if (strstr(s->code,"|sp|")) resol -= 20.0;
			if (resol<best) { best = resol; seq = j; }
		}
		if (seq<0) continue;
		if (best>50.0) continue;
		fprintf(aln,"%s\n", b->seq[seq]->code+5);
		printf("picked %s (%f)\n", b->seq[seq]->code,best);
		total = -1;
		if (strstr(b->seq[seq]->code,"SEED:")) total = 0;
		if (strstr(b->seq[seq]->code,"seed:")) total = 0;
		for ( j=0; j<b_memb; j++ )
		{ Seqs	*s = b->seq[j];
		  int	more;
			total++;
			if (strstr(s->code,"SEED:")) total--;
			if (strstr(s->code,"seed:")) total--;
			c = (char*)strstr(s->title,"[+");
			if (!c) continue;
			sscanf(c+2,"%d", &more);
			total += more;
		}
		fprintf(aln,"[+%d+] ", total);
		c = (char*)strstr(b->seq[seq]->title,"+]");
		if (c) fprintf(aln,"%s\n", c+3);
			else fprintf(aln,"%s\n", b->seq[seq]->title);
		for ( j = 1; j <= b->length; j++ )
		{	Seqs	*s = b->seq[seq];
			char	res = s->residue[j];
			if (res==' ') continue;
			fprintf(aln,"%c",res); 
		}
		fprintf(aln,"*\n");
		kept++;
	}
	Pi(kept) NL
}

matin(file,mat)
	char	*file;
	int	mat[NACID][NACID];
{
	int	i, j, mat_const;
	char	acid[NACID], c;
	FILE	*mat_file;

	mat_file = fopen(file,"r");
	while( c = getc(mat_file), c != '\n' ) putchar(c); NL
	fscanf(mat_file,"%s\n",acid);
	printf("%s\n",acid);
	fscanf(mat_file,"%d\n",&mat_const);
	printf("matrix constant = %d\n",mat_const);
	for( i = 0; acid[i]; i++ ) 
	{	int	ai = acid[i]-'A';
		for( j = 0; acid[j]; j++ ) 
		{	int aj = acid[j]-'A';
			fscanf(mat_file,"%d",&mat[ai][aj]);
			mat[ai][aj] += mat_const;
			matrix[ai][aj] = mat[ai][aj];
		}
	}
	matrix[ 9][ 9] = 100;   /* set JJ to score 100 */
	matrix[14][14] = 500;   /* set OO to score 500 */
	matrix[20][20] = 1000;  /* set UU to score 1000 */
	nres = max(nres, i+1);
 }
		
seqin ()
{
	int	i, members, status,
		locseq, locref, istat,
		maxs, mins, key,
		inseq, user_in, more_seqs = 1;
	char	user_file[40], user_name[40],
		*file_name = user_file,
		string[10];
	FILE	*user_seqs = 0, *user_list = 0;

/* READ PARAMETER FILE HEADER */
	run_data = fopen("skip.run","r");
	fscanf(run_data,"%s %d %d %d %d %d %d %d",
		user_file,&inseq,&mins,&maxs,&pep_len,&table,&pep_span,&soft);
	maxs = min(maxs,SEQLENMAX);
	mins = max(mins,SEQLENMIN);
	if (mins<pep_len) {
		mins = pep_len;
		printf("Min.seq.len set to pep.len\n");
	}
	if (inseq<0) { region = 1; inseq = -inseq; }
	printf("\n%d sequences to be read in length range %d -> %d\n",
		inseq,mins,maxs);
	if (pep_len && delayed<2) init_slots(inseq,pep_len,soft);
/*
	if (inseq>TABLE) {
		printf("More than %d sequences so no names printed\n",TABLE);
		table = 0;
	} else {
		if (!table) printf("No names will be printed\n");
	}
*/
	if (!delayed && inseq>DELAY) delayed = 1;
	if (delayed==1) fclose(run_data);
	if (delayed<2) {
		block = (Blocks**)malloc(sizeof(Blocks*)*inseq);
		sequence = (Seqs**)malloc(sizeof(Seqs*)*(inseq+1));
	}
	user_in = check_code(user_file,&key);
	switch(user_in) {
		case 0:
			printf("Reading user file of names: %s\n",file_name+1);
			user_list = fopen(file_name+1,"r");
			if (!user_list) 
				printf("Cant open file: %s\n",file_name+1);
			break;
		case 1 :
			printf("Reading user file of seq.s: %s\n",file_name);
			user_seqs = fopen(file_name,"r");
			if (!user_seqs) 
				printf("Cant open file: %s\n",file_name);
			user_list = 0;
			break;
		case 2 :
			printf("Extracting databank sequences %s\n",file_name);
		case 3 :
			printf("Reading user alignment file %s\n",file_name);
	}
/* READ SEQUENCES */
	nblocks = nseqs = nkey = 0;
	while (nseqs<inseq && more_seqs) {
Pi(user_in) NL
		switch(user_in) {
		case 1 : 
			while(get_user_seq(key,user_seqs,maxs,mins,inseq));
			printf("End of user sequences\n");
			break;
		case 2 :
			get_bank_seq(key,file_name,maxs,mins,inseq);
			break;
		case 3 :
			while(get_user_aln(key,user_seqs,maxs,mins,inseq));
			printf("End of user alignments\n");
		}
		if (!user_list) break;
		fscanf(user_list,"%s", user_name);
		if (next_line(user_list)<0 || feof(user_list)) {
			printf("End of user list\n");
			more_seqs = 0;
			break;
		}
		if (nseqs==inseq) break;
		printf("Reading list: %s\n",user_name);
		user_in = check_code(user_name,&key);
		if (user_in==2) {
			file_name = user_name;
		} else {
			file_name = user_name+5;
			if (user_seqs) fclose(user_seqs);
			user_seqs = fopen(file_name,"r");
			if (!user_seqs)
				printf("Cant open file: %s\n",file_name);
		}
	}
	printf("%d sequences read into %d blocks, (%d are keyed)\n", 
		nseqs,nblocks,nkey);
	if (user_list) fclose(user_list);
	if (user_seqs) fclose(user_seqs);
}

check_code (file,key)
char	*file;
int	*key;
{
	int	in = 1;
	*key = 0;
	if (file[0] == '@') return 0;
	if (!strncmp(file,"USER:",5)) { *key =  1; return  1; }
	if (!strncmp(file,"User:",5)) { *key =  1; return  1; }
	if (!strncmp(file,"user:",5)) { *key =  0; return  1; }
	if (in) in = strncmp(file,"PATCHX:",7);
	if (in) in = strncmp(file,"LEEDS:",6);
	if (in) in = strncmp(file,"SWISS:",6);
	if (in) in = strncmp(file,"NBRF:",5);
	if (in) in = strncmp(file,"NEW:",4);
	if (!in) { *key =  1; return  2; }
	if (in) in = strncmp(file,"Patchx:",7);
	if (in) in = strncmp(file,"Leeds:",6);
	if (in) in = strncmp(file,"Swiss:",6);
	if (in) in = strncmp(file,"Nbrf:",5);
	if (in) in = strncmp(file,"New:",4);
	if (!in) { *key =  1; return  2; }
	if (in) in = strncmp(file,"patchx:",7);
	if (in) in = strncmp(file,"leeds:",6);
	if (in) in = strncmp(file,"swiss:",6);
	if (in) in = strncmp(file,"nbrf:",5);
	if (in) in = strncmp(file,"new:",4);
	if (!in) { *key =  0; return  2; }
	if (!strncmp(file,"ALNS:",5)) { *key =  1; return  3; }
	if (!strncmp(file,"Alns:",5)) { *key =  1; return  3; }
	if (!strncmp(file,"alns:",5)) { *key =  0; return  3; }
	*key = 1;
	return 1;
}

get_bank_seq(key,code,maxs,mins,in)
char	*code;
int	key, maxs, mins, in;
{
/*
	int	n = 0;
	while (sqnext_(code,&seq)) {
		n++;
		if (seq.len < mins) continue;
		if (seq.len > maxs) continue;
		load_seq(key,mins,maxs);
		if (nseqs==in) {
			printf("*NB* list not exhausted\n");
			return;
		}
	}
	if (!n) printf("No names fit!  Check GCG is set up\n");
*/
printf("Check GCG is not set up\n");
}

load_seq(key,mins,maxs)
int	key, mins, maxs;
{	Blocks	*b;
	Seqs	*s;
	int	j, print = 0;
	if(seq.len > maxs) return 0;
	if(seq.len < mins) return 0;
	if (nseqs && !print && delayed!=1) {
		if (sequence[nseqs]->title[0] != seq.doc[0]) print = 1;
	}
	nseqs++;
	if (table) printf("%d %s %d %s\n", nseqs,seq.name,seq.len,seq.doc);
	if (delayed==2) b = block[nblocks];
	else {
		b = block[nblocks] = (Blocks*)malloc(sizeof(Blocks));
		TEST(b)
	}
	nblocks++;
	if (delayed==1) b->seq = 0;
	else {
		b->seq = (Seqs**)malloc(sizeof(Seqs*));
		TEST(b->seq)
		s = sequence[nseqs] = (Seqs*)malloc(sizeof(Seqs));
		TEST(s)
		s->residue = (char*)malloc(sizeof(char)*(seq.len+1));
		TEST(s->residue)
		seq.name[CODE] = '\0';
		strcpy(s->code,seq.name);
		seq.doc[TITLE] = '\0';
		strcpy(s->title,seq.doc);
		s->length = seq.len;
		s->tree[0] = '\0';
		s->id = nblocks;
		b->seq[0] = s;
	}
	b->length = seq.len;
	b->members = 1;
	b->grouped = 0;
	b->key = key;
	if (key) nkey++;
	if (delayed==1) copy_seq(key, s->code, seq.seq, seq.seq, seq.len);
	else	s->odd = (int)copy_seq(key, s->code, seq.seq, s->residue,  seq.len);
}

copy_seq (key,code,res_in,res_out,length) 
char	*code, *res_in, *res_out;
int	length;
int	key;
{
	int	j, odd = 0, glitch = 0, plen = abs(pep_len);
	for(j=1; j <= length; j++) 
	{	char	res = res_in[j-1],
			upr = UPPER(res);
		if (res=='*') {
			length = j-1;
			return;
		}
		res_out[j] = res;
		if (res==' ') {
			glitch = j;
			continue;
		}
		if ( upr < 'A' || upr > 'Z' ) {
			printf("*NB* bad res >%c< in %s\n", res,code);
			res_out[j] = '?';
			glitch = j;
			odd += 10;
			continue;
		}
		if (strchr("JOUX",upr)) odd++;
		if (strchr("BZ",upr)) odd++;
		if (!plen) continue;
		if ( region && isupper(res_out[j])) {
			glitch = j;
			continue;
		}
		if (j-plen<glitch) continue;
		if (delayed==2) continue;
		find_slot(key,plen,res_out+(j-plen+1),nblocks);
	}
	return odd;
}

get_user_seq (key,user,maxs,mins,in)	
FILE	*user;
int 	key, maxs, mins, in;
{
	char	c, line[10000];
	int	i, n=0,
		coff=5;
	if (nseqs==in) {
		printf("*NB* seq limit reading user file\n");
		return 0;
	}
	seq.len = 0;
	while (1) { int io = read_line(user,line);
		if (io==0) continue;
		if (io<0) break;
		if (*line != '>' && strncmp(line,"SEED:",5)) continue;
		if (key) strncpy(seq.name,"USER>",coff);
		   else  strncpy(seq.name,"user>",coff);
		seq.name[coff] = (char)0;
		strncat(seq.name,line,CODE-coff);
		read_line(user,line);
		strncpy(seq.doc,line,TITLE);
		for ( i=0; i<5000; i++ ) {
			while(isspace((c=getc(user))));
			if (c=='*') {
				seq.len = n;
				load_seq(key,mins,maxs);
				seq.len = 0;
				return 1;
			}
			if (isalpha(c)) seq.seq[n++] = c;
		}
	}
	if (seq.len>0) load_seq(key,mins,maxs);
	return 0;
}
	
get_user_aln(key,user,maxs,mins,in)	
FILE	*user;
int 	key,maxs, mins, in;
{
	Blocks	*b;
	int	i, j, n,
		members,
		coff=5,
		next_aln = 0;
	char	line[500];
	if (read_line(user,line) < 0) return 0;;
	key = 999;
	if (strstr(line,"BLOCK")) key = -1;
	if (strstr(line,"block")) key =  0;
	if (strstr(line,"Block")) key =  1;
	if (key==999) return 1;
	read_line(user,line);
	sscanf(line,"%d",&members);
	Pi(members) NL
	if (members<2) {
		printf("*NB* too few members in USER alignment\n");
		return 0;
	}
	b = block[nblocks] = (Blocks*)malloc(sizeof(Blocks));
		TEST(b)
	nblocks++;
	b->seq = (Seqs**)malloc(sizeof(Seqs*)*members);
		TEST(b->seq)
	b->members = members;
	b->grouped = 0;
	b->key = key;
	if (key) nkey++;
	for (i=0; i<members; i++)
	{	Seqs	*s;
		char	*c;
		if (nseqs==in) {
			printf("*NB* seq limit reading alignment\n");
			return 0;
		}
		nseqs++;
		s = sequence[nseqs] = (Seqs*)malloc(sizeof(Seqs));
		TEST(s)
		s->residue = (char*)malloc(sizeof(char)*(maxs+1));
		TEST(s->residue)
		s->residue[0] = 0;
		s->id = nblocks;
		s->tree[0] = '\0';
		b->seq[i] = s;
		read_line(user,line);
		n = 0;
		if (!n && strstr(line," = ")) {
			c = (char*)strstr(line," = ")+3;
			if (strlen(c) > TITLE) c[TITLE] = (char)0;
			strcpy(s->title,c);
			c = (char*)strstr(line," = ");
			*c = (char)0;
			n = 1;
		}
		if (!n && strstr(line," : ")) {
			c = (char*)strstr(line," : ")+3;
			if (strlen(c) > TITLE) c[TITLE] = (char)0;
			strcpy(s->title,c);
			c = (char*)strstr(line," : ");
			*c = (char)0;
			n = 1;
		}
		if (!n && strchr(line,' ')) {
			c = (char*)strstr(line," ")+1;
			if (strlen(c) > TITLE) c[TITLE] = (char)0;
			strcpy(s->title,c);
			c = (char*)strstr(line," ");
			*c = (char)0;
			n = 1;
		}
		if (!n) s->title[0] = (char)0;
		strncpy(s->code,"ALNS>",5);
		s->code[5] = 0;
		if (strlen(line) > CODE-5) line[CODE-5] = (char)0;
		strcat(s->code,line);
		printf("%d %s : %s\n",nseqs,s->code,s->title);
	}
	n = 1;
	while (!feof(user)) 
	{	char	line[MAXNALN];
		read_line(user,line);
		if (!strncmp(line,"****",3)) {
			printf("End of Alignment\n");
			next_aln = 1;
			break;
		}
		for ( i=0; i<members; i++)
		{	Seqs	*s = b->seq[i];
			s->residue[n] = line[i];
			if (s->residue[n]==GAP) s->residue[n] = ' ';
			if (s->residue[n] < ' ') { n--; break; }
			if ( n>maxs ) {
				printf("*NB* alignment too long\n");
				n = maxs;
				break;
			}
		}
		if (n==maxs) break;
		n++;
	}
	b->length = n-1;
	if (key) printf("Key ");
	printf("alignment length = %d\n",b->length);
	for (i=0; i<members; i++)
	{	Seqs	*s = b->seq[i];
		s->length = b->length;
		copy_seq(key,s->code,s->residue+1,s->residue,s->length);
		if (delayed==1) {
			free(s->residue);
			free(s);
		}
	}
	return next_aln;
}

get_param ()
{
	int	listlim,
		old_mat_wt = mat_wt;

	if ( fscanf(run_data,"%d%d%d%d%d%d%d%d%d%d", &mul_wt,&damp,&mat_wt,
		&gap_pen,&span,&window,&minscore,&pepscore,&listlim,&output)
		== -1 ) return 0;
	list_limit = MAXLISLEN;
	if (listlim) list_limit = abs(listlim);
	if ( list_limit > MAXLISLEN ) list_limit = MAXLISLEN;
	if ( list_limit < 2 ) list_limit = 2;
	circles = 0;
	if (listlim<0) circles = 1;
	outline = 255;
	if (output<0) outline = 80;
	output = abs(output);
	printf("\n\nPARAMETERS for next cycle\n");
	Pi(mul_wt) Pi(damp) Pi(mat_wt) NL
	Pi(gap_pen) Pi(span) Pi(window) NL
	Pi(minscore) Pi(pepscore) Pi(list_limit) NL
	Pi(circles) Pi(output) Pi(outline) NL
	if (window>MAXWINLEN) {
		window = MAXWINLEN;
		printf("*NB* window width has been reduced to %d\n",window);
	}
	if ( mat_wt != old_mat_wt ) revise_matrix();
	return 1;
}

pre_order(pep_order)
int	*pep_order;
{
	int	i, j, n, m;
	int	*pep_pair;
	int	*place;
	Blocks	**new;
	char	*dead;
printf("pre_order ");
	pep_pair = (int*)malloc(sizeof(int)*npairs);	TEST(pep_pair)
	place = (int*)malloc(sizeof(int)*npairs);	TEST(place)
	pair = (Pairs*)malloc(sizeof(Pairs)*npairs);	TEST(pair)
printf(" malloc OK\n");
	reorder(pep_pair);
	for (n=0; n<npairs; n++) {
		unpack(&i,&j,nblocks,pep_pair[n]);
		pair[n].block_a = i-1;
		pair[n].block_b = j-1;
		place[n] = n;
	}
printf("preorder ");
	free(pep_pair);
/* CLUSTER */
	new = (Blocks**)malloc(sizeof(Blocks*)*nblocks); TEST(new);
printf(" new block mallock OK\n");
	if (pep_order[0]) {
		for ( i = 0; i < nblocks; i++) new[pep_order[i+1]-1] = block[i];
		for ( i = 0; i < nblocks; i++) block[i] = new[i];
	}
	list_alloc = min(LISTALLOC,nblocks/2);
	list = (Lists*)malloc(sizeof(Lists)*(list_alloc+1));	TEST(list)
	printf("cluster_pairs ");
	circles = 0;
	list_limit = MAXLISLEN;
	cluster_pairs(place);
	free(place);
	free(pair);
	group_alloc = 0;
	for ( i = 1; i <= nlists; i++ ) if( !((list+i)->dead) ) group_alloc++;
	group = (Groups*)malloc(sizeof(Groups)*group_alloc);	TEST(group)
	printf("pack_lists ");
	pack_lists();
	free(list);
/* LOOP OVER THE GROUPS */
	n = m = 0;
	dead = (char*)malloc(sizeof(char)*ngroups);	TEST(dead);
	for ( i = 0; i < ngroups; i++ ) dead[i] = 0;
	for ( i = 0; i < nblocks; i++) {
		if ( block[i]->grouped ) 
		{	int	grp = block[i]->grouped-1;
			Groups	*g = group+grp;
			if ( dead[grp] ) continue;
			for ( j = 0; j < g->members; j++ ) {
				new[n] = g->blk[j];
				m++; n++;
			}
			dead[grp] = 1;
		} else  new[n++] = block[i];
	}
	if(n!=nblocks) printf("*NB* lost some sequences in sort!");
	printf("%d/%d blocks moved\n", m,n);
	free(dead);
	if (delayed) { 
		delayed = 2; 
		seqin();
	}
	printf("\nNew sequence order is:\n");
	for ( i=0; i < nblocks; i++)
	{	char c = ' ', d = '-';
		pep_order[i+1] = new[i]->seq[0]->id;
		block[i] = new[i];
		if (block[i]->grouped) d = '='; 
		block[i]->grouped = 0;
		if (block[i]->key) c = '*';
		for ( j=0; j<block[i]->members; j++)
		{	Seqs	*s = block[i]->seq[j];
			printf("%c Block %d <%c%c%c %d\t%s\t%5d   %s\n", c,
				i,d,d,d,s->id,s->code,(int)s->length,s->title);
		}
	}
	printf("\n('*' = keyed block)\n\n");
	free(new);
}	
		
align ()
{
	extern	sort();
	int	*place; 
	int	*scores,
		size, i;

	if (!get_param()) return 0;
/* SCORE */
	size = min(span,nblocks);
	pair_alloc = min(PAIRALLOC,(nblocks*size)-(size*size)/2+1);	
	pair = (Pairs*)malloc(sizeof(Pairs)*pair_alloc);	TEST(pair)
	scores = (int*)malloc(sizeof(int)*pair_alloc);		TEST(scores)
	printf("score_pairs ");
	i = score_pairs(&scores);
	if ( !i ) {
		free(pair);
		free(scores);
		return 1;
	}
	place  = (int*)malloc(sizeof(int)*pair_alloc);	TEST(place)
	printf("(%d pairs) ",npairs);
	sort(0, 0, scores, place, npairs, 1 );
	free(scores);
/* CLUSTER */
	list_alloc = min(LISTALLOC,nblocks/2);
	list = (Lists*)malloc(sizeof(Lists)*(list_alloc+1));	TEST(list)
	printf("cluster_pairs ");
	cluster_pairs(place);
	free(place);
	free(pair);
	group_alloc = 0;
	for ( i = 1; i <= nlists; i++ ) if( !((list+i)->dead) ) group_alloc++;
	group = (Groups*)malloc(sizeof(Groups)*group_alloc);	TEST(group)
	printf("pack_lists ");
	pack_lists();
	free(list);
	printf("\n");
/* ALIGN */
	aligned_pairs = (Alns*) malloc(sizeof(Alns)*list_limit);
			TEST(aligned_pairs)
	for ( i = 0; i < ngroups; i++ ) 
	{	Groups	*g = group+i;
		align_pairs(aligned_pairs,g);
		expand_seqs(aligned_pairs,g);
	}
	printf("update_blocks\n");
	update_blocks();
	free(group);
	free(aligned_pairs);
/* PRINT */
	return print_aln();
}

print_aln()
{
	int	live_blocks = 0,
		aligned = 0,
		i, j, k;

	for ( i = 0; i < nblocks; i++ ) 
	{ 	Blocks	*b = block[i];
		int	b_memb = b->members, n = 1;
		if ( b_memb ) live_blocks++;
		if ( b_memb > 1 ) aligned++;
		if (!b->grouped) continue;
		b->grouped = 0;
		if ( !output ) continue;
		if ( b_memb < output) continue;
		printf("block %d = %d seqs\n",i,b_memb);
		for ( j=0; j<b_memb; j++ ) 
		{	Seqs	*s = b->seq[j];
			printf("*%s :	%s\n",s->code,s->title);
		}
		NL
		while ( n < b->length && output > 0 ) 
		{	Seqs	*top = b->seq[0];
			for ( j = 0; j < b_memb; j++ )
			{	Seqs	*s = b->seq[j];
				int	m = min( n+outline, b->length+1 );
				for ( k = n; k < m; k++ )
				{	char	res = s->residue[k];
/*
					if (j && res!=' ' &&
						 res==top->residue[k]) res='.';
*/
					printf("%c",res); 
				}
				NL 
			} NL
			n += outline; 
		} NL
	}
	printf("%d groups ( %d aligned )\n\n",live_blocks,aligned);
	return live_blocks-1;
}

write_aln()
{	int	i, j, k;
	FILE	*aln, *lis, *tre;
	aln = fopen("final.aln","w");
	lis = fopen("final.lis","w");
	tre = fopen("final.tre","w");
	for ( i = 0; i < nblocks; i++ ) 
	{ 	Blocks	*b = block[i];
		int	b_memb = b->members, n = 1;
		for ( j=0; j<b_memb; j++ )
		{ Seqs	*s = b->seq[j];
			fprintf(lis,"%s =%d= %s\n", s->code,s->odd,s->title);
			fprintf(tre,"%s %s =%d= %s\n",
				s->tree,s->code,s->odd,s->title);
		}
		if ( !output ) continue;
		if ( b_memb < output) continue;
		fprintf(aln,"Block %d\n %d seqs\n",i,b_memb);
		for ( j=0; j<b_memb; j++ )
		{ Seqs	*s = b->seq[j];
			fprintf(aln,"%s = %s\n", s->code,s->title);
		}
		for (k=1; k <= b->length; k++) {
			for ( j = 0; j < b_memb; j++ )
			{	Seqs	*s = b->seq[j];
				char	res = s->residue[k];
				if (res==' ') res = GAP;
				fprintf(aln,"%c",res); 
			}
			fprintf(aln,"\n");
		}
		fprintf(aln,"\n");
	}
	fprintf(aln,"\n");
}

revise_matrix ()
{
	int	wt = abs(mat_wt),
		a_wt = 10-wt,
		b_wt = wt,
		i, j;
	if (mat_wt<0) printf("matrix is\n");
	for ( i=0; i < NACID; i++ ) {
		for ( j=0; j < NACID; j++ ) {
			matrix[i][j] =	a_wt*matrix_a[i][j] +
					b_wt*matrix_b[i][j];
			matrix[i][j] = matrix[i][j]/10;
			if (mat_wt<0) printf("%3d",matrix[i][j]);
		}
		if (mat_wt<0) NL
	}
	if (mat_wt<0) NL
}

/*	ALIGN SEQUENCE SEGMENTS PAIRWISE	*/

score_pairs (scores)
	int	**scores;
{
	int	i, j, ni, nj,
		pep_score = 0,
		sum_scores = 0,
		pscores = 0;
	float	ps = 0.0, sp = 0.0,
		ms = 0.0, mp = 0.0;
	npairs = 0;
	ni = -1;
	NL
	for (i = 0; i < nblocks; i++ ) 
	{	Blocks	*bi = block[i];
		int	nbi = bi->members,
			lbi = bi->length;
		if ( !nbi ) continue;
		nj = ++ni;
		if (table) printf("\n");
		for ( j = 0; j < nblocks; j++ )
		{	Blocks	*bj = block[j];
			int 	nbj = bj->members,
				lbj = bj->length, 
				s;
			if ( !nbj ) continue;
			if (j<i) {
				if (table) printf("    ");
				continue;
			}
			if (j==i) {
				if (table) printf("%5d=block",i);
				continue;
			}
			nj++;
			if (!bi->key && !bj->key) {
				if (table) printf("  ><");
				continue;
			}
			if ( nj-ni > span ) {
				if (table) printf("  >s");
				continue;
			}
			if ( abs(lbi-lbj) > window ) {
				if (table) printf("  >w");
				continue;
			}
			if (pep_len && pepscore)
			{	int	cut = pepscore, k, l;
				for (k=0; k<nbi; k++)
				{	int	sk = bi->seq[k]->id;
					for (l=0; l<nbj; l++)
					{	int	sl = bj->seq[l]->id;
						pep_score=pscore(sk,sl,nblocks);
						pep_score = abs(pep_score);
						if (pep_score<cut) pep_score=0;
						if (pep_score) break;
					}
					if (pep_score) break;
				}
				if (!pep_score) {
					if (table) printf("  <p");
					continue;
				}
			}
			s = score( bi, bj, 0 )*100/min(lbi,lbj);
			if (table) printf("%4d",s);
			if (s >= minscore)
			{	Pairs	*p = pair+npairs;
				(*scores)[npairs] = s;
				sum_scores += s;
				p->block_a = i;
				p->block_b = j;
				npairs++;
				if (pep_score) {
					pscores++;
					ms = ms + s;
					mp = mp + pep_score;
					sp = sp + (float)s/(float)pep_score;
				} else 	ps = ps + s;
			}
			if ( npairs >= pair_alloc ) {
				pair_alloc += PAIRALLOC;
				pair = (Pairs*) realloc(
					pair,sizeof(Pairs)*pair_alloc);
					TEST(pair)
				*scores = (int*) realloc(
					*scores,sizeof(int)*pair_alloc);
					TEST(*scores)
			}
		}
	}
	if ( !npairs ) {
		printf("\nNo pairs over cutoff\n\n");
		return 0;
	}
	sum_scores = sum_scores/npairs;
	printf("\naverage score over cutoff = %d\n",sum_scores);
	if (!pep_len) return npairs;
	if (pscores) {
		ms = ms/pscores; mp = mp/pscores;
		printf("mean score = %f, mean pep_count = %f\n", ms,mp);
		printf("score/peptides = %f", sp/pscores);
		printf("  (over %d pairs)\n", pscores);
		if (pscores<npairs) {
			i = npairs-pscores;
			printf("mean non-pept score = %f",ps/(float)i);
			printf("  (over %d pairs)\n", i);
		}
	}
	NL
	return npairs;
} 

cluster_pairs (place)
	int	*place;
{
	int	i, j;
	nlists = 0;
	for ( i = 0; i < npairs; i++)						/* LOOP OVER ORDERED PAIRS */
	{	int	a = (pair+place[i])->block_a,	
			b = (pair+place[i])->block_b,
			used_a = block[a]->grouped,
			used_b = block[b]->grouped;
		for ( j = 1; j <= nlists; j++)					/* LOOP OVER LISTS */
		{	Lists	*lis = list+j;
			if ( lis->closed || lis->dead ) continue;
			if (extend(j,a,b)) continue; 
			if (extend(j,b,a)) continue;
		}
		if ( used_a || used_b ) continue;
		nlists++;
		if ( nlists >= list_alloc ) 
		{	int	old_alloc = list_alloc;
			list_alloc += LISTALLOC;
			list = (Lists*) realloc(list,sizeof(Lists)*list_alloc);
			TEST(list)
		}
		{	Lists	*new_list = list+nlists;			/* CREATE A NEW LIST */
			Halfs	*l = new_list->left,
				*r = new_list->right;
			new_list->closed = FALSE;
			new_list->dead = FALSE;
			l->last = r->last = 1;
			l->member[1] = a;
			r->member[1] = b;
			if (list_limit == 2) new_list->closed = TRUE;
			block[a]->grouped = block[b]->grouped = nlists;		/* SET TO INDICATE BLOCK HAS BEEN ALLOCATED A GROUP */
		}
	}
	printf("total list_alloc = %d\n", list_alloc);
}

extend (current,blka,blkb)
	int	current, blka, blkb;
{
	Blocks	*block_a = block[blka],
		*block_b = block[blkb];
	int	used = block_b->grouped,
		new_side_len, new_list_len, i;
	Lists	*list_a = list+current;
	Halfs	*al = list_a->left,
		*ar = list_a->right,
		*side = 0;

	if (used<0) return 1;
	for (i=1; i<=al->last; i++) {
	int	al_edge = al->member[i];
		if ( blka == al_edge ) side = al;
	}
	for (i=1; i<=ar->last; i++) {
	int	ar_edge = ar->member[i];
		if ( blka == ar_edge ) side = ar;
	}
	if ( !side ) return 0;
	if ( !used ) {
		side->last++;
		if ( side->last >= MAXLISLEN ) return 1;
		new_list_len = al->last + ar->last;
		if (new_list_len == list_limit ) {
			list_a->closed = TRUE;
			block_b->grouped = -1;					/* BOTH list AND block closed TO AVOID EITHER JOINING LATER LISTS */
		}
		side->member[side->last] = blkb;
		block_b->grouped = current;
		return 1;
	} else {
	Lists	*list_b = list+used;
	Halfs	*bl = list_b->left,
		*br = list_b->right;
	int	bl_edge = bl->member[bl->last],
		br_edge = br->member[br->last];
		if ( list_b->dead || list_b->closed ) return 1;			/* THE joinING LIST IS closed */
		if ( used == current ) {
			if ( circles ) list_a->closed = TRUE;			/* CLOSE OFF THE LIST */
			return 1;
		}
		if ( blkb != bl_edge && blkb != br_edge ) return 1;		/* NOT AN EDGE IN THE OTHER LIST */ 
		new_side_len = bl->last + br->last + side->last;
		if ( new_side_len >= MAXLISLEN ) return 1;			/* NO ROOM TO MERGE LISTS */
		new_list_len = bl->last + br->last + al->last + ar->last;
		if ( new_list_len > list_limit ) return 1;
		if ( new_list_len == list_limit ) {
				list_a->closed = TRUE;
			list_b->closed = TRUE;					/* BOTH LISTS closed TO AVOID EITHER JOINING LATER LISTS */
		}
		list_b->dead = TRUE;
		if ( blkb == bl_edge ) {					/* COPY THE LIST ONTO THE LEFT */
			for ( i = bl->last; i > 0; i-- )
				side->member[++(side->last)] = bl->member[i];
			for ( i = 1; i <= br->last; i++)
				side->member[++(side->last)] = br->member[i];
		}
		if ( blkb == br_edge ) {					/* COPY THE LIST ONTO THE RIGHT */
			for ( i = br->last; i > 0; i-- ) 
				side->member[++(side->last)] = br->member[i];
			for ( i = 1; i <= bl->last; i++)
				side->member[++(side->last)] = bl->member[i];  
		}
		block[side->member[side->last]]->grouped = current; 		/* RESET EDGE POINTERS */
		return 1;
	}
}

/*	ALIGN SELECTED PAIRS AND SET UP POINTERS BETWEEN THEM	*/

pack_lists ()
{
	int	i, j;
	ngroups = 0;
	for ( i = 1; i <= nlists; i++ )						/* LOOK FOR A FREE END */
	{	Lists	*lis = list+i;
		Halfs	*l = lis->left,
			*r = lis->right;
		Groups	*g = group+ngroups;
		if ( lis->dead ) continue;
		ngroups++;
		g->members = 0;
		for ( j = l->last; j > 0; j-- ) 
		{	int	ml = l->member[j];	
			g->blk[(g->members)++] = block[ml];
			block[ml]->grouped = ngroups;
		}
		for ( j = 1; j <= r->last; j++)
		{	int	mr = r->member[j];
			g->blk[(g->members)++] = block[mr];
			block[mr]->grouped = ngroups;
		}
	}
}

align_pairs (aln,cluster)
	Alns	*aln;
	Groups	*cluster;
{
	int	nalns = cluster->members - 1,
		i, j, k;
	for ( k = 0; k < nalns; k++ )
	{ 	Alns	*al = aln+k,
			*ak = aln+k+1;
		int	*last = al->position_b,
			*next = ak->position_a;
		int	ncomp, minl;
		ak->block_a = cluster->blk[k];
		ak->block_b = cluster->blk[k+1];
		ak->score = score( ak->block_a, ak->block_b, ak );
		minl = min(ak->block_a->length,ak->block_b->length);
		ak->block_b->score_to_last = ak->score*100/minl;
		for ( i = 0; i <= ak->length+1; i++ )
			ak->last_link[i] = ak->next_link[i] = 0;		/* INITIALISE TOP ROW OF POINTERS */
		if ( !k ) continue;
		else for ( i = 0; i <= al->length+1; i++ ) 
			al->next_link[i] = 0;					/* INITIALISE BOTTOM ROW OF POINTERS */
		i = j = 1;							/* SET POINTERS BETWEEN ALIGNMENTS */
		while ( (i <= al->length) &&
			(j <= ak->length) ) {					/* LOOP TO FIND EQUIVALENT POSITIONS */
			if ( last[i] > next[j] ) { 
				ak->last_link[j++] = 0; 
				continue; }
			if ( last[i] < next[j] ) {
				al->next_link[i++] = 0; 
				continue; }
			if ( last[i] != next[j] ) continue;
			if ( last[i] && next[j] ) {
				al->next_link[i] = j; 
				ak->last_link[j] = i;
			} else  al->next_link[i] = ak->last_link[j] = 0;
			i++; j++;
		}
	}
}

/*	FILLS AN ARRAY WITH THE ALIGNMENTS	*/

expand_seqs (aln, cluster)
	Alns	*aln;
	Groups	*cluster;
{
	int	max_lines = 0, max_leng = MIN,
		grp_memb = cluster->members,
		new_len, 
		i, j, k, l;
	for ( i = 0; i < grp_memb; i++ )
	{	Blocks	*b = cluster->blk[i];
		int	blk_memb = b->members,
			blk_leng = b->length;
		max_leng = max(max_leng, blk_leng);
		max_lines += blk_memb;
	}
	for (i=0; i <= max_lines; i++) colmax[i]=0;
	aln_len = max_leng*10*max_lines+1;
	alignment = (char*)malloc(sizeof(char)*aln_len); TEST(alignment)
	for(i=0;i<aln_len;i++) alignment[i] = ' ';
	mincol = MAX;
	maxcol = MIN;
	enter( aln+1, max_leng/2, 0, max_leng*10, 1 );
	max_leng += max_leng;
	l = 0;
	for ( i = 0; i < grp_memb; i++ ) 
	{	Blocks	*b = cluster->blk[i];
		int	new_len = maxcol-mincol+1,
			size = sizeof(char)*(new_len+1),
			relocate = 0;
		if ( new_len > b->length ) {
			relocate = 1;
			b->length = new_len; 
		}
		for ( j = 0; j < b->members; j++)
		{	Seqs	*s = b->seq[j];
			int	n = 0;
			if ( relocate ) { 
				s->residue = (char*) realloc(s->residue,size);
				TEST(s->residue)
			}
			for ( k = mincol; k <= maxcol; k++)
				s->residue[++n] = alignment[k+l];
			l += max_leng;
		}
	}
	free(alignment);
}


enter ( aln, col, row, leng, pos)
	Alns	*aln;
	int	col, row, leng, pos;
{
	Blocks	*blk_a = aln->block_a,
		*blk_b = aln->block_b;
	int	memb_a = blk_a->members,
		memb_b = blk_b->members,
		a = aln->position_a[pos],
		b = aln->position_b[pos],
		down = aln->next_link[pos];
	int	cell;

	col = bump_check(row, col, memb_a, a, b);

	if (shift(aln,pos,LEFT,-1))
			enter(aln,   col-1, row,        leng,	pos-1);
	if (down) col = enter(aln+1, col,   row+memb_a, leng,	down );
	if (shift(aln,pos,RIGHT,aln->length))
			enter(aln,   col+1, row,        leng,	pos+1);

	cell = col + row*leng;
	if (a) 
	{	int	i;
		for ( i = 0; i < memb_a; i++)
		{	Seqs	*s = blk_a->seq[i];		if(alignment[cell] != ' ' && alignment[cell] != s->residue[a])
									printf("*** cell overwritten ***\n");
			alignment[cell]= s->residue[a];		if(cell >= aln_len || col < 0 || row < 0 ) 
									printf("*** alignment overflow *** %d %d\n",row,col);
			if ( col > colmax[row] ) colmax[row] = col;
			cell += leng;
			row++;
		}
	} else { 
		if ( col > colmax[row] ) colmax[row] = col;
		cell += memb_a*leng;
		row  += memb_a;
	}
	if (b) 
	{	int	i;
		for ( i = 0; i < memb_b; i++)
		{	Seqs	*s = blk_b->seq[i];		if(alignment[cell] != ' ' && alignment[cell] != s->residue[b])
									printf("*** cell overwritten ***\n");
			alignment[cell] = s->residue[b];	if(cell >= aln_len || col < 0 || row < 0 ) 
									printf("*** alignment overflow *** %d %d\n",row,col);
			if ( col > colmax[row] ) colmax[row] = col;
			cell += leng;
			row++;
		}
	} else if ( col > colmax[row] ) colmax[row] = col;
	maxcol = max(maxcol,col);				if ( maxcol > leng ) printf("*** maxcol > length *** (%d)\n",col);
	mincol = min(mincol,col);				if ( mincol < 0 ) printf("*** mincol < 0 *** (%d)\n",col);
	return col;
}

shift ( aln, pos, direction, edge)
	Alns	*aln;
	int	pos, direction, edge;
{
	int	new_pos = pos+direction;

	if (aln->last_link[new_pos]) return 0;
	if (new_pos*direction > edge) return 0;
	aln->last_link[pos] = -1;
	aln->last_link[new_pos] = -1;
	return 1;
}

bump_check (row, col, memb, a, b)
	int	col, memb, a, b;
{
	if (a) { 
		if (col <= colmax[row]) {
			col = colmax[row]+1; 
			return col;
		}
	}
	row += memb;
	if (b) { 
		if (col <= colmax[row]) {
			col = colmax[row]+1;
			return col;
		}
	}
	return col;
}

update_blocks ()
{
	int	i, j, k, l, len;
	char	temp[TREE+1];
	for ( i = 0; i < nblocks; i++ ) 
	{ 	Blocks	*b = block[i];
		if (!b->members) continue;
		if (b->grouped) continue;
		if (b->members == 1)
		{ char	*tree = b->seq[0]->tree;
			len = strlen(tree);
			strcpy(temp,tree);
			strncpy(tree+2,temp,len);
			tree[len+2] = 0;
			tree[0] = '-';
			tree[1] = '-';
			continue;
		}
		for ( j = 0; j < b->members; j++)
		{ char	*tree = b->seq[j]->tree;
			len = strlen(tree);
			strcpy(temp,tree);
			strncpy(tree+2,temp,len);
			tree[len+2] = 0;
			tree[0] = '.';
			tree[1] = '.';
		}
	}
	/* LOOP OVER THE GROUPS */
	for ( i = 0; i < ngroups; i++ )
	{	Groups	*g = group+i;
		Blocks	*core = g->blk[0];
		int	g_memb = g->members,
			c_memb = core->members,
			link, top, max_blk = MIN;
		core->grouped = 1;
		top = 999;
		link = c_memb - 1;
		for ( k = 0; k < c_memb; k++ )
		{ char	*tree = core->seq[k]->tree,
			*code = core->seq[k]->code;
			printf("     %s %s\n", tree,code);
			if (tree[0]=='b') { link = k; top = 0; }
			for (l=0; core->seq[link]->tree[l]=='.'; l++) {
				if (l==TREE) {
					printf("*NB* tree too big\n");
					break;
				}
				if (tree[l+1]=='b' && l<top) {
					link = k;
					top = l+1;
				}
			}
			len = strlen(tree);
			strcpy(temp,tree);
			strncpy(tree+2,temp,len);
			tree[len+2] = 0;
			tree[0] = '.';
			tree[1] = '.';
		}
		core->seq[link]->tree[0] = 'p';
		for (k=1; core->seq[link]->tree[k]=='.'; k++) {
			if (k==TREE) {
				printf("*NB* tree too big\n");
				break;
			}
			core->seq[link]->tree[k] = '-';
		}
		/* LOOP OVER CORE SEQUENCES */
		for ( j = 1; j < g_memb; j++ )
		{	Blocks	*b = g->blk[j];
			int	b_memb = b->members,
				new_size = 
				((core->members)+b_memb)*(sizeof(Seqs*));
			core->seq = (Seqs**) realloc(core->seq,new_size);
				TEST(core->seq)
			printf("%d\n",b->score_to_last);
			link = 0;
			top = 999;
			/* LOOP OVER BLOCK SEQUENCES */
			for ( k = 0; k < b_memb; k++ )
			{	int	m = (core->members)++;
				char	*tree = b->seq[k]->tree,
				    	*code = b->seq[k]->code;
				printf("     %s %s\n", tree,code);
				if (tree[0]=='p') { link = k; top = 0; }
				for (l=0; b->seq[link]->tree[l]=='.'; l++) {
					if (l==TREE) {
						printf("*NB* tree too big\n");
						break;
					}
					if (tree[l+1]=='p' && l<top) {
						link = k;
						top = l+1;
					}
				}
				len = strlen(tree);
				strcpy(temp,tree);
				strncpy(tree+2,temp,len);
				tree[len+2] = 0;
				tree[0] = '.';
				tree[1] = '.';
				core->seq[m] = b->seq[k];
				if (b->key) core->key = 1;
			}
			b->seq[link]->tree[0] = 'b';
			for (k=1; b->seq[link]->tree[k]=='.'; k++) {
				if (k==TREE) {
					printf("*NB* tree too big\n");
					break;
				}
				b->seq[link]->tree[k] = '-';
			}
			b->members = 0;
		}
		NL NL
		link = 0;
		for ( k = core->members-1; k >= 0; k-- )
		{ char	*tree = core->seq[k]->tree;
			if (tree[0]=='p') break;
			if (tree[0]=='b') { link = 1;}
			else { if (link) tree[0] = '|';}
		}
	}
}	
				
/*	CALCULATES THE MAXIMUM SCORE MATRIX	*/

score ( blka, blkb, aln) 
	Blocks	*blka, *blkb;
	Alns	*aln;
{
	Seqs	*seqa = blka->seq[0],
		*seqb = blkb->seq[0];
	int	*posa, *posb,
		lena = blka->length,
		lenb = blkb->length,
		block_aln = blka->members + blkb->members - 2,
		trace_back = (int)aln,
		ldif = lena-lenb,
		windo = window-1, winl, winr,
		wdisp, wstart, wstops,
		now = 0, last = 1,
		enter = 0, edge_set = 0,
		pati, patj,
		mati, matj, 
		i, j;
	int	maxrows[MAXWINLEN],
		toprows[MAXWINLEN],
		maxscore = MIN, 
		mat[2][MAXWINLEN],
		mul[2][MAXWINLEN];
	char	na, nb;
	for ( i = 0; i < window; i++ ) {
		mul[0][i] = 0;
		mul[1][i] = 0;
		mat[0][i] = WEE;
		mat[1][i] = WEE;
		pat[i][lenb] = 0;
		maxrows[i] = WEE;
		toprows[i] = lenb;
	}
	winl = max(0,min(windo,(windo+ldif)/2));
	winr = max(0,min(windo,(windo-ldif)/2));
	if (winl+winr < windo ) if ( lena < lenb ) winr++; else winl++;  
	for ( j = lenb; j > 0; j-- )
	{	int	wmid = lena-(lenb-j),
			maxcol = WEE,
			topcol = 0;
		if( --enter < 0 ) enter = window-1;
		wdisp = winl-wmid;
		wstart = max(   1, wmid-winl);
		wstops = min(lena, wmid+winr);
		if ( trace_back && wstops == lena ) pat[lena-wstart][j] = 0;
		nb = UPPER(seqb->residue[j])-'A';
		for ( i = wstops; i >= wstart; i-- ) 
		{	int	wi = i+wdisp,
				maxrow, toprow,
				diag, row, col,
				rat = enter+wi;
			if ( block_aln)
				mat[now][wi] = dist( blka, blkb, i, j );
			else {
				na = UPPER(seqa->residue[i])-'A';
				mat[now][wi] = matrix[na][nb];
			}
			if ( i == lena || j == lenb ) continue;
			if ( rat < 0 ) rat += window;
			if ( rat >= window ) rat -= window;
			if ( i == wstart && !edge_set ) {
				maxrows[rat] = WEE;
				toprows[rat] = j;
				if( i == 1 ) edge_set = TRUE;
			}
			maxrow = maxrows[rat];
			toprow = toprows[rat];
			diag = mat[last][wi];
			row = maxrow-gap_pen; 
			col = maxcol-gap_pen; 
			if (mul_wt) {
				mul[now][wi] = ((1+mat[now][wi])*(1+mul[last][wi]))/damp;
				diag += (mul[now][wi]*mul_wt)/100; 
			}
			if ( diag >= row && diag >= col ) {
				mat[now][wi] += diag;
/*
Pi(mul[last][wi]) Pi(mul[now][wi]) Pi(mat[last][wi]) Pi(mat[now][wi]) NL
*/
				if (trace_back) pat[wi][j] = 1; 
			} else {
				if ( col > row ) {
					mat[now][wi] += col;
					mul[now][wi] = 0;
					if (trace_back) 
						pat[wi][j] = topcol-i+1;
				} else {
					mat[now][wi] += row;
					mul[now][wi] = 0;
					if (trace_back)
						pat[wi][j] = -(toprow-j+1);
				}
			}
			if (diag > maxcol) {
				maxcol = diag;
				topcol = i;
			}
			if (diag > maxrow) {
				maxrows[rat] = diag;
				toprows[rat] = j;
			}
			if ( mat[now][wi] >= maxscore ) {
				maxscore = mat[now][wi];
				if ( trace_back ) {
					pati = wi; patj = j;
					mati = matj = 1;
					if ( i == 1 ) matj = j;
					if ( j == 1 ) mati = i;
				}
			}
		}
/*
if (trace_back && output) {
printf("\n");
printf("*%c",nb+'A');
for (i=1; i<wstart; i++)  printf("    ");
for (i=wstart; i<=wstops; i++) 
printf("%4d",max(-999,min(999,mul[now][i+wdisp])));
printf("\n");
printf("+%c",nb+'A');
for (i=1; i<wstart; i++) printf("    "); 
for (i=wstart; i<=wstops; i++) 
printf("%4d",max(-999,min(999,mat[now][i+wdisp])));
printf("\n");
}
*/		i = now; now = last; last = i;
	}
	if ( !trace_back ) return maxscore;
	posa = aln->position_a;
	posb = aln->position_b;
	aln->length = 0;
	trace( posa, posb, mati, matj, pati, patj, 0, 0, &aln->length );
	posa += aln->length; posb += aln->length;
	if ( *posa == lena ) {
		for ( i = *(posb)+1; i <= lenb; i++) { 
			*(++posa) = 0; *(++posb) = i; 
			(aln->length)++; 
		}
		if ( aln->length > MAXALNLEN ) 
			printf("*NB* aln_length exceeded\n");
		return maxscore;
	}
	if ( *posb == lenb )	
		for ( i = *(posa)+1; i <= lena; i++) { 
			*(++posb) = 0; *(++posa) = i; 
			(aln->length)++; 
		}
	if ( aln->length > MAXALNLEN ) 
		printf("*NB* aln_length exceeded\n");
	return maxscore;
}

trace ( posa, posb, mati, matj, pati, patj, lasti, lastj, n )
	int	*posa, *posb,
		mati,  matj, 
		pati,  patj,
		lasti, lastj, 
		*n;
{
	int	pij = pat[pati][patj],
		i, j;

	for ( i=lasti+1; i < mati; i++) 
		{ *(++posa) = i; *(++posb) = 0; (*n)++; } 
	for ( j=lastj+1; j < matj; j++) 
		{ *(++posa) = 0; *(++posb) = j; (*n)++; }
	*(++posa) = mati; *(++posb) = matj; (*n)++;

	if ( !pij ) return;

	if ( pij == 1 ) 
		trace ( posa, posb, mati+1,   matj+1,   pati,       patj+1,
			mati, matj, n );
	if ( pij <  1 ) 
		trace ( posa, posb, mati+1,   matj-pij, pati+pij+1, patj-pij,
			mati, matj, n );
	if ( pij >  1 ) 
		trace ( posa, posb, mati+pij, matj+1,   pati+pij-1, patj+1,
			mati, matj, n );
}

dist (a, b, pos_a,pos_b)
	Blocks	*a, *b;
	int	pos_a, pos_b;
{
	register	i, j, sum = 0;

	for ( i = 0; i < a->members; i++ ) 
	{	char	res_a = UPPER((a->seq[i])->residue[pos_a]);
		if ( res_a < 'A' || res_a > 'Z' ) continue;
		for ( j = 0; j < b->members; j++ )
		{	char	res_b = UPPER((b->seq[j])->residue[pos_b]);
			if ( res_b < 'A' || res_b > 'Z' ) continue;
			sum += matrix[res_a-'A'][res_b-'A'];
		}
	}
	sum = sum/(a->members*b->members);
	return sum;
}
