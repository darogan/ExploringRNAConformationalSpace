#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

main () {
FILE *file;
char	c;
int	i=0;
	file = fopen("test.txt","r");
	while(c=getc(file)) {
		printf("%d %d >%c<\n", i,c,c);
		if (c=='\n') i++;
	}
}
