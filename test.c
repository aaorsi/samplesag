#include <stdlib.h>
#include <stdio.h>

int main()
{
	int i, j;
	int k;
	const char *c[5];

	c[0] = "a";
	c[1] = "v";
	c[2] = "r";

	i = 0;
	j = 1;
	k = i++ + j;

	printf("i %d j %d k %d c0 %s c1 %s c2 %s \n",i,j,k, c[0], c[1], c[2]);

}
	
