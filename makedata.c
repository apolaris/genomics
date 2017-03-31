#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<memory.h>
#include<time.h>

#define LEN 150
#define NUM 100
char ran(char a, int k, int b)
{
    srand(time(NULL)+k+10000*b);
    int i = rand()%11;
    if(i<2)
        return '_';
    if(i<3)
        return 'A';
    if(i<4)
        return 'C';
    if(i<5)
        return 'T';
    if(i<6)
        return 'G';
    return a;
}
int main()
{
    FILE *f;
    f = fopen("seq.txt" , "w");
    char seq[LEN], data[LEN];
    int i, j, key;
    for(i = 0;i < LEN; i++)
    {
        srand((unsigned)time(NULL)+i);
        key = rand()%4;
        //printf("%d", key);
        switch(key){
            case 0: seq[i] = 'A'; break;
            case 1: seq[i] = 'C'; break;
            case 2: seq[i] = 'T'; break;
            case 3: seq[i] = 'G'; break;
        }
        fprintf(f,"%c", seq[i]);
    }
    fclose(f);

    printf("%s \n", seq);

    f = fopen("data.txt", "w");
    for(j = 0; j < NUM ; j++)
    {
        for(i = 0; i < LEN ; i++)
        {
            data[i] = ran(seq[i], i, j);
            fprintf(f, "%c", data[i]);
        }
        fprintf(f, "\n");
        //fprintf(f,"\n", seq[i]);
    }
    fclose(f);



    return 1 ;
}
