#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<memory.h>
#include<time.h>

#define LEN 200
#define NUM 30
int flag = 0;
char ran(char a, int k, int b)
{
    srand(time(NULL)+k*b);
    int i = rand()%100;
    if (flag == 0)
    {
        if(i<20)
            return 'A' + i;
        if(i == 20)
            return '_';
        return a;
    }
    else
    {
        if(i<20)
            return '_';
        if(i<40)
            return 'A'+i-20;
        return a;
    }
}
int main()
{
    FILE *f;
    f = fopen("seq.txt", "w");
    char seq[LEN], data[LEN];
    int i, j, key;
    for(i = 0; i < LEN; i++)
    {
        srand((unsigned)time(NULL)+i);
        key = rand()%4;
        //printf("%d", key);
        switch(key)
        {
        case 0:
            seq[i] = 'A';
            break;
        case 1:
            seq[i] = 'C';
            break;
        case 2:
            seq[i] = 'T';
            break;
        case 3:
            seq[i] = 'G';
            break;
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
            if (data[i] == '_')
            {
                flag = 1;
            }
            else
            {
                flag = 0;
            }
            fprintf(f, "%c", data[i]);
        }
        fprintf(f, "\n");
        //fprintf(f,"\n", seq[i]);
    }
    fclose(f);



    return 1 ;
}
