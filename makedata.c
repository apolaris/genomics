#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
int LEN;
int NUM;
typedef struct Tree{
        struct Tree *left;
        struct Tree *right;
        int seed;
        char *seq;
        int depth;
}tree;
FILE *fp;
int iii= 0, jjj = 0;
void addTree(tree *t, int se, int depth)
{
    if(se > t->seed)
    {
        if(t->right != NULL)
        {
            addTree(t->right, se, t->depth);
        }
        else
        {
            tree *node = malloc(sizeof(tree));
            node->left = node -> right = NULL;
            node->depth = depth+1;
            node->seed = se;
            t->right = node;
            return;
        }
    }
    else
    {
        if(t->left != NULL)
        {
            addTree(t->left, se, t->depth);
        }
        else
        {
            tree *node = malloc(sizeof(tree));
            node->left = node -> right = NULL;
            node->depth = depth+1;
            node->seed = se;
            t->left = node;
            return;
        }
    }
}
void freetree(tree* t)
{
    if(t->right != NULL)
    {
        freetree(t->right);
        free(t->right);
    }
    if(t->left != NULL)
    {
        freetree(t->left);
        free(t->left);
    }

    if(t->seq != NULL)
        free(t->seq);
}
void printtree(tree* t)
{
    if(t->left)
        printtree(t->left);
    printf("%d %d\n", t->seed, t->depth);
    if(t->right)
        printtree(t->right);
}
void makeseq(tree *t, char *seq)
{
    int i, key;
    t->seq = malloc(sizeof(char)* LEN);
    t->seq[0] = seq[0];
    for (i = 1; i< LEN; i++)
    {
        srand((unsigned)time(NULL)*i + t->seed);
        key = rand()%1000;
        if(key < 20)
            t->seq[i] = 'A' + key;
        else if(key < 30)
            t->seq[i] = '_';
        else if(key < 40 && t->seq[i-1] == '_')
            t->seq[i] = '_';
        else
            t->seq[i] = seq[i];
    }
    jjj++;
    if(t->right)
        makeseq(t->right, t->seq);
    if(t->left)
        makeseq(t->left, t->seq);

}
void printseq(tree *t)
{
    if(t->right == NULL && t->left == NULL)
    {
        fprintf(fp, "%s\n", t->seq);
        iii++;
        return;
    }
    if(t->left)
    {
        printseq(t->left);
    }
    if(t->right)
    {
        printseq(t->right);
    }

    return;
}
int main()
{
    printf("please input the length of the reference and the number of the individuals:\n");
    scanf("%d %d", &LEN, &NUM);
    fp = fopen("datatree2.txt", "w");
    int seed[NUM], i = NUM, key;
    printf("%d nodes\n", i);
    tree *T = malloc(sizeof(tree));
    T->seed = 5000;
    T->right = T->left = NULL;
    T->depth = 0;
    printf("begin ");
    for(i = 0; i<NUM;i++)
    {
        srand((unsigned)time(NULL)*i);
        seed[i] = rand()%10000;
        printf("%d ", seed[i]);
        addTree(T, seed[i], 0);
    }
    T->seq = malloc(sizeof(char) * LEN);
    for(i = 0;i < LEN; i++)
    {
        srand((unsigned)time(NULL)*i);
        key = rand()%23;
        T->seq[i] = 'A'+key;
    }
    printf("build OK \n");
    printtree(T);
    printf("printtree\n");
    makeseq(T->right, T->seq);
    makeseq(T->left, T->seq);
    printseq(T);
    freetree(T);
    free(T);
    free(fp);
    return 1;
}
