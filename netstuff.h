#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#ifndef RAN2_H

#include "ran2.h"

#endif


int compare (const void * a, const void * b) {
      return ( *(int*)a - *(int*)b );
}

typedef struct Lista list;
struct Lista {
    int allocated;
    int leng;
    int *list;
};

list initList(int size) {
    list lista;
    lista.allocated = size;
    lista.leng = 0;
    lista.list = (int*) malloc(size* sizeof(int));
    if (!lista.list) printf("MEMORY ERROR!\n");
    return lista;
}

void markVisited(bool *visited, list *Nvec, int index) {
    int i,neig;
    for (i=0;i<Nvec[index].leng;i++){
        neig = Nvec[index].list[i];
        if ( ! visited[neig] ) {
            visited[neig] = true;
            markVisited(visited, Nvec, neig);
        }
    }
}

int connectedComponenets(list *Nvec, int size) {
    bool visited[size];
    int i;
    int compCount = 0;
    for (i=0;i<size;i++)
        visited[i] = false;

    for (i=0;i<size;i++) {
        if ( ! visited[i] ) {
            compCount++;
            visited[i] = true;
            markVisited(visited, Nvec, i);
        }
    }
    return compCount;
}

void resizeList(list *lista, int size) {
    //printf("List resized\n");
    lista->allocated += size;
    int *new_list = (int*) realloc(lista->list, lista->allocated* sizeof(int));
    if (!new_list)
    {
        printf("NO MEMORYYYY\n\n\n");
    }
    lista->list = new_list;
}

//this implements a binary search to find the index to which the new value should be appended
int findApIndex(list *lista, int *value) {
    int low,high,mid;
    low = 0;
    high = lista->leng-1;
    mid = (high + low)/2;
    while (high > low)
    {
        if (lista->list[mid] > *value)
            high = mid;
        else
            low = mid+1;
        mid = (high + low)/2;
    }
    return mid;
}

int isInList(list *lista, int value) {
    char *pItem;
    pItem = (char*) bsearch (&value, lista->list, lista->leng, sizeof(int), compare);
    if (pItem == NULL) return 0;
    else return 1;
}

void appendList(list *lista, int value) {
    //checks if the allocated space is enough
    //printf("\nleng = %i allocated = %i\n",lista->leng, lista->allocated);
    if (lista->leng >= lista->allocated - 2) resizeList(lista, 10);
    //finds where the new value should be appended
    if (lista->leng > 0 && lista->list[lista->leng-1] > value)
    {
        int index = findApIndex(lista, &value);
        //shifts all the bigger values
        //printf("index = %i old %i new %i\n",index,lista->list[index],value);
        for (int i=lista->leng;i>index;i--) {
            lista->list[i] = lista->list[i-1];
        }
        lista->list[index] = value;
    }
    else
        lista->list[lista->leng] = value;

    lista->leng++;
    //qsort (lista->list, lista->leng, sizeof(int), compare);
}

void removeConnection(list *Nvec, int *node1, int *node2) {
    int j;
    int index1 = findApIndex(&Nvec[*node1], node2);
    int index2 = findApIndex(&Nvec[*node2], node1);
    for (j=index2;j<Nvec[*node1].leng;j++)
        Nvec[*node1].list[j] = Nvec[*node1].list[j+1];
    Nvec[*node1].leng--;
    for (j=index1;j<Nvec[*node2].leng;j++)
        Nvec[*node2].list[j] = Nvec[*node2].list[j+1];
    Nvec[*node2].leng--;
}

void removeNode(list *lista, int value) {
    //checks if the allocated space is enough
    //printf("\nleng = %i allocated = %i\n",lista->leng, lista->allocated);
    if (lista->leng >= lista->allocated - 1) resizeList(lista, 10);
    //finds where the new value should be appended
    if (lista->leng > 0 && lista->list[lista->leng-1] > value)
    {
        int index = findApIndex(lista, &value);
        //shifts all the bigger values
        //printf("index = %i old %i new %i\n",index,lista->list[index],value);
        for (int i=lista->leng;i>index;i--) {
            lista->list[i] = lista->list[i-1];
        }
        lista->list[index] = value;
    }
    else
        lista->list[lista->leng] = value;

    lista->leng++;
    //qsort (lista->list, lista->leng, sizeof(int), compare);
}

list* randomGraph(int N, float k, struct ran_state *state) {
    int i,j,node1,node2,count,nCount;
    list* Nvec;
    Nvec = (list*) malloc(sizeof(list)* N);
    int size = ((int) k)+1;
    //printf("Beginning %i\n",size);
    for (i=0;i<N;i++) Nvec[i] = initList(size);

    //printf("Beginning 0.5\n");
    int auxN[N];
    for (i=0;i<N;i++) auxN[i] = 0;

    nCount = 0;
    //printf("Beginning 2\n");
    int max = (int) N*k/2.;
    while (nCount< max)
    {
        node1 = ran2int(0,N,state);
        //if (node1>=N || node1<0) printf("node1 too big\n");
        //printf("node1 = %i\n",node1);
        //reinitializes the options
        count = 0;
        for (i=0;i<N;i++)
        {
            //check if connection is availiable
            if (i != node1 && (! isInList(&Nvec[node1], i)))
            {
                //if (count>=N || count<0) printf("count too big\n");
                //printf("before auxN.");
                auxN[count] = i;
                //printf("after auxN.");
                count++;
            }
        }
        if (count > 0)
        {
            int index2 = ran2int(0,count,state);
            //if (index2>=N || index2<0) printf("index2 too big\n");
            //printf("before index2 %i / %i",index2,count);
            node2 = auxN[index2];
            //printf("node1=%i node2=%i \n",node1,node2);
            appendList(&Nvec[node1],node2);
            //printf("between appendList.");
            appendList(&Nvec[node2],node1);
            nCount++;
        }
        //printf("Node loop ended. %i of %i\n",nCount, max);
    }
    return Nvec;
}

list* randomModularGraph(int N, int n1, float k, float mu, struct ran_state *state) {
    int i,j,node1,node2,count,nCount;
    list* Nvec;
    Nvec = (list*) malloc(sizeof(list)* N);
    for (i=0;i<N;i++) Nvec[i] = initList(((int) k)+1);
    int auxN[N];

    nCount = 0;
    while (nCount<(int) N*(1.-mu)*k/2.)
    {
        node1 = ran2int(0,N,state);
        if (node1 < n1)
        {
            //reinitializes the options
            for (i=0;i<n1;i++) auxN[i] = 0;
            count = 0;
            for (i=0;i<n1;i++)
            {
                //check if connection is availiable
                if (i != node1 && (! isInList(&Nvec[node1], i)))
                {
                    auxN[count] = i;
                    count++;
                }
            }
            if (count > 0)
            {
                node2 = auxN[ran2int(0,count,state)];
                appendList(&Nvec[node1],node2);
                appendList(&Nvec[node2],node1);
                nCount++;
            }
        }
        else
        {
            //reinitializes the options
            for (i=n1;i<N;i++) auxN[i] = 0;
            count = 0;
            for (i=n1;i<N;i++)
            {
                //check if connection is availiable
                if (i != node1 && ! isInList(&Nvec[node1], i))
                {
                    auxN[count] = i;
                    count++;
                }
            }
            if (count > 0)
            {
                node2 = auxN[ran2int(0,count,state)];
                appendList(&Nvec[node1],node2);
                appendList(&Nvec[node2],node1);
                nCount++;
            }
        }
    }
    nCount = 0;
    while (nCount<(int) N*mu*k/2.)
    {
        node1 = ran2int(0,N,state);
        if (node1 > n1)
        {
            //reinitializes the options
            for (i=0;i<n1;i++) auxN[i] = 0;
            count = 0;
            for (i=0;i<n1;i++)
            {
                if (! isInList(&Nvec[node1], i))
                {
                    auxN[count] = i;
                    count++;
                }
            }
            if (count > 0)
            {
                node2 = auxN[ran2int(0,count,state)];
                appendList(&Nvec[node1],node2);
                appendList(&Nvec[node2],node1);
                nCount++;
            }
        }
        else
        {
            //reinitializes the options
            for (i=n1;i<N;i++) auxN[i] = 0;
            count = 0;
            for (i=n1;i<N;i++)
            {
                if (! isInList(&Nvec[node1], i))
                {
                    auxN[count] = i;
                    count++;
                }
            }
            if (count > 0)
            {
                node2 = auxN[ran2int(0,count,state)];
                appendList(&Nvec[node1],node2);
                appendList(&Nvec[node2],node1);
                nCount++;
            }
        }
    }
    return Nvec;
}
