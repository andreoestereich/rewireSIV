#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include <omp.h>
#include "ran2.h"
#include "netstuff.h"

double contsIntOcurance (int *vector, int *size, int tobCounted) {
    double mean = 0.0;
    int i;
    for (i=0; i < *size; i++)
        mean += (double) (vector[i] == tobCounted);
    return mean/ *size;
}

double averageInt (int *vector, int *size) {
    double mean = 0.0;
    int i;
    for (i=0; i < *size; i++) mean += (double) vector[i];
    return mean/ *size;
}

double averageDouble (double *vector, int *size) {
    double mean = 0.0;
    int i;
    for (i=0; i < *size; i++) mean = mean + vector[i];
    mean = mean/ *size;
    return mean;
}

typedef struct inputs typeInput;
//structure to save all the inputs
struct inputs {
    int N;
    float k;
    int tmax;
    int saved_steps;
    int cnt;
    int samples;
    double D; // fraction of agents with opiniton in [0,1]
    double lambda; // infection rate
    double alpha; // recovery rate
    double phi; // vaccination rate
    double epsilon; //distance of separation for rewiring
    double r; // rate of rewiring
    double w; // stregth of the influence of infection on opinions
    int pHist;
    int netType;
};
typeInput input =
{10000,    //int N;
50.0,     //float k;
4000,     //int tmax;
2000,     //int saved_steps;
3,        //int cnt;
20,       //int samples;
0.4,      //double D; // fraction of agents with opiniton in [0,1]
1.0,      //double lambda; // infection rate
0.1,      //double alpha; // recovery rate
0.01,     //double phi; // vaccination rate
1.00,     //double epsilon; //distance of separation for rewiring
0.00,     //double r; // rate of rewiring
0.0,      //double w; // stregth of the influence of infection on opinions
0,        //int pHist;
0};       //int netType;

void printFile(char *text) {
    FILE * outputFile;
    char filename[200];

    //printf("Open the output file\n");
    sprintf(filename, "SIRrewire_l_%lf_w_%lf_r_%lf.dat",input.lambda,input.w,input.r);

    /* check if file exists if it does append to it if it does't generate a header */
    outputFile = fopen( filename, "r" );
    if (outputFile == NULL) {
        outputFile = fopen( filename, "w" );
        //making the file header
        fprintf(outputFile,"# coupled infection with knetic opinion\n");
        fprintf(outputFile,"# Order parameter time evolution\n#$N=%i \\, t_{max}=%i \\, t_{stored}=%i \\, k=%g \\, D=%g \\, r=%g \\, \\epsilon=%g \\, \\lambda=%g netType=%i$\n", input.N, input.tmax, input.saved_steps, input.k , input.D, input.r, input.epsilon, input.lambda, input.netType);
        fprintf(outputFile, "#mean_op, mean_I, mean_S, mean_R, mean_X, max_I, connected\n");
    }
    else
        fclose(outputFile);

    outputFile = fopen( filename, "a" );
    fprintf(outputFile,"%s\n",text);
    fclose(outputFile);
}

void inteEpid (double *opinions, list *Nvec, char *epid, struct ran_state *state) {
    int i, j, node1, node2;
    /* char auxState[input.N]; */
    /* for (i=0;i<input.N;i++) auxState[i] = 'C'; */
    /* printf("epid int set\n"); */
    for (i=0;i<input.N;i++)
    {
        node1 = ran2int(0,input.N,state);
    //printf("what?\n");
        if (Nvec[node1].leng > 0)
        {
            if (epid[node1] == 'S')
            {
        //printf("the\n"); (opinions[node1]+1.0)/2.0
                if (ran2(state) < ((opinions[node1]+1.0)/2.0))
                    epid[node1] = 'R';
                else
                {
                    //printf("node1 = %i node leng = %i\n",node1,Nvec[node1].leng);
                    node2 = Nvec[node1].list[ran2int(0,Nvec[node1].leng,state)];
                    //printf("node1 = %i node2 = %i\n",node1,node2);
                    if ( epid[node2] == 'I' )
                        epid[node1] = (ran2(state) < input.lambda ) ? 'I' : 'S';
                    else
                        epid[node1] = 'S';
                }
            }
            else
            {
                if (epid[node1] == 'I')
                    epid[node1] = (ran2(state) < input.alpha ) ? 'S' : 'I';
                else
                    epid[node1] = (ran2(state) < input.phi ) ? 'S' : 'R';;
            }
        }
        else
            epid[node1] = 'S';
    }
    /* for (i=0;i<input.N;i++) */
    /*     epid[i] = auxState[i]; */
}

void inteEpidCP (double *opinions, list *Nvec, char *epid, struct ran_state *state) {
    int i, j, node1, node2;
    double lamEff;
    /* char auxState[input.N]; */
    /* for (i=0;i<input.N;i++) auxState[i] = 'C'; */
    /* printf("epid int set\n"); */
    for (i=0;i<input.N;i++)
    {
        node1 = ran2int(0,input.N,state);
    //printf("what?\n");
        if (Nvec[node1].leng > 0)
        {
            if (epid[node1] == 'S')
            {
                if (ran2(state) < (opinions[node1]+1.0)/2.0)
                    epid[node1] = 'R';
                else
                {
                    node2 = Nvec[node1].list[ran2int(0,Nvec[node1].leng,state)];
                    lamEff = 0.0;
                    for (j=0; j< Nvec[node1].leng ; j++)
                        lamEff += (double) (epid[j] == 'I');
                    lamEff *= input.lambda/((double) Nvec[node1].leng);

                    epid[node1] = (ran2(state) < lamEff ) ? 'I' : 'S';
                }
            }
            else
            {
                if (epid[node1] == 'I')
                    epid[node1] = (ran2(state) < input.alpha ) ? 'S' : 'I';
                else
                    epid[node1] = (ran2(state) < input.phi ) ? 'S' : 'R';;
            }
        }
        else
            epid[node1] = 'S';
    }
    /* for (i=0;i<input.N;i++) */
    /*     epid[i] = auxState[i]; */
}

void inteEpidVar1a (double *opinions, list *Nvec, char *epid, struct ran_state *state) {
    int i, j, node1, node2;
    /* char auxState[input.N]; */
    /* for (i=0;i<input.N;i++) auxState[i] = 'C'; */
    /* printf("epid int set\n"); */
    for (i=0;i<input.N;i++)
    {
        node1 = ran2int(0,input.N,state);
        if (Nvec[node1].leng > 0)
        {
            if (epid[node1] == 'S')
            {
                if (ran2(state) < (opinions[node1]+1.0)/2.0)
                    epid[node1] = 'R';
                /* Unnecessary just to make clear what happens */
                /* else */
                /*     epid[node1] = 'S'; */
            }
            else
            {
                if (epid[node1] == 'I')
                {
                    if (ran2(state) < input.lambda/Nvec[node1].leng )
                        for(j=0;j<Nvec[node1].leng; j++)
                            epid[Nvec[node1].list[j]] = ( epid[Nvec[node1].list[j]] == 'R' ) ? 'R' : 'I' ;
                    else
                        epid[node1] = (ran2(state) < input.alpha ) ? 'S' : 'I';
                }
                else
                    epid[node1] = (ran2(state) < input.phi ) ? 'S' : 'R';;
            }
        }
        /* else */
        /*     epid[node1] = 'S'; */
    }
    /* for (i=0;i<input.N;i++) */
    /*     epid[i] = auxState[i]; */
}

//interaction of the opinions parallel
void inteOprefRw (double *opinions, char *epid, list *Nvec, struct ran_state *state) {
    int i, node1, node2, index1, index2, j, k, count, shuf1, shuf2;
    int auxN[input.N];
    for (j=0; j< input.N; j++)
        auxN[j] = j;
    double field;
    for (i=0; i< input.N; i++)
    {
        node1 = ran2int(0,input.N,state);
        //case in witch there is no neighbours has to be ignored
        if (Nvec[node1].leng > 0) {
                index2 = ran2int(0,Nvec[node1].leng,state);
                node2 = Nvec[node1].list[index2];
                //if (i == node2) //printf("%i -> %i\n",i,node2);
                //if (*N <= node2) //printf("%i -> %i\n",i,node2);
                if ((Nvec[node1].leng > 1) &&
                        (Nvec[node2].leng > 1) &&
                        (fabs(opinions[node1] - opinions[node2]) > input.epsilon) &&
                        (ran2(state) < input.r)
                            )
                {
                    //find candidates for new connection
                    /* count = 0; */
                    /* for (j=0; j< input.N; j++) */
                    /*     if (j != node1 && j != node2 && */
                    /*             (! isInList(&Nvec[node1], j)) && */
                    /*             (fabs(opinions[node1] - opinions[j]) < input.epsilon) */
                    /*             ) */
                    /*     { */
                    /*         auxN[count] = j; */
                    /*         count++; */
                    /*     } */

                    /* if (count > 0) */
                    /* { */
                    /*     removeConnection(Nvec, &node1, &node2); */

                    /*     node2 = auxN[ran2int(0,count,state)]; */
                    /*     appendList(&Nvec[node1],node2); */
                    /*     appendList(&Nvec[node2],node1); */
                    /* } */
                    /* else { */
                    /*     double probs[input.N]; */
                    /*     double totalProb = 0.0; */
                    /*     for (j=0; j< input.N; j++) { */
                    /*         if (j != node1 && j != node2 && */
                    /*                 (! isInList(&Nvec[node1], j))) */
                    /*         { */
                    /*             totalProb += 1./fabs(opinions[node1] - opinions[j]); */
                    /*             probs[count] = totalProb; */
                    /*             auxN[count] = j; */
                    /*             count++; */
                    /*         } */
                    /*     } */
                    /*     if (count > 0) { */
                    /*         double ranVal = ran2(state)*totalProb; */
                    /*         //now do a binary search for the selected index */
                    /*         int low,high,mid; */
                    /*         low = 0; */
                    /*         high = count; */
                    /*         mid = (high + low)/2; */
                    /*         while (high > low) { */
                    /*             if (probs[mid] > ranVal) */
                    /*                 high = mid; */
                    /*             else */
                    /*                 low = mid+1; */
                    /*             mid = (high + low)/2; */
                    /*         } */
                    /*         removeConnection(Nvec, &node1, &node2); */

                    /*         node2 = auxN[mid]; */
                    /*         appendList(&Nvec[node1],node2); */
                    /*         appendList(&Nvec[node2],node1); */
                    /*     } */
                    /* } */
                    /* partOpin[node1] = opinions[node1]; */

                    for (j=0; j< input.N; j++)
                    {
                        shuf1 = ran2int(0,input.N-j,state)+j;
                        shuf2 = auxN[j];
                        auxN[j] = auxN[shuf1];
                        auxN[shuf1] = shuf2;
                    }
                    for (j=0; j< input.N; j++)
                        if (auxN[j] != node1 && auxN[j] != node2 &&
                                (! isInList(&Nvec[node1], auxN[j])) &&
                                (fabs(opinions[node1] - opinions[auxN[j]]) < input.epsilon)
                                )
                        {
                            removeConnection(Nvec, &node1, &node2);
                            node2 = auxN[j];
                            appendList(&Nvec[node1],node2);
                            appendList(&Nvec[node2],node1);
                            break;
                        }
                }
                else
                {
                    field = 0.0;
                    for (j=0; j< Nvec[node1].leng ; j++)
                        field += (double) (epid[j] == 'I');
                    field *= input.w/((double) Nvec[node1].leng);

                    opinions[node1] += (ran2(state) * opinions[node2]) + field;
                    opinions[node1] = fmax(fmin(opinions[node1],1.0),-1.0);
                }
            }
        /* for (i=0; i< input.N; i++) */
        /*      opinions[node1] = partOpin[node1]; */
        //}
        //printf("%lf\n",opinions[*N+i]);
    }
}

void inteOpDisc (double *opinions, char *epid, list *Nvec, struct ran_state *state) {
    int i, node1, node2, index1, index2, j, shuf1, shuf2, count, k;
    int auxN[input.N];
    for (j=0; j< input.N; j++)
        auxN[j] = j;
    double field;
    for (i=0; i< input.N; i++)
    {
        node1 = ran2int(0,input.N,state);
        //case in witch there is no neighbours has to be ignored
        if (Nvec[node1].leng > 0) {
                index2 = ran2int(0,Nvec[node1].leng,state);
                node2 = Nvec[node1].list[index2];
                //if (i == node2) //printf("%i -> %i\n",i,node2);
                //if (*N <= node2) //printf("%i -> %i\n",i,node2);
                if ((Nvec[node1].leng > 1) &&
                        (Nvec[node2].leng > 1) &&
                        (fabs(opinions[node1] - opinions[node2]) > input.epsilon) &&
                        (ran2(state) < input.r)
                            )
                {
                    /* //find candidates for new connection */
                    /* count = 0; */
                    /* for (j=0; j< input.N; j++) */
                    /*     if (j != node1 && j != node2 && */
                    /*             (! isInList(&Nvec[node1], j)) && */
                    /*             (fabs(opinions[node1] - opinions[node2]) < input.epsilon) */
                    /*             ) */
                    /*     { */
                    /*         auxN[count] = j; */
                    /*         count++; */
                    /*     } */

                    /* if (count > 0) */
                    /* { */
                    /*     //removes the old connection */
                    /*     index1 = findApIndex(&Nvec[node1], &node2); */
                    /*     for (j=index2;j<Nvec[node1].leng;j++) */
                    /*         Nvec[node1].list[j] = Nvec[node1].list[j+1]; */
                    /*     Nvec[node1].leng--; */
                    /*     for (j=index1;j<Nvec[node2].leng;j++) */
                    /*         Nvec[node2].list[j] = Nvec[node2].list[j+1]; */
                    /*     Nvec[node2].leng--; */

                    /*     node2 = auxN[ran2int(0,count,state)]; */
                    /*     appendList(&Nvec[node1],node2); */
                    /*     appendList(&Nvec[node2],node1); */
                    /* } */
                    /* /1* partOpin[node1] = opinions[node1]; *1/ */

                    /* shuffle the list using fisher-yates */
                    for (j=0; j< input.N; j++)
                    {
                        shuf1 = ran2int(0,input.N-j,state)+j;
                        shuf2 = auxN[j];
                        auxN[j] = auxN[shuf1];
                        auxN[shuf1] = shuf2;
                    }
                    for (j=0; j< input.N; j++)
                    {
                        if (auxN[j] != node1 && auxN[j] != node2 &&
                                (! isInList(&Nvec[node1], auxN[j])) &&
                                (fabs(opinions[node1] - opinions[auxN[j]]) < input.epsilon)
                                )
                        {
                            /* removing old connection */
                            index1 = findApIndex(&Nvec[node1], &node2);
                            for (k=index2;k<Nvec[node1].leng;k++)
                                Nvec[node1].list[k] = Nvec[node1].list[k+1];
                            Nvec[node1].leng--;
                            for (k=index1;k<Nvec[node2].leng;k++)
                                 Nvec[node2].list[k] = Nvec[node2].list[k+1];
                            Nvec[node2].leng--;

                            /* adding the first valid reconection */
                            node2 = auxN[j];
                            appendList(&Nvec[node1],node2);
                            appendList(&Nvec[node2],node1);
                            break;
                        }
                    }
                }
                else
                {
                    field = 0.0;
                    for (j=0; j< Nvec[node1].leng ; j++)
                        field += (double) (epid[j] == 'I');
                    field *= input.w/((double) Nvec[node1].leng);

                    opinions[node1] += opinions[node2];
                    opinions[node1] += (double) (ran2(state) < field);
                    opinions[node1] = fmax(fmin(opinions[node1],1.0),-1.0);
                }
            }
        /* for (i=0; i< input.N; i++) */
        /*      opinions[node1] = partOpin[node1]; */
        //}
        //printf("%lf\n",opinions[*N+i]);
    }
}

void timeSeriesNW(void) {
    //starts writing the output file
    FILE * outputFile;
    char out_name[200];
    //printf("Starting timeseries for p=%f\n",input.p);
    sprintf(out_name, "SISrewire_l_%lf_w_%lf_r_%lf.dat",input.lambda,input.w,input.r);
    outputFile = fopen( out_name, "w" );

    struct ran_state state = {-time(NULL), 123456789, 0, {0}};

    list* Nvec;

    //allocates the vectors that carry the system data
    double *opinions;
    char *epid;
    //after the N-th vector component comes the the structure for time t+1
    opinions = (double *) malloc(sizeof(double)*input.N);
    //int *signs[input.N];
    epid = (char *)malloc(sizeof(char)*input.N);
    //printf("Getting memory\n");

    //check if the memory was allocated
    if (epid == NULL || opinions == NULL )
    {
        perror("error: malloc() \n");
        exit(EXIT_FAILURE);
    } //vectors sucessfully allocated

    //making the file header
    fprintf(outputFile,"#BCS with memory\n");
    fprintf(outputFile,"# Order parameter time evolution\n#$N=%i \\, t_{max}=%i \\, t_{stored}=%i \\, k=%g \\, D=%g \\, r=%g \\, \\epsilon=%g \\, \\lambda=%g netType=%i$\n", input.N, input.tmax, input.saved_steps, input.k, input.D, input.r, input.epsilon, input.lambda, input.netType);

    /* printf("Generating network\n"); */
    Nvec = randomGraph(input.N, input.k,&state);
    /* printf("Network generated\nInitial state\n"); */

    //just to get to the stationary state without kepping record
    int t,i,j;
    int count = 0;
    double mean_o, mean_I, mean_S, mean_R, mean_X;
    //printf("Serie temporal gerada.\n");
    //

    //reset them opinions
    for (i=0; i< input.N; i++)
    {
        opinions[i] = (ran2(&state) > input.D) ? -1.0 : 1.0;
        epid[i] = (ran2(&state) > 0.1) ? 'S' : 'I';

    }
    /*for (i=input.N/2; i< input.N; i++)
    {
        opinions[i] = ran2(&state);
        epid[i] = 'S';

    }*/
    for (i=0; i< input.N ; i++)
        count += Nvec[i].leng;

    for (t=0; t< input.tmax ; t++)
    {
        inteOprefRw (opinions, epid, Nvec, &state);
        inteEpid (opinions, Nvec, epid, &state);

        mean_o = 0.0;
        mean_I = 0.0;
        mean_S = 0.0;
        mean_R = 0.0;
        mean_X = 0.0;
        for (i=0; i< input.N ; i++)
        {
            mean_o += opinions[i];
            mean_I += (double) (epid[i] == 'I');
            mean_S += (double) (epid[i] == 'S');
            mean_R += (double) (epid[i] == 'R');
            for (j=0;j < Nvec[i].leng; j++)
                mean_X += fabs(opinions[i] - opinions[Nvec[i].list[j]]);
        }
        mean_o /= input.N;
        mean_I /= input.N;
        mean_S /= input.N;
        mean_R /= input.N;
        mean_X /= 2.*count;
        fprintf(outputFile, "%i\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\n",t, mean_o, mean_I , mean_S, mean_R, mean_X);

    }

        //
        //printf("Generating monte carlo\n");

    if (input.pHist == 1)
    {
        FILE * histFile;
        sprintf(out_name, "Hrewire_l_%lf.dat",input.lambda);
        histFile = fopen( out_name, "w" );

        //making the file header
        fprintf(histFile,"#BCS with memory histogram\n");
        for (i=0; i< input.N; i++)
            fprintf(histFile,"%f\n",opinions[i]);

        fclose(outputFile);
    }

    //printf("Fechando o arquivo de texto\n");
    fclose(outputFile);
    /* printf("Limpando lista de vizinhos\n"); */
    for (i=0; i< input.N; i++)
        free(Nvec[i].list);

    free(Nvec);
    free(opinions);
    free(epid);
}

void stationaryNW(double *opinions, char *epid, list *Nvec, struct ran_state *state) {
    //just to get to the stationary state without kepping record
    int t,i,j;
    double mean_o, mean_ep, mean_I, max_I, mean_S, mean_R;

    int count = 0;
    double mean_X = 0.0;
    double aux0, fneg;
    for (i=0; i< input.N ; i++)
        count += Nvec[i].leng;

    max_I  = 0.0;
    for (t=0; t< input.tmax ; t++)
    {
        mean_ep = 0.0;
        inteOprefRw (opinions, epid, Nvec, state);
        inteEpidCP (opinions, Nvec, epid, state);
        for (i=0; i< input.N ; i++)
            mean_ep += (double) (epid[i] == 'I');
        mean_ep /= input.N;
        max_I = fmax(mean_ep, max_I);
    }

    mean_o = 0.0;
    mean_I = 0.0;
    mean_S = 0.0;
    mean_R = 0.0;
    mean_X = 0.0;
    for (t=0; t< input.saved_steps; t++)
    {
        for (i=0;i < input.cnt; i++)
        {
            inteOprefRw (opinions, epid, Nvec, state);
            inteEpidCP (opinions, Nvec, epid, state);
        }

        mean_ep = 0.0;
        for (i=0; i< input.N ; i++)
        {
            mean_o += opinions[i];
            mean_ep += (double) (epid[i] == 'I');
            mean_S += (double) (epid[i] == 'S');
            mean_R += (double) (epid[i] == 'R');
            for (j=0;j < Nvec[i].leng; j++)
                mean_X += fabs(opinions[i] - opinions[Nvec[i].list[j]]);
        }
        mean_ep /= input.N;
        max_I = fmax(mean_ep, max_I);
        mean_I += mean_ep;
    }
    mean_o /= input.N*input.saved_steps;
    mean_I /= input.saved_steps;
    mean_S /= input.N*input.saved_steps;
    mean_R /= input.N*input.saved_steps;
    mean_X /= 2.0*count*input.saved_steps;

    #pragma omp critical
    {
        //starts writing the output file
        char output[200];
        sprintf(output, "%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%i",mean_o, mean_I , mean_S, mean_R,mean_X, max_I, connectedComponenets(Nvec, input.N));
        printFile(output);

    }
}

void setup_dcNet (void) {

    {
        struct ran_state state = {-(time(NULL)), 123456789, 0, {0} };

        list* Nvec;

        //allocates the vectors that carry the system data
        double *opinions;
        char *epid;
        int i;
        //after the N-th vector component comes the the structure for time t+1
        opinions = (double *) malloc(sizeof(double)*input.N);
        //int *signs[input.N];
        epid = (char *)malloc(sizeof(char)*input.N);
        //printf("Getting memory\n");

        //check if the memory was allocated
        if (epid == NULL || opinions == NULL )
        {
            perror("error: malloc() \n");
            exit(EXIT_FAILURE);
        } //vectors sucessfully allocated

        for (int samp=0; samp< input.samples; samp++)
        {
            //reset them opinions
            for (i=0; i< input.N; i++)
            {
                opinions[i] = (ran2(&state) > input.D) ? -1.0 : 1.0;
                epid[i] = (ran2(&state) > 0.1) ? 'S' : 'I';

            }
            /*for (i=input.N/2; i< input.N; i++)
            {
                opinions[i] = ran2(&state);
                epid[i] = 'S';

            }*/

            //printf("Generating graph\n");
            Nvec = randomGraph(input.N, input.k,&state);
            //printf("Network generated\nInitial state\n");

            //printf("Initial state\n");

            //printf("Starting dynamics\n");
            stationaryNW(opinions, epid, Nvec, &state);
            //timeSeriesNW(outputFile, opinions, signs, Nvec, &half);
            //samp = input.samples;
            /* printf("Limpando lista de vizinhos\n"); */
            for (i=0; i< input.N; i++)
                free(Nvec[i].list);

            free(Nvec);
        }
        /* printf("Limpando opiniões\n"); */
        free(opinions);
        free(epid);
    }
    /* printf("Fechando o arquivo de texto\n"); */
}

int main (int argc, char* argv[])
{
    int i;
    //loops over arguments changin the values of the given inputs
    for (i=1;i<argc;i++)
    {
        char text[20];
        strcpy(text,argv[i]);

        char *split = strtok(text,"=");

        //from here one the code will go throught the parameters changint the default value when necessary
        if (strcmp("N",split) == 0) input.N = (int) atof(strtok(NULL," "));
        else
        {
            if (strcmp("tmax",split) == 0) input.tmax = (int) atof(strtok(NULL," "));
            else {
                if (strcmp("saved_steps",split) == 0) input.saved_steps = (int) atof(strtok(NULL," "));
                else {
                    if (strcmp("cnt",split) == 0) input.cnt = (int) atof(strtok(NULL," "));
                    else {
                        if (strcmp("samples",split) == 0) input.samples = (int) atof(strtok(NULL," "));
                        else {
                            if (strcmp("N",split) == 0) input.N = (int) atof(strtok(NULL," "));
                            else {
                                if (strcmp("lambda",split) == 0) input.lambda = (double) atof(strtok(NULL," "));
                                else {
                                    if (strcmp("epsilon",split) == 0) input.epsilon = (double) atof(strtok(NULL," "));
                                    else {
                                        if (strcmp("k",split) == 0) input.k = (float) atof(strtok(NULL," "));
                                        else {
                                            if (strcmp("alpha",split) == 0) input.alpha = (double) atof(strtok(NULL," "));
                                            else {
                                                if (strcmp("phi",split) == 0) input.phi = (double) atof(strtok(NULL," "));
                                                else {
                                                if (strcmp("hist",split) == 0) input.pHist = (int) atof(strtok(NULL," "));
                                                    else {
                                                        if (strcmp("epsilon",split) == 0) input.epsilon = (double) atof(strtok(NULL," "));
                                                        else {
                                                            if (strcmp("r",split) == 0) input.r = (double) atof(strtok(NULL," "));
                                                            else {
                                                                if (strcmp("w",split) == 0) input.w = (double) atof(strtok(NULL," "));
                                                                else {
                                                                    if (strcmp("D",split) == 0) input.D = (double) atof(strtok(NULL," "));
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    if (input.samples == 0)
        timeSeriesNW();
    else
        setup_dcNet();

    return 0;
}
