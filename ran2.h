#define RAN2_H _

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)

/* ran2() - Return a random floating point value between 0.0 and
   1.0 exclusive.  If idum is negative, a new series starts (and
   idum is made positive so that subsequent teste1s using an unchanged
   idum will continue in the same sequence). */

struct ran_state {
    long idum;
    long idum2;
    long iy;
    long iv[NTAB];
};

float ran2(struct ran_state *state)
{
    int j;
    long k;
    /*static long idum2 = 123456789;
    static long iy = 0;
    static long iv[NTAB];*/
    float temp;

    if (state->idum <= 0)                               /* initialize */
    {
        if (-(state->idum) < 1)                           /* prevent idum == 0 */
            state->idum = 1;
        else
            state->idum = -(state->idum);                         /* make idum positive */
        state->idum2 = (state->idum);
        for (j = NTAB + 7; j >= 0; j--)             /* load the shuffle table */
        {
            k = (state->idum) / IQ1;
            state->idum = IA1 * (state->idum - k*IQ1) - k*IR1;
            if (state->idum < 0)
                state->idum += IM1;
            if (j < NTAB)
                state->iv[j] = state->idum;
        }
        state->iy = state->iv[0];
    }

    k = (state->idum) / IQ1;
    state->idum = IA1 * (state->idum - k*IQ1) - k*IR1;
    if (state->idum < 0)
        state->idum += IM1;
    k = state->idum2/IQ2;
    state->idum2 = IA2 * (state->idum2 - k*IQ2) - k*IR2;
    if (state->idum2 < 0)
        state->idum2 += IM2;
    j = state->iy / NDIV;
    state->iy = state->iv[j] - state->idum2;
    state->iv[j] = state->idum;
    if (state->iy < 1)
        state->iy += IMM1;
    if ((temp = AM * state->iy) > RNMX)
        return RNMX;                                /* avoid endpoint */
    else
        return temp;
}

int ran2int (int min, int max, struct ran_state *state) {
    int val = (int)(min + (ran2(state)*(max-min)));
    return val - (val == max);
}

void shuffle(int *array, int n, struct ran_state *state) {
    if (n > 1)
    {
        int i;
        for (i = 0; i < n - 1; i++)
        {
          int j = i + ran2int(0,32767,state) / (32767 / (n - i) + 1);
          int t = array[j];
          array[j] = array[i];
          array[i] = t;
        }
    }
}
