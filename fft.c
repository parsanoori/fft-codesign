#include <stdlib.h>
#include <stdio.h>
#include <math.h>

struct k_j {
    unsigned long long int k, j;
};

typedef struct k_j K_J;

int to_fixed(double input)
{
    return (int)(input * 16384);
}

double to_double(int input)
{
    return ((double)input / 268435456.0);
}


// bit_reversal_inner_loop memory map
volatile unsigned long long int *bit_reversal_inner_loop_in_k = (unsigned long long int *) 0x80000000;
volatile unsigned long long int *bit_reversal_inner_loop_in_j = (unsigned long long int *) 0x80000008;
volatile unsigned long long int *bit_reversal_inner_loop_out_k = (unsigned long long int *)0x80000010;
volatile unsigned long long int *bit_reversal_inner_loop_out_j = (unsigned long long int *)0x80000018;
volatile unsigned int *bit_reversal_inner_loop_load = (unsigned int *)0x80000020;
volatile unsigned int *bit_reversal_inner_loop_done = (unsigned int *)0x80000024;

// fixed_point_mult_add memory map
volatile unsigned long long int *fixed_point_mult_add_in_a = (unsigned long long int *)0x80000028;
volatile unsigned long long int *fixed_point_mult_add_in_b = (unsigned long long int *)0x80000030;
volatile unsigned long long int *fixed_point_mult_add_in_c = (unsigned long long int *)0x80000038;
volatile unsigned long long int *fixed_point_mult_add_in_d = (unsigned long long int *)0x80000040;
volatile unsigned long long int *fixed_point_mult_add_out = (unsigned long long int *)0x80000048;
volatile unsigned int *fixed_point_mult_add_sign = (unsigned int *)0x80000050;
volatile unsigned int *fixed_point_mult_add_load = (unsigned int *)0x80000054;
volatile unsigned int *fixed_point_mult_add_done = (unsigned int *)0x80000058;

// prototypes
K_J bit_reversal_inner_loop(unsigned long long int, unsigned long long int);
double fixed_point_mult_add(unsigned long long, unsigned long long, unsigned long long, unsigned long long);

/*
   This computes an in-place complex-to-complex FFT
   x and y are the real and imaginary arrays of 2^m points.
   dir =  1 gives forward transform
   dir = -1 gives reverse transform
*/
short FFT(double x[8], double y[8])
{
    int m = 3;
    int dir = 1;
    long n, i, i1, j, k, i2, l, l1, l2;
    double c1, c2, tx, ty, t1, t2, u1, u2, z;

    /* Calculate the number of points */

    /*
    n = 1;
    for (i = 0; i < m; i++)
        n *= 2;
    */
    n = 1 << m;

    /* Do the bit reversal */
    i2 = n >> 1;
    j = 0;
    for (i = 0; i < n - 1; i++)
    {
        if (i < j)
        {
            tx = x[i];
            ty = y[i];
            x[i] = x[j];
            y[i] = y[j];
            x[j] = tx;
            y[j] = ty;
        }

        /*
        k = i2;
        while (k <= j)
        {
            j -= k;
            k >>= 1;
        }
        j += k;
        */
        K_J kj = bit_reversal_inner_loop(i2, j);
        k = kj.k;
        j = kj.j;
    }

    /* Compute the FFT */
    c1 = -1.0;
    c2 = 0.0;
    l2 = 1;
    for (l = 0; l < m; l++)
    {
        l1 = l2;
        l2 <<= 1;
        u1 = 1.0;
        u2 = 0.0;
        for (j = 0; j < l1; j++)
        {
            for (i = j; i < n; i += l2)
            {
                i1 = i + l1;
                //t1 = u1 * x[i1] - u2 * y[i1];
                //t2 = u1 * y[i1] + u2 * x[i1];
                t1 = fixed_point_mult_add(to_fixed(u1), to_fixed(x[i1]), to_fixed(-u2), to_fixed(y[i1]));
                t2 = fixed_point_mult_add(to_fixed(u1), to_fixed(y[i1]), to_fixed(u2), to_fixed(x[i1]));
                x[i1] = x[i] - t1;
                y[i1] = y[i] - t2;
                x[i] += t1;
                y[i] += t2;
            }
            //z = u1 * c1 - u2 * c2;
            //u2 = u1 * c2 + u2 * c1;
            z = fixed_point_mult_add(to_fixed(u1), to_fixed(c1), to_fixed(-u2), to_fixed(c2));
            u2 = fixed_point_mult_add(to_fixed(u1), to_fixed(c2), to_fixed(u2), to_fixed(c1));
            u1 = z;
        }
        c2 = sqrt((1.0 - c1) / 2.0);
        if (dir == 1)
            c2 = -c2;
        c1 = sqrt((1.0 + c1) / 2.0);
    }

    /* Scaling for forward transform */
    if (dir == 1)
    {
        for (i = 0; i < n; i++)
        {
            x[i] /= n;
            y[i] /= n;
        }
    }

    return 0;
}

void show(const char *s, double bufr[], double bufi[])
{
    printf("%s", s);
    for (int i = 0; i < 8; i++)
        printf("%lf,%lf\n", bufr[i], bufi[i]);
}

int main()
{
    double bufr[] = {1, 1, 0, 0, 1, 1, 0, 0};
    double bufi[] = {0, 0, 0, 0, 0, 0, 0, 0};

    int ret_val = 0;

    show("\nData:\n", bufr, bufi);
    FFT(bufr, bufi);
    show("\nFFT:\n", bufr, bufi);

    getchar();

    return ret_val;
}

K_J bit_reversal_inner_loop(unsigned long long int k, unsigned long long int j)
{
    *bit_reversal_inner_loop_in_k = k;
    *bit_reversal_inner_loop_in_j = j;
    *bit_reversal_inner_loop_load = 1;
    while (*bit_reversal_inner_loop_done != 1)
        ;
    K_J result;
    result.k = *bit_reversal_inner_loop_out_k;
    result.j = *bit_reversal_inner_loop_out_j;
    *bit_reversal_inner_loop_load = 0;
    // printf("[debug] bit_reversal_inner_loop(%lld, %lld) = %lld\n", a, b, result);
    return result;
}


double fixed_point_mult_add(unsigned long long a, unsigned long long b, unsigned long long c, unsigned long long d)
{
    *fixed_point_mult_add_in_a = a;
    *fixed_point_mult_add_in_b = b;
    *fixed_point_mult_add_in_c = c;
    *fixed_point_mult_add_in_d = d;
    *fixed_point_mult_add_load = 1;
    while (*fixed_point_mult_add_done != 1)
        ;
    long long int result = *fixed_point_mult_add_out;
    if (*fixed_point_mult_add_sign == 1)
	result *= -1;
    *fixed_point_mult_add_load = 0;
    // printf("[debug] fixed_point_mult_add(%lld, %lld, %lld) = %lld\n", b, e, m, result);
    return to_double(result);
}


