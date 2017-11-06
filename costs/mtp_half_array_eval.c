// $Id: mtp_half_array_eval.c 1956 2017-10-26 12:41:10Z fabien $

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <strings.h>
#include <assert.h>
#include <math.h>
#include <sys/time.h> // srandom() on time
#include <pthread.h>

typedef unsigned long long int uint64;
typedef unsigned char uint8;

#define MAX_LL 20
#define MAX_PHIS (1 << MAX_LL)

static void init_PHI(int l, int PHI[l])
{
  for (int i = 1; i < l; i++)
  {
    double u = (double) random() / RAND_MAX;
    // quadratic bias
    double r = (1.0 - u * u);
    PHI[i] = (int) ((i - 1) * r);
    assert(0 <= PHI[i] && PHI[i] < i);
  }
}

static uint64 one_cost(int n, int i, uint8 used[], int l, int PHI[l])
{
  // if available, nothing to do
  if (i < n || i % 2 == 0 || used[i])
    return 0;
  // else compute element and dependencies
  used[i] = 1;
  int phi_i = PHI[i];
  uint64 cost = 1;
  switch (n)
  {
  case 6:
    cost += one_cost(n, (3 * phi_i + i) / 4, used, l, PHI);
  case 5:
    cost += one_cost(n, (phi_i + 3 * i) / 4, used, l, PHI);
  case 4:
    // cost += one_cost(n, 3 * (i-1) / 4, used);
    // BAD: stops easily if n is odd
    // cost += one_cost(n, i-n, used, l, PHI);
    // hmmm... special property, easy to stop with periodic clusters
    // cost += one_cost(n, i-n + n%2, used, l, PHI);
    cost += one_cost(n, 7 * i / 8, used, l, PHI);
  case 3:
    // cost += one_cost(n, 3*phi_i / 4, used);
    cost += one_cost(n, (phi_i + i) / 2, used, l, PHI);
  case 2:
    cost += one_cost(n, phi_i, used, l, PHI);
  case 1:
    cost += one_cost(n, i-1, used, l, PHI);
    break;
  default:
    fprintf(stderr, "unexpected n=%d\n", n);
    exit(2);
  }
  return cost;
}

static double average_cost(int n, int l)
{
  uint64 cost = 0;
  uint8 used[l];
  int PHI[l];
  init_PHI(l, PHI);
  for (int i = n; i < l; i++)
  {
    bzero(used, sizeof(used));
    cost += one_cost(n, i, used, l, PHI);
  }
  fputc('.', stderr);
  return (double) cost / l;
}

#define N_THREADS 4

typedef struct
{
  int n;
  int l;
  int N;
  double sum;
  double sum2;
} stuff_t;

static stuff_t * collect_stats(stuff_t * stuff)
{
  stuff->sum = 0.0;
  stuff->sum2 = 0.0;
  for (int i = 0; i < stuff->N; i++)
  {
    double cost = average_cost(stuff->n, stuff->l);
    stuff->sum += cost;
    stuff->sum2 += cost * cost;
  }
  return stuff;
}

int main(int argc, char* argv[])
{
  int n, ll, N, nthreads;
  if (argc != 4 && argc != 5)
  {
    fprintf(stderr, "usage: %s n ll N [nthreads]\n", argv[0]);
    return 1;
  }

  n = atoi(argv[1]);
  assert(2 <= n && n <= 6);

  ll = atoi(argv[2]);
  assert(5 <= ll && ll <= MAX_LL);

  N = atoi(argv[3]);
  assert(1 <= N);

  if (argc == 5)
    nthreads = atoi(argv[4]);
  else
    nthreads = N_THREADS;
  assert(nthreads >= 1);

  int l = 1 << ll;
  fprintf(stdout, "N=%d n=%d l=%d (2^%d)\n", N, n, l, ll);

  // random generator seed
  struct timeval tv;
  gettimeofday(&tv, NULL);
  srandom(tv.tv_sec * 1000000 + tv.tv_usec);

  double sum = 0.0, sum2 = 0.0;

  if (nthreads > 1)
  {
    // start some threads
    pthread_t threads[nthreads];
    stuff_t stuffs[nthreads];

    int NT = (N + nthreads - 1) / nthreads;
    N = NT * nthreads;

    for (int i = 0; i < nthreads; i++)
    {
      stuffs[i].n = n;
      stuffs[i].l = l;
      stuffs[i].N = NT;
      int r = pthread_create(&threads[i], NULL,
                             (void*(*)(void*))collect_stats, &stuffs[i]);
      assert(r == 0);
    }

    for (int i = 0; i < nthreads; i++)
    {
      int r = pthread_join(threads[i], NULL);
      assert(r == 0);
      sum += stuffs[i].sum;
      sum2 += stuffs[i].sum2;
    }
  }
  else // no threads
  {
    stuff_t stuff;
    stuff.n = n;
    stuff.l = l;
    stuff.N = N;
    collect_stats(&stuff);
    sum += stuff.sum;
    sum2 += stuff.sum2;
  }

  fputc('\n', stderr);

  double average = sum / N;
  double stddev = sqrt(( sum2 - sum * sum / N)/N);
  fprintf(stdout,
          "N=%d n=%d l=%d (2^%d) cost: %.3f +- %.3f\n",
          N, n, l, ll, average, stddev);

  return 0;
}
