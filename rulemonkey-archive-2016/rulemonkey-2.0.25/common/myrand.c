/*
** Copyright 2008 Northern Arizona University and
** The Translational Genomics Research Institute
**
**  This file is part of Rulemonkey.
**
**  Rulemonkey is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Rulemonkey is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Rulemonkey.  If not, see <http://www.gnu.org/licenses/>.
**
*/


#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0;
#endif

#include "myrand.h"
#include "output.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/(MBIG))

#define MT_ARRAY_SIZE (1024 * 2 * 2)

static int mt_array_size;
static int mt_array_avail;
static double *mt_array;

static int inext,inextp;
static long ma[56]; /* The value 56 (range ma[1..55]) is special and */
static int iff=0;   /* should not be modified; see Knuth. */
static int myseed=MSEED;
static long idum;

/*
** According to Knuth, any large MBIG, and any smaller (but still large)
** MSEED can be substituted  for the above values.
*/
void ran3_seed(int seed)
{
  myseed = seed;

  while (myseed < 1000000)
  {
    myseed += 100000;
  }

  while (myseed > MSEED)
  {
    myseed -= 1000;
  }

  // Initialized the sequence
  ran3(-1);
}

/*
** Returns a uniform random deviate between 0.0 and 1.0. Set idum to any
** negative value to initialize or reinitialize the sequence.
** From Numerical Recipies in C:
** http://www.library.cornell.edu/nr/bookcpdf/c7-1.pdf
*/
float ran3()
{
  long mj,mk;
  int i,ii,k;
  if (idum < 0 || iff == 0)
  {
    /* Initialization. */
    iff=1;
    /* Initialize ma[55] using the seed idum and the */
    /* large number MSEED. */
    mj=labs(myseed-labs(idum));
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++)
    {
      /* Now initialize the rest of the table, */
      /* in a slightly random order, */
      /* with numbers that are not especially random. */
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }

    /* We randomize them by warming up the generator */
    for (k=1;k<=4;k++)
    {
      for (i=1;i<=55;i++)
      {
        ma[i] -= ma[1+(i+30) % 55];
        if (ma[i] < MZ) ma[i] += MBIG;
      }
    }

    inext=0;   /* Prepare indices for our first generated number. */
    inextp=31; /* The constant 31 is special; see Knuth. */
    idum=1;
  }

  /* Here is where we start, except on initialization. */
  /* Increment inext and inextp, wrapping around */
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1; /* 56 to 1. */

  /* Generate a new random number subtractively. */
  mj=ma[inext]-ma[inextp];

  /* Be sure that it is in range. */
  if (mj < MZ) mj += MBIG;

  /* Store it, and output the derived uniform deviate. */
  ma[inext]=mj;
  return mj*FAC;
}

int
mysrand_mt(const uint32_t arg_seed)
{
  uint32_t seed = arg_seed;

  if (0 == seed)
  {
    // No seed specified, so use time to initialize.
    seed = time(NULL);
  }
  debug_printf(0, "Random number seed used: %d\n", seed);
  init_gen_rand(seed);

  mt_array_size = get_min_array_size();
  if (mt_array_size < MT_ARRAY_SIZE)
  {
    mt_array_size = MT_ARRAY_SIZE;
  }

  // Allocate aligned memory to contain generated random numbers.
  // Code taken from dSFMT:/html/howto-compile.html (sample2.c)
#if defined(__APPLE__) || \
    (defined(__FreeBSD__) && __FreeBSD__ >= 3 && __FreeBSD__ <= 6)
    mt_array = malloc(sizeof(double) * mt_array_size);
    if (mt_array == NULL) {
        error_printf("can't allocate memory.\n");
        return 0;
    }
#elif defined(_POSIX_C_SOURCE)
    if (posix_memalign((void **)&mt_array, 16,
                       sizeof(double) * mt_array_size) != 0) {
        error_printf("can't allocate memory.\n");
        return 0;
    }
#elif defined(__GNUC__) && (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 3))
    mt_array = memalign(16, sizeof(double) * mt_array_size);
    if (mt_array == NULL) {
        error_printf("can't allocate memory.\n");
        return 0;
    }
#else /* in this case, gcc doesn't support SSE2 */
    mt_array = malloc(sizeof(double) * mt_array_size);
    if (mt_array == NULL) {
        error_printf("can't allocate memory.\n");
        return 0;
    }
#endif

  mt_array_avail = 0;

  return 1;
}

/**
** Returns a random number in the range [0, 1)
** Example (Choose an index to an array with 10 elements):
** index = myrand_mt() * 10;
** Note that rounding is not needed if used in this way.
*/
double
myrand_mt()
{
  if (mt_array_avail == 0)
  {
    fill_array_close_open(mt_array, mt_array_size);
    mt_array_avail = mt_array_size;
  }

  return mt_array[--mt_array_avail];
}

void
myrand_end()
{
  free(mt_array);
}

