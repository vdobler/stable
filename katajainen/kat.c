#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

typedef struct {
  int a;
  int b;
} pair;



#define bool int
#define TRUE 1
#define FALSE 0


/*********************************************************************
 * Definitions counters, which are used for counting the number of
 * comparisons and assignments
 */

#define COUNTERS 1


#ifdef COUNTERS
#define add_comparisons(x) (comparisons  += (x))
#else
#define add_comparisons(x) ;
#endif
#ifdef COUNTERS
#define add_assignments(x) (assignments  += (x))
#else
#define add_assignments(x) ;
#endif

#define FUNCTIONS  1 

#define RIGHT 1 


/*********************************************************************
 * This type of items we are sorting
 */

#define TYPE int

/*********************************************************************
 * In my PC I have minimum macro, but not in UNIX-machine
 */

#ifndef __TURBOC__

#define min(a,b) ( ((a)<(b)) ? (a) : (b) )

#endif


/*********************************************************************
 * GLOBAALIEN muuttujien esittely
 */

long comparisons;
long assignments;

/*********************************************************************
 */

void swap2(pair A[], int i, int j)
{
  pair temp = A[i];
  add_assignments(3);
  A[i] = A[j];
  A[j] = temp;
}

/*********************************************************************
 */

void reverse2(pair A[], int low, int high)
{
  int             t;
  for (t = 0; t < ((high - low + 1) >> 1); ++t)
    swap2(A, low + t, high - t);
}

/*********************************************************************
 */

void move_right(pair A[], int low, int high, int k)
{
  if (k > 0) {
    reverse2(A, low, high);
    reverse2(A, high + 1, high + k);
    reverse2(A, low, high + k);
  }
}


/*********************************************************************
 */

void change_places(pair A[], int buf, int xl, int xr)
{
  pair temp;

  add_assignments(1);
  /* The hole is created */
  temp = A[buf];

  while (xl < xr) {
    add_assignments(2);
    A[buf++] = A[xl];
    A[xl++] = A[buf];
  }
  add_assignments(2);
  A[buf] = A[xl];
  A[xl] = temp;
}



/*********************************************************************
  We define logarithm of 2 by ourselves when compiling with Turbo C */

double 
log2(double x)
{
  double          i = 1, j = 0;
  while (i < x) {
    i *= 2;
    j += 1;
  }
  return (j);
}

/*********************************************************************
  We define logarithm of 4 by ourselves */
double 
log4(double x)
{
  double          i = 1, j = 0;
  while (i < x) {
    i *= 4;
    j += 1;
  }
  return (j);
}

/*********************************************************************
  We define logarithm of 8 by ourselves */
double 
log8(double x)
{
  double          i = 1, j = 0;
  while (i < x) {
    i *= 8;
    j += 1;
  }
  return (j);
}




/*********************************************************************
 */

void 
insertion_sort_f(pair A[], int low, int high, int (*f) (pair *, pair *))
{
  int  i, j;
  pair temp;
  for (j = low + 1; j <= high; ++j) {
    add_assignments(1);
    temp = A[j];
    i = j - 1;	/* i >= low && A[i] > temp */
    while (i >= low && (f(&A[i], &temp) > 0)) {
      add_comparisons(1);
      A[i + 1] = A[i];
      --i;
    }
    add_comparisons(1);
    add_assignments(1);
    A[i + 1] = temp;
  }
}




/*********************************************************************
	The procedure merges subarrays A[xl..xr] and A[yl..yr] into the
	subarray A[buf..]. The elements in the result area are swapped with
	the merged elements. */

void 
straight_merge_f(pair A[], int buf, int xl, int xr, int yl, int yr,
		 int (*f) (pair *, pair *))
{
  pair temp;
  add_assignments(1);
  /*
   * The hole is created by storing the first element (A[buf]) of
   * result area. By that we get rid of swap-operations.
   */
  temp = A[buf];
  while (xl < xr && yl < yr) {
    add_comparisons(1);
    add_assignments(2);	/* A[xl] <= A[yl] */
    if (f(&A[xl], &A[yl]) <= 0) {
      A[buf++] = A[xl];
      A[xl++] = A[buf];
    } else {
      A[buf++] = A[yl];
      A[yl++] = A[buf];
    }
  }
  while (xl <= xr && yl <= yr) {
    add_comparisons(1);
    add_assignments(2);	/* A[xl] <= A[yl] */
    if (f(&A[xl], &A[yl]) <= 0) {
      A[buf++] = A[xl];
      A[xl] = ((xl < xr) ? A[buf] : temp);
      ++xl;
    } else {
      A[buf++] = A[yl];
      A[yl] = ((yl < yr) ? A[buf] : temp);
      ++yl;
    }
  }
  if (buf == xl || buf == yl)
    return;
  add_assignments(1);
  /* The hole is created again */
  temp = A[buf];
  if (yl > yr) {
    while (xl < xr) {
      add_assignments(2);
      A[buf++] = A[xl];
      A[xl++] = A[buf];
    }
    add_assignments(2);
    A[buf] = A[xl];
    A[xl] = temp;
  } else {
    while (yl < yr) {
      add_assignments(2);
      A[buf++] = A[yl];
      A[yl++] = A[buf];
    }
    add_assignments(2);
    A[buf] = A[yl];
    A[yl] = temp;
  }
}


/********************************************************************
	The procedure search the highest index k, yl <= k <= yr, by
	binarysearch for which the relation A[k] <= z is true. */

int 
high_binary_search_f(pair A[], int yl, int yr, pair z, int (*f) (pair *, pair *))
{
  int             mid, low = yl, high = yr;
  if (yr < yl)
    return (yl - 1);
  while ((high - low) >= 1) {
    mid = low + ((high - low + 1) >> 1);
    add_comparisons(1);	/* A[mid] <= z */
    if (f(&A[mid], &z) <= 0)
      low = mid;
    else
      high = mid - 1;
  }
  add_comparisons(1);	/* A[low] <= z */
  if (f(&A[low], &z) <= 0)
    return (low);
  else
    return (low - 1);
}

/*********************************************************************
 * Ohjelma lomittaa lohkon A[xl..xr] lohkoon A[yl..yr]. Lohkot ovat
 * vierekkaisia, eli xr+1 == yl. Ensimmaisen lohkon alkioiden paikat
 * etsitaan kayttaen puolitushakua. Alkiot siirretaan paikalle
 * listakaannoilla.
 */
void 
pardo_forward_block_merge_f(pair A[], int xl, int xr, int yl, int yr,
			    int (*f) (pair *, pair *))
{
  int             h;
  while (xr - xl >= 0) {
    h = high_binary_search_f(A, xr + 1, yr, A[xl], f);
    move_right(A, xl, xr, h - xr);
    xl = xl + h - xr + 1;
    xr = h;
  }
}

/*********************************************************************
	The procedure merges subarrays A[xl..xr] and A[yl..yr] into the
	subarray A[buf..] by binary merge. The elements in the result area are
	swapped with the merged elements. */

void 
binary_merge_f(pair A[], int buf, int xl, int xr, int yl, int yr,
	       int (*f) (pair *, pair *))
{
  int i, m, n;
  pair temp;
  m = xr - xl + 1;
  n = yr - yl + 1;
  add_assignments(1);
  /*
   * The hole is created by storing the first element (A[buf]) of
   * result area. By that we get rid of swap-operations.
   */
  temp = A[buf];
  while (m > 0 && n > 0) {
    if (m <= n) {
      i = yl + n / m - 1;
      add_comparisons(1);	/* A[xl] >= A[i] */
      if (f(&A[xl], &A[i]) >= 0) {
	while (yl < i) {
	  add_assignments(2);
	  A[buf++] = A[yl];
	  A[yl++] = A[buf];
	}
	add_assignments(2);
	A[buf++] = A[yl];
	A[yl] = ((yl < yr) ? A[buf] : temp);
	++yl;
	n = yr - yl + 1;
      } else {
	if ((yl - 1) != (i = high_binary_search_f(A, yl, i - 1, A[xl], f))) {
	  while (yl < i) {
	    add_assignments(2);
	    A[buf++] = A[yl];
	    A[yl++] = A[buf];
	  }
	  add_assignments(2);
	  A[buf++] = A[yl];
	  A[yl] = ((yl < yr) ? A[buf] : temp);
	  ++yl;
	  n = yr - yl + 1;
	}
	add_assignments(2);
	A[buf++] = A[xl];
	A[xl] = ((xl < xr) ? A[buf] : temp);
	++xl;
	--m;
      }
    } else {
      i = xl + m / n - 1;
      add_comparisons(1);	/* A[yl] >= A[i] */
      if (f(&A[yl], &A[i]) >= 0) {
	while (xl < i) {
	  add_assignments(2);
	  A[buf++] = A[xl];
	  A[xl++] = A[buf];
	}
	add_assignments(2);
	A[buf++] = A[xl];
	A[xl] = ((xl < xr) ? A[buf] : temp);
	++xl;
	m = xr - xl + 1;
      } else {
	if ((xl - 1) != (i = high_binary_search_f(A, xl, i - 1, A[yl], f))) {
	  while (xl < i) {
	    add_assignments(2);
	    A[buf++] = A[xl];
	    A[xl++] = A[buf];
	  }
	  add_assignments(2);
	  A[buf++] = A[xl];
	  A[xl] = ((xl < xr) ? A[buf] : temp);
	  ++xl;
	  m = xr - xl + 1;
	}
	add_assignments(2);
	A[buf++] = A[yl];
	A[yl] = ((yl < yr) ? A[buf] : temp);
	++yl;
	--n;
      }
    }
  }
  if (buf == xl || buf == yl)
    return;
  add_assignments(1);
  /* The hole is created again */
  temp = A[buf];
  if (n <= 0) {
    while (xl < xr) {
      add_assignments(2);
      A[buf++] = A[xl];
      A[xl++] = A[buf];
    }
    add_assignments(2);
    A[buf] = A[xl];
    A[xl] = temp;
  } else {
    while (yl < yr) {
      add_assignments(2);
      A[buf++] = A[yl];
      A[yl++] = A[buf];
    }
    add_assignments(2);
    A[buf] = A[yl];
    A[yl] = temp;
  }
}

/*********************************************************************
	Non-recursive merge sort of the subarray A[left..rigth] by using
	locations A[buf..left1-1] or A[rigth+1..buf] as the extra space. If
	to_rigth is 0 then the buffer-area is in A[buf..left-1], else if
	to_rigth is 1 then it is in A[buf..buf+rigth-left]. The assumption is
	that |A[buffer-area]| = |A[left..rigth]|. */

void 
two_way_msort_f(pair A[], int left, int rigth, int buf,
		int to_rigth, int (*f) (pair *, pair *))
{
  int i, cursor, left1, rigth1, left2, rigth2, runs,  size, nelements;
  nelements = rigth - left + 1;
  size = 1;
  while (size < nelements) {
    runs = nelements / size;
    if (nelements % size)
      ++runs;
    cursor = left;
    for (i = 1; i <= (runs >> 1); ++i) {
      left1 = cursor;
      left2 = left1 + size;
      rigth1 = left2 - 1;
      rigth2 = (((rigth1 + size) > rigth) ? rigth : rigth1 + size);
      straight_merge_f(A, buf, left1, rigth1, left2, rigth2, f);
      buf += (rigth2 - left1 + 1);
      cursor += (rigth2 - left1 + 1);
    }
    if (runs % 2 == 1) {	/* Number of runs was odd. */
      change_places(A, buf, cursor, rigth);
    }
    buf = left;
    if (to_rigth) {
      left += nelements;
      rigth += nelements;
      to_rigth = 0;
    } else {
      left -= nelements;
      rigth -= nelements;
      to_rigth = 1;
    }
    size = size << 1;
  }
}



/*********************************************************************
	The procedure searches the smallest elements of two blocks
	A[left..cut] and A[cut+1..rigth]. */

void 
find_smallest_elements_f(pair A[], int left, int cut, int rigth,
			 int *last1, int *last2, int (*f) (pair *, pair *))
{
  int             bufsize, leftblocksize;

  bufsize = 0;
  leftblocksize = cut - left + 1;
  *last1 = left - 1;
  *last2 = cut;

  while (*last2 < rigth && bufsize < leftblocksize) {
    add_comparisons(1);
    /* A[*last1+1] <= A[*last2+1] */
    if (f(&A[*last1 + 1], &A[*last2 + 1]) <= 0) {
      *last1 = *last1 + 1;
      --leftblocksize;
    } else {
      *last2 = *last2 + 1;
    }
    ++bufsize;
  }
}


/*********************************************************************
 */

void 
simple_tkp_f(pair A[], int left, int rigth, int (*f) (pair *, pair *))
{
  int             start, end, buf, cut, nelements, half, nrounds_is_odd;
  cut = rigth;
  nelements = cut - left + 1;
  half = (nelements >> 1);
  if (nelements == 1)
    return;
  nrounds_is_odd = (int) ceil(log2((double) half)) % 2;
  if (nrounds_is_odd) {	/* Merging first to rigth. */
    buf = cut - half + 1;
    end = buf - 1;
    start = buf - half;
  } else {		/* Merging first to left. */
    end = cut;
    start = end - half + 1;
    buf = start - half;
  }
  two_way_msort_f(A, start, end, buf, nrounds_is_odd, f);
  cut = cut - half;
  nelements = cut - left + 1;
  half = (nelements >> 1);
  /* Handling of the	first n - O(sqrt(n)) elements */
  while (half >= (int) sqrt((double) rigth - left + 1)) {
    nrounds_is_odd = (int) ceil(log2((double) half)) % 2;
    if (nrounds_is_odd) {	/* Merging first to left. */
      buf = left;
      start = buf + half;
      end = start + half - 1;
      nrounds_is_odd = 0;
    } else {	/* Merging first to rigth. */
      start = left;
      end = start + half - 1;
      buf = end + 1;
      nrounds_is_odd = 1;
    }
    two_way_msort_f(A, start, end, buf, nrounds_is_odd, f);
    binary_merge_f(A, cut - half + 1, left, left + half - 1, cut + 1,
		   rigth, f);
    cut = cut - half;
    nelements = cut - left + 1;
    half = (nelements >> 1);
  }
  /* Handling of the first O(sqrt(n)) elements */
  insertion_sort_f(A, left, left + nelements - 1, f);
  pardo_forward_block_merge_f(A, left, left + nelements - 1,
			      left + nelements, rigth, f);
}

/*********************************************************************
	Sort the area A[left..rigth] with minimum space. First the second half
	is sorted by using first half as the auxiliary space. Then the first
	quarter is sorted by using the second quarter as the auxiliary space,
	and merged with the sorted second half, into locations starting from
	the second quarter. Now the first quarter is handled similarly,
	leaving 1/8 of the array unsorted, and so on.  The process stops after
	log N phases, leaving the first element unsorted. It is finally
	inserted to the correct position.  */

void 
tkp_f(pair A[], int left, int rigth, int (*f) (pair *, pair *))
{
  int             start, end, buf, cut, nelements, half, nrounds_is_odd;
  cut = rigth;
  nelements = cut - left + 1;
  half = (nelements >> 1);
  if (nelements == 1)
    return;
  nrounds_is_odd = (int) ceil(log2((double) half)) % 2;
  if (nrounds_is_odd) {	/* Merging first to rigth. */
    buf = cut - half + 1;
    end = buf - 1;
    start = buf - half;
  } else {		/* Merging first to left. */
    end = cut;
    start = end - half + 1;
    buf = start - half;
  }
  two_way_msort_f(A, start, end, buf, nrounds_is_odd, f);
  cut = cut - half;
  nelements = cut - left + 1;
  half = (nelements >> 1);
  /* Handling of the	first n - O(n / log n) elements */
  while (nelements >= (rigth - left + 1) / log2((double) rigth - left + 1)) {
    nrounds_is_odd = (int) ceil(log2((double) half)) % 2;
    if (nrounds_is_odd) {	/* Merging first to left. */
      buf = left;
      start = buf + half;
      end = start + half - 1;
      nrounds_is_odd = 0;
    } else {	/* Merging first to rigth. */
      start = left;
      end = start + half - 1;
      buf = end + 1;
      nrounds_is_odd = 1;
    }
    two_way_msort_f(A, start, end, buf, nrounds_is_odd, f);
    binary_merge_f(A, cut - half + 1, left, left + half - 1, cut + 1,
		   rigth, f);
    cut = cut - half;
    nelements = cut - left + 1;
    half = (nelements >> 1);
  }
  {
    int             last1, last2;
    /* Handling of the first O(n / log n) elements */
    simple_tkp_f(A, left, cut, f);
    /*
     * Searching of the smallest O(n / log n) elements of two
     * ordererd block, the first is left..cut and the other is
     * cut+1 ..rigth
     */
    find_smallest_elements_f(A, left, cut, rigth, &last1, &last2, f);
    /*
     * The smallest elements are in left..last1 and cut+1..last2
     * If the smallest elements are in front already, the job is
     * done.
     */
    if (last1 != cut)
      if (last2 == rigth)
	move_right(A, left, cut, rigth - cut);
      else {
	straight_merge_f(A, left, last1 + 1, cut, 1, 0, f);
	binary_merge_f(A, last2 - cut + last1 + 1, left,
		       left + cut - last1 - 1, last2 + 1, rigth, f);
	simple_tkp_f(A, left, last1 - left + 1 + last2 - cut, f);
      }
  }
}

int integerLess(int *a, int *b) {
  return (*a) - (*b);
}

int pairLess(pair *x, pair *y) {
  return  (x->a) - (y->a);
}




int main() {
  pair A[1000000];
  int N =1000000;
  
  int i;
  for (i=0; i<N; i++) {
    A[i].a = rand();
    A[i].b = i;
  }

  tkp_f(A, 0, N-1, pairLess);

  // sorted?
  for (i=1; i<N; i++) {
    if (A[i-1].a > A[i].a) {
      printf("Bad Sort: %d/%d  %d %d\n", i-1,i, A[i-1].a, A[i].a);
    }
  }

  // stable?
  int lastA, lastB; 
  lastA = -1;
  lastB = 0;
  for (i = 0; i < N; i++) {
    break;
    if (lastA != A[i].a) {
      lastA = A[i].a;
      lastB = A[i].b;
      continue;
    }
    if (A[i].b <= lastB) {
      printf("Bad Stable %d: %d %d\n", i, lastB, A[i].b);
    }
    lastB = A[i].b;
  }
  
  printf("n = %d   comparisons = %.2e  assignments = %.2e\n",
	 N, (double)comparisons, (double)assignments);
}
