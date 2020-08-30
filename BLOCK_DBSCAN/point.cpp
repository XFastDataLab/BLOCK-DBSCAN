/*
 *
 * Copyright (C) 2017-  Yewang Chen<ywchen@hqu.edu.cn;nalandoo@gmail.com;1611414017@hqu.edu.cn>
 * License: GPL v1
 * This software may be modified and distributed under the terms
 * of license.
 *
 */
#include "point.h"
#define NDEBUG
#include<assert.h>
#include<math.h>
#include <string.h>
#include <errno.h>

const int batch = 120;//must be a multiple of 8

int point_len = 0;
//typedef float v4sf __attribute__ ((mode(V4SF)));


int posix_memalign(void **mptr, size_t, size_t bytes) {
  *mptr = malloc(bytes);
  return *mptr ? 0 : ENOMEM;
}


//Assumption: points are a multiples of 8 long
float distance_1(float* p1, float* p2)
{
//Euclidean distance
    float sum = 0.0;
    for(int i = 0; i < DIME_NUM-1;i++)
    {
        float d1 =p1[i] -p2[i];
        d1 *= d1;
        sum = sum+d1;
    }
  return sqrt(sum);

//Chebyshev distance
//    float sum =0;
//    for(int i = 0; i < DIME_NUM-1;i++)
//    {
//        float d2 =fabsf(p2[i] -p1[i]);
//        if(d2 > sum)
//            sum = d2;
//    }
//    return sum;
}
/*
//Assumption: points are a multiples of 8 long
float sse_distance(point p1, point p2, float upper_bound)
{
  v4sf sum = {0.,0.,0.,0.};
  float *end = p1 + point_len;
  upper_bound *= upper_bound;
  for (float *batch_end = p1 + batch; batch_end <= end; batch_end += batch)
    {
      for (; p1 != batch_end; p1+=8, p2+=8)
	{
	  v4sf v1 = __builtin_ia32_loadaps(p1);
	  v4sf v2 = __builtin_ia32_loadaps(p2);
	  v4sf v3 = __builtin_ia32_loadaps(p1+4);
	  v4sf v4 = __builtin_ia32_loadaps(p2+4);
	  v1 = __builtin_ia32_subps(v1, v2);
	  v3 = __builtin_ia32_subps(v3, v4);
	  v1 = __builtin_ia32_mulps(v1, v1);
	  v3 = __builtin_ia32_mulps(v3, v3);
	  v1 = __builtin_ia32_addps(v1,v3);
	  sum = __builtin_ia32_addps(sum,v1);
	}
      v4sf temp = __builtin_ia32_addps(sum,__builtin_ia32_shufps(sum,sum,14));
      temp = __builtin_ia32_addss(temp,__builtin_ia32_shufps(temp,temp,1));
      if (((float *)&temp)[0] > upper_bound)
	{
	  temp = __builtin_ia32_sqrtss(temp);
	  return ((float *)&temp)[0];
	}
    }
  for (; p1 != end; p1+=8, p2+=8)
    {
      v4sf v1 = __builtin_ia32_loadaps(p1);
      v4sf v2 = __builtin_ia32_loadaps(p2);
      v4sf v3 = __builtin_ia32_loadaps(p1+4);
      v4sf v4 = __builtin_ia32_loadaps(p2+4);
      v1 = __builtin_ia32_subps(v1, v2);
      v3 = __builtin_ia32_subps(v3, v4);
      v1 = __builtin_ia32_mulps(v1, v1);
      v3 = __builtin_ia32_mulps(v3, v3);
      v1 = __builtin_ia32_addps(v1,v3);
      sum = __builtin_ia32_addps(sum,v1);
    }
  sum = __builtin_ia32_addps(sum,__builtin_ia32_shufps(sum,sum,14));
  sum = __builtin_ia32_addss(sum,__builtin_ia32_shufps(sum,sum,1));
  sum = __builtin_ia32_sqrtss(sum);
  return ((float *) & sum)[0];
}*/

/*
float distance(point p1, point p2, float upper_bound)
{
  return fabsf(p1 - p2);
}

v_array<point> parse_points(FILE *input)
{
  v_array<point> ret;
  for (int i = 0; i< 1000; i++)
    push(ret,(float) i);
  return ret;
}

void print(point &p)
{
  printf("%f ",p);
  printf("\n");
}

*/

v_array<point > parse_points(FILE *input)
{
  v_array<point > parsed;
  char c;
  v_array<float> p;
  while ( (c = getc(input)) != EOF )
    {
      ungetc(c,input);

      while ((c = getc(input)) != '\n' )
	{
	  while (c != '0' && c != '1' && c != '2' && c != '3'
		 && c != '4' && c != '5' && c != '6' && c != '7'
		 && c != '8' && c != '9' && c != '\n' && c != EOF && c != '-')
	    c = getc(input);
	  if (c != '\n' && c != EOF) {
	    ungetc(c,input);
	    float f;
	    fscanf(input, "%f",&f);
	    push(p,f);
	  }
	  else
	    if (c == '\n')
	      ungetc(c,input);
	}

      if (p.index %8 > 0)
	for (int i = 8 - p.index %8; i> 0; i--)
	  push(p,(float) 0.);
      float *new_p;
      posix_memalign((void **)&new_p, 16, p.index*sizeof(float));
      memcpy(new_p,p.elements,sizeof(float)*p.index);

      if (point_len > 0 && point_len != p.index)
	{
	  printf("Can't handle vectors of differing length, bailing\n");
	  exit(0);
	}

      point_len = p.index;
      p.index = 0;
      push(parsed,new_p);
    }
  return parsed;
}

void print(point &p)
{
  for (int i = 0; i<point_len; i++)
    printf("%f ",p[i]);
  printf("\n");
}


v_array<point > parse_points(float *input, int data_size, int dim)
{
  v_array<point > parsed;
  char c;
  v_array<float> p;

  for (int i=0;i<data_size;i++)
  {
      for (int j=0;j<dim;j++){
         push(p,*(input+i*dim+j));
      }
      if (p.index %8 > 0){
         for (int i = 8 - p.index %8; i> 0; i--)
            push(p,(float) 0.);
      }

      float *new_p;
      posix_memalign((void **)&new_p, 16, p.index*sizeof(float));
      memcpy(new_p,p.elements,sizeof(float)*p.index);

      if (point_len > 0 && point_len != p.index)
      {
          printf("Can't handle vectors of differing length, bailing\n");
          exit(0);
      }

      point_len = p.index;
      p.index = 0;
      push(parsed,new_p);
    }
  return parsed;
}

