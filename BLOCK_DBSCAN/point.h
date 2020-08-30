/*
 *
 * Copyright (C) 2017-  Yewang Chen<ywchen@hqu.edu.cn;nalandoo@gmail.com;1611414017@hqu.edu.cn>
 * License: GPL v1
 * This software may be modified and distributed under the terms
 * of license.
 *
 */
#include "stack1.h"
#include <stdio.h>
#include <stdlib.h>
#include "DataPoint.h"
#define NDEBUG
#include<assert.h>
#include <string.h>
#include <iostream>
#include <errno.h>
using namespace std;

typedef float* point;


//float complete_distance(float* v1, float* v2);
float distance_1(float* v1, float* v2);
v_array<point > parse_points(FILE *input);
void print(point &p);
int posix_memalign(void **mptr, size_t, size_t bytes);
v_array<point > parse_points(float *input, int data_size, int dim);
