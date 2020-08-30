/*
 *
 * Copyright (C) 2017-  Yewang Chen<ywchen@hqu.edu.cn;nalandoo@gmail.com;1611414017@hqu.edu.cn>
 * License: GPL v1
 * This software may be modified and distributed under the terms
 * of license.
 *
 */
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iosfwd>
#include "cover_tree.h"
#include <windows.h>
#include "cyw_timer.h"
using namespace std;

//Cluster analysis type.
class ClusterAnalysis
{
private:
	vector<DataPoint> dadaSets;        //the dataset.
	vector< vector< int > > TIC;
	vector< vector< int > > OCS;
	vector< vector< int > > BOR;
    int* types;
    int* scanned;
    int* out_core_id;
    v_array<point> data_set;
	unsigned int dimNum;
	float radius;
	unsigned int dataNum;
	unsigned int minPTs;
    int iti;
    int ito;
    int clusterId;
    node node_data;
public:
    void SetArrivalPoints(DataPoint& dp);
	ClusterAnalysis(){}                //Default constructor.
	float distance_2(point p1, point p2);
	bool DoDBSCANRecursive();
	bool Init(char* fileName, double radius, int minPTs);    //Initialization operation
	bool WriteToFile(char* fileName);    //save results
	bool Judge_density_reachable(vector<int>,vector<int>);
	void clusterId_adjust();
	void Handle_Inner_Cores();
    void Handle_Outer_Cores();
    void Retrieve_Borders();
    void Merge(int a,int b);
    int get_dim(char* s, char* delims);
    float*  get_data(char* s, int dim,char* delims);
    void read_data_dim_size(char* filename, int* data_dim, int* data_size, char* delims);
    float* read_data(char* filename, char* delims);
    float* read_data(char* filename, char* delims, int* dim, int* data_size);
};
