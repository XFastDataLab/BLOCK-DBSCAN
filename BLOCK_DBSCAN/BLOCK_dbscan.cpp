/*
 *
 * Copyright (C) 2017-  Yewang Chen<ywchen@hqu.edu.cn;nalandoo@gmail.com;1611414017@hqu.edu.cn>
 * License: GPL v1
 * This software may be modified and distributed under the terms
 * of license.
 *
 */
#include "ClusterAnalysis.h"
#include <cstdio>
#include <ctime>
#include <cstdlib>
#include "cyw_timer.h"




int main()
{
	ClusterAnalysis myClusterAnalysis;                        //Clustering algorithm object declaration.
	cout<<"loading the file..."<<endl;
	myClusterAnalysis.Init("data\\points_new.txt",0.05,55);        //Algorithm initialization eps=0.05 minPts=55.
	cout<<"clusting the data..."<<endl;
    CYW_TIMER build_timer;
    build_timer.start_my_timer();
	myClusterAnalysis.DoDBSCANRecursive();                    //Perform C1_DBSCAN.
	build_timer.stop_my_timer();
    printf("the running time is %f\n",build_timer.get_my_timer());
    myClusterAnalysis.WriteToFile("data\\result.txt");      //Save the result.
	system("pause");    //print the result

	return 0;
}
