/*
 *
 * Copyright (C) 2017-  Yewang Chen<ywchen@hqu.edu.cn;nalandoo@gmail.com;1611414017@hqu.edu.cn>
 * License: GPL v1
 * This software may be modified and distributed under the terms
 * of license.
 *
 */
#include "ClusterAnalysis.h"
#include <algorithm>


bool ClusterAnalysis::Init(char* fileName, double radius, int minPTs){
    this->radius = radius;        //set the radius
	this->minPTs = minPTs;        //set the minPTs
	this->dimNum = DIME_NUM;    //set the dim of the dataset
    int dim;
    int data_size;
    std::cout<<"reading data...\n";

    float *raw_data=read_data(fileName," ",&dim, &data_size);
    dataNum = data_size;
    types = (int*)malloc(sizeof(int)*data_size);
    scanned = (int*)malloc(sizeof(int)*data_size);
    out_core_id = (int*)malloc(sizeof(int)*data_size);
    for(int i  = 0; i <dataNum;i++){
	    out_core_id[i] = 0;
	    scanned[i]= 0;
	    types[i] = 0;
		DataPoint tempDP;
		double tempDimData[DIME_NUM];
		for (int j = 0; j<DIME_NUM; j++)
		{
			tempDimData[j] = raw_data[i*DIME_NUM+j];
		}
		tempDP.SetDimension(tempDimData);
		tempDP.SetDpId(i);
		tempDP.SetVisited(false);
		tempDP.SetClusterId(-1);
		dadaSets.push_back(tempDP);
	}
    std::cout<<"push to array...\n";
    data_set= parse_points(raw_data,dataNum,DIME_NUM);
    free(raw_data);
    cout<<"dataNum is: "<<dataNum<<endl;;

    CYW_TIMER build_timer;
    build_timer.start_my_timer();
    std::cout<<"building the trees...\n";
    node_data = batch_create(data_set);
    build_timer.stop_my_timer();

    printf("the tree building time is %f\n",build_timer.get_my_timer());
    return true;
}

bool ClusterAnalysis::WriteToFile(char* fileName){
	ofstream of1(fileName);
	for (unsigned long i = 0; i<dataNum; i++)
	{
        of1 << i << '\t';
		of1 << dadaSets[i].GetClusterId()<< '\t'<<endl;
	}
	of1.close();
	return true;
}

void ClusterAnalysis::SetArrivalPoints(DataPoint& dp){
    int query_num = 1; //set the query number as 1
    float *query_data=(float*)malloc(sizeof(float)*DIME_NUM);
	for(int i = 0;i<DIME_NUM;i++){
        query_data[i]  = dp.GetDimension()[i];
	}
    v_array<point> queries=parse_points(query_data,1,DIME_NUM);
    free(query_data);
    node node_query = batch_create(queries);
    v_array<v_array<point> > res;
    scanned[dp.GetDpId()] = 1;
    epsilon_nearest_neighbor(node_data, node_query,res, radius);
    int cout1 = res[0].index-1;
    int cout2 = 0;
    if(cout1 < minPTs)
        return;
    vector <int> tic;
    vector <int> ocs;
    vector <int> bor;
    for(int i = 0; i<res[0].index;i++){
        if(distance_1(queries[0],res[0][i])<= radius/2)
        {
            cout2++;
        }
    }
    if(cout1 > minPTs ){
        if( cout2>= minPTs )
        {
            tic.push_back(dp.GetDpId());
            types[dp.GetDpId()] = 2;
            for(int i = 0; i<res[0].index;i++)
            {
                if(distance_1(queries[0],res[0][i])<= radius/2)
                {
                    int id =res[0][i][DIME_NUM-1];
                    scanned[id] = 1;
                    if(types[id] != 1)
                        types[id] = 2;
                    tic.push_back(id);
                }
            }
            TIC.push_back(tic);
            tic.clear();
        }
        else
        {
            ocs.push_back(dp.GetDpId());
            types[dp.GetDpId()] = 1;
            for(int i = 0; i<res[0].index;i++)
            {
                int id =res[0][i][DIME_NUM-1];
                if(types[id] == 0)
                    types[id] = 3;
                ocs.push_back(id);
            }
            OCS.push_back(ocs);
            ocs.clear();
        }
    }
    else{
        types[dp.GetDpId()] = 3;
        bor.push_back(dp.GetDpId());
        for(int i = 0; i<res[0].index;i++)
        {
            int id =res[0][i][DIME_NUM-1];
            bor.push_back(id);
        }
        BOR.push_back(bor);
        bor.clear();
    }
}

bool ClusterAnalysis::DoDBSCANRecursive(){
//    int it1=0;
//    int it2=0;
//    int it3=0;
//    int it4=0;
    for (unsigned long i = 0; i<dataNum; i++){
	    if(scanned[i] == 0)
        {
//	        it1++;
            SetArrivalPoints(dadaSets[i]);
        }

	}
	cout<<"Handle Inner Cores..."<<endl;
    Handle_Inner_Cores();
    cout<<"Handle Outer Cores..."<<endl;
    Handle_Outer_Cores();
    cout<<"retrieve borders..."<<endl;
    Retrieve_Borders();
	cout << "the result contains " << clusterId << " clusters" << endl;
//    cout<<"the number of RangeQuery is "<<it1<<endl;
//    cout<<"the number of inner core blocks is "<<TIC.size()<<endl;
//    for(int i = 0;i<dataNum;i++){
//        if(types[i] == 2)
//            it2++;
//        if(types[i] == 1)
//            it3++;
//        if(types[i] == 0)
//            it4++;
//
//    }
//    cout<<"the number of inner-core points is "<<it2<<endl;
//    cout<<"the number of outer core points is "<<OCS.size()<<endl;
//    cout<<"the number of filtered points is "<<dataNum - it1<<endl;
//    cout<<"the number of noise points is "<<it4<<endl;
//    cout<<"the number of core points is "<<it2+it3<<endl;
//    ofstream of1("data\\center_points.txt");
//	 for (unsigned long i = 0; i<TIC.size(); i++){
//            for (int d = 0; d<DIME_NUM; d++)
//                of1 << dadaSets[TIC[i][0]].GetDimension()[d] << '\t';
//            of1<<types[TIC[i][0]]<<endl;
//	 }
//	of1.close();
//
//  ofstream of2("data\\TIC.txt");
//	for (unsigned long i = 0; i<TIC.size(); i++)
//	{
//        for(int j = 0;j<TIC[i].size();j++)
//        {
//            for (int d = 0; d<DIME_NUM; d++)
//                of2 << dadaSets[TIC[i][j]].GetDimension()[d] << '\t';
//            of2<<types[TIC[i][j]]<<endl;
//        }
//	}
//	of2.close();
//
//	ofstream of3("data\\OCS_center_points.txt");
//	for (unsigned long i = 0; i<OCS.size(); i++)
//	{
//        for (int d = 0; d<DIME_NUM; d++)
//            of3 << dadaSets[OCS[i][0]].GetDimension()[d] << '\t';
//        of3<<types[OCS[i][0]]<<endl;
//	}
//	of3.close();
//
//	ofstream of4("data\\Borders.txt");
//	for (unsigned long i = 0; i<dataNum; i++)
//	{
//	    if(types[i] == 3)
//	    {
//        for (int d = 0; d<DIME_NUM; d++)
//            of4 << dadaSets[i].GetDimension()[d] << '\t';
//        of4<<types[i]<<endl;
//	    }
//	}
//	of4.close();

//    ofstream of5("data\\After_merge.txt");
// 	for (int i = 0; i<dataNum; i++)
//	{
//	    if(types[i] == 1 ||types[i]==2)
//	    {
//        for (int d = 0; d<DIME_NUM; d++)
//            of5 << dadaSets[i].GetDimension()[d] << '\t';
//        of5<<dadaSets[i].GetClusterId()<<endl;
//	    }
//	}
//	of5.close();
	return true;
}

void ClusterAnalysis::Handle_Inner_Cores(){
    for(int i = 0;i <TIC.size();i++){
        if(dadaSets[TIC[i][0]].GetClusterId() == -1)
        {
            clusterId++;
            for(int m = 0; m <TIC[i].size();m++)
            {
                dadaSets[TIC[i][m]].SetClusterId(clusterId);
            }
            for(int j = i+1; j <TIC.size();j++)
            {
                if(dadaSets[TIC[j][0]].GetClusterId() == -1)
                {
                    float d = distance_1(data_set[TIC[i][0]],data_set[TIC[j][0]]);
                    if(d <= radius)
                    {
                        for(int m = 0; m <TIC[j].size();m++)
                        {
                            dadaSets[TIC[j][m]].SetClusterId(clusterId);
                        }

                    }
                    else if(d > radius && d < 2*radius)
                    {
                        if(Judge_density_reachable(TIC[i],TIC[j]) == true)
                        {
                            for(int m = 0; m <TIC[j].size();m++)
                            {
                                dadaSets[TIC[j][m]].SetClusterId(clusterId);
                            }
                        }
                    }
                }
            }
            for(int j = i+1; j <TIC.size();j++){
                int init_id = dadaSets[TIC[i][0]].GetClusterId();
                if(dadaSets[TIC[j][0]].GetClusterId() > -1 && dadaSets[TIC[j][0]].GetClusterId()!=init_id)
                {
                    float d = distance_1(data_set[TIC[i][0]],data_set[TIC[j][0]]);
                    if(d <= radius)
                    {
                        int temp_id = dadaSets[TIC[j][0]].GetClusterId();
                        for(int m = 0; m <TIC.size();m++)
                        {
                            if(dadaSets[TIC[m][0]].GetClusterId() == init_id)
                            {
                                for(int n = 0;n < TIC[m].size();n++)
                                {
                                    dadaSets[TIC[m][n]].SetClusterId(temp_id);
                                }
                            }
                        }
                        clusterId--;
                    }
                    else if(d > radius && d < 2*radius)
                    {
                        if(Judge_density_reachable(TIC[i],TIC[j]) == true)
                        {
                            int temp_id = dadaSets[TIC[j][0]].GetClusterId();
                            for(int m = 0; m <TIC.size();m++)
                            {
                                if(dadaSets[TIC[m][0]].GetClusterId() == init_id)
                                {
                                    for(int n = 0;n < TIC[m].size();n++)
                                    {
                                        dadaSets[TIC[m][n]].SetClusterId(temp_id);
                                    }
                                }
                            }
                        clusterId--;
                        }
                    }
                }
            }
        }
    }
}

bool ClusterAnalysis::Judge_density_reachable(vector<int> a,vector<int> b){
    if(a.size()>10*minPTs || b.size()>10*minPTs ){
        bool flag = true;
        int m = 0;
        int n = 0;
        float mindist1 = distance_1(data_set[a[0]] , data_set[b[0]]);
        float mindist2 = mindist1;
        for(int i = 0;i<10;i++)
        {
            Loop_Judge:;
            for(int i = 1; i<a.size();i++)
            {
                float dist = distance_1(data_set[b[n]] , data_set[a[i]]);
                if(dist < mindist1)
                {
                    m = i;
                    mindist1 = dist;
                }
            }
            int temp = n;
            for(int i = 1; i<b.size();i++)
            {
                float dist = distance_1(data_set[a[m]] , data_set[b[i]]);
                if(dist < mindist1)
                {
                    n = i;
                    mindist2 = dist;
                }
            }
            if(n == temp)
            {
                float mindistance = distance_1(data_set[a[m]] , data_set[b[n]]);
                if(mindistance <= radius)
                    return true;
                else
                    return false;
            }
            else
            {
                goto Loop_Judge;
            }
        }
        float mindistance = distance_1(data_set[a[m]] , data_set[b[n]]);
        if(mindistance <= radius)
            return true;
        else
            return false;
    }
    else{
        for(int i = 1;i<a.size();i++)
        {
            for(int j =1;j<b.size();j++)
            {
                if(distance_1(data_set[a[i]] , data_set[b[j]]) <= radius)
                {
                    return true;
                }
            }
        }
        return false;
    }
}

void ClusterAnalysis::Retrieve_Borders(){
    for(int i = 0; i<BOR.size();i++){
        for(int j = 0; j<BOR[i].size();j++){
            if(types[BOR[i][j]] == 2 || types[BOR[i][j]] == 1)
            {
                int id = dadaSets[BOR[i][j]].GetClusterId();
                dadaSets[BOR[i][0]].SetClusterId(id);
                break;
            }
        }
    }
    int it = 0;
    pair<int,int> a[10000];
    for(int i = 0;i<dataNum;i++){
        if (dadaSets[i].GetClusterId()!=-1)
        {
            for(int j = 0;j <10000;j++)
            {
                if(a[j].first == dadaSets[i].GetClusterId())
                {
                    dadaSets[i].SetClusterId(a[j].second);
                    goto LoopRetrieve;
                }
            }
            it++;
            a[it].first = dadaSets[i].GetClusterId();
            a[it].second = it;
            dadaSets[i].SetClusterId(a[it].second);
        }
        LoopRetrieve:;
    }
    clusterId = it;
}

void ClusterAnalysis::Handle_Outer_Cores(){
    for(int i = 0;i<OCS.size();i++){
        for(int j = 0;j<TIC.size();j++)
        {
            float d = distance_1(data_set[OCS[i][0]],data_set[TIC[j][0]]);
            if(d<0)
                continue;
            else
                break;
        }
    }
    for(int i = 0;i<OCS.size();i++){
        for(int j = 0;j<OCS[i].size();j++)
        {
            if(types[OCS[i][j]] == 2)
            {
                ito++;
                int c_id = dadaSets[OCS[i][j]].GetClusterId();
                for(int m =0;m<OCS[i].size();m++)
                {
                    dadaSets[OCS[i][m]].SetClusterId(c_id);
                    out_core_id[OCS[i][m]] = 1;
                }
                break;
            }
        }
    }
    for(int i = 0; i<OCS.size();i++){
        if(out_core_id[OCS[i][0]] == 0)
        {
            clusterId++;
            for(int k = 0;k<OCS[i].size();k++)
            {
                out_core_id[OCS[i][k]] = 1;
                dadaSets[OCS[i][k]].SetClusterId(clusterId);
            }
            for(int j = 0; j<OCS[i].size();j++)
            {
                if(types[OCS[i][j]] == 1 && dadaSets[OCS[i][j]].GetClusterId() !=-1)
                {
                    ito++;
                    int id1 = dadaSets[OCS[i][j]].GetClusterId();
                    for(int k = 0;k<OCS[i].size();k++)
                    {
                        dadaSets[OCS[i][k]].SetClusterId(id1);
                    }
                    goto Loop1;
                }
            }Loop1:;
        }
        else
        {
            int id = dadaSets[OCS[i][0]].GetClusterId();
            for(int k = 1;k<OCS[i].size();k++)
            {
                out_core_id[OCS[i][k]] = 1;
                dadaSets[OCS[i][k]].SetClusterId(id);
            }
        }
    }
}

int ClusterAnalysis::get_dim(char* s, char* delims){
    char *val_str = NULL;
    val_str = strtok( s, delims );
    int dim=0;
    while( val_str != NULL ) {
        dim++;
        val_str = strtok( NULL, delims );
    }
    return dim;
}

float*  ClusterAnalysis::get_data(char* s, int dim,char* delims){
    float* temp= (float*) malloc (dim*sizeof(float));
    char *val_str = NULL;
    val_str = strtok( s, delims );
    int counter=0;
    while( val_str != NULL ) {
        temp[counter]=atof(val_str);
        counter++;
        val_str = strtok( NULL, delims );
    }
    return temp;
}

void ClusterAnalysis::read_data_dim_size(char* filename, int* data_dim, int* data_size, char* delims){
    int n_size=0;
    int dim=0;
    char s[10000];
    freopen(filename,"r",stdin);
    while(gets(s))
    {
        if (dim==0)
           dim=get_dim(s,delims);
        n_size ++;
    }
    *data_dim=dim;
    *data_size=n_size;
    fclose(stdin);
}

float* ClusterAnalysis::read_data(char* filename, char* delims){
    int dim, n_size;
    read_data_dim_size(filename,&dim, &n_size, delims);

    float* data= (float*) malloc (n_size*dim*sizeof(float));
    freopen(filename,"r",stdin);
    int counter=0;
    char s[10000];
    while(gets(s))
    {
        float* tmp_data= get_data( s, dim,delims);
        memcpy(data+counter*dim,tmp_data,dim*sizeof(float));
        counter++;
        free(tmp_data);
    }
    fclose(stdin);

    return data;
}

float* ClusterAnalysis::read_data(char* filename, char* delims, int* dim, int* data_size){
    read_data_dim_size(filename,dim, data_size, delims);

    float* data= (float*) malloc ((*data_size)*(*dim)*sizeof(float));
    freopen(filename,"r",stdin);
    int counter=0;
    char s[10000];
    while(gets(s))
    {
        float* tmp_data= get_data( s,*dim,delims);
        memcpy(data+counter*(*dim),tmp_data,(*dim)*sizeof(float));
        counter++;
        free(tmp_data);
    }
    fclose(stdin);

    return data;
}
