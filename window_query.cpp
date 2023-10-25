#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <typeinfo> 
#include <tuple>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <map>
#include <set>
#include "Constants.h"
#include "utils.h"

using namespace std;


vector<int> get_cluster_candidate(double *ql, double *qh, double * grid_info, int ** grid_cluster, int * grid_block_nums, int &scale)
{
    int ql_xid = ((ql[0] - grid_info[4])<0?0:(ql[0] - grid_info[4])) / grid_info[2];
    int ql_yid = ((ql[1] - grid_info[5])<0?0:(ql[1] - grid_info[5])) / grid_info[3];
    int qh_xid = ((qh[0] - grid_info[4])<grid_info[8]?(qh[0] - grid_info[4]):grid_info[6]) / grid_info[2];
    int qh_yid = ((qh[1] - grid_info[5])<grid_info[9]?(qh[1] - grid_info[5]):grid_info[7]) / grid_info[3];
    int y_grid_num = grid_info[1];
    
    int count = 0;
    vector<int> cids;
    int i,j,m,gid;
    
    for (i = ql_xid; i < qh_xid+1; i++)
    {
        for (j = ql_yid; j < qh_yid + 1; j++)
        {
            gid = i*y_grid_num + j;
            for (m = 0; m < grid_block_nums[gid]; m++) 
            {
                cids.push_back(grid_cluster[gid][m]);
            }
        }  
    }

    vector<int> distinct_cids = remove_duplicates(cids);
    return distinct_cids;
}

tuple<int,int> convert2multiid(int cluster_long_id, int scale)
{
    return {cluster_long_id%scale, (cluster_long_id / scale) % scale};
}

int classify_cluster(double *ql, double *qh, int pid, int k, double *center, double R, double * qc, double qr)
{
    double dist = 0;
    for(int d = 0; d<2;d++)
    {
        if (center[d]>qh[d])
        {
            dist += pow((center[d]-qh[d]),2);
        }    
        else if(center[d]<ql[d])
        {
            dist += pow((center[d]-ql[d]),2);
        } 
    }
    if(pow(dist,0.5)>R)
    {
        return -1; 
    }
    else
    {
        if((center[0]-R >=ql[0])&&(center[1]-R>=ql[1])&&(center[1]+R<=qh[1])&&(center[0]+R<=qh[0]))
        {
            return 1;
        }
        else
        {
            return 0;
        }
    }
}

tuple<vector<double *>, int> range_query(double *query, double * grid_info,int ** grid_cluster, int * grid_block_nums,double*** centers,double** R,int ** cluster_block_nums,double *** split_points, int *** block_ids, double *** blocks, int * block_count, int &B, int &scale)
{
    double ql[2] = {query[0],query[1]};
    double qh[2] = {query[2],query[3]};
    double query_center[2] = {(ql[0]+qh[0])*0.5, (ql[1]+qh[1])*0.5};
    double query_r = cal_dist2(query_center,ql,2);

    vector<int> cid = get_cluster_candidate(ql, qh, grid_info, grid_cluster, grid_block_nums,scale);

    int block_visit = 0;
    int visit_point = 0;
    int block_cur_count,m,j,cluster_block_num;
    double * split_point;
    double ** select_block;
    int cid_num = cid.size();
    tuple<int,int> mids;
    int pid;
    int k;
    int tag;
    // vector<double *> result;
    int * select_bid;
    int * pid_list;
    int * k_list;
    pid_list = new int[cid_num];
    k_list = new int[cid_num];
    int total_blk = 0;
    for (int i = 0; i<cid_num; i++)
    {
        mids = convert2multiid(cid[i],scale);
        pid_list[i] = get<0>(mids);
        k_list[i] = get<1>(mids);
    }

    vector<double *> result;
    int result_len = 0;

    for (int i = 0; i<cid_num; i++)
    {
        pid = pid_list[i];
        k = k_list[i];
        double center[2] = {centers[pid][k][0],centers[pid][k][1]};
        tag = classify_cluster(ql, qh, pid, k, center, R[pid][k], query_center, query_r);
        select_bid = block_ids[pid][k];
        if(tag == -1){continue;}
        else if(tag == 1)
        {
            for (j = 0; j<cluster_block_nums[pid][k]; j++)
            {
                select_block = blocks[select_bid[j]];
                block_cur_count = block_count[select_bid[j]];
                block_visit++;
                for (m = 0; m<block_cur_count; m++)
                {
                    result.push_back(select_block[m]);
                }
            }
            continue;
        }
        else
        {
            split_point = split_points[pid][k];
            cluster_block_num = cluster_block_nums[pid][k];
            auto end_loc = lower_bound(split_point, split_point+cluster_block_num, qh[0]) - split_point;
            auto start_loc = (upper_bound(split_point, split_point+cluster_block_num, ql[0]) - split_point) -1;
            start_loc = start_loc<0?0:start_loc;
            for (j = start_loc; j<end_loc; j++)
            {
                select_block = blocks[select_bid[j]];
                block_cur_count = block_count[select_bid[j]];
                block_visit++;
                for (m = 0; m<block_cur_count; m++)
                {
                    if(select_block[m][0]>qh[0]){break;}
                    double point[2] = {select_block[m][0],select_block[m][1]};
                    if ((point[0]>=ql[0])&&(point[1]>=ql[1])&&(point[1]<=qh[1]))
                    {
                        result.push_back(select_block[m]);
                    }         
                }
            }

        }
    }

    return {result,block_visit};
}

int main(){
    sync();
    string dataset = Constants::DATASETS;
    string model_record_path = Constants::MODEL_R_PATH;

    int dim = Constants::DIM;
    int B = Constants::B;
    int scale = Constants::SCALE;
    string distribution = Constants::DISTRIBUTION;

    int level_num = Constants::LEVEL;
    int level_cluster[level_num] = {20,20,20,20};
    int partition = 1;
    for(int i = 0; i<level_num-1; i++)
    {
        partition *= level_cluster[i];
    }

    int data_size = getFileLine(Constants::INDEX_PATH);
    double ** dat =  get_2d_points_csv(Constants::INDEX_PATH,data_size);

    int query_num = getFileLine(Constants::WINDOW_WORKLOAD_PATH + "window_"+distribution+"_1.0_0.0001_.csv");
    double ** query = get_2d_points_csv(Constants::WINDOW_WORKLOAD_PATH + "window_"+distribution+"_1.0_0.0001_.csv",query_num);
    tuple<int *,int> block_info = get_block_count(Constants::BLOCK_INFO);
    int * block_count = get<0>(block_info);
    int total_block_num = get<1>(block_info);

    double *** blocks = get_blocks(dat,block_count,total_block_num,B);
    tuple<int *, int> base_info = get_block_count(Constants::WINDOW_WORKLOAD_PATH + "window_result_"+distribution+"_"+to_string(data_size)+"_1.0_0.0001_.csv");
    
    int * base_len = get<0>(base_info);

    tuple<int ***, int **> cluter_block = get_block_ids(partition, model_record_path, Constants::CLUSTER_BLOCK);
    int *** block_ids = get<0>(cluter_block);
    int ** cluster_block_nums = get<1>(cluter_block);
    double *** split_points = get_split_points(partition, model_record_path, Constants::SPLIT_PTS);

    tuple<double ***, double**, int*> cluster_infos = get_cluster_infos(partition, model_record_path, Constants::CLUSTER_INFO);
    double*** centers;
    double** R;
    int *sub_cluster_num;
    centers = get<0>(cluster_infos);
    R = get<1>(cluster_infos);
    sub_cluster_num = get<2>(cluster_infos);
    
    int grid_num = getFileLine(Constants::GRID_CLUSTER);
    tuple<int **, int*> grid_cluster_info = get_cluster_id(Constants::GRID_CLUSTER, grid_num);
    int ** grid_cluster = get<0>(grid_cluster_info);
    int * grid_block_nums  = get<1>(grid_cluster_info);

    int grid_param_num = getFileLine(Constants::GRID_INFO);
    double * grid_info = get_grid_info(Constants::GRID_INFO,grid_param_num);

    int prediction=0;
    long query_time=0;
    long long s_time = 0;
    double recall = 0;
    double total_block_access = 0;
    int result_len,block_access;
    tuple<vector<double *>, int> query_result;
    std::cout<<"query num: "<<query_num<<endl;
    for (int n = 0; n< query_num; n++)
    {
        auto start = chrono::high_resolution_clock::now();
        query_result = range_query(query[n], grid_info, grid_cluster, grid_block_nums, centers, R, cluster_block_nums,split_points, block_ids, blocks, block_count,B, scale);
        auto end = chrono::high_resolution_clock::now();
        int result_len = get<0>(query_result).size();
        int block_access = get<1>(query_result);
        recall += (double)result_len/(double)base_len[n];
        total_block_access += block_access;
        query_time += chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    }
    std::cout<<"time: "<<(query_time/query_num)<<endl;
    std::cout<<"recall: "<<(recall/query_num)<<endl;
    std::cout<<"access: "<<((double)total_block_access/(double)query_num)<<endl;
    return 0;
}

