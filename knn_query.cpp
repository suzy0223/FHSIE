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
#include <numeric> 
#include "Constants.h"
#include "utils.h"
  
using namespace std;




double cal_dist(double *point1, double *point2)
{
    double temp_dist = sqrt(pow((point1[0] - point2[0]), 2) + pow((point1[1] - point2[1]), 2));
    return temp_dist;
}

tuple<int **,int> get_cluster_candidate(double *query, double search_range, double * grid_info, int ** grid_cluster, int * grid_block_nums, int &scale)
{
    int ql_xid = ((query[0] - grid_info[4]-search_range)<0?0:(query[0] - grid_info[4]-search_range)) / grid_info[2];
    int ql_yid = ((query[1] - grid_info[5]-search_range)<0?0:(query[1] - grid_info[5]-search_range)) / grid_info[3];
    int qh_xid = ((query[0] - grid_info[4]+search_range)<(grid_info[6] - grid_info[4])?(query[0] - grid_info[4] + search_range):(grid_info[6] - grid_info[4])) / grid_info[2];
    int qh_yid = ((query[1] - grid_info[5]+search_range)<(grid_info[7] - grid_info[5])?(query[1] - grid_info[5] + search_range):(grid_info[7] - grid_info[5])) / grid_info[3];
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
    count = distinct_cids.size();
    int ** cluster_candidate = convert2multiclusterid(distinct_cids,count,scale);
    return {cluster_candidate,count};
}

tuple<int **,int> filter_clster(int ** cid, double *** centers, double ** R, double * query_center, double query_r,int &cid_num)
{
    int i_list[cid_num];
    int count = 0;
    for (int i = 0; i < cid_num; i++)
    {
        int tmp_pid = cid[i][0];
        int tmp_k = cid[i][1];
        // cout<<"i ="<<i<<"pid = "<<tmp_pid<<"k = "<<tmp_k<<endl;
        double center[2] = {centers[tmp_pid][tmp_k][0],centers[tmp_pid][tmp_k][1]};
        double r = R[tmp_pid][tmp_k];
        double dist = 0;
        for(int d = 0; d<2;d++)
        {
            if (center[d]>query_center[d]+query_r)
            {
                dist += pow((center[d]-(query_center[d]+query_r)),2);
            }    
            else if(center[d]<query_center[d]-query_r)
            {
                dist += pow((center[d]-(query_center[d]-query_r)),2);
            } 
        }
        if(pow(dist,0.5)<=r)
        {
            i_list[count] = i;
            count++;
        }
    }

    int ** ids;
    ids = new int*[count];
    for (int i = 0; i < count; i++)
    {
        ids[i] = cid[i_list[i]];
    }

    return {ids,count};
}

tuple<double **, int, int> range_query(double *query, double * grid_info,int ** grid_cluster, int * grid_block_nums,double*** centers,double** R,int ** cluster_block_nums,double *** split_points, int *** block_ids, double *** blocks, int * block_count,double &search_range, int &B, int &scale)
{ 
    tuple<int**,int> cid_info = get_cluster_candidate(query, search_range, grid_info, grid_cluster, grid_block_nums,scale);
    // double search_r = search_range * sqrt(2);
    cid_info = filter_clster(get<0>(cid_info), centers, R, query, search_range, get<1>(cid_info));

    int ** cid_selected = get<0>(cid_info);
    int cid_selected_num = get<1>(cid_info);

    int access_block = 0;
    for (int i = 0; i<cid_selected_num; i++)
    {
        int tmp_pid = cid_selected[i][0];
        int tmp_k = cid_selected[i][1];
        access_block += cluster_block_nums[tmp_pid][tmp_k];
    }

    double **result;
    result = new double*[access_block*B];
    int result_len = 0;
    int block_visit = 0;
    for (int i = 0; i<cid_selected_num; i++)
    {
        int tmp_pid = cid_selected[i][0];
        int tmp_k = cid_selected[i][1];
        double * split_point = split_points[tmp_pid][tmp_k];
        int cluster_block_num = cluster_block_nums[tmp_pid][tmp_k];
        auto end_loc = lower_bound(split_point, split_point+cluster_block_num, query[0]+search_range) - split_point;
        auto start_loc = (upper_bound(split_point, split_point+cluster_block_num, query[0]-search_range) - split_point) -1;
        start_loc = start_loc<0?0:start_loc;
        int * select_bid = block_ids[tmp_pid][tmp_k];

        for (int j = start_loc; j<end_loc; j++)
        {
            double ** select_block = blocks[select_bid[j]];
            int block_cur_count = block_count[select_bid[j]];

            block_visit++;
            for (int m = 0; m<block_cur_count; m++)
            {
                if((select_block[m][0]>(query[0]+search_range))){break;}
                double tmp_dist = cal_dist(select_block[m],query);
                if(tmp_dist<=search_range)
                {
                    result[result_len] = select_block[m];
                    result_len++;
                }
            }
        }
    }

    return {result,result_len, block_visit};
}

double cal_recall(double ** baseline, double ** predict, int &k, int &dim)
{
    double meet = 0;
    
    for(int i=0; i<k; i++)
    {
        double tmp_dist = sqrt(cal_dist2(predict[i],baseline[i],dim));
        if (tmp_dist<0.00000000000005){meet++;}
    }

    double recall = meet/k;
    
    return recall;
}

int main(){
    string dataset = Constants::DATASETS;
    string model_record_path = Constants::MODEL_R_PATH;
    int dim = Constants::DIM;
    int B = Constants::B;
    int scale = Constants::SCALE;
    string distribution = Constants::DISTRIBUTION;

    int level_num = Constants::LEVEL;
    int level_cluster[level_num] = {27,27,27};

    int partition = 1;
    cout<<"getting datasets and workload..."<<endl;

    int data_size = getFileLine(Constants::INDEX_PATH);
    double ** dat =  get_2d_points_csv(Constants::INDEX_PATH,data_size);
    
    int query_num = getFileLine(Constants::KNN_WORKLOAD_PATH + "knn_" + distribution+".csv");
    cout<<query_num<<endl;
    double ** query = get_2d_points_csv(Constants::KNN_WORKLOAD_PATH + "knn_" + distribution+".csv", query_num);
    

    cout<<"getting blocks..."<<endl;
    tuple<int *,int> block_info = get_block_count(Constants::BLOCK_INFO);
    int * block_count = get<0>(block_info);
    int total_block_num = get<1>(block_info);
    double *** blocks = get_blocks(dat,block_count,total_block_num,B);

    cout<<"getting model..."<<endl;
    double *** model;
    model = new double**[level_num-1];
    for(int i = 0; i<level_num-1; i++)
    {
        partition *= level_cluster[i];
        model[i] = get_2d_points_csv(Constants::K_MEANS_PATH + "/level"+to_string(i)+"model.csv",partition);
    }
    
    cout<<"getting parts and clusters info..."<<endl;
    tuple<int *,int> parts_info = get_block_count(Constants::PART_INFO);
    int * pid = get<0>(parts_info);
    int total_parts_num = get<1>(parts_info);
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

    cout<<"getting grid info..."<<endl;
    int grid_num = getFileLine(Constants::GRID_CLUSTER);
    tuple<int **, int*> grid_cluster_info = get_cluster_id(Constants::GRID_CLUSTER, grid_num);
    int ** grid_cluster = get<0>(grid_cluster_info);
    int * grid_block_nums  = get<1>(grid_cluster_info);

    int grid_param_num = getFileLine(Constants::GRID_INFO);
    double * grid_info = get_grid_info(Constants::GRID_INFO,grid_param_num);

    cout<<"getting knn base..."<<endl;
    double *** base = knn_base_result(Constants::KNN_WORKLOAD_PATH + "knn_result_" + distribution+ "_"+to_string(data_size)+".csv");
    
    cout<<"start query now"<<endl;
    int K_nums_list[5] = {1, 5, 25, 50, 75};
    for (int ks = 0; ks<5; ks++)
    {
        long query_time=0;
        double recall_record = 0.0;
        int prediction=0;
        double total_block_access = 0;
        int knn_k = K_nums_list[ks];
        cout<<"knn query k = "<<knn_k<<endl;
        int next_idx,mid_begin,mid_end;
        for (int i = 0; i< query_num; i++)
        {
        auto start = chrono::high_resolution_clock::now();
        double cur_query[2] = {query[i][0],query[i][1]};
        next_idx = 0;
        mid_begin = next_idx;
        
        for(int j = 0; j<level_num-1; j++)
        {
            mid_begin *= level_cluster[j];
            mid_end = mid_begin + level_cluster[j];
            next_idx = predict_position(cur_query, model[j],level_cluster[j],dim,mid_begin,mid_end);
            mid_begin += next_idx;
        }
        int leaf_id = mid_begin;

        int k = predict_position_arr(cur_query, centers[leaf_id],sub_cluster_num[leaf_id],dim);
        
        int * cluster_bid =  block_ids[leaf_id][k];
        int cluster_block_num = cluster_block_nums[leaf_id][k];
        double * cluster_split = split_points[leaf_id][k];
        int result_len = 0;
        double **result;
        result = new double*[B];
        int result_sort[B];
        double block_access = 0;
        auto left = (lower_bound(cluster_split, cluster_split+cluster_block_num, cur_query[0]) - cluster_split) - 1;
        left = left<0?0:left;
        int pts_num = block_count[cluster_bid[left]];
        double ** blk = blocks[cluster_bid[left]];
        block_access++;

        vector<double> result_dist;
        for (int n = 0; n<pts_num; n++)
        {
            result[result_len] = blk[n];
            result_dist.push_back(cal_dist(blk[n],cur_query));
            result_len++;
        }

        double dist_max = 0.0;
        if(result_len>=knn_k)
        {
            nth_element(result_dist.begin(),result_dist.begin()+knn_k-1,result_dist.end());
            dist_max = result_dist[knn_k-1];

        }
        else
        {
            nth_element(result_dist.begin(),result_dist.begin()+result_len-1,result_dist.end());
            dist_max = result_dist[result_len-1];
        }
        
        double search_range;
        if (dist_max != 0)
        {
            search_range = dist_max;
        }else
        {
            search_range = pow((double)knn_k/data_size, (double)1/dim) * 1;
        }
        tuple<double **, int, int> tmp_result;
       
        while(true)
        {
            tmp_result = range_query(cur_query, grid_info, grid_cluster, grid_block_nums, centers, R, cluster_block_nums,split_points, block_ids, blocks, block_count, search_range,B,scale);
            result_len = get<1>(tmp_result);
            block_access += get<2>(tmp_result);

            if (result_len>=knn_k)
            {
                result = get<0>(tmp_result);
                break;
            }else
            {
                search_range = sqrt(knn_k / (result_len + 0.5)) * search_range;
            }
        }

        double * final_result_dist = get_dist_list_arr(cur_query, result,result_len,dim);
        int final_result_sort[result_len];
        argsort(final_result_dist,result_len,final_result_sort);
        int knn_result_len = result_len<knn_k?result_len:knn_k;
        double **knn_result;
        knn_result = new double*[knn_result_len];
        for(int j=0; j<knn_result_len; j++)
        {
            knn_result[j] = result[final_result_sort[j]];
        }

        auto end = chrono::high_resolution_clock::now();
        query_time += chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        total_block_access += block_access;
        double recall = cal_recall(base[i], knn_result, knn_k, dim);
        recall_record += recall;
        
        }
    std::cout<<(query_time/query_num)<<endl;
    std::cout<<(recall_record/query_num)<<endl;
    std::cout<<(total_block_access/query_num)<<endl;

    }
    return 0;
}
