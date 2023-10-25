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
#include "Constants.h"
#include "utils.h"
  
using namespace std;


int main(){
    
    cout<<Constants::DATASETS<<endl;
    string dataset = Constants::DATASETS;
    string model_record_path = Constants::MODEL_R_PATH;

    int level_num = Constants::LEVEL;
    int level_cluster[level_num] = {27,27,27};

    int dim = Constants::DIM;
    int B = Constants::B;
    int partition = 1;

    int dat_line = getFileLine(Constants::INDEX_PATH);
    // loading all points
    double ** dat =  get_2d_points_csv(Constants::INDEX_PATH,dat_line);

    tuple<int *,int> block_info = get_block_count(Constants::BLOCK_INFO);
    
    int * block_count = get<0>(block_info);
    int total_block_num = get<1>(block_info);

    // Load cluster centres of K-Means model at each level
    double *** model;
    model = new double**[level_num-1];
    for(int i = 0; i<level_num-1; i++)
    {
        partition *= level_cluster[i];
        model[i] = get_2d_points_csv(Constants::K_MEANS_PATH + "/level"+to_string(i)+"model.csv",partition);
    }
    
    // allocate points to each blocks to build the index
    double *** blocks = get_blocks(dat,block_count,total_block_num, B);

    tuple<int ***, int **> cluter_block = get_block_ids(partition, model_record_path, Constants::CLUSTER_BLOCK);
    int *** block_ids = get<0>(cluter_block);
    int ** cluster_block_nums = get<1>(cluter_block);

    double *** split_points = get_split_points(partition, model_record_path, Constants::SPLIT_PTS);
    tuple<double ***, double**, int*> cluster_infos = get_cluster_infos(partition, model_record_path, Constants::CLUSTER_INFO);
    double*** centers;
    int *sub_cluster_num;
    centers = get<0>(cluster_infos);
    sub_cluster_num = get<2>(cluster_infos);
    
    int prediction=0; // record hit times
    long query_time=0;
    double total_block_access = 0;
    double ** query = dat; // all points as query
    int query_num = dat_line;
    int next_idx,mid_begin,mid_end;

    cout<<"query_num: "<<query_num<<endl;
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
        int k = predict_position(cur_query, centers[leaf_id],sub_cluster_num[leaf_id],dim,0,sub_cluster_num[leaf_id]);
        
        // get cluster info
        int * cluster_bid =  block_ids[leaf_id][k];
        double * cluster_split = split_points[leaf_id][k];
        int cluster_block_num = cluster_block_nums[leaf_id][k];

        // use split values to filter blocks
        auto left = (lower_bound(cluster_split, cluster_split+cluster_block_num, cur_query[0]) - cluster_split) - 1;
        auto right = upper_bound(cluster_split, cluster_split+cluster_block_num, cur_query[0]) - cluster_split;
        left = left<0?0:left;
        right = right+1;
        right = right>cluster_block_num?cluster_block_num:right;

        int c = 0;
        double block_access = 0;
        for (int j = left; j<right; j++)
        {
            int bid = cluster_bid[j];
            c = scan_block(cur_query, blocks[bid],block_count[bid],dim);
            block_access++;
            if (c==1){break;}
        }
        auto end = chrono::high_resolution_clock::now();
        query_time += chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        prediction += c;
        total_block_access += block_access;
    }
    std::cout<<"hit times: "<<prediction<<endl;
    std::cout<<(query_time/query_num)<<endl;
    std::cout<<(total_block_access/query_num)<<endl;

    return 0;
}
