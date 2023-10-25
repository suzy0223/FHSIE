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
#include <iomanip>
#include "utils.h"
#include "Constants.h"
  
using namespace std;

double cal_dist_arr(double *point1, double *point2)
{
    double temp_dist = sqrt(pow((*point1 - *point2), 2) + pow((*(point1+1) - *(point2+1)), 2));
    return temp_dist;
}


tuple<int **, int *, int *> get_block_id_insert(string file, int &cluster_num, double insert_ratio)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    int ** id;
    id = new int*[cluster_num];
    int * cluster_block_num = new int[cluster_num];
    int * max_cluster_block_num = new int[cluster_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        int allocate_size = (vec.size())*(1+insert_ratio*80);
        id[cur_line] = new int[allocate_size];
        cluster_block_num[cur_line] = vec.size();
        max_cluster_block_num[cur_line] = allocate_size;
        if (vec.size() > 0)
        {
            for (int i = 0; i < vec.size(); i++)
            {
                id[cur_line][i] = stoi(vec[i]);
            }
        }
        cur_line++;
    }
    infile.close();
    return {id,cluster_block_num,max_cluster_block_num};
}

tuple<int ***,int**,int**> get_block_ids_reloc(int model_number, double insert_ratio, string path, string filename)
{
    int *** ids;
    ids = new int**[model_number];
    int **  cluster_block_nums = new int*[model_number];
    int **  max_cluster_block_nums = new int*[model_number];
    for (int i=0; i<model_number; i++)
    {
        string file = path + to_string(i) + "/original/" + filename;
        int cluster_num = getFileLine(file);
        tuple<int**,int*,int*> block_info = get_block_id_insert(file, cluster_num, insert_ratio);
        ids[i] = get<0>(block_info);
        cluster_block_nums[i] = get<1>(block_info);
        max_cluster_block_nums[i] = get<2>(block_info);
    }
    
    return {ids, cluster_block_nums,max_cluster_block_nums};
}

void update_grid_cluster(int ** grid_cluster,int grid_num, int *grid_block_nums, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);
    outfile<<fixed;
    for(int i =0; i<grid_num; i++)
    {          
        for(int j=0;j<grid_block_nums[i];j++)
        {
            outfile<<to_string(grid_cluster[i][j]);
            if(j!=grid_block_nums[i]-1)
            {
                outfile<<",";
            }
            else
            {
                outfile<<endl;
            }
        }
    }
    outfile.close();
}

double ** get_split_point_insert(string file, int &cluster_num, double insert_ratio)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    double ** split_point;
    split_point = new double*[cluster_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        int allocate_size = (vec.size())*(1+insert_ratio*80);
        split_point[cur_line] = new double[allocate_size];
        if (vec.size() > 0)
        {
            for (int i = 0; i < vec.size(); i++)
            {
                split_point[cur_line][i] = stod(vec[i]);
            }
        }
        cur_line++;
    }
    infile.close();
    return split_point;
}

double *** get_split_points_insert(int model_number, double insert_ratio, string path, string filename)
{
    double *** split_points;
    split_points = new double**[model_number];
    for (int i=0; i<model_number; i++)
    {
        string file = path + to_string(i) + "/original/" + filename;
        int cluster_num = getFileLine(file);
        split_points[i] = get_split_point_insert(file, cluster_num,insert_ratio);
    }

    return split_points;
}


tuple<double ***, int *, double ***, int ***, int **, int ***, int **, int *, int, double **> insert_points(double *query, double *** model, int * level_cluster, int level_num, int partition, int B, int dim, int * pid, int total_parts_num, int *sub_cluster_num, double * grid_info,int ** grid_cluster, int * grid_block_nums,double*** centers,double** R,int ** cluster_block_nums,double *** split_points, int *** block_ids, double *** blocks, int * block_count, int &block_num, int *** cluster_intersect_grid_id,int max_blocks_num, double insert_ratio, int ** max_cluster_block_nums, int scale)
{
    int next_idx = 0;
    int mid_begin = next_idx;
    int mid_end;

    for(int j = 0; j<level_num-1; j++)
    {
        mid_begin *= level_cluster[j];
        mid_end = mid_begin + level_cluster[j];
        next_idx = predict_position(query, model[j],level_cluster[j],dim,mid_begin,mid_end);
        mid_begin += next_idx;
    }

    int leaf_id = mid_begin; 
    // if this part does not have a block
    if (!(sub_cluster_num[leaf_id]^0))
    {
        auto p_tmp_id = lower_bound(pid, pid+total_parts_num, leaf_id)- pid;
        leaf_id = pid[p_tmp_id];
    }

    int k = predict_position_arr(query, centers[leaf_id],sub_cluster_num[leaf_id],dim);
    // cout<<"leaf_id = "<<leaf_id<<" k = "<<k<<endl;

    int * cluster_bid =  block_ids[leaf_id][k];
    double * cluster_center = centers[leaf_id][k];
    double * cluster_split = split_points[leaf_id][k];
    double cluster_r = R[leaf_id][k];
    int cluster_block_num = cluster_block_nums[leaf_id][k];
    int max_cluster_block_num = max_cluster_block_nums[leaf_id][k];
    
    int block_loc = (upper_bound(cluster_split,cluster_split + cluster_block_num, query[0]) - cluster_split)-1;
    // cout<<"block_loc = "<<block_loc<<endl;
    block_loc = block_loc<0?0:block_loc;
    int block_bid = cluster_bid[block_loc];
    // cout<<"block_bid = "<<block_bid<<endl;
    // get block size before insert
    int block_size = block_count[block_bid];
    
    // insert the points(sorted) and update block size
    double tmp_dat_x[block_size];
    for(int i = 0; i<block_size; i++)
    {
        tmp_dat_x[i] = blocks[block_bid][i][0];
    }
    int point_loc = upper_bound(tmp_dat_x,tmp_dat_x+block_size,query[0])-tmp_dat_x;
    int next_pos = block_size-1;
    int cur_pos = block_size;
    while(cur_pos>=point_loc)
    {
        blocks[block_bid][cur_pos] = blocks[block_bid][next_pos];
        next_pos--;
        cur_pos--;
    }

    blocks[block_bid][point_loc] = query;
    block_size += 1;

    if(block_size<=B)
    {
        block_count[block_bid] = block_size;
        if(query[0]<cluster_split[block_loc])
        {
            split_points[leaf_id][k][block_loc] = query[0];
        }
        // check the R, if change update the R and the grid
        double points_dist = cal_dist_arr(query,cluster_center);
        if (points_dist > cluster_r)
        {
            R[leaf_id][k] = points_dist;
            //  compute grid query belongs to and update if necessary
            int xid = ((query[0] - grid_info[4])<0?0:(query[0] - grid_info[4])) / grid_info[2];
            int yid = ((query[1] - grid_info[5])<0?0:(query[1] - grid_info[5])) / grid_info[3];

            int cluster_x_grid_b = cluster_intersect_grid_id[leaf_id][k][0];
            int cluster_x_grid_e = cluster_intersect_grid_id[leaf_id][k][1];
            int cluster_y_grid_b = cluster_intersect_grid_id[leaf_id][k][2];
            int cluster_y_grid_e = cluster_intersect_grid_id[leaf_id][k][3];

            int tmp_x_begin = xid<cluster_x_grid_b?xid:cluster_x_grid_b;
            int tmp_x_end = xid<cluster_x_grid_e?cluster_x_grid_e:xid;
            int tmp_y_begin = yid<cluster_y_grid_b?yid:cluster_y_grid_b;
            int tmp_y_end = yid<cluster_y_grid_e?cluster_y_grid_e:yid;
            if ((tmp_x_begin != cluster_x_grid_b) || (tmp_x_end != cluster_x_grid_e) || (tmp_y_begin != cluster_y_grid_b) || (tmp_y_end != cluster_y_grid_e))
            {
                cluster_intersect_grid_id[leaf_id][k][0] = tmp_x_begin;
                cluster_intersect_grid_id[leaf_id][k][1] = tmp_x_end;
                cluster_intersect_grid_id[leaf_id][k][2] = tmp_y_begin;
                cluster_intersect_grid_id[leaf_id][k][3] = tmp_y_end;

                for(int m =tmp_x_begin;m<=tmp_x_end;m++)
                {
                    for(int n = tmp_y_begin;n<=tmp_y_end;n++)
                    {
                        int grid_id = int(m * grid_info[1] + n);
                        int cid = leaf_id+k*scale;
                        
                        auto tmp_loc = upper_bound(grid_cluster[grid_id],grid_cluster[grid_id]+grid_block_nums[grid_id],cid) - grid_cluster[grid_id];
                        if ((cid!=grid_cluster[grid_id][tmp_loc])&&(cid!=(grid_cluster[grid_id][(tmp_loc-1)<0?0:(tmp_loc-1)]))&&(cid!=(grid_cluster[grid_id][(tmp_loc+1)>(grid_block_nums[grid_id])?(grid_block_nums[grid_id]):(tmp_loc+1)])))
                        {
                            grid_block_nums[grid_id]+=1;
                            // resize the grid_cluster
                            int * new_grid_cluster = new int[grid_block_nums[grid_id]];
                            new_grid_cluster = grid_cluster[grid_id];
                            new_grid_cluster[grid_block_nums[grid_id]-1] = cid;
                            grid_cluster[grid_id] = new_grid_cluster;
                        }
                    }
                }       
            }       
        }
        return {blocks,block_count,split_points,block_ids,cluster_block_nums,cluster_intersect_grid_id,grid_cluster,grid_block_nums,block_num,R};
    }
    else
    {
        int sub_block_size = block_size*0.5;
        if(block_num+1<=max_blocks_num)
        {
            for(int i = 0;i<(block_size - sub_block_size);i++)
            {
                blocks[block_num][i] = blocks[block_bid][sub_block_size+i];
            }
        }
        else
        {
            // build a new blocks double *** blocks
            double *** new_blocks;
            new_blocks = new double**[block_num+1];
            new_blocks = blocks;

            new_blocks[block_num] = new double*[B+1];
            for(int i = 0;i<(block_size - sub_block_size);i++)
            {
                new_blocks[block_num][i] = new double[2];
                memcpy(new_blocks[block_num][i], blocks[block_bid][sub_block_size+i],2*sizeof(double));
            }
            blocks = new_blocks;
        }

        // upload block count
        block_count[block_bid] = sub_block_size;
        int * new_block_count = new int[block_num+1];
        memcpy(new_block_count,block_count,block_num*sizeof(int));
        delete [] block_count;
        new_block_count[block_num] = block_size - sub_block_size;

        // update block id, split points and cluster block numbers
        int leaf_cluster_num;

        for(int i = 0; i<partition; i++)
        {
            leaf_cluster_num = sub_cluster_num[i];
            
            // loacte the sub dataset(leafnode)
            if(i==leaf_id)
            {
                for(int j = 0; j<leaf_cluster_num; j++)
                {
                    if(j==k)
                    {
                        int *tmp_block_ids;
                        tmp_block_ids = new int[cluster_block_num+1];
                        if(cluster_block_num>=max_cluster_block_num){cout<<"block num not enough, leaf id = "<<leaf_id<<" k = "<<k<<endl;}
                        memcpy(tmp_block_ids,block_ids[i][j],cluster_block_num*sizeof(int));
                        block_ids[i][j][block_loc+1] = block_num;

                        double *tmp_split_points;
                        tmp_split_points = new double[cluster_block_num+1];
                        memcpy(tmp_split_points,split_points[i][j],cluster_block_num*sizeof(double));

                        split_points[i][j][block_loc] = blocks[block_bid][0][0];
                        split_points[i][j][block_loc+1] = blocks[block_bid][sub_block_size][0];
                        
                        for(int n = block_loc+2; n<=cluster_block_num; n++)
                        {
                            block_ids[i][j][n] = tmp_block_ids[n-1];
                            split_points[i][j][n] = tmp_split_points[n-1];
                        }

                        cluster_block_nums[i][j] += 1;

                        delete [] tmp_block_ids;
                        delete [] tmp_split_points;
                    }
                }
            }
        }
        // check the R, if change update the R and the grid
        double points_dist = cal_dist_arr(query,cluster_center);
        if (points_dist > cluster_r)
        {
            R[leaf_id][k] = points_dist;
            //  compute grid query belongs to and update if necessary
            int xid = ((query[0] - grid_info[4])<0?0:(query[0] - grid_info[4])) / grid_info[2];
            int yid = ((query[1] - grid_info[5])<0?0:(query[1] - grid_info[5])) / grid_info[3];

            int cluster_x_grid_b = cluster_intersect_grid_id[leaf_id][k][0];
            int cluster_x_grid_e = cluster_intersect_grid_id[leaf_id][k][1];
            int cluster_y_grid_b = cluster_intersect_grid_id[leaf_id][k][2];
            int cluster_y_grid_e = cluster_intersect_grid_id[leaf_id][k][3];

            int tmp_x_begin = xid<cluster_x_grid_b?xid:cluster_x_grid_b;
            int tmp_x_end = xid<cluster_x_grid_e?cluster_x_grid_e:xid;
            int tmp_y_begin = yid<cluster_y_grid_b?yid:cluster_y_grid_b;
            int tmp_y_end = yid<cluster_y_grid_e?cluster_y_grid_e:yid;
            if ((tmp_x_begin != cluster_x_grid_b) || (tmp_x_end != cluster_x_grid_e) || (tmp_y_begin != cluster_y_grid_b) || (tmp_y_end != cluster_y_grid_e))
            {
                cluster_intersect_grid_id[leaf_id][k][0] = tmp_x_begin;
                cluster_intersect_grid_id[leaf_id][k][1] = tmp_x_end;
                cluster_intersect_grid_id[leaf_id][k][2] = tmp_y_begin;
                cluster_intersect_grid_id[leaf_id][k][3] = tmp_y_end;

                for(int m =tmp_x_begin;m<=tmp_x_end;m++)
                {
                    for(int n = tmp_y_begin;n<=tmp_y_end;n++)
                    {
                        int grid_id = int(m * grid_info[1] + n);
                        int cid = leaf_id+k*scale;

                        // if(grid_block_nums[grid_id]==0)
                        // {
                        //     grid_block_nums[grid_id]+=1;
                        //     grid_cluster[grid_id][0] = cid;
                        //     continue;
                        // }
                        
                        auto tmp_loc = upper_bound(grid_cluster[grid_id],grid_cluster[grid_id]+grid_block_nums[grid_id],cid) - grid_cluster[grid_id];
                        if ((cid!=grid_cluster[grid_id][tmp_loc])&&(cid!=(grid_cluster[grid_id][(tmp_loc-1)<0?0:(tmp_loc-1)]))&&(cid!=(grid_cluster[grid_id][(tmp_loc+1)>(grid_block_nums[grid_id])?(grid_block_nums[grid_id]):(tmp_loc+1)])))
                        {
                            grid_block_nums[grid_id]+=1;
                            // resize the grid_cluster
                            int * new_grid_cluster = new int[grid_block_nums[grid_id]];
                            new_grid_cluster = grid_cluster[grid_id];
                            new_grid_cluster[grid_block_nums[grid_id]-1] = cid;
                            grid_cluster[grid_id] = new_grid_cluster;
                        }
                    }
                }       
            }       
        }
        
        return {blocks,new_block_count,split_points,block_ids,cluster_block_nums,cluster_intersect_grid_id,grid_cluster,grid_block_nums,block_num+1,R};

    }
}


int main(){
    string dataset = Constants::DATASETS;
    string model_record_path = Constants::MODEL_R_PATH;
    int dim = Constants::DIM;
    int B = Constants::B;
    int scale = Constants::SCALE;
    string distribution = Constants::DISTRIBUTION;
    string workload_path = Constants::INSERT_WORKLOAD_PATH;

    int level_num = Constants::LEVEL;
    int level_cluster[level_num] = {20,20,20,20};
    int partition = 1;
    string insert_distribution = Constants::INSERT_DISTRIBUTION;
    double insert_ratio = 1.0;
    double insert_points_num = 16000000;
    string insert_num = "16000000";
    cout<<"getting datasets and workload..."<<endl;
    int data_size = getFileLine(Constants::INDEX_PATH);
    double ** dat =  get_2d_points_csv(Constants::INDEX_PATH,data_size);
    
    cout<<Constants::INSERT_WORKLOAD_PATH +"insert_" + distribution+ "_"+insert_distribution+"_"+to_string(data_size)+ "_" + insert_num +"_.csv"<<endl;
    int query_num = getFileLine(Constants::INSERT_WORKLOAD_PATH +"insert_" + distribution+ "_"+insert_distribution+"_"+to_string(data_size)+ "_" + insert_num +"_.csv");
    double ** query = get_2d_points_csv(Constants::INSERT_WORKLOAD_PATH +"insert_" + distribution+ "_"+insert_distribution+"_"+to_string(data_size)+ "_" + insert_num +"_.csv",query_num);
    cout<<query[0][0]<<", "<<query[0][1]<<endl;
    cout<<"getting blocks..."<<endl;
    tuple<int *,int> block_info = get_block_count(Constants::BLOCK_INFO);
    int * block_count = get<0>(block_info);
    int total_block_num = get<1>(block_info);
    cout<<"total_block_num = "<<total_block_num<<endl;
    int max_blocks_num = total_block_num*(1+80*insert_ratio);
    cout<<"max_blocks_num = "<<max_blocks_num<<endl;
    double *** blocks = get_blocks_reloc(dat,block_count,total_block_num,B+1,max_blocks_num);

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

    tuple<int ***,int**,int**> cluter_block = get_block_ids_reloc(partition, insert_ratio, model_record_path, "cluster_block.csv");
    int *** block_ids = get<0>(cluter_block);
    int ** cluster_block_nums = get<1>(cluter_block);
    int ** max_cluster_block_nums = get<2>(cluter_block);

    double *** split_points = get_split_points_insert(partition, insert_ratio, model_record_path, "split_pts.csv");
    tuple<double ***, double**, int*> cluster_infos = get_cluster_infos(partition, model_record_path, "cluster_info.csv");
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
    double * grid_info = get_grid_info(Constants::GRID_INFO, grid_param_num);

    int *** cluster_intersect_grid_id = get_grid_ids(partition, model_record_path, "grid_ids.csv");
    
    cout<<"start query now"<<endl;
    double query_time=0;
    for (int i = 0; i< query_num; i++)
    {
        // cout<<i<<" = "<<query[i][0]<<","<<query[i][1]<<endl;
        auto start = chrono::high_resolution_clock::now();
        tuple<double ***, int *, double ***, int ***, int **, int ***, int **, int *, int, double **> update_result = insert_points(query[i], model, level_cluster, level_num, partition, B, dim, pid,total_parts_num, sub_cluster_num, grid_info, grid_cluster, grid_block_nums, centers, R, cluster_block_nums, split_points, block_ids, blocks, block_count, total_block_num, cluster_intersect_grid_id, max_blocks_num, insert_ratio, max_cluster_block_nums, scale);
        blocks = get<0>(update_result);
        block_count = get<1>(update_result);
        split_points = get<2>(update_result);
        block_ids = get<3>(update_result);
        cluster_block_nums = get<4>(update_result);
        cluster_intersect_grid_id = get<5>(update_result);
        grid_cluster= get<6>(update_result);
        grid_block_nums= get<7>(update_result);
        total_block_num = get<8>(update_result);
        R = get<9>(update_result);
        auto end = chrono::high_resolution_clock::now();
        query_time += chrono::duration_cast<chrono::nanoseconds>(end - start).count();
        if((i==1600000-1)||(i==4000000-1)||(i==8000000-1)||(i==12000000-1))
        {
            std::cout<<(query_time/(double)(i+1))<<endl;
            write_blocks(blocks, block_count,total_block_num, dim,  "dataset/" + dataset + "/insert_"+to_string(i+1)+"_"+insert_distribution+"_.csv");
            cout<<"writing results to file"<<endl;
            write_block_count(block_count,total_block_num,model_record_path + "block_count_insert_"+to_string(i+1)+"_"+insert_distribution+"_.csv");
            cout<<"writing block_count to file"<<endl;
            write_split_points(split_points,partition,sub_cluster_num,cluster_block_nums,model_record_path, "split_pts_insert_"+to_string(i+1)+"_"+insert_distribution+"_.csv");
            cout<<"writing split_points to file"<<endl;
            write_block_ids(block_ids,partition,sub_cluster_num,cluster_block_nums,model_record_path, "cluster_block_insert_"+to_string(i+1)+"_"+insert_distribution+"_.csv");
            cout<<"writing block_ids to file"<<endl;
            write_cluster_grid_info_ids(cluster_intersect_grid_id,partition,sub_cluster_num,4,model_record_path, "grid_ids_insert_"+to_string(i+1)+"_"+insert_distribution+"_.csv");
            cout<<"writing cluster_intersect_grid_id to file"<<endl;
            write_cluster_info(partition, centers, R, sub_cluster_num, dim, model_record_path, "cluster_info_insert_"+to_string(i+1)+"_"+insert_distribution+"_.csv");
            cout<<"writing cluster_info to file"<<endl;
            update_grid_cluster(grid_cluster,grid_num,grid_block_nums,model_record_path + "grid_cluster_insert_"+to_string(i+1)+"_"+insert_distribution+"_.csv");
            cout<<"writing cluster_info to file"<<endl;
        }
    }
    std::cout<<(query_time/query_num)<<endl;
    
    write_blocks(blocks, block_count,total_block_num, dim,  "dataset/" + dataset + "/insert_"+to_string(insert_points_num)+"_"+insert_distribution+"_.csv");
    cout<<"writing results to file"<<endl;
    write_block_count(block_count,total_block_num, model_record_path + "block_count_insert_"+to_string(insert_points_num)+"_"+insert_distribution+"_.csv");
    cout<<"writing block_count to file"<<endl;
    write_split_points(split_points,partition,sub_cluster_num,cluster_block_nums,model_record_path, "split_pts_insert_"+to_string(insert_points_num)+"_"+insert_distribution+"_.csv");
    cout<<"writing split_points to file"<<endl;
    write_block_ids(block_ids,partition,sub_cluster_num,cluster_block_nums,model_record_path, "cluster_block_insert_"+to_string(insert_points_num)+"_"+insert_distribution+"_.csv");
    cout<<"writing block_ids to file"<<endl;
    write_cluster_grid_info_ids(cluster_intersect_grid_id,partition,sub_cluster_num,4,model_record_path, "grid_ids_insert_"+to_string(insert_points_num)+"_"+insert_distribution+"_.csv");
    cout<<"writing cluster_intersect_grid_id to file"<<endl;
    write_cluster_info(partition, centers, R, sub_cluster_num, dim, model_record_path, "cluster_info_insert_"+to_string(insert_points_num)+"_"+insert_distribution+"_.csv");
    cout<<"writing cluster_info to file"<<endl;
    update_grid_cluster(grid_cluster,grid_num,grid_block_nums,model_record_path + "grid_cluster_insert_"+to_string(insert_points_num)+"_"+insert_distribution+"_.csv");
    cout<<"writing cluster_info to file"<<endl;
    return 0;
}
