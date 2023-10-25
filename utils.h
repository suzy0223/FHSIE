#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include <vector>
#include <typeinfo> 
#include <tuple>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <map>
#include <numeric> 
#include <iomanip>
#include <random>
using namespace std;

template<class T>
inline int argmin(T first, T last) 
{
    return distance(first, min_element(first, last));
}

template<class T>
inline int argmax(T first, T last) 
{
    return distance(first, max_element(first, last));
}

template<class T>
inline void argsort(T *array, int num, int *index) {
    const auto function = [array](int a, int b) noexcept -> bool {
        return array[a] < array[b];
    };
    assert(num < INT_MAX);
    int *temp = new int[num];
    std::iota(temp, temp + num, 0);
    std::sort(temp, temp + num, function);
    memcpy(index, temp, num * sizeof(int));
    delete[] temp;
}

struct Point
{
    double * pt;
    int cluster;
};

int check_dir(string dirpath);

int getFileLine(string file);
double ** get_2d_points_comma(string file, int line_num);
double ** get_2d_points_space(string file, int line_num);
double ** get_2d_points_csv(string file, int line_num);
tuple<int *,int> get_block_count(string file);
double *** get_model(double **centers,int model_number, int cluster_number);
double *** get_blocks(double **data,int *block_count,int &block_num, int &B);
double *** get_blocks_reloc(double **data,int *block_count,int &block_num, int B, int max_blocks_num);
tuple<int **, int *> get_block_id(string file, int &cluster_num);
tuple<int ***,int**> get_block_ids(int model_number, string path, string filename);
double ** get_split_point(string file, int &cluster_num);
double *** get_split_points(int model_number, string path, string filename);
tuple<double **,double *, int> get_cluster_info(string file,int &total_cluster_num);
tuple<double ***, double**, int*> get_cluster_infos(int model_number, string path, string filename);
double cal_dist2(double * point1, double * point2, int dim);

void write_points(Point * pts, int point_num, int dim, string filepath);
void write_double_points(double ** pts, int point_num, int dim, string filepath);
void write_grid_ids(int ** pts, int point_num, int dim, string filepath);
void write_split_points(double *** split_points,int model_number, int *sub_cluster_num, int ** cluster_block_nums, string savedir, string filepath);
void write_split_pts(double ** pts, int point_num, int * dim, string filepath);

void write_blocks(double *** blocks, int * block_count, int total_block_num, int dim, string filepath);
void write_block_ids(int *** block_bids,int model_number, int *sub_cluster_num, int ** cluster_block_nums, string savedir, string filepath);

void write_int_points(int ** pts, int point_num, int dim, string filepath);
void write_cluster_blk(int ** pts, int point_num, int * dim, string filepath);
void write_grid_cluster(vector<vector<int>> pts, string filepath);
void write_int_list(vector<int> pts, int point_num, string filepath);
void write_block_count(int * block_count, int total_block_num, string filepath);
void write_cluster_grid_info_ids(int *** block_bids,int model_number, int *sub_cluster_num, int grid_param_num, string savedir, string filepath);
void write_cluster_info(int model_number, double *** center, double ** R, int * sub_cluster_num, int dim, string savedir, string filepath);

void get_dist_list(double *dist_list, double *query, double **points, int &model_size,int dim, int m_begin, int m_end);
int predict_position(double *query,double **model,int &model_size, int dim, int m_begin, int m_end);
int predict_position_arr(double *query,double **model,int &model_size, int dim);
int scan_block(double *query, double **block, int &block_cur_count, int &dim);

double * get_grid_info(string file, int &grid_param_num);
tuple<int **, int *> get_cluster_id(string file, int &cluster_num);
vector<vector<double>> get_points(string file);
double *** knn_base_result(string file);
int ** get_grid_id(string file,int &cluster_num);
int *** get_grid_ids(int model_number, string path, string filename);
double * get_dist_list_arr(double *query, double **points, int &model_size,int dim);

vector<int> remove_duplicates(vector<int> vec);
int ** convert2multiclusterid(vector<int> cluster_long_id, int &id_nums, int scale);

#endif