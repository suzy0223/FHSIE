#include <iostream>
#include <fstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <iterator>
#include <string>
#include <algorithm>
#include <cassert>
#include <sstream>

#include "utils.h"
using namespace std;


int check_dir(string dirpath)
{
    ifstream fin(dirpath);
    if (!fin)
    {
        string commmand = "mkdir -p " + dirpath;
        return system(commmand.c_str());
    }
    return 0;
}

int getFileLine(string file)
{
    /*
    Getting total lines of files
    Return: number of file lines
    */
    FILE * fp=fopen(file.c_str(),"rt+");

    int i=0;
    fseek(fp,0,0);
    char line[256]={0};
    while(fgets(line,255,fp))
    {
        i++;
    }
    fclose(fp);
	return i;  
}

double ** get_2d_points_space(string file, int line_num)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());
    double** data_list;
    data_list = new double*[line_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(" "));
        data_list[cur_line] = new double[vec.size()];
        if (vec.size() > 0)
        {   
            for (int i = 0; i < vec.size(); i++)
            {
                data_list[cur_line][i] = stod(vec[i]);
            }
        }
        cur_line++;
    }
    infile.close();
    return data_list;
}

double ** get_2d_points_comma(string file, int line_num)
{
    ifstream infile; 
    infile.open(file.data()); 
    assert(infile.is_open());

    double** data_list;
    data_list = new double*[line_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        data_list[cur_line] = new double[vec.size()];
        if (vec.size() > 1)
        {   
            for (int i = 0; i < vec.size(); i++)
            {
                data_list[cur_line][i] = stod(vec[i]);
            }
        }
        cur_line++;
    }
    infile.close();
    return data_list;
}


double ** get_2d_points_csv(string file, int line_num)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    double** data_list;
    data_list = new double*[line_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        data_list[cur_line] = new double[vec.size()];
        if (vec.size() > 0)
        {   
            for (int i = 0; i < vec.size(); i++)
            {
                data_list[cur_line][i] = stod(vec[i]);
            }
        }
        cur_line++;
    }
    infile.close(); 
    return data_list;
}


tuple<int *,int> get_block_count(string file)
{
    /* 
    Return data_list,block_num
    data_list: number of points in each block; 
    block_num: total block_num for whole dataset
    */
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());
    int block_num = getFileLine(file);

    int * data_list = new int[block_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        if (vec.size() > 0)
        {
            data_list[cur_line] = stoi(vec[0]);
        }
        cur_line++;
    }
    infile.close();            
    return {data_list,block_num};
}

double *** get_model(double **centers,int model_number, int cluster_number)
{
    double *** model;
    model = new double**[model_number]; 
    for (int i=0; i < model_number; i++)
    {
        model[i] = new double*[cluster_number];
        vector<vector<double>> sub_model;
        for (int j=0; j<cluster_number; j++)
        {
            model[i][j] = centers[i*cluster_number+j];
        }
    }
    return model;
}

double *** get_blocks(double **data,int *block_count,int &block_num, int &B)
{
    double*** blocks;
    blocks = new double**[block_num];
    int loc = 0;
    for (int i=0; i < block_num; i++)
    {
        blocks[i] = new double*[B];
        for (int j=0; j<block_count[i]; j++)
        {
            blocks[i][j] = data[loc];
            loc++;
        }
    }
    return blocks;
}

double *** get_blocks_reloc(double **data,int *block_count,int &block_num, int B, int max_blocks_num)
{
    double*** blocks;
    blocks = new double**[max_blocks_num];
    int loc = 0;
    for (int i=0; i < block_num; i++)
    {
        blocks[i] = new double*[B];
        for (int j=0; j<block_count[i]; j++)
        {
            blocks[i][j] = new double[2];
            blocks[i][j] = data[loc];
            loc++;
        }
    }

    for (int i=block_num; i < max_blocks_num; i++)
    {
        blocks[i] = new double*[B];
        for (int j=0; j<B; j++)
        {
            blocks[i][j] = new double[2];
        }
    }
    return blocks;
}


tuple<int **, int *> get_block_id(string file, int &cluster_num)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    int ** id;
    id = new int*[cluster_num];
    int * cluster_block_num = new int[cluster_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        id[cur_line] = new int[vec.size()];
        cluster_block_num[cur_line] = vec.size();
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
    return {id,cluster_block_num};
}

tuple<int ***,int**> get_block_ids(int model_number, string path, string filename)
{
    /*
    model_number: number of penultimate level's partition
    path filename: last level's cluster info saved path and filename

    relationship of: clusters - block ids

    Returns: block ids of each cluster; number of block in each cluster
    */
    int *** ids;
    ids = new int**[model_number];
    int **  cluster_block_nums = new int*[model_number];
    for (int i=0; i<model_number; i++)
    {
        string file = path + to_string(i) + "/original/" + filename;
        int cluster_num = getFileLine(file);
        tuple<int**,int*> block_info = get_block_id(file, cluster_num);
        ids[i] = get<0>(block_info);
        cluster_block_nums[i] = get<1>(block_info);
    }

    return {ids, cluster_block_nums};
}

double ** get_split_point(string file, int &cluster_num)
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
        split_point[cur_line] = new double[vec.size()];
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

double *** get_split_points(int model_number, string path, string filename)
{
    /*
    model_number: number of penultimate level's partition
    path filename: last level's split pts info saved path and filename

    Returns: split values (each block's begin value)
    */

    double *** split_points;
    split_points = new double**[model_number];
    for (int i=0; i<model_number; i++)
    {
        string file = path + to_string(i) + "/original/" + filename;
        int cluster_num = getFileLine(file);
        split_points[i] = get_split_point(file, cluster_num);
    }

    return split_points;
}

tuple<double **,double *, int> get_cluster_info(string file,int &total_cluster_num)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    double** center;
    double * radiu = new double[total_cluster_num];
    center = new double*[total_cluster_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        radiu[cur_line] = stod(vec[2]);
        center[cur_line] = new double[vec.size()-1];
        if (vec.size()-1 > 0)
        {
            for (int i = 0; i < vec.size()-1; i++)
            {
                center[cur_line][i] = stod(vec[i]);
            }
        }
        cur_line++;
    }
    infile.close(); 
    return {center,radiu,cur_line};
}

tuple<double ***, double**, int*> get_cluster_infos(int model_number, string path, string filename)
{
    /*
    model_number: number of penultimate level's partition
    path filename: last level's cluster info saved path and filename

    Returns: centers and radius of each cluster; number of cluster in this partition
    */
    double *** centers;
    double ** R;
    int *sub_cluster_num= new int[model_number];
    centers = new double**[model_number];
    R = new double*[model_number];
    for (int i=0; i<model_number; i++)
    {
        // cout<<i<<endl;
        string file = path + to_string(i) + "/original/" + filename;
        int total_cluster_num = getFileLine(file);
        tuple<double **, double *, int> info;
        info = get_cluster_info(file,total_cluster_num);
        centers[i] = get<0>(info);
        R[i] = get<1>(info);
        sub_cluster_num[i] = get<2>(info);
    }

    return {centers,R,sub_cluster_num};
}


void get_dist_list(double *dist_list, double *query, double **points, int &model_size,int dim, int m_begin, int m_end)
{
    
    int cur_pointer = 0;
    for (int i = m_begin; i<m_end; i++)
    {
        dist_list[cur_pointer] = sqrt(cal_dist2(query, points[i],dim));
        cur_pointer++;
    }
    return;
}

void write_points(Point * pts, int point_num, int dim, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);

    for(int i =0; i<point_num; i++)
    {
        for(int d=0;d<dim;d++)
        {
            if(d!=dim-1)
            {
                outfile<<setprecision(16)<<pts[i].pt[d]<<",";
            }
            else
            {
                outfile<<setprecision(16)<<pts[i].pt[d]<<endl;
            }
        }
    }

    outfile.close();
}


void write_double_points(double ** pts, int point_num, int dim, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);

    for(int i =0; i<point_num; i++)
    {
        for(int d=0;d<dim;d++)
        {
            if(d!=dim-1)
            {
                outfile<<setprecision(16)<<pts[i][d]<<",";
            }
            else
            {
                outfile<<setprecision(16)<<pts[i][d]<<endl;
            }
        }
    }

    outfile.close();
}

void write_block_count(int * block_count, int total_block_num, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);
    outfile<<fixed;
    for(int i=0;i<total_block_num;i++)
    {
        outfile<<to_string(block_count[i])<<endl;  
    }
    outfile.close();
}


void write_grid_ids(int ** pts, int point_num, int dim, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);

    for(int i =0; i<point_num; i++)
    {
        for(int d=0;d<dim;d++)
        {
            if(d!=dim-1)
            {
                outfile<<pts[i][d]<<",";
            }
            else
            {
                outfile<<pts[i][d]<<endl;
            }
        }
    }

    outfile.close();
}

void write_cluster_grid_info_ids(int *** block_bids,int model_number, int *sub_cluster_num, int grid_param_num, string savedir, string filepath)
{
    for(int i =0; i<model_number; i++)
    {
        string savepath = savedir+to_string(i)+ "/insert/" +filepath;
        check_dir(savedir+to_string(i)+ "/insert/");
        ofstream outfile;
        outfile.open(savepath,ios::out);
        outfile<<fixed;
        for(int j=0;j<sub_cluster_num[i];j++)
        {
            for(int n=0; n<grid_param_num; n++)
            {
                if(n!=grid_param_num-1)
                {
                    outfile<<to_string(block_bids[i][j][n])<<",";
                }
                else
                {
                    outfile<<to_string(block_bids[i][j][n])<<endl;
                }
            }
        }
        outfile.close();
    } 
}

void write_blocks(double *** blocks, int * block_count, int total_block_num, int dim, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);
    outfile<<fixed;
    for(int i=0;i<total_block_num;i++)
    {
        for(int j = 0; j<block_count[i];j++)
        {
            for(int d = 0; d<dim; d++)
            {
                outfile<<setprecision(16)<<blocks[i][j][d];
                if(d!=dim-1)
                {
                    outfile<<",";
                }
                else
                {
                    outfile<<endl;
                }
            }
        }
        
    }
   
    outfile.close();
}

void write_block_ids(int *** block_bids,int model_number, int *sub_cluster_num, int ** cluster_block_nums, string savedir, string filepath)
{
    for(int i =0; i<model_number; i++)
    {
        string savepath = savedir+to_string(i)+ "/insert/" +filepath;
        check_dir(savedir+to_string(i)+ "/insert/");
        ofstream outfile;
        outfile.open(savepath,ios::out);
        outfile<<fixed;
        for(int j=0;j<sub_cluster_num[i];j++)
        {
            for(int n=0; n<cluster_block_nums[i][j];n++)
            {
                
                outfile<<to_string(block_bids[i][j][n]);
                if(n!=cluster_block_nums[i][j]-1)
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
}

void write_split_points(double *** split_points,int model_number, int *sub_cluster_num, int ** cluster_block_nums, string savedir, string filepath)
{
    for(int i =0; i<model_number; i++)
    {
        string savepath = savedir+to_string(i)+ "/insert/" +filepath;
        check_dir(savedir+to_string(i)+ "/insert/");
        ofstream outfile;
        outfile.open(savepath,ios::out);
        outfile<<fixed;
        for(int j=0;j<sub_cluster_num[i];j++)
        {
            for(int n=0; n<cluster_block_nums[i][j];n++)
            {
                outfile<<setprecision(16)<<split_points[i][j][n];
                if(n!=cluster_block_nums[i][j]-1)
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
}


void write_split_pts(double ** pts, int point_num, int * dim, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);

    for(int i =0; i<point_num; i++)
    {
        for(int d=0;d<dim[i];d++)
        {
            if(d!=dim[i]-1)
            {
                outfile<<setprecision(16)<<pts[i][d]<<",";
            }
            else
            {
                outfile<<setprecision(16)<<pts[i][d]<<endl;
            }
        }
    }

    outfile.close();
}

void write_int_points(int ** pts, int point_num, int dim, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);
    for(int i =0; i<point_num; i++)
    {
        for(int d=0;d<dim;d++)
        {
            if(d!=dim-1)
            {
                outfile<<pts[d]<<",";
            }
            else
            {
                outfile<<pts[d]<<endl;
            }
        }
    }

    outfile.close();
}

void write_cluster_blk(int ** pts, int point_num, int * dim, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);

    for(int i =0; i<point_num; i++)
    {
        for(int d=0;d<dim[i];d++)
        {
            if(d!=dim[i]-1)
            {
                outfile<<pts[i][d]<<",";
            }
            else
            {
                outfile<<pts[i][d]<<endl;
            }
        }
    }

    outfile.close();
}


void write_cluster_info(int model_number, double *** center, double ** R, int * sub_cluster_num, int dim, string savedir, string filepath)
{
    for(int i =0; i<model_number; i++)
    {
        string savepath = savedir+to_string(i)+ "/insert/" +filepath;
        check_dir(savedir+to_string(i)+ "/insert/");
        ofstream outfile;
        outfile.open(savepath,ios::out);
        outfile<<fixed;
        for(int j=0;j<sub_cluster_num[i];j++)
        {
            for(int n=0; n<dim; n++)
            {
                outfile<<setprecision(16)<<center[i][j][n]<<",";
            }
            outfile<<setprecision(16)<<R[i][j]<<endl;
        }
        outfile.close();
    } 
}


void write_grid_cluster(vector<vector<int>> pts, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);

    for(int i =0; i<pts.size(); i++)
    {
        for(int d=0;d<pts[i].size();d++)
        {
            if(d!=pts[i].size()-1)
            {
                outfile<<pts[i][d]<<",";
            }
            else
            {
                outfile<<pts[i][d]<<endl;
            }
        }
        if(pts[i].size()==0)
        {
            outfile<<to_string(-1)<<endl;
        }
    }

    outfile.close();
}


void write_int_list(vector<int> pts, int point_num, string filepath)
{
    ofstream outfile;
    outfile.open(filepath,ios::out);

    for(int i =0; i<point_num; i++)
    {
        outfile<<pts[i]<<endl;
    }

    outfile.close();
}


double cal_dist2(double * point1, double * point2, int dim)
{
    double temp_dist = 0.0;
    for(int i = 0; i<dim; i++)
    {
        temp_dist += pow((point1[i] - point2[i]), 2);
    }
    return temp_dist;
}


int predict_position(double *query,double **model,int &model_size, int dim, int m_begin, int m_end)
{
    /*
    predict next cluster at next level
    */
    int next_id;
    double *dist_list = new double[model_size];
    get_dist_list(dist_list,query,model,model_size,dim,m_begin,m_end);
    next_id = argmin(dist_list, dist_list+model_size);
    delete dist_list;
    return next_id;
}

int predict_position_arr(double *query,double **model,int &model_size, int dim)
{
    int next_id;
    double *dist_list = get_dist_list_arr(query,model,model_size,dim);
    next_id = argmin(dist_list, dist_list+model_size);
    delete dist_list;
    return next_id;
}

int scan_block(double *query, double **block, int &block_cur_count, int &dim)
{
    for (int m = 0; m < block_cur_count; m++)
    {
        bool flag = true;
        for (int d = 0; d<dim; d++)
        {
            if (query[d]!=block[m][d]){
                flag = false;
                break;
            }
        }

        if (flag){
            return 1;
        }
    }
    return 0;
}

double * get_grid_info(string file, int &grid_param_num)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    double *data_list = new double[grid_param_num+2];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        if (vec.size() > 0)
        {
            data_list[cur_line] = stod(vec[0]);
        }
        cur_line++;
    }
    data_list[cur_line] = data_list[6]-data_list[4];
    data_list[cur_line+1] = data_list[7]-data_list[5];
    infile.close();
    return data_list;
}

tuple<int **, int *> get_cluster_id(string file, int &cluster_num)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    int ** id;
    id = new int*[cluster_num];
    int * cluster_block_num = new int[cluster_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        id[cur_line] = new int[vec.size()];
        cluster_block_num[cur_line] = vec.size();
        if (vec.size() > 0)
        {
            for (int i = 0; i < vec.size(); i++)
            {
                id[cur_line][i] = stoi(vec[i]);
                if(stoi(vec[i])==-1){cluster_block_num[cur_line]=0;}
            }
        }
        cur_line++;
    }
    infile.close();
    return {id,cluster_block_num};
}


vector<vector<double>> get_points(string file)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());

    vector<vector<double>> data_list;

    string line;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(" "));
        vector<double> point;
        if (vec.size() > 1)
        {
            for (int i = 0; i < vec.size(); i++)
            {
                point.push_back(stod(vec[i]));
            }
            data_list.push_back(point);
        }
    }
    infile.close();
    return data_list;
}

double *** knn_base_result(string file)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open());
    int fileline = getFileLine(file);

    double ***base;
    base = new double**[fileline];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(" "));
        base[cur_line] = new double*[vec.size()];
        if (vec.size() > 0)
        {
            for(int i=0; i<vec.size(); i++)
            {
                vector<string> elemt;
                boost::algorithm::split(elemt, vec[i], boost::is_any_of(","));
                base[cur_line][i] = new double[elemt.size()];
                for(int j =0; j<elemt.size(); j++)
                {
                    base[cur_line][i][j] = stod(elemt[j]);
                }
            }
        }
        cur_line++;
    }
    infile.close();
    return base;
}

int ** get_grid_id(string file,int &cluster_num)
{
    ifstream infile; 
    infile.open(file.data());
    assert(infile.is_open()); 

    int ** id;
    id = new int*[cluster_num];

    string line;
    int cur_line = 0;
    while(getline(infile,line))
    {
        vector<string> vec;
        boost::algorithm::split(vec, line, boost::is_any_of(","));
        id[cur_line] = new int[vec.size()];
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
    return id;
}

int *** get_grid_ids(int model_number, string path, string filename)
{
    int *** ids;
    ids = new int**[model_number];
    for (int i=0; i<model_number; i++)
    {
        string file = path + to_string(i) + "/original/" + filename;
        int cluster_num = getFileLine(file);
        ids[i] = get_grid_id(file,cluster_num);
    }

    return ids;
}

double * get_dist_list_arr(double *query, double **points, int &model_size,int dim)
{
    // int points_num = points.size();
    // static double dist_list[36];
    double *dist_list = new double[model_size];
    for (int i = 0; i<model_size; i++)
    {
        dist_list[i] = sqrt(cal_dist2(query, points[i],dim));
    }
    return dist_list;
}

vector<int> remove_duplicates(vector<int> vec) 
{

	sort(vec.begin(), vec.end());
    vec.erase(unique(vec.begin(),vec.end()),vec.end());
    return vec;
}

int ** convert2multiclusterid(vector<int> cluster_long_id, int &id_nums, int scale)
{
    int ** cluster_2d_ids;
    cluster_2d_ids = new int*[id_nums];
    for (int i = 0; i<id_nums; i++)
    { 
        int id = cluster_long_id[i];
        cluster_2d_ids[i] = new int[2];
        cluster_2d_ids[i][0] = id%scale;
        cluster_2d_ids[i][1] = (id / scale) % scale;
    }

    return cluster_2d_ids;
}
