#include <iostream>
#include <fstream>
#include <cassert>
#include <string>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <vector>
#include <set>
#include <typeinfo> 
#include <tuple>
#include <math.h>
#include <algorithm>
#include <chrono>
#include <map>
#include <numeric> 
#include <iomanip>
#include <random>
#include "utils.h"
#include "Constants.h"

using namespace std;

void generate_points(Point * pts, double ** dat, int point_num, int dim)
{
    for(int i=0; i<point_num; i++)
    {
        memcpy(pts[i].pt,dat[i],sizeof(double)*dim);
        pts[i].pt = dat[i];
        pts[i].cluster = 0;
    }

    return;
}


int nearest(Point pt, Point * cent, int cluster_num, int dim)
{
    int cluster_idx,i;
    double dist, min_d;
    min_d = sqrt(cal_dist2(pt.pt, cent[0].pt, dim));
    cluster_idx = 0;
    for(i=1; i<cluster_num; i++)
    {
        dist = sqrt(cal_dist2(pt.pt,cent[i].pt,dim));
        if(dist<min_d)
        {
            min_d = dist;
            cluster_idx = i;
        }
    }
    return cluster_idx;
}

void kmeans_plusplus(Point * pts, int point_num, Point * centroids, int cluster_num, int dim)
{
	int j,cid;
	int selectedIndex;
	int cluster;
	double sum;
	double d;
	double random;	
	double * cumulativeDistances;
	double * shortestDistance;
 
	cumulativeDistances = new double[point_num];
    shortestDistance = new double[point_num];

	/* Pick the first cluster centroids at random. */
	selectedIndex = rand() % point_num;
    memcpy(centroids[0].pt,pts[selectedIndex].pt,sizeof(double)*dim);

    /* init the shortestdist*/
    for (j = 0; j < point_num; ++j)
    {
        shortestDistance[j] = HUGE_VAL;
    }

	/* Select the centroids for the remaining clusters. */
	for (cluster = 1; cluster < cluster_num; cluster++) {
 
		/* For each point find its closest distance to any of
		   the previous cluster centers */
		for ( j = 0; j < point_num; j++ ) 
        {
            d = sqrt(cal_dist2(pts[j].pt,centroids[cluster - 1].pt,dim));
			if (d < shortestDistance[j])
            {
                shortestDistance[j] = d;
            }
		}
 
		/* Create an array of the cumulative distances. */
		sum = 0.0;
		for (j = 0; j < point_num; j++) {
			sum += shortestDistance[j];
			cumulativeDistances[j] = sum;
		}

		/* Select a point at random. Those with greater distances
		   have a greater probability of being selected. */
		random = (float) rand() / (float) RAND_MAX * sum;
		selectedIndex = upper_bound(cumulativeDistances, cumulativeDistances+point_num, random) - cumulativeDistances;
 
		/* assign the selected point as the center */
        memcpy(centroids[cluster].pt,pts[selectedIndex].pt,sizeof(double)*dim);
	}

    for (j = 0; j < point_num; j++)
    {
        cid = nearest(pts[j], centroids, cluster_num, dim);
        pts[j].cluster = cid;
        centroids[cid].cluster++;
    }
 
	delete [] shortestDistance;
	delete [] cumulativeDistances;
 
	return;
}


tuple<Point *,Point *> k_means(Point * points, Point * clusters, int point_num, int cluster_num, int maxEpoch, int dim, double tol)
 {
    int cid, i, selectedIndex, tmp_idx;
    double center_shift;
    double ** center_old;
    center_old = new double*[cluster_num];

    for(int i=0;i<cluster_num;i++)
    {
        center_old[i] = new double[dim];
    }
    
    kmeans_plusplus(points, point_num, clusters, cluster_num, dim);

    int it = 0;
    int max_clusters_size = 0;
    int max_cluster_idx = -1;
    vector<int> backup_idx;

    // train model
    for (it =0; it<maxEpoch; it++)
    {
        center_shift = 0.0;
        max_clusters_size = 0;
        backup_idx.clear();
        for (i = 0; i < cluster_num; i++)
        {
            memcpy(center_old[i], clusters[i].pt,sizeof(double)*dim);
            memset(clusters[i].pt,0,sizeof(double)*dim);
            if(clusters[i].cluster>max_clusters_size)
            {
                max_clusters_size=clusters[i].cluster;
                max_cluster_idx = i;
            }
        }

        for (i = 0; i < point_num; i++)
        {
            cid = points[i].cluster;
            if(cid==max_cluster_idx)
            {
                backup_idx.push_back(i);
            }
            // clusters[cid].cluster++;
            for(int d=0; d<dim; d++)
            {
                clusters[cid].pt[d] += points[i].pt[d];
            }
        }

        // check if there is a cluster is null
        for(i = 0; i < cluster_num; i++)
        {
            if(clusters[i].cluster==0)
            {
                /* Pick the empty cluster centroids from the cluster has the largest size. */
                tmp_idx = rand() % max_clusters_size;
	            selectedIndex = backup_idx[tmp_idx];
                memcpy(clusters[i].pt, points[selectedIndex].pt, sizeof(double)*dim);
                clusters[i].cluster++;
                clusters[points[selectedIndex].cluster].cluster--;
                backup_idx.erase(backup_idx.begin()+tmp_idx);
                max_clusters_size--;
            }
        }

        for(i = 0; i < cluster_num; i++)
        {
            // compute the new center
            for(int d=0; d<dim; d++)
            {
                clusters[i].pt[d] /= clusters[i].cluster;
            }
            clusters[i].cluster = 0;
            // compute center shift
            center_shift += sqrt(cal_dist2(center_old[i],clusters[i].pt,dim));
        }

        // clustering results according to latest centers and update points number in a cluster
        for (i = 0; i < point_num; i++)
        {
            cid = nearest(points[i], clusters, cluster_num, dim);
            if(cid!=points[i].cluster)
            {
                points[i].cluster = cid;
            }
            clusters[cid].cluster++;
        }

        if(center_shift<tol)
        {
           break;
        }
    }

    // check if there is a cluster is null
    int flag = 1;
    for(i = 0; i < cluster_num; i++)
    {
        if(clusters[i].cluster==0)
        {
            flag = 0;
        }
        delete [] center_old[i];
    }

    if(flag==0)
    {
        for(i = 0; i < cluster_num; i++)
        {
            clusters[i].cluster=0;
        }
        kmeans_plusplus(points, point_num, clusters, cluster_num, dim);
    }

    delete center_old;
    return {clusters,points};
}

int main()
{ 
    
    /*
    Init
    dataset: distribution_size
    distribution: skewed, uniform, normal, OSM, tiger, default skewed
    B: block size, default 100
    */

    string dataset = Constants::DATASETS;
    string distribution = Constants::DISTRIBUTION;
    int skewness = Constants::SKEWNESS;
    int data_size = Constants::DATASIZE;
    int scale = Constants::SCALE;
    string model_path = Constants::MODEL_PATH;
    check_dir(model_path);
    string kmeans_model_path = Constants::K_MEANS_PATH;
    check_dir(model_path + "/kmeans/");
    check_dir(kmeans_model_path);
    string record_model_path = Constants::MODEL_R_PATH;
    check_dir(model_path + "/record/");
    check_dir(record_model_path);
    string data_path = dataset + "/data/";
    srand(1);
    
    int level_num = Constants::LEVEL;
    int level_cluster[level_num] = {20,20,20,20};
    int dim = Constants::DIM;
    int B = Constants::B;
    int maxEpoch = 300;
    double tolerance = 0.0001;
    int cur_level = 0;
    int partition,i,j,d,cluster_num,tol_cid;
    double x_min = 0.0;
    double y_min = 0.0;
    double x_max = 1.0;
    double y_max = 1.0;
    int min_points;

    cout<<model_path<<endl;
    cout<<record_model_path<<endl;
    cout<<kmeans_model_path<<endl;
    cout<<"scale: "<<scale<<endl;
    cout<<"data_size: "<<data_size<<endl;

    /* Allocate the interior nodes center space in advance */
    Point ** index_centers;
    index_centers = new Point*[level_num-1];
    cluster_num = 1;
    for(i=0; i<level_num-1; i++)
    {
        cluster_num *= level_cluster[i];
        index_centers[i] = new Point[cluster_num];
        for(j=0;j<cluster_num;j++)
        {
            index_centers[i][j].pt = new double[dim];
        }
    }

    /* reading dataset */
    cout<<"reading dataset"<<endl;
    cout<<Constants::INDEX_DIR + distribution +"_"+to_string(data_size)+"_"+to_string(skewness)+"_"+to_string(dim)+"_.csv"<<endl;
    int point_num = getFileLine(Constants::INDEX_DIR + distribution +"_"+to_string(data_size)+"_"+to_string(skewness)+"_"+to_string(dim)+"_.csv");
    double ** dat =  get_2d_points_comma(Constants::INDEX_DIR + distribution +"_"+to_string(data_size)+"_"+to_string(skewness)+"_"+to_string(dim)+"_.csv", data_size);
    Point * pts;
    
    pts = new Point[point_num];
    for(int i=0; i<point_num; i++)
    {
       pts[i].pt = new double[dim];
    }
    cout<<"convert to Point"<<endl;
    generate_points(pts, dat, point_num, dim);
    printf("Dataset size = %d\n",point_num);
    
    /* for root node, do kmeans on whole dataset */
    auto start = chrono::high_resolution_clock::now();
    printf("Training Root node\n");
    cluster_num = level_cluster[cur_level];
    min_points = point_num/cluster_num;
    tuple<Point *,Point*> fit_result = k_means(pts, index_centers[cur_level], point_num, cluster_num, maxEpoch, dim, tolerance);
    index_centers[cur_level] = get<0>(fit_result);
    pts = get<1>(fit_result);
    /* train data in each partition */
    printf("Writing Root node model...\n");
    write_points(index_centers[cur_level], level_cluster[cur_level], dim, kmeans_model_path+"level"+to_string(cur_level)+"model.csv");
    cur_level++;

    /* temporary record*/
    Point * new_pts;
    new_pts = new Point[point_num];
    for(int i=0; i<point_num; i++)
    {
       new_pts[i].pt = new double[dim];
    }

    Point * sub_center;
    Point * sub_pts;
    int sub_size, cur_pos, cur_sub_point;
    sub_size = 0;
    partition = 1;
    tuple<Point *,Point*> sub_fit_result;

    /* for interior node, sub model */
    while(cur_level<level_num-1)
    {
        /* compute submodel number at each level*/
        cur_pos = 0;
        partition *= level_cluster[cur_level-1]; 
        cluster_num = level_cluster[cur_level];
        
        /* train data in each partition */
        for(i=0;i<partition;i++)
        {
            /* gain points of each sub partition first */
            sub_size = index_centers[cur_level-1][i].cluster;
            if(sub_size<=0){continue;}
            printf("partition id = %d, partition size = %d\n",i, sub_size);

            sub_pts = new Point[sub_size];
            sub_center = new Point[cluster_num];
            cur_sub_point = 0;
            
            for(j=0;j<point_num;j++)
            {
                if(pts[j].cluster == i)
                {
                    sub_pts[cur_sub_point].pt = new double[dim];
                    memcpy(sub_pts[cur_sub_point].pt, pts[j].pt, sizeof(double)*dim);
                    cur_sub_point++;
                }
            }

            for(int j=0;j<cluster_num;j++)
            {
                sub_center[j].pt = new double[dim];
            }

            printf("Finish gain points in the partition: sub_pts\n");

            /* train sub model, saving pts result into new pts and new centers, saving cluster result to new cluster */
            sub_fit_result = k_means(sub_pts, sub_center, sub_size, cluster_num, maxEpoch, dim, tolerance);
            printf("Finish K-Means, cluster number = %d\n",cluster_num);
            sub_center = get<0>(sub_fit_result);
            sub_pts = get<1>(sub_fit_result);
            for(j=0;j<sub_size;j++)
            {
                memcpy(new_pts[cur_pos].pt,sub_pts[j].pt,sizeof(double)*dim);
                new_pts[cur_pos].cluster = sub_pts[j].cluster + i*level_cluster[cur_level];
                cur_pos++;
                delete [] sub_pts[j].pt;
            }
            delete sub_pts;

            for(j=0;j<cluster_num;j++)
            {
                tol_cid = i*level_cluster[cur_level]+j;
                index_centers[cur_level][tol_cid].cluster = sub_center[j].cluster;
                memcpy(index_centers[cur_level][tol_cid].pt, sub_center[j].pt, sizeof(double)*dim);
                delete [] sub_center[j].pt;
            }
            delete sub_center;
        }

        /* After complete training a interior level, exchange the new pts and pts result and write the new center result*/
        printf("writing partition result...\n");
        write_points(index_centers[cur_level], partition*level_cluster[cur_level], dim, kmeans_model_path+"level"+to_string(cur_level)+"model.csv");
        for(i=0; i<point_num; i++)
        {
            memcpy(pts[i].pt, new_pts[i].pt, sizeof(double)*dim);
            pts[i].cluster = new_pts[i].cluster;
        }
        printf("Complete the level = %d\n", cur_level);
        printf("Recording the points result for interior level\n");
        cur_level++;
    }

    /* Traing and record index info at leaf node */
    partition *= level_cluster[cur_level-1];
    cur_pos = 0;
    int max_cluster_num = level_cluster[cur_level];

    /* record info */
    vector<int> block_pts;
    vector<int> part_info;
    int bid = 0;
    double ** split_pts;
    int ** bid_list;
    double *** cluster_info;
    cluster_info = new double**[partition];

    /* temporary various */
    int p, pos,block_num,cluster_pts;
    Point single_cluster;
    double * cluster_center;
    double R;
    int * block_num_list;
    vector<double> R_list;
    int * cluster_num_list;
    cluster_num_list = new int[partition];

    for(i=0;i<partition;i++)
    {
        /* gain points of each sub partition first */
        sub_size = index_centers[cur_level-1][i].cluster;
        if(sub_size!=0){part_info.push_back(i);}
        if(sub_size<=0){continue;}
        cluster_num = (sub_size/B+1)<max_cluster_num?(sub_size/B+1):max_cluster_num;
        cluster_num_list[i] = cluster_num;

        sub_pts = new Point[sub_size];
        sub_center = new Point[cluster_num];
        split_pts = new double*[cluster_num];
        bid_list = new int*[cluster_num];
        cluster_info[i] = new double*[cluster_num];
        block_num_list = new int[cluster_num];
        cur_sub_point = 0;
        
        for(j=0;j<point_num;j++)
        {
            if(pts[j].cluster == i)
            {
                sub_pts[cur_sub_point].pt = new double[dim];
                memcpy(sub_pts[cur_sub_point].pt, pts[j].pt, sizeof(double)*dim);
                cur_sub_point++;
            }
        }

        for(int j=0;j<cluster_num;j++)
        {
            sub_center[j].pt = new double[dim];
        }

        /* train sub model, saving pts result into new pts and new centers, saving cluster result to new cluster */
        sub_fit_result = k_means(sub_pts, sub_center, sub_size, cluster_num, maxEpoch, dim, tolerance);
        sub_center = get<0>(sub_fit_result);
        sub_pts = get<1>(sub_fit_result);

        /*
        Check each cluster
        Compute the R, record the cluster info(center R)
        Devided into blocks if cluster is larger than B
        Record the block id, split point and block size
        */
       for(j=0; j<cluster_num; j++)
       {
           single_cluster = sub_center[j];
           cluster_pts = single_cluster.cluster;
           cluster_center = single_cluster.pt;
           Point * ptsinc; // saving points in this cluster
           double * dist_list; // saving dist between pts and center
           int * idx;
           ptsinc = new Point[cluster_pts];
           dist_list = new double[cluster_pts];
           double x[cluster_pts];
           idx = new int[cluster_pts];
           pos = 0;
           cluster_info[i][j] = new double[dim+1];
           memcpy(cluster_info[i][j],cluster_center,sizeof(double)*dim);

           for(p = 0; p<sub_size;p++)
           {
               if(sub_pts[p].cluster == j)
               {
                   ptsinc[pos].pt = new double[dim];
                   memcpy(ptsinc[pos].pt, sub_pts[p].pt, sizeof(double)*dim);
                   x[pos] = sub_pts[p].pt[0];
                   idx[pos] = pos;
                   dist_list[pos] = sqrt(cal_dist2(sub_pts[p].pt,cluster_center,2));
                   pos++;
               }
           }

           R = dist_list[argmax(dist_list,dist_list+pos)];
           cluster_info[i][j][dim] = R;
           R_list.push_back(R);
           argsort(x,cluster_pts,idx);

           if(pos%100==0)
           {
               block_num = pos/B;
           }
           else
           {
               block_num = pos/B+1;
           }
           
           bid_list[j] = new int[block_num];
           split_pts[j] = new double[block_num];
           block_num_list[j] = block_num;
           for(p=0; p<block_num; p++)
           {
               bid_list[j][p] = bid;
               if(p==block_num-1){        
                   for(int n = 0; n<(cluster_pts - p*B); n++)
                   {
                       memcpy(new_pts[cur_pos].pt,ptsinc[idx[p*B+n]].pt,sizeof(double)*dim);
                       if(n==0){split_pts[j][p] = new_pts[cur_pos].pt[0];}
                       new_pts[cur_pos].cluster = bid;
                       cur_pos++;
                       delete [] ptsinc[idx[p*B+n]].pt;
                   }
                   block_pts.push_back(cluster_pts - p*B);
               }
               else{
                   for(int n = 0; n<B; n++)
                   {
                       memcpy(new_pts[cur_pos].pt,ptsinc[idx[p*B+n]].pt,sizeof(double)*dim);
                       if(n==0){split_pts[j][p] = new_pts[cur_pos].pt[0];}
                       new_pts[cur_pos].cluster = bid;
                       cur_pos++;
                       delete [] ptsinc[idx[p*B+n]].pt;
                   }
                   block_pts.push_back(B);
               }
               bid++;
           }

           delete ptsinc;
           delete [] dist_list;
           delete [] idx;
        }

        /*
        Write result to docments:cluster info, cluster block, split points
        */
       check_dir(record_model_path + to_string(i) + "/");
       check_dir(record_model_path + to_string(i) + "/" + "original/");

       write_double_points(cluster_info[i], cluster_num, dim+1, record_model_path + to_string(i) + "/" + "original/cluster_info.csv");
       write_split_pts(split_pts, cluster_num, block_num_list, record_model_path + to_string(i) + "/" + "original/split_pts.csv");
       write_cluster_blk(bid_list, cluster_num, block_num_list, record_model_path + to_string(i) + "/" + "original/cluster_block.csv");

        for(j=0;j<sub_size;j++)
        {
            delete [] sub_pts[j].pt;
        }
        delete sub_pts;

        for(j=0;j<cluster_num;j++)
        {
            delete [] sub_center[j].pt;
            delete [] split_pts[j];
            delete [] bid_list[j];
        }

        delete sub_center;
        delete split_pts;
        delete bid_list;
        delete [] block_num_list;
    }

    /* Write blocks(new pts) and block counts*/
    write_points(new_pts, point_num, dim, Constants::INDEX_DIR+"index_points"+to_string(cur_level)+".csv");
    write_int_list(block_pts, block_pts.size(), record_model_path+"block_count.csv");
    write_int_list(part_info, part_info.size(), record_model_path+"part_info.csv");

    /*
    Design the grid
    -R_list used to compute the grid edge size
    -Scan the partition again to record grid information
    */
    int c_num = R_list.size();
    sort(R_list.begin(),R_list.end());
    double edge_len;
    cout<<c_num*0.95<<endl;
    edge_len = R_list[c_num*0.95];
    double x_range,y_range,gap_x,gap_y;
    x_range = x_max-x_min+0.001;
    y_range = y_max-y_min+0.001;
    int split_num_x,split_num_y;
    split_num_x = x_range/edge_len;
    split_num_y = y_range/edge_len;
    gap_x = x_range/split_num_x;
    gap_y = y_range/split_num_y;

    string grid_info_path;
    grid_info_path = record_model_path+"grid_info.csv";
    ofstream outfile;
    outfile.open(grid_info_path,ios::out);
    outfile<<split_num_x<<endl;
    outfile<<split_num_y<<endl;
    outfile<<setprecision(16)<<gap_x<<endl;
    outfile<<setprecision(16)<<gap_y<<endl;
    outfile<<setprecision(16)<<x_min<<endl;
    outfile<<setprecision(16)<<y_min<<endl;
    outfile<<setprecision(16)<<x_max<<endl;
    outfile<<setprecision(16)<<y_max<<endl;
    outfile.close();

    vector<vector<int>> grid_cluster;
    int cell_num;
    cell_num = split_num_x*split_num_y;
    grid_cluster.resize(cell_num);

    double x_begin_val,x_end_val,y_begin_val,y_end_val;
    int x_begin_gid,x_end_gid,y_begin_gid,y_end_gid,grid_id;
    double center_x;
    double center_y;
    int ** grid_id_range; //xbegin,xend,ybegin,yend
    

    for(i=0; i<partition; i++)
    {
        grid_id_range = new int*[cluster_num_list[i]];
        for(j=0; j<cluster_num_list[i];j++)
        {
            grid_id_range[j] = new int[dim*2];
            center_x = cluster_info[i][j][0];
            center_y = cluster_info[i][j][1];
            R = cluster_info[i][j][2];
            x_begin_val = center_x - R - x_min;
            x_begin_val = x_begin_val<0?0:x_begin_val;
            x_end_val = center_x + R - x_min;
            x_end_val = x_end_val<(x_range-0.001)?x_end_val:(x_range-0.001);
            x_begin_gid = x_begin_val/gap_x;
            x_end_gid = x_end_val/gap_x;
            grid_id_range[j][0] = x_begin_gid;
            grid_id_range[j][1] = x_end_gid;

            y_begin_val = center_y - R - y_min;
            y_begin_val = y_begin_val<0?0:y_begin_val;
            y_end_val = center_y + R - y_min;
            y_end_val = y_end_val<(y_range-0.001)?y_end_val:(y_range-0.001);
            y_begin_gid = y_begin_val/gap_y;
            y_end_gid = y_end_val/gap_y;
            grid_id_range[j][2] = y_begin_gid;
            grid_id_range[j][3] = y_end_gid;
            
            for(int m = x_begin_gid; m<x_end_gid+1; m++)
            {
                for(int n = y_begin_gid; n<y_end_gid+1; n++)
                {
                    grid_id = m * split_num_y + n;
                    grid_cluster[grid_id].push_back(i+j*scale);
                }
            }
        }
        // writing grid id range
        write_grid_ids(grid_id_range, cluster_num_list[i], dim*2, record_model_path + to_string(i) + "/" + "original/grid_ids.csv");

        for(j=0; j<cluster_num_list[i];j++)
        {
            delete [] grid_id_range[j];
        }
        delete grid_id_range;
    }

    //writing grid_cluster
    write_grid_cluster(grid_cluster, record_model_path+"grid_cluster.csv");

    auto end = chrono::high_resolution_clock::now();
    long fit_time;
    fit_time = chrono::duration_cast<chrono::nanoseconds>(end - start).count();
    cout<<fit_time<<endl;
}
