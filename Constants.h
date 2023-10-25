#ifndef CONSTANTS_H
#define CONSTANTS_H
#include <string>
using namespace std;
class Constants
{
public:
    static const int DIM = 2;
    static const int B = 100;
    static const int LEVEL = 4;
    static const int SCALE = 10000;
    static const int SKEWNESS = 4;
    static const int DATASIZE = 16000000;

    static const string DATASETS;
    static const string DISTRIBUTION;
    static const string INSERT_DISTRIBUTION;

    static const string MODEL_PATH;
    static const string DATA_PATH;
    static const string MODEL_R_PATH;

    static const string WORKLOAD_PATH;
    static const string WINDOW_WORKLOAD_PATH;
    static const string KNN_WORKLOAD_PATH;
    static const string INSERT_WORKLOAD_PATH;

    static const string INDEX_DIR;
    static const string INDEX_PATH;

    static const string K_MEANS_PATH;
    static const string BLOCK_INFO;
    static const string CLUSTER_BLOCK;
    static const string SPLIT_PTS;
    static const string CLUSTER_INFO;

    static const string GRID_CLUSTER;
    static const string GRID_INFO;
    static const string PART_INFO;

    Constants();
};

#endif
