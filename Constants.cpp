#include "Constants.h"

const string Constants::DATASETS = "skewed_16m";
const string Constants::DISTRIBUTION = "skewed";
const string Constants::INSERT_DISTRIBUTION = "normal";

const string Constants::MODEL_PATH = "model/";
const string Constants::DATA_PATH = Constants::DATASETS + "/data/";
const string Constants::MODEL_R_PATH = "model/record/"+Constants::DATASETS+"/";

const string Constants::INDEX_DIR = "dataset/" + Constants::DATASETS +  "/";
const string Constants::INDEX_PATH = Constants::INDEX_DIR + "index_points"+to_string(Constants::LEVEL-1)+".csv";

const string Constants::WORKLOAD_PATH = "workload/";
const string Constants::WINDOW_WORKLOAD_PATH = Constants::WORKLOAD_PATH + "window/";
const string Constants::KNN_WORKLOAD_PATH = Constants::WORKLOAD_PATH + "knn/";
const string Constants::INSERT_WORKLOAD_PATH = Constants::WORKLOAD_PATH + "insert/";


const string Constants::K_MEANS_PATH = "model/kmeans/" + Constants::DATASETS + "/";
const string Constants::BLOCK_INFO = Constants::MODEL_R_PATH + "block_count.csv";
const string Constants::CLUSTER_BLOCK = "cluster_block.csv";
const string Constants::SPLIT_PTS = "split_pts.csv";
const string Constants::CLUSTER_INFO = "cluster_info.csv";

const string Constants::GRID_CLUSTER = Constants::MODEL_R_PATH + "grid_cluster.csv";
const string Constants::GRID_INFO = Constants::MODEL_R_PATH + "grid_info.csv";
const string Constants::PART_INFO = Constants::MODEL_R_PATH + "part_info.csv";

Constants::Constants()
{
}
