# FHSIE
Official code for the paper 'A Fast Hybrid Spatial Index with External Memory Support (FHSIE)'
https://ieeexplore.ieee.org/abstract/document/10148115

## How to use
1. boost
homepage: https://www.boost.org/

2. Change path
change {boost_1_78_0/include's root} to your own path


## Parameter Settings 
16M: level_num=4; level_custer[level_num] = {20,20,20,20};
256M: level_num=4; level_custer[level_num] = {40,40,40,40};

## Running code
{reaplce /home/suzy/boost_1_78_0/include with your own boost_1_78_0/include's root}
{reaplce /home/suzy/boost_1_78_0/lib with your own boost_1_78_0/lib's root}
1. Runing the utils.cpp and Constant.cpp
g++ -O3 -std=c++14 -o utils.o -c utils.cpp -g -I/home/suzy/boost_1_78_0/include
g++ -O3 -std=c++14 -o Constants.o -c Constants.cpp -g -I/home/suzy/boost_1_78_0/include

2. Running AKS to build the embedding, we upload an index on skewed_16m, if you just want to check the query, you could skip this step
g++ -O3 -std=c++14 -o AKS.o -c AKS.cpp -g -I/home/suzy/boost_1_78_0/include
g++ -O3 -std=c++14 -o AKS AKS.o Constants.o utils.o -I/home/suzy/boost_1_78_0/include -lpthread -L/home/suzy/boost_1_78_0/lib -lboost_filesystem
./AKS

3. Running {query type}_query.cpp to conduct query
An example of point query:
g++ -O3 -std=c++14 -o point_query.o -c point_query.cpp -g -I/home/suzy/boost_1_78_0/include
g++ -O3 -std=c++14 -o point_query point_query.o Constants.o utils.o -I/home/suzy/boost_1_78_0/include -lpthread -L/home/suzy/boost_1_78_0/lib -lboost_filesystem
./point_query

For window and knn query, we replace the point_query with window_query and knn_query, repectively.

## Citation
@inproceedings{FHSIE,
  title={A Fast Hybrid Spatial Index with External Memory Support},
  author={Su, Xinyu and Qi, Jianzhong and Tanin, Egemen},
  booktitle={2023 IEEE 39th International Conference on Data Engineering Workshops (ICDEW)},
  pages={67--73},
  year={2023},
  organization={IEEE}
}
