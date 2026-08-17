#include <vector>
#include <string>
#include <random>
#include <memory>
#include <unordered_map>
#include <cmath>
#include <Rcpp.h>
// #include <RcppParallel.h>

using namespace Rcpp;
// using namespace RcppParallel;


const double EM_CONSTANT = 0.5772156649; // Euler-Mascheroni constant
std::mt19937 gen(42); // Fixed random seed for reproducibility

// Define a structure for an Isolation Tree Node
struct IsolationTreeNode {
  bool is_leaf;
  int size;
  std::string feature;
  double value;
  bool is_missing_split;
  std::unique_ptr<IsolationTreeNode> left;
  std::unique_ptr<IsolationTreeNode> right;
  
  IsolationTreeNode(int s) : is_leaf(true), size(s), value(0), is_missing_split(false) {}
  IsolationTreeNode(
    std::string feat, 
    double val, 
    bool missing_split, 
    std::unique_ptr<IsolationTreeNode> l, 
    std::unique_ptr<IsolationTreeNode> r)
    : is_leaf(false), size(0), 
      feature(std::move(feat)), 
      value(val), is_missing_split(missing_split), 
      left(std::move(l)), right(std::move(r)) {}
};

// Convert R DataFrame to C++ data structure
std::vector<std::unordered_map<std::string, double>> convert_dataframe(const DataFrame& df) {
    std::vector<std::unordered_map<std::string, double>> data;
    int n_rows = df.nrows();
    CharacterVector col_names = df.names();
    int n_cols = col_names.size();
    
    for (int i = 0; i < n_rows; ++i) {
      std::unordered_map<std::string, double> row;
      for (int j = 0; j < n_cols; ++j) {
        NumericVector col = df[j];
        row[as<std::string>(col_names[j])] = col[i];
      }
      data.push_back(row);
    }
    return data;
}

// Function to create an isolation tree
std::unique_ptr<IsolationTreeNode> isolation_tree(
    const std::vector<std::unordered_map<std::string, double>>& data, 
    int depth = 0, int max_depth = -1) {
  
  int n = data.size();
  if (max_depth == -1) {
    max_depth = static_cast<int>(log2(n));
  }
  if (n <= 1 || depth >= max_depth) {
    return std::make_unique<IsolationTreeNode>(n);
  }
  
  // Create vector of features to split on
  std::vector<std::string> features;
  for (const auto& [key, _] : data[0]) {
    features.push_back(key);
  }
  
  std::uniform_int_distribution<> feature_dist(0, features.size() - 1);
  std::string split_feature = features[feature_dist(gen)];
  
  // Can split on numeric or missing value
  double min_val = 0.0;
  double max_val = 0.0;
  bool has_missing = false;
  bool has_valid = false;
  for (const auto& row : data) {
    double val = row.at(split_feature);
    if (std::isnan(val)) {
      has_missing = true;
      continue;
    }
    if (!has_valid) {
      min_val = val;
      max_val = val;
      has_valid = true;
    } else {
      min_val = std::min(min_val, val);
      max_val = std::max(max_val, val);
    }
  }

  if (!has_valid || (min_val == max_val && !has_missing)) {
    return std::make_unique<IsolationTreeNode>(n);
  }

  // TODO: Chance to chose missing is 50/50. Could make less likely. Test
  bool is_missing_split = false;
  double split_value = 0.0;
  if (has_missing && (std::bernoulli_distribution(0.5)(gen) || min_val == max_val)) {
    is_missing_split = true;
  } else {
    std::uniform_real_distribution<> split_dist(min_val, max_val);
    split_value = split_dist(gen);
  }
  
  std::vector<std::unordered_map<std::string, double>> left_data, right_data;
  for (const auto& row : data) {
    if (is_missing_split) {
      if (std::isnan(row.at(split_feature))) {
        left_data.push_back(row);
      } else {
        right_data.push_back(row);
      }
    } else {
      if (row.at(split_feature) < split_value) {
        left_data.push_back(row);
      } else {
        right_data.push_back(row);
      }
    }
  }
  
  return std::make_unique<IsolationTreeNode>(
    split_feature, split_value, is_missing_split,
    isolation_tree(left_data, depth + 1, max_depth),
    isolation_tree(right_data, depth + 1, max_depth)
  );
}

// Train an isolation forest
std::vector<std::unique_ptr<IsolationTreeNode>> isolation_forest(
        const std::vector<std::unordered_map<std::string, double>>& data, 
        int n_trees, 
        int max_depth
) {
    std::vector<std::unique_ptr<IsolationTreeNode>> forest;
    forest.reserve(n_trees);
    
    for (int i = 0; i < n_trees; ++i) {
        forest.push_back(isolation_tree(data, 0, max_depth));
    }
    
    return forest;
}

// Compute path length of a single observation in a tree
int path_length(const IsolationTreeNode* node, const std::unordered_map<std::string, double>& obs, int depth = 0) {
    if (node->is_leaf) {
        return depth + node->size;
    }
    
    auto it = obs.find(node->feature);
    if (it == obs.end()) {
        return depth; // If feature is missing, return current depth
    }
    
    if (it->second < node->value) {
        return path_length(node->left.get(), obs, depth + 1);
    } else {
        return path_length(node->right.get(), obs, depth + 1);
    }
}

// Compute anomaly scores for data
// [[Rcpp::export]]
DataFrame calculate_anomaly_score(
        DataFrame df, 
        int n_trees, 
        int max_depth
) {
    
    std::vector<std::unordered_map<std::string, double>> data = convert_dataframe(df);
    
    int n = data.size();
    
    // Generate the forest
    auto forest = isolation_forest(data, n_trees, max_depth);
    
    // Compute average path lengths
    std::vector<double> avg_path_length(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double total_path_length = 0.0;
        for (const auto& tree : forest) {
            total_path_length += path_length(tree.get(), data[i]);
        }
        avg_path_length[i] = total_path_length / n_trees;
    }
    
    // Compute anomaly scores
    double c_n = 2 * (log(n - 1) + EM_CONSTANT) - (2 * (n - 1) / n);
    std::vector<double> scores(n);
    for (int i = 0; i < n; ++i) {
        scores[i] = std::pow(2, -(avg_path_length[i] / c_n));
    }
    
    // Return DataFrame with results
    return DataFrame::create(
        _["avg_depth"] = avg_path_length,
        _["anomaly_score"] = scores
    );
}

// #include <vector>
// #include <string>
// #include <random>
// #include <memory>
// #include <unordered_map>
// #include <cmath>
// #include <Rcpp.h>

// using namespace Rcpp;

// const double EM_CONSTANT = 0.5772156649; // Euler-Mascheroni constant
// std::mt19937 gen(42); // Fixed random seed for reproducibility

// struct IsolationTreeNode {
//   bool is_leaf;
//   int size;
//   std::string feature;
//   double value;
//   bool is_missing_split;
//   std::unique_ptr<IsolationTreeNode> left;
//   std::unique_ptr<IsolationTreeNode> right;

//   IsolationTreeNode(int s)
//     : is_leaf(true),
//       size(s),
//       feature(""),
//       value(0.0),
//       is_missing_split(false),
//       left(nullptr),
//       right(nullptr) {}

//   IsolationTreeNode(
//     std::string feat,
//     double val,
//     bool missing_split,
//     std::unique_ptr<IsolationTreeNode> l,
//     std::unique_ptr<IsolationTreeNode> r
//   )
//     : is_leaf(false),
//       size(0),
//       feature(std::move(feat)),
//       value(val),
//       is_missing_split(missing_split),
//       left(std::move(l)),
//       right(std::move(r)) {}
// };

// std::vector<std::unordered_map<std::string, double>>
// convert_dataframe(const DataFrame& df) {
//   std::vector<std::unordered_map<std::string, double>> data;

//   int n_rows = df.nrows();
//   CharacterVector col_names = df.names();
//   int n_cols = col_names.size();

//   data.reserve(n_rows);

//   for (int i = 0; i < n_rows; ++i) {
//     std::unordered_map<std::string, double> row;
//     row.reserve(n_cols);

//     for (int j = 0; j < n_cols; ++j) {
//       NumericVector col = df[j];
//       row[as<std::string>(col_names[j])] = col[i];
//     }
//     data.push_back(std::move(row));
//   }

//   return data;
// }

// std::unique_ptr<IsolationTreeNode>
// isolation_tree(
//   const std::vector<std::unordered_map<std::string, double>>& data,
//   int depth = 0,
//   int max_depth = -1
// ) {
//   int n = static_cast<int>(data.size());

//   if (max_depth == -1) {
//     max_depth = static_cast<int>(std::log2(n));
//   }

//   if (n <= 1 || depth >= max_depth) {
//     return std::make_unique<IsolationTreeNode>(n);
//   }

//   std::vector<std::string> features;
//   features.reserve(data[0].size());

//   for (const auto& kv : data[0]) {
//     features.push_back(kv.first);
//   }

//   std::uniform_int_distribution<> feat_dist(0, features.size() - 1);
//   std::string split_feature = features[feat_dist(gen)];

//   double min_val = data[0].at(split_feature);
//   double max_val = min_val;
//   int missing_count = 0;

//   for (const auto& row : data) {
//     double val = row.at(split_feature);

//     if (std::isnan(val)) {
//       missing_count++;
//       continue;
//     }

//     min_val = std::min(min_val, val);
//     max_val = std::max(max_val, val);
//   }

//   bool has_missing = (missing_count > 0);

//   double missing_fraction =
//     static_cast<double>(missing_count) / static_cast<double>(n);

//   if (min_val == max_val && !has_missing) {
//     return std::make_unique<IsolationTreeNode>(n);
//   }

//   // Calibrated missing split probability:
//   // - rare missing -> rare split (more anomalous)
//   // - common missing -> more missing splits (less anomalous)
//   double p_miss = std::min(0.5, missing_fraction);

//   bool is_missing_split = false;
//   double split_value = 0.0;

//   if (has_missing &&
//       (std::bernoulli_distribution(p_miss)(gen) || min_val == max_val)) {
//     is_missing_split = true;
//   } else {
//     std::uniform_real_distribution<> split_dist(min_val, max_val);
//     split_value = split_dist(gen);
//   }

//   std::vector<std::unordered_map<std::string, double>> left_data;
//   std::vector<std::unordered_map<std::string, double>> right_data;
//   left_data.reserve(n);
//   right_data.reserve(n);

//   for (const auto& row : data) {
//     double v = row.at(split_feature);

//     if (is_missing_split) {
//       if (std::isnan(v)) {
//         left_data.push_back(row);
//       } else {
//         right_data.push_back(row);
//       }
//     } else {
//       if (v < split_value) {
//         left_data.push_back(row);
//       } else {
//         right_data.push_back(row);
//       }
//     }
//   }

//   return std::make_unique<IsolationTreeNode>(
//     split_feature,
//     split_value,
//     is_missing_split,
//     isolation_tree(left_data, depth + 1, max_depth),
//     isolation_tree(right_data, depth + 1, max_depth)
//   );
// }

// std::vector<std::unique_ptr<IsolationTreeNode>>
// isolation_forest(
//   const std::vector<std::unordered_map<std::string, double>>& data,
//   int n_trees,
//   int max_depth
// ) {
//   std::vector<std::unique_ptr<IsolationTreeNode>> forest;
//   forest.reserve(n_trees);

//   for (int i = 0; i < n_trees; ++i) {
//     forest.push_back(isolation_tree(data, 0, max_depth));
//   }

//   return forest;
// }

// // IMPORTANT FIX: respect is_missing_split during traversal
// int path_length(
//   const IsolationTreeNode* node,
//   const std::unordered_map<std::string, double>& obs,
//   int depth = 0
// ) {
//   if (node->is_leaf) {
//     return depth + node->size;
//   }

//   auto it = obs.find(node->feature);
//   if (it == obs.end()) {
//     // Feature absent from obs map. Fall back: stop at current depth.
//     return depth;
//   }

//   double x = it->second;

//   if (node->is_missing_split) {
//     if (std::isnan(x)) {
//       return path_length(node->left.get(), obs, depth + 1);
//     }
//     return path_length(node->right.get(), obs, depth + 1);
//   }

//   if (std::isnan(x)) {
//     // Numeric split but x missing: choose a side deterministically.
//     // Here: treat missing as going left.
//     return path_length(node->left.get(), obs, depth + 1);
//   }

//   if (x < node->value) {
//     return path_length(node->left.get(), obs, depth + 1);
//   }

//   return path_length(node->right.get(), obs, depth + 1);
// }

// // [[Rcpp::export]]
// DataFrame calculate_anomaly_score(DataFrame df, int n_trees, int max_depth) {
//   std::vector<std::unordered_map<std::string, double>> data =
//     convert_dataframe(df);

//   int n = static_cast<int>(data.size());

//   auto forest = isolation_forest(data, n_trees, max_depth);

//   std::vector<double> avg_path_length(n, 0.0);

//   for (int i = 0; i < n; ++i) {
//     double total_path_length = 0.0;

//     for (const auto& tree : forest) {
//       total_path_length += path_length(tree.get(), data[i]);
//     }

//     avg_path_length[i] = total_path_length / static_cast<double>(n_trees);
//   }

//   double c_n = 2.0 * (std::log(n - 1.0) + EM_CONSTANT) -
//     (2.0 * (n - 1.0) / n);

//   std::vector<double> scores(n, 0.0);

//   for (int i = 0; i < n; ++i) {
//     scores[i] = std::pow(2.0, -(avg_path_length[i] / c_n));
//   }

//   return DataFrame::create(
//     _["avg_depth"] = avg_path_length,
//     _["anomaly_score"] = scores
//   );
// }


// TODO: Parallelization would be faster within RcppParallel framework
// but requires setup and testing.
// Parallel task struct to process each data frame
// struct ParallelTask : public Worker {
//     // Input list of data frames
//     const List& df_list;
//     int n_trees;
//     int max_depth;
//     // Output list to store results
//     List& result;
// 
//     // Constructor to initialize the worker
//     ParallelTask(const Rcpp::List& df_list, int n_trees, int max_depth, Rcpp::List& result)
//         : df_list(df_list), n_trees(n_trees), max_depth(max_depth), result(result) {}
//     
//     
//     // Function to process each individual task in parallel
//     void operator()(std::size_t begin, std::size_t end) {
//         for (std::size_t i = begin; i < end; ++i) {
//             DataFrame df = Rcpp::as<DataFrame>(df_list[i]);
//             std::vector<std::unordered_map<std::string, double>> data = convert_dataframe(df);
//             Rcpp::DataFrame scores = calculate_anomaly_score(data, 
//                                                              n_trees, 
//                                                              max_depth);
//             result[i] = scores;
//         }
//     }
// };
// Parallelized function to process a vector of data frames
// List train_anomaly_model(List df_list, int n_trees = 100, 
//                          int max_depth = -1, int n_cores = 1) {
//     int n = df_list.size();
//     List result(n);
//     
//     // Create the ParallelTask object
//     ParallelTask task(df_list, n_trees, max_depth, result);
//     
//     // Run the parallel task using parallelFor
//     parallelFor(0, n, task);
//     
//     return result;
// }