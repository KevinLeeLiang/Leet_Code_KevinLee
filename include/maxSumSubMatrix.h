//
// Created by garen_lee on 2024/5/13.
/**
  ******************************************************************************
  * @file           : maxSumSubMatrix.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/13
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_MAXSUMSUBMATRIX_H
#define CTEST_PRO_MAXSUMSUBMATRIX_H

#include <vector>
#include <climits>
#include <algorithm>
#include <set>

using namespace std;
namespace maxSumSubMatrix {
    class Solution {
    public:
        int maxSumSubmatrix(vector<vector<int>> &matrix, int k) {
            int ans = INT_MIN;
            int m = matrix.size(), n = matrix[0].size();
            for (int i = 0; i < m; ++i) {
                vector<int> sum(n);
                for (int j = i; j < m; ++j) {
                    for (int c = 0; c < n; ++c) {
                        sum[c] += matrix[j][c];
                    }
                    set<int> sum_set{0};
                    int s = 0;
                    for (int v : sum) {
                        s += v;
                        auto lb = sum_set.lower_bound(s - k);
                        if (lb != sum_set.end()) {
                            ans = max(ans, s - *lb);
                        }
                        sum_set.insert(s);
                    }
                }
            }
            return ans;
        }
    };
}

#endif//CTEST_PRO_MAXSUMSUBMATRIX_H
