//
// Created by garen_lee on 2024/5/14.
/**
  ******************************************************************************
  * @file           : largestDivisibleSubset.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/14
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_LARGESTDIVISIBLESUBSET_H
#define CTEST_PRO_LARGESTDIVISIBLESUBSET_H

#include <algorithm>
#include <vector>

using namespace std;
namespace largestDivisibleSubset {
    class solution {
    public:
        vector<int> largestDivisibleSubset(vector<int> &nums) {
            int n = nums.size();
            if (n == 0) return {};

            sort(nums.begin(), nums.end());

            vector<int> dp(n, 1); // dp[i] 表示以 nums[i] 结尾的最大整除子集的长度
            vector<int> prev(n, -1); // prev[i] 表示在以 nums[i] 结尾的最大整除子集中，nums[i] 前面的元素的下标

            int maxLen = 1, maxIdx = 0;
            for (int i = 1; i < n; ++i) {
                for (int j = 0; j < i; ++j) {
                    if (nums[i] % nums[j] == 0 && dp[i] < dp[j] + 1) {
                        dp[i] = dp[j] + 1;
                        prev[i] = j;
                    }
                }
                if (dp[i] > maxLen) {
                    maxLen = dp[i];
                    maxIdx = i;
                }
            }

            vector<int> result;
            while (maxIdx != -1) {
                result.push_back(nums[maxIdx]);
                maxIdx = prev[maxIdx];
            }

            reverse(result.begin(), result.end()); // 最大整除子集是递增序列，需要反转
            return result;
        }
    };
}

#endif//CTEST_PRO_LARGESTDIVISIBLESUBSET_H
