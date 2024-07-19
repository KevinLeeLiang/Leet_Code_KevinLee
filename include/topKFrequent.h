//
// Created by garen_lee on 2024/5/9.
/**
  ******************************************************************************
  * @file           : topKFrequent_solution.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/9
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_TOPKFREQUENT_SOLUTION_H
#define CTEST_PRO_TOPKFREQUENT_SOLUTION_H

#include<vector>
#include<unordered_map>
#include<queue>
using namespace std;
class topKFrequent_solution {
public:
  vector<int> topKFrequent(vector<int>& nums, int k) {
    unordered_map<int, int> map;
    for (auto num : nums){
      map[num]++;
    }
    auto cmp = [](pair<int, int>&a, pair<int, int>& b) {
      return a.second > b.second;
    };
    priority_queue<pair<int, int>, vector<pair<int, int>>, decltype(cmp)> pq(cmp);
    for (auto pair : map) {
      pq.push(pair);
      if (pq.size() > k){
        pq.pop();
      }
    }

    vector<int>result;
    while (!pq.empty()) {
      result.push_back(pq.top().first);
      pq.pop();
    }
    return result;
  }
};


#endif//CTEST_PRO_TOPKFREQUENT_SOLUTION_H
