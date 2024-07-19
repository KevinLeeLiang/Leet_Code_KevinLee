//
// Created by garen_lee on 2024/5/11.
/**
  ******************************************************************************
  * @file           : SummaryRanges.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/11
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_SUMMARYRANGES_H
#define CTEST_PRO_SUMMARYRANGES_H

#include <vector>
#include <map>

using namespace std;

class SummaryRanges {
private:
    map<int, int> intervals; // key为区间left.value为区间right
public:
    SummaryRanges() {

    }

    void addNum(int value) {
        int left = value, right = value;
        auto it = intervals.upper_bound(value);
        if (it != intervals.begin() && (--it)->second + 1 < value) {
            ++it;
        }
        while (it != intervals.end() && value + 1 >= it->first) {
            left = min(left, it->first);
            right = max(right, it->second);
            it = intervals.erase(it);
        }
        intervals[left] = right;
    }

    vector<vector<int>> getIntervals() {
        vector<vector<int>> result;
        for (const auto &pair : intervals) {
            result.push_back({pair.first, pair.second});
        }
        return result;
    }
};


#endif//CTEST_PRO_SUMMARYRANGES_H
