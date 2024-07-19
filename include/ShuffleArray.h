//
// Created by garen_lee on 2024/5/21.
/**
  ******************************************************************************
  * @file           : ShuffleArray.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/21
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_SHUFFLEARRAY_H
#define CTEST_PRO_SHUFFLEARRAY_H
#include <vector>
#include <random>
namespace Shuffle {
  using namespace std;
  class ShuffleArray {
  private:
    vector<int>original_;
  public:
    ShuffleArray(vector<int>& nums) {
      original_ = nums;
    }

    vector<int> reset() {
      return original_;
    }

    vector<int> shuffle() {
      std::vector<int> shuffled = original_;
      std::random_device rd;
      std::mt19937 g(rd());
      for (int i = shuffled.size() - 1; i > 0; --i) {
        std::uniform_int_distribution<int> distribution(0, i);
        int j = distribution(g);
        std::swap(shuffled[i], shuffled[j]);
      }
      return shuffled;
    }
  };
}

#endif//CTEST_PRO_SHUFFLEARRAY_H
