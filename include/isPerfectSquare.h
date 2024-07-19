//
// Created by garen_lee on 2024/5/14.
/**
  ******************************************************************************
  * @file           : isPerfectSquare.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/14
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_ISPERFECTSQUARE_H
#define CTEST_PRO_ISPERFECTSQUARE_H

namespace isPerfectSquare {
    class solution {
    public:
        bool isPerfectSquare(int num) {
            // 特殊情况处理
            if (num < 2) {
                return true;
            }

            // 使用二分查找来寻找平方根
            long left = 2, right = num / 2;
            while (left <= right) {
                long mid = left + (right - left) / 2;
                long guess = mid * mid;
                if (guess == num) {
                    return true;
                } else if (guess < num) {
                    left = mid + 1;
                } else {
                    right = mid - 1;
                }
            }

            return false;
        }
    };
}

#endif//CTEST_PRO_ISPERFECTSQUARE_H
