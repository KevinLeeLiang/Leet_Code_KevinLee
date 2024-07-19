//
// Created by garen_lee on 2024/5/9.
/**
  ******************************************************************************
  * @file           : reverseVowels.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/9
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_REVERSEVOWELS_H
#define CTEST_PRO_REVERSEVOWELS_H

#include "algorithm"
#include <unordered_set>
using namespace std;
class reverseVowels_solution {
private:
  const unordered_set<char> vowels = {'a', 'e', 'i', 'o', 'u', 'A', 'E', 'I', 'O', 'U'};
  bool isVowel(char c) {
    return vowels.count(c) != 0;
  }
public:
  string reverseVowels(string s){
    int slow, fast;
    slow = 0;
    fast = s.size() - 1;
    while (slow < fast) {
      if (!isVowel(s[slow])) {
        slow++;
        continue;
      }
      if (!isVowel(s[fast])) {
        fast--;
        continue;
      }
      char tmp = s[slow];
      s[slow] = s[fast];
      s[fast] = tmp;
      slow ++;
      fast --;
    }
    return s;
  }

};


#endif//CTEST_PRO_REVERSEVOWELS_H
