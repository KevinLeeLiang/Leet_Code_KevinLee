//
// Created by garen_lee on 2024/5/7.
/**
  ******************************************************************************
  * @file           : NestedIterator.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/7
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_NESTEDITERATOR_H
#define CTEST_PRO_NESTEDITERATOR_H

class NestedInteger {
public:
  // Return true if this NestedInteger holds a single integer, rather than a nested list.
  bool isInteger() const;
  // Return the single integer that this NestedInteger holds, if it holds a single integer
  // The result is undefined if this NestedInteger holds a nested list
  int getInteger() const;
  // Return the nested list that this NestedInteger holds, if it holds a nested list
  // The result is undefined if this NestedInteger holds a single integer
  const vector<NestedInteger> &getList() const;
};

class NestedIterator {
public:
  NestedIterator(vector<NestedInteger> &nestedList) {

  }

  int next() {

  }

  bool hasNext() {

  }
};


#endif//CTEST_PRO_NESTEDITERATOR_H
