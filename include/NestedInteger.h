//
// Created by garen_lee on 2024/5/21.
/**
  ******************************************************************************
  * @file           : NestedInteger.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/5/21
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_NESTEDINTEGER_H
#define CTEST_PRO_NESTEDINTEGER_H

#include <iostream>
#include <vector>
#include <string>

using namespace std;
namespace NestedInteger {
    class NestedInteger {
    public:
        // Constructor initializes an empty nested list.
        NestedInteger() : isInt(false) {}

        // Constructor initializes a single integer.
        NestedInteger(int value) : isInt(true), integer(value) {}

        // Return true if this NestedInteger holds a single integer, rather than a nested list.
        bool isInteger() const { return isInt; }

        // Return the single integer that this NestedInteger holds, if it holds a single integer
        // The result is undefined if this NestedInteger holds a nested list
        int getInteger() const { return integer; }

        // Set this NestedInteger to hold a single integer.
        void setInteger(int value) {
            isInt = true;
            integer = value;
        }

        // Set this NestedInteger to hold a nested list and adds a nested integer to it.
        void add(const NestedInteger &ni) {
            isInt = false;
            list.push_back(ni);
        }

        // Return the nested list that this NestedInteger holds, if it holds a nested list
        // The result is undefined if this NestedInteger holds a single integer
        std::vector<NestedInteger> &getList() { return list; }

    private:
        bool isInt;
        int integer;
        vector<NestedInteger> list;
    };


    class Solution {
    public:
        NestedInteger deserialize(string s) {
            int index = 0;
            return parseNestedInteger(s, index);
        }

    private:
        NestedInteger parseNestedInteger(const string &s, int &index) {
            if (s[index] == '[') {
                index++;
                NestedInteger ni;
                while (s[index] != ']') {
                    ni.add(parseNestedInteger(s, index));
                    if (s[index] == ',') {
                        index++;
                    }
                }
                index++;
                return ni;
            } else {
                return parseInteger(s, index);
            }
        }

        NestedInteger parseInteger(const string &s, int &index) {
            int sign = 1;
            if (s[index] == '-') {
                index++;
                sign = -1;
            }
            int num = 0;
            while (index < s.size() && isdigit(s[index])) {
                num = num * 10 + (s[index] - '0');
                index++;
            }
            return NestedInteger(sign * num);
        }
    };
}

#endif//CTEST_PRO_NESTEDINTEGER_H
