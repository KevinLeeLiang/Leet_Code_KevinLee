//
// Created by garen_lee on 2024/6/13.
/**
  ******************************************************************************
  * @file           : createTree.h
  * @author         : garen_lee
  * @brief          : None
  * @attention      : None
  * @date           : 2024/6/13
  ******************************************************************************
  */
//

#ifndef CTEST_PRO_CREATETREE_H
#define CTEST_PRO_CREATETREE_H

#include<vector>
#include<queue>
#include <iostream>
#include <string>

using namespace std;
namespace TreeNode {
    struct TreeNode {
        int val;
        TreeNode *left;
        TreeNode *right;

        TreeNode() : val(0), left(nullptr), right(nullptr) {}

        TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}

        TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
    };

    void bfs(TreeNode *root, int index, vector<int> tree_vals) {
        if (tree_vals.size() == 0)
            return;
        queue<TreeNode *> q;
        root = new TreeNode(tree_vals[0]);
        q.push(root);
        index++;
        while (!q.empty()) {
            if (index >= tree_vals.size())
                break;
            auto p = q.front();
            q.pop();
            if (tree_vals[index] != 0) {
                p->left = new TreeNode(tree_vals[index]);
                q.push(p->left);
            }
            index++;
            if (tree_vals[index] != 0) {
                p->right = new TreeNode(tree_vals[index]);
                q.push(p->right);
            }
            index++;
        }
    }

    TreeNode *createTree(vector<int> tree_vals) {
        if (tree_vals.size() == 0)
            return nullptr;
        queue<TreeNode *> q;
        TreeNode *root = new TreeNode(tree_vals[0]);
        int index = 0;
        q.push(root);
        index++;
        while (!q.empty()) {
            if (index >= tree_vals.size())
                break;
            auto p = q.front();
            q.pop();
            if (index < tree_vals.size() && tree_vals[index] != 0) {
                p->left = new TreeNode(tree_vals[index]);
                q.push(p->left);
            }
            index++;
            if (index < tree_vals.size() && tree_vals[index] != 0) {

                p->right = new TreeNode(tree_vals[index]);
                q.push(p->right);
            }
            index++;
        }
        return root;
    }

    TreeNode *createTree2(vector<int> tree_vals) {
        if (tree_vals.size() == 0)
            return nullptr;
        queue<TreeNode *> q;
        TreeNode *root = new TreeNode(tree_vals[0]);
        int index = 0;
        q.push(root);
        index++;
        while (!q.empty()) {
            if (index >= tree_vals.size())
                break;
            auto p = q.front();
            q.pop();
            if (index < tree_vals.size() && tree_vals[index] >= 0) {
                p->left = new TreeNode(tree_vals[index]);
                q.push(p->left);
            }
            index++;
            if (index < tree_vals.size() && tree_vals[index] >= 0) {

                p->right = new TreeNode(tree_vals[index]);
                q.push(p->right);
            }
            index++;
        }
        return root;
    }

    string print_tree(TreeNode *root) {
        string s;
        queue<TreeNode *> q;
        q.push(root);
        while (!q.empty()) {
            auto size = q.size();
            for (int i = 0; i < size; ++i) {
                auto node = q.front();
                q.pop();
                if (node == nullptr) {
                    s += "0";
                } else {
                    s += to_string(node->val);
                    q.push(node->left);
                    q.push(node->right);
                }
            }
        }
        for (int i = s.size() - 1; i >=0; --i) {
            if (s[i] == '0') {
                s.pop_back();
            } else {
                break;
            }
        }
        return s;
    }
};

#endif//CTEST_PRO_CREATETREE_H
