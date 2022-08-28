//
// Created by Tahsinul Haque Abir on 7/4/22.
//

#ifndef ABIR_RIGHTSIDEVIEW_H
#define ABIR_RIGHTSIDEVIEW_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>

using namespace std;

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode() : val(0), left(nullptr), right(nullptr) {}
    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
};

// Good article // must view
// https://leetcode.com/problems/binary-tree-right-side-view/solution/

// The difference is how to find the end of the level, i.e. the rightmost element:
//
// 1. Two queues, one for the previous level and one for the current. (memory intensive)
//
// 2. One queue with sentinel to mark the end of the level.
//
// 3. One queue + level size measurement.
// Third approach is the best

class rightSideView {
    vector<int> rightSideViewF(TreeNode* root) {
        vector<int> ans;
        if(root == nullptr){
            return ans;
        }

        queue<TreeNode*> q;
        q.push(root);

        while(!q.empty()){
            int lenLevel = q.size();
            for(int i = 0; i < lenLevel; ++i){
                TreeNode *cur = q.front();
                q.pop(); // popping cur

                if(i == lenLevel - 1){
                    ans.push_back(cur->val);
                }
                // pushing cur's child
                if(cur->left != nullptr){
                    q.push(cur->left);
                }
                if(cur->right != nullptr){
                    q.push(cur->right);
                }
            }
        }
        return ans;
    }
};


#endif //ABIR_RIGHTSIDEVIEW_H
