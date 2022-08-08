//
// Created by Tahsinul Haque Abir on 7/3/22.
//

#ifndef ABIR_LEVELORDERBOTTOM_H
#define ABIR_LEVELORDERBOTTOM_H

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

class levelOrderBottom {
public:
    // A traditional BFS but we will insert the values in stack when we visit them

    // This is the best structure of BFS
    vector<vector<int>> levelOrderBottomf(TreeNode* root) {
        vector<vector<int>> ans;
        if(root == nullptr){
            return ans;
        }

        queue<TreeNode*> q;

        // push it to queue
        q.push(root);

        while(!q.empty()) {
            vector<int> level; // val of each level will be separate in the level arr
            int qLen = q.size(); // qLen is equal to the numNodes in a particular level
            for(int i = 0; i < qLen; ++i){
                TreeNode* curNode = q.front();
                level.push_back(curNode->val);
                if(curNode->left != nullptr){
                    q.push(curNode->left);
                }
                if(curNode->right != nullptr){
                    q.push(curNode->right);
                }
                q.pop(); // popping curNode
            }

            if(!level.empty()){
                ans.push_back(level);
            }
        }
        reverse(ans.begin(), ans.end());

        return ans;
    }

    vector<vector<int>> zigzagLevelOrder(TreeNode* root) {
        vector<vector<int>> ans;
        if(root == nullptr){
            return ans;
        }

        queue<TreeNode*> q;

        // push it to queue
        q.push(root);

        int levelNum = 0;
        while(!q.empty()) {
            vector<int> level; // val of each level will be separate in the level arr
            int qLen = q.size(); // qLen is equal to the numNodes in a particular level
            levelNum++;
            for(int i = 0; i < qLen; ++i){
                TreeNode* curNode = q.front();
                level.push_back(curNode->val);
                if(curNode->left != nullptr){
                    q.push(curNode->left);
                }
                if(curNode->right != nullptr){
                    q.push(curNode->right);
                }
                q.pop(); // popping curNode
            }

            if(levelNum % 2 == 0){ //even level gula reversed hbe
                reverse(level.begin(), level.end());
            }

            if(!level.empty()){
                ans.push_back(level);
            }
        }
        return ans;
    }
};


#endif //ABIR_LEVELORDERBOTTOM_H
