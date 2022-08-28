//
// Created by Tahsinul Haque Abir on 7/5/22.
//

#ifndef ABIR_DISTANCEK_H
#define ABIR_DISTANCEK_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
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

class DistanceK {
public:
    void createGraph(unordered_map<int, vector<int>>& treeHash, TreeNode *child, TreeNode *parent){
        //visited array lagbe na cz no cycle. It is a tree
        if(parent != nullptr){
            treeHash[child->val].push_back(parent->val);
            treeHash[parent->val].push_back(child->val);
        }
        if(child->left != nullptr){
            createGraph(treeHash, child->left, child);
        }
        if(child->right != nullptr){
            createGraph(treeHash, child->right, child);
        }
    }
    vector<int> distanceK(TreeNode* root, TreeNode* target, int k) {
        // First make an adjacency list and convert the tree into graph
        // this will be O(V)
        unordered_map<int, vector<int>> treeHash;
        vector<int> ans;
        if(root == nullptr){
            return ans;
        }
        createGraph(treeHash, root, nullptr);

        // Now we will do BFS
        queue<int> q;
        unordered_set<int> visited;
        q.push(target->val);

        int levelNum = -1;
        while(levelNum != k){
            int levelLen = q.size();
            levelNum++;

            for(int i = 0; i < levelLen; ++i){
                int cur = q.front();
                q.pop();
                visited.insert(cur);
                if(levelNum == k){
                    ans.push_back(cur);
                }
                for(auto &neighbor: treeHash[cur]){
                    if(visited.count(neighbor) == 0){
                        q.push(neighbor);
                    }
                }
            }
        }
        return ans;
    }
};


#endif //ABIR_DISTANCEK_H
