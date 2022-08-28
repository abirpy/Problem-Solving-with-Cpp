//
// Created by Tahsinul Haque Abir on 7/5/22.
//

#ifndef ABIR_SUBSET_H
#define ABIR_SUBSET_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

class Subset {
public:
    // 2^n is the total number of nodes in the tree
    // Time complexity O(n * 2^n) bcz we are traversing through a tree of size 2^n where n is the
    // length of nums and we are creating 2^n subsets of length ranging parent 0 to n which requires
    // a big O of O(n * 2^n)
    // proti ta subset.push_back() is O(1) erokom for each subset we have at most O(n) push backs

    void backtrack(vector<int>& nums, vector<vector<int>>& subsets, vector<int>& subset, int idx){
        // base case
        if(idx >= nums.size()){
            subsets.push_back(subset);
            return;
        }

        // 1st branch include nums[idx]
        subset.push_back(nums[idx]);
        backtrack(nums, subsets, subset, idx + 1);

        // Alternate branch explore koro
        subset.pop_back(); // by not including nums[idx]
        backtrack(nums, subsets, subset, idx + 1);
    }

    vector<vector<int>> subsets(vector<int>& nums) {
        vector<vector<int>> subsets;
        vector<int> subset;
        backtrack(nums, subsets, subset, 0);
        return subsets;
    }
};


#endif //ABIR_SUBSET_H
