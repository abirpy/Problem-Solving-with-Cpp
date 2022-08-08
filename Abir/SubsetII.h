//
// Created by Tahsinul Haque Abir on 7/6/22.
//

#ifndef ABIR_SUBSETII_H
#define ABIR_SUBSETII_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

class SubsetII {
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
        int curLast = subset.back();
        int nextElemIdx = idx + 1;
        int next = nextElemIdx < nums.size() ? nums[nextElemIdx] : -1;
        while(nextElemIdx < nums.size() && next == curLast){  // Finding the next elem and
            // excluding
            // repetition
            nextElemIdx++;
            if(nextElemIdx < nums.size()){
                next = nums[nextElemIdx];
            }
        }
        subset.pop_back(); // by not including nums[idx]
        backtrack(nums, subsets, subset, nextElemIdx);
    }

    vector<vector<int>> subsetsWithDup(vector<int>& nums) {
        sort(nums.begin(), nums.end());
        vector<vector<int>> subsets;
        vector<int> subset;
        backtrack(nums, subsets, subset, 0);
        return subsets;
    }
};


#endif //ABIR_SUBSETII_H
