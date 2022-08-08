//
// Created by Tahsinul Haque Abir on 7/7/22.
//

#ifndef ABIR_COMBINATIONSUM_H
#define ABIR_COMBINATIONSUM_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

// Complexity is very nicely explained in the article. It's a n-ary tree
// https://leetcode.com/problems/combination-sum/solution/
// Our basic strategy of avoiding repeting combination will be same. We'll continue from an
// element (including it) forward only
// if the minimal elem is m and t is the target sum the height of the tree is gonna be T/m.
// therefore time complexity is O(n^T/m)

// We sort the array and do pruning to make the dfs efficient
class CombinationSum {
public:
    bool backtrack(vector<int>& candidates, int remainingSum, int idx, vector<vector<int>>&
    combinations, vector<int>& combination){
        if(remainingSum <= 0){
            if(remainingSum == 0){
                combinations.push_back(combination);
            }
            return false;
        }

        for(int i = idx; i < candidates.size(); ++i){
            combination.push_back(candidates[i]);

            if(!backtrack(candidates, remainingSum - candidates[i],  i, combinations, combination)){
                combination.pop_back(); // loop r majkhan theke ber hoye ashtesi so pop_back()
                return true; // this return is kind of like a break and a signal to the calling
                // func that don't inside this if condition cz that would create an infinite loop
            }
            // we are backtracking but only in the forward direction(including current cz
            // repetition allowed) that's why i not i + 1
            combination.pop_back(); // natural pop_back()
        }

        return true; // doesn't matter tho
    }

    vector<vector<int>> combinationSum(vector<int>& candidates, int target) {
        // must be sorted
        sort(candidates.begin(), candidates.end());
        // the set is [1,..n] size n
        vector<vector<int>> combinations;
        vector<int> combination;

        backtrack(candidates, target, 0, combinations,combination);
        return combinations;
    }
};


#endif //ABIR_COMBINATIONSUM_H
