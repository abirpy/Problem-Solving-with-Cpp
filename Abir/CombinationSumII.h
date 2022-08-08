//
// Created by Tahsinul Haque Abir on 7/7/22.
//

#ifndef ABIR_COMBINATIONSUMII_H
#define ABIR_COMBINATIONSUMII_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

// one digit can be more than once in the array but amra combination e repeat korte parbo na
// For example, [10,1,2,7,6,1,5] => [1, 1, 2, 5, 6, 7, 10]
// ekhane we can't have combinations like [2,2,....2] because je koy bar ase or theke beshi na
// we also can;t have [1, 2] and [2, 1] at the same time even thought we have two 1
// the way to tackle this would be to modify our tree
//      1   2 5 6 7 10 (amra ekhane second 1 skip kortesi)
//     / \
//    1   2

// Time complexity: O(n * 2^n) because the total number of possible combinations without
// repetition is actually the number of possible subsets which is 2^n and we are forming the
// outputs of max size n, so n * 2^n
// Shb shomoy tree r dike takaye time complexity ber kora jabe na for example, ei problem r tree
// ta n-ary tree cilo but the number of possible combinations without rep is 2^n which means
// there 2^n branches in our tree and thus the result
class CombinationSumII {
public:
    bool backtrack(vector<int>& candidates, int remainingSum, int idx, vector<vector<int>>&
    combinations, vector<int>& combination){
        if(remainingSum <= 0){
            if(remainingSum == 0){
                combinations.push_back(combination);
            }
            return false; // mane need to trim this branch
        }

        for(int i = idx; i < candidates.size(); ++i){
            combination.push_back(candidates[i]);
            if(!backtrack(candidates, remainingSum - candidates[i],  i + 1, combinations,
                          combination)){
                combination.pop_back(); // loop r majkhan theke ber hoye ashtesi so pop_back()
                return true; // this return is kind of like a break and a signal to the calling
                // func that don't inside this if condition cz that would create an infinite loop
            }
            // we are backtracking but only in the forward direction(including current cz
            // repetition allowed) that's why i not i + 1
            combination.pop_back(); // natural pop_back()

            //ekhn repetition handle koraa lagbe
            // I need to take i  to the last instance of repeptition if repetition exists
            int cur = candidates[i];
            while(i < candidates.size() && candidates[i] == cur){
                i++;
            }
            i--; // i is now at the last instance of repeptition if repetition exists
        }

        return true; // doesn't matter tho
    }

    vector<vector<int>> combinationSum2(vector<int>& candidates, int target) {
        // must be sorted
        sort(candidates.begin(), candidates.end());
        // the set is [1,..n] size n
        vector<vector<int>> combinations;
        vector<int> combination;

        backtrack(candidates, target, 0, combinations,combination);
        return combinations;
    }
};


#endif //ABIR_COMBINATIONSUMII_H
