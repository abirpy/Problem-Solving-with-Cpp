//
// Created by Tahsinul Haque Abir on 7/7/22.
//

#ifndef ABIR_COMBINATIONSUMIII_H
#define ABIR_COMBINATIONSUMIII_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

class CombinationSumIII {
public:
    int outputSize;
    bool backtrack(int remainingSum, int idx, vector<vector<int>>&
    combinations, vector<int>& combination){
        if(remainingSum <= 0 || combination.size() > outputSize){
            if(remainingSum == 0 && combination.size() == outputSize){
                combinations.push_back(combination);
            }
            return false;
        }

        for(int i = idx; i <= 9; ++i){
            combination.push_back(i);

            if(!backtrack(remainingSum - i,  i + 1, combinations, combination)){
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

    // Time complexity: O(k * nCk)
    vector<vector<int>> combinationSum3(int k, int n) {
        vector<vector<int>> combinations;
        vector<int> combination;
        outputSize = k;
        backtrack(n, 1, combinations,combination);
        return combinations;
    }
};


#endif //ABIR_COMBINATIONSUMIII_H
