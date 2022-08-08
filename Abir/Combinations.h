//
// Created by Tahsinul Haque Abir on 7/7/22.
//

#ifndef ABIR_COMBINATIONS_H
#define ABIR_COMBINATIONS_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

// Time: O(k * nCk!) n! is the number of branches and each time we are forming a combination of
// size n
// Additional space: O(k)

// Ei type r backtracking problem is we either count the total number of branches or nodes jeta
// tight hoy and then multiple by the output size(n)

class Combinations {
public:
    void backtrack(int n, int k, int idx, vector<vector<int>>& combinations, vector<int>&
            combination){
        if(combination.size() == k){
            combinations.push_back(combination);
            return;
        }

        for(int i = idx; i <= n; ++i){
            combination.push_back(i);
            backtrack(n, k, i + 1, combinations, combination); // we are backtracking but only
            // in the forward direction that's why idx + 1 and knon visited array lagtese na
            combination.pop_back();
        }
    }

    vector<vector<int>> combine(int n, int k) {
        // the set is [1,..n] size n
        vector<vector<int>> combinations;
        vector<int> combination;

        backtrack(n, k, 1, combinations,combination);
        return combinations;
    }
};


#endif //ABIR_COMBINATIONS_H
