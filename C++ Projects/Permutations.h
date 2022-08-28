//
// Created by Tahsinul Haque Abir on 7/7/22.
//

#ifndef ABIR_PERMUTATIONS_H
#define ABIR_PERMUTATIONS_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

// Time: O(n * n!) n! is the number of branches and each time we are forming a permutaion of size n
// Additional space: O(n)

// Ei type r backtracking problem is we either count the total number of branches or nodes jeta
// tight hoy and then multiple by the output size(n)

class Permutations {
    void backtrack(vector<int>& nums, vector<int>& visited, vector<vector<int>>& permutations,
                   vector<int>& permutation){
        if(permutation.size() == nums.size()){
            permutations.push_back(permutation);
            return;
        }
        for(int i = 0; i < nums.size(); ++i){
            if(visited[i] == 0){
                permutation.push_back(nums[i]);
                visited[i] = 1;

                // recursive call
                backtrack(nums, visited, permutations, permutation);

                // backtrack
                permutation.pop_back();
                visited[i] = 0;
            }
        }

    }
    vector<vector<int>> permute(vector<int>& nums) {
        vector<vector<int>> permutations;
        vector<int> permutation;
        vector<int> visited(nums.size(), 0);
        backtrack(nums, visited, permutations,permutation);
        return permutations;
    }
};


#endif //ABIR_PERMUTATIONS_H
