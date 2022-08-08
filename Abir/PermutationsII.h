//
// Created by Tahsinul Haque Abir on 7/7/22.
//

#ifndef ABIR_PERMUTATIONSII_H
#define ABIR_PERMUTATIONSII_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

class PermutationsII {
    int numSize;
    void backtrack(unordered_map<int, int>& hashmap, vector<vector<int>>& permutations,
              vector<int>& permutation){
        if(permutation.size() == numSize){
            permutations.push_back(permutation);
            return;
        }
        for(auto &numCountPair: hashmap){
            if(numCountPair.second != 0){
                permutation.push_back(numCountPair.first);
                numCountPair.second--;

                // recursive call
                backtrack(hashmap, permutations, permutation);

                // backtrack
                permutation.pop_back();
                numCountPair.second++;
            }
        }
    }

    vector<vector<int>> permuteUnique(vector<int>& nums) {
        numSize = nums.size();
        vector<vector<int>> permutations;
        vector<int> permutation;

        // we need to convert nums into a hashmap with val and count and use it as a visited
        // while backracking
        unordered_map<int, int> hashmap;
        for(auto &num: nums){
            hashmap[num]++;
        }

        backtrack(hashmap, permutations,permutation);
        return permutations;
    }
};


#endif //ABIR_PERMUTATIONSII_H
