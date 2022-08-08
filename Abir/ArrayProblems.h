//
// Created by Tahsinul Haque Abir on 8/4/22.
//

#ifndef ABIR_ARRAYPROBLEMS_H
#define ABIR_ARRAYPROBLEMS_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;

class ArrayProblems {
    // use a hash table for the most optimzied solution
    bool containsDuplicate(vector<int>& nums) {
        unordered_set<int> hash;
        for(int &i: nums){
            if(hash.count(i) == 1){ // mane eta ekbar dhukano ase
                return true;
            }
            hash.insert(i);
        }
        return false;
    }

    int missingNumber(vector<int>& nums) {
        int n = nums.size();
        int sum = 0;
        for(int i: nums){
            sum += i;
        }
        return ((n * (n + 1)) / 2) - sum;
    }

    // all numbers are in the range [1, n] which means you can use the array indices as hash
    // n is the size so valid indices are  0 to n-1
    // nums[i] - 1 minus 1 karon valid idx upto n - 1
    vector<int> findDisappearedNumbers(vector<int>& nums) {
        vector<int> ans;
        for(int i = 0; i < nums.size(); ++i){
            if(nums[abs(nums[i]) - 1] > 0){
                nums[abs(nums[i]) - 1] = -nums[abs(nums[i]) - 1]; // negative marking
            }
        }

        for(int i = 0; i < nums.size(); ++i){
            if(nums[i] > 0){ // mane i + 1 nai thakle ith idx - hoto
                ans.push_back(i + 1);
            }
        }
        return ans;
    }

    // a ^ a = 0 and a ^ 0 = a
    int singleNumber(vector<int>& nums) {
        int res = 0;
        for(int &i: nums){
            res = res ^ i; // XOR
        }
        return res;
    }

    vector<vector<int>> construct2DArray(vector<int>& original, int m, int n) {
        if(m * n != original.size()){
            return vector<vector<int>>();
        }
        vector<vector<int>> ans(m, vector<int>(n));

        for(int i = 0; i < original.size(); ++i){
            ans[i / n][i % n] = original[i];
        }
        return ans;
    }
};


#endif //ABIR_ARRAYPROBLEMS_H
