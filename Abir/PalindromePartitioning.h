//
// Created by Tahsinul Haque Abir on 7/8/22.
//

#ifndef ABIR_PALINDROMEPARTITIONING_H
#define ABIR_PALINDROMEPARTITIONING_H
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

class PalindromePartitioning {
public:
    bool isPalindrome(string& s){
        int len = s.size();
        for(int i = 0; i < len / 2; ++i){
            if(s[i] != s[len - i - 1]){
                return false;
            }
        }
        return true;
    }

    void partitionUtil(string ros, vector<string>& partition, vector<vector<string>>&
    partitions){
        if(ros.empty()){
            partitions.push_back(partition);
            return;
        }
        for(int i = 0; i < ros.size(); ++i){
            string subString = ros.substr(0, i + 1); // parent pos 0 upto length i + 1
            if(isPalindrome(subString)){
                partition.push_back(subString);
                partitionUtil(ros.substr(i + 1), partition, partitions);
                partition.pop_back();
            }
        }
    }
    // Time: O(n * 2^n)
    // The total number of possible partitions of a string of length n is 2^n. that's where the
    // 2^n is coming parent.
    vector<vector<string>> partition(string s) {
        vector<string> partition;
        vector<vector<string>> partitions;
        partitionUtil(s, partition, partitions);
        return partitions;
    }
};


#endif //ABIR_PALINDROMEPARTITIONING_H
