//
// Created by Tahsinul Haque Abir on 7/16/22.
//

#ifndef ABIR_DYNAMICPROGRAMMING_H
#define ABIR_DYNAMICPROGRAMMING_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cmath>

using namespace std;

class DynamicProgramming {
public:
    // maximum profit you can achieve by picking n items with a bag size of weightCapacity
    // our relevant weight and profit array [1....n]
    int knapsack01(int weight[], int profit[], int weightCapacity, int n, vector<vector<int>>& mem){
        if(mem[n][weightCapacity] != -1){
            return mem[n][weightCapacity];
        }
        if(n == 0 || weightCapacity == 0){
            return mem[n][weightCapacity] = 0;
        }
        if(weight[n] > weightCapacity){
            return mem[n][weightCapacity] = knapsack01(weight, profit, weightCapacity, n - 1); // not picking up the last
            // item
        }
        return mem[n][weightCapacity] = max(profit[n] + knapsack01(weight, profit, weightCapacity
        - weight[n], n - 1, mem), // took nth item
                   knapsack01(weight, profit, weightCapacity, n - 1, mem)); // not taking it
    }

    int knapsack01(int weight[], int profit[], int weightCapacity, int n){
        int dp[n + 1][weightCapacity + 1];

        for(int i = 0; i <= n; ++i){
            for(int j = 0;  j <= weightCapacity; ++j){
                if(i == 0 || j == 0){
                    dp[i][j] = 0;
                }
                else if(weight[i] > j){ // processing the ith item // j is the current
                    // weightCapacity
                    // if weight array is [1...n] then weight[i] otherwise i - 1
                    dp[i][j] = dp[i - 1][j]; // not including current item weight capacity
                    // remains same
                }
                else{
                    dp[i][j] = max(profit[i] + dp[i - 1][j - weight[i]], dp[i - 1][j]);
                }
            }
        }
        return dp[n][weightCapacity];
    }

    // if the set contains 0, then specially treat kora lagbe. You have to avoid taking 0 at the
    // beginning and then finally add 2^(number of zeroes in the set) with final ans
    int countNumSubsetSum(int set[], int target, int n, unordered_map<int, unordered_map<int, int>>
    mem){
        // target can be negative in some base cases that's why can't take array
        if(mem[target][n] != -1){
            return mem[target][n];
        }

        if(target != 0 && n == 0){ // no
            return 0;
        }
        else if(target == 0 && n >= 0){
            return 1;
        }
        else if(set[n - 1] > target){ // indexing r jonne n - 1
            // we cannot take it
            return mem[target][n] = countNumSubsetSum(set, target, n - 1, mem);
        }
        else{
            return mem[target][n] = countNumSubsetSum(set, target - set[n - 1], n - 1, mem) +
                    countNumSubsetSum(set, target, n - 1, mem);
        }
    }

    string make_pair(int i, int j){
        return to_string(i) + "*" + to_string(j);
    }

    // if the set contains 0, then specially treat kora lagbe. You have to avoid taking 0 at the
    // beginning and then finally add 2^(number of zeroes in the set) with final ans

    int countNumSubsetSumTab(vector<int>& set, int target, int n){
        int mem[target + 1][n + 1];

        int count0 = 0, sum = 0;
        for(int j = 0; j < n; ++j){
            if(set[j] == 0){
                count0++;
            }
            sum += set[j];
        }

        if(target > sum){
            return 0;
        }

        for(int i = 0; i <= target; ++i){
            for(int j = 0; j <= n; ++j){
                if(i == 0 && j >= 0){
                    mem[i][j] = 1;
                }
                else if(i != 0 && j == 0){
                    mem[i][j] = 0;
                }
                else if(set[j - 1] > i || set[j - 1] == 0){ // i represents cur target and j
                    // represents cur n
                    mem[i][j] = mem[i][j - 1];
                }
                else{
                    mem[i][j] = mem[i - set[j - 1]][j - 1] + mem[i][j - 1]; // take it and don't
                    // take it
                }
            }
        }

        return mem[target][n] * pow(2, count0);
    }

    // min number of coins required to form target where each coin can be taken infinite number
    // of times
    // if not possible return -1
    // unbounded knapsack
    int coinChangeI(vector<int>& set, int target, int n, vector<vector<int>>& mem){
        if(mem[n][target] != -1){
            return mem[n][target] >= INT16_MAX ? -1 : mem[n][target];
        }
        if(target == 0){
            return 0;
        }
        else if(n == 0){ // && target != 0
            return INT16_MAX;
        }
        else if(set[n - 1] > target){
            return mem[n][target] = coinChangeI(set, target, n - 1, mem);
        }
        else{
            return mem[n][target] = min((1 + coinChangeI(set, target - set[n - 1], n, mem)), //
                                        // took it but can take it more so n unchanged
                    coinChangeI(set, target, n - 1, mem)); // didn't take it
        }
    }

    // unbounded knapsack
    // 1 + mem[i][j - set[i - 1]] we are counting the number of coins it take to form target and
    // we want to minimize that number
    // so we have to count the number of things we include in our bag that's why 1 + mem[i][j -
    // set[i - 1]]
    int coinChangeITab(vector<int>& set, int target, int n){
        int mem[n + 1][target + 1];

        for(int i = 0; i <= n; ++i){
            for(int j = 0; j <= target; ++j){
                if(j == 0){
                    mem[i][j] = 0;
                }
                else if(i == 0){ // j!= 0
                    mem[i][j] = INT16_MAX; // this is a technique to make the math work when we
                    // take min
                }
                else if(set[i - 1] > j){
                    mem[i][j] = mem[i - 1][j];
                }
                else{
                    mem[i][j] = min(1 + mem[i][j - set[i - 1]], mem[i - 1][j]);
                }
            }
        }

        return mem[n][target] >= INT16_MAX ? -1 : mem[n][target];
    }

    // count the total number of ways we can form target where each coin can be taken infinite
    // number of times
    int coinChangeII(vector<int>& set, int target, int n, vector<vector<int>>& mem){
        if(mem[n][target] != -1){
            return mem[n][target];
        }
        if(target == 0){
            return 1; // found one way to form target
        }
        else if(n == 0){ // target != 0
            return 0;
        }
        else if(set[n - 1] > target){ // if set contains 0 then 0 skip kora lagbe and last 2^
            // (number of zero) ad kora lagbe
            return mem[n][target] = coinChangeII(set, target, n - 1, mem); // not including it
        }
        else{
            return mem[n][target] = coinChangeII(set, target, n - 1, mem) + // not including it
                                    coinChangeII(set, target - set[n - 1], n, mem); // including
                                    // it but can take n more times
        }
    }

    int coinChangeIITab(vector<int>& set, int target, int n){
        int mem[n + 1][target + 1];

        for(int i = 0; i <= n; ++i){
            for(int j = 0; j <= target; ++j){
                if(j == 0){
                    mem[i][j] = 1;
                }
                else if(i == 0){ // j!= 0
                    mem[i][j] = 0;
                }
                else if(set[i - 1] > j){
                    mem[i][j] = mem[i - 1][j];
                }
                else{
                    mem[i][j] = mem[i][j - set[i - 1]] + mem[i - 1][j];
                }
            }
        }

        return mem[n][target];
    }

    // evabe rod cato so thaat you find max profit
    // profit by selling a rod of len n is price[n - 1]
    // nahole len array die dite pare
    // unbounded knapsack
    int rodCutting(int price[], int totalLen, int n, vector<vector<int>>& mem){
        if(mem[n][totalLen] != -1){
            return mem[n][totalLen];
        }

        if(totalLen == 0 || n == 0){
            return 0;
        }
        else if(n > totalLen){
            return mem[n][totalLen] = rodCutting(price, totalLen, n - 1, mem);
        }
        else {
            return mem[n][totalLen] = max(price[n - 1] + rodCutting(price, totalLen - n, n, mem),
                       rodCutting(price, totalLen, n - 1, mem));
        }
    }

    int rodCuttingTab(vector<int>& price, int totalLen, int n){
        vector<vector<int>> mem(n + 1, vector<int>(totalLen + 1));

        for(int i = 0; i <= n; ++i){
            for(int j = 0; j <= totalLen; ++j){
                if(i == 0 || j == 0){
                    mem[i][j] = 0;
                }
                else if(i > j){
                    mem[i][j] = mem[i - 1][j];
                }
                else{
                    mem[i][j] = max(price[i - 1] + mem[j - i][i], mem[i - 1][j]);
                }
            }
        }
        return mem[n][totalLen];
    }

    // you can either remove an elem parent left side of array or right side
    // find the min number of steps to reduce X to 0
    // in minimization problems you can assign infinity to a branch which makes the situation
    // impossible
    // giving an infinite weight is essentially equal to eliminating the branch
    string make_pair_3(int l, int r, int target){
        return to_string(l) + "*" + to_string(r) + "*" + to_string(target);
    }

    // This problem can be optimally solved by converting into subArray sum problem
    int minOpToReduceXto0(vector<int>& set, int l, int r, int target, unordered_map<string, int>&
            mem){
        string key = make_pair_3(l, r, target);
        if(mem.count(key) == 1){
            return mem[key] >= INT16_MAX ? -1 : mem[key];
        }

        if(target == 0){
            return 0;
        }
        else if(l > r){
            return mem[key] = INT16_MAX; // infinite weight means not possible
        }
        else{
            return mem[key] = min(
                    1 + minOpToReduceXto0(set, l + 1, r, target - set[l], mem),
                    1 + minOpToReduceXto0(set, l, r - 1, target - set[r], mem)
                    );
        }
    }

    // p1 and p2 will point to 0 of s1 and s2 initalliy
    bool interlearvingString(string& s1, string& s2, string& s3, string& resultStr, int p1, int
    p2, unordered_map<string, bool>& dp){
        string key = make_pair(p1, p2);
        if(dp.count(key) != 0){
            return dp[key];
        }

        if(p1 == s1.length() && p2 == s2.length()){
            return resultStr == s3;
        }
        else{
            if(p1 < s1.length()){
                resultStr.push_back(s1[p1]);
                if(interlearvingString(s1, s2, s3, resultStr, p1 + 1, p2, dp)){
                    return true;
                }
                else{
                    dp[key] = false;
                }
                resultStr.pop_back();
            }
            if(p2 < s2.length()){
                resultStr.push_back(s2[p1]);
                if(interlearvingString(s1, s2, s3, resultStr, p1, p2 + 1, dp)){
                    return true;
                }
                else{
                    dp[key] = false;
                }
                resultStr.pop_back();
            }
        }
    }

    int LCS(string &s1, string &s2, string& lcs){
        // create a dp table
        int len1 = s1.length();
        int len2 = s2.length();
        int dp[len1 + 1][len2 + 1];

        for(int i = 0; i <= len1; ++i){
            int flag = true;
            for(int j = 0; j <= len2; ++j){
                if(i == 0 || j == 0){
                    dp[i][j] = 0;
                }
                else if(s1[i - 1] == s2[j - 1]){ // matching char
                    if(flag){
                        lcs.push_back(s1[i - 1]);
                        flag = false;
                    }
                    dp[i][j] = 1 + dp[i - 1][j - 1];
                }
                else{ // not matching char
                    dp[i][j] = max(dp[i - 1][j], dp[i][j - 1]);
                }
            }
        }
        return dp[len1][len2];
    }

    bool wordBreakUtil(string& s, unordered_set<string>& dictSet, unordered_map<int, bool>& map,
                       int i) {
        if(map.count(i) == 1){
            return map[i];
        }
        if(i == s.length()){
            return true;
        }
        for(int j = i + 1; j <= s.length(); ++j){
            if(dictSet.count(s.substr(i, j - i)) && wordBreakUtil(s, dictSet, map, j)){
                return map[i] = true;
            }
        }
        return map[i] = false;
    }

    // recursion + memoization
    bool wordBreak(string s, vector<string>& wordDict) {
        unordered_set<string> dictSet;
        for(string &s : wordDict){
            dictSet.insert(s);
        }
        unordered_map<int, bool> map;
        return wordBreakUtil(s, dictSet, map, 0);
    }

    // unbounded knapsack
    int climbStairsUtil(int n, vector<int>& dp) {
        if(n == 0){
            return 1;
        }
        else if(n < 0){
            return 0;
        }

        if(dp[n] != -1){
            return dp[n];
        }
        else{
            return dp[n] = climbStairsUtil(n - 1, dp) + climbStairsUtil(n - 2, dp);
        }
    }

    // unbounded knapsack
    int climbStairs(int n) {
        vector<int> dp(n + 1, -1);
        return climbStairsUtil(n, dp);
    }

    int climbStairsTab(int n){
        vector<int> dp(n + 1);

        for(int i = 0; i <= n; ++i){
            if(i == 0 || i == 1){
                dp[i] = 1;
            }
            else{
                dp[i] = dp[i - 1] + dp[i - 2];
            }
        }
        return dp[n];
    }

    // Kadanes
    int maxSubArray(vector<int>& nums) {
        int muh = nums[0], msf = nums[0]; // max upto here and max so far

        for(int i = 1; i < nums.size(); ++i){
            muh = max(nums[i], muh + nums[i]);
            msf = max(msf, muh);
        }
        return msf;
    }

    class NumArray {
    private:
        vector<int> cumSum;
        int len;
    public:
        NumArray(vector<int>& nums) {
            len = nums.size();
            int sum = 0;
            for(int i = 0; i < len; ++i){
                sum += nums[i];
                cumSum.push_back(sum);
            }
        }

        int sumRange(int left, int right) {
            return left == 0 ? cumSum[right] : cumSum[right] - cumSum[left - 1];
        }
    };

    vector<int> countBits(int n) {
        vector<int> dp(n + 1);
        dp[0] = 0;

        // dividing a number by 2 means right shifting it by 1. if the number is even the no of
        // set bits remain same. Otherwise it increases by 1
        for(int i = 1; i <= n; ++i){
            dp[i] = dp[i / 2] + i%2;
        }
        return dp;
    }

    // sum = sum of all elem
    // s1 = sum of + and s2 = sum of -
    // sum = s1 + s2
    // target = s1 - s2 = s1 - (sum - s1) = 2s1 - sum
    // => s1 = (sum + target) / 2;
    // so this problem is about the number of ways we can form of sum value (sum + target) / 2
    int findTargetSumWays(vector<int>& nums, int target) {
        int sum = 0;
        int n =  nums.size();
        for(int i = 0; i < n; ++i){
            sum += nums[i];
        }

        int newTarget;
        if(target > sum || (sum + target) % 2 != 0){
            return 0;
        }
        else{
            newTarget = (sum + target) / 2;
        }

        unordered_map<int, unordered_map<int, int>> dp;
        int count0 = 0;

        for(int &i: nums){
            if(i == 0){
                count0++;
            }
        }

        for(int i = 0; i <= n; ++i){
            bool flag = true;
            for(int j = 0; j <= newTarget; ++j){
                if(j == 0){
                    dp[i][j] = 1;
                }
                else if(i == 0){ // j != 0
                    dp[i][j] = 0;
                }
                else if(nums[i - 1] > j || nums[i - 1] == 0){
                    dp[i][j] = dp[i - 1][j];
                }
                else{
                    dp[i][j] = dp[i - 1][j - nums[i - 1]] + dp[i - 1][j];
                    // including and not including it
                }
            }
        }

        return dp[n][newTarget] * pow(2, count0);
    }
    int robUtil(vector<int>& nums, vector<vector<int>>& mem, int s, int e){
        if(s == e){
            return nums[s];
        }
        else if(s > e){
            return 0;
        }

        if(mem[s][e] != -1){
            return mem[s][e];
        }
        else{
            return mem[s][e] = max(nums[s] + robUtil(nums, mem, s + 2, e),
                                    robUtil(nums, mem, s + 1, e));
        }
    }


    int rob(vector<int>& nums){
        vector<vector<int>> mem(nums.size(), vector<int>(nums.size(), -1));
        return robUtil(nums, mem, 0, nums.size() - 1);
    }

    int robTab(vector<int>& nums){
        int len = nums.size();
        vector<vector<int>> dp(len, vector<int>(len));

        for(int j = 0; j < len; ++j){ // i -> start
            for(int i = 0; i < len; ++i){ // j -> end
                if(i == j - 1){
                    dp[i][j] = max(nums[j], dp[i][j - 1]);
                }
                else if(i == j){
                    dp[i][j] = nums[j];
                }
                else if(i > j){
                    dp[i][j] = 0;
                }
                else{
                    dp[i][j] = max(nums[j] + dp[i][j - 2], dp[i][j - 1]);
                }
            }
        }

        return dp[0][len - 1];
    }

    int robTab2(vector<int>& nums){
        int len = nums.size();
        vector<vector<int>> dp(len, vector<int>(len));

        for(int j = 0; j < len; ++j){ // i -> start
            for(int i = 0; i < len; ++i){ // j -> end
                if(i == j - 1){
                    dp[i][j] = max(nums[j], dp[i][j - 1]);
                }
                else if(i == j){
                    dp[i][j] = nums[j];
                }
                else if(i > j){
                    dp[i][j] = 0;
                }
                else if(j ==  len - 1 && i == 0){
                    dp[i][j] = max(nums[j] + dp[i + 1][j - 2], dp[i][j - 1]);
                    // when includes last don't include first
                }
                else{
                    dp[i][j] = max(nums[j] + dp[i][j - 2], dp[i][j - 1]);
                }
            }
        }

        return dp[0][len - 1];
    }

    int coinChange(vector<int>& coins, int amount) {
        int n = coins.size();
        vector<vector<int>> dp(n + 1, vector<int>(amount + 1));

        for(int i = 0; i <= n; ++i){
            for(int j = 0; j <= amount; ++j){
                if(j == 0){
                    dp[i][j] = 0;
                }
                else if(i == 0){
                    dp[i][j] = INT16_MAX;
                }
                else if(coins[i - 1] > j){
                    dp[i][j] = dp[i - 1][j];
                }
                else{
                    // unbounded knapsack
                    dp[i][j] = min(1 + dp[i][j - coins[i - 1]], dp[i - 1][j]);
                }
            }
        }

        return dp[n][amount] >= INT16_MAX ? -1 : dp[n][amount];
    }

    int maxProduct(vector<int>& nums) {
        int globalMax = nums[0];
        int maxIncludingIdx = nums[0];
        int minIncludingIdx = nums[0];

        bool flag = false;
        for(int i = 1; i < nums.size(); ++i){
            if(nums[i] == 0){
                flag = true;
                maxIncludingIdx = 1;
                minIncludingIdx = 1;
                continue;
            }
            int tmp = maxIncludingIdx;
            maxIncludingIdx = max(maxIncludingIdx * nums[i], max(minIncludingIdx * nums[i],nums[i]));
            minIncludingIdx = min(tmp * nums[i], min(minIncludingIdx * nums[i], nums[i]));
            globalMax = max(globalMax, maxIncludingIdx);
        }

        return flag ? max(globalMax, 0) : globalMax;
    }

    int lengthOfLIS(vector<int>& nums) {
        int maxL = 1;

        vector<int> LIS(nums.size(), 1);

        for(int i = 1; i < nums.size(); ++i){
            for(int j = 0; j < i; ++j){
                if(nums[i] > nums[j] && LIS[j] + 1 > LIS[i]){
                    LIS[i] = LIS[j] + 1;
                }
                maxL = max(maxL, LIS[i]);
            }
        }
        return maxL;
    }

    // longest palindromic substring
    string longestPalindrome(string s) {
        pair<int, int> longestRange(0, 0);
        int len = s.size();
        bool dp[len][len];

        // Age first col fillup kora lagbe then second col
        // substring problem gula generally col by col fillup kora lage
        // for a proper substring i <= j
        // i -> row j -> col
        for(int j = 0; j < len; ++j){
            for(int i = 0; i < len; ++i){
                if(i >= j){ // empty subs and one char sub are palindrome
                    dp[i][j] = true;
                }
                else{
                    dp[i][j] = s[i] == s[j] && dp[i + 1][j - 1];

                    if(dp[i][j] && j - i > longestRange.second - longestRange.first){
                        longestRange.first = i;
                        longestRange.second = j;
                    }
                }
            }
        }

        return s.substr(longestRange.first, longestRange.second - longestRange.first + 1);
    }

    // O(n) solution cz mem is O(n)
    int numDecodingsUtil(string &s, unordered_map<string, char>& encodings, int idx,
                          unordered_map<int, int>& mem){
        if(idx == s.size()){
            return 1;
        }
        if(mem.count(idx) == 1){
            return mem[idx];
        }

        int count = 0;
        for(int i = idx + 1; i <= s.size(); ++i){
            if(encodings.count(s.substr(idx, i - idx))){
                count += numDecodingsUtil(s, encodings, i, mem);
            }
            else{
                break;
            }
        }
        return mem[idx] = count;
    }


    int numDecodings(string s) {
        unordered_map<string, char> encodings;
        unordered_map<int, int> mem;

        encodings["1"] = 'A';
        for(int i = 1; i <= 25; ++i){
            encodings[to_string(i + 1)] = 'A' + i;
        }

        return numDecodingsUtil(s, encodings, 0, mem);
    }

    int uniquePaths(int m, int n) {
        int dp[m + 1][n + 1];

        for(int i = 0; i <= m; ++i){
            for(int j = 0; j <= n; ++j){
                if(i == 0 || j == 0){
                    dp[i][j] = 0;
                }
                else if(i == 1 && j == 1){
                    dp[i][j] = 1;
                }
                else{
                    dp[i][j] = dp[i][j - 1] + dp[i - 1][j];
                }
            }
        }

        return dp[m][n];
    }

    bool canJumpUtil(vector<int>& nums, int idx, unordered_map<int, bool>& mem){
        if(idx == nums.size() - 1){
            return true;
        }
        else if(nums[idx] == 0){
            return false;
        }

        if(mem.count(idx) == 1){
            return mem[idx];
        }
        else{
            for(int i = 1; i <= nums[idx]; ++i){
                if(canJumpUtil(nums, idx + i, mem)){
                    return mem[idx] = true; // ekta branch true pailei return
                }
            }
            return mem[idx] = false; // kno branch true na hoile then false
        }

    }

    bool canJump(vector<int>& nums) {
        unordered_map<int, bool> mem;
        return canJumpUtil(nums, 0, mem);
    }

    // Greedy sol
    bool canJumpG(vector<int>& nums) {
        int post = nums.size() - 1; // starting parent the end and bringing the post to the front

        for(int i = nums.size() - 2; i >= 0; --i){
            if(i + nums[i] > post){
                post = i;
            }
        }

        // post can be now either at 0 or some other values. 0 hole possible
        return post == 0;
    }

    // Greedy sol
    // you are assured that it is possible to reach the end
    bool numOfMinJumpG(vector<int>& nums) {
        // a greedy bfs solu with a window
        int res = 0;
        int l = 0, r = 0;

        while(r < nums.size()){
            int farthest = 0;
            for(int i = l; i <= r; ++i){
                farthest = max(farthest, i + nums[i]); // being greedy and calculating farthest
            }
            // moving the window
            l = r + 1;
            r = farthest;
            res++;
        }
        return res;
    }

    // just all col % 2 kore dile memory optimized hoy onek
    int countSubstrings(string s) {
        int count = 0;
        int len = s.size();
        bool dp[len][2];

        for(int j = 0; j < len; ++j){ // col
            for(int i = 0; i < len; ++i){ // row
                if(i >= j){
                    dp[i][j % 2] = true;
                    if(i == j){
                        count++;
                    }
                }
                else{
                    dp[i][j % 2] = s[i] == s[j] && dp[i + 1][(j - 1) %2];
                    if(dp[i][j % 2]){
                        count++;
                    }
                }
            }
        }
        return count;
    }

    // Variation of the LIS problem
    int findNumberOfLIS(vector<int>& nums) {
        int numLIS = 1;
        int maxLen = 1;

        vector<int> LIS(nums.size(), 1);
        vector<int> LISLen(nums.size(), 1);

        for(int i = 1; i < nums.size(); ++i){
            for(int j = 0; j < i; ++j){
                if(nums[i] > nums[j]){
                    if(LIS[j] + 1 > LIS[i]){
                        LIS[i] = LIS[j] + 1;
                        LISLen[i] = LISLen[j];
                    }
                    else if(LIS[j] + 1 == LIS[i]){ // mane already ekta ase
                        LISLen[i] += LISLen[j];
                    }
                }
            }
            maxLen = max(maxLen, LIS[i]);
        }

        int result = 0;
        // there can be more than one maxLen like for the array [2,2,2,2,2] the LIS array will be
        // [1,1,1,1,1] so you canhave to add all LISLen[i] which are at maxLen idxes
        for(int i = 0; i < LIS.size(); ++i){
            if(LIS[i] == maxLen){
                result += LISLen[i];
            }
        }
        return result;
    }

    string canPartitionKSubsetsMakeKey(int i, int j, int l){
        return to_string(i) + "*" + to_string(j) + "*" + to_string(l);
    }

    // O(k * 2^n) solu hard backtracing problem watch Neetcode
    bool canPartitionKSubsetsUtil(vector<int>& nums, int idx, int target, int& initalTarget, int k,
                                  unordered_set<int>& visited, unordered_map<string, bool>& mem){
        string key = canPartitionKSubsetsMakeKey(idx, target, k);
        if(k == 0){
            return true;
        }
        else if(idx < 0 && target != 0){ // && target != 0
            return false;
        }
        else if(target == 0){
            return mem[key] = canPartitionKSubsetsUtil(nums, nums.size() - 1, initalTarget,
                                                 initalTarget, k - 1, visited, mem);
        }
        else if(mem.count(key) == 1){
            return mem[key];
        }
        else if(nums[idx] > target){
            return mem[key] = canPartitionKSubsetsUtil(nums, idx - 1, target, initalTarget, k,
                                                     visited, mem); // not
            // including it
        }
        else{
            if(visited.count(idx) == 0){
                visited.insert(idx);
                // if nums[idx] hasn't been taken yet in any subset which summed up target
                if(canPartitionKSubsetsUtil(nums, idx - 1, target - nums[idx], initalTarget, k,
                                            visited, mem)){ //
                    // including it
                    return mem[key] = true;
                }
                visited.erase(idx);
            }
            // has been visited before so let's try by not including it
            return mem[key] = canPartitionKSubsetsUtil(nums, idx - 1, target, initalTarget, k,
                                                       visited, mem); // not including it
        }
    }

    // O(k * 2^n) solu hard backtracing problem watch Neetcode
    bool canPartitionKSubsets(vector<int>& nums, int k) {
        unordered_set<int> visited;
        unordered_map<string, bool> mem;
        int sum = 0;
        for(int i = 0; i < nums.size(); ++i){
            sum += nums[i];
        }

        int target;
        if(sum % k == 0){
            target = sum / k;
            return canPartitionKSubsetsUtil(nums, nums.size() - 1, target, target, k,
                                            visited, mem);
        }
        else{
            return false;
        }
    }

    string maxProfitKey(int i, int j){
        return to_string(i) + "*" + to_string(j);
    }
    // O(n) solution cz the number of possible keys is n * 2 = 2n which is O(n)
    int maxProductutil(vector<int>& prices, unordered_map<string, int>& mem, int idx, int canBuy){
        if(idx >= prices.size()){
            return 0;
        }
        string key = maxProfitKey(idx, canBuy);
        if(mem.count(key) == 1){
            return mem[key];
        }
        if(canBuy == 1){ // in the state of buying
            int buy = maxProductutil(prices, mem, idx + 1, 0) - prices[idx]; // state changed
            int coolDown = maxProductutil(prices, mem, idx + 1, 1); // state unchanged
            return mem[key] = max(buy, coolDown);
        }
        else{ // canBuy == 0
            int sell = maxProductutil(prices, mem, idx + 2, 1) + prices[idx]; // state changed
            int coolDown = maxProductutil(prices, mem, idx + 1, 0); // state unchanged
            return mem[key] = max(sell, coolDown);
        }
    }


    // interesting state DP problem
    // there can be only two state for this problem: buying state or not buying => selling state
    // and you can always have cooldowns
    // after selling once you must have one cooldown
    int maxProfit(vector<int>& prices) {
        unordered_map<string, int> mem;
        return maxProductutil(prices, mem, 0, 1);
    }

    // No cooldown constraint
    // you can only buy and sell once
    // buy low sell high
    // Sliding window problem
    int maxProfit2(vector<int>& prices) {
        int l = 0, r = 1;
        int maxProfit = 0;
        while(r < prices.size()){
            maxProfit = max(maxProfit, prices[r] - prices[l]);
            if(prices[r] < prices[l]){
                l = r;
            }
            r++;
        }

        return maxProfit; // peak - valley
    }
};




#endif //ABIR_DYNAMICPROGRAMMING_H
