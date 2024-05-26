//
// Created by Tahsinul Haque Abir on 8/12/22.
//

#ifndef BLIND75_BLIND_H
#define BLIND75_BLIND_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode() : val(0), left(nullptr), right(nullptr) {}
    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
};

const vector<vector<int>> DIR({
      {1, 0},
      {-1, 0},
      {0, 1},
      {0, -1}
      });

class Blind {
public:
    int maxProfit(vector<int>& prices){
        int l = 0, h = 1;
        int maxP = 0;

        while(h < prices.size()){
            int profit = prices[h] - prices[l];
            maxP = max(maxP, profit);

            if(prices[h] < prices[l]){
                l = h;
            }
            h++;
        }
        return maxP;
    }

    // unsorted so hashtable most optimized
    bool containsDuplicate(vector<int>& nums){
        unordered_set<int> hash;

        for(int &i: nums){
            if(hash.count(i) == 1){
                return true;
            }
            hash.insert(i);
        }
        return false;
    }

    vector<int> productExceptSelf(vector<int>& nums) {
        vector<int> ans(nums.size());
        int prefix = 1;

        for(int i = 0; i < nums.size(); ++i){
            ans[i] = prefix;
            prefix *= nums[i];
        }

        int postfix = 1;
        for(int i = nums.size() - 1; i >= 0; --i){
            ans[i] *= postfix;
            postfix *= nums[i];
        }

        return ans;
    }

    int maxProduct(vector<int>& nums){
        if(nums.size() == 1){
            return nums[0];
        }
        int globalMax, lmax, lmin;
        globalMax = nums[0];
        lmin = lmax = nums[0] == 0 ? 1 : nums[0];

        for(int i = 1; i < nums.size(); ++i){
            if(nums[i] == 0){
                lmin = lmax = 1;
            }
            int tmp = lmin;
            lmin = min(lmax * nums[i], min(tmp * nums[i], nums[i]));
            lmax = max(lmax * nums[i], max(tmp * nums[i], nums[i]));
            globalMax = max(globalMax, lmax);
        }
        return globalMax;
    }

    int findMin(vector<int>& nums) {
        int l = 0, r = nums.size();

        if(l <= r){
            return nums[l];
        }

        while(l <= r){
            int mid = (l + r) / 2;

            if(mid - 1 >= 0 && nums[mid] < nums[mid - 1]){
                return nums[mid];
            }
            else if(mid + 1 < nums.size() && nums[mid] > nums[mid + 1]){
                return nums[mid + 1];
            }
            else if(nums[mid] > nums[l]){
                l = mid + 1;
            }
            else{
                r = mid - 1;
            }
        }
        return -1;
    }

    int search(vector<int>& nums, int target) {
        int l = 0, r = nums.size() - 1;

        while(l <= r){
            int mid = (l + r) / 2;

            if(nums[mid] == target){
                return mid;
            }
            else if(nums[mid] >= nums[l]){
                if(target >= nums[l] && target < nums[mid]){
                    r = mid - 1;
                }
                else{
                    l = mid + 1;
                }
            }
            else{
                if(target > nums[mid] && target <= nums[r]){
                    l = mid + 1;
                }
                else{
                    r = mid - 1;
                }
            }
        }

        return -1;
    }

    int findFirst(vector<int>& nums, int target){
        int l = 0, r = nums.size();

        while(l < r){
            int mid = (l + r) / 2;
            if(target <= nums[mid]){
                r = mid;
            }
            else{
                l = mid + 1;
            }
        }
        // l == r ekhn correct position e
        return l < nums.size() && nums[l] == target ? l : -1;
    }

    int findLast(vector<int>& nums, int target){
        int l = 0, r = nums.size() - 1;

        while(l <= r){
            int mid = (l + r) / 2;

            if(nums[mid] == target && (mid + 1 == nums.size() || nums[mid] != nums[mid + 1])){
                return mid;
            }
            else if(nums[mid] > target){
                r = mid - 1;
            }
            else{
                l = mid + 1;
            }
        }
        return -1;
    }


    vector<int> searchRange(vector<int>& nums, int target) {
        vector<int> ans(2, -1);
        if(nums.size() == 0){
            return ans;
        }

        ans[0] = findFirst(nums, target);
        ans[1] = findLast(nums, target);
        cout << ans[0] << ans[1];
        return ans;
    }


    // [1....n]
//    int firstBadVersion(int n) {
//        int l = 1, r = n;
//
//        while (l < r){
//            int mid = (l + r) / 2;
//            if(isBadVersion(mid)){
//                r = mid;
//            }
//            else{
//                l = mid + 1;
//            }
//        }
//        return l;
//    }

    int maxArea(vector<int>& height) {
        int l = 0, r = height.size() - 1;
        int maxArea = 0;

        while(l < r){
            int minHeight = min(height[l], height[r]);
            int area = (r - l) * minHeight;
            maxArea = max(maxArea, area);

            if(height[l] < height[r]){
                l++;
            }
            else{
                r--;
            }
        }
        return maxArea;
    }

    vector<vector<int>> threeSum(vector<int>& nums) {
        sort(nums.begin(), nums.end());
        vector<vector<int>> ans;

        for(int i = 0; i < nums.size(); ++i){
            if(i != 0 && nums[i] == nums[i - 1]){ // can't take a value as first which has
                // already been taken before
                continue;
            }

            int first = nums[i]; // first
            // now you have to find to find l and r such that l + r = -f
            int l = i + 1, r = nums.size() - 1;

            while(l < r){
                int threeSum = first + nums[l] + nums[r];
                if(threeSum < 0){
                    l++;
                }
                else if(threeSum > 0){
                    r--;
                }
                else{
                    ans.push_back(vector<int>({first, nums[l], nums[r]}));
                    l++;
                    // ekta  l++ tow amra noramlly kori
                    // the next loop is for time when when you have similar value
                    while(l < r && nums[l] == nums[l - 1]){
                        l++;
                    }
                }
            }

        }

        return ans;
    }

    // Dividing by 2 is the same as >> 1 (right shifted by 1)
    vector<int> countBits(int n){
        vector<int> dp(n + 1);
        dp[0] = 0;

        for(int i = 1; i <= n; ++i){
            dp[i] = dp[i / 2] + i % 2;
        }
        return dp;
    }

    uint32_t reverseBits(uint32_t n) {
        uint32_t res = 0; // 32 0's

        // n & (1 << i)) == 0 mane oi pos e set bit nai != 0 mane set bit ase
        for(int i = 0; i < 32; ++i){
            if((n & (1 << i)) != 0){ // != 0 dewa lagbe cz it can be anything like 00000..100
                res = res | (1 << (31 - i)); // firs er dike e replace korte chacci with set bit
            }
        }
        return res;
    }

    // addition is simlar to xor and bitwise left shifting
    // Java use 2's complement to represent neg numbers and take care of negative numbers for you
    // neg nie alada chinta kora lage na
    int getSum(int a, int b) {
        while(b != 0){
            int tmp = a;
            a = a ^ b;
            b = (tmp & b) << 1; // b is the carry representative
        }
        return a;
    }

    int climbStairs(int n){
        int dp[n + 1];
        dp[0] = dp[1] = 1;

        for(int i = 2; i <= n; ++i){
            dp[i] = dp[i - 1] + dp[i - 2];
        }
        return dp[n];
    }

    int coinChange(vector<int>& coins, int amount){
        int len = coins.size();
        vector<vector<int>> dp(len + 1, vector<int>(amount + 1));

        for(int i = 0; i <= len; ++i){
            for(int j = 0; j <= amount; ++j){
                if(j == 0){
                    dp[i][j] = 0;
                }
                else if(i == 0){
                    dp[i][j] = INT16_MAX;
                }
                else if(coins[i - 1] > j){
                    dp[i][j] = dp[i - 1][j]; // not taking the ith coin because ile amount axceed
                    // kore
                }
                else{
                    dp[i][j] = min(dp[i - 1][j], 1 + dp[i][j - coins[i - 1]]); // not taking &
                    // taking
                }
            }
        }

        return dp[len][amount] >= INT16_MAX ? -1 : dp[len][amount];
    }


    int lengthOfLIS(vector<int>& nums) {
        int maxLen = 1;
        vector<int> LIS(nums.size(), 1);

        for(int i = 1; i < nums.size(); ++i){
            for(int j = 0; j < i; ++j){
                if(nums[i] > nums[j] && LIS[j] + 1 > LIS[i]){
                    LIS[i] = LIS[j] + 1;
                    maxLen = max(maxLen, LIS[i]);
                }
            }
        }
        return maxLen;
    }

    // Find ALL LIS
    // find all increasing subsets/subsequences
    void findAllLISDriver(vector<int>& nums, int idx, vector<vector<int>>& ans,
                                   vector<int>& tmp){
        // only increasing branch ke grow korte dibo
        if(tmp.size() >= 2 && tmp[tmp.size() - 1] <= tmp[tmp.size() - 2]){
            return; // cut the branch
        }
        if(idx == nums.size()){
            ans.push_back(tmp);
            return;
        }

        tmp.push_back(nums[idx]);
        findAllLISDriver(nums, idx + 1, ans, tmp);
        tmp.pop_back();
        findAllLISDriver(nums, idx + 1, ans, tmp);
    }

    vector<vector<int>> findAllLIS(vector<int>& nums){
        vector<vector<int>> ans;
        vector<int> tmp;
        findAllLISDriver(nums, 0, ans, tmp);
        return ans;
    }

    string make_key(int i, int j){
        return to_string(i) + "*" + to_string(j);
    }
    int longestCommonSubsequenceUtil(string text1, string text2, int i, int j,
                                     unordered_map<string, int>& mem) {
        if(i == text1.size() || j == text2.size()){
            return 0;
        }
        string key = make_key(i, j);
        if(mem.count(key) == 1){
            return mem[key];
        }
        else if(text1[i] == text2[j]){
            return mem[key] = 1 + longestCommonSubsequenceUtil(text1, text2, i + 1, j + 1, mem);
        }
        else{
            return mem[key] = max(longestCommonSubsequenceUtil(text1, text2, i + 1, j, mem),
                       longestCommonSubsequenceUtil(text1, text2, i, j + 1, mem));
        }
    }

    int longestCommonSubsequence(string text1, string text2) {
        unordered_map<string, int> mem;
        return longestCommonSubsequenceUtil(text1, text2, 0, 0, mem);
    }


    class Solution {
    private:
        string str;
        vector<string> ans;
        unordered_set<string> hash;

        void wordBreakUtil(int i, string tmp){
            if(i == str.size()){
                tmp.pop_back(); // popping an extra space
                ans.push_back(tmp);
                return;
            }

            for(int j = i + 1; j <= str.size(); ++j){
                string subS = str.substr(i, j - i);
                if(hash.count(subS) == 1){
                    wordBreakUtil(j, tmp + subS + " ");
                }
            }
        }
    public:
        vector<string> wordBreak(string s, vector<string>& wordDict) {
            str = s;
            for(auto &s: wordDict){
                hash.insert(s);
            }

            wordBreakUtil(0, "");
            return ans;
        }
    };


    void combinationSumUtil(vector<vector<int>>& ans, vector<int>& candidates, int i, int target,
                            vector<int>& tmp){
        if(target == 0){
            ans.push_back(tmp);
            return;
        }
        else if(i == candidates.size()){
            return;
        }

        if(candidates[i] <= target){
            tmp.push_back(candidates[i]);
            combinationSumUtil(ans, candidates, i, target - candidates[i], tmp); // including
            tmp.pop_back();
        }
        combinationSumUtil(ans, candidates, i + 1, target, tmp); // not including
    }

    vector<vector<int>> combinationSum(vector<int>& candidates, int target){
        vector<vector<int>> ans;
        vector<int> tmp;
        combinationSumUtil(ans, candidates, 0, target, tmp);
        return ans;
    }

    void combinationSum2Util(vector<vector<int>>& ans, vector<int>& candidates, int i, int target,
                            vector<int>& tmp){
        if(target == 0){
            ans.push_back(tmp);
            return;
        }
        else if(i == candidates.size()){
            return;
        }

        if(i < candidates.size() && candidates[i] <= target){
            tmp.push_back(candidates[i]);
            combinationSum2Util(ans, candidates, i + 1, target - candidates[i], tmp); //
            // including only once
            tmp.pop_back();
        }

        i++;
        while(candidates[i] == candidates[i - 1]){
            i++;
        }
        combinationSum2Util(ans, candidates, i, target, tmp); // not including
    }

    vector<vector<int>> combinationSum2(vector<int>& candidates, int target){
        sort(candidates.begin(), candidates.end());
        vector<vector<int>> ans;
        vector<int> tmp;
        combinationSum2Util(ans, candidates, 0, target, tmp);
        return ans;
    }

    int rob(vector<int>& nums){
        // top down dp approach
        int len = nums.size();
        vector<int> dp(len);
        dp[len - 1] = nums[len - 1]; // max profit possible to make starting from last house is
        // just to rob it
        if(len == 1){
            return dp[0];
        }
        dp[len - 2] = max(nums[len - 2], nums[len - 1]);

        for(int i = len - 3; i >= 0; ++i){
            dp[i] = max(nums[i] + dp[i + 2], dp[i + 1]); // rob cur and add max profit possible
            // starting from i + 2 // Don't rob cur
        }
        return dp[0];
    }



    int robIIUtil(vector<int>& nums, int i, int flag, unordered_map<string, int>& mem){
        if((i == nums.size() - 1 && flag == 1) || i >= nums.size()){ // flag turned on so first elem new hoise can't take
            // last one
            return 0;
        }
        else if(i == nums.size() - 1 && flag == 0){
            return nums[i];
        }

        string key = make_key(i, flag);
        if(mem.count(key) == 1){
            return mem[key];
        }
        else if(i == 0){ // i == 0 first case e amra flag nie deal korbo so eta exceptional case
            return mem[key] = max(nums[i] + robIIUtil(nums, i + 2, 1, mem), robIIUtil(nums, i +
                                                                                            1, 0, mem));
            // taking first and not taking first
            // first nile turning flag = 1
        }
        else{ // normal case flag r je value ase thakbe we will just pass it forward
            return mem[key] = max(nums[i] + robIIUtil(nums, i + 2, flag, mem), robIIUtil(nums, i +
            1, flag, mem)); // taking first and not taking first
        }
    }


    int robII(vector<int>& nums){
        unordered_map<string, int> mem;
        return robIIUtil(nums, 0, 0, mem);
    }

    int robIITab(vector<int>& nums){
        // top down dp approach
        int len = nums.size();
        vector<int> dp(len);
        dp[len - 1] = nums[len - 1]; // max profit possible to make starting from last house is
        // just to rob it
        if(len == 1){
            return dp[0];
        }
        dp[len - 2] = max(nums[len - 2], nums[len - 1]);

        for(int i = len - 3; i >= 0; ++i){
            dp[i] = max(nums[i] + dp[i + 2], dp[i + 1]); // rob cur and add max profit possible
            // starting from i + 2 // Don't rob cur
        }
        return dp[0];
    }

    // O(n) solution cz ekhane string size max can be 2
    // so substr O(1)
    // and bairer loop o O(1)
    // Word break e we went through all substrings that's why O(n^2 * n) = O(n^3)
    int numDecodingsUtil(string& s, int i, unordered_map<string, char>& hashMap, vector<int>&
            mem){
        if(i == s.size()){
            return 1;
        }
        else if(i > s.size()){
            return 0;
        }

        if(mem[i] != -1){
            return mem[i];
        }

        int count = 0;
        for(int j = i + 1; j <= i + 2; ++j){ // "226" j can only take 2 places front of i // max
            // string size limit restricted to 2
            if(hashMap.count(s.substr(i, j - i)) == 1){
                count += numDecodingsUtil(s, j, hashMap, mem);
            }
        }
        return mem[i] = count;
    }

    int numDecodings(string s){
        vector<int> mem(s.size(), -1);
        unordered_map<string, char> hashMap;

        for(int i = 0; i <= 25; ++i){
            hashMap[to_string(i + 1)] = 'A' + i;
        }

        return numDecodingsUtil(s, 0, hashMap, mem);
    }

    bool canJump(vector<int>& nums){
        int postIdx = nums.size() - 1;

        for(int i = nums.size() - 2; i >= 0; --i){
            if(i + nums[i] >= postIdx){ // postIdx can be reached
                postIdx = i + nums[i];
            }
        }

        return postIdx == 0; // the final value must be 0 for it to be possible to reach the end
    }

    // visited can be 0, 1, 2
    // 0 unprocessed, 1 processing, 2 processed
    bool isCyclicDirected(vector<vector<int>>& adj, vector<int> visited, int cur){
        visited[cur] = 1; // processing
        for(auto &n: adj[cur]){
            if(visited[n] == 1){ // processing obosthai pewe gesi mane cycle ase
                return true;
            }
            else if(visited[n] == 0){
                if(isCyclicDirected(adj, visited, n)){
                    return true;
                }
            }
        }
        visited[cur] = 2; // processed
        return false;
    }

    // visited 0 & 1
    bool isCyclicUndirected(vector<vector<int>>& adj, vector<int>& visited, int cur, int prev){
        visited[cur] = 1;

        for(auto &n: adj[cur]){
            if(visited[n] == 1 && n != prev){
                return true;
            }
            else if(visited[n] == 0){
                if(isCyclicUndirected(adj, visited, n, cur)){
                    return true;
                }
            }
        }
        return false; // no cycle
    }



    bool canFinish(int numCourses, vector<vector<int>>& prerequisites) {
        vector<vector<int>> adj(numCourses);

        for(auto &e: prerequisites){
            adj[e[1]].push_back(e[0]); // directed
        }

        vector<int> visited(numCourses, 0); // unprocessed
        for(int i = 0; i < numCourses; ++i){
            if(visited[i] == 0){
                if(isCyclicDirected(adj, visited, i)){
                    return false; // cyclic so not possible
                }
            }
        }
        return true; // no cycle
    }

    // 1 2 3 4 100 200
    int longestConsecutive(vector<int>& nums) {
        unordered_set<int> s;
        for(int &i: nums){
            s.insert(i);
        }

        int maxLen = 1;
        for(int &i: nums){
            if(s.count(i - 1) == 0){ // i is the beginning of a new seq
                int len = 1;
                while(s.count(i + 1) == 1){
                    i++;
                    len++;
                }
                maxLen = max(maxLen, len);
            }
        }
        return maxLen;
    }

    vector<vector<int>> insert(vector<vector<int>>& intervals, vector<int>& newInterval){
        vector<vector<int>> ans;
        bool flag = true;
        for(int i = 0; i < intervals.size(); ++i){
            if(newInterval[1] < intervals[i][0]){ // newEnd < start
                ans.push_back(newInterval); // will happen once
                flag = false;
                // after the new interval is pushed mane copy paste rest of the list into ans as
                // it is
                while(i < intervals.size()){
                    ans.push_back(intervals[i]);
                    i++;
                }
                return ans;
            }
            else if(newInterval[0] > intervals[i][1]){
                ans.push_back(intervals[i]);
            } // newStart > end
            else{
                // merge and change newInterval
                newInterval = vector<int> ({min(newInterval[0], intervals[i][0]), max
                                            (newInterval[1], intervals[i][1])});
            }
        }

        if(flag){
            ans.push_back(newInterval);
        }
        return ans;
    }

    vector<vector<int>> merge(vector<vector<int>>& intervals){
        sort(intervals.begin(), intervals.end()); // must be sorted for the merging process
        vector<vector<int>> ans;
        ans.push_back(intervals[0]); // first r ta dhukalam

        for(int i = 1; i < intervals.size(); ++i){
            if(intervals[i][0] <= ans.back()[1]){ // is the cur one starting before prev one ends
                // this one condition is enough since we have the array sorted
                // left boundary already thik ase since array sorted
                ans.back()[1] = max(ans.back()[1], intervals[i][1]);
            }
            else{
                ans.push_back(intervals[i]);
            }
        }
        return ans;
    }

    int eraseOverlapIntervals(vector<vector<int>>& intervals) {
        int minErase = 0;
        sort(intervals.begin(), intervals.end());

        // we will be greedy whenever we found a overlap to choose which interval to erase to
        // minimize our chance of future overlapping
        // the interval which ends later needs to be removed as that increases the chance of
        // future overlapping
        int endOfLast = intervals[0][1];
        for(int i = 1; i < intervals.size(); ++i){
            if(intervals[i][0] < endOfLast){ // cur one starts before the prev one ends
                // overlapping condition here
                minErase++; // ekta erase kora lagbe at least
                endOfLast = min(endOfLast, intervals[i][1]); // greedy decision to optimize
                // future erasings
            }
            else{
                endOfLast = intervals[i][1];
            }
        }

        return minErase;
    }

    bool canAttendMeetings(vector<vector<int>>& intervals){
        if(intervals.size() == 0){
            return true;
        }
        sort(intervals.begin(), intervals.end());
        int lastOfPrev = intervals[0][1];

        for(int i = 1; i < intervals.size(); ++i){
            if(intervals[i][0] < lastOfPrev){ // overlapping
                return false;
            }
            else{
                lastOfPrev = intervals[i][1];
            }
        }
        return true; // no overlapping
    }

    // the min number of conference room will the maximum number of overlapping intervals at any
    // one point in time
    // sorting + heap(based on ending time)
    int minMeetingRooms(vector<vector<int>>& intervals){
        sort(intervals.begin(), intervals.end());
        priority_queue<int, vector<int>, greater<int>> minHeap;
        minHeap.push(intervals[0][1]);
        int minRoom = 1;

        for(int i = 1; i < intervals.size(); ++i){
            if(intervals[i][0] < minHeap.top()){
                // overlapping so need to assign a new room
                minRoom++;
                // whenever assigned a new room push its ending time to minHeap
                minHeap.push(intervals[i][1]);
            }
            else{ // latest free room use korbe intervals[i] meeting
                minHeap.pop(); // older meeting end time ber koro
                minHeap.push(intervals[i][1]); //notun ta dhukao
            }
        }
        return minRoom;
    }

    // sliding window
    int lengthOfLongestSubstring(string s){
        if(s.empty()){
            return 0;
        }
        int l = 0, r = 1;
        int maxLen = 1;

        unordered_set<char> hashSet;
        hashSet.insert(s[l]);

        while(r < s.size()){
            while(hashSet.count(s[r]) == 1){ // already ase set e
                // start pushing left boundary
                hashSet.erase(s[l]);
                l++;
            }
            hashSet.insert(s[r]);
            maxLen = max(maxLen, r - l + 1);
            r++;
        }
        return maxLen;
    }

    // return the max occurance
    int maxFreq(unordered_map<char, int>& hashMap){
        int maxF = 0;
        for(auto &p: hashMap){
            maxF = max(maxF, p.second);
        }
        return maxF;
    }

    // in a particular substring, you want to make things like the maximum freq char
    // that why you need to replace windowLen - maxFreq to get the number of required replacements
    int characterReplacement(string s, int k) {
        int maxLen = 0;
        int l = 0, r = 0;
        unordered_map<char, int> hashMap;
        hashMap[s[l]]++; // entered the first char

        while(r < s.size()){
            int windowLen = r - l + 1;
            if(windowLen - maxFreq(hashMap) <= k){ // valid window
                maxLen = max(maxLen, windowLen);
                // push the right boundary to the right
                r++;
                hashMap[s[r]]++;
            }
            else{
                // push left boundary
                hashMap[s[l]]--;
                l++;
            }
        }
        return maxLen;
    }

    bool isValid(unordered_map<char, int>& map1, unordered_map<char, int>& map2){
        for(auto &p: map1){
            char key = p.first;
            if(map2[key] < map1[key]){
                return false;
            }
        }
        return true;
    }

    string minWindow(string s, string t) {
        if(s.size() < t.size()){
            return "";
        }
        else if(s.size() == 1){
            if(s == t){
                return s;
            }
            else{
                return "";
            }
        }

        unordered_map<char, int> map1;
        unordered_map<char, int> map2;
        for(char &c: t){
            map1[c]++;
        }

        int l = 0, r = 0;
        int minLen = INT16_MAX;
        int minL = 0;
        if(map1.count(s[r]) == 1){
            map2[s[r]]++;
        }

        while(r < s.size()){
            if(!isValid(map1, map2)){ // push right
                r++;
                if(map1.count(s[r]) == 1){
                    map2[s[r]]++;
                }
            }
            else{
                // push left
                int window = r - l + 1;
                if(window < minLen){
                    minLen = window;
                    minL = l;
                }
                if(map1.count(s[l]) == 1){
                    map2[s[l]]--;
                }
                l++;
            }
        }
        return minLen == INT16_MAX ? "" : s.substr(minL, minLen);
    }

    struct hashFunction
    {
        size_t operator()(const vector<int>
                          &myVector) const
        {
            std::hash<int> hasher;
            size_t answer = 0;

            for (int i : myVector)
            {
                answer ^= hasher(i) + 0x9e3779b9 +
                          (answer << 6) + (answer >> 2);
            }
            return answer;
        }
    };


    vector<vector<string>> groupAnagrams(vector<string>& strs) {
        vector<vector<string>> ans;
        // we can actually use a vector as a key
        // keys will be like 10201331..3329 (26) ta count thakbe
        unordered_map<vector<int>, vector<string>, hashFunction> hashMap;

        for(auto &s: strs){
            vector<int> count(26, 0); // for each 26 letters
            for(auto &c: s){
                count[c - 'a']++; // 'a' -'a' = 0
            }
            hashMap[count].push_back(s); // hashmap[count] is a list of strings
        }

        for(auto &p: hashMap){
            ans.push_back(p.second);
        }
        return ans;
    }

    // string parsing
    bool isPalindrome(string s) {
        int l = 0, r = s.size() - 1;
        while(l <= r && r < s.size()){
            while(l < s.size() && !isalnum(s[l])){
                l++;
            }
            while(r < s.size() && !isalnum(s[r])){
                r--;
            }

            if(l < s.size() && r < s.size() && toupper(s[l]) != toupper(s[r])){
                return false;
            }
            else{
                l++;
                r--;
            }
        }
        return true;
    }

    // naive solution is O(n^3)
    // dp is O(n^2)
    // dp of all possible substrings
    // s[1:3] = s[2:2] && s[1] == s[3]
    // column wise fillup kore jete hbe
    /*   0 1 2 3
     * 0
     * 1
     * 2
     * 3
     */
    string longestPalindrome(string s){
        int maxLen = 1;
        int maxL;
        maxL = -1;
        int rows, cols;
        rows = cols = s.size();
        vector<vector<bool>> dp(rows, vector<bool>(2));

        // column by column traversal
        // we always need the column before so 2 ta coumn die kaaj chalano jay
        // dp[i : j] i -> row j -> col
        for(int j = 0; j < cols; ++j){ // right boundary // for each col
            for(int i = 0; i < rows; ++i){ // left boundary // for each row
                if(i >= j){
                    dp[i][j % 2] = true;
                }
                else{
                    dp[i][j % 2] = s[i] == s[j] && dp[i + 1][(j - 1) % 2];
                }

                int windowLen = j - i + 1;
                if(dp[i][j % 2] && windowLen >= maxLen){
                    maxLen = windowLen;
                    maxL = i;
                }
            }
        }

        return s.substr(maxL, maxLen);
    }

    // count the number of palinromic substrings
    int countSubstrings(string s){
        int count = 0;
        int rows, cols;
        rows = cols = s.size();
        vector<vector<bool>> dp(rows, vector<bool>(2));

        // column by column traversal
        // we always need the column before so 2 ta coumn die kaaj chalano jay
        // dp[i : j] i -> row j -> col
        for(int j = 0; j < cols; ++j){ // right boundary // for each col
            for(int i = 0; i < rows; ++i){ // left boundary // for each row
                if(i >= j){
                    dp[i][j % 2] = true;
                }
                else{
                    dp[i][j % 2] = s[i] == s[j] && dp[i + 1][(j - 1) % 2];
                }

                if(j >= i && dp[i][j % 2]){ // j >= i for a valid substring
                    count++;
                }
            }
        }

        return count;
    }


    // Bucket sort
    vector<int> topKFrequent(vector<int>& nums, int k) {
        vector<int> ans;
        unordered_map<int, int> hashMap;
        for(int &i: nums){
            hashMap[i]++;
        }

        // make a count array for bucket sort
        // the max value of k can be nums.size()
        vector<vector<int>> count(nums.size() + 1);
        // will do count -> num mapping (ulta)
        for(auto &p: hashMap){
            count[p.second].push_back(p.first);
        }

        // iterate from the end of count array
        for(int i = nums.size(); k > 0 && i >= 1; --i){
            if(!count[i].empty()){
                for(int j = 0; k > 0 && j < count[i].size(); ++j){
                    ans.push_back(count[i][j]);
                    k--;
                }
            }
        }
        return ans;
    }
};




// ROGUE CLASSES

class MedianFinder {
private:
    vector<int> v;
public:
    // naive solution would be inserting like insertion sort
    MedianFinder() {}

    void addNum(int num) {
        v.push_back(num);
        int i = v.size() - 1;
        while(i >= 1 && v[i] < v[i - 1]){
            swap(v[i], v[i - 1]);
            i--;
        }
    }

    double findMedian() {
        int len = v.size();
        if(len % 2 == 0){
            int f = v[(len / 2) - 1];
            int s = v[len / 2];
            return (double)(f + s) / 2;
        }
        else{
            return v[len / 2];
        }
    }
};

class MedianFinderOp {
private:
    priority_queue<int> firstHeap; // max heap
    priority_queue<int, vector<int>, greater<int>> secondHeap; // min heap
public:
    // naive solution would be inserting like insertion sort
    MedianFinderOp() {}

    void addNum(int num) {
        firstHeap.push(num); // by default

        int diff = (int)(firstHeap.size() - secondHeap.size());
        int balance = abs(diff);

        while(balance > 1){
            secondHeap.push(firstHeap.top());
            firstHeap.pop();

            diff = (int)(firstHeap.size() - secondHeap.size());
            balance = abs(diff);
        }

        int firstHalfMax = firstHeap.empty() ? INT32_MIN : firstHeap.top();
        int secondHalfMin = secondHeap.empty() ? INT32_MAX : secondHeap.top();
        while(firstHalfMax > secondHalfMin){
            firstHeap.pop();
            secondHeap.push(firstHalfMax);

            firstHalfMax = firstHeap.top();
            secondHalfMin = secondHeap.top();
        }
    }

    double findMedian() {
        if(firstHeap.size() > secondHeap.size()){
            return firstHeap.top();
        }
        else if(firstHeap.size() < secondHeap.size()){
            return secondHeap.top();
        }
        else{
            return double(firstHeap.top() + secondHeap.top()) / 2;
        }
    }
};

class TreeProblems {
public:
    int maxDepth(TreeNode* root) {
        if(root == nullptr){
            return 0;
        }
        return 1 + max(maxDepth(root->left) , maxDepth(root->right));
    }

    bool isSameTree(TreeNode* p, TreeNode* q){
        if(p == nullptr && q == nullptr){
            return true;
        }
        else if(p == nullptr || q == nullptr){ // we are sure if one is null the other is not at
            // this point
            return false;
        }
        else{
            return p->val == q->val && isSameTree(p->left, q->left) &&
            isSameTree(p->right, q->right);
        }
    }

    TreeNode* invertTree(TreeNode* root){
        if(root == nullptr){
            return nullptr;
        }
        // root != nullptr
        TreeNode* newRoot = new TreeNode(root->val);
        newRoot->right = invertTree(root->left);
        newRoot->left = invertTree(root->right);
        return newRoot;
    }

    // the idea is we will return the maximum possible from each node without diverging at that
    // point i.e. choosing either left or right option
    int maxPathSumUtil(TreeNode* root, int& globalMaxima){
        if(root == nullptr){
            return 0;
        }
        // root can't be nullptr
        int maxLeftWithoutDiverging = maxPathSumUtil(root->left, globalMaxima);
        int maxRightWithoutDiverging = maxPathSumUtil(root->right, globalMaxima);

        // calculating global maxima considering the path can diverge at the current node
        if(maxLeftWithoutDiverging < 0){
            globalMaxima = max(globalMaxima, root->val + maxRightWithoutDiverging);
        }
        vector<int> tmp({globalMaxima, root->val, maxLeftWithoutDiverging + root->val,
                           root->val + maxRightWithoutDiverging, maxLeftWithoutDiverging +
                           root->val + maxRightWithoutDiverging});
        globalMaxima = *max_element(tmp.begin(), tmp.end()); // returns the iterator to the max
        // element

        // but returning the maximum path possible without divergence at the current node
        return root->val + max(maxLeftWithoutDiverging, max(maxRightWithoutDiverging, 0));
    }

    int maxPathSum(TreeNode* root) {
        int globalMaxima = INT16_MIN;
        maxPathSumUtil(root, globalMaxima);
        return globalMaxima;
    }

    vector<vector<int>> levelOrder(TreeNode* root){
        vector<vector<int>> ans;
        if(root == nullptr){
            return ans;
        }
        queue<TreeNode*> q;
        q.push(root);

        while(!q.empty()){
            int levelLen = q.size();
            vector<int> tmp;
            for(int i = 0; i < levelLen; ++i){
                TreeNode* cur = q.front();
                q.pop();
                tmp.push_back(cur->val);
                if(cur->left){
                    q.push(cur->left);
                }
                if(cur->right){
                    q.push(cur->right);
                }
            }
            ans.push_back(tmp);
        }
        return ans;
    }

    // gonna be a pre-order traversal
    // Encodes a tree to a single string.
    string serialize(TreeNode* root) {
        if(root == nullptr){
            return "n,";
        }
        string ans;
        ans += to_string(root->val) + "," + serialize(root->left) + serialize(root->right);
        return ans;
    }

    TreeNode* deserializeUtil(string data, int& i){
        string num;
        while(data[i] != ','){
            num += data[i];
            i++;
        }
        i++; // taking i to correct position for the next call

        if(num == "n"){
            return nullptr;
        }

        int numInt = stoi(num);
        // preorder building of tree
        TreeNode* root = new TreeNode(numInt);
        root->left = deserializeUtil(data, i);
        root->right = deserializeUtil(data, i);
        return root;
    }


    // Decodes your encoded data to tree.
    TreeNode* deserialize(string data) { // , delimitter
        int i = 0;
        return deserializeUtil(data, i);
    }

    bool isSubtree(TreeNode* root, TreeNode* subRoot){
        if(root == nullptr){
            return subRoot == nullptr;
        }

        // post order dfs so no extra work // O(n) visiting ever node only once
        return isSubtree(root->left, subRoot) || isSubtree(root->right, subRoot) ||
                isSameTree(root, subRoot);
    }

    TreeNode* buildTreeUtil(vector<int>& preorder, vector<int>& inorder, int p1, int p2, int i1,
                            int i2){
        if(p1 > p2 || i1 > i2){
            return nullptr;
        }
        TreeNode* root = new TreeNode(preorder[p1]);

        // build a hashmap to optimize this tep so that we can get the index of preorder[p1] in
        // inorder in O(1) times
        int i;
        for(i = i1; i <= i2; ++i){
            if(inorder[i] == preorder[p1]){
                break;
            }
        }
        // now i is at the correct idx of inorder
        int leftTreeSize = i - i1;
        int rightTreeSize = i2 - i;

        root->left = buildTreeUtil(preorder, inorder, p1 + 1, p1 + leftTreeSize, i1, i - 1);
        root->right = buildTreeUtil(preorder, inorder, p1 + leftTreeSize + 1, p2, i + 1, i2);
        return root;
    }

    TreeNode* buildTree(vector<int>& preorder, vector<int>& inorder) {
        return buildTreeUtil(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1);
    }

    bool isValidBSTUtil(TreeNode* root, int l, int h){
        if(root == nullptr){
            return true;
        }
        return l < root->val && root->val < h && isValidBSTUtil(root->left, l, root->val) &&
                isValidBSTUtil(root->right, root->val, h);
    }

    bool isValidBST(TreeNode* root) {
        return isValidBSTUtil(root, INT32_MIN, INT32_MAX);
    }

    // O(logn) solution for a balanced tree
    TreeNode* lowestCommonAncestor(TreeNode* root, TreeNode* p, TreeNode* q){
        if(root == p || root == q){
            return root;
        }
        else if(p->val < root->val && root->val < q->val ||
                q->val < root->val && root->val < p->val){
            return root;
        }
        else if(p->val < root->val && q->val < root->val){
            // reduced search space to half of the tree
            return lowestCommonAncestor(root->left, p, q);
        }
        else{
            return lowestCommonAncestor(root->right, p, q);
        }
    }

    // O(n)
    // it's like searching the nodes
    TreeNode* lowestCommonAncestorBinTree(TreeNode* root, TreeNode* p, TreeNode* q){
        if(root == nullptr){
            return nullptr;
        }
        else if(root == p || root == q){
            return root;
        }

        TreeNode* leftN = lowestCommonAncestorBinTree(root->left, p, q); // like a search
        TreeNode* rightN = lowestCommonAncestorBinTree(root->right, p, q); // like a search

        if(leftN && rightN){
            return root;
        }
        else{
            return leftN == nullptr ? rightN : leftN;
        }
    }
};

class TrieNode{
public:
    bool isEnd;
    unordered_map<char, TrieNode*> hashMap;
    TrieNode(bool flag): isEnd(false){}
};

class Trie {
    TrieNode* root;
public:
    Trie() {
        root = new TrieNode(false);
    }

    void insert(string word) {
        TrieNode* tmp(root);
        int i;
        for(i = 0; i < word.size(); ++i){
            if((tmp->hashMap).count(word[i]) == 0){ // nei
                (tmp->hashMap)[word[i]] = new TrieNode(false);
            }
            tmp = (tmp->hashMap)[word[i]];
        }
        tmp->isEnd = true;
    }

    bool search(string word) {
        TrieNode* tmp(root);
        int i;
        for(i = 0; i < word.size(); ++i){
            if((tmp->hashMap).count(word[i]) == 0){ // nei
                return false; // khuje paini
            }
            tmp = (tmp->hashMap)[word[i]];
        }
        return tmp->isEnd;
    }

    bool startsWith(string prefix) {
        TrieNode* tmp(root);
        int i;
        for(i = 0; i < prefix.size(); ++i){
            if((tmp->hashMap).count(prefix[i]) == 0){ // nei
                return false; // khuje paini
            }
            tmp = (tmp->hashMap)[prefix[i]];
        }
        return true; // isEnd doesn't matter prefix r shb word thake manei true
    }
};

class WordDictionary {
    TrieNode* root;
public:
    WordDictionary() {
        root = new TrieNode(false);
    }

    void addWord(string word) {
        TrieNode* tmp(root);
        int i;
        for(i = 0; i < word.size(); ++i){
            if((tmp->hashMap).count(word[i]) == 0){ // nei
                (tmp->hashMap)[word[i]] = new TrieNode(false);
            }
            tmp = (tmp->hashMap)[word[i]];
        }
        tmp->isEnd = true;
    }

    bool searchUtil(string word, int s, TrieNode* root){
        TrieNode* tmp(root);
        int i;
        for(i = s; i < word.size(); ++i){
            if(word[i] == '.'){
                for(auto &p: tmp->hashMap){
                    if(searchUtil(word, i + 1, p.second)){
                        return true;
                    }
                }
                return false;
            }
            else{
                if((tmp->hashMap).count(word[i]) == 0){ // nei
                    return false; // khuje paini
                }
                tmp = (tmp->hashMap)[word[i]];
            }
        }
        return tmp->isEnd;
    }

    bool search(string word) {
        return searchUtil(word, 0, root);
    }
};

class WordSearchII {
private:
    int rows;
    int cols;
    Trie t;
    vector<string> ans;
    unordered_set<string> ansSet;
    vector<vector<bool>> boardVisited;
    void backtracking(vector<vector<char>>& board, string tmp, int r, int c){
        if(t.search(tmp)){
            ansSet.insert(tmp);
        }
        boardVisited[r][c] = true;
        for(auto &p: DIR){
            int newRow = r + p[0];
            int newCol = c + p[1];
            if(newRow >= 0 && newRow < rows && newCol >= 0 && newCol < cols &&
            !boardVisited[newRow][newCol] && t.startsWith(tmp + board[newRow][newCol])){
                backtracking(board, tmp + board[newRow][newCol], newRow, newCol);
            }
        }
        boardVisited[r][c] = false;
    }
public:
    vector<string> findWords(vector<vector<char>>& board, vector<string>& words) {
        for(string &s: words){
            t.insert(s);
        }
        rows = board.size();
        cols = board[0].size();
        boardVisited.resize(rows, vector<bool>(cols, false));

        for(int r = 0; r < rows; ++r){
            for(int c = 0; c < cols; ++c){
                if(t.startsWith(string(1, board[r][c]))){
                    backtracking(board, string(1, board[r][c]), r, c);
                }
            }
        }
        for(auto &e: ansSet){
            ans.push_back(e);
        }
        return ans;
    }
};

class Codec {
public:
    // Encodes a list of strings to a single string.
    string encode(vector<string>& strs) {
        string encodedStr = "";
        for(string &s: strs){
            encodedStr += to_string(s.size()) + "#" + s;
        }
        return encodedStr;
    }

    // Decodes a single string to a list of strings.
    vector<string> decode(string s) {
        vector<string> output;

        int i = 0;
        while(i < s.size()){
            string len = "";
            while(i < s.size() && s[i] != '#'){
                len += s[i];
                ++i;
            }
            // s[i] == '#' which is delimeter
            int lenInt = stoi(len);
            string str = s.substr(i + 1, lenInt); // i is at delimitter
            output.push_back(str);
            i = i + lenInt + 1;
        }
        return output;
    }
};

class WordSearch {
private:
    int rows;
    int cols;
public:
    // you can use the board itself as visited by using non-alphabetic chars to mark visited places
    // then after backtracking undoing what you have done
    bool existUtil(vector<vector<char>>& board, vector<vector<int>>& boardVisit, string& word, int
    i, int r, int c){
        if(i == word.size() - 1){
            return true;
        }

        // mark cur as visited
        boardVisit[r][c] = 1;
        for(auto &d: DIR){
            int newR = r + d[0];
            int newC = c + d[1];
            if(newR >= 0 && newR < rows && newC >= 0 && newC < cols && boardVisit[newR][newC] ==
            0 && board[newR][newC] == word[i + 1]){
                if(existUtil(board, boardVisit, word, i + 1, newR, newC)){
                    return true;
                }
            }
        }
        boardVisit[r][c] = 0; // ei func theke ber howar shomoy homoy unvisited howe ber hbe shb
        return false; // kono direction e jawei match pay nai
    }

    // Word Search I
    bool exist(vector<vector<char>>& board, string word) {
        rows = board.size();
        cols = board[0].size();

        vector<vector<int>> boardVisit(rows, vector<int>(cols, 0));
        for(int r = 0; r < rows; ++r){
            for(int c = 0; c < cols; ++c){
                if(board[r][c] == word[0]){
                    if(existUtil(board, boardVisit, word, 0, r, c)){
                        return true;
                    }
                }
            }
        }
        return false;
    }
};

class AlienOrder {
public:
    // visited can be 0, 1, 2
    // 0 unprocessed, 1 processing, 2 processed // this will be mapped here
    bool isCyclicDirected(unordered_map<char, vector<char>>& adj, unordered_map<char, int>&
            visited, char cur){
        visited[cur] = 1; // processing
        for(auto &n: adj[cur]){
            if(visited[n] == 1){ // processing obosthai pewe gesi mane cycle ase
                return true;
            }
            else if(visited[n] == 0){
                if(isCyclicDirected(adj, visited, n)){
                    return true;
                }
            }
        }
        visited[cur] = 2; // processed
        return false;
    }

    void topoSort(unordered_map<char, vector<char>>& adj, unordered_map<char, int>& visited,
                  stack<char>& s, char cur){
        visited[cur] = 1; // visited
        for(auto &n: adj[cur]){
            if(visited[n] == 0){
                topoSort(adj, visited, s, n);
            }
        }
        s.push(cur);
    }

    string alienOrder(vector<string>& words) {
        unordered_map<char, vector<char>> adj;
        unordered_map<char, int> visited1;
        unordered_map<char, int> visited2;

        for(auto &word: words){
            for(auto &c: word){
                adj.insert(make_pair(c, vector<char>()));
            }
        }

        for(int i = 0; i < words.size() - 1; ++i){
            int m = 0;
            bool flag = true;
            while(m < words[i].size() && m < words[i + 1].size()){
                if(words[i][m] == words[i + 1][m]){
                    m++;
                    continue;
                }
                else{
                    adj[words[i][m]].push_back(words[i + 1][m]);
                    flag = false;
                    break;
                }
            }
            if(flag){ // abc and ab situation hoise // abc must come after ab
                if(words[i + 1].length() < words[i].length()){ // impossible to form
                    return "";
                }
            }
        }

        for(auto &p: adj){
            if(visited1[p.first] == 0){
                if(isCyclicDirected(adj, visited1, p.first)){
                    return ""; // cycle paisi topo sort not possible
                }
            }
        }

        // no cycle
        stack<char> s;
        for(auto &p: adj){
            if(visited2[p.first] == 0){
                topoSort(adj, visited2, s, p.first);
            }
        }

        string output = "";
        while(!s.empty()){
            output.push_back(s.top());
            s.pop();
        }
        return output;
    }
};

class NumIslands {
    int rows;
    int cols;
public:
    void dfs(vector<vector<char>>& grid, int r, int c){
        // mark as visited
        grid[r][c] = 0;
        for(auto &d: DIR){
            int newR = r + d[0];
            int newC = c + d[1];
            if(newR >= 0 && newR < rows && newC >= 0 && newC < cols &&
                grid[newR][newC] == '1'){ // not visited // 0 mane water or already visited
                dfs(grid, newR, newC);
            }
        }
    }

    // find numComponents of directed graph
    int numIslands(vector<vector<char>>& grid) {
        int numIslands = 0;
        rows = grid.size();
        cols = grid[0].size();

        // use the grid itself as a visited set
        // whenever you visit a land(1) make it 0 so it won't be visited again
        for(int r = 0; r < rows; ++r){
            for(int c = 0; c < cols; ++c){
                if(grid[r][c] == '1'){
                    numIslands++;
                    dfs(grid, r, c);
                }
            }
        }
        return numIslands;
    }
};

class PacificAtlantic{
private:
    typedef pair<int, int> Pair;

    int rows;
    int cols;
    struct hashFunction
    {
        size_t operator()(const pair<int ,
                int> &x) const
        {
            return x.first ^ x.second;
        }
    };

    bool isValid(int r, int c){
        return r >= 0 && r < rows && c >= 0 && c < cols;
    }

    void dfs(vector<vector<int>>& heights, int r, int c, int prevHeight, unordered_set<Pair, hashFunction>&
    visited){
        visited.insert(Pair(r, c));

        for(auto &d: DIR){
            int newR = r + d[0];
            int newC = c + d[1];
            if(isValid(newR, newC) && heights[newR][newC] >= prevHeight && visited.count(Pair(newR, newC)) == 0){
                dfs(heights, newR, newC, heights[r][c], visited);
            }
        }
    }

public:
    vector<vector<int>> pacificAtlantic(vector<vector<int>>& heights) {
        rows = heights.size();
        cols = heights[0].size();
        vector<vector<int>> ans;
        unordered_set<Pair, hashFunction> setP;
        unordered_set<Pair, hashFunction> setA;

        for(int c = 0; c < cols; ++c){
            if(setP.count(Pair(0, c)) == 0){
                dfs(heights, 0, c, 0, setP);
            }
            if(setA.count(Pair(rows - 1, c)) == 0){
                dfs(heights, rows - 1, 0, c, setA);
            }
        }

        for(int r = 0; r < rows; ++r){
            if(setP.count(Pair(r, 0)) == 0){
                dfs(heights, r, 0, 0, setP);
            }
            if(setA.count(Pair(r, cols - 1)) == 0){
                dfs(heights, r, cols - 1, 0, setA);
            }
        }

        for(auto &p: setP){
            if(setA.count(p) == 1){
                ans.push_back(vector<int>({p.first, p.second}));
            }
        }
        return ans;
    }
};

class CloneGraph {
    class Node {
    public:
        int val;
        vector<Node*> neighbors;
        Node() {
            val = 0;
            neighbors = vector<Node*>();
        }
        Node(int _val) {
            val = _val;
            neighbors = vector<Node*>();
        }
        Node(int _val, vector<Node*> _neighbors) {
            val = _val;
            neighbors = _neighbors;
        }
    };

    Node* DFS(Node* oldNode, unordered_map<Node*, Node*>& visited){
        // copy oldNode
        Node* newNode = new Node(oldNode->val);
        visited[oldNode] = newNode; // mark visited

        for(auto &n: oldNode->neighbors){
            if(visited.count(n) == 0){
                (newNode->neighbors).push_back(DFS(n, visited));
            }
            else{ // neighbor already has been made
                (newNode->neighbors).push_back(visited[n]);
            }
        }
        return newNode;
    }
public:
    Node* cloneGraph(Node* node) {
        if(node == nullptr){
            return nullptr;
        }
        unordered_map<Node*, Node*> visited;
        return DFS(node, visited);
    }
};




#endif //BLIND75_BLIND_H
