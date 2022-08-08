//
// Created by Tahsinul Haque Abir on 8/1/22.
//

// Also Contains Sliding Window problems

#ifndef ABIR_TWOPOINTERS_H
#define ABIR_TWOPOINTERS_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;
typedef pair<int, int> Pair;

class TwoPointers {
public:
    // One pass oslution
    vector<int> twoSum(vector<int>& nums, int target) {
        unordered_map<int, int> hash; // (k, v) will be (nums[i], i)
        for(int i = 0; i < nums.size(); ++i){
           int comp = target - nums[i];
           if(hash.count(comp) == 1){
               return vector<int>({i, hash[comp]}); // hash[comp] contains idx of the complement
           }
           hash[nums[i]] = i; // (nums[i], i)
        }
        return vector<int>();
    }

    // O(n) is a variation of the merging algorithm
    vector<int> sortedSquares(vector<int>& nums) {
        vector<int> ans;
        int i, j;
        for(i = 0; i < nums.size(); ++i){
            if(nums[i] >= 0){
                break;
            }
        }
        // now i is at the first positive index
        j = i - 1;
        // j is at the idx before

        // now merge the positive array going froward and the negative one going backwards
        while (i < nums.size() && j >= 0){
            if(nums[i] < abs(nums[j])){
                ans.push_back(pow(nums[i], 2));
                i++;
            }
            else if(nums[i] > abs(nums[j])){
                ans.push_back(pow(nums[j], 2));
                j--;
            }
            else{
                // equal // so push twice
                ans.push_back(pow(nums[i], 2));
                ans.push_back(pow(nums[i], 2));
                i++;
                j--;
            }
        }

        while(i < nums.size()){
            ans.push_back(pow(nums[i], 2));
            i++;
        }

        while(j >= 0){
            ans.push_back(pow(nums[j], 2));
            j--;
        }

        return ans;
    }

    // build the two output string and compare
    // time: O(M + N)
    bool backspaceCompare(string s, string t) {
        string s1, t1;
        // string builder always uses push and pop as they are O(1) never concatanate
        // O(n)
        for(char &c: s){
            if(c != '#'){
                s1.push_back(c);
            }
            else{
                if(!s1.empty()){
                    s1.pop_back();
                }
                // empty hole just skip
            }
        }

        // O(n)
        for(char &c: t){
            if(c != '#'){
                t1.push_back(c);
            }
            else{
                if(!t1.empty()){
                    t1.pop_back();
                }
                // empty hole just skip
            }
        }

        return s1 == t1; // O(n)
    }

    // O(n) sol just find the sum and do math
    // only one duplicate can occur any number of times
    // Floyd cycle detection algoirthm
    // find the entrance of the cycle
    int findDuplicate(vector<int>& nums) {
        int s = nums[0];
        int f = nums[0];

        do{
            f = nums[nums[f]];
            s = nums[s];
        } while (s != f);
        // s = f now
        // find the entrance of the cycle
        f = nums[0]; // keep s in its place
        while(s != f){
            s = nums[s];
            f = nums[f];
        }
        return s;
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
                    l++; // ekta  l++ tow amra noramlly kori
                    // the next loop is for time when when you have similar value
                    while(l < r && nums[l] == nums[l - 1]){
                        l++;
                    }
                }
            }
        }

        return ans;
    }

    // eta te same value nie pera nai
    // like same set duibar dhukleo pera nai
    int threeSumClosest(vector<int>& nums, int target) {
        sort(nums.begin(), nums.end());
        int minDiff = INT32_MAX; // minimise korte bolse so eta always max possible value dite hbe
        int closestSum = target;

        for(int i = 0; i < nums.size(); ++i){
            int first = nums[i]; // first
            // now you have to find to find l and r such that l + r = -f
            int l = i + 1, r = nums.size() - 1;

            while(l < r){
                int threeSum = first + nums[l] + nums[r];
                if(threeSum == target){
                    return target; // the closest value to target is it itself
                }
                else{
                    if(abs(threeSum - target) < minDiff){
                        minDiff = abs(threeSum - target);
                        closestSum = threeSum;
                    }
                }

                if(threeSum < target){
                    l++;
                }
                else if(threeSum > target){
                    r--;
                }
            }
        }
        return closestSum;
    }


    // all positive in nums bolei ei formula work kore
    int numSubarrayProductLessThanK(vector<int>& nums, int k) {
        if(k == 0){
            return 0;
        }
        int prod = 1;
        int result = 0;
        int l, r;
        l = r = 0;

        while(r < nums.size()){
            prod = prod * nums[r];
            while(l <= r && prod >= k){
                prod = prod / nums[l];
                l++;
            }
            if(l <= r){
                result += r - l + 1; // this is the formula to count it works try by using an example
            }
            r++;
        }
        return result;
    }

    // O(n^2) solution would use finding all the sub-arrays
    int minSubArrayLen(int target, vector<int>& nums) {
        // O(n) solution optimized
        int l, r;
        l = r = 0;

        int minLen = INT32_MAX;
        int sum = nums[0];
        while(l <= r && r < nums.size()){
            while(l <= r && sum >= target){
                minLen = min(minLen, r - l + 1);
                sum -= nums[l];
                l++;
            }

            if(sum < target){
                r++;
                if(r < nums.size()){
                    sum += nums[r];
                }
            }
        }
        return minLen == INT32_MAX ? 0 : minLen;
    }

    int totalFruit(vector<int>& fruits) {
        unordered_map<int, int> hash;
        int maxFruit = 1;
        int l, r;
        l = r = 0;
        hash[fruits[0]]++;

        while(l <= r && r < fruits.size()){
            while(hash.size() > 2){
                hash[fruits[l]]--;
                if(hash[fruits[l]] == 0){
                    hash.erase(fruits[l]);
                }
                l++;
            }

            maxFruit = max(maxFruit, r - l + 1);

            if(hash.size() <= 2){
                r++;
                if(r < fruits.size()){
                    hash[fruits[r]]++;
                }
            }
        }
        return maxFruit;
    }

    // O(26n) solution
    bool checkInclusion(string s1, string s2) {
        if(s1.size() > s2.size()){
            return false;
        }

        unordered_map<char, int> hash1; // for s1
        unordered_map<char, int> hash2; // for subs of s2

        // make a hash of s1
        for(char &c: s1){
            hash1[c]++;
        }

        int l, r;
        int windowLen = s1.size();
        l = 0;
        r = l + windowLen - 1;

        // put the first s2 sub in hash2
        for(int i = l; i <= r; ++i){
            hash2[s2[i]]++;
        }


        while(hash1 != hash2 && r < s2.size()){
            hash2[s2[l]]--;

            if(hash2[s2[l]] == 0){
                hash2.erase(s2[l]);
            }

            l++;
            r++;

            if(r < s2.size()){
                hash2[s2[r]]++;
            }
        }

        return hash1 == hash2;
    }

    int maxFreq(unordered_map<char, int>& map){
        int maxF = 0;
        for(auto &p: map){
            maxF = max(maxF, p.second);
        }
        return maxF;
    }


    int characterReplacement(string s, int k) {
        unordered_map<char, int> map; // as we slide through subset it will produce hash to give
        // us maxFreq
        map[s[0]]++;
        int l, r;
        l = r = 0;
        int maxLenSub = 0;

        // we can apply a trick to optimize the O(26n) solution by keeping maxFreq variable and
        // updating it everytime we do map[s[r]]++ and overestimating by not caring about map[s[l]]--
        while(r < s.size()){ // O(n)
            int windowLen = r - l + 1;
            int maxFq = maxFreq(map); // O(26)

            if(windowLen - maxFq <= k){ // valid window
                maxLenSub = max(maxLenSub, windowLen);
                r++;
                if(r < s.size()){
                    map[s[r]]++;
                }
            }
            else{
                map[s[l]]--;
                l++;
            }
        }
        return maxLenSub;
    }

    // the most optimized solu to this problem will be using a montonously decreasing queue.
    // Actually we need to use a deque for O(n)
    // but jehetu max bolse the first thought would be using a max heap keeping track of (idx, val)
    // and that's O(nlogn) pretty optimal
//    vector<int> maxSlidingWindow(vector<int>& nums, int k) {
//
//    }

    // Given a string s, find the length of the longest substring without repeating characters.
    int lengthOfLongestSubstring(string s) {
        // using a set seems intuitive
        unordered_set<char> hashSet;
        int maxLen = 0;
        int l, r;
        l = r = 0;

        while(r < s.size()){

            if(hashSet.count(s[r]) == 0){
                maxLen = max(maxLen, r - l + 1);
                hashSet.insert(s[r]);
                r++;
            }

            while(hashSet.count(s[r]) == 1){
                hashSet.erase(s[l]);
                l++;
            }
        }
        return maxLen;
    }
};


#endif //ABIR_TWOPOINTERS_H
