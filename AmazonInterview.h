//
// Created by Tahsinul Haque Abir on 8/19/22.
//

#ifndef BLIND75_AMAZONINTERVIEW_H
#define BLIND75_AMAZONINTERVIEW_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>
#include "Blind.h"

using namespace std;

class Node {
public:
    int val;
    Node* next;
    Node* random;

    Node(int _val) {
        val = _val;
        next = NULL;
        random = NULL;
    }
};

class AmazonInterview {
public:
    string intToRomanUtil(int num, unordered_map<int, string>& dp){
        if(dp.count(num) == 1){
            return dp[num];
        }

        if(1 < num && num < 4 || 5 < num && num < 9){
            return dp[num] = intToRomanUtil(num - 1, dp) + "I";
        }
        else if(num % 10 == 0 && (10 < num && num < 40 || 50 < num && num < 90)){
            return dp[num] = intToRomanUtil(num - 10, dp) + "X";
        }
        else if(num % 100 == 0 && (100 < num && num < 400 || 500 < num && num < 900)){
            return dp[num] = intToRomanUtil(num - 100, dp) + "C";
        }
        else if(num % 1000 == 0 && (1000 < num && num < 4000)){
            return dp[num] = intToRomanUtil(num - 1000, dp) + "M";
        }
        else{
            int placeVal = 0;
            string str = "";
            while(num != 0){
                str = intToRomanUtil((num % 10) * pow(10, placeVal), dp) + str;
                placeVal++;
                num /= 10;
            }
            return dp[num] = str;
        }
    }
    // 1 <= num <= 3999
    string intToRoman(int num) {
        unordered_map<int, string> dp;

        dp[1] = "I";
        dp[4] = "IV";
        dp[5] = "V";
        dp[9] = "IX";
        dp[10] = "X";
        dp[40] = "XL";
        dp[50] = "L";
        dp[90] = "XC";
        dp[100] = "C";
        dp[400] = "CD";
        dp[500] = "D";
        dp[900] = "CM";
        dp[1000] = "M";

        return intToRomanUtil(num, dp);
    }

    int variance(unordered_map<char, int>& hashMap){ // O(26)
        int maxF = 0, minF = INT32_MAX;
        for(auto &p: hashMap){
            maxF = max(maxF, p.second);
            minF = min(minF, p.second);
        }
        return maxF - minF;
    }

    // leetcode hard. this will not work
    int largestVariance(string s) {
        int maxVar = 0;
        int l = 0, r = 0;
        unordered_map<char, int> hashMap;
        hashMap[s[l]]++;

        while(l <= r && r < s.size()){
            int curVar = variance(hashMap);

            if(l == r || curVar >= maxVar){
                maxVar = max(maxVar, curVar);
                r++;
                hashMap[s[r]]++;
            }
            else{
                hashMap[s[l]]--;
                if(hashMap[s[l]] == 0){
                    hashMap.erase(s[l]);
                }
                l++;
            }
        }

        return maxVar;
    }

    int sumSubarrayMins(vector<int>& arr) {
        int sum = 0;
        // sort kora jabe na to hold the arr property
        for(int i = 0; i < arr.size(); ++i){
            int minElem = arr[i];
            for(int j = i; j < arr.size(); ++j){
                minElem = min(minElem, arr[j]);
                sum += minElem;
            }
        }

        return sum;
    }

    int totalStrength(vector<int>& strength) {
        unsigned long long int globalSum = 0;
        int modder = pow(10, 9) + 7;

        for(int i = 0; i < strength.size(); ++i){
            int localSum = 0;
            int minima = strength[i];
            for(int j = i; j < strength.size(); ++j){
                localSum += strength[j];
                minima = min(minima, strength[j]);
                globalSum += localSum * minima;
            }
        }
        return globalSum % modder;
    }

//    static bool comparator(string& s1, string& s2){
//        int idx1 = s1.find(" "); // first occurance idx
//        int idx2 = s2.find(" ");
//
//        if(isalpha(s1[idx1 + 1]) && isalpha(s2[idx2 + 1])){ // letter logs
//            string s1Data = s1.substr(idx1 + 1);
//            string s2Data = s2.substr(idx2 + 1);
//            if(s1Data == s2Data){
//                return s1[idx1 - 1] < s2[idx2 - 1]; //identifier serial
//            }
//            return s1Data < s2Data; //increasing order
//        }
//        else if(isalpha(s1[idx1 + 1]) && !isalpha(s2[idx2 + 1])){
//            return true; // the order is correct
//        }
//        else if(!isalpha(s2[idx2 + 1]) && isalpha(s1[idx1 + 1])){
//            return false; // we want it to be flipped
//        }
//        else{
//            return false; // means we don't care about the situation and we want to keep the
//            // initial relative order intact
//        }
//    }

    struct {
        bool operator()(string& s1, string& s2){
            int idx1 = s1.find(" "); // first occurance idx
            int idx2 = s2.find(" ");

            if(isalpha(s1[idx1 + 1]) && isalpha(s2[idx2 + 1])){ // letter logs
                string s1Data = s1.substr(idx1 + 1);
                string s2Data = s2.substr(idx2 + 1);
                if(s1Data == s2Data){
                    return s1[idx1 - 1] < s2[idx2 - 1]; //identifier serial
                }
                return s1Data < s2Data; //increasing order
            }
            else if(isalpha(s1[idx1 + 1]) && !isalpha(s2[idx2 + 1])){
                return true; // the order is correct
            }
            else if(isalpha(s2[idx2 + 1]) && !isalpha(s1[idx1 + 1])){
                return false; // we want it to be flipped
            }
            else{
                return false; // means we don't care about the situation and we want to keep the
                // initial relative order intact
            }
        }
    } myCmp;

    vector<string> reorderLogFiles(vector<string>& logs) {
        sort(logs.begin(), logs.end(), myCmp);
        return logs;
    }

    // fractional knapsack
    // must be a static func or a struct
    static bool comparator(vector<int>& v1, vector<int>& v2){
        return v1[1] > v2[1]; // descending order
    }
    int maximumUnits(vector<vector<int>>& boxTypes, int truckSize) {
        int totalUnits = 0;
        sort(boxTypes.begin(), boxTypes.end(), comparator);

        for(int i = 0; truckSize > 0 && i < boxTypes.size(); ++i){
            int minima = min(boxTypes[i][0], truckSize);
            totalUnits += minima * boxTypes[i][1];
            truckSize -= minima;
        }
        return totalUnits;
    }

    Node* copyRandomList(Node* head) {
        unordered_map<Node*, Node*> hashMap;
        Node* newHead = new Node(head->val);
        hashMap[head] = newHead;

        Node* oldNode(head);
        Node* newNode(newHead);

        while(oldNode != nullptr){
            if(oldNode->next == nullptr){
                newNode->next = nullptr;
            }
            else{
                if(hashMap.count(oldNode->next)){
                    newNode->next = hashMap[oldNode->next];
                }
                else{
                    newNode->next = new Node(oldNode->next->val);
                    hashMap[oldNode->next] = newNode->next;
                }
            }

            if(oldNode->random == nullptr){
                newNode->random = nullptr;
            }
            else{
                if(hashMap.count(oldNode->random)){
                    newNode->random = hashMap[oldNode->random];
                }
                else{
                    newNode->random = new Node(oldNode->random->val);
                    hashMap[oldNode->random] = newNode->random;
                }
            }

            oldNode = oldNode->next;
            newNode = newNode->next;
        }

        return newHead;
    }

    void bfs(unordered_map<int, vector<int>>& adj, int cur, vector<int>& ans, int k){
        queue<int> q;
        unordered_set<int> visited;
        q.push(cur);

        while(!q.empty() && k + 1 != 0){
            int levelLen = q.size();
            for(int i = 0; i < levelLen; ++i){
                int cur = q.front();
                q.pop();
                visited.insert(cur);
                if(k == 0){
                    ans.push_back(cur);
                }
                for(auto &n: adj[cur]){
                    if(visited.count(n) == 0){
                        q.push(n);
                    }
                }
            }
            k--; // one level done
        }
    }

    void turnTreeToGraph(TreeNode* root, unordered_map<int, vector<int>>& adj){
        if(root == nullptr){
            return;
        }
        else if(root->left == nullptr && root->right == nullptr){
            return;
        }
        else if(root->left != nullptr && root->right == nullptr){
            adj[root->val].push_back(root->left->val);
            adj[root->left->val].push_back(root->val);
            turnTreeToGraph(root->left, adj);
        }
        else if(root->right != nullptr && root->left == nullptr){
            adj[root->val].push_back(root->right->val);
            adj[root->right->val].push_back(root->val);
            turnTreeToGraph(root->right, adj);
        }
        else {
            adj[root->val].push_back(root->left->val);
            adj[root->left->val].push_back(root->val);
            turnTreeToGraph(root->left, adj);

            adj[root->val].push_back(root->right->val);
            adj[root->right->val].push_back(root->val);
            turnTreeToGraph(root->right, adj);
        }
    }

    // Graph problem
    vector<int> distanceK(TreeNode* root, TreeNode* target, int k) {
        unordered_map<int, vector<int>> adj;
        vector<int> ans;
        turnTreeToGraph(root, adj);
        bfs(adj, target->val, ans, k);
        return ans;
    }

    int appealSum(string s) {
        unordered_set<char> hashSet;
        int count = 0;

        for(int i = 0; i < s.size(); ++i){
            hashSet.clear();
            for(int j = i; j < s.size(); ++j){
                if(hashSet.count(s[j]) == 0){
                    hashSet.insert(s[j]);
                }
                count += hashSet.size();
            }
        }
        return count;
    }

    int triangularSum(vector<int>& nums) {
        if(nums.size() == 1){
            return nums[0];
        }
        vector<int> tmp;
        for(int i = 0; i< nums.size() - 1; ++i){
            tmp.push_back((nums[i] + nums[i + 1]) % 10);
        }
        return triangularSum(tmp);
    }

    int minSwaps(string s) {
        string tmp1("1");
        string tmp2("0");

        for(int i = 0; i < s.size(); ++i){
            tmp1 += tmp1.back() == '1' ? '0' : '1';
            tmp2 += tmp2.back() == '1' ? '0' : '1';
        }

        int num1tmp1 = 0;
        int num0tmp1 = 0;
        int num1tmp2 = 0;
        int num0tmp2 = 0;
        for(int i = 0; i < s.size(); ++i){
            if(s[i] != tmp1[i]){
                if(s[i] == '0'){
                    num1tmp1++;
                }
                else{
                    num0tmp1++;
                }
            }
            if(s[i] != tmp2[i]){
                if(s[i] == '0'){
                    num1tmp2++;
                }
                else{
                    num0tmp2++;
                }
            }
        }

        if(num1tmp1 == num0tmp1 && num1tmp2 == num0tmp2){
            return min(num1tmp1, num1tmp2);
        }
        else if(num1tmp1 != num0tmp1 && num1tmp2 == num0tmp2){
            return num1tmp2;
        }
        else if(num1tmp1 == num0tmp1 && num1tmp2 != num0tmp2){
            return num1tmp1;
        }
        else{
            return -1;
        }
    }


    int numberOfWaysUtil(string& s, int i, string& tmp){
        if(tmp.size() == 3){
            if(tmp == "101" || tmp == "010"){
                return 1;
            }
            else{
                return 0;
            }
        }
        if(i == s.size()){
            return 0;
        }

        int totalWays = 0;
        tmp.push_back(s[i]);
        totalWays += numberOfWaysUtil(s, i + 1, tmp);
        tmp.pop_back();
        totalWays += numberOfWaysUtil(s, i + 1, tmp);
        return totalWays;
    }

    int numberOfWays(string s) {
        string tmp = "";
        vector<int> mem(s.size(), -1);
        return numberOfWaysUtil(s, 0, tmp);
    }
};

class KthLargest {
private:
    priority_queue<int, vector<int>, greater<int>> minHeap;
    int heapSize;
public:
    // constructor
    // O(k + (n - k)logk) = O(nlogk)
    KthLargest(int k, vector<int>& nums): heapSize(k) {
        for(int i = 0; i < nums.size(); ++i){
            add(nums[i]);
        }
    }

    // O(logn) for each query
    int add(int val) {
        if(minHeap.size() < heapSize){
            minHeap.push(val);
        }
        else{
            if(val > minHeap.top()){
                minHeap.pop();
                minHeap.push(val);
            }
        }

        return minHeap.top();
    }
};

class ParkingSystem {
private:
    int arr[4];
public:
    ParkingSystem(int big, int medium, int small) {
        arr[0] = big;
        arr[1] = medium;
        arr[2] = small;
    }

    bool addCar(int carType) {
        if(arr[carType] == 0){
            return false;
        }
        else{
            arr[carType]--;
            return true;
        }
    }
};

// This is a DFS solution
// But this problem can be more optimally solved with a BFS
class GetFood {
private:
    int rows;
    int cols;
    int getFoodUtil(vector<vector<char>>& grid, vector<vector<int>>& mem, int r, int c){
        if(grid[r][c] == '#'){
            return 0;
        }

        if(mem[r][c] != -2){
            return mem[r][c];
        }

        grid[r][c] = 'X'; // marked as visited basically

        int minDist = INT16_MAX;
        for(auto &d: DIR){
            int newR = r + d[0];
            int newC = c + d[1];

            if(newR >= 0 && newR < rows && newC >= 0 && newC < cols &&
               grid[newR][newC] != 'X'){
                int dist = getFoodUtil(grid, mem, newR, newC);
                if(dist != -1){
                    minDist = min(minDist, 1 + dist);
                }
            }
        }
        grid[r][c] = 'O';
        return mem[r][c] = minDist >= INT16_MAX ? -1 : minDist;
    }
public:
    int getFood(vector<vector<char>>& grid) {
        rows = grid.size();
        cols = grid[0].size();
        vector<vector<int>> mem(rows, vector<int>(cols, -2));

        int sr, sc;
        for(int r = 0; r < rows; ++r){
            for(int c = 0; c < cols; ++c){
                if(grid[r][c] == '*'){
                    sr = r;
                    sc = c;
                    break;
                }
            }
        }

        return getFoodUtil(grid, mem, sr, sc);
    }
};

class RottingOrange {
private:
    int rows;
    int cols;
public:
    int orangesRotting(vector<vector<int>>& grid) {
        int mins = 0;
        int freshOranges = 0;
        rows = grid.size();
        cols = grid[0].size();

        queue<vector<int>> q;
        for(int r = 0 ; r < rows; ++r){
            for(int c = 0; c < cols; ++c){
                if(grid[r][c] == 2){
                    q.push(vector<int>({r, c}));
                }
                else if(grid[r][c] == 1){
                    freshOranges++;
                }
            }
        }

        if(q.empty() && freshOranges == 0){
            return 0;
        }

        // multi source BFS
        while(!q.empty()){
            int levelLen = q.size();
            mins++;
            for(int i = 0; i < levelLen; ++i){
                int curRow = q.front()[0];
                int curCol = q.front()[1];
                q.pop();
                for(auto &d: DIR){
                    int nRow = curRow + d[0];
                    int nCol = curCol + d[1];
                    if(nRow >= 0 && nRow < rows && nCol >= 0 && nCol < cols &&
                        grid[nRow][nCol] == 1){
                        // has to be fresh
                        // make it rotten now
                        freshOranges--;
                        grid[nRow][nCol] = 2;
                        q.push(vector<int>({nRow, nCol}));
                    }
                }
            }
        }

        return freshOranges == 0 ? mins : -1;
    }

};

class ConcatanatedWords {
private:
    vector<string> ans;
    unordered_set<string> hashSet;
    void util(vector<string>& words, string& tmp, int& numWords){
        if(tmp.size() > 30){
            return;
        }
        if(hashSet.count(tmp) == 1 && numWords >= 2){
            ans.push_back(tmp);
        }

        // protibar we are starting from the beginning since there's no order during
        // concatanation and one word can be added multiple times
        for(int i = 0; i < words.size(); ++i){
            int idx = tmp.size();
            tmp += words[i];
            numWords++;
            util(words, tmp, numWords);
            tmp.erase(idx);
            numWords--;
        }
    }
public:
    vector<string> findAllConcatenatedWordsInADict(vector<string>& words)     {
        for(auto &s: words){
            hashSet.insert(s);
        }
        string tmp = "";
        int numWords = 0;
        util(words, tmp, numWords);
        return ans;
    }
};


#endif //BLIND75_AMAZONINTERVIEW_H
