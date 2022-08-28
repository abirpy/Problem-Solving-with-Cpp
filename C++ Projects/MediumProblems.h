//
// Created by Tahsinul Haque Abir on 8/26/22.
//

#ifndef ABIR_MEDIUMPROBLEMS_H
#define ABIR_MEDIUMPROBLEMS_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;

const vector<vector<int>> DIR({
      {1, 0},
      {-1, 0},
      {0, 1},
      {0, -1}
});

struct ListNode {
    int val;
    ListNode *next;
    ListNode() : val(0), next(nullptr) {}
    ListNode(int x) : val(x), next(nullptr) {}
    ListNode(int x, ListNode *next) : val(x), next(next) {}
};

class MediumProblems{
public:
    void permuteUtil(vector<int>& nums, vector<int>& output, vector<vector<int>>& ans, vector<int>&
            visited){
        if(output.size() == nums.size()){
            ans.push_back(output);
            return;
        }

        for(int i = 0; i < nums.size(); ++i){
            if(visited[i] == 0){
                visited[i] = 1;
                output.push_back(nums[i]);
                permuteUtil(nums, output, ans, visited);
                output.pop_back();
                visited[i] = 0;
            }
        }
    }
    vector<vector<int>> permute(vector<int>& nums) {
        vector<int> visited(nums.size(), 0);
        vector<vector<int>> ans;
        vector<int> output;
        permuteUtil(nums, output, ans, visited);
        return ans;
    }

    vector<int> twoSum(vector<int>& numbers, int target) {
        vector<int> ans;
        int l = 0, r = numbers.size() - 1;
        while(l < r){
            int twoSum = numbers[l] + numbers[r];
            if(twoSum == target){
                ans.push_back(numbers[l]);
                ans.push_back(numbers[r]);
                l++;
                r--;
            }
            else if(twoSum < target){
                l++;
            }
            else{
                r--;
            }
        }
        return ans;
    }

    int lengthOfLIS(vector<int>& nums) {
        vector<int> LIS(nums.size(), 1);
        int lenLIS = 1;

        for(int i = 0; i < nums.size(); ++i){
            for(int j = 0; j < i; ++j){
                if(nums[i] > nums[j] && LIS[j] + 1 > LIS[i]){
                    LIS[i] = LIS[j] + 1;
                }
            }
            lenLIS = max(lenLIS, LIS[i]);
        }
        return lenLIS;
    }

    void transpose(vector<vector<int>>& matrix){
        for(int i = 0; i < matrix.size(); ++i){
            for(int j = i + 1; j < matrix[0].size(); ++j){
                swap(matrix[i][j], matrix[j][i]);
            }
        }
    }

    void reverseRows(vector<vector<int>>& matrix){
        int col = matrix[0].size();
        for(int i = 0; i < matrix.size(); ++i){
            for(int j = 0; j < col / 2; ++j){
                swap(matrix[i][j], matrix[i][col -j - 1]);
            }
        }
    }

    void rotate(vector<vector<int>>& matrix) {
        transpose(matrix);
        reverseRows(matrix);
    }

    int partition(vector<int>& nums, int l, int r){
        // l is the pivot idx
        int pivot = nums[l];
        int sl = l, sr = r;

        while(sl < sr){
            while(nums[sr] > pivot){
                sr--;
            }
            while(sl < sr && nums[sl] <= pivot){
                sl++;
            }
            if(sl < sr){
                swap(nums[sl], nums[sr]);
            }
        }
        // now sr is at the correct idx of the pivot
        swap(nums[l], nums[sr]);
        return sr;
    }

    int findKthLargest(vector<int>& nums, int k) {
        // find kth largest means (n - k + 1)th smallest elem which is at the index n - k
        int n = nums.size();
        int l = 0, r = n - 1;
        int targetIdx = n - k;
        int pivotIdx = partition(nums, l, r);

        while(pivotIdx != targetIdx){
            if(pivotIdx < targetIdx){
                pivotIdx = partition(nums, pivotIdx + 1, r);
            }
            else{
                pivotIdx = partition(nums, l, pivotIdx - 1);
            }
        }
        return nums[targetIdx];
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

    void setZeroes(vector<vector<int>>& matrix) {
        int rows = matrix.size();
        int cols = matrix[0].size();
        vector<int> rowMarker(rows, 1);
        vector<int> colMarker(cols, 1);

        for(int r = 0; r < rows; ++r){
            for(int c = 0; c < cols; ++c){
                if(matrix[r][c] == 0){
                    rowMarker[r] = 0;
                    colMarker[c] = 0;
                }
            }
        }

        for(int r = 0; r < rows; ++r){
            if(rowMarker[r] == 0){
                for(int c = 0; c < cols; ++c){
                    matrix[r][c] = 0;
                }
            }
        }

        for(int c = 0; c < cols; ++c){
            if(colMarker[c] == 0){
                for(int r = 0; r < rows; ++r){
                    matrix[r][c] = 0;
                }
            }
        }
    }

    bool canJumpUtil(vector<int>& nums, int i, unordered_map<int, bool>& mem){
        if(i == 0){
            return true;
        }
        if(mem.count(i) == 1){
            return mem[i];
        }
        for(int k = i - 1; k >= 0; --k){
            if(k + nums[k] >= i && canJumpUtil(nums, k, mem)){
                return mem[i] = true;
            }
        }
        return mem[i] = false;
    }

    bool canJump(vector<int>& nums){
        unordered_map<int, bool> mem;
        return canJumpUtil(nums, nums.size() - 1, mem);
    }

    vector<int> spiralOrder(vector<vector<int>>& matrix) {
        vector<int> output;
        int rows = matrix.size();
        int cols = matrix[0].size();
        int u = 0, d = rows - 1, l = 0, r = cols - 1;
        int dir = 0;
        while(u <= d && l <= r){
            if(dir == 0){
                // go right
                for(int i = l; i <= r; ++i){
                    output.push_back(matrix[u][i]);
                }
                u++;
            }
            else if(dir == 1){
                for(int i = u; i <= d; ++i){
                    output.push_back(matrix[i][r]);
                }
                r--;
            }
            else if(dir == 2){
                for(int i = r; i >= l; --i){
                    output.push_back(matrix[d][i]);
                }
                d--;
            }
            else{
                for(int i = d; i >= u; --i){
                    output.push_back(matrix[i][l]);
                }
                l++;
            }
            dir = (dir + 1) % 4;
        }
        return output;
    }

    int findDuplicate(vector<int>& nums){
        int s = 0;
        int f = 0;

        while(s != f){
            s = nums[s];
            f = nums[nums[f]];
        }
        // now s = f
        // find the entrance now
        f = nums[0];
        while(s != f){
            s = nums[s];
            f = nums[f];
        }
        return s;
    }

    bool isRobotBounded(string instructions) {
        int dirX = 0, dirY = 1; // north (0, 1)
        int x = 0, y = 0;

        for(int i = 1; i <= 4; ++i){
            for(int k = 0; k < instructions.size(); ++k){
                if(instructions[k] == 'G'){
                    x = x + dirX;
                    y = y + dirY;
                }
                else if(instructions[k] == 'L'){
                    // rotate DIR
                    swap(dirX, dirY);
                    dirX = -dirX;
                }
                else{
                    swap(dirX, dirY);
                    dirY = -dirY;
                }
            }
        }
        return (x == 0 && y == 0) && (dirX == 0 && dirY == 1);
    }
};

class LL {
    void insertLast(ListNode* &head, ListNode* &tail, ListNode* &node){
        if(head == tail && head == nullptr){
            head = tail = node;
        }
        else{
            tail->next = node;
            tail = node;
        }
    }
    int lenList(ListNode* head){
        int count = 0;
        for(ListNode* tmp = head; tmp != nullptr; tmp = tmp->next){
            count++;
        }
        return count;
    }
public:
    ListNode* removeNthFromEnd(ListNode* head, int n) {
        ListNode* newHead(nullptr);
        ListNode* newTail(nullptr);

        int len = lenList(head);

        int count = 0;
        while(head != nullptr){
            count++;
            ListNode* tmp = head;
            head = head->next;
            tmp->next = nullptr;
            if(count != len - n + 1){
                insertLast(newHead, newTail, tmp);
            }
        }
        return newHead;
    }
};

class PacificAtlanticWaterFlow {
private:
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
    unordered_set<vector<int>, hashFunction> setP;
    unordered_set<vector<int>, hashFunction> setA;
    int rows;
    int cols;
    void dfs(vector<vector<int>>& heights, unordered_set<vector<int>, hashFunction>& visited, int
    r, int c, int prevHeight){
        if(r < 0 || r >= rows || c < 0 || c >= cols || visited.count(vector<int>({r, c})) == 1 ||
        heights[r][c] < prevHeight){
            return;
        }
        visited.insert(vector<int>({r, c}));
        for(auto &d: DIR){
            int newR = r + d[0];
            int newC = c + d[1];
            dfs(heights, visited, newR, newC, heights[r][c]);
        }
    }
public:
    vector<vector<int>> pacificAtlantic(vector<vector<int>>& heights) {
        vector<vector<int>> ans;
        rows = heights.size();
        cols = heights[0].size();

        for(int c = 0; c < cols; ++c){
            if(setP.count(vector<int>({0, c})) == 0){
                dfs(heights, setP, 0, c, 0);
            }
            if(setA.count(vector<int>({rows - 1, c})) == 0){
                dfs(heights, setA, rows - 1, c, 0);
            }
        }

        for(int r = 0; r < rows; ++r){
            if(setP.count(vector<int>({r, 0})) == 0){
                dfs(heights, setP, r, 0, 0);
            }
            if(setA.count(vector<int>({r, cols - 1})) == 0){
                dfs(heights, setA, r, cols - 1, 0);
            }
        }

        for(auto &v: setP){
            if(setA.count(v) == 1){
                ans.push_back(v);
            }
        }
    }
};


#endif //ABIR_MEDIUMPROBLEMS_H
