//
// Created by Tahsinul Haque Abir on 8/7/22.
//

#ifndef ABIR_AMAZONPROBLEMS_H
#define ABIR_AMAZONPROBLEMS_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;
struct ListNode {
    int val;
    ListNode *next;
    ListNode(int x) : val(x), next(NULL) {}
};

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;

    TreeNode() : val(0), left(nullptr), right(nullptr) {}

    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}

    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
};

class AmazonProblems {
public:
    void letterCombinationsUtil(string& digits, int idx, unordered_map<int, vector<char>>& map,
                                string& output, vector<string>& ans){
        if(idx == digits.size()){
            ans.push_back(output);
            return;
        }

        for(auto &c: map[digits[idx] - '0']){
            output.push_back(c);
            letterCombinationsUtil(digits, idx + 1, map, output, ans);
            output.pop_back();
        }
    }

    vector<string> letterCombinations(string digits) {
        vector<string> ans;
        if(digits.size() == 0){
            return ans;
        }
        unordered_map<int, vector<char>> map;
        for(int i = 2; i <= 6; ++i){
            int j = (i - 2) * 3;
            int counter = 0;
            while(counter != 3){
                map[i].push_back('a' + j);
                j++;
                counter++;
            }
        }
        map[7] = vector<char>({'p','q', 'r', 's'});
        map[8] = vector<char>({'t', 'u', 'v'});
        map[9] = vector<char>({'w', 'x', 'y', 'z'});

        string output("");
        letterCombinationsUtil(digits, 0, map, output, ans);
        return ans;
    }

    // valid parenthesis
    bool isValid(string s) {
        stack<int> stk;

        for(char &c: s){
            if(c == '(' || c == '[' || c == '{'){
                stk.push(c);
            }
            else if(c == ')' || c == ']' || c == '}'){
                if(stk.empty()){
                    return false;
                }
                else{
                    if(c == ')' && stk.top() == '('){
                        stk.pop();
                        continue;
                    }
                    if(c == ')' && stk.top() != '('){
                        return false;
                    }
                    if(c == ']' && stk.top() == '['){
                        stk.pop();
                        continue;
                    }
                    if(c == ']' && stk.top() != '['){
                        return false;
                    }
                    if(c == '}' && stk.top() == '{'){
                        stk.pop();
                        continue;
                    }
                    if(c == '}' && stk.top() != '{'){
                        return false;
                    }
                }
            }
        }

        if(!stk.empty()){
            return false;
        }
        return true;
    }

    void generateParenthesisUtil(int n, int numOpen, int numClosed, vector<string>& ans, string&
    output) {
        if(numOpen == n && numOpen == numClosed){
            ans.push_back(output);
            return;
        }

        if(numOpen < n){
            output.push_back('(');
            generateParenthesisUtil(n, numOpen + 1, numClosed, ans, output);
            output.pop_back();
        }

        if(numOpen > numClosed && numClosed < n){
            output.push_back(')');
            generateParenthesisUtil(n, numOpen, numClosed + 1, ans, output);
            output.pop_back();
        }
    }
    vector<string> generateParenthesis(int n) {
        vector<string> ans;
        string output = "(";
        generateParenthesisUtil(n, 1, 0, ans, output);
        return ans;
    }

    void insertLast(ListNode* &head, ListNode* &tail, ListNode* &node){
        if(head == tail && head == nullptr){
            head = tail = node;
        }
        else{
            tail->next = node;
            tail = tail->next;
        }
        tail->next = nullptr;
    }

    ListNode* merge(ListNode* l1, ListNode* l2){
        if(l1 == nullptr){
            return l2;
        }
        else if(l2 == nullptr){
            return l1;
        }

        ListNode *nh, *nt;
        nh = nt = nullptr;
        while(l1 != nullptr && l2 != nullptr){
            if(l1->val < l2->val){
                ListNode* tmp1(l1);
                l1 = l1->next;
                insertLast(nh, nt, tmp1);
            }
            else if(l1->val > l2->val){
                ListNode* tmp2(l2);
                l2 = l2->next;
                insertLast(nh, nt, tmp2);
            }
            else{
                ListNode* tmp1(l1);
                l1 = l1->next;
                insertLast(nh, nt, tmp1);
                ListNode* tmp2(l2);
                l2 = l2->next;
                insertLast(nh, nt, tmp2);
            }
        }

        while(l1 != nullptr){
            ListNode* tmp1(l1);
            l1 = l1->next;
            insertLast(nh, nt, tmp1);
        }

        while(l2 != nullptr){
            ListNode* tmp2(l2);
            l2 = l2->next;
            insertLast(nh, nt, tmp2);
        }
        return nh;
    }

    ListNode* mergeKLists(vector<ListNode*>& lists) {
        if(lists.size() == 0){
            return nullptr;
        }
        else if(lists.size() == 1){
            return lists[0];
        }

        while(lists.size() != 1){
            vector<ListNode*> tmp;
            for(int i = 0; i < lists.size(); i = i + 2){
                ListNode* f = lists[i];
                ListNode* s = i == lists.size() - 1 ? nullptr : lists[i + 1];
                tmp.push_back(merge(f, s));
            }
            lists = tmp;
        }
        return lists[0];
    }

    int uniquePathsWithObstacles(vector<vector<int>>& obstacleGrid)     {
        int r = obstacleGrid.size();
        int c = obstacleGrid[0].size();
        vector<vector<int>> numWays(r, vector<int>(c));

        for(int i = 0; i < c; ++i){
            if(obstacleGrid[0][i] == 1){
                while(i < c){
                    numWays[0][i] = 0;
                    i++;
                }
                break;
            }
            numWays[0][i] = 1;
        }

        for(int i = 0; i < r; ++i){
            if(obstacleGrid[i][0] == 1){
                while(i < r){
                    numWays[i][0] = 0;
                    i++;
                }
                break;
            }
            numWays[i][0] = 1;
        }

        for(int i = 1; i < r; ++i){
            for(int j = 1; j < c; ++j){
                if(obstacleGrid[i][j] == 1){
                    numWays[i][j] = 0;
                }
                else{
                    numWays[i][j] = numWays[i - 1][j] + numWays[i][j - 1];
                }
            }
        }

        return numWays[r - 1][c - 1];
    }

    vector<vector<int>> levelOrder(TreeNode* root) {
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
                TreeNode* curNode = q.front();
                q.pop();
                tmp.push_back(curNode->val);
                if(curNode->left != nullptr){
                    q.push(curNode->left);
                }
                if(curNode->right != nullptr){
                    q.push(curNode->right);
                }
            }
            ans.push_back(tmp);
        }
    }

    vector<int> rightSideView(TreeNode* root) {
        vector<int> ans;
        if(root == nullptr){
            return ans;
        }
        queue<TreeNode*> q;
        q.push(root);

        while(!q.empty()){
            int levelLen = q.size();
            int rightMost;
            for(int i = 0; i < levelLen; ++i){
                TreeNode* curNode = q.front();
                q.pop();
                rightMost = curNode->val; // ei loop theke ber howar shomoy rightMost will have
                // the last value of the level
                if(curNode->left != nullptr){
                    q.push(curNode->left);
                }
                if(curNode->right != nullptr){
                    q.push(curNode->right);
                }
            }
            ans.push_back(rightMost);
        }
        return ans;
    }

    // nth row has n + 1 values
    vector<int> getRow(int rowIndex) {
        if(rowIndex == 0){
            return vector<int>({1});
        }
        else if(rowIndex == 1){
            return vector<int>({1, 1});
        }

        vector<int> ans(rowIndex + 1);
        vector<int> tmp = getRow(rowIndex - 1);
        ans[0] = 1;
        ans[rowIndex] = 1;
        for(int i = 0; i < tmp.size() - 1; ++i){
            ans[i + 1] = tmp[i] + tmp[i + 1];
        }
        return ans;
    }

    vector<int> twoSum(vector<int>& numbers, int target) {
        vector<int> res(2, -1);
        int l = 0, r = numbers.size() - 1;

        while(l < r){
            int sum = numbers[l] + numbers[r];
            if(sum == target){
                res[0] = l + 1;
                res[1] = r + 1;
                return res;
            }
            else if(sum < target){ // underestimate kore felso
                l++;
            }
            else{ // overestimate
                r--;
            }
        }
        return res;
    }

    // time: O(logn)
    // space: O(n)
    // a better solu would use no extra space by using floyd's cycle detection algo
    bool isHappy(int n) {
        unordered_set<int> hash;

        while(n != 1){
            int sum = 0;
            while(n > 0) {
                sum += pow(n % 10, 2);
                n = n / 10;
            }
            if(hash.count(sum) == 1){ // already found so cycle
                return false;
            }
            else{
                hash.insert(sum);
                n = sum;
            }
        }
        return true; // false return hoyni mane true
    }

    // sieve of eratosthenes
    int countPrimes(int n) {
        vector<bool> prime(n + 1, true);
        int count = 0;

        for(int i = 2; i <= sqrt(n); ++i){
            if(prime[i]){ // true
                for(int j = i*i; j <= n; j = j + i){ // mark all multiple of i as false
                    prime[j] = false;
                }
            }
        }

        for(int i = 2; i < n; i++){
            if(prime[i]){
                count++;
            }
        }

        return count;
    }

    // neg marking
    vector<int> findErrorNums(vector<int>& nums) {
        vector<int> ans;
        for(int i = 0; i < nums.size(); ++i){
            if(nums[abs(nums[i]) - 1] < 0){ // already less than zero mane duplicate ashce
                ans.push_back(abs(nums[i]));
            }
            else{
                nums[abs(nums[i]) - 1] = -nums[abs(nums[i]) - 1];
            }
        }

        for(int i = 0; i < nums.size(); ++i){
            if(nums[i] > 0){
                ans.push_back(i + 1);
            }
        }
        return ans;
    }
};


#endif //ABIR_AMAZONPROBLEMS_H
