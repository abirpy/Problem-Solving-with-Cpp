//
// Created by Tahsinul Haque Abir on 7/8/22.
//

#ifndef ABIR_GENERATEPARENTHESIS_H
#define ABIR_GENERATEPARENTHESIS_H
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

// Time complexity: O(n * 2^2n) binary tree with height 2n
class GenerateParenthesis {
    void generateParenthesisUtil(int n, int numOpen, int numClosed, vector<string>& ans,
                                           string& output) {
        if(numOpen == numClosed && numOpen == n){
            ans.push_back(output);
            return;
        }
        if(numOpen < n){
            output.push_back('(');
            generateParenthesisUtil(n, numOpen + 1, numClosed, ans, output);
            output.pop_back();
        }
        if(numOpen > numClosed){
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
};


#endif //ABIR_GENERATEPARENTHESIS_H
