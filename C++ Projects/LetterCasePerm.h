//
// Created by Tahsinul Haque Abir on 7/5/22.
//

#ifndef ABIR_LETTERCASEPERM_H
#define ABIR_LETTERCASEPERM_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

class LetterCasePerm {
public:
    void permUtil(string& s, string& target, vector<string>& ans, int idx){
        if(target.length() == s.length()){
            ans.push_back(target);
            return;
        }

        if(isalpha(s[idx])){
            vector<int> tmp({tolower(s[idx]), toupper(s[idx])});
            for(int k = 0; k < 2; k++){
                target.push_back(tmp[k]);
                permUtil(s, target, ans, idx + 1);
                target.pop_back();
            }
        }
        else{
            target.push_back(s[idx]);
            permUtil(s, target, ans, idx + 1);
            target.pop_back();
        }
    }
    vector<string> letterCasePermutation(string s) {
        vector<string> ans;
        string target = "";
        permUtil(s, target, ans, 0);
        return ans;
    }
};


#endif //ABIR_LETTERCASEPERM_H
