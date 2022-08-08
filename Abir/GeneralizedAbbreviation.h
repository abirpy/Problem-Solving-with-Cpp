//
// Created by Tahsinul Haque Abir on 7/8/22.
//

#ifndef ABIR_GENERALIZEDABBREVIATION_H
#define ABIR_GENERALIZEDABBREVIATION_H
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

class GeneralizedAbbreviation {
public:
    void util(string& word, int idx, string abb, int count, vector<string>& ans){
        if(idx == word.size()){
            if(count != 0){
                char c = '0' + count;
                ans.push_back(abb + c);
            }
            else{
                ans.push_back(abb);
            }
            return;
        }

        // Now take word[idx] into abb
        if(count != 0){
            char c = '0' + count;
            util(word, idx + 1, abb + c + word[idx], 0, ans);
        }
        else{
            util(word, idx + 1, abb + word[idx], count, ans);
        }

        util(word, idx + 1, abb, count + 1, ans); // now explore the possible options
    }

    vector<string> generateAbbreviations(string word) {
        vector<string> ans;
        string abb("");
        util(word, 0, abb, 0, ans);
        return ans;
    }
};


#endif //ABIR_GENERALIZEDABBREVIATION_H
