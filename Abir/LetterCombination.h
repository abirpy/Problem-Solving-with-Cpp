//
// Created by Tahsinul Haque Abir on 7/8/22.
//

#ifndef ABIR_LETTERCOMBINATION_H
#define ABIR_LETTERCOMBINATION_H
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;


class LetterCombination {
public:
    void letterCombinationsUtil(string& digits, int idx, unordered_map<int, vector<char>>& map,
                                string& output, vector<string>& ans){
        if(output.size() == digits.size()){
            ans.push_back(output);
            return;
        }


        for(char& ch: map[digits[idx] - '0']){
            output.push_back(ch);
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
};


#endif //ABIR_LETTERCOMBINATION_H
