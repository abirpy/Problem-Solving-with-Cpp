//
// Created by Tahsinul Haque Abir on 7/5/22.
//

#ifndef ABIR_BACKTRACKING_H
#define ABIR_BACKTRACKING_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <algorithm>

using namespace std;

const vector<vector<int>> DIR{
        {0, 1}, //north
        {0, -1}, //s
        {1, 0}, //e
        {-1, 0} //w
};

class Backtracking {
private:
    int row;
    int col;
public:
    bool isValid(int r, int c, vector<vector<char>>& board){
        if(r < 0 || c < 0 || r >= row || c >= col || board[r][c] == '#'){
            return false;
        }
        return true;
    }

    bool isPrefix(string& word, string& target, int idx){ // Optimizing. Instead of iterating
        // over the whole string we are just checking the last char // O(1)
        if(word[idx] == target.back()){
            return true;
        }
        return false;
    }

    bool existUtil(vector<vector<char>>& board, string& word, string& target, int r, int c, int
    idx){
        target.push_back(board[r][c]);
        board[r][c] = '#';

        if(target == word){
            return true;
        }
        else if(isPrefix(word, target, idx)){ // if target is not equal to word it might be a prefix
            for(auto &dirPair: DIR){
                if(isValid(r + dirPair[0], c + dirPair[1], board)){
                    if(existUtil(board, word, target, r + dirPair[0], c + dirPair[1], idx + 1)){
                        return true;
                    }
                }
            }
        }

        // If none of the conditions are true then backtrack by undoing everything
        board[r][c] = target.back();
        target.pop_back();
        return false;
    }

    bool exist(vector<vector<char>>& board, string word){
        row = board.size();
        col = board[0].size();
//        vector<vector<int>> visited(row, vector<int>(col, 0));
        // Instead of using visited we can use the board itself as a flag
        string target = "";

        for(int r = 0; r < row; ++r){
            for(int c = 0; c < col; ++c){
                if(existUtil(board, word, target, r, c, -1)){
                    return true;
                }
            }
        }
        return false;
    }
};


#endif //ABIR_BACKTRACKING_H
