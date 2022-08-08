//
// Created by Tahsinul Haque Abir on 7/1/22.
//

#ifndef ABIR_ISLANDS_H
#define ABIR_ISLANDS_H

#include <iostream>
#include <vector>

using namespace std;

const vector<vector<int>> DIR({
      {0, 1},
      {0, -1},
      {1, 0},
      {-1, 0}
});

class Islands {
private:
    int rows;
    int cols;
public:
    bool isValid(vector<vector<char>>& grid, int r, int c){
        if(r < 0 || r >= rows || c < 0 || c >= cols || grid[r][c] == '0'){ // '0' mane water
            return false;
        }
        return true;
    }
    void dfs(vector<vector<char>>& grid, int r, int c){
        grid[r][c] = '2'; // mark as visited 2 mane visited

        for(auto &dirPair : DIR){
            if(isValid(grid, r + dirPair[0], c + dirPair[1]) && grid[r + dirPair[0]][c + dirPair[1]] == '1'){
                dfs(grid, r + dirPair[0], c + dirPair[1]);
            }
        }
    }

    int numIslands(vector<vector<char>>& grid) {
        rows = grid.size();
        cols = grid[0].size();
        int count = 0;

        for(int r = 0; r < rows; ++r){
            for(int c = 0; c < cols; ++c){
                if(grid[r][c] == '1'){ // Unvisited
                    count++;
                    dfs(grid, r, c);
                }
            }
        }
        return count;
    }
};


#endif //ABIR_ISLANDS_H
