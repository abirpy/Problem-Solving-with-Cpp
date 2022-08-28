//
// Created by Tahsinul Haque Abir on 7/1/22.
//

#ifndef ABIR_COMPONENTS_H
#define ABIR_COMPONENTS_H

#include <iostream>
#include <vector>
#include <unordered_set>

using namespace std;

class Components {
    void dfs(vector<vector<int>> adj, unordered_set<int> visited, int cur){
        visited.insert(cur); //cur node is visited

        for(auto &neighbor : adj[cur]){
            if(visited.count(neighbor) == 0){
                dfs(adj, visited, neighbor);
            }
        }
    }
    int countComponents(int n, vector<vector<int>>& edges) {
        // Create an adjacency list
        vector<vector<int>> adj(n);
        for(auto &edge : edges){
            adj[edge[0]].push_back(edge[1]);
            adj[edge[1]].push_back(edge[0]);
        }
        unordered_set<int> visited; // visited set
//        int cur;
//        for(int i = 0; i < n; ++i){
//            if(!adj[i].empty()){
//                cur = i;
//                break;
//            }
//        }

        int count = 0;
        for(int i = 0; i < n; ++i){
            if(visited.count(i) == 0){
                count++;
                dfs(adj, visited, i);
            }
        }

        return count;
    }
};


#endif //ABIR_COMPONENTS_H
