//
// Created by Tahsinul Haque Abir on 7/1/22.
//

#ifndef ABIR_VALIDTREE_H
#define ABIR_VALIDTREE_H

#include <iostream>
#include <vector>
#include <unordered_set>

using namespace std;

class ValidTree {
private:
    int numNodes;
public:
    bool isCyclic(vector<vector<int>>& adj, unordered_set<int>& visited, int prev, int cur){
        visited.insert(cur); //cur node is visited

        for(auto &neighbor : adj[cur]){
            if(visited.count(neighbor) == 1 && neighbor != prev){
                return true;
            }
            else if(visited.count(neighbor) != 1){
                 if(isCyclic(adj, visited, cur, neighbor)){ // if cyclic return kore felo
                     // but not cyclic hole don't return pura graph traverse kore lagbe
                     return true;
                 }
            }
        }
        // pura graph traverse shesh but no cycle
        return false;
    }

    bool validTree(int n, vector<vector<int>>& edges) {
        if(edges.size() == 0 && n == 1){
            return true;
        }

        if(edges.size() == 0 && n > 1){
            return true;
        }

        // it has n nodes. So for it to be a tree it must have n - 1 edges
        if(edges.size() != n - 1){
            return false;
        }

        // Create an adjacency list
        vector<vector<int>> adj(n);
        numNodes = n;
        for(auto &edge : edges){
            adj[edge[0]].push_back(edge[1]);
            adj[edge[1]].push_back(edge[0]);
        }
        unordered_set<int> visited; // visited set
        int cur;
        for(int i = 0; i < numNodes; ++i){
            if(!adj[i].empty()){
                cur = i;
                break;
            }
        }

        if(isCyclic(adj, visited, -1, cur)){ // cyclic
            return false;
        }
        else{ // not cyclic. So check the next condition that the graph is connected or not
            return visited.size() == numNodes;
        }
    }
};


#endif //ABIR_VALIDTREE_H
