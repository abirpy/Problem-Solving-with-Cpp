//
// Created by Tahsinul Haque Abir on 7/3/22.
//

#ifndef ABIR_MHT_H
#define ABIR_MHT_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <queue>

using namespace std;

// BFS problem
class MHT {
public:
    vector<int> findMinHeightTrees(int n, vector<vector<int>>& edges) {
        if(n  == 1){
            return vector<int>({0});
        }
        if(n == 2){
            return vector<int>({0, 1});
        }
        // undirected graph
        vector<vector<int>> adj(n);
        vector<int> degreeArr(n, 0);
        for(auto &edge: edges){
            adj[edge[0]].push_back(edge[1]);
            adj[edge[1]].push_back(edge[0]);
            degreeArr[edge[0]]++;
            degreeArr[edge[1]]++;
        }

        queue<int> q;
        // Fill the queue for the first time with degreeArr
        for(int i = 0; i < n; ++i){
            if(degreeArr[i] == 1){
                q.push(i);
            }
        }
        int count = n;
        vector<int> ans;
        while(count > 2){
            // if the queue becomes empty and count > 2 then hocce we have to
            // push things from ans to queue
            if(q.empty()){
                while(!ans.empty()){
                    q.push(ans.back());
                    ans.pop_back();
                }
            }

            int cur = q.front();
            q.pop(); // removing cur

            count--;
            degreeArr[cur]--; // decreasing degree of cur by 1
            for(auto &neighbor: adj[cur]){ // ei loop ta ek bar choleb bcz
                // degree of cur = 1
                degreeArr[neighbor]--;
                if(degreeArr[neighbor] == 1){
                    ans.push_back(neighbor);
                }
            }
        }
        return ans;
    }
};


#endif //ABIR_MHT_H
