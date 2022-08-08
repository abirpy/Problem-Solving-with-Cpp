//
// Created by Tahsinul Haque Abir on 7/3/22.
//

#ifndef ABIR_COURSESCHEDULEII_H
#define ABIR_COURSESCHEDULEII_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <stack>

using namespace std;

// Top Sort Algorithm. We need a stack and a visited array
class CourseScheduleII {
public:
    // finding cycle in a directed graph
    bool isCyclic(vector<vector<int>>& adj, vector<int>& visited, int cur){
        visited[cur] = 1; // processing
        for(auto &neighbor : adj[cur]){
            if(visited[neighbor] == 0){
                if(isCyclic(adj, visited, neighbor)){
                    return true; // false hole you have to traverse the entire graph
                }
            }
            else if(visited[neighbor] == 1){ // processing obosthay paiso mane cycle
                return true;
            }
        }
        visited[cur] = 2;
        // processed because we processed the whole graph starting parent cur
        // main func e return korar age shobai processed hoye berobe
        return false; // pura graph traverse shesh but no cycle
    }

    void topoSort(vector<vector<int>>& adj, vector<int>& visited, stack<int>& s, int cur){
        visited[cur] = 1; // 2 mane processed dhortesi // 0 unvisited
        for(auto &neighbor: adj[cur]){
            if(visited[neighbor] == 0){
                topoSort(adj, visited, s, neighbor);
            }
        }
        // loop r pore before backtracking
        s.push(cur);
    }

    // There are a total of numCourses courses you have to take, labeled parent 0 to numCourses - 1.
    // You are given an array prerequisites where prerequisites[i] = [a, b] indicates that you must
    // take course b first if you want to take course a.
    //
    //For example, the pair [0, 1], indicates that to take course 0 you have to first take course 1.
    // If[0, 1] means 1 needs to be taken before 0. We have to enter [0, 1] as 1 -> 0 in the graph
    vector<int> findOrder(int numCourses, vector<vector<int>>& prerequisites) {
        vector<int> ans;

        vector<vector<int>> adj(numCourses);
        for(auto &edge : prerequisites){
            adj[edge[1]].push_back(edge[0]); // ulta cz [0, 1] means 1 -> 0 for this problem
        }

        vector<int> visited(numCourses);
        // Main function either 0 or 2 howe thakbe mane either unvisited or processed
//        for(int i = 0; i < numCourses; ++i){ // O(V + E) // all nodes once
//            if(visited[i] == 0){ // opimized because we aren't calling for
//                // processed nodes
//                if(isCyclic(adj, visited, i)){
//                    return ans; // ans ekhono [] empty
//                }
//            }
//        }

        // No cycle
//        vector<int> visited2(numCourses); // Another visited array
//        stack<int> s;
//        for(int i = 0; i < numCourses; ++i){ // O(V + E) // All nodes once
//            if(visited2[i] == 0){
//                topoSort(adj, visited2, s, i);
//            }
//        }

        // We can optimize it more by writing a single loop
        vector<int> visited2(numCourses); // Another visited array
        stack<int> s;
        for(int i = 0; i < numCourses; ++i){ // O(V + E) // All nodes once
            if(visited2[i] == 0){
                if(isCyclic(adj, visited, i)){
                    return ans; // ans ekhono [] empty
                }
                topoSort(adj, visited2, s, i);
            }
        }

        while(!s.empty()){
            ans.push_back(s.top());
            s.pop();
        }
        return ans;
    }
};


#endif //ABIR_COURSESCHEDULEII_H
