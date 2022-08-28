//
// Created by Tahsinul Haque Abir on 7/2/22.
//

#ifndef ABIR_COURSESCHEHDULE_H
#define ABIR_COURSESCHEHDULE_H

#include <iostream>
#include <vector>
#include <unordered_set>

using namespace std;

class CourseSchehdule {
public:
    // isCyclic for directed graph
    bool isCyclic(vector<vector<int>>& adj, vector<int>& visited, int cur){
//        visited.insert(cur);
//        for(auto &neighbor : adj[cur]){
//            if(visited.count(neighbor) == 0){
//                if(isCyclic(adj, visited, neighbor)){
//                    return true; // false hole you have to traverse the entire graph
//                }
//            }
//            else{
//                return true;
//            }
//        }
//        visited.erase(cur);
//        return false; // pura graph traverse shesh but no cycle

        // Time complexity O(|E| + |V|^2) because for each V we are calling the isCyclic function
        // unlike the case of undirected graph where the visited array remains untouched
        // throughout the whole process and we only call unconnected components. That's that
        // was O(V+ E). For deatiled time compleixity check leetcode

        // Optimized approach uses unvisited - 0, processing - 1, processed - 2

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

    bool canFinish(int numCourses, vector<vector<int>>& prerequisites) {
//        // make an edge list
//        // Create an adjacency list
//        vector<vector<int>> adj(numCourses);
//        for(auto &edge : prerequisites){
//            adj[edge[0]].push_back(edge[1]);
//        }
//        unordered_set<int> visited; // visited set
//
//        for(int i = 0; i < numCourses; ++i){
//            if(visited.count(i) == 0){
//                if(isCyclic(adj, visited, i)){
//                    return false;
//                }
//            }
//        }
//        return true; // ekhono false return hoyni mane canFinish no cycles

        // Optimized approach uses unvisited - 0, processing - 1, processed - 2
        //We have to use an array as hashmap or a hashmap for this approach
        // make an edge list
        // Create an adjacency list
        vector<vector<int>> adj(numCourses);
        for(auto &edge : prerequisites){
            adj[edge[0]].push_back(edge[1]);
        }

        vector<int> visited(numCourses);

        for(int i = 0; i < numCourses; ++i){
            if(visited[i] == 0){ // opimized because we aren't calling for
                // processed nodes
                if(isCyclic(adj, visited, i)){
                    return false;
                }
            }
        }
        return true; // ekhono false return hoyni mane canFinish no cycles
    }
};


#endif //ABIR_COURSESCHEHDULE_H
