//
// Created by Tahsinul Haque Abir on 7/28/22.
//

#ifndef ABIR_GRAPHPROBLEMS_H
#define ABIR_GRAPHPROBLEMS_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;

class GraphProblems {
public:
    // visited can have three states 0, 1, 2
    // main function e return r age any node will be either unprocessed 0 or processed 2
    bool isCyclicDirected(vector<vector<int>>& adj, vector<int>& visited, int cur){
        visited[cur] = 1; // processsing
        for(auto &neighbor : adj[cur]){
            if(visited[neighbor] == 1){
                return true;
            } // processing stage e paisi mane loop ase
            else if(visited[neighbor] == 0){ // unprocessed
                if(isCyclicDirected(adj, visited, neighbor)){ // false hole you have to traverse
                    // the entire graph
                    return true;
                }
            }
        }
        visited[cur] = 2; // cur already processed now
        return false; // pura graph traverse shesh but true return hoyni
    }

    // visited has 2 states 0 and 1
    bool isCyclicUndirected(vector<vector<int>>& adj, vector<int>& visited, int prev, int
    cur){
        visited[cur] = 1;
        for(auto &neighbor : adj[cur]){
            if(visited[neighbor] == 1 && neighbor != prev){
                return true;
            }
            else if(visited[neighbor] == 0){
                if(isCyclicUndirected(adj, visited, cur, neighbor)){
                    // false hole you have to traverse the entire graph
                    return true;
                }
            }
        }
        return false;
    }

    // regular DFS
    void dfs(vector<vector<int>>& adj, unordered_set<int>& visited, int cur){
        visited.insert(cur);
        for(auto &neighbor : adj[cur]){
            if(visited.count(cur) == 0){
                dfs(adj, visited, cur);
            }
        }
    }

    /// regular BFS
    // nodes numbered from 0 to adj.size() -1
    void bfs(vector<vector<int>>& adj){
        unordered_set<int> visited;
        queue<int> q;
        q.push(0); // starting from cur 0

        while(!q.empty()){
            int cur = q.front();
            q.pop();
            // queue te ek jinish multiple times thakte pare but this checks ensure our BFS
            // remains optimzied
            if(visited.count(cur) != 1){ // not visited yet
                visited.insert(cur); // mark as visited
                cout << cur << " "; // do the work
                for(int &n: adj[cur]){
                    if(visited.count(n) != 1){ // not visited yet
                        q.push(n);
                    }
                }
            }
        }
    }


    // num of componenets in a graph
    int numComponents(vector<vector<int>>& adj){
        int numComp = 0;
        unordered_set<int> visited;

        // We have nodes parent 0 to adj.size() - 1 in our graph
        for(int i = 0; i < adj.size(); ++i){
            if(visited.count(i) == 0){
                numComp++;
                dfs(adj, visited, i);
            }
        }

        return numComp;
    }

    // this is a modified dfs where we use graph coloring
    // for a graph to be bipartite there must be no odd edge cycle
    // to fulfill the above condition, if we start coloring the graaph with two colors, there
    // must not be 2 adjacent nodes with the same color. If there are 2 adjacent nodes with same
    // color, it means the graph has an odd edge cycle.
    // visited can be 3 values: -1 (unvisited), 0 (color 1), 1(color 2)
    // initially main func theke color = 1 pass hbe
    bool isBipartiteUtil(vector<vector<int>>& adj, vector<int>& visited, int cur, int color){
        visited[cur] = color; // color cur
        for(auto &neighbor: adj[cur]){
            int otherColor = 1 - visited[cur];
            if(visited[neighbor] == color){ // mane cur neighbor r color eki
                // odd edge cycle ase
                return false;
            }
            else if(visited[neighbor] == -1){ // not visited
                if(!isBipartiteUtil(adj, visited, neighbor, otherColor)){
                    return false;
                }
            }
        }
        // ekhono false return hoyni mane must be true
        return true;
    }

    // driver func
    bool isBipartite(vector<vector<int>>& graph) {
        vector<vector<int>> adj(graph.size());
        vector<int> visited(graph.size(), -1); // initialising with -1
        for(int i = 0; i < adj.size(); ++i){
            adj[i] = graph[i];
        }

        // eta kora lagtese for multi components
        for(int i = 0; i < adj.size(); ++i){
            if(visited[i] == -1){
                if(!isBipartiteUtil(adj, visited, i, 1)){ // everytime you start you send color 1
                    return false;
                }
            }
        }
        // ekhono false return hoyni mane must be true
        return true;
    }

    // A modified DFS with Stack where you have to store the cur vertex in the stack after doing
    // everything during backtracking
    // TopoSort algorithm
    void courseScheduleIIDFS(vector<vector<int>>& adj, vector<int>& visited, stack<int>& s, int
    cur){
        visited[cur] = 1;

        for(int neighbor: adj[cur]){
            if(visited[neighbor] == 0){
                courseScheduleIIDFS(adj, visited, s, neighbor);
            }
        }
        // after the recursion during backtracking stack e push korbo
        s.push(cur);
    }

    // driver
    vector<int> findOrder(int numCourses, vector<vector<int>>& prerequisites) {
        vector<vector<int>> adj(numCourses);
        vector<int> visited(numCourses, 0);
        vector<int> visited2(numCourses, 0);
        stack<int> s;
        vector<int> ans;

        for(auto &v : prerequisites){
            adj[v[1]].push_back(v[0]);
        }

        for(int i = 0; i < numCourses; ++i){
            if(visited[i] == 0){
                if(isCyclicDirected(adj, visited, i)){
                    return ans;
                }
            }
        }

        for(int i = 0; i < numCourses; ++i){
            if(visited2[i] == 0){
                courseScheduleIIDFS(adj, visited2, s, i); // initial value bepar na
                // it's not like tomake indegree 0 emon vertex diyei shuru kora lagbe
            }
        }

        while(!s.empty()){
            ans.push_back(s.top());
            s.pop();
        }
        return ans;
    }

    struct GraphN{
        int parent;
        int child;
        int weight;
        GraphN(int parent, int child, int cost): parent(parent), child(child), weight(cost) {}
    };

    struct childWeight{
        int child;
        int weight;
        childWeight(int c, int w): child(c), weight(w){}
    };

    struct myCmp{
        bool operator() (GraphN &n1, GraphN &n2){
            return n1.weight > n2.weight; // min heap tai grater than
            // max heap lagle <
        }
    };

    // prim's algorithm O(ElogV) as there will be at most E edges in the pq and logE is same
    //    // as logV^2 which is just logV
    // Regular BFS with prority queue
    // same idea like normal BFS: remove, see if it is visited or not, if not mark as visited,
    // work/print, add all
    // adjacent vertices which aren't visited
    // returns the min cost for the possible MST
    // 1 -> (2, 4), (4, 1) means 1 is connected to 2 and the edge weight is 4
    int Prim(unordered_map<int, vector<childWeight>>& adj, GraphN cur){
        priority_queue<GraphN, vector<GraphN>, myCmp> pq;
        unordered_set<int> visited;
        pq.push(cur);
        int minCost = 0;

        while (!pq.empty()){
            GraphN cur = pq.top();
            pq.pop();
            if(visited.count(cur.child) == 0){ // not visited yet
                visited.insert(cur.child); // mark as visited
                minCost += cur.weight;
                // done the work
                for(auto &n: adj[cur.child]){
                    // cur.child is the parent to n.child
                    if(visited.count(n.child) == 0){
                        pq.push(GraphN(cur.parent, n.child, n.weight));
                    }
                }
            }
        }
        return minCost;
    }

    // Undirected graph single component
    int primDriver(vector<vector<int>> &edgeList){
        unordered_map<int, vector<childWeight>> adj;
        for(auto &v: edgeList){
            adj[v[0]].push_back(childWeight(v[1], v[2]));
            adj[v[1]].push_back(childWeight(v[0], v[2]));
        }

        return Prim(adj, GraphN(-1, 0, 0)); // first parent -1, first to hocce any node you can
        // decide and first weight non-existent
    }

    // djikstra's algorithm O(ElogV) as there will be at most E edges in the pq and logE is same
    // as logV^2 which is just logV
    // Modified BFS with prority queue
    // 1 -> (2, 4), (4, 1) means 1 is connected to 2 and the edge weight is 4
    void Djikstra(unordered_map<int, vector<childWeight>>& adj, unordered_set<int>& visited,
    GraphN cur){
        // visited: 1 -> (3, 2) will be the format. It means
        // parent of 1 is 2 and the cost of this path is 2
        priority_queue<GraphN, vector<GraphN>, myCmp> pq;
        pq.push(cur);

        while (!pq.empty()){
            GraphN cur = pq.top();
            pq.pop();
            if(visited.count(cur.child) == 0){ // not visited yet
                visited.insert(cur.child); // mark as visited
                // done the work
                for(auto &n: adj[cur.child]){
                    // cur.child is the parent to n.child
                    if(visited.count(n.child) == 0){
                        pq.push(GraphN(cur.parent, n.child, cur.weight + n.weight)); // only diff
                        // with prim you add cur.weight with n.weight
                    }
                }
            }
        }
    }

     void djikstraDriver(vector<vector<int>> &edgeList){
        unordered_map<int, vector<childWeight>> adj;
        for(auto &v: edgeList){
            adj[v[0]].push_back(childWeight(v[1], v[2]));
            adj[v[1]].push_back(childWeight(v[0], v[2]));
        }
        unordered_set<int> visited;
        Djikstra(adj, visited, GraphN(-1, 0, 0)); // first parent -1 and we want the
        // djikstra path starting from 0

        // You can sent another visited structure if you want to print everything
     }

    // O(V^2logV)
    int primWithAdjMatrix(vector<vector<int>>& adj){
        // nodes running from 0 to adj.size() - 1
        // let's 0 our cur and start
        priority_queue<GraphN, vector<GraphN>, myCmp> pq;
        unordered_set<int> visited;
        pq.push(GraphN(-1, 0, 0)); // out starting cur Node
        int minCost = 0;

        while (!pq.empty()){
            GraphN cur = pq.top();
            pq.pop();
            if(visited.count(cur.child) == 0){ // not visited yet
                visited.insert(cur.child); // mark as visited
                // done the work
                for(int neigh = 0; neigh < adj.size(); ++neigh){
                    if(adj[cur.child][neigh] != INT32_MAX && visited.count(neigh) == 0){
                        pq.push(GraphN(cur.child, neigh, adj[cur.child][neigh]));
                    }
                }
            }
        }
        return minCost;
    }

    // Floyd Warshall Algorithm
    // All pairs shortest path
    // we use adjacency matrix
    // all the nodes are from 0 to adj.size() - 1
    // adj[i][j] contains the minimum cost/weight for traversing from i to j
    // if any cell has INT_MAX it means it is impossible to go from i to j
    // O(n^3)
    void floydWarshall(vector<vector<int>>& adj){
        for(int k = 0; k < adj.size(); ++k){
            for(int i = 0; i < adj.size(); ++i){
                for(int j = 0; j < adj.size(); ++j){
                    if(adj[i][k] == INT32_MAX || adj[k][j] == INT32_MAX){
                        continue;
                    }
                    else if(adj[i][k] + adj[k][j] < adj[i][j]){
                        adj[i][j] = adj[i][k] + adj[k][j];
                    }
                }
            }
        }
    }

    void floydWarshallDriver(vector<vector<int>>& adj){
        floydWarshall(adj);

        // can detect negative edge weight cycle
        for(int i = 0 ; i < adj.size(); ++i){
            if(adj[i][i] < 0){ // adj[i][i] should always be 0 if there's no negative edge weight
                // cycle
                cout << "There's a negative edge weight cycle" << endl;
            }
        }

        // now adj contains the shortest distance between all pairs
        // it's like a djikstra but for all the possible starting points
    }

    // Kosaraju's algorithm
    // O(V + E)
    // do DFS (store val in stack while backtracking), transpose the graph, DFS again and find
    // the strongly connected components while doing so
    // all the nodes are from 0 to adj.size() - 1

    // the idea is when we transpose a graph SCC r moddhe connection still thake properly and you
    // can traverse them but SCC and ar non-SCC r moddhe connection invert hoye jay which means
    // we have to initiate DFS multiple times for each SCC.
    // this can't be done norammly like the way we do in undirected graph bcz strongly connected
    // mane khali one way connection isn't enough
    void transpose(vector<vector<int>>& adj, vector<vector<int>>& rev){
        for(int i = 0; i < adj.size(); ++i){
            for(int &j: adj[i]){ //(i , j)
                rev[j].push_back(i); //(j, i)
            }
        }
    }

    // DFS + Stack
    void kosarjauDFS1(vector<vector<int>>& adj, unordered_set<int>& visited, stack<int>& s, int
    cur){
        visited.insert(cur);
        for(int &j: adj[cur]){
            if(!visited.count(j)){
                kosarjauDFS1(adj, visited, s, j);
            }
        }
        // while backtracking push to stack
        s.push(cur);
    }

    // a normal DFS
    // visit cur, mark as visited. For each vertex adjacent if it's not visited call dfs on them
    void kosarjauDFS2(vector<vector<int>>& adj, unordered_set<int>& visited, int cur){
        cout << cur << " ";
        visited.insert(cur);
        for(int &j: adj[cur]){
            if(!visited.count(j)){
                kosarjauDFS2(adj, visited,j);
            }
        }
    }

    void kosrajuDriver(vector<vector<int>>& edgeList){
        vector<vector<int>> adj(8);
        vector<vector<int>> rev(8);
        unordered_set<int> visited;
        stack<int> s;
        // directed graph
        // this algo is mainly used for finding the strongly connected components of a directed
        // graph
        for(auto &edge: edgeList){
            adj[edge[0]].push_back(edge[1]);
        }

        kosarjauDFS1(adj, visited, s, 0); // cur is 0
        transpose(adj, rev); // rev is the transpose

        visited.clear(); // making visited empty to use it again
        // using rev for this dfs
        int numSCC = 0;

        // this part of calling and counting is like the way we call and count to find the number
        // of components in an undirected graph but here we use stack to get the information
        while(!s.empty()) {
            int cur = s.top();
            s.pop();
            if(!visited.count(cur)){
                // each time you call the dfs2 there's a new component
                numSCC++;
                kosarjauDFS2(rev, visited, cur);
                cout << endl; // going to next line after printing a component
            }
        }
    }

    // a modified BFS where you continuously remove nodes of degree 1 from the graph until the
    // queue size is <= 2
    // no visited array needed for this bfs since you are removing visited edges
    // you are guaranted that edges are n - 1 and they will form a tree
    vector<int> findMinHeightTrees(int n, vector<vector<int>>& edges) {
        vector<int> ans;
        if(n == 1){
            ans.push_back(0);
            return ans;
        }

        vector<unordered_set<int>> adj(n);
        // make adjacency list
        for(auto &e: edges){ // undirected edge
            adj[e[0]].insert(e[1]);
            adj[e[1]].insert(e[0]);
        }

        queue<int> q;
        for(int i = 0; i < adj.size(); ++i){
            if(adj[i].size() == 1){ // nodes of degree 1
                q.push(i);
            }
        }

        while(n > 2){
            int levelLen = q.size();
            n = n - levelLen;
            for(int i = 0; i < levelLen; ++i){
                int cur = q.front();
                q.pop();
                for(auto &n: adj[cur]){ // adj[cur] will always have only 1 node
                    adj[n].erase(cur);
                    if(adj[n].size() == 1){
                        q.push(n);
                    }
                }
            }
        }

        while(!q.empty()){
            ans.push_back(q.front());
            q.pop();
        }

        return ans;
    }

    void djikstraNetworkDelay(vector<vector<GraphN>>& adj, vector<int>& timeTakenToReach, GraphN
    cur){
        priority_queue<GraphN, vector<GraphN>, myCmp> pq;
        unordered_set<int> visited;
        pq.push(cur);

        while(!pq.empty()){
            GraphN cur = pq.top();
            pq.pop();
            if(visited.count(cur.child) == 0){
                visited.insert(cur.child);
                timeTakenToReach[cur.child] = cur.weight;

                for(auto &n: adj[cur.child]){
                    if(visited.count(n.child) == 0){
                        pq.push(GraphN(n.parent, n.child, cur.weight + n.weight));
                    }
                }
            }
        }
    }


    // Djikstra algorithm
    // directed graph
    int networkDelayTime(vector<vector<int>>& times, int n, int k) {
        vector<vector<GraphN>> adj(n + 1);
        vector<int> timeTakenToReach(n + 1, INT32_MAX);
        for(auto &time: times){
            adj[time[0]].push_back(GraphN(time[0], time[1], time[2]));
        }

        djikstraNetworkDelay(adj, timeTakenToReach, GraphN(-1, k, 0));

        int minTime(0);
        for(int i = 1; i <= n; ++i){
            minTime = max(minTime, timeTakenToReach[i]);
        }

        return minTime == INT32_MAX ? -1 : minTime;
    }
};

class SolutionFindItinerary {
private:
    int numEdges;
public:
    // there is no idea of visited in this dfs as each node can be visited more than once but
    // each edge can only be visited once. so whenever we visit an edge we remove it
    // euler path with fixed starting
    // Backtracking + DFS
    // O(d^|E|) d is the max number of branches and E is the height
    bool dfsFindItenary(unordered_map<string, multiset<string>>& adj, vector<string>& ans, string cur){
        if(ans.size() == numEdges + 1){
            return true;
        }

        for(auto &neighbor: adj[cur]){ // we are using a set instead of a min heap bcz we can
            // iterate over a set which we need here but not a min heap
            ans.push_back(neighbor);
            adj[cur].erase(neighbor);
            if(dfsFindItenary(adj, ans, neighbor)){
                return true;
            }
            // backtracking
            adj[cur].insert(neighbor);
            ans.pop_back();
        }
        return false;
    }

    vector<string> findItinerary(vector<vector<string>>& tickets) {
        numEdges = tickets.size();
        unordered_map<string, multiset<string>> adj;
        vector<string> ans;
        ans.push_back("JFK"); // starting

        for(auto &v: tickets){
            adj[v[0]].insert(v[1]); // directed
        }

        dfsFindItenary(adj, ans, "JFK");
        return ans;
    }
};

// isCycleUndirected most optimized with Disjoint set data structure
// must be used for undirected graph
class disjointSet{
private:
    struct node{
        int parent;
        int rank;
        node(int p, int r) : parent(p), rank(r) {}
    };
    // make dsuf private variable so that it is available to all class funcs
    vector<node> dsuf; // ekhane vector nicci cz all nodes are 0, 1, 2, .. otherwise
    // unordered_map lagbe
    // child -> (parent, rank) structure thakbe


    int find(int v){ // return v's absolute parent
        if(dsuf[v].parent == -1){
            return v;
        }
        else{
            return dsuf[v].parent = find(dsuf[v].parent); // optimization we only care about
            // absolute parent
        }
    }

    void union_op(int fromP, int toP){
///     basic union
///     dsuf[fromP].parent = toP; // dsuf[fromP] storing a node

        // union by rank
        // rank jar less she point korbe to the other one
        // no change in rank
        if(dsuf[fromP].rank < dsuf[toP].rank){
            dsuf[fromP].parent = toP;
        }
        else if(dsuf[fromP].rank > dsuf[toP].rank){
            dsuf[toP].parent = fromP;
        }
        else{ // equal
            dsuf[fromP].parent = toP;
            // toP got pointed by fromP so rank of toP will increase
            dsuf[toP].rank++;
        }
    }

public:
    // time: O(E)
    // time for find and union is amortized O(1)
    // there's an advanced proof that time for find and union is a function is alpha(V). It's a
    // very slow growing func and the max value possible for our universe is 4
    bool isCyclicUndirected(vector<vector<int>>& edgeList, int numNodes){
        dsuf.resize(numNodes, node(-1, 0)); // additional space
        for(auto &e: edgeList){
            int fromP = find(e[0]);
            int toP = find(e[1]);

            if(fromP == toP){
                return true; // cycle peye gesi
                // already ekta connection ase and arekta connection banate boltese
            }

            union_op(fromP, toP); // ekta edge ase eta amader data structure r maddhome express
            // is a union_op
        }
        // gone through all edges ekhono cycle paini
        return false;
    }
};


#endif //ABIR_GRAPHPROBLEMS_H
