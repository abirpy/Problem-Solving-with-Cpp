//
// Created by Tahsinul Haque Abir on 7/31/22.
//

#ifndef ABIR_HEAPPROBLEMS_H
#define ABIR_HEAPPROBLEMS_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;

class HeapProblems {
public:
    // like the problem finding kth smallest element from two sorted linked lists
    // uses merging algorithm
    // now we basically have r sorted linked lists and we need to use a min heap to keep track of
    // the min element
    // the cols are also so sorted which means even in the worst case we have to go through the
    // first k rows to get the kth smallest element. Actually we have to go through min(colNum,
    // k) rows
    int kthSmallest(vector<vector<int>>& matrix, int k) {
        priority_queue<int, vector<int>, greater<int>> pq;
        int rows, cols;
        rows = matrix.size();
        cols = matrix[0].size();

        int minK;
        for(int c = 0; c < cols; ++c){
            for(int r = 0; r < min(rows, k); ++r){
                pq.push(matrix[r][c]);
            }
            minK = pq.top();
            pq.pop();
            k--;
            if(k == 0){  //optimization
                return minK;
            }
        }

        while(k != 0){
            minK = pq.top();
            pq.pop();
            k--;
        }

        return minK;
    }

    struct myCmpKSmallest{
        bool operator()(vector<int>& v1, vector<int>& v2){
            return v1[0] + v1[1] >= v2[0] + v2[1]; // min heap lagbe so ulta >=(descending)
        }
    };

    vector<vector<int>> kSmallestPairs(vector<int>& nums1, vector<int>& nums2, int k) {
        vector<vector<int>> ans;
        priority_queue<vector<int>, vector<vector<int>>, myCmpKSmallest> pq;

        for(int i = 0; i < nums1.size(); ++i){
            for(int j = 0; j < nums2.size(); ++j){
                pq.push(vector<int>({nums1[i], nums2[j]}));
            }
        }
        int i = 0;
        while(i != k && !pq.empty()){
            ans.push_back(pq.top());
            pq.pop();
            i++;
        }
        return ans;
    }

    // VVI INFO about heap formation
    // https://leetcode.com/problems/top-k-frequent-elements/solution/
    // time required to build a heap of size k from a hash map of size n is O(nlogk)
    // we push the first k elem into the heap which is klogk and then for the rest n - k elem, we
    // check one by one if they are less/greater than top of the heap and then push/pop which is
    // (n -k)logk in total it's O(nlogk) where k < n. It's much more optimal thank nlogn


    // You can't sort a hashmap directly. You can copy all pairs into a vector and then try
    // sorting the vector.
    // You can also build a heap directly through priority queue constructor in O(n). This is
    // more optimal
};


#endif //ABIR_HEAPPROBLEMS_H
