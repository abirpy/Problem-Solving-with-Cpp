//
// Created by Tahsinul Haque Abir on 7/30/22.
//

#ifndef ABIR_GREEDYPROBLEMS_H
#define ABIR_GREEDYPROBLEMS_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;

class GreedyProblems {
public:
/////     Buy and Sell Stock II (A greedy problem actually)
    // DP solution (not optimal)
    int maxProfitUtil(vector<int>& prices, vector<vector<int>>& mem, int idx, int buy){ // int buy
        // is a bool 1 mane buying
        // mode nahole 0 mane selling mode
        if(idx == prices.size()){
            return 0;
        }
        if(mem[idx][buy] != -1){
            return mem[idx][buy];
        }

        if(buy == 1){
            return mem[idx][buy] = max(maxProfitUtil(prices, mem, idx + 1, buy - 1) - prices[idx],
                       maxProfitUtil(prices, mem, idx + 1, buy));
        }
        else { // buy == 0
            return mem[idx][buy] = max(maxProfitUtil(prices, mem, idx + 2, buy + 1) + prices[idx],
                       maxProfitUtil(prices, mem, idx + 1, buy));
        }

    }
    int maxProfit(vector<int>& prices) {
        vector<vector<int>> mem(prices.size(), vector<int>(2, -1));
        return maxProfitUtil(prices, mem, 0, 1);
    }

    // can do any number of transaction I like
    // so the strategy is always buy when prices[i] < prices[i + 1] and sell on the next day
    // greedy strategy
    // jehetu any amount of transaction kora jabe tai shb gula low-high theke profit korle beshi
    // laav hbe
    int maxProfitGreedy(vector<int>& prices){
        int maxProfit = 0;
        for(int i = 0; i < prices.size() - 1; ++i){
            if(prices[i] < prices[i + 1]){
                maxProfit += prices[i + 1] - prices[i];
            }
        }
        return maxProfit;
    }

    ///// Best Time to Buy and Sell Stock III
    // At most 2 transactions
    // this is a DP problem
    int maxProfitIIIUtil(vector<int>& prices, vector<vector<vector<int>>>& mem, int idx, int buy,
                         int numTransactions) {
        if(idx >= prices.size() || numTransactions == 2){
            return 0;
        }
        if(mem[idx][buy][numTransactions] != -1){
            return mem[idx][buy][numTransactions];
        }

        if(buy == 1){
            return mem[idx][buy][numTransactions] = max(maxProfitIIIUtil(prices, mem, idx + 1, buy -
            1, numTransactions) - prices[idx], maxProfitIIIUtil(prices, mem, idx + 1, buy, numTransactions));
        }
        else { // buy == 0 // whenever sold numTransaction++
            return mem[idx][buy][numTransactions] = max(maxProfitIIIUtil(prices, mem, idx + 1, buy
            + 1, numTransactions + 1) + prices[idx], maxProfitIIIUtil(prices, mem, idx + 1, buy,
                                                                      numTransactions));
        }
    }
    int maxProfitIII(vector<int>& prices) {
        vector<vector<vector<int>>> mem(prices.size(), vector<vector<int>>(2, vector<int>(3, -1)));
        return maxProfitIIIUtil(prices, mem, 0, 1, 0);
    }

    ////Jump GameII
    // min number of jumps required to reach the end
    // A greedy bfs
    int jump(vector<int>& nums) {
        int l, r;
        l = r = 0;
        int numJumps = 0;

        while (r < nums.size() - 1){
            int farthest = r + 1; // min possible value of farthest
            for(int i = l; i <= r; ++i){
                farthest = max(farthest, i + nums[i]);
            }
            l = r + 1; // next window
            r = farthest;
            numJumps++;
        }
        return numJumps;
    }

    //// Jump Game III
    //// A weird DFS
    //Given an array of non-negative integers arr, you are initially positioned at start index of the array. When you are at index i, you can jump to i + arr[i] or i - arr[i], check if you can reach to any index with value 0.
    // Notice that you can not jump outside of the array at any time.
    // This is actually a dfs problem where you have to mark things as visited
    // not a dp problem, a regular recursive dfs
    // non negative value in arr. so negative marking kore visited rakha jabe
    bool canReach(vector<int>& arr, int start) {
        if(start >= 0 && start < arr.size() && arr[start] >= 0){
            if(arr[start] == 0){
                return true;
            }
            arr[start] = -arr[start]; // marking as visited
            return canReach(arr, start + arr[start]) || canReach(arr, start - arr[start]);
        }
        // ekhono true return hoini mane false you got trapped in a cycle
        return false;
    }

    // weird greedy problem
    //best explanation: https://www.youtube.com/watch?v=Z2Plc8o1ld4
    // O(n) where n is the size of tasks
    int leastInterval(vector<char>& tasks, int n) {
        unordered_map<char, size_t> hash;
        size_t maxFreq = 0;
        size_t countNumCharWithMaxFreq = 0;
        for(char &c: tasks){
            hash[c]++;
            maxFreq = max(maxFreq, hash[c]);
        }

        for(auto &p: hash){ //iterating through hash pairs
            if(p.second == maxFreq){
                countNumCharWithMaxFreq++;
            }
        }

        // formula: number of groups * group_size + countNumCharWithMaxFreq
        return (int)max(tasks.size(), ((maxFreq - 1) * (n + 1) +countNumCharWithMaxFreq));
    }
//
//    struct myCmp{
//
//    };
    // the comparator should always be static
    static bool myCmp(vector<int> v1, vector<int> v2){
        return v1[1] <= v2[1]; // return true if second value is in increasing
        // order and then first value in increasing order
    }
    // a greedy problem which requires sorting first based on end points
    // O(nlogn)
    int findMinArrowShots(vector<vector<int>>& points) {
        sort(points.begin(), points.end(), myCmp);
        // second value in increasing order and then first value in increasing orer
        int e = points[0][1];
        int arrow = 1; // at least 1 ta arrow to lagbei
        for(auto &point: points){
            if(point[0] <= e){ // all starting points which are <= current e will have
                // overlapping or
                // touching
                continue;
            }
            // ekhane ashbo for a point whose s > e
            e = point[1]; // new ending point
            arrow++; // notun ending mark mane one more arrow lagbe
        }

        return arrow;
    }

    // pq e comparator always struct r moddhe lekhte hoy
    struct myCmpReorganizeString{
    // max heap requires ascending order
    // min heap requires descending order
        bool operator()(pair<char, int>& p1, pair<char, int>& p2){
            return p1.second <= p2.second; //we need a max heap mane ascending order lagbe that's
            // why ulta
            // <= sign. Sorting chaile >= ditam
        }
    };
    typedef pair<char, int> Pair;
    // greedy + heap problem related to using frequency map
    string reorganizeString(string s) {
        unordered_map<char, int> freq;
        for(auto &c: s){
            freq[c]++;
        }

        // now make a max heap
        priority_queue<Pair, vector<Pair>, myCmpReorganizeString> maxHeap
        (freq.begin(), freq.end()); // making maxHeap through heapify O(n)

        string ans = "";
        Pair prev({' ', 0}); // there's nothing like '' empty char bujhate gele ' '
        while (!maxHeap.empty()){
            Pair cur = maxHeap.top();
            maxHeap.pop();
            if(prev.second > 0){
                maxHeap.push(prev);
            }

            ans += cur.first; // max_freq char add korlam
            cur.second--; // ekbar use holo
            prev = cur;
        }

        if(prev.second != 0){ // ekhono elem baki ase but queue empty hoye gese mane not possible
            return "";
        }
        return ans;
    }

    // eta ashole default way so eta lekha lagbe na
    static bool myCmpMerge(vector<int>& v1, vector<int>& v2){
        return v1[0] <= v2[0]; // ascending order e sorting hbe
    }

    // Merge Intervals
    // Greedy approach sort on the basis of first element of the interval
    vector<vector<int>> merge(vector<vector<int>>&  intervals) {
        vector<vector<int>> ans;
        sort(intervals.begin(), intervals.end(), myCmpMerge);
        ans.push_back(intervals[0]);
        for(int i = 1; i < intervals.size(); ++i){
            int lastOfPrev = ans.back()[1];
            if(intervals[i][0] <= lastOfPrev){
                ans.back()[1] = (int)max(ans.back()[1], intervals[i][1]);
            }
            else{
                ans.push_back(intervals[i]);
            }
        }
        return ans;
    }

    // variant of the merge interval problem
    vector<vector<int>> employeeFreeTime(vector<vector<int>>& schedule) {
        // merge intervals first
        vector<vector<int>> mergedSchedule(merge(schedule));
        vector<vector<int>> ans;
        vector<int> tmp;
        for(int i = 0; i < mergedSchedule.size(); ++i){
            for(int j = 0; j < mergedSchedule[i].size(); ++j){
                if((i == 0 && j == 0) || (i == mergedSchedule.size() - 1 && j == 1)){
                    continue;
                }
                tmp.push_back(mergedSchedule[i][j]);
                if(tmp.size() == 2){
                    ans.push_back(tmp);
                    tmp.clear();
                }
            }
        }
        return ans;
    }

    // a greedy problem which uses sorting and heap
    // the idea is simple. We will choose the course first whose deadline is before.
    // while taking this decision and moving forward if we cross the deadline anytime, we will
    // remove the course requiring max number of hours from the heap as that would allow us to
    // keep more NUMBER OF COURSES
    // return the heap size at the end
    static bool myCmpCourses(vector<int>& v1, vector<int>& v2){
        return v1[1] <= v2[1]; // ascending order e sorting hbe
    }
    //time: O(nlogn)
    int scheduleCourse(vector<vector<int>>& courses) {
        sort(courses.begin(), courses.end(), myCmpCourses);
        priority_queue<int> pq; //maxHeap
        int time = 0;

        for(auto& course: courses){
            time += course[0]; // adding time required to finish the course
            pq.push(course[0]);
            while(time > course[1]){ // deadline is course[1]
                int max = pq.top();
                pq.pop(); // removing the max
                time -= max;
            }
        }
        return pq.size();
    }
};


#endif //ABIR_GREEDYPROBLEMS_H
