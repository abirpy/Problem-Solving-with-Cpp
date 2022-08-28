//
// Created by Tahsinul Haque Abir on 8/1/22.
//

#ifndef ABIR_INTERVALPROBLEMS_H
#define ABIR_INTERVALPROBLEMS_H
#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;

class IntervalProblems {
public:
    // Merge Intervals
    // Greedy approach sort on the basis of first element of the interval
    vector<vector<int>> merge(vector<vector<int>>&  intervals) {
        vector<vector<int>> ans;
        sort(intervals.begin(), intervals.end());
        ans.push_back(intervals[0]); // this will give us info about the prev later
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

    // Merge Intervals (Already sorted)
    // Greedy approach sort on the basis of first element of the interval
    vector<vector<int>> mergeWithoutSort(vector<vector<int>>&  intervals) {
        vector<vector<int>> ans;
        ans.push_back(intervals[0]); // this will give us info about the prev later
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

    // interval problem
    // meeting room I
    bool canAttendMeetings(vector<vector<int>>& intervals) {
        if(intervals.size() == 0){
            return true;
        }
        // a person can attend all the meeting if there is no overlapping meeting
        // you can sort the intervals first and use a prev to find if there is any overlap while
        // iterating through the intervals
        sort(intervals.begin(), intervals.end());//default sort will sort based on first value
        int prev = intervals[0][1]; // we'll track the ending time of the last meeting
        for(int i = 1; i < intervals.size(); ++i){
            if(intervals[i][0] < prev){// overlapping hocce
                return false;
            }
            prev = intervals[i][1];
        }
        return true; // ekhono overlap payni
    }

    // interval problem
    // another min heap problem
    // first sort on the basis of starting time and then keep the meetings in min heap based on
    // ending times
    int minMeetingRooms(vector<vector<int>>& intervals) {
        sort(intervals.begin(), intervals.end()); // first value r basisi e default sort hbe
        priority_queue<int, vector<int>, greater<int>> pq;
        pq.push(intervals[0][1]); // ending time of the first meeting which started the first
        for(int i = 1; i < intervals.size(); ++i){
            // check if the top of the heap which is supposed to be free at first is fee or not
            if(intervals[i][0] >= pq.top()) { // mane free
                // ei room e khn intervals[i] meeting hbe
                pq.pop();
                pq.push(intervals[i][1]);
            }
            else{ // not free
                // assign new room to intervals[i] which is gonna end at intervals[i][1]
                pq.push(intervals[i][1]);
            }
        }
        return pq.size();
    }

    // neetcode
    int eraseOverlapIntervals(vector<vector<int>>& intervals) {
        // goal: find the maxNumberOfOverlapping
        int numErase = 0;
        int prev = intervals[0][1]; // just need to trak the end of the last one
        for(int i = 1; i < intervals.size(); ++i){
            if(intervals[i][0] < prev){ // condition of overlapping
                // mane you have to remove one now to avoid overlapping
                // you remove the one which ends later and keep the one which ends earlier; common
                // sense this way you provide less chance of collision
                // we actually don;t remove we simulate removal through prev
                numErase++;
                prev = min(prev, intervals[i][1]);
            }
            else{ // no collision
                // move prev to the end of current interval
                prev = intervals[i][1];
            }
        }
        return numErase;
    }

    // it's like the way merging algorithm find intersection
    vector<vector<int>> intervalIntersection(vector<vector<int>>& firstList, vector<vector<int>>& secondList) {
        vector<vector<int>> ans;
        // you take two pointers
        int i = 0, j = 0;
        while(i < firstList.size() && j < secondList.size()){
            if(secondList[j][0] <= firstList[i][1] && secondList[j][1] >= firstList[i][0]){
                vector<int> tmp(2);
                tmp[0] = max(firstList[i][0], secondList[j][0]);
                tmp[1] = min(firstList[i][1], secondList[j][1]);
                ans.push_back(tmp);
            }
            if(firstList[i][1] < secondList[j][1]){
                i++;
            }
            else if(firstList[i][1] > secondList[j][1]){
                j++;
            }
            else{
                i++;
                j++;
            }
        }
        return ans;
    }

    vector<vector<int>> insert(vector<vector<int>>& intervals, vector<int>& newInterval) {
        if(intervals.empty()){
            intervals.push_back(newInterval);
            return intervals;
        }
        for(auto it = intervals.begin(); it != intervals.end() + 1; ++it){
            if(it != intervals.end() && newInterval[0] <= (*it)[1] && newInterval[1] >= (*it)[0]){
                //overlapping
                // merge the inervals
                (*it)[0] = min((*it)[0], newInterval[0]);
                (*it)[1] = max((*it)[1], newInterval[1]);
                return mergeWithoutSort(intervals);
            }
            else if((it == intervals.begin() || (*(it - 1))[1] < newInterval[0]) &&
                    (it == intervals.end() || newInterval[1] < (*it)[0])){
                // no overlapping
                intervals.insert(it, newInterval);
                break;
            }
        }
        return intervals;
    }
};


#endif //ABIR_INTERVALPROBLEMS_H
