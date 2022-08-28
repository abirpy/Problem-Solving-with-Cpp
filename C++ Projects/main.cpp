#include <iostream>
#include <vector>
//#include "Islands.h"
//#include "ValidTree.h"
//#include "CourseSchehdule.h"
//#include "CourseScheduleII.h"
//#include "MHT.h"
//#include "DistanceK.h"
//#include "Backtracking.h"
//#include "LetterCasePerm.h"
//#include "Subset.h"
//#include "Combinations.h"
//#include "CombinationSum.h"
//#include "PalindromePartitioning.h"
//#include "LetterCombination.h"
//#include "GeneralizedAbbreviation.h"
//#include "TreeProblems.h"
//#include "DynamicProgramming.h"
//#include "GreedyProblems.h"
//#include "HeapProblems.h"
//#include "IntervalProblems.h"
//#include "TwoPointers.h"
//#include "LinkedList.h"
//#include "AmazonProblems.h"
//#include "GraphProblems.h"
#include "MediumProblems.h"

using namespace std;

int findUniqueValues(vector<int> experience) {
    sort(experience.begin(), experience.end());
    unordered_set<double> hashSet;

    int l = 0, r = experience.size() - 1;

    while(l < r){
        cout << experience[l] + experience[r] << endl;
        double average = double((experience[l] + experience[r])) / 2;
        if(hashSet.count(average) == 0){
            hashSet.insert(average);
        }
        l++;
        r--;
    }
    return hashSet.size();
}

int maxLength(vector<int> a, int k) {
    int maxLen = 0;
    int l = 0, r = 0;
    int sum = a[l];
    while(r < a.size()){
        while(sum <= k){
            maxLen = max(maxLen, r - l + 1);
            r++;
            if(r < a.size()){
                sum += a[r];
            }
        }
        while(l <= r && sum > k){
            sum -= a[l];
            l++;
        }
    }
    return maxLen;
}

// regular DFS
void dfs(vector<vector<int>>& adj, unordered_set<int>& visited, int cur){
    visited.insert(cur);
    for(auto &neighbor : adj[cur]){
        if(visited.count(neighbor) == 0){
            dfs(adj, visited, neighbor);
        }
    }
}


vector<int> getTheGroups(int n, vector<string> queryType, vector<int> students1, vector<int> students2) {
    vector<vector<int>> adj(n + 1);
    vector<int> ans;
    unordered_set<int> visited;
    int sumNumGroups = 0;

    for(int i = 0; i < queryType.size(); ++i){
        if(queryType[i] == "Total"){
            int student1 = students1[i];
            int student2 = students2[i];
            visited.clear();
            dfs(adj, visited, student1);
            sumNumGroups += visited.size();

            dfs(adj, visited, student2);
            sumNumGroups += visited.size() - sumNumGroups;

            ans.push_back(sumNumGroups);
        }
        else {
            adj[students1[i]].push_back(students2[i]);
            adj[students2[i]].push_back(students1[i]);
        }
    }

    return ans;
}


int getPalindromesCount(string s) {
    int len = s.size();
    vector<vector<int>> dp(len + 1, vector<int>(len + 1, 0));

    for(int i = 0; i < len; ++i){
        dp[i][i] = 1;
    }

    for(int l = 2; l <= len; ++l){
        for (int i = 0; i <= len - l; i++) {
            int k = l + i - 1;
            if (s[i] == s[k])
                dp[i][k] = dp[i][k - 1] + dp[i + 1][k] + 1;
            else
                dp[i][k] = dp[i][k - 1] + dp[i + 1][k]
                           - dp[i + 1][k - 1];
        }
    }

    return dp[0][len - 1];
}



int main() {
    string s = "010110";
    cout << getPalindromesCount(s);
//    cout << 9 / 2;
//
//    cout << findUniqueValues(experience);

//    vector<vector<char>> board({
//        {'A'}
//    });

//    Islands i;
//    cout << i.numIslands(v) << endl;
//    ValidTree t;
//    cout << t.validTree(5, v) << endl;
//    CourseSchehdule c;
//    cout << c.canFinish(2, v) << endl;
//    CourseScheduleII c2;
//    c2.findOrder(2, v);
//    MHT t;
//    t.findMinHeightTrees(6, v);
//
//    TreeNode* one = new TreeNode(5);
//    TreeNode* two = new TreeNode(4);
//    TreeNode* three = new TreeNode(6);
//    TreeNode* four = new TreeNode(3);
//    TreeNode* five = new TreeNode(7);
////    TreeNode* six = new TreeNode(15);
////    TreeNode* seven = new TreeNode(7);
//    one->left = two;
//    one->right = three;
//    two->left = nullptr;
//    two->right = nullptr;
//    three->left = four;
//    three->right = five;
//
//    TreeProblems t;
////    cout << t.kthSmallest(one, 1);
////    t.pathSumII(one, -1);
////    bool flag = one == t.lowestCommonAncestorBinaryTree(one, three, two);
////    cout << flag;
//    t.isValidBST(one);


//
//    levelOrderBottom l;
//    l.levelOrderBottomf(one);
//    DistanceK d;
//    d.distanceK(one, four, 4);

//    Backtracking b;
//    cout << b.exist(board, "A");
//    LetterCasePerm l;
//    l.letterCasePermutation(string("a1b2"));
//    Subset s;
//    vector<int> v({1, 2, 2});
//    s.subsets(v);
/*    Combinations c;
    c.combine(4, 2);*/
//    vector<int> candidates({1, 2, 3});
//    CombinationSum c;
//    c.combinationSum(candidates, 3);
//    PalindromePartitioning p;
//    p.partition("aab");
//    LetterCombination l;
//    l.letterCombinations("7");
//    GeneralizedAbbreviation g;
//    g.generateAbbreviations("wor");
//    TreeProblems t;
////    vector<int> nums({1, 1, 1, 1, 1});
////    t.findTargetSumWays(nums, 3);
////    DynamicProgramming d;
////    vector<int> v({7,1,5,3,6,4});
////    cout << d.countNumSubsetSumTab(v, 3, 4) << endl;
////    cout << d.coinChangeITab(v, 4, 0) << endl;
////    cout << d.coinChangeIITab(v, 3, 4) << endl;
////    cout << d.rodCuttingTab(v, 4, 4) << endl;
////    string s1 = "AAA", s2 = "", lcs = "";
////    int lcsLen = d.LCS(s1, s2, lcs);
////    cout << lcs.substr(0, lcsLen);
////
////    cout << d.findTargetSumWays(v, -200);
////    cout << d.robTab2(v) << endl;
////    cout << d.coinChange(v, 11);
////    cout << d.lengthOfLIS(v) << endl;
////    cout << d.longestPalindrome("babad") << endl;
////    cout << d.numDecodings("11106");
////    cout << d.canJump(v);
////    cout << d.findNumberOfLIS(v);
////    cout << d.canPartitionKSubsets(v, 4);
////    cout << d.maxProfit(v) << endl;
////    cout << d.maxProfit2(v);
//
//    vector<vector<int>> points({
//                                       {1, 2}, {1, 3}, {5, 6}, {4, 10}
//    });
////    GraphProblems g;
////    g.isBipartite(v);
////    g.findOrder(2, v);
////    vector<vector<string>> tickets({
////        {"MUC","LHR"},{"JFK","MUC"},{"SFO","SJC"},{"LHR","SFO"}
////    });
////    SolutionFindItinerary s;
////    s.findItinerary(tickets);
////    g.djikstraDriver(v);
////    g.kosrajuDriver(edgeList);
////    g.findMinHeightTrees(6, edgeList);
//    vector<int> v({3, 0, 2, 1, 2});
//    GreedyProblems g;
//    g.canReach(v, 2);
//    g.findMinArrowShots(points);
//    string s = "aab";
////    g.reorganizeString(s);
////    g.merge(points);
//    g.employeeFreeTime(points);

//    HeapProblems h;
//    vector<vector<int>> matrix({
//       {1, 2}, {3, 3}
//    });
//
//    h.kthSmallest(matrix, 2);
//    IntervalProblems i;
//    vector<vector<int>> v1({
//        {0,2},{5,10},{13,23},{24,25}
//    });
//    vector<vector<int>> v2({
//        {1,5},{8,12},{15,24},{25,26}
//    });
//    i.intervalIntersection(v1, v2);

//    vector<int> nums({0, 1, 2, 2});
//    TwoPointers t;
////    t.sortedSquares(nums);
////    t.findDuplicate(nums);
////    t.threeSum(nums);
////    t.threeSumClosest(nums, 1);
////    t.totalFruit(nums);
//
//    LinkedList l;
//    ListNode* head = new ListNode(1);
//    head->next = new ListNode(4);
//    head->next->next = new ListNode(5);
//
//    ListNode* head2 = new ListNode(1);
//    head->next = new ListNode(3);
//    head->next->next = new ListNode(4);
//
//    ListNode* head3 = new ListNode(2);
//    head->next = new ListNode(6);
//
//    vector<ListNode*> lists({head, head2, head3});

//    l.reverseBetween(head, 2, 3);
//    l.swapPairs(head);
//    l.sortList(head);
//    vector<vector<int>> p({
//        {1, 0}
//    });
//    AmazonProblems a;
//    a.uniquePathsWithObstacles(p);
//    a.mergeKLists(lists);
//
//
//    vector<vector<int>> edges({
//      {2,1,1},
//      {2,3,1},
//      {3,4,1}
//    });
//
//    GraphProblems g;
//    g.networkDelayTime(edges, 4, 2);
//    MediumProblems m;
////    vector<vector<int>> v({
////                                  {1, 2, 3},
////                                  {2, 5, 6},
////                                  {8, 9, 5}
////    });
////    m.rotate(v);
//    LL l;
//    ListNode* head = new ListNode(1);
//    head->next = new ListNode(2);
//    l.removeNthFromEnd(head, 1);
}
