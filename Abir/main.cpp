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
//#include "GraphProblems.h"
//#include "GreedyProblems.h"
//#include "HeapProblems.h"
//#include "IntervalProblems.h"
//#include "TwoPointers.h"
//#include "LinkedList.h"
#include "AmazonProblems.h"

using namespace std;

int main() {
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
}
