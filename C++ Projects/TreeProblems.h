//
// Created by Tahsinul Haque Abir on 7/15/22.
//

#ifndef ABIR_TREEPROBLEMS_H
#define ABIR_TREEPROBLEMS_H

#include <iostream>
#include "string"
#include <vector>
#include <queue>
#include <stack>
#include <unordered_map>
#include <queue>
#include <algorithm>

using namespace std;

struct TreeNode {
    int val;
    TreeNode *left;
    TreeNode *right;
    TreeNode() : val(0), left(nullptr), right(nullptr) {}
    TreeNode(int x) : val(x), left(nullptr), right(nullptr) {}
    TreeNode(int x, TreeNode *left, TreeNode *right) : val(x), left(left), right(right) {}
};

class TreeProblems {
public:
    bool isSameTreeF(TreeNode* p, TreeNode* q) {
        if(p == nullptr && q == nullptr){
            return true;
        }
        else if(p == nullptr && q != nullptr){
            return false;
        }
        else if(p != nullptr && q == nullptr){
            return false;
        }
        else{
            return (p->val == q->val) && isSameTreeF(p->left, q->left) && isSameTreeF(p->right,
                                                                                      q->right);
        }
    }
    bool hasPathSum(TreeNode* root, int targetSum) {
        if(root == nullptr){
            return false;
        }
        else if(root->left == nullptr && root->right == nullptr){
            return root->val == targetSum;
        }
        else if(root->left != nullptr && root->right == nullptr){
            return hasPathSum(root->left, targetSum - root->val);
        }
        else if(root->left == nullptr && root->right != nullptr){
            return hasPathSum(root->right, targetSum - root->val);
        }
        else{
            return hasPathSum(root->left, targetSum - root->val) || hasPathSum(root->right, targetSum - root->val);
        }
    }

    int heightOfBinaryTree(TreeNode* root, int& maxDiameter) {
        if(root == nullptr){
            return 0;
        }
        else{
            int heightOfLeftTree = heightOfBinaryTree(root->left, maxDiameter);
            int heightOfRightTree = heightOfBinaryTree(root->right, maxDiameter);

            // calculate maxDiameter
            // heightOfLeftTree + heightOfRightTree is the diameter of the tree with the current
            // node as the centre
            maxDiameter = max(maxDiameter, heightOfLeftTree + heightOfRightTree);

            return 1 + max(heightOfLeftTree, heightOfRightTree);
        }
    }

    int diameterOfBinaryTree(TreeNode* root) {
        int maxDiameter = 0;
        heightOfBinaryTree(root, maxDiameter);
        return maxDiameter;
    }

    TreeNode* mergeTrees(TreeNode* root1, TreeNode* root2) {
        TreeNode* root3 = nullptr;
        if(root1 == nullptr && root2 == nullptr){
            return root3;
        }
        else if(root1 != nullptr && root2 == nullptr){
            root3 = root1;
        }
        else if(root1 == nullptr && root2 != nullptr){
            root3 = root2;
        }
        else{
            root3 = new TreeNode(root1->val + root2->val);
            root3->left = mergeTrees(root1->left, root2->left);
            root3->right = mergeTrees(root1->right, root2->right);
        }
        return root3;
    }

    // of a Binary Search Tree
    TreeNode* lowestCommonAncestor(TreeNode* root, TreeNode* p, TreeNode* q) {
        if(p->val <= root->val && q->val >= root->val){
            return root;
        }
        else if(p->val < root->val && q->val < root->val){
            return lowestCommonAncestor(root->left, p, q);
        }
        else{
            return lowestCommonAncestor(root->right, p, q);
        }
    }

    // time: O(Size of root tree * size of subRoot tree)
    // in the worst case for each node of main tree you have to call isSameTreeF and traverse the
    // entire sub-tree
    bool isSubtree(TreeNode* root, TreeNode* subRoot) {
        // if subRoot is a sub tree of root tree, it must exist as a same tree with either the
        // main tree or one of its childs
        if(root == nullptr && subRoot == nullptr){
            return true;
        }
        else if(root == nullptr && subRoot != nullptr){
            return false;
        }
        else if(root != nullptr && subRoot == nullptr){
            return true;
        }
        else{
            if(isSameTreeF(root, subRoot)){
                return true;
            }
            else{
                return isSubtree(root->left, subRoot) || isSubtree(root->right, subRoot);
            }
        }
    }

    TreeNode* invertTree(TreeNode* root) {
        if(root == nullptr){
            return root;
        }
        // someone will invert the left and right sub-tree for me
        // then I will just horizontally flip( or swap) those trees keeping the root at the centre
        TreeNode* tmp = root->right;
        root->right = invertTree(root->left);
        root->left = invertTree(tmp);
        return root;
    }

    void kthSmallestUtil(TreeNode* root, int k, vector<int>& val) {
        //inorder traversal
        if(root == nullptr){
            return;
        }

        kthSmallestUtil(root->left, k, val); // left
        if(val.size() < k){
            val.push_back(root->val); // process
        }
        kthSmallestUtil(root->right, k, val); // right
    }

    int kthSmallest(TreeNode* root, int k) {
        int count = 0;
        vector<int> val; // you have to store the entire inorder traversal in the val
        kthSmallestUtil(root, k, val);
        return val[k - 1];
    }

    // this algo is O(n^2)
    // The number of possible path parent root to leaf is at most N/2 (the number of leaves).
    // N/2 * N = N^2 since you have to push_back N/2 tmp arrays each of length N(skewed binary tree)
    // in the
    // worst case
    void pathSumUtil(TreeNode* root, int targetSum, vector<int>& tmp, vector<vector<int>>& output) {
        if(root->left == nullptr && root->right == nullptr && targetSum == 0){ // reached the leaf
            output.push_back(tmp);
            // this line is O(n). we have to copy tmp to output
            return;
        }

        // include left child
        // root->left->val <= targetSum ei optimization kora jeto when we have only positive val
        // in the tree. Tree te negative val thakle eta kora jabe na
        if(root->left){
            tmp.push_back(root->left->val);
            pathSumUtil(root->left, targetSum - root->left->val, tmp, output);
            tmp.pop_back();
        }


        if(root->right){
            tmp.push_back(root->right->val);
            pathSumUtil(root->right, targetSum - root->right->val, tmp, output);
            tmp.pop_back();
        }
    }

    vector<vector<int>> pathSum(TreeNode* root, int targetSum) {
        vector<vector<int>> output;
        if(root == nullptr){
            return output;
        }
        vector<int> tmp;
        tmp.push_back(root->val);
        pathSumUtil(root, targetSum - root->val, tmp, output);
        return output;
    }


    // return num of subarray with sum k in O(n)
    // ei structure ta important otherwise edge cases where k = 0 pass hbe na
    int subarraySum(vector<int>& nums, int k) {
        unordered_map<int, int> hash;
        int count = 0;
        int prefixSum = 0;
        hash[0] = 1; // mane 0 always exists

        for(int i = 0; i < nums.size(); ++i){
            prefixSum += nums[i];

            if(hash.count(prefixSum - k) == 1){
                count += hash[prefixSum - k];
            }
            hash[prefixSum]++; // hash e value change pore kora lagbe after checking
        }
        return count;
    }

    void pathSumIIUtil(TreeNode* root, int targetSum, long int& prefixSum, int& count,
                      unordered_map<long int, int>& hash) {
        if(root == nullptr){
            return;
        }
        // counting at the beginning
        prefixSum += root->val;
        if(hash.count(prefixSum - targetSum) == 1){
            count += hash[prefixSum - targetSum];
        }
        hash[prefixSum]++;

        // recursive calls
        pathSumIIUtil(root->left, targetSum, prefixSum, count, hash);
        pathSumIIUtil(root->right, targetSum, prefixSum, count, hash);

        // undoing everything done so far
        hash[prefixSum]--;
        prefixSum -= root->val;
    }

    // Similar to the problem number of subArray with sum k
    int pathSumII(TreeNode* root, int targetSum) {
        long int prefixSum = 0;
        int count = 0;
        unordered_map<long int, int> hash;
        hash[0] = 1; // mane 0 always exists
        pathSumIIUtil(root, targetSum, prefixSum, count, hash);
        return count;
    }

    bool isNodePresent(TreeNode* root, TreeNode* p){
        if(root == nullptr){
            return false;
        }
        else if(root == p){
            return true;
        }
        return isNodePresent(root->left, p) || isNodePresent(root->right, p);
    }

    // Not a binary search tree
    // O(n^2) solution we are trying all the nodes as a potential LCA
    TreeNode* lowestCommonAncestorBinaryTree(TreeNode* root, TreeNode* p, TreeNode* q) {
        if(root == p || root == q){ // p & q are present in the tree and p != q
            return root;
        }
        // have to try both options
        else if((isNodePresent(root->left, p) && isNodePresent(root->right, q)) ||
                (isNodePresent(root->left, q) && isNodePresent(root->right, p))){
            return root;
        }
        // tar mane p and q are in the same side of the root now
        // so root->left is a potential candidate
        else if(isNodePresent(root->left, p)){ //  && isNodePresent(root->left, q) ei call ta
            // kora lagbe na cz we know both are in the same side
            return lowestCommonAncestorBinaryTree(root->left, p, q);
        }
        else{ // root->right is a potential candidate
            return lowestCommonAncestorBinaryTree(root->right, p, q);
        }
    }

    // Optimal solution
    // we can track the entire path parent root to p and parent to q for both and then find the most
    // right common node but that would require extra memory. the below soln uses the returned
    // value(like a true/false) to find the LCA

    TreeNode* lowestCommonAncestorBinaryTreeO(TreeNode* root, TreeNode* p, TreeNode* q){
        if(root == p || root == q){
            return root; // p or q jekono ekta khuje paisi so returning that
        }
        if(root == nullptr){
            return nullptr;
        }

        TreeNode *leftN(nullptr), *rightN(nullptr);

        if(root->left != nullptr){
            leftN = lowestCommonAncestorBinaryTree(root->left, p, q); // what would be the LCA of
            // the binary tree considering root->left as the new root
            // LCA na thakle null return korbe
        }
        if(root->right != nullptr){
            rightN = lowestCommonAncestorBinaryTree(root->right, p, q); // what would be the LCA of
            // the binary tree considering root->right as the new root
            // LCA na thakle null return korbe
        }
        if(leftN != nullptr && rightN != nullptr){
            return root;
        }
        return leftN == nullptr ? rightN : leftN;
    }


    // BFS problem
    // we have to use the property of a tree being stored in an array
    // this we put left child of ith idx in 2*i + 1 and right child of i in 2*i + 2
    int widthOfBinaryTree(TreeNode* root) {
        queue<pair<TreeNode*, int>> q;
        q.push(pair(root, 0));
        int maxWidth = 0;

        while(!q.empty()){
            int levelLen = q.size();
            int l = q.front().second; // left boundary for this level

            TreeNode* curNode;
            int curIdx;
            for(int i = 0; i < levelLen; ++i){
                curNode = q.front().first;
                curIdx = q.front().second;
                q.pop();
                if(curNode->left){ // will not enter nullptr in queues
                    q.push(pair(curNode->left, 2 * curIdx + 1));
                }
                if(curNode->right){
                    q.push(pair(curNode->right, 2 * curIdx + 2));
                }
            }
            // curIdx is actually at the right boundary of our level now

            maxWidth = max(maxWidth, curIdx - l + 1);
        }

        return maxWidth;
    }

    // Very interesting and important problem
    // O(n ^2)
    TreeNode* buildTreeUtil(vector<int>& preorder, vector<int>& inorder, int pl, int pr, int il,
                            int ir) {
        if(pl > pr || il > ir){
            return nullptr;
        }
        TreeNode* root = new TreeNode(preorder[pl]);
        // find the idx of preorder[0] in inorder // O(n)
        ///// Optimize it using a hashmap so that the idx can be found in O(1)
        int idx = -1;
        for(int i = il; i <= ir; ++i){
            if(inorder[i] == preorder[pl]){
                idx = i;
                break;
            }
        }

        // build the left and right subtree
        root->left = buildTreeUtil(preorder, inorder, pl + 1, idx - il + pl, il, idx - 1);
        root->right = buildTreeUtil(preorder, inorder, idx - il + pl + 1, pr, idx + 1, ir);
        return root;
    }


    TreeNode* buildTree(vector<int>& preorder, vector<int>& inorder) {
        return buildTreeUtil(preorder, inorder, 0, preorder.size() - 1, 0, inorder.size() - 1);
    }

    int smallestValBST(TreeNode* root){
        TreeNode* tmp = root;
        while(tmp->left != nullptr){
            tmp = tmp->left;
        }
        return tmp->val;
    }

    int largestValBST(TreeNode* root){
        TreeNode* tmp = root;
        while(tmp->right != nullptr){
            tmp = tmp->right;
        }
        return tmp->val;
    }

    bool isValidBST(TreeNode* root) {
        if(root->left == nullptr && root->right == nullptr){
            return true;
        }
        else if(root->left != nullptr && root->right == nullptr){
            return root->left->val < root->val && largestValBST(root->left) < root->val &&
            isValidBST(root->left);
        }
        else if(root->left == nullptr && root->right != nullptr){
            return root->right->val > root->val && root->val < smallestValBST(root->right) && isValidBST
            (root->right);
        }
        else{
            return root->left->val < root->val && root->right->val > root->val &&
                    largestValBST(root->left) < root->val && root->val < smallestValBST
                    (root->right) && isValidBST(root->left) && isValidBST(root->right);
        }
    }

    bool isValidBSTOptimal(TreeNode* root, long long low, long long high) {
        if(root == nullptr){
            return true;
        }
        else if(root->left == nullptr && root->right == nullptr){
            return low < root->val && root->val < high;
        }
        else if(root->left != nullptr && root->right == nullptr){
            return root->left->val < root->val && isValidBSTOptimal(root->left, low, root->val);
        }
        else if(root->left == nullptr && root->right != nullptr){
            return root->right->val > root->val && isValidBSTOptimal(root->right, root->val, high);
        }
        else{
            return root->left->val < root->val && root->right->val > root->val &&
                    isValidBSTOptimal(root->left, low, root->val) &&
                    isValidBSTOptimal(root->right, root->val, high);
        }
    }

    // Not optimal O(n) solution
    bool isValidBSTOp(TreeNode* root) {
        return isValidBSTOptimal(root, LONG_LONG_MIN, LONG_LONG_MAX);
    }


};

#endif //ABIR_TREEPROBLEMS_H0
