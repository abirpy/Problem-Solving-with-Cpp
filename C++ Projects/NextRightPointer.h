//
// Created by Tahsinul Haque Abir on 7/3/22.
//

#ifndef ABIR_NEXTRIGHTPOINTER_H
#define ABIR_NEXTRIGHTPOINTER_H

#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <algorithm>

using namespace std;

class Node {
public:
    int val;
    Node* left;
    Node* right;
    Node* next;

    Node() : val(0), left(NULL), right(NULL), next(NULL) {}

    Node(int _val) : val(_val), left(NULL), right(NULL), next(NULL) {}

    Node(int _val, Node* _left, Node* _right, Node* _next)
            : val(_val), left(_left), right(_right), next(_next) {}
};

class NextRightPointer {
    // Traditional BFS solution
    // Time O(V) and Space O(V)
    Node* connect(Node* root) {
        if(root == nullptr){
            return root;
        }

        queue<Node*> q;
        // push it to queue
        q.push(root);

        while(!q.empty()) {
            int qLen = q.size(); // qLen is equal to the numNodes in a particular level
            for(int i = 0; i < qLen; ++i){
                Node* curNode = q.front();
                q.pop(); // popping curNode
                Node* nextNode = i == qLen - 1 ? nullptr : q.front();
                curNode->next = nextNode;

                // Pushing next level nodes
                if(curNode->left != nullptr){
                    q.push(curNode->left);
                }
                if(curNode->right != nullptr){
                    q.push(curNode->right);
                }
            }
        }
        return root;
    }

    // Optimazed solution doesn't use a queue
    // time: O(V) Space: O(1)
    // Initially all next are set to nullptr
    Node* connectOpt(Node* root){
        if(root == nullptr){
            return root;
        }

        Node *leftMost = root;
        while(leftMost->left != nullptr){
            Node* head = leftMost;
            while(head != nullptr){
                // C1
                head->left->next = head->right;

                //C2
                if(head->next != nullptr){
                    head->right->next = head->next->left;
                }
                head = head->next;
            }
            leftMost = leftMost->left;
        }
        return root;
    }

    // Optimazed solution doesn't use a queue
    // time: O(V) Space: O(1)
    // Initially all next are set to nullptr
    // For any size tree not only perfect trees
    Node* connectOptII(Node* root){
        if(root == nullptr){
            return root;
        }

        // Very cringe solution
        Node *parent(root), *child(nullptr), *childHead(nullptr);

        while(parent != nullptr){// Outer loop parent top to bottom
            while(parent != nullptr){// Inner loop parent left to right
                if(parent->left != nullptr){
                    if(childHead == nullptr){
                        childHead = parent->left;
                    }
                    else{
                        child->next = parent->left;
                    }
                    child = parent->left;
                }
                if(parent->right != nullptr){
                    if(childHead == nullptr){
                        childHead = parent->right;
                    }
                    else{
                        child->next = parent->right;
                    }
                    child = parent->right;
                }
                parent = parent->next;
            }
            parent = childHead;
            child = childHead = nullptr;
        }
    }
};


#endif //ABIR_NEXTRIGHTPOINTER_H
