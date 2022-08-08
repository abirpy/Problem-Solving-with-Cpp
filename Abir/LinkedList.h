//
// Created by Tahsinul Haque Abir on 8/3/22.
//

#ifndef ABIR_LINKEDLIST_H
#define ABIR_LINKEDLIST_H

#include <iostream>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <queue>
#include <stack>
#include <set>

using namespace std;
struct ListNode {
    int val;
    ListNode *next;
    ListNode(int x) : val(x), next(NULL) {}
};

class LinkedList {
public:
    typedef ListNode* ListNodePtr;

    // slow fast
    bool hasCycle(ListNode *head) {
        if(head == nullptr || head->next == nullptr){
            return false;
        }

        ListNodePtr s, f;
        s = f = head;

        while(f != nullptr){
            s = s->next;
            f = f->next;
            if(f != nullptr){
                f = f->next;
            }
            if(s == f){
                return true;
            }
        }
        return false;
    }

    ListNode *detectCycle(ListNode *head) {
        if(head == nullptr || head->next == nullptr){
            return nullptr;
        }
        ListNodePtr s, f;
        s = f = head;

        while(f != nullptr){
            s = s->next;
            f = f->next;
            if(f != nullptr){
                f = f->next;
            }
            if(s == f){
                break;
            }
        }
        if(s != f){
            return nullptr; // no cycle present
        }
        // ekhane asha mane cycle ase and s == f
        s = head;
        while(s != f){
            s = s->next;
            f = f->next;
            if(s == f){ // entry point of cycle
                return s;
            }
        }
        // ekhane s == f and s is the head
        return s;
    }

    // return the node prev to mid
    // if even return first mid
    ListNode* middleNode(ListNode* head) {
        if(head->next == nullptr){
            return head;
        }
        ListNodePtr s, f, prev;
        s = f = head;

        while(f != nullptr && f->next != nullptr){
            prev = s;
            s = s->next;
            f = f->next;
            if(f != nullptr){
                f = f->next;
            }
        }
        return prev;
    }

    // by reference cz egula main function theke pathacci to return later
    void insertFirst(ListNodePtr& head, ListNodePtr& tail, ListNodePtr& node){
        if(head == nullptr){
            head = tail = node;
            tail->next = nullptr;
        }
        else{
            node->next = head;
            head = node;
        }
    }

    // by reference cz egula main function theke pathacci to return later
    void insertLast(ListNodePtr& head, ListNodePtr& tail, ListNodePtr& node){
        if(head == tail && head == nullptr){
            head = tail = node;
        }
        else{
            tail->next = node;
            tail = node;
        }
        tail->next = nullptr;
    }

    ListNodePtr reverseList(ListNodePtr head){
        ListNodePtr newHead(nullptr), newTail(nullptr);

        while(head != nullptr){
            ListNodePtr tmp = head;
            head = head->next;
            insertFirst(newHead, newTail, tmp);
        }
        return newHead;
    }

    // Sapce: O(n) but recursion is cool
    ListNodePtr reverseListRecursive(ListNodePtr head){
        if(head == nullptr || head->next == nullptr){
            return head;
        }
        ListNodePtr newHead = reverseList(head->next);
        head->next->next = head;
        head->next = nullptr;
        return newHead;
    }

    ListNode* reverseBetween(ListNode* head, int left, int right) {
        if(left == right || head->next == nullptr){
            return head;
        }

        int i = 1;
        ListNodePtr newHead, newTail;
        newHead = newTail = nullptr;

        while(i < left && head != nullptr){
            ListNodePtr tmp = head;
            head = head->next;
            insertLast(newHead, newTail, tmp);
            i++;
        }
        // i = left
        ListNodePtr newTail2(nullptr);
        while(i <= right && head != nullptr){
            ListNodePtr tmp = head;
            head = head->next;
            if(newTail != nullptr){
                insertFirst(newTail->next, newTail2, tmp); // newTail->next will act as head of the
            }
            else{
                insertFirst(newHead, newTail2, tmp); // newTail->next will act as head of the
            }

            // reversed part
            i++;
        }

        // now i is after right; don't need i anymore
        // head is after right too
        ListNodePtr newTail3(nullptr);
        while(head != nullptr){
            ListNodePtr tmp = head;
            head = head->next;
            insertLast(newTail2->next, newTail3, tmp);
        }

        return newHead;
    }

    ListNode* rotateRight(ListNode* head, int k) {
        if(head == nullptr || head->next == nullptr){
            return head;
        }
        int len = 1;
        ListNodePtr tmp(head);
        while(tmp->next != nullptr){
            tmp = tmp->next;
            len++;
        }

        if(k % len == 0){
            return head;
        }

        // tmp is at tail now
        tmp->next = head; // made it circular

        ListNodePtr tmp2(head);
        for(int i = 1; i < len - k % len; ++i){
            tmp2 = tmp2->next;
        }
        ListNodePtr tmp3 = tmp2;
        tmp2 = tmp2->next;
        tmp3->next = nullptr;
        return tmp2;
    }


    ListNode* swapPairs(ListNode* head) {
        if(head == nullptr || head->next == nullptr){
            return head;
        }

        ListNodePtr nh, nt, prev;
        nh = nt = prev = nullptr;

        int i = 1;
        bool flag = true;
        while(head != nullptr){
            ListNodePtr tmp(head);
            head = head->next;

            if(i % 2 == 1){
                prev = nt;
                insertLast(nh, nt, tmp);
            }
            else{
                if(flag){
                    insertFirst(nh, nt, tmp);
                    flag = false;
                }
                else{
                    insertFirst(prev->next, nt, tmp);
                }
            }
            i++;
        }
        return nh;
    }

    ListNode* reverseKGroup(ListNode* head, int k) {
        int len = 0;
        ListNodePtr tmp(head);
        while(tmp != nullptr){
            tmp = tmp->next;
            len++;
        }

        int times = len / k; // eto ta sublist reverse kora lagbe
        if(times == 0){
            return head;
        }
        // at least 1 ta sublist reverse kora lagbe
        ListNodePtr nh, nt;
        nh = nt = nullptr;

        int i = 0;
        while(i < k){
            ListNodePtr tmp = head;
            head = head->next;
            insertFirst(nh, nt, tmp);
            i++;
        }
        times--;
        // first sublist reverse shesh
        while(times != 0){
           int i = 0;
           ListNodePtr nt2(nullptr); // a dummy required for the sublists
           while(i < k){
               ListNodePtr tmp = head;
               head = head->next;
               insertFirst(nt->next, nt2, tmp); // nt->next acts as a head for the sublists
               i++;
           }
           nt = nt2; // continuously changing nt our original new tail
           times--;
       }

       // badbaki list if exists
        while(head != nullptr){
            ListNodePtr tmp = head;
            head = head->next;
            insertLast(nh, nt, tmp);
        }
        return nh;
    }


    ListNode* deleteDuplicates(ListNode* head) {
        if(head == nullptr || head->next == nullptr){
            return head;
        }
        ListNodePtr nh, nt;
        nh = nt = nullptr;
        while(head != nullptr){
            ListNodePtr tmp = head;
            head = head->next;

            while(head != nullptr && head->val == tmp->val){ // avoiding same value
                head = head->next;
            }
            insertLast(nh, nt, tmp);
        }
        return nh;
    }

    ListNode* addTwoNumbers(ListNode* l1, ListNode* l2) {
        ListNodePtr nh, nt;
        nh = nt = nullptr;

        ListNodePtr t1(l1), t2(l2);
        int carry = 0;
        while(t1 != nullptr && t2 != nullptr){
            ListNodePtr newNode = new ListNode((t1->val + t2->val + carry) % 10);
            carry = (t1->val + t2->val+ carry) / 10;
            insertLast(nh, nt, newNode);
            t1 = t1->next;
            t2 = t2->next;
        }

        while(t1 != nullptr){ // t2 == nullptr
            ListNodePtr newNode = new ListNode((t1->val + carry) % 10);
            carry = (t1->val + carry) / 10;
            insertLast(nh, nt, newNode);
            t1 = t1->next;
        }

        while(t2 != nullptr){ // t1 == nullptr
            ListNodePtr newNode = new ListNode((t2->val + carry) % 10);
            carry = (t2->val+ carry) / 10;
            insertLast(nh, nt, newNode);
            t2 = t2->next;
        }

        if(carry != 0){
            ListNodePtr newNode = new ListNode(carry);
            insertLast(nh, nt, newNode);
        }

        return nh;
    }

    // merge in O(m + n) with O(1) space
    ListNodePtr mergeSortedLists(ListNodePtr l1, ListNodePtr l2){
        ListNodePtr nh, nt;
        nh = nt = nullptr;
        while(l1 != nullptr && l2 != nullptr){
            if(l1->val < l2->val){
                ListNodePtr tmp1(l1);
                l1 = l1->next;
                insertLast(nh, nt, tmp1);
            }
            else if(l1->val > l2->val){
                ListNodePtr tmp2(l2);
                l2 = l2->next;
                insertLast(nh, nt, tmp2);
            }
            else{
                ListNodePtr tmp1(l1);
                l1 = l1->next;
                insertLast(nh, nt, tmp1);
                ListNodePtr tmp2(l2);
                l2 = l2->next;
                insertLast(nh, nt, tmp2);
            }
        }

        while(l1 != nullptr){
            ListNodePtr tmp1(l1);
            l1 = l1->next;
            insertLast(nh, nt, tmp1);
        }

        while(l2 != nullptr){
            ListNodePtr tmp2(l2);
            l2 = l2->next;
            insertLast(nh, nt, tmp2);
        }
        return nh;
    }

    // mergeSort
    ListNode* sortList(ListNode* head) {
        if(head == nullptr || head->next == nullptr){
            return head;
        }
        ListNodePtr prevMid = middleNode(head);
        ListNodePtr mid = prevMid->next;
        prevMid->next = nullptr;
        ListNodePtr h1 = sortList(head);
        ListNodePtr h2 = sortList(mid);
        return mergeSortedLists(h1, h2);
    }

// The idea is to create two sublist one with odd and one with even and then add them
//    ListNode* oddEvenList(ListNode* head) {
//
//    }

    // The solution with O(1) space is pretty tough
    // First you find the midpoint // actually you need to return the pointer to the node below
    // the midpoint. So that we can separate first ans second list
    // then you reverse the second half
    // then you match the first and reversed second half
//    bool isPalindrome(ListNode* head) {
//
//    }

};


#endif //ABIR_LINKEDLIST_H
