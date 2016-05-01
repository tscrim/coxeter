/**
 * Coxeter test suite.
 *
 * Coxeter version 3.1  Copyright (C) 2015 Travis Scrimshaw
 * See file main.cpp for full copyright notice
 */

#include <iostream>
#include <stdio.h>
#include <time.h>

#include "tests.h"
#include "coxgroup.h"
#include "coxtypes.h"
#include "interactive.h"

namespace coxeter {

using namespace std;

void testA1() {
  cout<<"Testing A1"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("A"), (Rank) 1);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testA2() {
  cout<<"Testing A2"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("A"), (Rank) 2);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testA3() {
  cout<<"Testing A3"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("A"), (Rank) 3);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testB3() {
  cout<<"Testing B3"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("B"), (Rank) 3);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testH3() {
  cout<<"Testing H3"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("H"), (Rank) 3);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testD4() {
  cout<<"Testing D4"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("D"), (Rank) 4);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testE6() {
  cout<<"Testing E6"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("E"), (Rank) 6);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testG2() {
  cout<<"Testing G2"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("G"), (Rank) 2);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testAffineA2() {
  cout<<"Testing a2"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("a"), (Rank) 2);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testAffineC3() {
  cout<<"Testing c3"<<endl;
  clock_t start = clock();
  CoxGroup* W = interactive::coxeterGroup(Type("c"), (Rank) 3);
  delete W;
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testAffineD4() {
  cout<<"Testing d4"<<endl;
  CoxGroup* W = interactive::coxeterGroup(Type("d"), (Rank) 4);
  delete W;
  clock_t start = clock();
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testStar3() {
  cout<<"Testing star3,3,3"<<endl;
  clock_t start = clock();
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testHowlett1() {
  cout<<"Testing Howlett1"<<endl;
  clock_t start = clock();
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testPerformace() {
  cout<<"Testing performance"<<endl;
  clock_t start = clock();
  printf("Time taken: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

void testAll() {
  clock_t start = clock();
  testA1();
  testA2();
  testA3();
  testB3();
  testD4();
  testE6();
  testG2();
  testAffineA2();
  testAffineC3();
  testAffineD4();
  testStar3();
  testHowlett1();
  testPerformace();
  printf("Total time for all tests: %.2fs\n", (double)(clock() - start)/CLOCKS_PER_SEC);
}

} // namespace coxeter
