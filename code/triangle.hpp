#ifndef TRIANGLE_HPP
#define TRIANGLE_HPP 23
#include "node.hpp"
#include <cassert>
class Triangle{
public:
  Triangle(int id,Node node0,Node node1,Node node2)
  :id_(id)
  {
    nodes_[0] = node0;
    nodes_[1] = node1;
    nodes_[2] = node2;
  }
  Node& node(int k){
    assert(k>=0&&k<3);
    return nodes_[k];
  }
private:
  int id_;
  Node nodes_[3];
};

#endif
