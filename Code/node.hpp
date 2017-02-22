#ifndef NODE_HPP
#define NODE_HPP 23
class Node{
public:
  Node(){};
  Node(int id,double r,double z)
  :id_(id),r_(r),z_(z)
  {};
  Node(Node const& node)
  :id_(node.id_),r_(node.r_),z_(node.z_)
  {};
  int id(){
    return id_;
  }
  void operator= (Node const& n)
  {
    id_= n.id_;
    r_ = n.r_;
    z_ = n.z_;
  }

private:
  int id_;
  double r_;
  double z_;

};
#endif
