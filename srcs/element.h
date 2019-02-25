#ifndef HEADER_ELEMENT
#define HEADER_ELEMENT

#include <vector>


class element{

    public:
    element();
    element(int Tag, int numberNodes, std::vector<int> nodes);
    ~element();
    int m_Tag;
    int m_numberNodes;
    std::vector<int> m_nodes;
};

element::element(){

    m_Tag = 0;
    m_numberNodes = 0;
    m_nodes = {};
}

element::element(int Tag, int numberNodes, std::vector<int> nodes){

    m_Tag = Tag;
    m_numberNodes = numberNodes;
    m_nodes = nodes;
}

element::~element(){
}

#endif