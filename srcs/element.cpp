#include "element.h"

element :: element(){

    m_numberNodes = 0;
    m_nodes = {};
}

element :: element(int numberNodes, std::vector<int> nodes){

    m_numberNodes = numberNodes;
    m_nodes = nodes;
}