#include <cstdio>
#include <gmsh.h>
#include <errno.h>
#include <iostream>
#include <utils.h>

#include "Mesh.h"
#include "solver.h"
#include "configParser.h"

// This is test unit with a very simple test case
namespace factory = gmsh::model::geo;

bool isEqual(std::vector<double> &a, std::vector<double> &b){
    return std::equal(a.begin(), a.end(), b.begin(),
               [](double value1, double value2)
               {
                   constexpr double epsilon = 0.001;
                   return std::fabs(value1 - value2) < epsilon;
               });
}

int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 1);

    // Create a rectangle [1,1] with 2 elements
    gmsh::model::addDiscreteEntity(2, 1);
    gmsh::model::mesh::setNodes(2, 1, {1, 2, 3, 4},
                                {0., 0., 0.,
                                 1., 0., 0.,
                                 1., 1., 0.,
                                 0., 1., 0.});
    gmsh::model::mesh::setElements(2, 1, {2}, {{1, 2}},
                                   {{1, 2, 3, 1, 3, 4}});
    gmsh::model::mesh::generate(2);
    gmsh::write("testUnit.msh");

    // Init test
    Config config = config::parseConfig("doc/config/config.conf");
    Mesh mesh("test", config);

    std::cout << "------------------------------------" << std::endl;
    std::cout << "          Begin Unit Test" << std::endl;
    std::cout << "------------------------------------" << std::endl;

    // Normals
    std::vector<double> e0n0_exact = {0, -1, 0};
    std::vector<double> e0n1_exact = {1, 0, 0};
    std::vector<double> e0n2_exact = {-sqrt(2) / 2., sqrt(2) / 2., 0};
    std::vector<double> e1n0_exact = {sqrt(2) / 2., -sqrt(2) / 2., 0};
    std::vector<double> e1n1_exact = {0, 1, 0};
    std::vector<double> e1n2_exact = {-1, 0, 0};
    std::vector<double> e0n0(&mesh.fNormal(mesh.elFId(0, 0)), &mesh.fNormal(mesh.elFId(0, 0)) + 3);
    std::vector<double> e0n1(&mesh.fNormal(mesh.elFId(0, 1)), &mesh.fNormal(mesh.elFId(0, 1)) + 3);
    std::vector<double> e0n2(&mesh.fNormal(mesh.elFId(0, 2)), &mesh.fNormal(mesh.elFId(0, 2)) + 3);
    std::vector<double> e1n0(&mesh.fNormal(mesh.elFId(1, 0)), &mesh.fNormal(mesh.elFId(1, 0)) + 3);
    std::vector<double> e1n1(&mesh.fNormal(mesh.elFId(1, 1)), &mesh.fNormal(mesh.elFId(1, 1)) + 3);
    std::vector<double> e1n2(&mesh.fNormal(mesh.elFId(1, 2)), &mesh.fNormal(mesh.elFId(1, 2)) + 3);
    std::transform(e0n0.begin(), e0n0.end(), e0n0.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, mesh.elFOrientation(0, 0)));
    std::transform(e0n1.begin(), e0n1.end(), e0n1.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, mesh.elFOrientation(0, 1)));
    std::transform(e0n2.begin(), e0n2.end(), e0n2.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, mesh.elFOrientation(0, 2)));
    std::transform(e1n0.begin(), e1n0.end(), e1n0.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, mesh.elFOrientation(1, 0)));
    std::transform(e1n1.begin(), e1n1.end(), e1n1.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, mesh.elFOrientation(1, 1)));
    std::transform(e1n2.begin(), e1n2.end(), e1n2.begin(),
                   std::bind(std::multiplies<double>(), std::placeholders::_1, mesh.elFOrientation(1, 2)));
    assert(isEqual(e0n0, e0n0_exact));
    assert(isEqual(e0n1, e0n1_exact));
    assert(isEqual(e0n2, e0n2_exact));
    assert(isEqual(e1n0, e1n0_exact));
    assert(isEqual(e1n1, e1n1_exact));
    assert(isEqual(e1n2, e1n2_exact));
    std::cout << "[1] All normals are correctly calculated and oriented" << std::endl;

    // Mass matrix
    std::vector<int> eleTypes;
    gmsh::model::mesh::getElementTypes(eleTypes, 2);

    int eleType2D = eleTypes[0];
    std::string name;
    int dim, order, numNodes;
    std::vector<double> paramCoord;
    gmsh::model::mesh::getElementProperties(eleType2D, name, dim, order, numNodes, paramCoord);

    std::vector<int> elementTags, nodeTags;
    gmsh::model::mesh::getElementsByType(eleType2D, elementTags, nodeTags, -1);

    std::vector<double> intpts, bf;
    int numComp;
    gmsh::model::mesh::getBasisFunctions(eleType2D, "Gauss3", "Lagrange", intpts, numComp, bf);
    int numIntPts = (int) intpts.size() / 4.;

    std::vector<double> jac, det, pts;
    gmsh::model::mesh::getJacobians(eleType2D, "Gauss3", jac, det, pts);

    std::vector<double> e1MassMatrix_exact(numNodes * numNodes, 0);
    std::vector<double> e2MassMatrix_exact(numNodes * numNodes, 0);
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            for (int g = 0; g < numIntPts; g++) {
                e1MassMatrix_exact[i * numNodes + j] +=
                        intpts[g * 4 + 3] * bf[g * numNodes + i] * bf[g * numNodes + j] * det[0 * numIntPts + g];
                e2MassMatrix_exact[i * numNodes + j] +=
                        intpts[g * 4 + 3] * bf[g * numNodes + i] * bf[g * numNodes + j] * det[1 * numIntPts + g];
            }
        }
    }

    std::vector<double> e1MassMatrix(numNodes * numNodes, 0);
    std::vector<double> e2MassMatrix(numNodes * numNodes, 0);
    mesh.getElMassMatrix(0, false, e1MassMatrix.data());
    mesh.getElMassMatrix(1, false, e2MassMatrix.data());

    assert(isEqual(e1MassMatrix, e1MassMatrix_exact));
    assert(isEqual(e2MassMatrix, e2MassMatrix_exact));

    std::cout << "[2] All mass matrices have been correctly calculated" << std::endl;

    // Derivative of shape function along physical direction
    std::vector<double> intptsGrad, Gradbf;
    gmsh::model::mesh::getBasisFunctions(eleType2D, "Gauss3", "GradLagrange", intptsGrad, numComp, Gradbf);

    int N = 3;
    int el = 0;
    std::vector<double> eljac(N * N);
    for (int g = 0; g < numIntPts; ++g) {
        for (int f = 0; f < numNodes; f++) {
            std::copy(&jac[el * 9 * numIntPts + g * 9], &jac[el * 9 * numIntPts + g * 9] + 9, eljac.begin());
            lapack::solve(eljac.data(), &Gradbf[g * 3 * numNodes + f * 3], N);
        }
    }

    std::vector<double> a = {1, 0, 0};
    std::vector<double> u(nodeTags.size() * numNodes, 1);

    std::vector<double> e1StiffVector_exact(numNodes, 0);
    std::vector<double> e2StiffVector_exact(numNodes, 0);
    for (int i = 0; i < numNodes; i++) {
        for (int j = 0; j < numNodes; j++) {
            for (int g = 0; g < numIntPts; g++) {
                for (int x = 0; x < 3; ++x) {
                    e1StiffVector_exact[i] +=
                            intptsGrad[g * 4 + 3] * Gradbf[g * numNodes * 3 + i * 3 + x] * a[x] * bf[g * numNodes + j] *
                            u[0 * numNodes + j] * det[0 * numIntPts + g];
                }
            }
        }
    }

    std::vector<double> e1StiffVector(numNodes, 0);
    std::vector<double> e2StiffVector(numNodes, 0);
    mesh.getElStiffVector(0, a.data(), u.data(), e1StiffVector.data());
    mesh.getElStiffVector(1, a.data(), u.data(), e2StiffVector.data());

    assert(isEqual(e1StiffVector, e1StiffVector_exact));

    std::cout << "[3] Stiffness vector correctly calculated." << std::endl;

    // Flux
    std::vector<int> nodes;
    gmsh::model::mesh::getElementEdgeNodes(eleType2D, nodes);
    int c = gmsh::model::addDiscreteEntity(1);
    int eleType1D = gmsh::model::mesh::getElementType("line", order);
    gmsh::model::mesh::setElementsByType(1, c, eleType1D, {}, nodes);
    std::vector<double> fintpts, fbf;
    gmsh::model::mesh::getBasisFunctions(eleType1D, "Gauss4", "Lagrange", fintpts, numComp, fbf);
    std::vector<double> fjac, fdet, fpts;
    gmsh::model::mesh::getJacobians(eleType1D, "Gauss4", fjac, fdet, fpts, c);

    {
        std::vector<double> flux(2);
        std::vector<double> a = {1, 0, 0};
        std::vector<double> u = {1,1,1,2,2,2};
        mesh.getFlux(2, a.data(), u.data(), flux.data());
        std::vector<double> flux_exact ={-0.75,-0.75};
        assert(isEqual(flux, flux_exact));
    }

    std::cout << "[4] Flux correctly computed." << std::endl;

    gmsh::finalize();

    return EXIT_SUCCESS;

    // Small code to vizualize the normals.
    /*std::vector<double> viewNormals;
    for(int f=0; f<m_fNum; f++) {
    std::vector<double> fNodeCoord, fNodeCoord2, paramCoords;
    gmsh::model::mesh::getNode(fNodeTag(f), fNodeCoord, paramCoords);
    gmsh::model::mesh::getNode(fNodeTag(f, 1), fNodeCoord2, paramCoords);
    for(int x=0; x<3; ++x)
    viewNormals.push_back((fNodeCoord[x]+fNodeCoord2[x])/2);
    for(int x=0; x<3; ++x)
    viewNormals.push_back(fNormal(f, x));
    }
    int normalTag;
    gmsh::view::add("normals", normalTag);
    gmsh::view::addListData(normalTag, "VP", m_elNum, viewNormals);
    gmsh::view::write(normalTag, "normal.pos");*/

}
