// Markus Buchholz

#include <iostream>
#include <vector>
#include <Eigen/Dense>
#include <math.h>
#include <tuple>

int nodes = 3;
int DOF = 2 * nodes;
float factor = std::pow(10, 4);
float inch = 2.56;

//---------------------------------------------------------------------------------------

std::tuple<Eigen::Matrix<float, 6, 6>, std::vector<Eigen::Matrix<float, 4, 1>>> computeStiffness(std::vector<std::pair<float, float>> nodes, std::vector<float> L, std::vector<float> EA, std::vector<std::vector<int>> nx)
{
    // global stiffness
    Eigen::Matrix<float, 6, 6> K = Eigen::Matrix<float, 6, 6>::Zero();
    std::vector<Eigen::Matrix<float, 4, 1>> A;

    std::vector<Eigen::Matrix<float, 4, 4>> stiffnessLocal;
    nodes.push_back(nodes[nodes.size()]);

    for (int ii = 0; ii < nodes.size() - 1; ii++)
    {

        // trasnforamtion matrix from local to global coordinative system
        Eigen::Matrix<float, 4, 1> a;
        a << -(nodes[ii + 1].first - nodes[ii].first) / L[ii],
            -(nodes[ii + 1].second - nodes[ii].second) / L[ii],
            (nodes[ii + 1].first - nodes[ii].first) / L[ii],
            (nodes[ii + 1].second - nodes[ii].second) / L[ii];
        A.push_back(a);
        Eigen::Matrix<float, 1, 4> at = a.transpose();
        Eigen::Matrix<float, 4, 4> k = a * (EA[ii] / L[ii]) * at;

        Eigen::Matrix<float, 6, 6> ki = Eigen::Matrix<float, 6, 6>::Zero();

        // A(row, column)

        for (int i = 0; i < nx[ii].size(); i++)
        {

            for (int j = 0; j < nx[ii].size(); j++)
            {

                ki(nx[ii][i] - 1, nx[ii][j] - 1) = k(i, j);
            }
        }

        K += ki;
    }
    return std::make_tuple(K, A);
}
//---------------------------------------------------------------------------------------

std::tuple<std::vector<float>, Eigen::Matrix<float, 1, 4>> computeInternalForces(Eigen::Matrix<float, 6, 6> K, std::vector<Eigen::Matrix<float, 4, 1>> A, Eigen::Matrix<float, 6, 1> P, std::vector<float> EA, std::vector<float> L, std::vector<std::vector<int>> dofDir)
{

    Eigen::Matrix<float, 2, 2> kff = K.block<2, 2>(0, 0);
    Eigen::Matrix<float, 4, 2> krf = K.block<4, 2>(2, 0);
    Eigen::Matrix<float, 2, 1> Pf = P.block<2, 1>(0, 0);

    std::vector<float> N;

    // displacement
    Eigen::Matrix<float, 6, 1> U = Eigen::Matrix<float, 6, 1>::Zero();

    Eigen::Matrix<float, 2, 1> Uf;

    // displacement Uf (U1, U2)
    Uf = kff.inverse() * Pf;

    U(0, 0) = Uf(0, 0);
    U(1, 0) = Uf(1, 0);

    // Axial fores

    for (int ii = 0; ii < A.size(); ii++)
    {
        Eigen::Matrix<float, 4, 1> u;
        u << U(dofDir[ii][0] - 1, 0), U(dofDir[ii][1] - 1, 0), U(dofDir[ii][2] - 1, 0), U(dofDir[ii][3] - 1, 0);
        float ni = (EA[ii] / L[ii]) * A[ii].dot(u);
        N.push_back(ni);
    }

    // Raction forces
    Eigen::Matrix<float, 1, 4> R = krf * Uf;

    return std::make_tuple(N, R);
}
//---------------------------------------------------------------------------------------

int main()
{

    float p2 = 10;
    float u1 = {0}; //?
    float u2 = {0}; // ?

    // std::vector<std::pair<int, int>> dofDir = {{3, 4}, {1, 5}, {6, 2}};
    std::vector<std::vector<int>> dofDir = {{3, 4, 1, 5}, {1, 5, 6, 2}, {3, 4, 6, 2}};

    // external forces
    Eigen::Matrix<float, 6, 1> P;
    P << 0, p2, 0, 0, 0, 0;

    // postion of nodes
    std::vector<std::pair<float, float>> nodes = {{0, 0}, {120, 0}, {72, 96}};

    // length of each element
    float la = 120;
    float lb = std::sqrt(std::pow(48, 2) + std::pow(96, 2));
    float lc = std::sqrt(std::pow(72, 2) + std::pow(96, 2));
    std::vector<float> L = {la, lb, lc};

    // stiffness of each element
    std::vector<float> EA = {10000, 20000, 30000};

    auto KA = computeStiffness(nodes, L, EA, dofDir);

    auto NR = computeInternalForces(std::get<0>(KA), std::get<1>(KA), P, EA, L, dofDir);

    std::cout << "------axial forces------"
              << "\n";
    for (auto &ni : std::get<0>(NR))
    {
        std::cout << ni << "\n";
    }
    std::cout << "------reaction forces------"
              << "\n";
    std::cout << std::get<1>(NR) << "\n";

    return 0;
}
