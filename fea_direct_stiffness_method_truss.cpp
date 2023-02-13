// Markus Buchhholz

#include <iostream>
#include <vector>
#include <math.h>
#include <tuple>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <algorithm>

double E = 10000; // modulus of elasticity
double A = 0.111; //cross-section of the bar
double EA = E * A;

//---------------------------------------------------------------------------------------

std::tuple<Eigen::Matrix<double, 16, 16>, std::vector<Eigen::Matrix<double, 4, 1>>, std::vector<double>, Eigen::Matrix<double, 16, 16>> structuralAnalysis(std::vector<std::pair<double, double>> nodes, std::vector<std::pair<int, int>> bars, std::vector<std::vector<int>> nx)
{
    // length
    std::vector<double> L;
    // transformation from local to global coordinative system
    std::vector<Eigen::Matrix<double, 4, 1>> A;
    // NDOF
    int DOF = 2;
    constexpr int NDOF = 16; 

    // global stiffness
    Eigen::Matrix<double, 16, 16> K = Eigen::Matrix<double, 16, 16>::Zero();
    Eigen::Matrix<double, 16, 16> MM = Eigen::Matrix<double, 16, 16>::Zero();
    Eigen::Matrix<double, 4, 4> I;
    I << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;

    for (int ii = 0; ii < bars.size(); ii++)
    {
        Eigen::Matrix<double, 4, 1> a;

        double x0 = nodes[bars[ii].first].first;
        double y0 = nodes[bars[ii].first].second;

        double x1 = nodes[bars[ii].second].first;
        double y1 = nodes[bars[ii].second].second;

        double li = std::sqrt(std::pow(x1 - x0, 2) + std::pow(y1 - y0, 2));
        L.push_back(li);

        a << -(x1 - x0) / li,
            -(y1 - y0) / li,
            (x1 - x0) / li,
            (y1 - y0) / li;

        A.push_back(a);

        // global stiffness

        Eigen::Matrix<double, 16, 16> ki = Eigen::Matrix<double, 16, 16>::Zero();
        Eigen::Matrix<double, 1, 4> at = a.transpose();
        Eigen::Matrix<double, 4, 4> k = a * (EA / li) * at;

        // A(row, column)

        for (int i = 0; i < nx[ii].size(); i++)
        {

            for (int j = 0; j < nx[ii].size(); j++)
            {

                ki(nx[ii][i] - 0, nx[ii][j] - 0) = k(i, j);
            }
        }
        K += ki;

        // natural frequency
        double DN = 0.1;
        double Aa = 0.111;
        Eigen::Matrix<double, 16, 16> mi = Eigen::Matrix<double, 16, 16>::Zero();

        Eigen::Matrix<double, 4, 4> m = DN * Aa * li * (0.001 / 386.09) * 0.5 * I;

        for (int i = 0; i < nx[ii].size(); i++)
        {

            for (int j = 0; j < nx[ii].size(); j++)
            {

                mi(nx[ii][i] - 0, nx[ii][j] - 0) = m(i, j);
            }
        }
        MM += mi;
    }


    return std::make_tuple(K, A, L, MM);
}
// https://stackoverflow.com/questions/13290395/how-to-remove-a-certain-row-or-column-while-using-eigen-library-c
//---------------------------------------------------------------------------------------
void removeRow(Eigen::MatrixXd &matrix, unsigned int rowToRemove)
{
    unsigned int numRows = matrix.rows() - 1;
    unsigned int numCols = matrix.cols();

    if (rowToRemove < numRows)
        matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

    matrix.conservativeResize(numRows, numCols);
}
//---------------------------------------------------------------------------------------

void removeColumn(Eigen::MatrixXd &matrix, unsigned int colToRemove)
{
    unsigned int numRows = matrix.rows();
    unsigned int numCols = matrix.cols() - 1;

    if (colToRemove < numCols)
        matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

    matrix.conservativeResize(numRows, numCols);
}

//---------------------------------------------------------------------------------------

std::tuple<std::vector<double>, Eigen::Matrix<double, 1, 4>, Eigen::Matrix<double, 16, 1>, std::vector<double>> computeInternalForces(Eigen::Matrix<double, 16, 16> Ki, std::vector<Eigen::Matrix<double, 4, 1>> A, Eigen::Matrix<double, 16, 1> P, std::vector<double> L, std::vector<std::vector<int>> dofDir, Eigen::Matrix<double, 16, 16> MMi)
{


    Eigen::MatrixXd K(16, 16);
    Eigen::MatrixXd K2(16, 16);
    Eigen::MatrixXd MM(16, 16);
    K = Ki;
    K2 = Ki;
    MM = MMi;

    // slicing for kff
    removeColumn(K, 8);
    removeColumn(K, 8);
    removeRow(K, 8);
    removeRow(K, 8);
    removeColumn(K, 0);
    removeColumn(K, 0);
    removeRow(K, 0);
    removeRow(K, 0);

    // slicing for krf
    removeColumn(K2, 8);
    removeColumn(K2, 8);
    for (int ii = 0; ii < 6; ii++)
    {

        removeRow(K2, 11);
    }
    removeColumn(K2, 0);
    removeColumn(K2, 0);

    for (int ii = 0; ii < 6; ii++)
    {

        removeRow(K2, 2);
    }

    // slicing for mff
    removeColumn(MM, 8);
    removeColumn(MM, 8);
    removeRow(MM, 8);
    removeRow(MM, 8);
    removeColumn(MM, 0);
    removeColumn(MM, 0);
    removeRow(MM, 0);
    removeRow(MM, 0);

    Eigen::Matrix<double, 12, 12> kff = K;
    Eigen::Matrix<double, 4, 12> krf = K2;

    //-------Natural frequency-------------

    Eigen::Matrix<double, 12, 12> mff = MM;

    Eigen::MatrixXd matrix;
    matrix = mff.inverse() * kff;
    Eigen::EigenSolver<Eigen::MatrixXd> es(matrix);
    Eigen::MatrixXcd eigValeues = es.eigenvalues();

    std::vector<double> nfreq;

    for (int ii = 0; ii < eigValeues.rows(); ii++)
    {

        double w = std::sqrt((real(eigValeues(ii))));
        nfreq.push_back(w / (2 * M_PI));
    }
    std::sort(nfreq.begin(), nfreq.end());

    //-------Stress forces-------------

    Eigen::Matrix<double, 12, 1> Pf = P.block<12, 1>(4, 0);
    std::vector<double> N;

    // displacement
    Eigen::Matrix<double, 16, 1> U;
    Eigen::Matrix<double, 12, 1> Uf;
    Uf = kff.inverse() * Pf;
    U << 0, 0, Uf(0, 0), Uf(1, 0), Uf(2, 0), Uf(3, 0), Uf(4, 0), Uf(5, 0), 0, 0, Uf(6, 0), Uf(7, 0), Uf(8, 0), Uf(9, 0), Uf(10, 0), Uf(11, 0);

    // Axial fores

    for (int ii = 0; ii < A.size(); ii++)
    {
        Eigen::Matrix<double, 1, 4> u;
        u << U(dofDir[ii][0] - 0, 0), U(dofDir[ii][1] - 0, 0), U(dofDir[ii][2] - 0, 0), U(dofDir[ii][3] - 0, 0);
        double ni = (EA / L[ii]) * A[ii].dot(u);
        N.push_back(ni);
    }

    // Raction forces
    Eigen::Matrix<double, 1, 4> R = krf * Uf;

    return std::make_tuple(N, R, U, nfreq);
}

//---------------------------------------------------------------------------------------

int main()
{

    std::vector<std::pair<double, double>> nodes = {{0, 120}, {120, 120}, {240, 120}, {360, 120}, {0, 0}, {120, 0}, {240, 0}, {360, 0}};
    std::vector<std::pair<int, int>> bars = {{0, 1}, {1, 2}, {2, 3}, {4, 5}, {5, 6}, {6, 7}, {5, 1}, {6, 2}, {7, 3}, {0, 5}, {4, 1}, {1, 6}, {5, 2}, {2, 7}, {6, 3}};

    std::vector<std::vector<int>> dofDir = {{0, 1, 2, 3}, {2, 3, 4, 5}, {4, 5, 6, 7}, {8, 9, 10, 11}, {10, 11, 12, 13}, {12, 13, 14, 15}, {10, 11, 2, 3}, {12, 13, 4, 5}, {14, 15, 6, 7}, {0, 1, 10, 11}, {8, 9, 2, 3}, {2, 3, 12, 13}, {10, 11, 4, 5}, {4, 5, 14, 15}, {12, 13, 6, 7}};

    // external forces
    double load = -10;
    Eigen::Matrix<double, 16, 1> P = Eigen::Matrix<double, 16, 1>::Zero();
    P(15, 0) = load;

    // support displacement
    Eigen::Matrix<double, 1, 4> Ur = Eigen::Matrix<double, 1, 4>::Zero();

    // conditions of DOF
    Eigen::Matrix<int, 8, 2> DC = Eigen::Matrix<int, 8, 2>::Ones();
    DC(0, 0) = DC(0, 1) = 0; // support does not move
    DC(4, 0) = DC(4, 1) = 0; // support does not move

    auto KALMM = structuralAnalysis(nodes, bars, dofDir);

    auto NRUF = computeInternalForces(std::get<0>(KALMM), std::get<1>(KALMM), P, std::get<2>(KALMM), dofDir, std::get<3>(KALMM));

    std::cout << "------axial forces------"
              << "\n";
    for (auto &n : std::get<0>(NRUF))
    {
        std::cout << n << "\n";
    }

    std::cout << "------reaction forces------"
              << "\n";
    std::cout << std::get<1>(NRUF).transpose() << "\n";

    std::cout << "------displacement------"
              << "\n";
    for (int ii = 0; ii < std::get<2>(NRUF).rows() - 1; ii = ii + 2)
    {

        std::cout << "x= " << std::get<2>(NRUF)(ii) << ", y= " << std::get<2>(NRUF)(ii + 1) << "\n";
    }

    std::cout << "------natural frequencies------"
            << "\n";
    for (auto &f : std::get<3>(NRUF))
    {
        std::cout << f << "\n";
    }
}
