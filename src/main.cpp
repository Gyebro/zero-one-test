#include <iostream>
#include <vector>
#include "microchaos.h"
#include "gnuplot_i.h"

double correlationCoefficient(std::vector<double>& X, std::vector<double>& Y) {
    if (X.size() != Y.size()) {
        std::cout << "Error: different input vector lengths during correlation calculation\n";
    }
    double sum_X = 0, sum_Y = 0, sum_XY = 0;
    double squareSum_X = 0, squareSum_Y = 0;
    size_t N = X.size();
    for (size_t i = 0; i < N; i++) {
        sum_X += + X[i];
        sum_Y += + Y[i];
        sum_XY += + X[i] * Y[i];
        squareSum_X += X[i] * X[i];
        squareSum_Y += Y[i] * Y[i];
    }
    double n = static_cast<double>(N);
    // Correlation coefficient.
    double corr = (n * sum_XY - sum_X * sum_Y) / sqrt((n * squareSum_X - sum_X * sum_X)*(n * squareSum_Y - sum_Y * sum_Y));
    return corr;
}

double average(const std::vector<double>& vec) {
    double avg = 0;
    for (const double& d : vec) avg += d;
    avg /= static_cast<double>(vec.size());
    return avg;
}

double zeroOneTest(double c, size_t idx, size_t steps, const MicroChaosMapStaticRoundIn& system, bool show_figure = false) {
    vec2 y = {10,10};
    vec2 y1;
    std::vector<double> tra;
    for (size_t k = 0; k<steps; k++) {
        y1 = system.step(y);
        if (y1[idx] != y1[idx]) {
            std::cout << "Error: NaN or InF occured!\n";
            break;
        }
        y = y1; // Advance
        tra.push_back(y[idx]); // TODO: Instead of taking one element of the state vector, use other representations, like norm
    }
    std::vector<double> p = tra;
    std::vector<double> q = tra;
    size_t NN = tra.size();
    for (size_t n=0; n<NN-1; n++) {
        p[n+1] = p[n] + tra[n]*cos(c*static_cast<double>(n+1));
        q[n+1] = q[n] + tra[n]*sin(c*static_cast<double>(n+1));
    }
    size_t nd = NN/10;
    // diffusive behaviour
    double EPhi = average(tra);
    if (EPhi != EPhi) {
        return 0;
    }
    std::vector<double> Mc;
    std::vector<double> Vosc;
    std::vector<double> Dc;
    std::vector<double> Idx;
    Mc.resize(nd);
    Vosc.resize(nd);
    Dc.resize(nd);
    Idx.resize(nd);
    for (size_t n=0; n<nd; n++) {
        double mcn = 0;
        for (size_t j=0; j<NN-n; j++) {
            mcn += (pow(p[j+n] - p[j],2.0) + pow(q[j+n] - q[j],2.0));
        }
        Mc[n] = mcn/static_cast<double>(NN-n);
        Vosc[n] = EPhi*EPhi*(1.0 - cos(c*static_cast<double>(n+1)))/(1.0 - cos(c));
        if (Vosc[n] != Vosc[n]) {
            std::cout << "Error, Vosc is NAN\n";
        }
        Dc[n] = Mc[n]-Vosc[n];
        Idx[n] = static_cast<double>(n+1);
    }
    if (show_figure) {
        Gnuplot gp;
        gp.cmd("set terminal png size 1200,800");
        gp.cmd("set output 'test_for_chaos_snapshot_c_"+std::to_string(c)+"_Mc.png'");
        gp.set_pointsize(0.5);
        gp.plot_xy(Idx, Mc, "Mc(n)");
        gp.cmd("set output 'test_for_chaos_snapshot_c_"+std::to_string(c)+"_Dc.png'");
        gp.set_pointsize(0.5);
        gp.plot_xy(Idx, Dc, "Dc(n)");
    }
    // In the basic version, K is the slope of the line fitted on Log(Mc)/Lon(Idx)
    // In this version Kc is the correlation coeff between the index vector and Dc
    double Kc = correlationCoefficient(Idx, Dc);
    return Kc;
}


int main() {

    Gnuplot gp;

    double alpha = 0.06;
    double delta = 0;

    // System and it's control parameter domain
    MicroChaosParamDomain dom(alpha, delta);
    MicroChaosMapStaticRoundIn system(0, 0, alpha, delta);

    system.setP(dom.getBestP());
    system.setD(dom.getBestD());

    std::cout << "Using control parameters [P, D] = " << system.getP() << ", " << system.getD() <<std::endl;

    std::vector<double> X;
    std::vector<double> corrList;

    size_t take_snapshot_idx = 117; // Save a sample plot of Kc versus c

    // Execute zero-one test (for different c values) in parallel
#pragma omp parallel
    {
        std::vector<double> corrList_private;
        std::vector<double> X_private;
#pragma omp for nowait // fill _private arrays in parallel
        for (size_t C = 1; C<500; C++) {
            double c = M_PI*static_cast<double>(C)/500.0;
            double corr = zeroOneTest(c, 1, 10000, system, take_snapshot_idx == C);
            corrList_private.push_back(corr);
            X_private.push_back(c);
        }
#pragma omp critical // collect values
        corrList.insert(corrList.end(), corrList_private.begin(), corrList_private.end());
        X.insert(X.end(), X_private.begin(), X_private.end());
    }

    double avgKc = average(corrList);
    std::cout << "Alpha = " << alpha << ", avg(Kc) = " << avgKc << std::endl;
    gp.cmd("set terminal png size 1200,800");
    gp.cmd("set output 'test_for_chaos_alpha_"+std::to_string(alpha)+"_Kc-c_multithread.png'");
    gp.set_pointsize(0.5);
    gp.plot_xy(X, corrList, "Kc - c");
    gp.reset_all();


    // Here's the non-parallel version for reference:
    X.resize(0);
    corrList.resize(0);
    for (size_t C = 1; C<500; C++) {
        double c = M_PI*static_cast<double>(C)/500.0;
        double corr = zeroOneTest(c, 1, 10000, system, take_snapshot_idx == C);
        corrList.push_back(corr);
        X.push_back(c);
    }

    avgKc = average(corrList);
    std::cout << "Alpha = " << alpha << ", avg(Kc) = " << avgKc << std::endl;
    gp.cmd("set terminal png size 1200,800");
    gp.cmd("set output 'test_for_chaos_alpha_"+std::to_string(alpha)+"_Kc-c_singlethread.png'");
    gp.set_pointsize(0.5);
    gp.plot_xy(X, corrList, "Kc - c");
    gp.reset_all();

    std::cout << "Press [ENTER] to clean up and exit...\n";
    getchar();
    gp.remove_tmpfiles();

    return 0;
}
