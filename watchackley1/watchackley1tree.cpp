/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   watchackley1tree.cpp
 * Author: alusov
 *
 * Created on December 10, 2016, 6:42 PM
 */
#define TESTBBTREE
#include <dirent.h>
#include <string.h>
#include <sstream>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <vector>
#include <bags/bfsbag.hpp>
#include <bags/idgen.hpp>
#include <cutfact/lbcutfact/recordsupp.hpp>
#include <cutfact/lbcutfact/lbcutfactory.hpp>
#include <applycut/basic/serialcutapp.hpp>
#include <decomp/bisectdecomp.hpp>
#include <oneobj/contboxconstr/ackley1.hpp>
#include <solver/treewatchsolver.hpp>
#include <solver/basesolver.hpp>
#include <box/boxutils.hpp>
#include <common/interval.hpp>
#include <api.h>
#include <node.h>

#include "ackley1bndsupp.hpp"
#define EPSILON 0.001

void clear()
{
    //note current folder must be TESTBBTREE. Run system("exec pwd") to check it.
    system("exec rm ./../bbtree/data/*.*");
}

int is_equal(double x, double y) 
{
    return ::abs(x - y) < EPSILON ? 1 : 0;
}

bool stopper(const NUC::BaseSolver<double>& solver) {
    static int cnt = 0;
    const int maxCnt = 20000;
    if (cnt++ > maxCnt) {
        return true;
    } else {
        return false;
    }
}

std::string GetInfo(const NUC::Sub<double>& sub, NUC::RecordSupplier<double>& rs)
{
   std::stringstream ss;
   ss << "Box: " <<  snowgoose::BoxUtils::toString(sub.mBox) << '\n' \
      << "LB: " << sub.mScore << '\n' \
      << "UB: " << rs.getBound(sub.mBox);
   return ss.str();
}

int main() {
    clear();
    
    const int n = 2;

    // Setup problem
    OPTITEST::Ackley1ProblemFactory2 fact(0.9, 1.2, -0.1, 0.2);
    COMPI::MPProblem<double> *mpp = fact.getProblem();

    //Setup bag of sub problems
    NUC::Sub<double> sub(NUC::IdGen::Instance().Id(), 0, std::numeric_limits<double>::max(), *(mpp->mBox));   
    NUC::BFSBag<double> bag;
    bag.putSub(sub);
    //Setup Cut Factory
    NUC::RecordSupplier<double> rs(std::numeric_limits<double>::max());
    COMPI::Functor<double>* pf = mpp->mObjectives.at(0);
    TESTNUC::Ackley1BoundSupplier<double> ibs(n);
    NUC::LBCutFactory<double> cf(EPSILON, rs, ibs);

    // Setup decomposer
    NUC::BisectDecomp<double> bisdec;
    // Setup cut applicator 
    NUC::SerialCutApplicator<double> cutapp;
    // Setup solver
    NUC::TreeWatchSolver<double> solver(bag, bisdec, cf, cutapp);

    // Set stopper for solver
    solver.setStopper(stopper);
    //Setup step watchers
    BBTree::create_node(0, sub.mId, BBTree::NodeState::NewBorn, 1, GetInfo(sub, rs));
    auto wtf = [&](const NUC::Sub<double>& sub,
            const std::vector<std::shared_ptr <NUC::Cut <double> > >& cutv,
            const std::vector< NUC::Sub<double> >& subv,
            const NUC::BaseSolver<double>& slv) {
        if(subv.size())
            BBTree::update_node(sub.mId, sub.mIdEnd, BBTree::NodeState::GiveBirth);
        else
            BBTree::update_node(sub.mId, sub.mIdEnd, BBTree::NodeState::EndUp);
    
        for(int i=0; i < subv.size(); i++)
            BBTree::create_node(sub.mId, subv[i].mId, BBTree::NodeState::NewBorn, i+1, GetInfo(subv[i], rs));       
    };
    solver.addStepWatcher(wtf);

    double x[n];
    //Setup sub evaluators
    auto sf = [&](NUC::Sub<double>& s) {
        snowgoose::BoxUtils::getCenter(s.mBox, x);
        double v = pf->func(x);
        rs.updateRv(v);
        s.mScore = ibs.getBound(s.mBox);
    };
    solver.addSubEval(sf);

    // Run solver
    solver.solve(); 

    std::cout << "Best value found : " << rs.getBound(sub.mBox) << "\n";
    
    SG_ASSERT(is_equal(rs.getBound(sub.mBox), 2.57993));
    
    //BBTree::storage.SaveBuf();
    
    return 0;
}

