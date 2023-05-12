#include <vector>

#include "eigen-3.4.0/Eigen/Dense"

#include "solver.h"

class Optimizer {
    public:
        bool run();

        Optimizer(Molecule &mol, int solver_type, double solver_threshold);
        ~Optimizer();
    protected:
        // objective function
        virtual double objective();

        // callback function
        virtual void callback();

        // function to print results
        virtual void results();

        // gradient
        Eigen::VectorXd g;

        // optimization parameters
        std::vector<double *> para;

        // steepest_descent algorithm
        bool steepest_descent();
        void grad();
        double stepsize;

        Molecule &mol;

        // ability to dynamically choose correct solver
        RHFSolver rhf;
        UHFSolver uhf;
        Solver *solver;
};

class GeomOptimizer : public Optimizer {
    public:
        GeomOptimizer(Molecule &mol, int solver_type, double solver_threshold);
        ~GeomOptimizer();
    protected:
        // override objective function
        double objective();

        void callback();

        void results();
};

class BasisOptimizer : public Optimizer {
    public:
        BasisOptimizer(Molecule &mol, int solver_type, double solver_threshold);
        ~BasisOptimizer();
    protected:
        // override objective function
        double objective();

        void callback();

        void results();
};