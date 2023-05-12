#include <vector>
#include <map>
#include <set>

#include "eigen-3.4.0/Eigen/Dense"

#include "molecule.h"
#include "basis.h"

#pragma once

class Solver {
    public:
        // function to run the respective SCF procedure
        virtual double run(bool verbose);

        // pass pointers to basis function zeta values for optimizer
        std::vector<double *> get_basis();

        // reinitialize basis functions
        void reinit();

        Solver(const Molecule &mol, double threshold);
        ~Solver();

    protected:
        const Molecule &mol;

        // variables for SCF
        Eigen::MatrixXd hcore, s, ortho, s_diag, s_evecs;
        Eigen::VectorXd g_abcd;
        std::map<std::set<std::set<int>>, int> abcd_map;
        double threshold;

        // storage for basis functions
        std::vector<BasisFunction *> basis;
        int nbas;

        // calculates the nuclear-nuclear interactions
        void nuclear();
        double v_nn;

        // calculates 1e integrals/2e integrals
        void oneint();
        void twoint();

        // calculate the integrals over the primitives
        std::array<double, 3> oneint_prim(int a, int b);
        double twoint_prim(int a, int b, int c, int d);
};

/* 
takes a molecule as constructor input and runs a restricted HF calculation with it
*/
class RHFSolver : public Solver {
    public:
        RHFSolver(const Molecule &mol, double threshold);
        // run the SCF procedure
        double run(bool verbose);
    private:
        // variables for RHF SCF
        Eigen::MatrixXd fock, fock_ortho, p, ao_coeff, mo_energies, coeff;

        // occupation number matrix
        Eigen::MatrixXd nocc;

        // RHF specific methods
        void dens_init();
        void dens();
        void guess();
        double energy();
        void new_fock();
        void mulliken();
};

class UHFSolver : public Solver {

    public:
        UHFSolver(const Molecule &mol, double threshold);
        // run the SCF procedure
        double run(bool verbose);

    private:
        // variables for RHF SCF
        Eigen::MatrixXd fock_alp, fock_bet, fock_alp_ortho, 
        fock_bet_ortho, p_alp, p_bet, ao_coeff_alp, 
        ao_coeff_bet, mo_energies_alp, mo_energies_bet,
        coeff_alp, coeff_bet;

        // occupation number matrices
        Eigen::MatrixXd nocc_alp, nocc_bet;

        // UHF specific methods
        void dens_init();
        void dens();
        void guess();
        double energy();
        void new_fock();
        double contamination();
};