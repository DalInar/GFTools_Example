#include <iostream>
#include "math.h"
#include <alps/gf/tail.hpp>
#include "alps/hdf5.hpp"

typedef alps::gf::omega_k_sigma_gf_with_tail cluster_matsubara_kspace_gf_with_tail;
typedef alps::gf::omega_k_sigma_gf cluster_matsubara_kspace_gf;
typedef alps::gf::two_index_gf<double, alps::gf::momentum_index_mesh, alps::gf::index_mesh> cluster_gf_tail;

int main() {
    double beta=1;
    int nfreq = 512;
    std::cout << "Example of creating, saving, loading, and manipulating GFTools Greens functions" << std::endl;
    std::cout << "First, create a matsubara greens function with Beta="<<beta<<", NFreq="<<nfreq<<" on a 2x2 square lattice (4 sites)" <<std::endl;
    alps::gf::matsubara_positive_mesh freq_mesh(beta, nfreq);

    std::cout << "Now create the 2x2 momentum mesh"<<std::endl;
    alps::gf::momentum_index_mesh::container_type points(boost::extents[4][2]);
    points[0][0] = 0.0; points[0][1] = 0.0;
    points[1][0] = 0.0; points[1][1] = M_PI;
    points[2][0] = M_PI; points[2][1] = 0.0;
    points[3][0] = M_PI; points[3][1] = M_PI;
    alps::gf::momentum_index_mesh mom_mesh(points);

    std::cout << "Now create the spin (flavor) mesh"<<std::endl;
    alps::gf::index_mesh spins(2);

    std::cout << "Create the example Greens function, G(k,iw)"<<std::endl;
    cluster_matsubara_kspace_gf_with_tail G(
            cluster_matsubara_kspace_gf(freq_mesh, mom_mesh, spins)
    );
    std::cout<< "G is indexed by G(frequency, momentum, spin)"<<std::endl;

    std::cout << "Create dummy high frequency tail coefficients" <<std::endl;
    cluster_gf_tail c1(G.mesh2(), G.mesh3());   //The tail is a function of the same momentum and spin mesh as G
    cluster_gf_tail c2(G.mesh2(), G.mesh3());

    for(alps::gf::index f(0); f<2; f++){
        for(alps::gf::momentum_index k(0); k<4; k++){
            c1(k,f) = 1.;
            c2(k,f) = -2*(std::cos(mom_mesh.points()[k()][0])+std::cos(mom_mesh.points()[k()][1]));
        }
    }

    std::cout << "Now set the tails for the Greens function"<<std::endl;
    G.set_tail(1,c1);
    G.set_tail(2,c2);

    std::cout << "Set the elements of G to arbitrary values" << std::endl;
    for(alps::gf::matsubara_index w(0); w<nfreq; w++){
        for(alps::gf::momentum_index k(0); k<4; k++){
            for(alps::gf::index f(0); f<2; f++){
                G(w,k,f) = std::complex<double>(w()*M_PI*k()*f(), k()+1.);
            }
        }
    }

    std::cout << "The element G(4,2,1) is "<<G(alps::gf::matsubara_index(4), alps::gf::momentum_index(2), alps::gf::index(1)) << std::endl;
    std::cout << "The element (0,0) of the first order tail of G (c1) is "<<G.tail(1)(alps::gf::momentum_index(0), alps::gf::index(0)) << std::endl;

    std::cout << "Save G to h5 file (G.h5)" <<std::endl;
    alps::hdf5::archive ar("G.h5", alps::hdf5::archive::WRITE);
    G.save(ar, "/G");
    ar.close();

    std::cout << "Read to G_new from h5 file (G.h5)"<<std::endl;
    alps::hdf5::archive ar_read("G.h5", alps::hdf5::archive::READ);
    cluster_matsubara_kspace_gf_with_tail G_new(alps::gf::omega_k_sigma_gf(
            alps::gf::matsubara_positive_mesh(5, 2),
            alps::gf::momentum_index_mesh(9,2),
            alps::gf::index_mesh(2)));  //This new Greens function is initialized with the wrong meshes, but this will be corrected during load

    G_new.load(ar_read, "/G");

    std::cout << "The element G_new(4,2,1) is "<<G_new(alps::gf::matsubara_index(4), alps::gf::momentum_index(2), alps::gf::index(1)) << std::endl;
    std::cout << "The element (0,0) of the first order tail of G_new (c1) is "<<G_new.tail(1)(alps::gf::momentum_index(0), alps::gf::index(0)) << std::endl;

    std::cout << "By the way, what frequency, momentum, and spin values are at G(4,2,1)?" <<std::endl;
    std::cout << "Frequency 4 is w="<<freq_mesh.points()[4]<<std::endl;
    std::cout << "Momentum 2 is k=("<< mom_mesh.points()[2][0] <<", "<< mom_mesh.points()[2][1]<<")"<<std::endl;
    std::cout << "Spin 1 is s="<<spins.points()[1]<<std::endl;


    return 0;
}