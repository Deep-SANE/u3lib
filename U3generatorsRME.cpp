#include <chrono>
#include <iostream>

#include "src/U3basis.h"
#include "src/U3generatorsRME.h"

int main(){

	int lambda,mu;
    std::cout << "Please enter SU(3) quantum numbers (lambda,mu): ";
    std::cin >> lambda >> mu;

    u3::irrep u3top(lambda,mu);

    //building the physical basis excluding those which have zero multiplicity
    std::vector<u3::physicalvector> physbasis = u3::physicalbasis(u3top);

    std::cout << "The physical basis states are: " << std::endl;
    for (int i=0; i<physbasis.size(); i++){
        std::cout << "(" << i << "): |" << "(" << physbasis[i].toprow.lam() << "," << physbasis[i].toprow.mu() << "), kappa_max=" << physbasis[i].kappa << ", L=" << physbasis[i].L << ">" << std::endl;
    }

    std::cout << "This program computes reduced matrix elements of SU(3) generator: <(lambda,mu) kappa' L'|| T || (lambda,mu) kappa L>" << std::endl;

    int option;
    std::cout << "Enter 0 to choose L or 1 to choose Q: ";
    std::cin >> option;

    int i_bra;
    std::cout << "Enter the index of the bra that you are interested: ";
    std::cin >> i_bra;

    int i_ket;
    std::cout << "Enter the index of the ket that you are interested: ";
    std::cin >> i_ket;

    std::vector<int> bra_kappa, ket_kappa;
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> RME;

    // Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

    if(option==0){
    	RME = u3::RME_L(physbasis[i_bra], physbasis[i_ket], bra_kappa, ket_kappa);
    }
    else if(option==1){
    	RME = u3::RME_Q(physbasis[i_bra], physbasis[i_ket], bra_kappa, ket_kappa);
    }

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "\tkappa = ";
    for(int i=0; i<ket_kappa.size(); i++){
    	std::cout << ket_kappa[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "kappa'" << std::endl;
    for(int i=0; i<bra_kappa.size(); i++){
    	std::cout << bra_kappa[i] << "\t\t" << RME.row(i) << std::endl;
    }

    // Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or " << duration.count()/60.0 << " minutes" << std::endl;

	return 0;
}
