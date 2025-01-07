#include <iostream>
#include <iomanip>
#include <vector>
#include <chrono>

#include <U3basis.h>
#include <U3SO3_to_U3U2U1.h>
#include <utils/nullspaceSolver.h>


int main() {

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

    int index;
    std::cout << "Please enter the index of the state that you are interested: ";
    std::cin >> index;

    // Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

    //Construct the canonical basis associated with the physical basis state that we are looking at
    std::vector<u3::gelfand> canobasis = u3::canonicalbasis(physbasis[index]);

    //Building the projection matrix 
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Lproj = Lprojection(physbasis[index],canobasis);
    //and solving its null space
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> nullspace = phys_to_cano(physbasis[index],Lproj); 
    //after solving the nullspace, the dimension of the nullspace is the multiplicity of L, and the multiplicity within the physical basis state labeling is modified  

    // End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();  
          
    std::cout << "The Gelfand states and the expansion coefficients for kappa = 1,...," << physbasis[index].kappa << " are:" << std::endl;

    for (int i=0; i<canobasis.size(); i++){
        std::cout << "|" << canobasis[i] << "> \t" << nullspace.row(i) << std::endl;
    } 

    // Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or " << duration.count()/60.0 << " minutes" << std::endl;


    return 0;
}
