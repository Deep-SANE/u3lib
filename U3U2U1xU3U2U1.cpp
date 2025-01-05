#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

#include "src/U3basis.h"
#include "src/decomposition.h"
#include "src/U3U2U1xU3U2U1.h"


int main() {

	int n13a, n23a, n33a;
	std::cout << "Please enter the first irrep [n13a,n23a,n33a]: ";
	std::cin >> n13a >> n23a >> n33a;
	u3::irrep u3a(n13a,n23a,n33a);
	if (n33a!=0){
		u3a.normalize();
		std::cout << "Your irrep is normalized to " << u3a << std::endl;
	}

	int n13b, n23b, n33b;
	std::cout << "Please enter the second irrep [n13b,n23b,n33b]: ";
	std::cin >> n13b >> n23b >> n33b;
	u3::irrep u3b(n13b,n23b,n33b);
	if (n33b!=0){
		u3b.normalize();
		std::cout << "Your irrep is normalized to " << u3b << std::endl;
	}

	std::vector<int> rhomaxs;
	std::vector<u3::irrep> u3cs = decomposition(u3a,u3b,rhomaxs);

	std::cout << "The normalized coupled irreps are: " << std::endl;
	for(int i=0; i<u3cs.size(); i++){
		std::cout << "(" << i << "): [" << u3cs[i].n13-u3cs[i].n33 << "," << u3cs[i].n23-u3cs[i].n33 << ",0] \t rho = 1 ... " << rhomaxs[i] << std::endl;
	}

	int index;
	std::cout << "Please enter the index of the coupled irrep you are interested: ";
	std::cin >> index;

	// Start time measurement
    auto t0 = std::chrono::high_resolution_clock::now();

	std::vector<u3::gelfand> gArr1;
	std::vector<u3::gelfand> gArr2;
	std::vector<u3::gelfand> gArr3;

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> coeffs;

	for(int w3=u3cs[index].n33; w3<=u3cs[index].n13; w3++){
		for(int w1=u3cs[index].n13; w1>=u3cs[index].n33; w1--){
			if( w3==u3cs[index].n33 && w1==u3cs[index].n13){
				u3::gelfand HWc(u3cs[index].n13,u3cs[index].n23,u3cs[index].n33,u3cs[index].n13,u3cs[index].n23,u3cs[index].n13); 
				u3::HighestWeightCoeffs(u3a,u3b,HWc,gArr1,gArr2,gArr3,coeffs);
			}
			else{
				u3::LowerWeightCoeffs(u3a,u3b,u3cs[index],w3,w1,gArr1,gArr2,gArr3,coeffs);
			}
		}
	}

	// End time measurement
    auto t1 = std::chrono::high_resolution_clock::now();


	std::cout << "< Ga; Gb | Gc > \t\t rho = 1 ... " << rhomaxs[index] << std::endl;
	for(int j=0; j<gArr1.size(); j++){
		gArr3[j].normalize();
		std::cout << "< " << gArr1[j] << " ; " << gArr2[j] << " | " << gArr3[j] << " >" << "\t" << coeffs.row(j) << std::endl;
	}
	
	// Calculate the duration
    std::chrono::duration<double> duration = t1 - t0;

    // Output the time in seconds
    std::cout << "Elapsed time: " << duration.count() << " seconds or " << duration.count()/60.0 << " minutes" << std::endl;

	return 0;
}