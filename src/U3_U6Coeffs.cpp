#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

#include <U3basis.h>
#include <decomposition.h>
#include <U3_U6Z6Coeffs.h>


int main() {

	std::cout << "This program compute U6-recoupling coefficients of three U(3) irreps" << std::endl;

	int n13a, n23a, n33a;
	std::cout << "Please enter the I. irrep [n13a,n23a,n33a]: ";
	std::cin >> n13a >> n23a >> n33a;
	u3::irrep u3a(n13a,n23a,n33a);
	if (n33a!=0){
		u3a.normalize();
		std::cout << "Your irrep is normalized to " << u3a << std::endl;
	}

	int n13b, n23b, n33b;
	std::cout << "Please enter the II. irrep [n13b,n23b,n33b]: ";
	std::cin >> n13b >> n23b >> n33b;
	u3::irrep u3b(n13b,n23b,n33b);
	if (n33b!=0){
		u3b.normalize();
		std::cout << "Your irrep is normalized to " << u3b << std::endl;
	}

	int n13c, n23c, n33c;
	std::cout << "Please enter the III. irrep [n13c,n23c,n33c]: ";
	std::cin >> n13c >> n23c >> n33c;
	u3::irrep u3c(n13c,n23c,n33c);
	if (n33c!=0){
		u3c.normalize();
		std::cout << "Your irrep is normalized to " << u3c << std::endl;
	}

	//couple the I. and II. irreps
	std::vector<int> rhomaxsab;
	std::vector<u3::irrep> u3ab = u3::decomposition(u3a,u3b,rhomaxsab);

	std::cout << "The normalized coupled irreps of the I and II are: " << std::endl;
	for(int i=0; i<u3ab.size(); i++){
		std::cout << "(" << i << "): [" << u3ab[i].n13-u3ab[i].n33 << "," << u3ab[i].n23-u3ab[i].n33 << ",0] \t rho = 1 ... " << rhomaxsab[i] << std::endl;
	}

	//couple the II. and III. irreps
	std::vector<int> rhomaxsbc;
	std::vector<u3::irrep> u3bc = u3::decomposition(u3b,u3c,rhomaxsbc);

	std::cout << "The normalized coupled irreps of the II and III are: " << std::endl;
	for(int i=0; i<u3bc.size(); i++){
		std::cout << "(" << i << "): [" << u3bc[i].n13-u3bc[i].n33 << "," << u3bc[i].n23-u3bc[i].n33 << ",0] \t rho = 1 ... " << rhomaxsbc[i] << std::endl;
	}

	int i_ab, i_bc, i_abc;
	std::cout << "Please enter the index of the I-II coupled irrep you are interested: ";
	std::cin >> i_ab;
	std::cout << "Please enter the index of the II-III coupled irrep you are interested: ";
	std::cin >> i_bc;

	//couple the three irreps in scheme (I-II)-III
	std::vector<int> rhomaxsab_c;
	std::vector<u3::irrep> u3ab_c = u3::decomposition(u3ab[i_ab],u3c,rhomaxsab_c);

	std::cout << "The normalized coupled irreps of the I, II and III are: " << std::endl;
	for(int i=0; i<u3ab_c.size(); i++){
		std::cout << "(" << i << "): [" << u3ab_c[i].n13-u3ab_c[i].n33 << "," << u3ab_c[i].n23-u3ab_c[i].n33 << ",0] \t rho = 1 ... " << rhomaxsab_c[i] << std::endl;
	}
	std::cout << "Please enter the index of the I-II-III coupled irrep you are interested: ";
	std::cin >> i_abc;

	//couple the three irreps in scheme I-(II-III)
	std::vector<int> rhomaxsa_bc;
	std::vector<u3::irrep> u3a_bc = u3::decomposition(u3a,u3bc[i_bc],rhomaxsa_bc);

	int flag=0;
	for(int j=0; j<u3a_bc.size(); j++){
		if( u3ab_c[i_abc].n13==u3a_bc[j].n13 && u3ab_c[i_abc].n23==u3a_bc[j].n23 && u3ab_c[i_abc].n33==u3a_bc[j].n33 ){
			flag++;
			break;
		}
	}

	if(flag>0){
		std::vector<int> rhoab, rhoab_c, rhobc, rhoa_bc;

		// Start time measurement
    	auto t0 = std::chrono::high_resolution_clock::now();

		std::vector<double> coeffs = u3::U6(u3a,u3b,u3ab_c[i_abc],u3c,u3ab[i_ab],rhoab,rhoab_c,u3bc[i_bc],rhobc,rhoa_bc);

		// End time measurement
    	auto t1 = std::chrono::high_resolution_clock::now();

		u3ab[i_ab].normalize();
		u3bc[i_bc].normalize();
		u3ab_c[i_abc].normalize();

		std::cout << "U([ni3a], [ni3b], [ni3],  [ni3c] ;[ni3ab],rab,rab_c,[ni3bc],rbc,ra_bc)" << std::endl;
		std::cout << "  " << u3a << "," << u3b << "," << u3ab_c[i_abc] << "," << u3c << ";" << u3ab[i_ab] << ",\t    " << u3bc[i_bc] << std::endl;

		for(int i=0; i<coeffs.size(); i++){
			std::cout << "\t\t\t\t\t  " << rhoab[i] << "   " << rhoab_c[i] << "\t\t     " << rhobc[i] << "  " << rhoa_bc[i] << "\t -> " << coeffs[i] << std::endl;
		}

		// Calculate the duration
    	std::chrono::duration<double> duration = t1 - t0;

    	// Output the time in seconds
    	std::cout << "Elapsed time: " << duration.count() << " seconds or " << duration.count()/60.0 << " minutes" << std::endl;
	}
	else{
		std::cout << "Recoupling does not exist" << std::endl;
	}

	return 0;
}
