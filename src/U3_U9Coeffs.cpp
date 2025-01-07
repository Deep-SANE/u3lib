#include <iostream>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>

#include <U3basis.h>
#include <decomposition.h>
#include <U3_U9Coeffs.h>


int main() {

	std::cout << "This program compute U9-recoupling coefficients of four U(3) irreps" << std::endl;

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

	int n13d, n23d, n33d;
	std::cout << "Please enter the IV. irrep [n13d,n23d,n33d]: ";
	std::cin >> n13d >> n23d >> n33d;
	u3::irrep u3d(n13d,n23d,n33d);
	if (n33d!=0){
		u3d.normalize();
		std::cout << "Your irrep is normalized to " << u3d << std::endl;
	}


	//couple the I. and II. irreps
	std::vector<int> rhomaxsab;
	std::vector<u3::irrep> u3ab = u3::decomposition(u3a,u3b,rhomaxsab);

	std::cout << "The normalized coupled irreps of the I and II are: " << std::endl;
	for(int i=0; i<u3ab.size(); i++){
		std::cout << "(" << i << "): [" << u3ab[i].n13-u3ab[i].n33 << "," << u3ab[i].n23-u3ab[i].n33 << ",0] \t rho = 1 ... " << rhomaxsab[i] << std::endl;
	}

	//couple the I. and III. irreps
	std::vector<int> rhomaxsac;
	std::vector<u3::irrep> u3ac = u3::decomposition(u3a,u3c,rhomaxsac);

	std::cout << "The normalized coupled irreps of the I and III are: " << std::endl;
	for(int i=0; i<u3ac.size(); i++){
		std::cout << "(" << i << "): [" << u3ac[i].n13-u3ac[i].n33 << "," << u3ac[i].n23-u3ac[i].n33 << ",0] \t rho = 1 ... " << rhomaxsac[i] << std::endl;
	}

	//couple the II. and IV. irreps
	std::vector<int> rhomaxsbd;
	std::vector<u3::irrep> u3bd = u3::decomposition(u3b,u3d,rhomaxsbd);

	std::cout << "The normalized coupled irreps of the II and IV are: " << std::endl;
	for(int i=0; i<u3bd.size(); i++){
		std::cout << "(" << i << "): [" << u3bd[i].n13-u3bd[i].n33 << "," << u3bd[i].n23-u3bd[i].n33 << ",0] \t rho = 1 ... " << rhomaxsbd[i] << std::endl;
	}

	//couple the III. and IV. irreps
	std::vector<int> rhomaxscd;
	std::vector<u3::irrep> u3cd = u3::decomposition(u3c,u3d,rhomaxscd);

	std::cout << "The normalized coupled irreps of the III and IV are: " << std::endl;
	for(int i=0; i<u3cd.size(); i++){
		std::cout << "(" << i << "): [" << u3cd[i].n13-u3cd[i].n33 << "," << u3cd[i].n23-u3cd[i].n33 << ",0] \t rho = 1 ... " << rhomaxscd[i] << std::endl;
	}

	int i_ab, i_ac, i_bd, i_cd, i_ab_cd, i_ac_bd;
	std::cout << "Please enter the index of the I-II coupled irrep you are interested: ";
	std::cin >> i_ab;
	std::cout << "Please enter the index of the I-III coupled irrep you are interested: ";
	std::cin >> i_ac;
	std::cout << "Please enter the index of the II-IV coupled irrep you are interested: ";
	std::cin >> i_bd;
	std::cout << "Please enter the index of the III-IV coupled irrep you are interested: ";
	std::cin >> i_cd;

	//couple the 4 irreps in the scheme (I-II)-(III-IV)
	std::vector<int> rhomaxsab_cd;
	std::vector<u3::irrep> u3ab_cd = u3::decomposition(u3ab[i_ab],u3cd[i_cd],rhomaxsab_cd);

	std::cout << "The coupled irreps of the I, II, III and IV are: " << std::endl;
	for(int i=0; i<u3ab_cd.size(); i++){
		std::cout << "(" << i << "): [" << u3ab_cd[i].n13-u3ab_cd[i].n33 << "," << u3ab_cd[i].n23-u3ab_cd[i].n33 << ",0] \t rho = 1 ... " << rhomaxsab_cd[i] << std::endl;
	}
	std::cout << "Please enter the index of the I-II-III-IV coupled irrep you are interested: ";
	std::cin >> i_ab_cd;

	//couple the 4 irreps in scheme (I-III)-(II-IV)
	std::vector<int> rhomaxsac_bd;
	std::vector<u3::irrep> u3ac_bd = u3::decomposition(u3ac[i_ac],u3bd[i_bd],rhomaxsac_bd);

	int flag1=0, flag2=0;
	for(int j=0; j<u3ac_bd.size(); j++){
		if( u3ab_cd[i_ab_cd].n13==u3ac_bd[j].n13 && u3ab_cd[i_ab_cd].n23==u3ac_bd[j].n23 && u3ab_cd[i_ab_cd].n33==u3ac_bd[j].n33 ){
			i_ac_bd = j;
			flag1++;
			break;
		}
	}

	//couple the 4 irreps in scheme ((I-II)-III)-IV
	std::vector<int> rhomaxsab_c;
	std::vector<u3::irrep> u3ab_c = u3::decomposition(u3ab[i_ab],u3c,rhomaxsab_c);

	//couple the 4 irreps in scheme ((I-III)-II)-IV
	std::vector<int> rhomaxsac_b;
	std::vector<u3::irrep> u3ac_b = u3::decomposition(u3ac[i_ac],u3b,rhomaxsac_b);

	//an array to store the couple I-II-III irreps
	std::vector<int> rhomaxsabc;
	std::vector<u3::irrep> u3abc;

	for(int i=0; i<u3ab_c.size(); i++){
		for(int j=0; j<u3ac_b.size(); j++){
			if( u3ab_c[i].n13==u3ac_b[j].n13 && u3ab_c[i].n23==u3ac_b[j].n23 && u3ab_c[i].n33==u3ac_b[j].n33 ){
				u3abc.push_back(u3ab_c[i]);
				rhomaxsabc.push_back(rhomaxsab_c[i]);
				break;
			}
		}
	}

	for(int i=0; i<u3abc.size(); i++){
		//std::cout << u3abc[i] << std::endl;
		std::vector<int> rhomaxsabc_d;
		std::vector<u3::irrep> u3abc_d = u3::decomposition(u3abc[i],u3d,rhomaxsabc_d);
		for(int j=0; j<u3abc_d.size(); j++){
			if( u3ab_cd[i_ab_cd].n13==u3abc_d[j].n13 && u3ab_cd[i_ab_cd].n23==u3abc_d[j].n23 && u3ab_cd[i_ab_cd].n33==u3abc_d[j].n33 ){
				flag2++;
				break;
			}
		}
	}

	if(flag1>0 && flag2==u3abc.size()){
		std::vector<int> rhoab, rhocd, rhoac, rhobd, rhoab_cd, rhoac_bd;

		// Start time measurement
    	auto t0 = std::chrono::high_resolution_clock::now();

		Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> coeffs = u3::U9(u3a,u3b,u3c,u3d,u3ab[i_ab],u3cd[i_cd],u3ac[i_ac],u3bd[i_bd],u3ab_cd[i_ab_cd], 
																rhoab,rhomaxsab[i_ab], rhocd,rhomaxscd[i_cd], rhoab_cd,rhomaxsab_cd[i_ab_cd], 
																rhoac,rhomaxsac[i_ac], rhobd,rhomaxsbd[i_bd], rhoac_bd,rhomaxsac_bd[i_ac_bd]);
		
		// End time measurement
    	auto t1 = std::chrono::high_resolution_clock::now();

		u3ab[i_ab].normalize();
		u3cd[i_cd].normalize();
		u3ac[i_ac].normalize();
		u3bd[i_bd].normalize();
		u3ab_cd[i_ab_cd].normalize();

		std::cout << u3a << "\t" << u3b << "\t" << u3ab[i_ab] << "\trhoab\t\t";
		for(int i=0; i<rhoab.size(); i++){
			std::cout << rhoab[i] << "\t";
		}
		std::cout << std::endl;
		std::cout << u3c << "\t" << u3d << "\t" << u3cd[i_cd] << "\trhocd\t\t";
		for(int i=0; i<rhocd.size(); i++){
			std::cout << rhocd[i] << "\t";
		}
		std::cout << std::endl;
		std::cout << u3ac[i_ac] << "\t" << u3bd[i_bd] << "\t" << u3ab_cd[i_ab_cd] << "\trhoac_bd\t";
		for(int i=0; i<rhoac_bd.size(); i++){
			std::cout << rhoac_bd[i] << "\t";
		}
		std::cout << std::endl;
		
		std::cout << "rhoac\trhobd\trhoab_cd\t" << std::endl;

		for(int i=0; i<rhoac.size(); i++){
			std::cout << rhoac[i] << "\t" << rhobd[i] << "\t" << rhoab_cd[i] << "\t\t\t" << coeffs.row(i) << std::endl;
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
