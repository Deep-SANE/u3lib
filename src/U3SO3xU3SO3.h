//protection against multiple header file inclusion
#ifndef U3SO3xU3SO3_H
#define U3SO3xU3SO3_H

#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "U3basis.h"
#include "U3SO3_to_U3U2U1.h"
#include "U3U2U1xU3U2U1.h"

extern "C"
{
#include <wigxjpf.h>
}

namespace u3{

//usual SO(3) Clebsh-Gordan coefficients <L1 M1; L2 M2 | L3 M3>
double CGC_SO3(int l1, int l2, int l3, int m1, int m2, int m3) {
   	int phase = ((-l1 + l2 - m3) % 2 == 0) ? 1 : -1;

   	wig_table_init(2*100,9);
  	wig_temp_init(2*100);

  	double val3j=wig3jj(2*l1, 2*l2, 2*l3, 2*m1, 2*m2, 2*(-m3));

  	wig_temp_free();
  	wig_table_free();

   	return phase * std::sqrt(2*l3 + 1) * val3j;
}

//generate all Wigner U3->SO(3) coefficients
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> WigU3SO3(u3::irrep &u3A, u3::irrep &u3B, u3::irrep &u3C, std::vector<u3::physicalvector> &physvecA, std::vector<u3::physicalvector> &physvecB, std::vector<u3::physicalvector> &physvecC){

	//find all the canonical coupling coefficients
	std::vector<u3::gelfand> gArrA;
	std::vector<u3::gelfand> gArrB;
	std::vector<u3::gelfand> gArrC;

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> canoCoupleCoeffs; //each column corresponds to one rho

	for(int w3=u3C.n33; w3<=u3C.n13; w3++){
		for(int w1=u3C.n13; w1>=u3C.n33; w1--){
			if( w3==u3C.n33 && w1==u3C.n13){
				u3::gelfand HW3(u3C.n13,u3C.n23,u3C.n33,u3C.n13,u3C.n23,u3C.n13); 
				u3::HighestWeightCoeffs(u3A,u3B,HW3,gArrA,gArrB,gArrC,canoCoupleCoeffs);
			}
			else{
				u3::LowerWeightCoeffs(u3A,u3B,u3C,w3,w1,gArrA,gArrB,gArrC,canoCoupleCoeffs);
			}
		}
	}
	
	//find expansion for the second irrep
    std::vector<u3::physicalvector> pbB = u3::physicalbasis(u3B);
    std::vector<int> kappaB;
    std::vector< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > allCanoCoeffsB; //each column corresponds to one kappa
    std::vector< std::vector<u3::gelfand> > allCanoBasisB;
    for (int i=0; i<pbB.size(); i++){
     
        //Construct the canonical basis associated with the physical basis state that we are looking at
        std::vector<u3::gelfand> canobasis = u3::canonicalbasis(pbB[i]);

        //Building the projection matrix 
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Lproj = Lprojection(pbB[i],canobasis);
        //and solving its null space
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> nullspace = phys_to_cano(pbB[i],Lproj); //after solving the nullspace, the dimension of the nullspace is the inner multiplicity of L

   		kappaB.push_back(pbB[i].kappa); //this multiplicity was computed after solving the nullspace
   		allCanoCoeffsB.push_back(nullspace);
        allCanoBasisB.push_back(canobasis);
    }


    //find expansion for the third irrep
    //first normalize the irrep to obtain physical basis states |(lamda,mu),kappa,L>
    int n33 = u3C.n33, lam3=u3C.n13-u3C.n23, mu3=u3C.n23-u3C.n33;
    u3::irrep u3c(lam3,mu3);
    std::vector<u3::physicalvector> pbC = u3::physicalbasis(u3c);
    std::vector<int> kappaC;
    std::vector< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > allCanoCoeffsC; //each column corresponds to one kappa
    std::vector< std::vector<u3::gelfand> > allCanoBasisC;
    for (int i=0; i<pbC.size(); i++){

        //Construct the canonical basis associated with the physical basis state that we are looking at
        std::vector<u3::gelfand> canobasis = u3::canonicalbasis(pbC[i]);

        //Building the projection matrix 
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Lproj = Lprojection(pbC[i],canobasis);
        //and solving its null space
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> nullspace = phys_to_cano(pbC[i],Lproj); //after solving the nullspace, the dimension of the nullspace is the inner multiplicity of L

   		kappaC.push_back(pbC[i].kappa); //this multiplicity was computed after solving the nullspace
   		allCanoCoeffsC.push_back(nullspace);

   		//std::cout << pbC[i] << std::endl;
   		//add back n33 to the canonical basis
   		for(int j=0; j<canobasis.size(); j++){
   			canobasis[j].n13 += n33;
   			canobasis[j].n23 += n33;
   			canobasis[j].n33 += n33;
   			canobasis[j].n12 += n33;
   			canobasis[j].n22 += n33;
   			canobasis[j].n11 += n33;

   			//std::cout << canobasis[j] << std::endl;
   		}
        allCanoBasisC.push_back(canobasis);
    }
    

    //compute Wigner coefficients: < (lamda1,mu1),kappa1,L1 ; (lamda2,mu2),kappa2,L2 | (lamda3,mu3),kappa3,L3 >
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> coeffs;
    int rhomax = canoCoupleCoeffs.cols();

    
    //find expansion for the first irrep
    std::vector<u3::physicalvector> pbA = u3::physicalbasis(u3A);
    for (int i=0; i<pbA.size(); i++){
     	
    	std::vector< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > allCanoCoeffsA; //each column corresponds to one kappa -- the inner multiplicity
    	std::vector< std::vector<u3::gelfand> > allCanoBasisA; //to store all canonical basis states associated with all physical basis states (L,M) with a given L

        //Construct the canonical basis associated with the heighest weight physical basis state (L=M) 
        allCanoBasisA.push_back(u3::canonicalbasis(pbA[i]));
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> Lproj = Lprojection(pbA[i],allCanoBasisA[0]);
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> nullspace = phys_to_cano(pbA[i],Lproj);
   		allCanoCoeffsA.push_back(nullspace);

        std::vector<double> gamma_coefficients; //to store sqrt((2L-gamma)! 2^gamma / (2L)! / gamma!)
        gamma_coefficients.push_back(1.0); //for gamma=0, i.e., L=M

        //now going lower from the HW state, remember that gamma = L-M
        for(int gamma=1; gamma<=2*pbA[i].L; gamma++){
        	//extract the expansion of the previous state
        	std::vector<u3::gelfand> cano_previous = allCanoBasisA[gamma-1];
        	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> canocoeffs_prev = allCanoCoeffsA[gamma-1];

        	//define a new array to store new expansion, in which all canonical basis states are obtained by acting with two generators E31 and E23
        	std::vector<u3::gelfand> cano_now; 

        	for(int j=0; j<cano_previous.size(); j++){
        		//action of E31
        		u3::gelfand g1 = cano_previous[j];
				g1.n12--;
				g1.n11--;
				if(g1.check()){
					u3::newstate(g1,cano_now);
				}

				u3::gelfand g2 = cano_previous[j];
				g2.n22--;
				g2.n11--;
				if(g2.check()){
					u3::newstate(g2,cano_now);
				}

				//determine the action of E23
				u3::gelfand g3 = cano_previous[j];
				g3.n12++;
				if(g3.check()){
					u3::newstate(g3,cano_now);
				}
				u3::gelfand g4 = cano_previous[j];
				g4.n22++;
				if(g4.check()){
					u3::newstate(g4,cano_now);
				}
        	}

        	//calculate expansion coefficients for the new state
        	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> canocoeffs_now(cano_now.size(),pbA[i].kappa);
        	for(int j=0; j<canocoeffs_now.rows(); j++){
				for(int k=0; k<canocoeffs_now.cols(); k++){
					canocoeffs_now(j,k) = 0.0;
				}
			}

			for(int j=0; j<cano_now.size(); j++){
				for(int m=0; m<cano_previous.size(); m++){
					//action of E31
					if( cano_now[j].n12==(cano_previous[m].n12-1) && cano_now[j].n22==cano_previous[m].n22 && cano_now[j].n11==(cano_previous[m].n11-1) ){
						for(int k=0; k<pbA[i].kappa; k++){
							canocoeffs_now(j,k) += canocoeffs_prev(m,k)*cano_previous[m].e31_12();
						}
					}
					else if( cano_now[j].n12==cano_previous[m].n12 && cano_now[j].n22==(cano_previous[m].n22-1) && cano_now[j].n11==(cano_previous[m].n11-1) ){
						for(int k=0; k<pbA[i].kappa; k++){
							canocoeffs_now(j,k) += canocoeffs_prev(m,k)*cano_previous[m].e31_22();
						}
					}

					//action of E23
					else if( cano_now[j].n12==(cano_previous[m].n12+1) && cano_now[j].n22==cano_previous[m].n22 && cano_now[j].n11==cano_previous[m].n11 ){
						for(int k=0; k<pbA[i].kappa; k++){
							canocoeffs_now(j,k) += canocoeffs_prev(m,k)*cano_previous[m].e23_12();
						}
					}
					else if( cano_now[j].n12==cano_previous[m].n12 && cano_now[j].n22==(cano_previous[m].n22+1) && cano_now[j].n11==cano_previous[m].n11 ){
						for(int k=0; k<pbA[i].kappa; k++){
							canocoeffs_now(j,k) += canocoeffs_prev(m,k)*cano_previous[m].e23_22();
						}
					}
				}
			}

  	     	allCanoBasisA.push_back(cano_now);
        	allCanoCoeffsA.push_back(canocoeffs_now);

        	//compute gamma coefficients
        	double LL=1.0, LLgamma, gammaFac;
        	for(int j=1; j<=2*pbA[i].L; j++){
        		LL *= j;
        		if(j==gamma){
        			gammaFac = LL;
        		}
        		if(j==(2*pbA[i].L-gamma)){
        			LLgamma = LL;
        		}
        	}
        	gamma_coefficients.push_back( std::sqrt( LLgamma*std::pow(2,gamma)/LL/gammaFac ) );
        }

        //last step is to compute Wigner coefficients < (lamda1,mu1)kappa1 L1; (lamda2,mu2)kappa2 L2 | (lamda3,mu3)kappa3 L3 >
        for(int k_A=1; k_A<=pbA[i].kappa; k_A++){
        	for(int j=0; j<pbB.size(); j++){
        		for(int k_B=1; k_B<=pbB[j].kappa; k_B++){
        			for(int k=0; k<pbC.size(); k++){
        				for(int k_C=1; k_C<=pbC[k].kappa; k_C++){

        					
        					double so3cgc = CGC_SO3(pbA[i].L, pbB[j].L, pbC[k].L, pbC[k].L-pbB[j].L, pbB[j].L, pbC[k].L);

        					if(so3cgc!=0.0){
        						int gamma = pbA[i].L + pbB[j].L - pbC[k].L;

        						u3::physicalvector pb1 = pbA[i], pb2 = pbB[j], pb3 = pbC[k];
        						pb1.kappa = k_A;
        						physvecA.push_back(pb1);
        						pb2.kappa = k_B;
        						physvecB.push_back(pb2);
        						pb3.kappa = k_C;
        						physvecC.push_back(pb3);

        						std::vector<u3::gelfand> cb1 = allCanoBasisA[gamma], cb2 = allCanoBasisB[j], cb3 = allCanoBasisC[k];
        						Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> canocoeffs1 = allCanoCoeffsA[gamma], canocoeffs2 = allCanoCoeffsB[j], canocoeffs3 = allCanoCoeffsC[k];

        						int old_row = coeffs.rows();

        						coeffs.conservativeResize(old_row+1,rhomax);

        						for(int rho=1; rho<=rhomax; rho++){
        							coeffs(old_row,rho-1) = 0.0;

        							for(int m=0; m<cb1.size(); m++){
        								for(int n=0; n<cb2.size(); n++){
        									for(int l=0; l<cb3.size(); l++){
        										double canocoupling = 0.0;
        										for(int p=0; p<canoCoupleCoeffs.rows(); p++){
        											if( u3::compare(cb1[m],gArrA[p]) && u3::compare(cb2[n],gArrB[p]) && u3::compare(cb3[l],gArrC[p]) ){
        												canocoupling = canoCoupleCoeffs(p,rho-1);
        											}
        										}
        										coeffs(old_row,rho-1) += canocoeffs1(m,k_A-1)*canocoeffs2(n,k_B-1)*canocoeffs3(l,k_C-1)*canocoupling;
        									}
        								}
        							}

        							coeffs(old_row,rho-1) *= (1/so3cgc)*gamma_coefficients[gamma];
        						}

        					}
        					else{
        						continue;
        					}
        				}
        			}
        		}
        	}
        }

        
    }

    return coeffs;
}

}

#endif