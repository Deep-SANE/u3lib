//protection against multiple header file inclusion
#ifndef U3U2U1XU3U2U1_H
#define U3U2U1XU3U2U1_H

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "U3basis.h"
#include "decomposition.h"
#include "utils/nullspaceSolver.h"
#include "utils/inverseSolver.h"


namespace u3 {

//expansion of a coupled gelfand basis state in terms of those uncoupled ones based on the weight selection rule; hence, this expansion is the same for those other coupled states that have the same weight with the state being expanded
void expansion(u3::irrep &u3a, u3::irrep &u3b, u3::gelfand &u3c, std::vector<u3::gelfand> &gArr1, std::vector<u3::gelfand> &gArr2){
	for(int n12a=u3a.n13; n12a>=u3a.n23; n12a--){
		for(int n22a=u3a.n23; n22a>=u3a.n33; n22a--){
			for(int n11a=n12a; n11a>=n22a; n11a--){
				for(int n12b=u3b.n13; n12b>=u3b.n23; n12b--){
					for(int n22b=u3b.n23; n22b>=u3b.n33; n22b--){
						for(int n11b=n12b; n11b>=n22b; n11b--){
							if ( (n12a+n22a-n11a+n12b+n22b-n11b)==u3c.e22() && (n11a+n11b)==u3c.e11() ){
								u3::gelfand g1(u3a.n13,u3a.n23,u3a.n33,n12a,n22a,n11a);		
								u3::gelfand g2(u3b.n13,u3b.n23,u3b.n33,n12b,n22b,n11b);

								//if(g1.check() && g2.check()){
									gArr1.push_back(g1);
									gArr2.push_back(g2);
								//}
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
	return;
}

void newstates(u3::gelfand &temp_gelfand, std::vector<u3::gelfand> &gArr1, std::vector<u3::gelfand> &gArr2, std::vector<u3::gelfand> &newgArr1, std::vector<u3::gelfand> &newgArr2, int index){
	if(newgArr1.size()==0){
		newgArr1.push_back(temp_gelfand);
		newgArr2.push_back(gArr2[index]);
	}
	else{
		int flag=0;
		for(int j=0; j<newgArr1.size(); j++){
			if( u3::compare(temp_gelfand,newgArr1[j]) && u3::compare(gArr2[index],newgArr2[j]) ){
				flag=1;
				break;
			}
		}
		if(flag==0){
			newgArr1.push_back(temp_gelfand);
			newgArr2.push_back(gArr2[index]);
		}
	}
	return;
}

//build the matrix representating a system of linear equations resulting from the action two generators E23 and E12 that annihiliate the highest weight state 
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> E23E12(std::vector<u3::gelfand> &gArr1, std::vector<u3::gelfand> &gArr2){
	
	int col_num = gArr1.size();

	std::vector<u3::gelfand> newgArr1;
	std::vector<u3::gelfand> newgArr2;

	for(int i=0; i<col_num; i++){
		
		//E23 and E12 act on the first state
		if( (gArr1[i].n12+1)<=gArr1[i].n13 ){
			u3::gelfand temp_gelfand = gArr1[i];
			temp_gelfand.n12++;
			newstates(temp_gelfand,gArr1,gArr2,newgArr1,newgArr2,i);			
		}
		if( (gArr1[i].n22+1)<=gArr1[i].n23 ){
			u3::gelfand temp_gelfand = gArr1[i];
			temp_gelfand.n22++;
			newstates(temp_gelfand,gArr1,gArr2,newgArr1,newgArr2,i);			
		}
		if( (gArr1[i].n11+1)<=gArr1[i].n12 ){
			u3::gelfand temp_gelfand = gArr1[i];
			temp_gelfand.n11++;
			newstates(temp_gelfand,gArr1,gArr2,newgArr1,newgArr2,i);			
		}

		//E23 and E12 act on the second state
		if( (gArr2[i].n12+1)<=gArr2[i].n13 ){
			u3::gelfand temp_gelfand = gArr2[i];
			temp_gelfand.n12++;
			newstates(temp_gelfand,gArr2,gArr1,newgArr2,newgArr1,i);			
		}
		if( (gArr2[i].n22+1)<=gArr2[i].n23 ){
			u3::gelfand temp_gelfand = gArr2[i];
			temp_gelfand.n22++;
			newstates(temp_gelfand,gArr2,gArr1,newgArr2,newgArr1,i);			
		}
		if( (gArr2[i].n11+1)<=gArr2[i].n12 ){
			u3::gelfand temp_gelfand = gArr2[i];
			temp_gelfand.n11++;
			newstates(temp_gelfand,gArr2,gArr1,newgArr2,newgArr1,i);			
		}
	}

	int row_num = newgArr1.size();
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mat(row_num,col_num);

	for(int i=0; i<row_num; i++){
		for(int j=0; j<col_num; j++){

			//action of E23
			if( newgArr1[i].n12==(gArr1[j].n12+1) && newgArr1[i].n22==gArr1[j].n22 && newgArr1[i].n11==gArr1[j].n11 && u3::compare(newgArr2[i],gArr2[j]) ){
				mat(i,j) = gArr1[j].e23_12();
			}
			else if( newgArr1[i].n12==gArr1[j].n12 && newgArr1[i].n22==(gArr1[j].n22+1) && newgArr1[i].n11==gArr1[j].n11 && u3::compare(newgArr2[i],gArr2[j]) ){
				mat(i,j) = gArr1[j].e23_22();
			}
			else if( u3::compare(newgArr1[i],gArr1[j]) && newgArr2[i].n12==(gArr2[j].n12+1) && newgArr2[i].n22==gArr2[j].n22 && newgArr2[i].n11==gArr2[j].n11 ){
				mat(i,j) = gArr2[j].e23_12();
			}
			else if( u3::compare(newgArr1[i],gArr1[j]) && newgArr2[i].n12==gArr2[j].n12 && newgArr2[i].n22==(gArr2[j].n22+1) && newgArr2[i].n11==gArr2[j].n11 ){
				mat(i,j) = gArr2[j].e23_22();
			}

			//action of E12
			else if( newgArr1[i].n12==gArr1[j].n12 && newgArr1[i].n22==gArr1[j].n22 && newgArr1[i].n11==(gArr1[j].n11+1) && u3::compare(newgArr2[i],gArr2[j]) ){
				mat(i,j) = gArr1[j].e12();
			}
			else if( u3::compare(newgArr1[i],gArr1[j]) && newgArr2[i].n12==gArr2[j].n12 && newgArr2[i].n22==gArr2[j].n22 && newgArr2[i].n11==(gArr2[j].n11+1) ){
				mat(i,j) = gArr2[j].e12();
			}

			else {
				mat(i,j) = 0.0;
			}
		}
	}

	return mat;

}


//highest-weight coefficients < G1; G2| HW3 >
void HighestWeightCoeffs(u3::irrep &u3a, u3::irrep &u3b, u3::gelfand &HWc, std::vector<u3::gelfand> &gArr1, std::vector<u3::gelfand> &gArr2, std::vector<u3::gelfand> &gArr3, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &coeffs){

	expansion(u3a,u3b,HWc,gArr1,gArr2);

	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mat = E23E12(gArr1,gArr2);

	coeffs = null_space(mat);

	for(int i=0; i<coeffs.rows(); i++){
		gArr3.push_back(HWc);
	}

	return;
}

//build gelfand basis states from a given weight 
std::vector<u3::gelfand> SameWeight(u3::irrep &u3, int w3, int w1){
	std::vector<u3::gelfand> sameweight;

	for(int n12=u3.n13; n12>=u3.n23; n12--){
		for(int n22=u3.n23; n22>=u3.n33; n22--){
			for(int n11=n12; n11>=n22; n11--){
				if( n11==w1 && (u3.n13+u3.n23+u3.n33-n12-n22)==w3 ){
					u3::gelfand g(u3.n13,u3.n23,u3.n33,n12,n22,n11);
					sameweight.push_back(g);
				}
				else{
					continue;
				}
			}
		}
	}

	return sameweight;
}


//lower-weight coefficients < G1; G2 | G3 >
void LowerWeightCoeffs(u3::irrep &u3a, u3::irrep &u3b, u3::irrep &u3c, int w3, int w1, std::vector<u3::gelfand> &gArr1, std::vector<u3::gelfand> &gArr2, std::vector<u3::gelfand> &gArr3, Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> &coeffs){

	//build all coupled states that have the same weight
	std::vector<u3::gelfand> gArr3_temp = SameWeight(u3c,w3,w1);
	

	if(gArr3_temp.size()!=0){

		//std::cout << gArr3_temp.size() << std::endl;
		std::vector<u3::gelfand> gArr1_temp;
		std::vector<u3::gelfand> gArr2_temp;

		//the coupled states having the same weight have the same expansion of uncoupled states
		expansion(u3a,u3b,gArr3_temp[0],gArr1_temp,gArr2_temp);
		//std::cout << gArr1_temp.size() << std::endl;

		coeffs.conservativeResize(coeffs.rows()+gArr3_temp.size()*gArr1_temp.size(),coeffs.cols());

		std::vector< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> > coeffs_temp;

		for(int outer_multi=0; outer_multi<coeffs.cols(); outer_multi++){

			//build the b-matrix
			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> b(gArr3_temp.size(),gArr3_temp.size());

			//build the lowering matrix from the action of E21 and E32
			Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> E21E32(gArr3_temp.size(),gArr1_temp.size());

			//set matrix elements to zero
			for(int i=0; i<gArr3_temp.size(); i++){
				for(int j=0; j<gArr3_temp.size(); j++){
					b(i,j) = 0.0;
				}
				for(int j=0; j<gArr1_temp.size(); j++){
					E21E32(i,j) = 0.0;
				}
			}

			for(int i=0; i<gArr3_temp.size(); i++){

				int flag = 0;

				//build with generator E21
				u3::gelfand g3_n11p1 = gArr3_temp[i];
				g3_n11p1.n11++; //this is a parent state
				if(g3_n11p1.check()){
					for(int k=0; k<coeffs.rows(); k++){
						if( u3::compare(g3_n11p1,gArr3[k]) ){
							for(int j=0; j<gArr3_temp.size(); j++){
								if( gArr3[k].n12==gArr3_temp[j].n12 && gArr3[k].n22==gArr3_temp[j].n22 && (gArr3[k].n11-1)==gArr3_temp[j].n11 ){
									b(i,j) = gArr3[k].e21();
								}
							}
							for(int j=0; j<gArr1_temp.size(); j++){
								for(int l=k; l<coeffs.rows(); l++){
									if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr2_temp[j],gArr2[l]) && gArr1_temp[j].n12==gArr1[l].n12 && gArr1_temp[j].n22==gArr1[l].n22 && gArr1_temp[j].n11==(gArr1[l].n11-1) ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr1[l].e21();
										flag++;
									}
									else if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr1_temp[j],gArr1[l]) && gArr2_temp[j].n12==gArr2[l].n12 && gArr2_temp[j].n22==gArr2[l].n22 && gArr2_temp[j].n11==(gArr2[l].n11-1) ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr2[l].e21();
										flag++;
									}
								}
							}
							break;
						}
					}
				}

				if(flag>0){
					continue;
				}


				//build with generator E32

				u3::gelfand g3_n12p1 = gArr3_temp[i];
				g3_n12p1.n12++; //this is a parent state
				if(g3_n12p1.check()){
					for(int k=0; k<coeffs.rows(); k++){
						if( u3::compare(g3_n12p1,gArr3[k]) ){
							for(int j=0; j<gArr3_temp.size(); j++){
								if( (gArr3[k].n12-1)==gArr3_temp[j].n12 && gArr3[k].n22==gArr3_temp[j].n22 && gArr3[k].n11==gArr3_temp[j].n11 ){
									b(i,j) = gArr3[k].e32_12();
								}
								else if( gArr3[k].n12==gArr3_temp[j].n12 && (gArr3[k].n22-1)==gArr3_temp[j].n22 && gArr3[k].n11==gArr3_temp[j].n11 ){
									b(i,j) = gArr3[k].e32_22();
								}
							}
							for(int j=0; j<gArr1_temp.size(); j++){
								for(int l=k; l<coeffs.rows(); l++){
									if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr2_temp[j],gArr2[l]) && gArr1_temp[j].n12==(gArr1[l].n12-1) && gArr1_temp[j].n22==gArr1[l].n22 && gArr1_temp[j].n11==gArr1[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr1[l].e32_12();
										flag++;
									}
									else if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr2_temp[j],gArr2[l]) && gArr1_temp[j].n12==gArr1[l].n12 && gArr1_temp[j].n22==(gArr1[l].n22-1) && gArr1_temp[j].n11==gArr1[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr1[l].e32_22();
										flag++;
									}
									else if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr1_temp[j],gArr1[l]) && gArr2_temp[j].n12==(gArr2[l].n12-1) && gArr2_temp[j].n22==gArr2[l].n22 && gArr2_temp[j].n11==gArr2[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr2[l].e32_12();
										flag++;
									}
									else if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr1_temp[j],gArr1[l]) && gArr2_temp[j].n12==gArr2[l].n12 && gArr2_temp[j].n22==(gArr2[l].n22-1) && gArr2_temp[j].n11==gArr2[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr2[l].e32_22();
										flag++;
									}
								}
							}
							break;
						}
					}
				}

				if(flag>0){
					continue;
				}

				u3::gelfand g3_n22p1 = gArr3_temp[i];
				g3_n22p1.n22++; //this is a parent state
				if(g3_n22p1.check()){
					for(int k=0; k<coeffs.rows(); k++){
						if( u3::compare(g3_n22p1,gArr3[k]) ){
							for(int j=0; j<gArr3_temp.size(); j++){
								if( (gArr3[k].n12-1)==gArr3_temp[j].n12 && gArr3[k].n22==gArr3_temp[j].n22 && gArr3[k].n11==gArr3_temp[j].n11 ){
									b(i,j) = gArr3[k].e32_12();
								}
								else if( gArr3[k].n12==gArr3_temp[j].n12 && (gArr3[k].n22-1)==gArr3_temp[j].n22 && gArr3[k].n11==gArr3_temp[j].n11 ){
									b(i,j) = gArr3[k].e32_22();
								}
							}
							for(int j=0; j<gArr1_temp.size(); j++){
								for(int l=k; l<coeffs.rows(); l++){
									if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr2_temp[j],gArr2[l]) && gArr1_temp[j].n12==(gArr1[l].n12-1) && gArr1_temp[j].n22==gArr1[l].n22 && gArr1_temp[j].n11==gArr1[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr1[l].e32_12();
										flag++;
									}
									else if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr2_temp[j],gArr2[l]) && gArr1_temp[j].n12==gArr1[l].n12 && gArr1_temp[j].n22==(gArr1[l].n22-1) && gArr1_temp[j].n11==gArr1[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr1[l].e32_22();
										flag++;
									}
									else if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr1_temp[j],gArr1[l]) && gArr2_temp[j].n12==(gArr2[l].n12-1) && gArr2_temp[j].n22==gArr2[l].n22 && gArr2_temp[j].n11==gArr2[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr2[l].e32_12();
										flag++;
									}
									else if( u3::compare(gArr3[k],gArr3[l]) && u3::compare(gArr1_temp[j],gArr1[l]) && gArr2_temp[j].n12==gArr2[l].n12 && gArr2_temp[j].n22==(gArr2[l].n22-1) && gArr2_temp[j].n11==gArr2[l].n11 ){
										E21E32(i,j) += coeffs(l,outer_multi)*gArr2[l].e32_22();
										flag++;
									}
								}
							}
							break;
						}
					}
				}
			}

			coeffs_temp.push_back(inverse(b)*E21E32);
			//std::cout << b.inverse() << std::endl << E21E32 << std::endl << coeffs_temp << std::endl << std::endl;
		}

		for(int i=0; i<gArr3_temp.size(); i++){
				for(int j=0; j<gArr1_temp.size(); j++){
					gArr1.push_back(gArr1_temp[j]);
					gArr2.push_back(gArr2_temp[j]);
					gArr3.push_back(gArr3_temp[i]);
					for(int outer_multi=0; outer_multi<coeffs.cols(); outer_multi++){
						Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> mat = coeffs_temp[outer_multi];
						if(abs(mat(i,j))<1e-10){
							mat(i,j)=0.0;
						}
						coeffs(gArr1.size()-1,outer_multi) = mat(i,j);
					}
				}
		}
	}
	return; 
}


}

#endif
