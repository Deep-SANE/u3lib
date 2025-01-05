//protection against multiple header file inclusion
#ifndef U3_U9COEFFS_H
#define U3_U9COEFFS_H

#include <iostream>
#include <vector>
#include <cmath>
#include <Eigen/Dense>

#include "U3basis.h"
#include "decomposition.h"
#include "U3_U6Z6Coeffs.h"


namespace u3 {


//input by users: u3a, u3b, u3c, u3d, u3ab, u3cd, u3ac, u3bd, u3abcd
Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> U9(u3::irrep &u3a, u3::irrep &u3b, u3::irrep &u3c, u3::irrep &u3d, u3::irrep &u3ab, u3::irrep &u3cd, u3::irrep &u3ac, u3::irrep &u3bd, u3::irrep &u3abcd,
																		std::vector<int> &rhoab, int rhoab_max, std::vector<int> &rhocd, int rhocd_max, std::vector<int> &rhoab_cd, int rhoab_cd_max, 
																		std::vector<int> &rhoac, int rhoac_max, std::vector<int> &rhobd, int rhobd_max, std::vector<int> &rhoac_bd, int rhoac_bd_max){
	//a matrix to store U9-coefficients
	Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic> U9coeffs(rhoac_max*rhobd_max*rhoab_cd_max, rhoab_max*rhocd_max*rhoac_bd_max);

	//couple the I-III and II irreps
	std::vector<int> rhomaxsac_b;
	std::vector<u3::irrep> u3ac_b = u3::decomposition(u3ac,u3b,rhomaxsac_b);

	//couple the I-II and III irreps
	std::vector<int> rhomaxsab_c;
	std::vector<u3::irrep> u3ab_c = u3::decomposition(u3ab,u3c,rhomaxsab_c);

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

	std::vector< std::vector<int> > rhoab_arr, rhocd_arr, rhoac_arr, rhobd_arr, rhoab_cd_arr, rhoac_bd_arr, rhoac_b_uarr, rhoac_b_zarr, rhoab_c_zarr, rhoab_c_uarr, rhoabc_d_arr1, rhoabc_d_arr2;
	std::vector< std::vector<double> > U6coeffs_arr1, Z6coeffs_arr, U6coeffs_arr2;

	for(int i=0; i<u3abc.size(); i++){
		std::vector<int> rab, rcd, rac, rbd, rab_cd, rac_bd, rac_b_u, rac_b_z, rab_c_z, rab_c_u, rabc_d_1, rabc_d_2;
		std::vector<double> U6coeffs_1, Z6coeffs, U6coeffs_2;

		U6coeffs_1 = u3::U6(u3ac,u3b,u3abcd,u3d,u3abc[i],rac_b_u,rabc_d_1,u3bd,rbd,rac_bd);
		Z6coeffs = u3::Z6(u3b,u3a,u3abc[i],u3c,u3ab,rab,rab_c_z,u3ac,rac,rac_b_z);
		U6coeffs_2 = u3::U6(u3ab,u3c,u3abcd,u3d,u3abc[i],rab_c_u,rabc_d_2,u3cd,rcd,rab_cd);

		rhoab_arr.push_back(rab);
		rhocd_arr.push_back(rcd);
		rhoac_arr.push_back(rac); 
		rhobd_arr.push_back(rbd); 
		rhoab_cd_arr.push_back(rab_cd); 
		rhoac_bd_arr.push_back(rac_bd); 
		rhoac_b_uarr.push_back(rac_b_u); 
		rhoac_b_zarr.push_back(rac_b_z); 
		rhoab_c_zarr.push_back(rab_c_z);
		rhoab_c_uarr.push_back(rab_c_u); 
		rhoabc_d_arr1.push_back(rabc_d_1); 
		rhoabc_d_arr2.push_back(rabc_d_2);

		U6coeffs_arr1.push_back(U6coeffs_1);
		U6coeffs_arr2.push_back(U6coeffs_2);
		Z6coeffs_arr.push_back(Z6coeffs);
/*
		for(int j=0; j<U6coeffs_1.size(); j++){
			std::cout << U6coeffs_1[j] << "\t";
		}
		std::cout << std::endl;
		for(int j=0; j<Z6coeffs.size(); j++){
			std::cout << Z6coeffs[j] << "\t";
		}
		std::cout << std::endl;
		for(int j=0; j<U6coeffs_2.size(); j++){
			std::cout << U6coeffs_2[j] << "\t";
		}
		std::cout << std::endl;
*/
	}

	int row=0, col=0;
	for(int rac=0; rac<rhoac_max; rac++){
		for(int rbd=0; rbd<rhobd_max; rbd++){
			for(int rab_cd=0; rab_cd<rhoab_cd_max; rab_cd++){
				rhoac.push_back(rac+1);
				rhobd.push_back(rbd+1);
				rhoab_cd.push_back(rab_cd+1);
			}
		}
	}

	for(int rab=0; rab<rhoab_max; rab++){
		for(int rcd=0; rcd<rhocd_max; rcd++){
			for(int rac_bd=0; rac_bd<rhoac_bd_max; rac_bd++){
				rhoab.push_back(rab+1);
				rhocd.push_back(rcd+1);
				rhoac_bd.push_back(rac_bd+1);
			}
		}
	}

	for(int m=0; m<rhoab_cd.size(); m++){ //m+1 is the multiplicity label for rhoac, rhobd, rhoab_cd
		for(int n=0; n<rhoac_bd.size(); n++){ //n+1 is the multiplicity label for rhoab, rhocd, rhoac_bd
			U9coeffs(m,n)=0.0;

			for(int i=0; i<u3abc.size(); i++){
				for(int j=0; j<U6coeffs_arr1[i].size(); j++){
					for(int k=0; k<Z6coeffs_arr[i].size(); k++){
						for(int l=0; l<U6coeffs_arr2[i].size(); l++){
							if( rhoab_arr[i][k]==rhoab[n] && rhocd_arr[i][l]==rhocd[n] && rhoac_bd_arr[i][j]==rhoac_bd[n] 
								&& rhoac_arr[i][k]==rhoac[m] && rhobd_arr[i][j]==rhobd[m] && rhoab_cd_arr[i][l]==rhoab_cd[m]
								&& rhoac_b_uarr[i][j]==rhoac_b_zarr[i][k] && rhoabc_d_arr1[i][j]==rhoabc_d_arr2[i][l] && rhoab_c_zarr[i][k]==rhoab_c_uarr[i][l] ){
								U9coeffs(m,n) += U6coeffs_arr1[i][j]*Z6coeffs_arr[i][k]*U6coeffs_arr2[i][l];
							}
						}
					}
				}
			}
		}
	}

	return U9coeffs;
}


}

#endif
