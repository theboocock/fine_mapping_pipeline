/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions using MKL library
 *
 * 2013 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

/////////////////
// data functions

void gcta::make_XMat_mkl(float *X)
{
    if(_mu.empty()) calcu_mu();
	
	cout<<"Recoding genotypes (individual major mode) ..."<<endl;
	unsigned long i=0, j=0, k=0, n=_keep.size(), m=_include.size();

    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        if(_dosage_flag){
            for(j=0; j<m; j++){
                if(_geno_dose[_keep[i]][_include[j]]<1e5){
                    if(_allele1[_include[j]]==_ref_A[_include[j]]) X[i*m+j]=_geno_dose[_keep[i]][_include[j]];
                    else X[i*m+j]=2.0-_geno_dose[_keep[i]][_include[j]];
                }
                else X[i*m+j]=1e6;
            }
            _geno_dose[i].clear();
        }
        else{
            for(j=0; j<_include.size(); j++){
                if(!_snp_1[_include[j]][_keep[i]] || _snp_2[_include[j]][_keep[i]]){
                    if(_allele1[_include[j]]==_ref_A[_include[j]]) X[i*m+j]=_snp_1[_include[j]][_keep[i]]+_snp_2[_include[j]][_keep[i]];
                    else X[i*m+j]=2.0-(_snp_1[_include[j]][_keep[i]]+_snp_2[_include[j]][_keep[i]]);
                }
                else X[i*m+j]=1e6;
            }
        }
    }
}

void gcta::std_XMat_mkl(float *X, vector<double> &sd_SNP, bool grm_xchr_flag, bool miss_with_mu, bool divid_by_std)
{
	if(_mu.empty()) calcu_mu();
	
    unsigned long i=0, j=0, k=0, n=_keep.size(), m=_include.size();
	sd_SNP.clear();
    sd_SNP.resize(m);
    if(_dosage_flag){
        #pragma omp parallel for private(i)
        for(j=0; j<m; j++){
            for(i=0; i<n; i++){
                double d_buf=X[i*m+j]-_mu[_include[j]];
                sd_SNP[j]+=d_buf*d_buf;
            }
            sd_SNP[j]/=(n-1.0);
        }
    }
    else{
        for(j=0; j<m; j++) sd_SNP[j]=_mu[_include[j]]*(1.0-0.5*_mu[_include[j]]);
    }
    if(divid_by_std){
        for(j=0; j<m; j++){
            if(fabs(sd_SNP[j])<1.0e-50) sd_SNP[j]=0.0;
            else sd_SNP[j]=sqrt(1.0/sd_SNP[j]);
        }        
    }
    
	#pragma omp parallel for private(j, k)
    for(i=0; i<n; i++){
		for(j=0; j<m; j++){
			k=i*m+j;
            if(X[k]<1e5){
                X[k]-=_mu[_include[j]];
                if(divid_by_std) X[k]*=sd_SNP[j];
            }
            else if(miss_with_mu) X[k]=0.0;
		}
	}
	
	if(!grm_xchr_flag) return;
	// for the X-chromosome
	check_sex();
	double f_buf=sqrt(0.5);
	
	#pragma omp parallel for private(j, k)
    for(i=0; i<n; i++){
        if(_sex[_keep[i]]==1){
            for(j=0; j<m; j++){
				k=i*m+j;
                if(X[k]<1e5) X[k]*=f_buf;
                else if(miss_with_mu) X[k]=0.0;
            }
        }
	}
}

/////////////////
// grm functions

void gcta::make_grm_mkl(bool grm_xchr_flag, bool inbred, bool output_bin, int grm_mtd, bool mlmassoc, bool diag_f3_flag)
{
	if(grm_xchr_flag) check_chrX();
	else check_autosome();
	
	unsigned long i=0, j=0, k=0, l=0, n=_keep.size(), m=_include.size();
	_geno_mkl=new float[n*m]; // alloc memory to X matrix
	
	make_XMat_mkl(_geno_mkl);
    vector<double> sd_SNP, sd_SNP_buf;
	if(grm_mtd==0) std_XMat_mkl(_geno_mkl, sd_SNP, grm_xchr_flag, false, true);
    else std_XMat_mkl(_geno_mkl, sd_SNP, grm_xchr_flag, false, false);
	
    if(!mlmassoc) cout<<"\nCalculating the genetic relationship matrix (GRM)"<<(grm_xchr_flag?" for the X chromosome":"")<<(_dosage_flag?" using imputed dosage data":"")<<" ... (Note: default speed-optimized mode, may use huge RAM)"<<endl;
    else cout<<"\nCalculating the genetic relationship matrix (GRM) ... "<<endl;
    
    // count the number of missing genotypes
    vector< vector<int> > miss_pos(n);
	bool * X_bool = new bool[n*m];
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            k=i*m+j;
            if(_geno_mkl[k]<1e5) X_bool[k]=true;
            else{
                _geno_mkl[k]=0.0;
                miss_pos[i].push_back(j);
                X_bool[k]=false;
            }
        }
    }
    
    // Calculate A_N matrix
	vector< vector<float> > A_N(n);
	for(i=0; i<n; i++) A_N[i].resize(n);
    #pragma omp parallel for private(j, k)
	for(i=0; i<n; i++){
		for(j=0; j<=i; j++){
			int comm=0;
			for(k=0; k<miss_pos[j].size(); k++) comm+=(int)X_bool[i*m+miss_pos[j][k]];
			A_N[i][j]=m-miss_pos[i].size()-comm;
		}
	}
    
    // Calculate sum of LD weights
    long double sum_wt=0.0, d_m=(double)m;
    if(grm_mtd==0) sum_wt=d_m;
    else if(grm_mtd==1){
        for(j=0; j<m; j++) sum_wt+=sd_SNP[j];
    }
    if(CommFunc::FloatEqual(sum_wt, 0.0)) throw("Error: the sum of the weights is zero!");
    
    #pragma omp parallel for private(j, k)
    for(i=0; i<n; i++){
        for(j=0; j<=i; j++) A_N[i][j]=A_N[i][j]*sum_wt/d_m;
    }
	
    // Calcuate WW'
	_grm_mkl=new float[n*n]; // alloc memory to A
	cblas_sgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, n, m, 1.0, _geno_mkl, m, _geno_mkl, m, 0.0, _grm_mkl, n);
    
    // re-calcuate the diagonals (Fhat3+1)
    if(diag_f3_flag){
        #pragma omp parallel for private(j,k,l)
        for(i=0; i<n; i++){
            l=i*n+i;
            _grm_mkl[l]=0.0;
            for(j=0; j<m; j++){
                k=i*m+j;
                _grm_mkl[l]+=_geno_mkl[k]*(_geno_mkl[k]+(_mu[_include[j]]-1.0)*sd_SNP[j]);
            }
        }
    }
    
    // Calculate A matrix
    #pragma omp parallel for private(j)
    for(i=0; i<n; i++){
        for(j=0; j<=i; j++){
            if(A_N[i][j]>0.0) _grm_mkl[i*n+j]/=A_N[i][j];
            else _grm_mkl[i*n+j]=0.0;
        }
    }
      
    if(inbred){
        #pragma omp parallel for private(j)
        for(i=0; i<n; i++){
            for(j=0; j<=i; j++) _grm_mkl[i*n+j]*=0.5;
        }
    }
    
    if(mlmassoc && grm_mtd==0){
        for(j=0; j<m; j++){
            if(fabs(sd_SNP[j])<1.0e-50) sd_SNP[j]=0.0;
            else sd_SNP[j]=1.0/sd_SNP[j];
        }
        #pragma omp parallel for private(j, k)
        for(i=0; i<n; i++){
            for(j=0; j<m; j++){
                k=i*m+j;
                if(_geno_mkl[k]<1e5) _geno_mkl[k]*=sd_SNP[j];
                else _geno_mkl[k]=0.0;
            }
        }
        delete[] X_bool;
    }
	else{
        // Output A_N and A
        string out_buf=_out;
        output_grm_mkl(_grm_mkl, A_N, output_bin);
        _out=out_buf;
        
        // free memory
        delete[] _geno_mkl;
        delete[] X_bool;
        delete[] _grm_mkl;
    }
}

void gcta::output_grm_mkl(float* A, vector< vector<float> > &A_N, bool output_grm_bin)
{
    unsigned long i=0, j=0, n=_keep.size();
 	string grm_file;
    
    if(output_grm_bin){
        // Save matrix A in binary file
        grm_file=_out+".grm.bin";
        fstream A_Bin(grm_file.c_str(), ios::out|ios::binary);
        if(!A_Bin) throw("Error: can not open the file ["+grm_file+"] to write.");
		int size=sizeof(float);
        for(i=0; i<n; i++){
            for(j=0; j<=i; j++) A_Bin.write((char*)&(A[i*n+j]), size);
        }
        A_Bin.close();
		cout<<"GRM of "<<n<<" individuals has been saved in the file ["+grm_file+"] (in binary format)."<<endl;
        
        string grm_N_file=_out+".grm.N.bin";
        fstream N_Bin(grm_N_file.c_str(), ios::out|ios::binary);
        if(!N_Bin) throw("Error: can not open the file ["+grm_N_file+"] to write.");
        size=sizeof(float);
        for(i=0; i<n; i++){
            for(j=0; j<=i; j++) N_Bin.write((char*)&(A_N[i][j]), size);
        }
        N_Bin.close();
        cout<<"Number of SNPs to calcuate the genetic relationship between each pair of individuals has been saved in the file ["+grm_N_file+"] (in binary format)."<<endl;
    }
    else{
        // Save A matrix in txt format
        grm_file=_out+".grm.gz";
        gzofstream zoutf;
        zoutf.open( grm_file.c_str() );
        if(!zoutf.is_open()) throw("Error: can not open the file ["+grm_file+"] to write.");
        cout<<"Saving the genetic relationship matrix to the file ["+grm_file+"] (in compressed text format)."<<endl;
        zoutf.setf(ios::scientific);
        zoutf.precision(6);
        for(i=0; i<n; i++){
            for(j=0; j<=i; j++) zoutf<<i+1<<'\t'<<j+1<<'\t'<<A_N[i][j]<<'\t'<<A[i*n+j]<<endl;
        }
        zoutf.close();
        cout<<"The genetic relationship matrix has been saved in the file ["+grm_file+"] (in compressed text format)."<<endl;
    }
	
	string famfile=_out+".grm.id";
	ofstream Fam(famfile.c_str());
	if(!Fam) throw("Error: can not open the file ["+famfile+"] to write.");
	for(i=0; i<n; i++) Fam<<_fid[_keep[i]]+"\t"+_pid[_keep[i]]<<endl;
	Fam.close();
	cout<<"IDs for the GRM file ["+grm_file+"] have been saved in the file ["+famfile+"]."<<endl;
}


///////////
// reml

bool gcta::comput_inverse_logdet_LDLT_mkl(eigenMatrix &Vi, double &logdet)
{
	unsigned long i=0, j=0, n=Vi.cols();
	double* Vi_mkl=new double[n*n];
	//float* Vi_mkl=new float[n*n];
	
#pragma omp parallel for private(j)
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Vi_mkl[i*n+j]=Vi(i,j);
		}
    }
	
	// MKL's Cholesky decomposition
	int info=0, int_n=(int)n;
	char uplo='L';
	dpotrf( &uplo, &int_n, Vi_mkl, &int_n, &info );
	//spotrf( &uplo, &n, Vi_mkl, &n, &info );
	if(info<0) throw("Error: Cholesky decomposition failed. Invalid values found in the matrix.\n");
	else if (info>0) return(false); //Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
	else{
		logdet=0.0;
		for(i=0; i<n; i++) {
			double d_buf=Vi_mkl[i*n+i];
			logdet+=log(d_buf*d_buf);
		}
		
		// Calcualte V inverse
		dpotri(&uplo, &int_n, Vi_mkl, &int_n, &info);
		//spotri( &uplo, &n, Vi_mkl, &n, &info );
		if(info<0) throw("Error: invalid values found in the varaince-covaraince (V) matrix.\n");
		else if (info>0) return(false); // Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
		else{
#pragma omp parallel for private(j)
			for(j=0; j<n; j++){
				for(i=0; i<=j; i++) Vi(i,j)=Vi(j,i)=Vi_mkl[i*n+j];
			}
		}
	}
	
	// free memory
	delete[] Vi_mkl;
	
	return true;
	
}

bool gcta::comput_inverse_logdet_LU_mkl(eigenMatrix &Vi, double &logdet)
{
	unsigned long i=0, j=0, n=Vi.cols();
	double* Vi_mkl=new double[n*n];
	
    #pragma omp parallel for private(j)
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Vi_mkl[i*n+j]=Vi(i,j);
		}
    }
    
    int N=(int)n;
    int *IPIV = new int[n+1];
    int LWORK = N*N;
    double *WORK = new double[n*n];
    int INFO;
    dgetrf(&N,&N,Vi_mkl,&N,IPIV,&INFO);
	if(INFO<0) throw("Error: LU decomposition failed. Invalid values found in the matrix.\n");
	else if (INFO>0){
        delete[] Vi_mkl;
        return(false); //Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
    }
	else{
		logdet=0.0;
		for(i=0; i<n; i++) {
			double d_buf=Vi_mkl[i*n+i];
			logdet+=log(fabs(d_buf));
		}
		
		// Calcualte V inverse
        dgetri(&N,Vi_mkl,&N,IPIV,WORK,&LWORK,&INFO);
		if(INFO<0) throw("Error: invalid values found in the varaince-covaraince (V) matrix.\n");
		else if (INFO>0) return(false); // Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
		else{
            #pragma omp parallel for private(j)
			for(j=0; j<n; j++){
				for(i=0; i<=j; i++) Vi(i,j)=Vi(j,i)=Vi_mkl[i*n+j];
			}
		}
	}
	
	// free memory
    delete[] Vi_mkl;
    delete[] IPIV;
    delete[] WORK;
	
	return true;
	
}

bool gcta::comput_inverse_logdet_LU_mkl_array(int n, float *Vi, double &logdet)
{
	unsigned long i=0, j=0;
	double* Vi_mkl=new double[n*n];
	
    #pragma omp parallel for private(j)
	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			Vi_mkl[i*n+j]=Vi[i*n+j];
		}
    }
    
    int N=(int)n;
    int *IPIV = new int[n+1];
    int LWORK = N*N;
    double *WORK = new double[n*n];
    int INFO;
    dgetrf(&N,&N,Vi_mkl,&N,IPIV,&INFO);
	if(INFO<0) throw("Error: LU decomposition failed. Invalid values found in the matrix.\n");
	else if (INFO>0){
        delete[] Vi_mkl;
        return(false); //Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
    }
	else{
		logdet=0.0;
		for(i=0; i<n; i++) {
			double d_buf=Vi_mkl[i*n+i];
			logdet+=log(fabs(d_buf));
		}
		
		// Calcualte V inverse
        dgetri(&N,Vi_mkl,&N,IPIV,WORK,&LWORK,&INFO);
		if(INFO<0) throw("Error: invalid values found in the varaince-covaraince (V) matrix.\n");
		else if (INFO>0) return(false); // Vi.diagonal()=Vi.diagonal().array()+Vi.diagonal().mean()*1e-3;
		else{
            #pragma omp parallel for private(j)
			for(j=0; j<n; j++){
				for(i=0; i<n; i++) Vi[i*n+j]=Vi_mkl[i*n+j];
			}
		}
	}
	
	// free memory
    delete[] Vi_mkl;
    delete[] IPIV;
    delete[] WORK;
	
	return true;
	
}
