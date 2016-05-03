/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * Implementations of functions for estimating the genetic relationship matrix
 *
 * 2010 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include "gcta.h"

void gcta::enable_grm_bin_flag()
{
	_grm_bin_flag=true;
}

void gcta::check_autosome()
{
    for(int i=0; i<_include.size(); i++){
        if(_chr[_include[i]]>_autosome_num) throw("Error: this option is for the autosomal SNPs only. Please check the option --autosome.");
    }
}

void gcta::check_chrX()
{
   for(int i=0; i<_include.size(); i++){
        if(_chr[_include[i]]!=(_autosome_num+1)) throw("Error: this option is for SNPs on the X chromosome only.");
    }
}

void gcta::check_sex()
{
    for(int i=0; i<_keep.size(); i++){
        if(_sex[_keep[i]]!=1 && _sex[_keep[i]]!=2) throw("Error: Sex information of the individual \""+_fid[_keep[i]]+" "+_pid[_keep[i]]+"\" is missing.\nUse --update-sex option to update the sex information of the individuals.");
    }
}

void gcta::output_grm_MatrixXf(bool output_grm_bin)
{
    int i=0, j=0;
    string grm_file;
    if(output_grm_bin){
        // Save matrix A in binary file
        grm_file=_out+".grm.bin";
        fstream A_Bin(grm_file.c_str(), ios::out|ios::binary);
        if(!A_Bin) throw("Error: can not open the file ["+grm_file+"] to write.");
		float f_buf=0.0;
		int size=sizeof(float);
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++){
				f_buf=(float)(_grm(i,j));
				A_Bin.write((char*)&f_buf, size);
			}
        }
        A_Bin.close();
        cout<<"GRM of "<<_keep.size()<<" individuals has been saved in the file ["+grm_file+"] (in binary format)."<<endl;

        string grm_N_file=_out+".grm.N.bin";
        fstream N_Bin(grm_N_file.c_str(), ios::out|ios::binary);
        if(!N_Bin) throw("Error: can not open the file ["+grm_N_file+"] to write.");
        f_buf=0.0;
        size=sizeof(int);
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++){
                f_buf=(float)(_grm_N(i,j));
                N_Bin.write((char*)&f_buf, size);
            }
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
        for(i=0; i<_keep.size(); i++){
            if(_grm_N.rows()>0) for(j=0; j<=i; j++) zoutf<<i+1<<'\t'<<j+1<<'\t'<<(int)_grm_N(i,j)<<'\t'<<_grm(i,j)<<endl;
            else for(j=0; j<=i; j++) zoutf<<i+1<<'\t'<<j+1<<"\t0\t"<<_grm(i,j)<<endl;
        }
        zoutf.close();
        cout<<"The genetic relationship matrix has been saved in the file ["+grm_file+"] (in compressed text format)."<<endl;
    }

    string famfile=_out+".grm.id";
	ofstream Fam(famfile.c_str());
	if(!Fam) throw("Error: can not open the file ["+famfile+"] to write.");
	for(i=0; i<_keep.size(); i++) Fam<<_fid[_keep[i]]+"\t"+_pid[_keep[i]]<<endl;
	Fam.close();
	cout<<"IDs for the GRM file ["+grm_file+"] have been saved in the file ["+famfile+"]."<<endl;
}

int gcta::read_grm_id(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only)
{
    // read GRM IDs
    string grm_id_file=grm_file+".grm.id";
    if(out_id_log) cout<<"Reading IDs of the GRM from ["+grm_id_file+"]."<<endl;
    ifstream i_grm_id(grm_id_file.c_str());
	if(!i_grm_id) throw("Error: can not open the file ["+grm_id_file+"] to read.");
	string str_buf, id_buf;
	vector<string> fid, pid;
	grm_id.clear();
	while(i_grm_id){
	    i_grm_id>>str_buf;
	    if(i_grm_id.eof()) break;
	    fid.push_back(str_buf);
		id_buf=str_buf+":";
	    i_grm_id>>str_buf;
	    pid.push_back(str_buf);
		id_buf+=str_buf;
		grm_id.push_back(id_buf);
		getline(i_grm_id, str_buf);
	}
	i_grm_id.close();
	int n=grm_id.size();
    if(out_id_log) cout<<n<<" IDs read from ["+grm_id_file+"]."<<endl;
    
    if(_id_map.empty()){
        _fid=fid;
        _pid=pid;
        _indi_num=_fid.size();
        _sex.resize(_fid.size());
        init_keep();
    }
    
    return(n);
}

void gcta::read_grm(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only)
{
	if(_grm_bin_flag) read_grm_bin(grm_file, grm_id, out_id_log, read_id_only);
	else read_grm_gz(grm_file, grm_id, out_id_log, read_id_only);
}

void gcta::read_grm_gz(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only)
{
    int n=read_grm_id(grm_file, grm_id, out_id_log, read_id_only);

    if(read_id_only) return;

    string grm_gzfile=grm_file+".grm.gz", str_buf;
    const int MAX_LINE_LENGTH = 1000;
    char buf[MAX_LINE_LENGTH];
    gzifstream zinf;
    zinf.open(grm_gzfile.c_str());
    if(!zinf.is_open()) throw("Error: can not open the file ["+grm_gzfile+"] to read.");

    int indx1=0, indx2=0, nline=0;
    double grm_buf=0.0, grm_N_buf;
    string errmsg="Error: failed to read ["+grm_gzfile+"]. The format of the GRM file has been changed?\nError occurs in line:\n";
    cout<<"Reading the GRM from ["+grm_gzfile+"]."<<endl;
    _grm.resize(n,n);
    _grm_N.resize(n,n);
    while(1){
        zinf.getline(buf, MAX_LINE_LENGTH, '\n');
        if(zinf.fail() || !zinf.good()) break;
        stringstream ss(buf);
        if(!(ss>>indx1)) throw(errmsg+buf);
        if(!(ss>>indx2)) throw(errmsg+buf);
        if(!(ss>>grm_N_buf)) throw(errmsg+buf);
        if(!(ss>>grm_buf)) throw(errmsg+buf);
		if(indx1 < indx2 || indx1>n || indx2>n) throw(errmsg+buf);
		if(grm_N_buf==0) cout<<"Warning: "<<buf<<endl;
		_grm_N(indx1-1,indx2-1)=_grm_N(indx2-1,indx1-1)=grm_N_buf;
		_grm(indx1-1,indx2-1)=_grm(indx2-1,indx1-1)=grm_buf;
		nline++;
        if(ss>>str_buf) throw(errmsg+buf);
    }
    zinf.close();
    cout<<"Pairwise genetic relationships between "<<n<<" individuals are included from ["+grm_gzfile+"]."<<endl;
}

void gcta::read_grm_bin(string grm_file, vector<string> &grm_id, bool out_id_log, bool read_id_only)
{
    int i=0, j=0, n=read_grm_id(grm_file, grm_id, out_id_log, read_id_only);
	
    if(read_id_only) return;
	
    string grm_binfile=grm_file+".grm.bin";
	ifstream A_bin(grm_binfile.c_str(), ios::in|ios::binary);
	if(!A_bin.is_open()) throw("Error: can not open the file ["+grm_binfile+"] to read.");
	_grm.resize(n,n);
    cout<<"Reading the GRM from ["+grm_binfile+"]."<<endl;
	int size=sizeof(float);
	float f_buf=0.0;
	for(i=0; i<n; i++){
		for(j=0; j<=i; j++){
			if(!(A_bin.read((char*)&f_buf, size))) throw("Error: the size of the ["+grm_binfile+"] file is incomplete?");
			_grm(j,i)=_grm(i,j)=f_buf;
		}
	}
	A_bin.close();
    
    string grm_Nfile=grm_file+".grm.N.bin";
    ifstream N_bin(grm_Nfile.c_str(), ios::in|ios::binary);
    if(!N_bin.is_open()) throw("Error: can not open the file ["+grm_Nfile+"] to read.");
    _grm_N.resize(n,n);
    cout<<"Reading the number of SNPs for the GRM from ["+grm_Nfile+"]."<<endl;
    size=sizeof(float);
    f_buf=0.0;
    for(i=0; i<n; i++){
        for(j=0; j<=i; j++){
            if(!(N_bin.read((char*)&f_buf, size))) throw("Error: the size of the ["+grm_Nfile+"] file is incomplete?");
            _grm_N(j,i)=_grm_N(i,j)=f_buf;
        }
    }
    N_bin.close();
	
    cout<<"Pairwise genetic relationships between "<<n<<" individuals are included from ["+grm_binfile+"]."<<endl;
}

void gcta::rm_cor_indi(double grm_cutoff)
{
    cout<<"Pruning the GRM with a cutoff of "<<grm_cutoff<<" ..."<<endl;

    int i=0, j=0, i_buf=0;

    vector<int> rm_grm_ID1, rm_grm_ID2;
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<i; j++){
            if(_grm(_keep[i],_keep[j])>grm_cutoff){
                rm_grm_ID1.push_back(_keep[i]);
                rm_grm_ID2.push_back(_keep[j]);
            }
        }
    }
    vector<int> rm_uni_ID(rm_grm_ID1);
	rm_uni_ID.insert(rm_uni_ID.end(), rm_grm_ID2.begin(), rm_grm_ID2.end());
    stable_sort(rm_uni_ID.begin(), rm_uni_ID.end());
    rm_uni_ID.erase(unique(rm_uni_ID.begin(), rm_uni_ID.end()), rm_uni_ID.end());
    map<int,int> rm_uni_ID_count;
    for(i=0; i<rm_uni_ID.size(); i++){
        i_buf=count(rm_grm_ID1.begin(), rm_grm_ID1.end(), rm_uni_ID[i])+count(rm_grm_ID2.begin(), rm_grm_ID2.end(), rm_uni_ID[i]);
        rm_uni_ID_count.insert(pair<int,int>(rm_uni_ID[i], i_buf));
    }
    map<int,int>::iterator iter1, iter2;
    for(i=0; i<rm_grm_ID1.size(); i++){
        iter1=rm_uni_ID_count.find(rm_grm_ID1[i]);
        iter2=rm_uni_ID_count.find(rm_grm_ID2[i]);
        if(iter1->second < iter2->second){
            i_buf=rm_grm_ID1[i];
            rm_grm_ID1[i]=rm_grm_ID2[i];
            rm_grm_ID2[i]=i_buf;
        }
    }
    stable_sort(rm_grm_ID1.begin(), rm_grm_ID1.end());
    rm_grm_ID1.erase(unique(rm_grm_ID1.begin(), rm_grm_ID1.end()), rm_grm_ID1.end());
    vector<string> removed_ID;
    for(i=0; i<rm_grm_ID1.size(); i++) removed_ID.push_back(_fid[rm_grm_ID1[i]]+":"+_pid[rm_grm_ID1[i]]);

    // update _keep and _id_map
    update_id_map_rm(removed_ID, _id_map, _keep);

    cout<<"After pruning the GRM, there are "<<_keep.size()<<" individuals ("<<removed_ID.size()<<" individuals removed)."<<endl;
}

void gcta::adj_grm(double adj_grm_fac)
{
    cout<<"Adjusting the GRM for sampling errors ..."<<endl;
    int i=0, j=0, n=_keep.size();
    double off_mean=0.0, diag_mean=0.0, off_var=0.0, diag_var=0.0, d_buf=0.0;
    for(i=0; i<n; i++){
        diag_mean+=_grm(_keep[i],_keep[i]);
        for(j=0; j<i; j++) off_mean+=_grm(_keep[i],_keep[j]);
    }
    diag_mean/=n;
    off_mean/=0.5*n*(n-1.0);
    for(i=0; i<n; i++){
        d_buf=_grm(_keep[i],_keep[i])-diag_mean;
        diag_var+=d_buf*d_buf;
        for(j=0; j<i; j++){
            d_buf=_grm(_keep[i],_keep[j])-off_mean;
            off_var+=d_buf*d_buf;
        }
    }
    diag_var/=n-1.0;
    off_var/=0.5*n*(n-1.0)-1.0;
    for(i=0; i<_keep.size(); i++){
        d_buf=1.0-(adj_grm_fac+1.0/_grm_N(_keep[i],_keep[i]))/diag_var;
        if(_grm(_keep[i],_keep[i])>0) _grm(_keep[i],_keep[i])=1.0+d_buf*(_grm(_keep[i],_keep[i])-1.0);
        for(j=0; j<i; j++){
            if(_grm_N(_keep[i],_keep[j])>0) _grm(_keep[i],_keep[j])*=1.0-(adj_grm_fac+1.0/_grm_N(_keep[i],_keep[j]))/off_var;
        }
    }
}

void gcta::dc(int dosage_compen)
{
    cout<<"Parameterizing the GRM under the assumption of ";
    if(dosage_compen==1) cout<<"full dosage compensation ..."<<endl;
    else if(dosage_compen==0) cout<<"no dosage compensation ..."<<endl;

    int i=0, j=0, i_buf=0;
    double c1=1.0, c2=1.0;
    if(dosage_compen==1){ c1=2.0; c2=sqrt(2.0); } // full dosage compensation
    else if(dosage_compen==0){ c1=0.5; c2=sqrt(0.5); } // on dosage compensation
    for(i=0; i<_keep.size(); i++){
        for(j=0; j<=i; j++){
            i_buf=_sex[_keep[i]]*_sex[_keep[j]];
            if(i_buf==1) _grm(i,j)*=c1;
            else if(i_buf==2) _grm(i,j)*=c2;
        }
    }
}

void gcta::manipulate_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag)
{
    int i=0, j=0;

    vector<string> grm_id;
    if(merge_grm_flag) merge_grm(grm_file);
    else read_grm(grm_file, grm_id);

    if(!keep_indi_file.empty()) keep_indi(keep_indi_file);
    if(!remove_indi_file.empty()) remove_indi(remove_indi_file);
    if(grm_cutoff>-1.0) rm_cor_indi(grm_cutoff);
    if(!sex_file.empty()) update_sex(sex_file);
    if(adj_grm_fac>-1.0) adj_grm(adj_grm_fac);
    if(dosage_compen>-1) dc(dosage_compen);
    if(grm_cutoff>-1.0 || !keep_indi_file.empty() || !remove_indi_file.empty()){
        eigenMatrix tmp(_grm);
        _grm.resize(_keep.size(),_keep.size());
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++) _grm(i,j)=tmp(_keep[i],_keep[j]);
        }
        tmp=_grm_N;
        _grm_N.resize(_keep.size(),_keep.size());
        for(i=0; i<_keep.size(); i++){
            for(j=0; j<=i; j++) _grm_N(i,j)=tmp(_keep[i],_keep[j]);
        }
    }
}

void gcta::save_grm(string grm_file, string keep_indi_file, string remove_indi_file, string sex_file, double grm_cutoff, double adj_grm_fac, int dosage_compen, bool merge_grm_flag, bool output_grm_bin)
{
    if(dosage_compen>-1) check_sex();
    manipulate_grm(grm_file, keep_indi_file, remove_indi_file, sex_file, grm_cutoff, adj_grm_fac, dosage_compen, merge_grm_flag);
    output_grm_MatrixXf(output_grm_bin);
}

void gcta::pca(string grm_file, string keep_indi_file, string remove_indi_file, double grm_cutoff, bool merge_grm_flag, int out_pc_num)
{
    manipulate_grm(grm_file, keep_indi_file, remove_indi_file, "", grm_cutoff, -2.0, -2, merge_grm_flag);
    _grm_N.resize(0,0);
    int i=0, j=0, n=_keep.size();
    cout<<"\nPerforming principal component analysis ..."<<endl;
    
    SelfAdjointEigenSolver<MatrixXd> eigensolver(_grm.cast<double>());
    MatrixXd evec=(eigensolver.eigenvectors());
    VectorXd eval=eigensolver.eigenvalues();

    string eval_file=_out+".eigenval";
    ofstream o_eval(eval_file.c_str());
    if(!o_eval) throw("Error: can not open the file ["+eval_file+"] to read.");
    for(i=_keep.size()-1; i>=0; i--) o_eval<<eval(i)<<endl;
    o_eval.close();
    cout<<"Eigenvalues of "<<n<<" individuals have been saved in ["+eval_file+"]."<<endl;
    string evec_file=_out+".eigenvec";
    ofstream o_evec(evec_file.c_str());
    if(!o_evec) throw("Error: can not open the file ["+evec_file+"] to read.");
    if(out_pc_num>n) out_pc_num=n;
    for(i=0; i<n; i++){
        o_evec<<_fid[_keep[i]]<<" "<<_pid[_keep[i]];
        for(j=n-1; j>=(n-out_pc_num); j--) o_evec<<" "<<evec(i,j);
        o_evec<<endl;
    }
    o_evec.close();
    cout<<"The first "<<out_pc_num<<" eigenvectors of "<<n<<" individuals have been saved in ["+evec_file+"]."<<endl;
}

void gcta::merge_grm(string merge_grm_file)
{
    vector<string> grm_files, grm_id;
    read_grm_filenames(merge_grm_file, grm_files);

    int f=0, i=0, j=0;
    for(f=0; f<grm_files.size(); f++){
        read_grm(grm_files[f], grm_id, false, true);
        update_id_map_kp(grm_id, _id_map, _keep);
    }
    vector<string> uni_id;
	for(i=0; i<_keep.size(); i++) uni_id.push_back(_fid[_keep[i]]+":"+_pid[_keep[i]]);
	_n=uni_id.size();
	if(_n==0) throw("Error: no individual is in common in the GRM files.");
	else cout<<_n<<" individuals in common in the GRM files."<<endl;

    vector<int> kp;
    eigenMatrix grm=eigenMatrix::Zero(_n, _n);
    eigenMatrix grm_N=eigenMatrix::Zero(_n, _n);
    for(f=0; f<grm_files.size(); f++){
        cout<<"Reading the GRM from the "<<f+1<<"th file ..."<<endl;
        read_grm(grm_files[f], grm_id);
        StrFunc::match(uni_id, grm_id, kp);
        for(i=0; i<_n; i++){
            for(j=0; j<=i; j++){
                if(kp[i]>=kp[j]){
                    grm(i,j)+=_grm(kp[i],kp[j])*_grm_N(kp[i],kp[j]);
                    grm_N(i,j)+=_grm_N(kp[i],kp[j]);
                }
                else{
                    grm(i,j)+=_grm(kp[j],kp[i])*_grm_N(kp[j],kp[i]);
                    grm_N(i,j)+=_grm_N(kp[j],kp[i]);
                }
            }
        }        
    }
    for(i=0; i<_n; i++){
        for(j=0; j<=i; j++){
            if(grm_N(i,j)==0) _grm(i,j)=0;
            else _grm(i,j)=grm(i,j)/grm_N(i,j);
            _grm_N(i,j)=grm_N(i,j);
        }
    }
    grm.resize(0,0);
    grm_N.resize(0,0);
    cout<<"\n"<<grm_files.size()<<" GRMs have been merged together."<<endl;
}

void gcta::read_grm_filenames(string merge_grm_file, vector<string> &grm_files, bool out_log)
{
    ifstream merge_grm(merge_grm_file.c_str());
    if(!merge_grm) throw("Error: can not open the file ["+merge_grm_file+"] to read.");
    string str_buf;
    grm_files.clear();
    vector<string> vs_buf;
    while(getline(merge_grm, str_buf)){
        if(!str_buf.empty()){
            if(StrFunc::split_string(str_buf, vs_buf)==1) grm_files.push_back(vs_buf[0]);
        }
    }
    if(out_log) cout<<"There are "<<grm_files.size()<<" GRM file names specified in ["+merge_grm_file+"]."<<endl;
    if(grm_files.size()>1000) throw("Error: too many GRM file names specified in ["+merge_grm_file+"]. Maximum is 1000.");
    if(grm_files.size()<1) throw("Error: no GRM file name is found in ["+merge_grm_file+"].");
}
