/*
 * GCTA: a tool for Genome-wide Complex Trait Analysis
 *
 * GCTA options
 *
 * 2010 by Jian Yang <jian.yang@uq.edu.au>
 *
 * This file is distributed under the GNU General Public
 * License, Version 2.  Please see the file COPYING for more
 * details
 */

#include <stdlib.h>
#include "gcta.h"

void option(int option_num, char* option_str[]);

int main(int argc, char* argv[])
{   
	cout<<"*******************************************************************"<<endl;
	cout<<"* Genome-wide Complex Trait Analysis (GCTA)"<<endl;
	cout<<"* version 1.21"<<endl;
	cout<<"* (C) 2010 Jian Yang, Hong Lee, Michael Goddard and Peter Visscher"<<endl;
	cout<<"* The University of Queensland"<<endl;
	cout<<"*******************************************************************"<<endl;
    
	long int time_used=0, start=time(NULL);
	time_t curr=time(0);
	cout<<"Analysis started: "<<ctime(&curr)<<endl;
	cout<<"Options:"<<endl;
	try{ option(argc, argv); }
	catch(const string &err_msg){ cerr<<"\n"<<err_msg<<endl; }
	catch(const char *err_msg){ cerr<<"\n"<<err_msg<<endl; }
	curr=time(0);
	cout<<"\nAnalysis finished: "<<ctime(&curr);
	time_used=time(NULL)-start;
	cout<<"Computational time: "<<time_used/3600<<":"<<(time_used%3600)/60<<":"<<time_used%60<<endl;

	return 0;
}

void option(int option_num, char* option_str[])
{
	int i=0, j=0;

	// OpenMP
    bool thread_flag=false;
	int thread_num=1;
	
    // raw genotype data
    string RG_fname_file="", RG_summary_file="";
    double GC_cutoff=0.7;

	// data management
	string bfile="", bfile2="", update_sex_file="", update_freq_file="", update_refA_file="", kp_indi_file="", rm_indi_file="", extract_snp_file="", exclude_snp_file="", extract_snp_name="", exclude_snp_name="", out="gcta";
	bool SNP_major=false, bfile_flag=false, make_bed_flag=false, dose_mach_flag=false, dose_beagle_flag=false, bfile2_flag=false, out_freq_flag=false, out_ssq_flag=false;
	bool ref_A=false, recode=false, recode_nomiss=false, save_ram=false, autosome_flag=false;
	int autosome_num=22, extract_chr_start=0, extract_chr_end=0;
	string dose_file="", dose_info_file="", update_impRsq_file="";
	double maf=0.0, max_maf=0.0, dose_Rsq_cutoff=0.0;

	// GRM
	bool ibc=false, ibc_all=false, grm_flag=false, grm_bin_flag=true, m_grm_flag=false, m_grm_bin_flag=true, make_grm_flag=false, make_grm_inbred_flag=false, make_grm_xchar_flag=false, grm_out_bin_flag=true, make_grm_f3_flag=false;
	bool pca_flag=false;
	double grm_adj_fac=-2.0, grm_cutoff=-2.0;
	int dosage_compen=-2, out_pc_num=20, make_grm_mtd=0;
	string grm_file="", paa_file="";

	// LD
	string LD_file="", i_ld_file="";
	bool LD=false, LD_search=false, LD_i=false;
	int LD_step=10;
	double LD_wind=2e7, LD_sig=0.05;

	// initialize paramters for simulation based on real genotype data
	bool simu_qt_flag=false, simu_cc=false, simu_emb_flag=false, simu_output_causal=false;
	int simu_rep=1, simu_case_num=0, simu_control_num=0;
	double simu_h2=0.1, simu_K=0.1, simu_gener=100, simu_seed=-CommFunc::rand_seed();
	string simu_causal="";
    
    // simulate unlinked SNPs
    bool simu_unlinked_flag=false;
    int simu_unlinked_n=1, simu_unlinked_m=1;
    double simu_unlinked_maf=0.0;

	// estimate genetic distance based on hapmap_data
	bool hapmap_genet_dst=false;
	string hapmap_genet_dst_file="";

	// REML analysis
	int mphen=1, mphen2=2, reml_mtd=0, MaxIter=100;
	double prevalence=-2.0, prevalence2=-2.0;
	bool reml_flag=false, pred_rand_eff=false, est_fix_eff=false, blup_snp_flag=false, no_constrain=false, reml_lrt_flag=false, no_lrt=false, bivar_reml_flag=false, ignore_Ce=false, within_family=false, reml_bending=false, HE_reg_flag=false, reml_diag_one=false, bivar_no_constrain=false;
	string phen_file="", qcovar_file="", covar_file="", qgxe_file="", gxe_file="", blup_indi_file="";
	vector<double> reml_priors, reml_priors_var, fixed_rg_val;
	vector<int> reml_drop;
	reml_drop.push_back(1);

	// Joint analysis of GWAS MA
	string massoc_file="", massoc_init_snplist="", massoc_cond_snplist="";
	int massoc_wind=1e7, massoc_top_SNPs=-1;
	double massoc_p=5e-8, massoc_collinear=0.9, massoc_sblup_fac=-1, massoc_gc_val=-1;
	bool massoc_slct_flag=false, massoc_joint_flag=false, massoc_sblup_flag=false, massoc_gc_flag=false, massoc_actual_geno_flag=false, massoc_backward_flag=false;
    
    // mixed linear model association 
    bool mlma_flag=false, mlma_loco_flag=false, mlma_no_adj_covar=false;
    
    int argc=option_num;
    vector<char *> argv(option_num+2);
    for(i=0; i<option_num; i++) argv[i]=option_str[i];
    argv[option_num]="gcta"; argv[option_num+1]="gcta";
	for(i=1; i<argc; i++){
		if(strcmp(argv[i],"--thread-num")==0){
			thread_num=atoi(argv[++i]);
			cout<<"--thread-num "<<thread_num<<endl;
			if(thread_num<1 || thread_num>1000) throw("\nError: --thread-num should be from 1 to 1000.\n");
		}
        // raw genotype data
		else if(strcmp(argv[i],"--raw-files")==0){
			RG_fname_file=argv[++i];
			cout<<"--raw-files "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--raw-summary")==0){
			RG_summary_file=argv[++i];
			cout<<"--raw-summary "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--gencall")==0){
			GC_cutoff=atof(argv[++i]);
			cout<<"--gencall "<<GC_cutoff<<endl;
			if(GC_cutoff<0.0 || GC_cutoff>1.0) throw("\nError: --gencall should be within the range from 0 to 1.\n");
		}

		// data management
		else if(strcmp(argv[i],"--bfile")==0){
			bfile_flag=true;
			bfile=argv[++i];
			cout<<"--bfile "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--make-bed")==0){
			make_bed_flag=true;
			cout<<"--make-bed "<<endl;
		}
		else if(strcmp(argv[i],"--bfile2")==0){
			bfile2_flag=true;
			bfile2=argv[++i];
			cout<<"--bfile2 "<<argv[i]<<endl;
		}
		else if(strcmp(argv[i],"--dosage-mach")==0){
			dose_mach_flag=true;
			dose_beagle_flag=false;
			dose_file=argv[++i];
			dose_info_file=argv[++i];
			cout<<"--dosage-mach "<<dose_file<<" "<<dose_info_file<<endl;
		}
		else if(strcmp(argv[i],"--dosage-beagle")==0){
			dose_beagle_flag=true;
			dose_mach_flag=false;
			dose_file=argv[++i];
			dose_info_file=argv[++i];
			cout<<"--dosage-beagle "<<dose_file<<" "<<dose_info_file<<endl;
		}
		else if(strcmp(argv[i],"--imput-rsq")==0){
			dose_Rsq_cutoff=atof(argv[++i]);
			cout<<"--imput-rsq "<<dose_Rsq_cutoff<<endl;
			if(dose_Rsq_cutoff<0.0 || dose_Rsq_cutoff>1.0) throw("\nError: --imput-rsq should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--update-imput-rsq")==0){
			update_impRsq_file=argv[++i];
			cout<<"--update-imput-rsq "<<update_impRsq_file<<endl;
            CommFunc::FileExist(update_impRsq_file);
		}
		else if(strcmp(argv[i],"--update-freq")==0){
			update_freq_file=argv[++i];
			cout<<"--update-freq "<<update_freq_file<<endl;
            CommFunc::FileExist(update_freq_file);
		}
		else if(strcmp(argv[i],"--update-ref-allele")==0){
			update_refA_file=argv[++i];
			cout<<"--update-ref-allele "<<update_refA_file<<endl;
            CommFunc::FileExist(update_refA_file);
		}
		else if(strcmp(argv[i],"--keep")==0){
			kp_indi_file=argv[++i];
			cout<<"--keep "<<kp_indi_file<<endl;
            CommFunc::FileExist(kp_indi_file);
		}
		else if(strcmp(argv[i],"--remove")==0){
			rm_indi_file=argv[++i];
			cout<<"--remove "<<rm_indi_file<<endl;
            CommFunc::FileExist(rm_indi_file);
		}
		else if(strcmp(argv[i],"--update-sex")==0){
			update_sex_file=argv[++i];
			cout<<"--update-sex "<<update_sex_file<<endl;
            CommFunc::FileExist(update_sex_file);
		}
		else if(strcmp(argv[i],"--chr")==0){
			extract_chr_start=extract_chr_end=atoi(argv[++i]);
			cout<<"--chr "<<extract_chr_start<<endl;
			if(extract_chr_start<1 || extract_chr_start>100) throw("\nError: --chr should be within the range from 1 to 100.\n");
		}
		else if(strcmp(argv[i],"--autosome-num")==0){
			autosome_num=atoi(argv[++i]);
			cout<<"--autosome-num "<<autosome_num<<endl;
            if(autosome_num<1 || autosome_num>100) throw("\nError: invalid number specified after the option --autosome-num.\n");
		}
		else if(strcmp(argv[i],"--autosome")==0){
		    autosome_flag=true;
			cout<<"--autosome"<<endl;
		}
		else if(strcmp(argv[i],"--extract")==0){
			extract_snp_file=argv[++i];
			cout<<"--extract "<<extract_snp_file<<endl;
            CommFunc::FileExist(extract_snp_file);
		}
		else if(strcmp(argv[i],"--exclude")==0){
			exclude_snp_file=argv[++i];
			cout<<"--exclude "<<exclude_snp_file<<endl;
            CommFunc::FileExist(exclude_snp_file);
		}
		else if(strcmp(argv[i],"--extract-snp")==0){
			extract_snp_name=argv[++i];
			cout<<"--extract-snp "<<extract_snp_name<<endl;
		}
		else if(strcmp(argv[i],"--exclude-snp")==0){
			exclude_snp_name=argv[++i];
			cout<<"--exclude-snp "<<exclude_snp_name<<endl;
		}		
        else if(strcmp(argv[i],"--maf")==0){
			maf=atof(argv[++i]);
			cout<<"--maf "<<maf<<endl;
			if(maf<0 || maf>0.5) throw("\nError: --maf should be within the range from 0 to 0.5.\n");
		}
		else if(strcmp(argv[i],"--max-maf")==0){
			max_maf=atof(argv[++i]);
			cout<<"--max-maf "<<max_maf<<endl;
			if(max_maf<=0) throw("\nError: --max-maf should be > 0.\n");
		}
		else if(strcmp(argv[i],"--out")==0){
			out=argv[++i];
			cout<<"--out "<<out<<endl;
		}
		else if(strcmp(argv[i],"--freq")==0){
			out_freq_flag=true;
            thread_flag=true;
			cout<<"--freq"<<endl;
		}
		else if(strcmp(argv[i],"--ssq")==0){
			out_ssq_flag=true;
			cout<<"--ssq"<<endl;
		}
		else if(strcmp(argv[i],"--recode")==0){
			recode=true;
            thread_flag=true;
			cout<<"--recode"<<endl;
		}
		else if(strcmp(argv[i],"--recode-nomiss")==0){
			recode_nomiss=true;
            thread_flag=true;
			cout<<"--recode-nomiss"<<endl;
		}
		else if(strcmp(argv[i],"--save-ram")==0){
			save_ram=true;
			cout<<"--save-ram"<<endl;
		}
		// GRM
		else if(strcmp(argv[i],"--paa")==0){
			paa_file=argv[++i];
			cout<<"--paa "<<paa_file<<endl;
            CommFunc::FileExist(paa_file);
		}
		else if(strcmp(argv[i],"--ibc")==0){
			ibc=true;
			cout<<"--ibc"<<endl;
		}
		else if(strcmp(argv[i],"--ibc-all")==0){
		    ibc=ibc_all=true;
			cout<<"--ibc-all"<<endl;
		}
		else if(strcmp(argv[i],"--mgrm")==0 || strcmp(argv[i],"--mgrm-bin")==0){
			m_grm_flag=true;
			grm_file=argv[++i];
			cout<<argv[i-1]<<" "<<grm_file<<endl;
		}
		else if(strcmp(argv[i],"--mgrm-gz")==0){
			m_grm_flag=true;
			m_grm_bin_flag=false;
            grm_bin_flag=false;
			grm_file=argv[++i];
			cout<<"--mgrm-gz "<<grm_file<<endl;
		}
		else if(strcmp(argv[i],"--grm")==0 || strcmp(argv[i],"--grm-bin")==0){
			grm_flag=true;
			grm_file=argv[++i];
			cout<<argv[i-1]<<" "<<grm_file<<endl;
		}
		else if(strcmp(argv[i],"--grm-gz")==0){
			grm_flag=true;
			m_grm_bin_flag=false;
            grm_bin_flag=false;
			grm_file=argv[++i];
			cout<<"--grm-gz "<<grm_file<<endl;
		}
		else if(strcmp(argv[i],"--make-grm")==0 || strcmp(argv[i],"--make-grm-bin")==0){
			make_grm_flag=true;
            thread_flag=true;
			cout<<argv[i]<<endl;
		}
        else if(strcmp(argv[i],"--make-grm-gz")==0){
		    make_grm_flag=true;
			grm_out_bin_flag=false;
            thread_flag=true;
			cout<<"--make-grm-gz"<<endl;
		}
		else if(strcmp(argv[i],"--make-grm-alg")==0){
            make_grm_flag=true;
			make_grm_mtd=atoi(argv[++i]);
            thread_flag=true;
			cout<<"--make-grm-alg "<<make_grm_mtd<<endl;
			if(make_grm_mtd<0 || make_grm_mtd>1) throw("\nError: --make-grm-alg should be 0 or 1.\n");
		}
        else if(strcmp(argv[i],"--make-grm-f3")==0){
		    make_grm_f3_flag=true;
			grm_out_bin_flag=true;
            thread_flag=true;
			cout<<"--make-grm-f3"<<endl;
		}
		else if(strcmp(argv[i],"--make-grm-xchr")==0 || strcmp(argv[i],"--make-grm-xchr-bin")==0){
            make_grm_flag=true;
			make_grm_xchar_flag=true;
            thread_flag=true;
			cout<<argv[i]<<endl;
		}
        else if(strcmp(argv[i],"--make-grm-xchr-gz")==0){
            make_grm_flag=true;
			make_grm_xchar_flag=true;
			grm_out_bin_flag=false;
            thread_flag=true;
			cout<<"--make-grm-xchr-gz"<<endl;
		}
        else if(strcmp(argv[i],"--make-grm-inbred")==0 || strcmp(argv[i],"--make-grm-inbred-bin")==0){
            make_grm_flag=true;
			make_grm_inbred_flag=true;
            thread_flag=true;
			cout<<argv[i]<<endl;
		}
        else if(strcmp(argv[i],"--make-grm-inbred-gz")==0){
            make_grm_flag=true;
			grm_out_bin_flag=false;
			make_grm_inbred_flag=true;
            thread_flag=true;
			cout<<"--make-grm-inbred-gz"<<endl;
		}
		else if(strcmp(argv[i],"--grm-adj")==0){
			grm_adj_fac=atof(argv[++i]);
			cout<<"--grm-adj "<<grm_adj_fac<<endl;
			if(grm_adj_fac<0 || grm_adj_fac>1) throw("\nError: the value to be specified after --grm-adj should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--dc")==0){
			dosage_compen=atoi(argv[++i]);
			cout<<"--dc "<<dosage_compen<<endl;
			if(dosage_compen!=0 && dosage_compen!=1) throw("\nError: the value to be specified after --dc should be 0 or 1.\n");
		}
		else if(strcmp(argv[i],"--grm-cutoff")==0){
			grm_cutoff=atof(argv[++i]);
			if(grm_cutoff>=-1 && grm_cutoff<=2) cout<<"--grm-cutoff "<<grm_cutoff<<endl;
            else grm_cutoff=-2;
		}
		else if(strcmp(argv[i],"--pca")==0){
		    pca_flag=true;
            thread_flag=true;
            i++;
            if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) { out_pc_num=20; i--; }
            else out_pc_num=atoi(argv[i]);
			cout<<"--pca "<<out_pc_num<<endl;
			if(out_pc_num<1) throw("\nError: the value to be specified after --pca should be positive.\n");
		}
		// estimation of LD structure
		else if(strcmp(argv[i],"--ld")==0){
			LD=true;
			LD_file=argv[++i];
			cout<<"--ld "<<LD_file<<endl;
            CommFunc::FileExist(LD_file);
		}
		else if(strcmp(argv[i],"--ld-step")==0){
			LD_search=true;
			LD_step=atoi(argv[++i]);
			cout<<"--ld-step "<<LD_step<<endl;
			if(LD_step<1 || LD_step>20) throw("\nError: --ld-step should be within the range from 1 to 20.\n");
		}
		else if(strcmp(argv[i],"--ld-wind")==0 || strcmp(argv[i],"--ld-pruning-wind")==0 || strcmp(argv[i],"--make-grm-wt-wind")==0){
			LD_wind=atof(argv[++i]);
			cout<<argv[i-1]<<" "<<LD_wind<<endl;
			LD_wind*=1000;
			if(LD_wind<1e3 || LD_wind>2e7){
                stringstream err_msg;
                err_msg<<"\nError: "<<argv[i-1]<<" should be 1Kb or 20Mb.\n";
                throw(err_msg.str());
            }
		}
		else if(strcmp(argv[i],"--ld-sig")==0){
			LD_sig=atof(argv[++i]);
			cout<<"--ld-sig "<<LD_sig<<endl;
			if(LD_sig<=0) throw("\nError: --ld-sig should be > 0.\n");
		}
		else if(strcmp(argv[i],"--ld-i")==0){
			LD_i=true;
			cout<<"--ld-i"<<endl;
		}
        // simulation based on real genotype data
		else if(strcmp(argv[i],"--simu-qt")==0){
			simu_qt_flag=true;
			cout<<"--simu-qt"<<endl;
		}
		else if(strcmp(argv[i],"--simu-cc")==0){
			simu_cc=true;
			simu_case_num=atoi(argv[++i]);
			simu_control_num=atoi(argv[++i]);
			cout<<"--simu-cc "<<simu_case_num<<" "<<simu_control_num<<endl;
            if(simu_case_num<10) throw("Error: --simu-cc, Invalid number of cases. Minimun number 10.");
            if(simu_control_num<10) throw("Error: --simu-cc, Invalid number of controls. Minimum number 10.");
		}
		else if(strcmp(argv[i],"--simu-rep")==0){
			simu_rep=atoi(argv[++i]);
			cout<<"--simu-rep "<<simu_rep<<endl;
            if(simu_rep<1 || simu_rep>10000) throw("Error: --simu-rep should be within the range from 1 to 10000.");
		}
		else if(strcmp(argv[i],"--simu-hsq")==0){
			simu_h2=atof(argv[++i]);
			cout<<"--simu-hsq "<<simu_h2<<endl;
            if(simu_h2>1.0 || simu_h2<0.0) throw("Error: --simu-h2 should be within the range from 0 to 1.");
		}
		else if(strcmp(argv[i],"--simu-k")==0){
			simu_K=atof(argv[++i]);
			cout<<"--simu-k "<<simu_K<<endl;
            if(simu_K>0.5 || simu_K<0.0001) throw("Error: --simu-K should be within the range from 0.0001 to 0.5.");
		}
		else if(strcmp(argv[i],"--simu-causal-loci")==0){
			simu_causal=argv[++i];
			cout<<"--simu-causal-loci "<<simu_causal<<endl;
            CommFunc::FileExist(simu_causal);
		}
		else if(strcmp(argv[i],"--simu-embayesb")==0){ // internal
			simu_emb_flag=true;
			cout<<"--simu-embayesb"<<endl;
		}
		else if(strcmp(argv[i],"--simu-ouput-causal")==0){ // internal
			simu_output_causal=true;
			cout<<"--simu-output-causal"<<endl;
		}
   		else if(strcmp(argv[i],"--simu-seed")==0){
			simu_seed=atof(argv[++i]);
			cout<<"--simu-seed "<<simu_seed<<endl;
            if(simu_seed<=100) throw("Error: --simu-seed should be >100.");
		}     
/*		else if(strcmp(argv[i],"--simu-gener")==0){
			simu_gener=atoi(argv[++i]);
            if(simu_gener<0 || simu_gener>1e5) throw("Error: --simu-gener should be within the range from 0 to 100000.");
			cout<<"--simu-gener "<<simu_gener<<endl;
		}
		*/
        // simulate unlinked SNPs
        /*else if(strcmp(argv[i],"--simu-unlinked")==0){
			simu_unlinked_flag=true;
			cout<<"--simu-unlinked"<<endl;
		}
        else if(strcmp(argv[i],"--simu-unlinked-n")==0){
			simu_unlinked_n=atoi(argv[++i]);
			cout<<"--simu-unlinked-n "<<simu_unlinked_n<<endl;
            if(simu_unlinked_n<1 || simu_unlinked_n>100000) throw("Error: --simu-unlinked-n should be within the range from 1 to 100,000.");
		}
        else if(strcmp(argv[i],"--simu-unlinked-m")==0){
			simu_unlinked_m=atoi(argv[++i]);
			cout<<"--simu-unlinked-m "<<simu_unlinked_m<<endl;
            if(simu_unlinked_m<1 || simu_unlinked_m>10000000) throw("Error: --simu-unlinked-m should be within the range from 1 to 10,000,000.");
		}
        else if(strcmp(argv[i],"--simu-unlinked-maf")==0){
			simu_unlinked_maf=atof(argv[++i]);
			cout<<"--simu-unlinked-maf "<<simu_unlinked_maf<<endl;
            if(simu_unlinked_maf<0 || simu_unlinked_maf>0.5) throw("Error: --simu-unlinked-maf should be within the range from 0 to 0.5.");
		}*/
        
		// calculate genetic dst based on HapMap data
		else if(strcmp(argv[i],"--hapmap-genet-dst")==0){
			hapmap_genet_dst=true;
			hapmap_genet_dst_file=argv[++i];
			cout<<"--hapmap-genet-dst "<<hapmap_genet_dst_file<<endl;
		}
		// estimate variance explained by all SNPs
		else if(strcmp(argv[i],"--HEreg")==0){
			HE_reg_flag=true;
            thread_flag=true;
			cout<<"--HEreg"<<endl;
		}
		else if(strcmp(argv[i],"--reml")==0){
			reml_flag=true;
            thread_flag=true;
			cout<<"--reml"<<endl;
            if(m_grm_flag) no_lrt=true;
		}
		else if(strcmp(argv[i],"--prevalence")==0){
			prevalence=atof(argv[++i]);
			cout<<"--prevalence "<<prevalence<<endl;
			if(prevalence<=0 || prevalence>=1) throw("\nError: --prevalence should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--reml-pred-rand")==0){
			pred_rand_eff=true;
			cout<<"--reml-pred-rand"<<endl;
		}
		else if(strcmp(argv[i],"--reml-est-fix")==0){
			est_fix_eff=true;
			cout<<"--reml-est-fix"<<endl;
		}
		else if(strcmp(argv[i],"--reml-alg")==0){
			reml_mtd=atoi(argv[++i]);
			cout<<"--reml-alg "<<reml_mtd<<endl;
			if(reml_mtd<0 || reml_mtd>2) throw("\nError: --reml-alg should be 0, 1 or 2.\n");
		}
		else if(strcmp(argv[i],"--reml-no-constrain")==0){
            reml_flag=true;
			no_constrain=true;
			cout<<"--reml-no-constrain"<<endl;
		}
		else if(strcmp(argv[i],"--reml-priors")==0){
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        reml_priors.push_back(atof(argv[i]));
		    }
		    i--;
			cout<<"--reml-priors ";
			bool err_flag=false;
			for(j=0; j<reml_priors.size(); j++){
			    cout<<reml_priors[j]<<" ";
			    if(reml_priors[j]>1.0 || reml_priors[j]<0.0) err_flag=true;
            }
			cout<<endl;
			if(err_flag || reml_priors.empty()) throw("\nError: --reml-priors. Prior values should be within the range from 0 to 1.\n");
		}
		else if(strcmp(argv[i],"--reml-priors-var")==0){
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        reml_priors_var.push_back(atof(argv[i]));
		    }
		    i--;
			cout<<"--reml-priors-var ";
			bool err_flag=false;
			for(j=0; j<reml_priors_var.size(); j++){
			    cout<<reml_priors_var[j]<<" ";
			    if(reml_priors_var[j]<0.0) err_flag=true;
            }
			cout<<endl;
			if(err_flag || reml_priors_var.empty()) throw("\nError: --reml-priors-var. Prior values should be positive.\n");
		}
		else if(strcmp(argv[i],"--reml-no-lrt")==0){
			no_lrt=true;
			cout<<"--reml-no-lrt"<<endl;
		}
		else if(strcmp(argv[i],"--reml-lrt")==0){
            no_lrt=false;
            reml_lrt_flag=true;
		    reml_drop.clear();
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        reml_drop.push_back(atoi(argv[i]));
		    }
		    i--;
			cout<<"--reml-lrt ";
			bool err_flag=false;
			for(j=0; j<reml_drop.size(); j++){
			    cout<<reml_drop[j]<<" ";
			    if(reml_drop[j]<1) err_flag=true;
            }
			cout<<endl;
			if(err_flag || reml_drop.empty()) throw("\nError: invalid values specified after --reml-lrt.\n");
		}
		else if(strcmp(argv[i],"--reml-maxit")==0){
			MaxIter=atoi(argv[++i]);
			cout<<"--reml-maxit "<<MaxIter<<endl;
			if(MaxIter<1 || MaxIter>10000) throw("\nError: --reml-maxit should be within the range from 1 to 10000.\n");
		}
		else if(strcmp(argv[i],"--reml-bending")==0){
			reml_bending=true;
			cout<<"--reml-bending "<<endl;
		}
		else if(strcmp(argv[i],"--reml-diag-one")==0){
			reml_diag_one=true;
			cout<<"--reml-diag-one "<<endl;
		}
		else if(strcmp(argv[i],"--pheno")==0){
			phen_file=argv[++i];
			cout<<"--pheno "<<phen_file<<endl;
            CommFunc::FileExist(phen_file);
		}
		else if(strcmp(argv[i],"--mpheno")==0){
			mphen=atoi(argv[++i]);
			cout<<"--mpheno "<<mphen<<endl;
			if(mphen<1) throw("Error: --mpheno should be > 0.");
		}
		else if(strcmp(argv[i],"--qcovar")==0){
			qcovar_file=argv[++i];
			cout<<"--qcovar "<<qcovar_file<<endl;
            CommFunc::FileExist(qcovar_file);
		}
		else if(strcmp(argv[i],"--covar")==0){
			covar_file=argv[++i];
			cout<<"--covar "<<covar_file<<endl;
            CommFunc::FileExist(covar_file);
		}
		else if(strcmp(argv[i],"--gxqe")==0){
			qgxe_file=argv[++i];
			cout<<"--gxqe "<<qgxe_file<<endl;
            CommFunc::FileExist(qgxe_file);
		}
		else if(strcmp(argv[i],"--gxe")==0){
			gxe_file=argv[++i];
			cout<<"--gxe "<<gxe_file<<endl;
            CommFunc::FileExist(gxe_file);
		}
		else if(strcmp(argv[i],"--blup-snp")==0){
		    blup_snp_flag=true;
			blup_indi_file=argv[++i];
			cout<<"--blup-snp "<<blup_indi_file<<endl;
            CommFunc::FileExist(blup_indi_file);
		}
		else if(strcmp(argv[i],"--reml-wfam")==0){
            reml_flag=true;
		    within_family=true;
			cout<<"--reml-within-family "<<endl;
		}
        else if(strcmp(argv[i],"--reml-bivar")==0){
			bivar_reml_flag=true;
            thread_flag=true;
            vector<int> mphen_buf;
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        mphen_buf.push_back(atoi(argv[i]));
		    }
		    i--;
            if(mphen_buf.size()<2 && mphen_buf.size()>0) throw("\nError: --reml-bivar. Please specify two traits for the bivariate REML analysis.");
            if(mphen_buf.size()==0){ mphen=1; mphen2=2; }
            else { mphen=mphen_buf[0]; mphen2=mphen_buf[1]; }
            if(mphen<1 || mphen2<1 || mphen==mphen2) throw("\nError: --reml-bivar. Invalid input parameters.");
			cout<<"--reml-bivar "<<mphen<<" "<<mphen2<<endl;
		} 
        else if(strcmp(argv[i],"--reml-bivar-prevalence")==0){
            vector<double> K_buf;
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        K_buf.push_back(atof(argv[i]));
		    }
		    i--;
            if(K_buf.size()<1 || K_buf.size()>2) throw("\nError: --reml-bivar-prevalence. Please specify the prevalences of the two diseases.");
            if(K_buf.size()==2){
                if(K_buf[0]<0.0 || K_buf[0]>1.0 || K_buf[1]<0.0 || K_buf[1]>1.0) throw("\nError: --reml-bivar-prevalence. Disease prevalence should be betwen 0 and 1.");
                cout<<"--reml-bivar-prevalence "<<K_buf[0]<<" "<<K_buf[1]<<endl;
                prevalence=K_buf[0];
                prevalence2=K_buf[1];
            }
            else{
                if(K_buf[0]<0.0 || K_buf[0]>1.0) throw("\nError: --reml-bivar-prevalence. Disease prevalence should be betwen 0 and 1.");
                cout<<"--reml-bivar-prevalence "<<K_buf[0]<<endl;
                prevalence=prevalence2=K_buf[0];
            }
		}
		else if(strcmp(argv[i],"--reml-bivar-nocove")==0){
			ignore_Ce=true;
			cout<<"--reml-bivar-nocove"<<endl;
		}
        else if(strcmp(argv[i],"--reml-bivar-lrt-rg")==0){
		    while(1){
		        i++;
		        if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0) break;
		        fixed_rg_val.push_back(atof(argv[i]));
		    }
		    i--;
			cout<<"--reml-bivar-lrt-rg ";
			bool err_flag=false;
			for(j=0; j<fixed_rg_val.size(); j++){
			    cout<<fixed_rg_val[j]<<" ";
			    if(fixed_rg_val[j]>1.0 || fixed_rg_val[j]<-1.0) err_flag=true;
            }
			cout<<endl;
			if(err_flag || fixed_rg_val.empty()) throw("\nError: --reml-bivar-lrt-rg. Any input paramter should be within the range from -1 to 1.\n");
            bool haveZero=false;
            if(CommFunc::FloatEqual(fixed_rg_val[0], 0.0)) haveZero=true;
            for(j=1; j<fixed_rg_val.size(); j++){
                if((CommFunc::FloatNotEqual(fixed_rg_val[0], 0.0) && haveZero) || (CommFunc::FloatEqual(fixed_rg_val[0], 0.0) && !haveZero)) throw("\nError: --reml-bivar-lrt-rg. Input paramters should be all zero or all non-zero values.\n");
            }
		}
        else if(strcmp(argv[i],"--reml-bivar-no-constrain")==0){
			bivar_no_constrain=true;
			cout<<"--reml-bivar-no-constrain"<<endl;
		}
		else if(strcmp(argv[i],"--massoc-file")==0){
			massoc_file=argv[++i];
			cout<<"--massoc-file "<<massoc_file<<endl;
            CommFunc::FileExist(massoc_file);
		}
		else if(strcmp(argv[i],"--massoc-slct")==0){
			massoc_slct_flag=true;
			cout<<"--massoc-slct"<<endl;
		}
        else if(strcmp(argv[i],"--massoc-top-SNPs")==0){
			massoc_slct_flag=true;
			massoc_top_SNPs=atoi(argv[++i]);
			cout<<"--massoc-top-SNPs "<<massoc_top_SNPs<<endl;
			if(massoc_top_SNPs<1 || massoc_top_SNPs>10000) throw("\nError: --massoc-top-SNPs should be within the range from 1 to 10000.\n");
		}
		else if(strcmp(argv[i],"--massoc-actual-geno")==0){
			massoc_actual_geno_flag=true;
			cout<<"--massoc-actual-geno"<<endl;
		}
		else if(strcmp(argv[i],"--massoc-p")==0){
			massoc_p=atof(argv[++i]);
			cout<<"--massoc-p "<<massoc_p<<endl;
			if(massoc_p>0.05 || massoc_p<=0) throw("\nError: --massoc-p should be within the range from 0 to 0.05.\n");
		}
		else if(strcmp(argv[i],"--massoc-collinear")==0){
			massoc_collinear=atof(argv[++i]);
			cout<<"--massoc-collinear "<<massoc_collinear<<endl;
			if(massoc_collinear>0.99 || massoc_collinear<0.01) throw("\nError: --massoc-collinear should be within the ragne from 0.01 to 0.99.\n");
		}
		else if(strcmp(argv[i],"--massoc-wind")==0){
			massoc_wind=atoi(argv[++i]);
			cout<<"--massoc-wind "<<massoc_wind<<endl;
			if(massoc_wind<100 || massoc_wind>100000) throw("\nError: invalid value for --massoc-wind. Valid range: 100 ~ 100000\n");
			massoc_wind*=1000;
		}
		else if(strcmp(argv[i],"--massoc-joint")==0){
			massoc_joint_flag=true;
			cout<<"--massoc-joint"<<endl;
		}
		else if(strcmp(argv[i],"--massoc-backward")==0){
            massoc_backward_flag=true;
			cout<<"--massoc-backward"<<endl;
		}
        else if(strcmp(argv[i],"--massoc-cond")==0){
			massoc_cond_snplist=argv[++i];
			cout<<"--massoc-cond "<<massoc_cond_snplist<<endl;
		}
        else if(strcmp(argv[i],"--massoc-gc")==0){
		    massoc_gc_flag=true;
            i++;
            if(strcmp(argv[i],"gcta")==0 || strncmp(argv[i], "--", 2)==0){
                massoc_gc_val=-1;
                i--;
            }
            else{
                massoc_gc_val=atof(argv[i]);
                if(massoc_gc_val<1 || massoc_gc_val>10) throw("\nError: invalid value specified after --massoc-gc.\n");
            }
			cout<<"--massoc-gc "<<((massoc_gc_val<0)?"":argv[i])<<endl;
		}
		else if(strcmp(argv[i],"--massoc-sblup")==0){
		    massoc_sblup_flag=true;
		    massoc_sblup_fac=atof(argv[++i]);
			cout<<"--massoc-sblup "<<massoc_sblup_fac<<endl;
			if(massoc_sblup_fac < 0) throw("\nError: invalid value for --massoc-sblup.\n");
		}
        else if(strcmp(argv[i],"--mlma")==0){
            reml_flag=false;
		    mlma_flag=true;
            thread_flag=true;
			cout<<"--mlma"<<endl;
		}
        else if(strcmp(argv[i],"--mlma-loco")==0){
            reml_flag=false;
		    mlma_loco_flag=true;
            thread_flag=true;
			cout<<"--mlma-loco "<<endl;
		}
        else if(strcmp(argv[i],"--mlma-no-adj-covar")==0){
            mlma_no_adj_covar=true;
			cout<<"--mlma-no-adj-covar "<<endl;
		}
		else if(strcmp(argv[i],"gcta")==0) break;
		else{ stringstream errmsg; errmsg<<"\nError: invalid option \""<<argv[i]<<"\".\n"; throw(errmsg.str()); }
		// genome() function
/*		else if(strcmp(argv[i],"--genome-pop")==0){
			genome_flag=true;
			genome_nSubPOP=atoi(argv[i+1]);
			genome_nSubSample.clear();
			genome_nSubSample.resize(genome_nSubPOP);
			for(j=0; j<genome_nSubPOP; j++) genome_nSubSample[j]=atoi(argv[i+2+j]);
		}
		else if(strcmp(argv[i],"--genome-N")==0) genome_popSize.assign(argv[i+1]);
		else if(strcmp(argv[i],"--genome-c")==0) genome_numIndepRegion=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-pieces")==0) genome_numPieces=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-len")==0) genome_pieceLen=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-s")==0) genome_SNP=atoi(argv[i+1]);
		else if(strcmp(argv[i],"--genome-rec")==0) genome_rec.assign(argv[i+1]);
		else if(strcmp(argv[i],"--genome-mut")==0) genome_mut=atof(argv[i+1]);
		else if(strcmp(argv[i],"--genome-mig")==0) genome_mig=atof(argv[i+1]);
*/
	}
	// conflicted options
	cout<<endl;
	if(bfile2_flag && !bfile_flag) throw("Error: the option --bfile2 should always go with the option --bfile.");
	if(m_grm_flag){
	    if(grm_flag){ grm_flag=false; cout<<"Warning: --grm option suppressed by the --mgrm option."<<endl; }
	    if(grm_cutoff>-1.0){ grm_cutoff=-2.0; cout<<"Warning: --grm-cutoff option suppressed by the --mgrm option."<<endl; }
	}
	if(pca_flag){
	    if(grm_adj_fac>-1.0){ grm_adj_fac=-2.0; cout<<"Warning: --grm-adj option suppressed by the --pca option."<<endl; }
	    else if(dosage_compen>-1){ grm_adj_fac=-2; cout<<"Warning: --dosage-compen option suppressed by the --pca option."<<endl; }
	}
    if(!gxe_file.empty() && !grm_flag && !m_grm_flag){
        cout<<"Warning: --gxe option is ignored because there is no --grm or --mgrm option specified."<<endl;
        gxe_file="";
    }
    if(pred_rand_eff && !grm_flag && !m_grm_flag){
        cout<<"Warning: --reml-pred-rand option is ignored because there is no --grm or --mgrm option specified."<<endl;
        pred_rand_eff=false;
    }
    if(dosage_compen>-1 && update_sex_file.empty()) throw("Error: you need to specify the sex information for the individuals by the option --update-sex because of the option --dc.");
    if(bfile2_flag && update_freq_file.empty()) throw("Error: you need to update the allele frequency by the option --update-freq because there are two datasets.");
    if(mlma_flag || mlma_loco_flag){
        if(!gxe_file.empty()) cout<<"Warning: the option --gxe option is disabled in this analysis."<<endl;
        if(!update_sex_file.empty()) cout<<"Warning: the option --update-sex option is disabled in this analysis."<<endl;
        if(grm_adj_fac>-1.0) cout<<"Warning: the option --grm-adj option is disabled in this analysis."<<endl;
        if(dosage_compen>-1.0) cout<<"Warning: the option --dc option is disabled in this analysis."<<endl;   
        if(est_fix_eff) cout<<"Warning: the option --reml-est-fix option is disabled in this analysis."<<endl; 
        if(pred_rand_eff) cout<<"Warning: the option --reml-pred-rand option is disabled in this analysis."<<endl; 
        if(reml_mtd!=0) cout<<"Warning: the option --reml-alg option is disabled in this analysis. The default algorithm AI-REML is used."<<endl;
        if(reml_lrt_flag) cout<<"Warning: the option --reml-lrt option is disabled in this analysis."<<endl; 
    }
	
	// OpenMP
    stringstream ss;
    ss<<thread_num;
    setenv("OMP_NUM_THREADS", ss.str().c_str(), 1);
	omp_set_num_threads(thread_num);
    if(thread_flag){
        if(thread_num==1) cout<<"Note: This is a multi-thread program. You could specify the number of threads by the --thread-num option to speed up the computation if there are multiple processors in your machine."<<endl;
        else cout<<"Note: the program will be running on "<<thread_num<<" threads."<<endl;
    }
		
    // set autosome
    if(autosome_flag){
        extract_chr_start=1;
        extract_chr_end=autosome_num;
    }
    if(make_grm_xchar_flag) extract_chr_start=extract_chr_end=(autosome_num+1);

	// Implement
	cout<<endl;
    gcta *pter_gcta=new gcta(autosome_num, out);//, *pter_gcta2=new gcta(autosome_num, rm_high_ld_cutoff, out);
	if(grm_bin_flag || m_grm_bin_flag) pter_gcta->enable_grm_bin_flag();
    //if(simu_unlinked_flag) pter_gcta->simu_geno_unlinked(simu_unlinked_n, simu_unlinked_m, simu_unlinked_maf);
    if(!RG_fname_file.empty()){
		if(RG_summary_file.empty()) throw("Error: please input the summary information for the raw data files by the option --raw-summary.");
        pter_gcta->read_IRG_fnames(RG_summary_file, RG_fname_file, GC_cutoff);
    }
    else if(bfile_flag || mlma_flag){
		if(hapmap_genet_dst) pter_gcta->genet_dst(bfile, hapmap_genet_dst_file);
		else{
		    if(bfile2_flag){
		        cout<<"There are two datasets specified (in PLINK binary PED format).\nReading dataset 1 ..."<<endl;
		        if(update_freq_file.empty()) throw("Error: since there are two dataset, you should update the allele frequencies that are calculated in the combined dataset.");
		    }
			pter_gcta->read_famfile(bfile+".fam");
            if(!kp_indi_file.empty()) pter_gcta->keep_indi(kp_indi_file);
			if(!rm_indi_file.empty()) pter_gcta->remove_indi(rm_indi_file);
			if(!update_sex_file.empty()) pter_gcta->update_sex(update_sex_file);
			if(!blup_indi_file.empty()) pter_gcta->read_indi_blup(blup_indi_file);
			pter_gcta->read_bimfile(bfile+".bim");
			if(!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
			if(!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
			if(extract_chr_start>0) pter_gcta->extract_chr(extract_chr_start, extract_chr_end);
			if(!extract_snp_name.empty()) pter_gcta->extract_single_snp(extract_snp_name);
			if(!exclude_snp_name.empty()) pter_gcta->exclude_single_snp(exclude_snp_name);
			if(!update_refA_file.empty()) pter_gcta->update_ref_A(update_refA_file);
			if(LD) pter_gcta->read_LD_target_SNPs(LD_file);
			pter_gcta->read_bedfile(bfile+".bed");
			if(!update_impRsq_file.empty()) pter_gcta->update_impRsq(update_impRsq_file);
			if(!update_freq_file.empty()) pter_gcta->update_freq(update_freq_file);
            if(dose_Rsq_cutoff>0.0) pter_gcta->filter_impRsq(dose_Rsq_cutoff);
			if(maf>0) pter_gcta->filter_snp_maf(maf);
			if(max_maf>0.0) pter_gcta->filter_snp_max_maf(max_maf);
			if(out_freq_flag) pter_gcta->save_freq(out_ssq_flag);
			else if(!paa_file.empty()) pter_gcta->paa(paa_file);
			else if(ibc) pter_gcta->ibc(ibc_all);
            else if(make_grm_flag) pter_gcta->make_grm_mkl(make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, false, make_grm_f3_flag);
			else if(recode || recode_nomiss) pter_gcta->save_XMat(recode_nomiss);
			else if(LD) pter_gcta->LD_Blocks(LD_step, LD_wind, LD_sig, LD_i, save_ram);
			else if(blup_snp_flag) pter_gcta->blup_snp_geno();
            else if(mlma_flag) pter_gcta->mlma(grm_file, phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
            else if(mlma_loco_flag) pter_gcta->mlma_loco(phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
			else if(massoc_slct_flag | massoc_joint_flag | massoc_backward_flag) pter_gcta->run_massoc_slct(massoc_file, massoc_wind, massoc_p, massoc_collinear, massoc_top_SNPs, massoc_joint_flag, massoc_gc_flag, massoc_gc_val, massoc_actual_geno_flag, massoc_backward_flag);
			else if(!massoc_cond_snplist.empty()) pter_gcta->run_massoc_cond(massoc_file, massoc_cond_snplist, massoc_wind, massoc_collinear, massoc_gc_flag, massoc_gc_val, massoc_actual_geno_flag);
			else if(massoc_sblup_flag) pter_gcta->run_massoc_sblup(massoc_file, massoc_wind, massoc_sblup_fac);
            else if(simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_seed, simu_output_causal, simu_emb_flag);
			else if(make_bed_flag) pter_gcta->save_plink();
		}
	}
	else if(dose_beagle_flag || dose_mach_flag){
        if(massoc_slct_flag | massoc_joint_flag | !massoc_cond_snplist.empty()) throw("Error: the --dosage option can't be used in combined with the --massoc options.");
		if(dose_mach_flag) pter_gcta->read_imp_info_mach(dose_info_file);
		else if(dose_beagle_flag) pter_gcta->read_imp_info_beagle(dose_info_file);
        if(!extract_snp_file.empty()) pter_gcta->extract_snp(extract_snp_file);
        if(!exclude_snp_file.empty()) pter_gcta->exclude_snp(exclude_snp_file);
        if(!extract_snp_name.empty()) pter_gcta->extract_single_snp(extract_snp_name);
        if(!exclude_snp_name.empty()) pter_gcta->exclude_single_snp(exclude_snp_name);
        if(extract_chr_start>0) cout<<"Warning: the option --chr, --autosome or --nonautosome is inactive for dosage data."<<endl;
        if(!update_refA_file.empty()) pter_gcta->update_ref_A(update_refA_file);
        if(dose_mach_flag) pter_gcta->read_imp_dose_mach(dose_file, kp_indi_file, rm_indi_file, blup_indi_file);
		else if(dose_beagle_flag) pter_gcta->read_imp_dose_beagle(dose_file, kp_indi_file, rm_indi_file, blup_indi_file);
		if(!update_sex_file.empty()) pter_gcta->update_sex(update_sex_file);
		if(!update_impRsq_file.empty()) pter_gcta->update_impRsq(update_impRsq_file);
		if(!update_freq_file.empty()) pter_gcta->update_freq(update_freq_file);
		if(dose_Rsq_cutoff>0.0) pter_gcta->filter_impRsq(dose_Rsq_cutoff);
		if(maf>0.0) pter_gcta->filter_snp_maf(maf);
		if(max_maf>0.0) pter_gcta->filter_snp_max_maf(max_maf);
		if(out_freq_flag) pter_gcta->save_freq(out_ssq_flag);
        else if(make_grm_flag) pter_gcta->make_grm_mkl( make_grm_xchar_flag, make_grm_inbred_flag, grm_out_bin_flag, make_grm_mtd, false, make_grm_f3_flag);
		else if(recode || recode_nomiss) pter_gcta->save_XMat(recode_nomiss);
		else if(blup_snp_flag) pter_gcta->blup_snp_dosage();
        else if(massoc_sblup_flag) pter_gcta->run_massoc_sblup(massoc_file, massoc_wind, massoc_sblup_fac);
        else if(simu_qt_flag || simu_cc) pter_gcta->GWAS_simu(bfile, simu_rep, simu_causal, simu_case_num, simu_control_num, simu_h2, simu_K, simu_seed, simu_output_causal, simu_emb_flag);
		else if(make_bed_flag) pter_gcta->save_plink();        
        else if(mlma_flag) pter_gcta->mlma(grm_file, phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
        else if(mlma_loco_flag) pter_gcta->mlma_loco(phen_file, qcovar_file, covar_file, mphen, MaxIter, reml_priors, reml_priors_var, no_constrain, make_grm_inbred_flag, mlma_no_adj_covar);
	}
    else if(HE_reg_flag) pter_gcta->HE_reg(grm_file, phen_file, kp_indi_file, rm_indi_file, mphen);
	else if((reml_flag || bivar_reml_flag) && phen_file.empty()) throw("\nError: phenotype file is required for reml analysis.\n");
    else if(bivar_reml_flag){
		pter_gcta->fit_bivar_reml(grm_file, phen_file, qcovar_file, covar_file, kp_indi_file, rm_indi_file, update_sex_file, mphen, mphen2, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, pred_rand_eff, est_fix_eff, reml_mtd, MaxIter, reml_priors, reml_priors_var, reml_drop, no_lrt, prevalence, prevalence2, no_constrain, ignore_Ce, fixed_rg_val, bivar_no_constrain);
    }
	else if(reml_flag){
		pter_gcta->fit_reml(grm_file, phen_file, qcovar_file, covar_file, qgxe_file, gxe_file, kp_indi_file, rm_indi_file, update_sex_file, mphen, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, pred_rand_eff, est_fix_eff, reml_mtd, MaxIter, reml_priors, reml_priors_var, reml_drop, no_lrt, prevalence, no_constrain, mlma_flag, within_family, reml_bending, reml_diag_one);
	}
	else if(grm_flag || m_grm_flag){
	    if(pca_flag) pter_gcta->pca(grm_file, kp_indi_file, rm_indi_file, grm_cutoff, m_grm_flag, out_pc_num);
	    else if(make_grm_flag) pter_gcta->save_grm(grm_file, kp_indi_file, rm_indi_file, update_sex_file, grm_cutoff, grm_adj_fac, dosage_compen, m_grm_flag, grm_out_bin_flag);
	}
	else throw("Error: no analysis has been launched by the option(s).\n");
	/*      if(genome_flag) pter_gcta->simu_genome(genome_popSize, genome_nSubPOP, genome_nSubSample,
	genome_numPieces, genome_pieceLen, genome_numIndepRegion,
	genome_SNP, genome_rec, genome_mut, genome_mig);
	*/

	delete pter_gcta;
}
