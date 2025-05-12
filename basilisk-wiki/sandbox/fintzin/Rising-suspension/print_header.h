void print_header_drop(FILE * f){
    fprintf(f,"i,t,tagmax,tracer");
    fprintf(f,",tag");
    fprintf(f,",realtag");
    fprintf(f,",tagmin");
    fprintf(f,",age");
    fprintf(f,",class");
    fprintf(f,",vol");
    fprintf(f,",S");
    fprintf(f,",Sc");
    fprintf(f,",ed");
    fprintf(f,",Mnorm");
    fprintf(f,",Pnorm");
    fprintf(f,",p");
    fprintf(f,",W2");
    #if dimension ==2
    fprintf(f,",Y_x");
    fprintf(f,",Y_y");
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",D_nbr_x");
    fprintf(f,",D_nbr_y");
    fprintf(f,",mv_x");
    fprintf(f,",mv_y");
    fprintf(f,",Fh_x");
    fprintf(f,",Fh_y");
    fprintf(f,",Fv_x");
    fprintf(f,",Fv_y");
    fprintf(f,",Fp_x");
    fprintf(f,",Fp_y");
    fprintf(f,",nc_x");
    fprintf(f,",nc_y");
    // I J
    fprintf(f,",M_x_x");
    fprintf(f,",M_x_y");
    fprintf(f,",M_y_x");
    fprintf(f,",M_y_y");

    fprintf(f,",P_x_x");
    fprintf(f,",P_x_y");
    fprintf(f,",P_y_x");
    fprintf(f,",P_y_y");
    
    fprintf(f,",WW_x_x");
    fprintf(f,",WW_x_y");
    fprintf(f,",WW_y_x");
    fprintf(f,",WW_y_y");
    
    fprintf(f,",T_x_x");
    fprintf(f,",T_x_y");
    fprintf(f,",T_y_x");
    fprintf(f,",T_y_y");

    fprintf(f,",rdT_x_x");
    fprintf(f,",rdT_x_y");
    fprintf(f,",rdT_y_x");
    fprintf(f,",rdT_y_y");

    fprintf(f,",Ms_x_x");
    fprintf(f,",Ms_x_y");
    fprintf(f,",Ms_y_x");
    fprintf(f,",Ms_y_y");
    
    #elif dimension ==3
    // X
    fprintf(f,",Y_x");
    fprintf(f,",Y_y");
    fprintf(f,",Y_z");
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",U_z");
    fprintf(f,",D_nbr_x");
    fprintf(f,",D_nbr_y");
    fprintf(f,",D_nbr_z");
    fprintf(f,",mv_x");
    fprintf(f,",mv_y");
    fprintf(f,",mv_z");
    fprintf(f,",Fh_x");
    fprintf(f,",Fh_y");
    fprintf(f,",Fh_z");
    fprintf(f,",Fv_x");
    fprintf(f,",Fv_y");
    fprintf(f,",Fv_z");
    fprintf(f,",Fp_x");
    fprintf(f,",Fp_y");
    fprintf(f,",Fp_z");
    fprintf(f,",nc_x");
    fprintf(f,",nc_y");
    fprintf(f,",nc_z");

    // I J
    fprintf(f,",M_x_x");
    fprintf(f,",M_x_y");
    fprintf(f,",M_x_z");
    fprintf(f,",M_y_x");
    fprintf(f,",M_y_y");
    fprintf(f,",M_y_z");
    fprintf(f,",M_z_x");
    fprintf(f,",M_z_y");
    fprintf(f,",M_z_z");

    fprintf(f,",P_x_x");
    fprintf(f,",P_x_y");
    fprintf(f,",P_x_z");
    fprintf(f,",P_y_x");
    fprintf(f,",P_y_y");
    fprintf(f,",P_y_z");
    fprintf(f,",P_z_x");
    fprintf(f,",P_z_y");
    fprintf(f,",P_z_z");

    fprintf(f,",WW_x_x");
    fprintf(f,",WW_x_y");
    fprintf(f,",WW_x_z");
    fprintf(f,",WW_y_x");
    fprintf(f,",WW_y_y");
    fprintf(f,",WW_y_z");
    fprintf(f,",WW_z_x");
    fprintf(f,",WW_z_y");
    fprintf(f,",WW_z_z");

    fprintf(f,",T_x_x");
    fprintf(f,",T_x_y");
    fprintf(f,",T_x_z");
    fprintf(f,",T_y_x");
    fprintf(f,",T_y_y");
    fprintf(f,",T_y_z");
    fprintf(f,",T_z_x");
    fprintf(f,",T_z_y");
    fprintf(f,",T_z_z");

    fprintf(f,",rdT_x_x");
    fprintf(f,",rdT_x_y");
    fprintf(f,",rdT_x_z");
    fprintf(f,",rdT_y_x");
    fprintf(f,",rdT_y_y");
    fprintf(f,",rdT_y_z");
    fprintf(f,",rdT_z_x");
    fprintf(f,",rdT_z_y");
    fprintf(f,",rdT_z_z");

    fprintf(f,",Ms_x_x");
    fprintf(f,",Ms_x_y");
    fprintf(f,",Ms_x_z");
    fprintf(f,",Ms_y_x");
    fprintf(f,",Ms_y_y");
    fprintf(f,",Ms_y_z");
    fprintf(f,",Ms_z_x");
    fprintf(f,",Ms_z_y");
    fprintf(f,",Ms_z_z");
    #endif
    fprintf(f,"\n");
}

void print_header_CA(FILE * f){
    fprintf(f,"i,t,NVOF,PHI,S");
    fprintf(f,",Vd,Vf,dissd,dissf,pd,pf");
    #if dimension == 2
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",Us_x");
    fprintf(f,",Us_y");
    fprintf(f,",Ud_x");
    fprintf(f,",Ud_y");
    fprintf(f,",Uf_x");
    fprintf(f,",Uf_y");
    fprintf(f,",RhoU_x");
    fprintf(f,",RhoU_y");
    fprintf(f,",A_x");
    fprintf(f,",A_y");

    fprintf(f,",UpUpf_x_x");
    fprintf(f,",UpUpf_x_y");
    fprintf(f,",UpUpf_y_x");
    fprintf(f,",UpUpf_y_y");
    
    fprintf(f,",UpUpd_x_x");
    fprintf(f,",UpUpd_x_y");
    fprintf(f,",UpUpd_y_x");
    fprintf(f,",UpUpd_y_y");

    fprintf(f,",UpUpf2_x_x");
    fprintf(f,",UpUpf2_x_y");
    fprintf(f,",UpUpf2_y_x");
    fprintf(f,",UpUpf2_y_y");
    
    fprintf(f,",UpUpd2_x_x");
    fprintf(f,",UpUpd2_x_y");
    fprintf(f,",UpUpd2_y_x");
    fprintf(f,",UpUpd2_y_y");

    fprintf(f,",UpUpUpf_x");
    fprintf(f,",UpUpUpf_x");
    
    fprintf(f,",UpUpUpd_x");
    fprintf(f,",UpUpUpd_y");
    #elif dimension == 3
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",U_z");
    fprintf(f,",Us_x");
    fprintf(f,",Us_y");
    fprintf(f,",Us_z");
    fprintf(f,",Ud_x");
    fprintf(f,",Ud_y");
    fprintf(f,",Ud_z");
    fprintf(f,",Uf_x");
    fprintf(f,",Uf_y");
    fprintf(f,",Uf_z");
    fprintf(f,",RhoU_x");
    fprintf(f,",RhoU_y");
    fprintf(f,",RhoU_z");
    fprintf(f,",A_x");
    fprintf(f,",A_y");
    fprintf(f,",A_z");

    fprintf(f,",UpUpf_x_x");
    fprintf(f,",UpUpf_x_y");
    fprintf(f,",UpUpf_x_z");
    fprintf(f,",UpUpf_y_x");
    fprintf(f,",UpUpf_y_y");
    fprintf(f,",UpUpf_y_z");
    fprintf(f,",UpUpf_z_x");
    fprintf(f,",UpUpf_z_y");
    fprintf(f,",UpUpf_z_z");

    fprintf(f,",UpUpd_x_x");
    fprintf(f,",UpUpd_x_y");
    fprintf(f,",UpUpd_x_z");
    fprintf(f,",UpUpd_y_x");
    fprintf(f,",UpUpd_y_y");
    fprintf(f,",UpUpd_y_z");
    fprintf(f,",UpUpd_z_x");
    fprintf(f,",UpUpd_z_y");
    fprintf(f,",UpUpd_z_z");

    fprintf(f,",UpUpf2_x_x");
    fprintf(f,",UpUpf2_x_y");
    fprintf(f,",UpUpf2_x_z");
    fprintf(f,",UpUpf2_y_x");
    fprintf(f,",UpUpf2_y_y");
    fprintf(f,",UpUpf2_y_z");
    fprintf(f,",UpUpf2_z_x");
    fprintf(f,",UpUpf2_z_y");
    fprintf(f,",UpUpf2_z_z");

    fprintf(f,",UpUpd2_x_x");
    fprintf(f,",UpUpd2_x_y");
    fprintf(f,",UpUpd2_x_z");
    fprintf(f,",UpUpd2_y_x");
    fprintf(f,",UpUpd2_y_y");
    fprintf(f,",UpUpd2_y_z");
    fprintf(f,",UpUpd2_z_x");
    fprintf(f,",UpUpd2_z_y");
    fprintf(f,",UpUpd2_z_z");

    fprintf(f,",Tf_x");
    fprintf(f,",Tf_y");
    fprintf(f,",Tf_z");

    fprintf(f,",Td_x");
    fprintf(f,",Td_y");
    fprintf(f,",Td_z");
    #endif
    fprintf(f,"\n");
}

void print_header_CAnst(FILE * f){
    fprintf(f,"i,t");
    #if dimension == 2
    fprintf(f,",D_x");
    fprintf(f,",D_y");
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",Ur_x");
    fprintf(f,",Ur_y");
    fprintf(f,",a_x");
    fprintf(f,",a_y");
    fprintf(f,",f");
    fprintf(f,",p");

    #elif dimension == 3
    fprintf(f,",D_x");
    fprintf(f,",D_y");
    fprintf(f,",D_z");
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",U_z");
    fprintf(f,",Ur_x");
    fprintf(f,",Ur_y");
    fprintf(f,",Ur_z");
    fprintf(f,",a_x");
    fprintf(f,",a_y");
    fprintf(f,",a_z");
    fprintf(f,",f");
    fprintf(f,",p");

    #endif
    fprintf(f,"\n");
}

void print_header_PA(FILE * f){
    fprintf(f,"i,t,Nb,V,S,Sc,p,ed");
    fprintf(f,",Mnorm");
    fprintf(f,",Pnorm");
    fprintf(f,",W2");
    #if dimension == 2
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",mv_x");
    fprintf(f,",mv_y");

    fprintf(f,",M_x_x");
    fprintf(f,",M_x_y");
    fprintf(f,",M_y_x");
    fprintf(f,",M_y_y");

    fprintf(f,",P_x_x");
    fprintf(f,",P_x_y");
    fprintf(f,",P_y_x");
    fprintf(f,",P_y_y");

    fprintf(f,",UpUp_x_x");
    fprintf(f,",UpUp_x_y");
    fprintf(f,",UpUp_y_x");
    fprintf(f,",UpUp_y_y");
    
    fprintf(f,",T_x_x");
    fprintf(f,",T_x_y");
    fprintf(f,",T_y_x");
    fprintf(f,",T_y_y");
    
    fprintf(f,",rdT_x_x");
    fprintf(f,",rdT_x_y");
    fprintf(f,",rdT_y_x");
    fprintf(f,",rdT_y_y");
    
    fprintf(f,",Ms_x_x");
    fprintf(f,",Ms_x_y");
    fprintf(f,",Ms_y_x");
    fprintf(f,",Ms_y_y");

    fprintf(f,",WW_x_x");
    fprintf(f,",WW_x_y");
    fprintf(f,",WW_y_x");
    fprintf(f,",WW_y_y");

    fprintf(f,",R_x");
    fprintf(f,",R_x");

    fprintf(f,",PFP_x_x");
    fprintf(f,",PFP_x_y");
    fprintf(f,",PFP_y_x");
    fprintf(f,",PFP_y_y");

    fprintf(f,",UR_x_x");
    fprintf(f,",UR_x_y");
    fprintf(f,",UR_y_x");
    fprintf(f,",UR_y_y");
    
    fprintf(f,",RR_x_x");
    fprintf(f,",RR_x_y");
    fprintf(f,",RR_y_x");
    fprintf(f,",RR_y_y");
    
    #elif dimension == 3
    fprintf(f,",U_x");
    fprintf(f,",U_y");
    fprintf(f,",U_z");
    fprintf(f,",mv_x");
    fprintf(f,",mv_y");
    fprintf(f,",mv_z");

    fprintf(f,",M_x_x");
    fprintf(f,",M_x_y");
    fprintf(f,",M_x_z");
    fprintf(f,",M_y_x");
    fprintf(f,",M_y_y");
    fprintf(f,",M_y_z");
    fprintf(f,",M_z_x");
    fprintf(f,",M_z_y");
    fprintf(f,",M_z_z");

    fprintf(f,",P_x_x");
    fprintf(f,",P_x_y");
    fprintf(f,",P_x_z");
    fprintf(f,",P_y_x");
    fprintf(f,",P_y_y");
    fprintf(f,",P_y_z");
    fprintf(f,",P_z_x");
    fprintf(f,",P_z_y");
    fprintf(f,",P_z_z");

    fprintf(f,",UpUp_x_x");
    fprintf(f,",UpUp_x_y");
    fprintf(f,",UpUp_x_z");
    fprintf(f,",UpUp_y_x");
    fprintf(f,",UpUp_y_y");
    fprintf(f,",UpUp_y_z");
    fprintf(f,",UpUp_z_x");
    fprintf(f,",UpUp_z_y");
    fprintf(f,",UpUp_z_z");
        
    fprintf(f,",T_x_x");
    fprintf(f,",T_x_y");
    fprintf(f,",T_x_z");
    fprintf(f,",T_y_x");
    fprintf(f,",T_y_y");
    fprintf(f,",T_y_z");
    fprintf(f,",T_z_x");
    fprintf(f,",T_z_y");
    fprintf(f,",T_z_z");
    
    fprintf(f,",rdT_x_x");
    fprintf(f,",rdT_x_y");
    fprintf(f,",rdT_x_z");
    fprintf(f,",rdT_y_x");
    fprintf(f,",rdT_y_y");
    fprintf(f,",rdT_y_z");
    fprintf(f,",rdT_z_x");
    fprintf(f,",rdT_z_y");
    fprintf(f,",rdT_z_z");

    fprintf(f,",Ms_x_x");
    fprintf(f,",Ms_x_y");
    fprintf(f,",Ms_x_z");
    fprintf(f,",Ms_y_x");
    fprintf(f,",Ms_y_y");
    fprintf(f,",Ms_y_z");
    fprintf(f,",Ms_z_x");
    fprintf(f,",Ms_z_y");
    fprintf(f,",Ms_z_z");

    fprintf(f,",WW_x_x");
    fprintf(f,",WW_x_y");
    fprintf(f,",WW_x_z");
    fprintf(f,",WW_y_x");
    fprintf(f,",WW_y_y");
    fprintf(f,",WW_y_z");
    fprintf(f,",WW_z_x");
    fprintf(f,",WW_z_y");
    fprintf(f,",WW_z_z");

    fprintf(f,",R_x");
    fprintf(f,",R_y");
    fprintf(f,",R_z");

    fprintf(f,",PFP_x_x");
    fprintf(f,",PFP_x_y");
    fprintf(f,",PFP_x_z");
    fprintf(f,",PFP_y_x");
    fprintf(f,",PFP_y_y");
    fprintf(f,",PFP_y_z");
    fprintf(f,",PFP_z_x");
    fprintf(f,",PFP_z_y");
    fprintf(f,",PFP_z_z");

    fprintf(f,",UR_x_x");
    fprintf(f,",UR_x_y");
    fprintf(f,",UR_x_z");
    fprintf(f,",UR_y_x");
    fprintf(f,",UR_y_y");
    fprintf(f,",UR_y_z");
    fprintf(f,",UR_z_x");
    fprintf(f,",UR_z_y");
    fprintf(f,",UR_z_z");

    fprintf(f,",RR_x_x");
    fprintf(f,",RR_x_y");
    fprintf(f,",RR_x_z");
    fprintf(f,",RR_y_x");
    fprintf(f,",RR_y_y");
    fprintf(f,",RR_y_z");
    fprintf(f,",RR_z_x");
    fprintf(f,",RR_z_y");
    fprintf(f,",RR_z_z");
    #endif
    fprintf(f,"\n");
}

void print_header_layer_files(FILE * f){
    int nL = round(Ls/(D/LDn));
    fprintf(f,"t");
    for (int i = 0; i < nL; i++)
      fprintf(f,",L%d",i);
    fprintf(f,"\n");   
}
