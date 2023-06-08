function val= lossfunction_growth_full(x,C,t_vector,a_vector,W)
    rho_u = x(1);
    rho_v = x(2);
    rho_w = x(3);
    sigma2_v = x(4);
    sigma2_w = x(5);
    sigma2_e = x(6);
    sigma2_z = x(7);
    psi = x(8);
 
    
    % theoretical moments
    cov_mat = covariance_growth_full(t_vector,a_vector,rho_u,rho_v,rho_w,sigma2_v,sigma2_w,sigma2_e,sigma2_z,psi);
        
    % distance matrix
    e = (cov_mat - C);
    
    %e(isnan(moments(:))) = [];
        val = e'*W*e;
% 
%     if exist('W','var')==1
%         val = e'*W*e;
%     else
%         val = e'*e;
%     end


%     e_phi_w = e_phi.*W_phi;
%     e_var_w = e_var.*W_var;
%     e_w = [e_var_w(:); e_phi_w(:)];
%     e_w(isnan(e_w)) = [];

%    val = sum(e.^2);
end