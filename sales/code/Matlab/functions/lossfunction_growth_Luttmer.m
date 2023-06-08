function val = lossfunction_growth_Luttmer(x,C,t_vector,a_vector,W)
%lossfunction_pss_autocov_level

    in_alpha = x(1);
    in_delta = x(2);
    rho = x(3);
    sigma2_q = x(4);
    sigma2_e = x(5);
    lambda = x(6);
    
    % theoretical moments
    cov_mat = covariance_growth_Luttmer(t_vector,a_vector,in_alpha,in_delta,rho,sigma2_q,sigma2_e,lambda);
        
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