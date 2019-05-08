
function [u_opt, w_opt]=D_opt(pars,Omega,u_vec,FIM_vec,bounds)

    N_u=length(u_vec);
    Np=length(pars);
    
    cvx_begin quiet
      %variable A(Np,Np) symmetric
      variable lambda(N_u)
      FIM_tot=zeros(Np,Np);
      for i=1:N_u
          FIM_tot=FIM_tot+FIM_vec(:,:,i)*lambda(i);
      end
      maximize (det_rootn(FIM_tot) )%det_rootn( FIM_tot )*det_inv(FIM_tot(N_th+1:end,N_th+1:end)))
      subject to
        sum(lambda) == 1;
        lambda >= 0;
    cvx_end
    lam = lambda; % save the solution for confidence ellipsoids

    lam(lam<1e-3)=0;
    boolLam=lam>0;
    nxtBool=logical(boolLam(1:end-1).*boolLam(2:end));

    flg=sum(nxtBool)~=0;
    
    cuttCnt=7;
    cntr=1;
    while flg&&cntr<cuttCnt
        new_u=[(u_vec([false; nxtBool])+u_vec([nxtBool; false]))/2];
        [u_vec,ind]=sort([u_vec new_u]);
        N_u=length(u_vec);
        FIM_tmp=cat(3,FIM_vec ,FIM_comp(pars,Omega,new_u,bounds));
        FIM_vec=FIM_tmp(:,:,ind);
        
        cvx_begin quiet
          %variable A(Np,Np) symmetric
          variable lambda(N_u)
          FIM_tot=zeros(Np,Np);
          for i=1:N_u
              FIM_tot=FIM_tot+FIM_vec(:,:,i)*lambda(i);
          end
          maximize (det_rootn(FIM_tot) )%det_rootn( FIM_tot )*det_inv(FIM_tot(N_th+1:end,N_th+1:end)))
          subject to
            sum(lambda) == 1;
            lambda >= 0;
        cvx_end
        lam=lambda;
        
        lam(lam<1e-3)=0;
        boolLam=lam>0;
        nxtBool=logical(boolLam(1:end-1).*boolLam(2:end));
        
        flg=sum(nxtBool)~=0;
        cntr=cntr+1;
    end
    
    if ~(cntr<cuttCnt)
        display('Warning iteration count reached!')
    end
    
    w_opt=lam(boolLam);
    u_opt=u_vec(boolLam);


