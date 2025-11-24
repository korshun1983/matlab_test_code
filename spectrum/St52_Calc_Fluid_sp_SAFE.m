%===============================================================================
% Fuild: Calculation ux, uy, uz, ur, uf, T energy, U energy, T-P, P
% and angle of P, poynting vector (px,py,pz,pr,pf)
%===============================================================================
function [ResFluid] = St52_Calc_Fluid_sp_SAFE(CompStructA,FEMatricesA,PlotAniso)
    ii_d=PlotAniso.num_layer;
    cone=complex(0.,1.);

    %===============================================================================
    % Loading the results for the specific point of eigenvalue variables
    %===============================================================================
    ii_f=PlotAniso.num_freq;
    %prev_folder=cd(char(PlotAniso.dir_input));
    dat_file = strcat(char(PlotAniso.Full_Dir_Name),'\Results-',num2str(CompStructA.Model.f_array(ii_f)),'.mat');
    load(dat_file);
    %cd(prev_folder);

    rho_fluid=FEMatricesA.PhysProp{ii_d}.rho*1000.;
    lambda_fluid=FEMatricesA.PhysProp{ii_d}.lambda*1.e9;

    %===============================================================================
    % Find potentials in fluid
    %===============================================================================
    
    DNodes=FEMatricesA.DNodesComp{ii_d};
    num_nodes=length(DNodes);

    eivec_var=Results.REig_vecs(1:length(Results.REig_vecs)/2,PlotAniso.ieig);
    
    ivar_st=PlotAniso.stfn_layers(ii_d,1);
    ivar_fn=PlotAniso.stfn_layers(ii_d,2);
    eivec_fi=eivec_var(ivar_st:ivar_fn);
    
    if ii_d>1 && ...
       ( strcmpi(CompStructA.Model.DomainType(ii_d-1),'fluid') && strcmpi(CompStructA.Model.DomainType(ii_d),'fluid') )
        ivar_st_prev=PlotAniso.stfn_layers(ii_d-1,1);
        ivar_fn_prev=PlotAniso.stfn_layers(ii_d-1,2);
        eivec_fi_prev=eivec_var(ivar_st_prev:ivar_fn_prev);
    end
    
    % Calculate pressures in all nodes
    ResFluid.xx=FEMatricesA.MeshNodes(1,DNodes)';
    ResFluid.yy=FEMatricesA.MeshNodes(2,DNodes)';
    ResFluid.PP=(-cone.*Results.omega_val.*CompStructA.Misc.F_conv).*rho_fluid.*eivec_fi; % pressure
    ResFluid.fiPP=angle(ResFluid.PP); % angle of pressure
    %fid=fopen('aaa.txt','w');
    %for iii=1:length(ResFluid.fiPP)
    %    fprintf(fid,'%f\n',ResFluid.fiPP(iii));
    %end
    %fclose(fid);
%         
%     figure(50); hold on; box on;
%     plot(ResFluid.xx,ResFluid.yy,'LineStyle','o','Markersize',5,'LineWidth',2,'Color','b');
%
%     % test ~ dp/dx for given points: (x1,y1)=(0.1016,0);
%     % (x2,y2)=(0.08909,-0.001514); (x3,y3)=(0.09481,002191)
%     % find nodes
%     for inode=1:num_nodes
%         if (abs(ResFluid.xx(inode)-0.1016)<1e-4) && (abs(ResFluid.yy(inode))<1e-4)
%            INode1=inode;
%            break;
%         end
%     end
%     for inode=1:num_nodes
%         if (abs(ResFluid.xx(inode)-0.08909)<1e-4) && (abs(ResFluid.yy(inode)-(-0.001514))<1e-4)
%            INode2=inode;
%            break;
%         end
%     end
%     for inode=1:num_nodes
%         if (abs(ResFluid.xx(inode)-0.09481)<1e-4) && (abs(ResFluid.yy(inode)-(0.002191))<1e-4)
%            INode3=inode;
%            break;
%         end
%     end
% 
%     xx2=ResFluid.xx(INode3); yy2=ResFluid.yy(INode3);
%     xx1=ResFluid.xx(INode2); yy1=ResFluid.yy(INode2);
%     aa=(yy2-yy1)/(xx2-xx1);
%     bb=yy1-aa*xx1;
%     xvar=-bb/aa; yvar=0;
%     ff_var=(eivec_fi(INode3)*(xx2-xvar)+eivec_fi(INode2)*(xvar-xx1))/(xx2-xx1);
%    
%     %var_PP=-(eivec_fi(INode1)-eivec_fi(INode2))/(ResFluid.xx(INode1)-ResFluid.xx(INode2));
%     var_PP=-(eivec_fi(INode1)-ff_var)/(ResFluid.xx(INode1)-xvar);
%     var_ux=var_PP/(-cone.*Results.omega_val.*CompStructA.Misc.F_conv);
%     
%     zzz=1;
   
    %===============================================================================
    % Calculate of ux, uy, uz, ur, ufi, fi in 10 node
    %===============================================================================
    ind_tria_fluid=find(FEMatricesA.MeshTri(11,:)==ii_d); % number of triangles in fluid
    tria_fluid=FEMatricesA.MeshTri(:,ind_tria_fluid);  % triangles with nodes
    meshprops_b=FEMatricesA.MeshProps.b(:,ind_tria_fluid); % coeff b of triangles
    meshprops_c=FEMatricesA.MeshProps.c(:,ind_tria_fluid); % coeff c of triangles
    meshprops_delta=FEMatricesA.MeshProps.delta(ind_tria_fluid); % delta of triangles

    num_tria_fluid=length(ind_tria_fluid);
      
    var_ind=tria_fluid(CompStructA.Advanced.N_nodes,:);
    ResFluid.xx_10=FEMatricesA.MeshNodes(1,var_ind)';
    ResFluid.yy_10=FEMatricesA.MeshNodes(2,var_ind)';
    
    ResFluid.ux_10=zeros(num_tria_fluid,1);
    ResFluid.uy_10=zeros(num_tria_fluid,1);
    ResFluid.uz_10=zeros(num_tria_fluid,1);
    ResFluid.ur_10=zeros(num_tria_fluid,1);
    ResFluid.uf_10=zeros(num_tria_fluid,1);
    ResFluid.TE_10=zeros(num_tria_fluid,1);
    ResFluid.UE_10=zeros(num_tria_fluid,1);
    ResFluid.TmU_10=zeros(num_tria_fluid,1);
    ResFluid.PP_10=zeros(num_tria_fluid,1);
    ResFluid.fiPP_10=zeros(num_tria_fluid,1);
    ResFluid.PV_10=zeros(num_tria_fluid,5);

    % calculate dff/dx and dff/dy at 10 node (L_1=L_2=L_3=1/3) by interpolation polynome
    % ff(x,y)=sum_{i=1}^{10} ff_i*N_i, where N_i is shape function
    % N_i=1/2*(3L_i-1)*(3L_i-2)*L_i, i=1,2,3 => dN_i/dx=-(1/2)*dL_i/dx=-(1/2)*b_i/(2*Del)
    CdNdx=zeros(CompStructA.Advanced.N_nodes,1);
    CdNdx(1:3)=-0.5;
    % N_4=(9/2)*L_1*L_2*(3*L_1-1) => dN_4/dx=(3/2)*dL_1/dx=(3/2)*b_1/(2*Del)
    CdNdx(4)=1.5;
    % N_5=(9/2)*L_1*L_2*(3*L_2-1) => dN_5/dx=(3/2)*dL_2/dx=(3/2)*b_2/(2*Del)
    CdNdx(5)=1.5;
    % N_6=(9/2)*L_2*L_3*(3*L_2-1) => dN_6/dx=(3/2)*dL_2/dx=(3/2)*b_2/(2*Del)
    CdNdx(6)=1.5;
    % N_7=(9/2)*L_2*L_3*(3*L_3-1) => dN_7/dx=(3/2)*dL_3/dx=(3/2)*b_3/(2*Del)
    CdNdx(7)=1.5;
    % N_8=(9/2)*L_1*L_3*(3*L_3-1) => dN_8/dx=(3/2)*dL_3/dx=(3/2)*b_3/(2*Del)
    CdNdx(8)=1.5;
    % N_9=(9/2)*L_1*L_3*(3*L_1-1) => dN_9/dx=(3/2)*dL_1/dx=(3/2)*b_1/(2*Del)
    CdNdx(9)=1.5;
    % N_10=27*L_1*L_2*L_3 => dN_10/d(x or y)=3*(dL_1/dx+dL_2/dx+dL_3/dx)
    %                                       =3*(b_1+b_2+b_3)/(2*Del)
    CdNdx(10)=3.;
    
    % the same for y-differentiation (b_i->c_i)
    CdNdy=zeros(CompStructA.Advanced.N_nodes,1);
    CdNdy(1:3)=-0.5; CdNdy(4:9)=1.5; CdNdy(10)=3.;
    
    kk_10=Results.REig_vals(PlotAniso.ieig,PlotAniso.ieig);
    for in_10=1:num_tria_fluid
        % find potentials for all nodes (1-10) in elementary triangle
        et_nod=tria_fluid(1:CompStructA.Advanced.N_nodes,in_10);
        ff=zeros(CompStructA.Advanced.N_nodes,1);
        for inn=1:CompStructA.Advanced.N_nodes
            ii_m = (DNodes==et_nod(inn));
            if sum(ii_m)==1
               ff(inn) = eivec_fi(ii_m);
            end
            if (length(CompStructA.Model.DomainType)==ii_d) && sum(ii_m)==0
               ff(inn)=0.;
            end
            if ii_d>1 && ...
               ( strcmpi(CompStructA.Model.DomainType(ii_d-1),'fluid') && strcmpi(CompStructA.Model.DomainType(ii_d),'fluid') )
               ii_m1=(FEMatricesA.DNodesComp{ii_d-1}==et_nod(inn)); 
               if sum(ii_m1)==1
                  ff(inn) = eivec_fi_prev(ii_m1); 
               end
            end
        end
        
        % ux displacement
        bb_2d=meshprops_b(:,in_10)./(2.*meshprops_delta(in_10));
        Cbb=[bb_2d(1); bb_2d(2); bb_2d(3); ...
             bb_2d(1); bb_2d(2); bb_2d(2); bb_2d(3); bb_2d(3); bb_2d(1); ...
             (bb_2d(1)+bb_2d(2)+bb_2d(3) )];
        ResFluid.ux_10(in_10)=-sum(CdNdx.*Cbb.*ff)./(-cone.*Results.omega_val.*CompStructA.Misc.F_conv);
        
        % uy displacement
        cc_2d=meshprops_c(:,in_10)./(2.*meshprops_delta(in_10));
        Ccc=[cc_2d(1); cc_2d(2); cc_2d(3); ...
             cc_2d(1); cc_2d(2); cc_2d(2); cc_2d(3); cc_2d(3); cc_2d(1); ...
             (cc_2d(1)+cc_2d(2)+cc_2d(3) )];
        ResFluid.uy_10(in_10)=-sum(CdNdy.*Ccc.*ff)./(-cone.*Results.omega_val.*CompStructA.Misc.F_conv);
        
        % uz displacement
        ResFluid.uz_10(in_10)=(kk_10./(Results.omega_val.*CompStructA.Misc.F_conv)).*ff(CompStructA.Advanced.N_nodes);
        
        % ur & uf displacement
        fi_10=atan(ResFluid.yy_10(in_10)/abs(ResFluid.xx_10(in_10)));
        if (ResFluid.xx_10(in_10)<0. && ResFluid.yy_10(in_10)>0.), fi_10= pi-fi_10; end
        if (ResFluid.xx_10(in_10)<0. && ResFluid.yy_10(in_10)<0.), fi_10=-pi-fi_10; end
        
        %ResFluid.ur_10(in_10)=sqrt(abs(ResFluid.ux(in_10)).^2+abs(ResFluid.uy(in_10)).^2);
        ResFluid.ur_10(in_10)= cos(fi_10)*ResFluid.ux_10(in_10)+sin(fi_10)*ResFluid.uy_10(in_10);
        
        % uf displacement
        ResFluid.uf_10(in_10)=-sin(fi_10)*ResFluid.ux_10(in_10)+cos(fi_10)*ResFluid.uy_10(in_10);
        
        % T energy
        varT=abs(ResFluid.ux_10(in_10)).^2+abs(ResFluid.uy_10(in_10)).^2+abs(ResFluid.uz_10(in_10)).^2;
        ResFluid.TE_10(in_10)=0.5*rho_fluid*varT.*(Results.omega_val.*CompStructA.Misc.F_conv).^2;
        
        % U energy
        varU=(rho_fluid./lambda_fluid).*abs(ff(CompStructA.Advanced.N_nodes)).^2;
        ResFluid.UE_10(in_10)=0.5*rho_fluid.*varU.*(Results.omega_val.*CompStructA.Misc.F_conv).^2;
        
        % T-U energy
        ResFluid.TmU_10(in_10)=ResFluid.TE_10(in_10)-ResFluid.UE_10(in_10);
        
        % Pressure at 10 node
        ResFluid.PP_10(in_10)=-cone.*Results.omega_val.*rho_fluid.*ff(CompStructA.Advanced.N_nodes); % pressure
        ResFluid.fiPP_10(in_10)=angle(ResFluid.PP_10(in_10)); % angle of pressure
        
        %                                   | sxx sxy sxz |                 | p 0 0 |
        % Poynting Vector, PV=-0.5*conj(v)* | sxy syy syz | =  0.5*conj(v)* | 0 p 0 |
        %                                   | sxz syz szz |                 | 0 0 p |
        
        var_vv=[ResFluid.ux_10(in_10); ResFluid.uy_10(in_10); ResFluid.uz_10(in_10)];
        var_vv=var_vv.*(-cone*Results.omega_val.*CompStructA.Misc.F_conv);
        conj_vv=conj(var_vv).';
        ResFluid.PV_10(in_10,1:3)=(conj_vv.*ResFluid.PP_10(in_10));
        
        % Poyinting Vector in Polar Coordinate (4 and 5 indexes)
        fi_10=atan(ResFluid.yy_10(in_10)/abs(ResFluid.xx_10(in_10)));
        if (ResFluid.xx_10(in_10)<0. && ResFluid.yy_10(in_10)>0.), fi_10= pi-fi_10; end
        if (ResFluid.xx_10(in_10)<0. && ResFluid.yy_10(in_10)<0.), fi_10=-pi-fi_10; end
        ResFluid.PV_10(in_10,4)= cos(fi_10)*ResFluid.PV_10(in_10,1)+sin(fi_10)*ResFluid.PV_10(in_10,2);
        ResFluid.PV_10(in_10,5)=-sin(fi_10)*ResFluid.PV_10(in_10,1)+cos(fi_10)*ResFluid.PV_10(in_10,2);
        
    end
    
    zzz=1;
end % program
    
    
    
    
    
    
