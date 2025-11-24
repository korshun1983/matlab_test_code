%===============================================================================
% Processing Fluid Results of Kinetic Energy
% for specific frequency number and all eigs
%===============================================================================
function [ResFluid] = St61_Proc_Fluid_TE_sp_SAFE(CompStructA,FEMatricesA,ProcTE)
    ii_l=ProcTE.ilayer;
    cone=complex(0.,1.);
   
    %===============================================================================
    % Build Polar Grid
    %===============================================================================
    % At first, we choose max_r of layer in r direction and find appropriate
    % triangles. Then among these triangles we find the vertex of
    % a triangle which is close to the origin and choose appropriate node at 10 node for RR grid.
    % Using this vertex as next start point we have cycle until RR(end)
    % reaches min r of this layer.
    
    rr=[]; 
    rr_min_prev_var=FEMatricesA.DomainRx(ii_l);
    
    % find solid triangles and their nodes
    tri_fluid_ind=(FEMatricesA.MeshTri(CompStructA.Advanced.N_nodes+1,:)==ii_l);
    tri_fluid=FEMatricesA.MeshTri(1:CompStructA.Advanced.N_nodes,tri_fluid_ind);
    if ii_l==1 % the smallest r of mesh
        rr_min_mesh=min(sqrt(FEMatricesA.MeshNodes(1,tri_fluid(1:3,:)).^2+FEMatricesA.MeshNodes(2,tri_fluid(1:3,:)).^2));
    else
        rr_min_mesh=FEMatricesA.DomainRx(ii_l-1);
    end
    
    % find all triangles lying on external circle
    rr_circ=[];
    for iel=1:size(tri_fluid,2)
        xx=FEMatricesA.MeshNodes(1,tri_fluid(:,iel));
        yy=FEMatricesA.MeshNodes(2,tri_fluid(:,iel));
        RR2=(xx./FEMatricesA.DomainRx(ii_l)).^2+(yy./FEMatricesA.DomainRy(ii_l)).^2;
        [RR21,~]=sort(RR2,'descend');
        if RR21(1)>0.99 && RR21(2)>0.99
           RR_10=sqrt(xx(10).^2+yy(10).^2);
           rr_circ(1:end+1)=[rr_circ(1:end) RR_10]; 
        end        
    end

    % set initial node with r=r_max and z=0
    r_find=FEMatricesA.DomainRx(ii_l); z_find=0.;
    
    key_count='y'; r_end_corr=1.; %r_end_corr=0.98; 
        
    i_count=0; 
    while (1)
        i_count=i_count+1;
        %if i_count>2, break; end
        node_find=find(FEMatricesA.MeshNodes(1,:)==r_find); %& FEMatricesA.MeshNodes(2,:)==z_find);
        ivar=size(node_find,2);
        if ivar>1,  node_find(ivar:end)=[]; end
        
        num_tri_node_find=[]; % find number of triangles with this node_find
        for iel=1:size(tri_fluid,2)
            if sum((tri_fluid(1:CompStructA.Advanced.N_nodes,iel)==node_find))==1
               num_tri_node_find(1:end+1)=[num_tri_node_find(1:end) iel];
            end
        end
        
        tri_nodes=tri_fluid(1:CompStructA.Advanced.N_nodes,num_tri_node_find);
        
        rr_min_var=rr_min_prev_var; tri_node_min=1;
        for iet=1:size(tri_nodes,2)
            tri_node_xx=FEMatricesA.MeshNodes(1,tri_nodes(:,iet))';
            tri_node_yy=FEMatricesA.MeshNodes(2,tri_nodes(:,iet))';
            rr_var3=sqrt(tri_node_xx(1:3).^2+tri_node_yy(1:3).^2); % choose only 1,2,3 nodes 
            for ia=1:3
                if rr_var3(ia)<rr_min_var; rr_min_var=rr_var3(ia); tri_node_min=iet; end
            end
        end
                
        tri_node_xx=FEMatricesA.MeshNodes(1,tri_nodes(:,tri_node_min))';
        tri_node_yy=FEMatricesA.MeshNodes(2,tri_nodes(:,tri_node_min))';
        rr_var=sqrt(tri_node_xx.^2+tri_node_yy.^2);
        if tri_nodes(1,tri_node_min)==node_find % 1 node
           if rr_var(2)>rr_var(3)
              ind_rr_1=9; ind_rr_2=8; ind_rr=3;
           else
              ind_rr_1=4; ind_rr_2=5; ind_rr=2;
           end
        end
        if tri_nodes(2,tri_node_min)==node_find  % 2 node
           if rr_var(1)>rr_var(3)
              ind_rr_1=6; ind_rr_2=7; ind_rr=3;
           else
              ind_rr_1=5; ind_rr_2=4; ind_rr=1;
           end
        end
        if tri_nodes(3,tri_node_min)==node_find % 3 node
           if rr_var(1)>rr_var(2)
              ind_rr_1=7; ind_rr_2=6; ind_rr=2;
           else
              ind_rr_1=8; ind_rr_2=9; ind_rr=1;
           end
        end
        ind_rr_all=[ind_rr_1; ind_rr_2; ind_rr];

        for iet=1:size(tri_nodes,2) % choose 10 node
            rr_add=sqrt( FEMatricesA.MeshNodes(1,tri_nodes(CompStructA.Advanced.N_nodes,iet)).^2 ...
                        +FEMatricesA.MeshNodes(2,tri_nodes(CompStructA.Advanced.N_nodes,iet)).^2);              
            rr(1:end+1)=[rr(1:end) rr_add];
        end
        rr_min_prev_var=rr_var(ind_rr);

        if rr_min_prev_var<=(rr_min_mesh*1.01), break; end
        
        r_find=tri_node_xx(ind_rr); z_find=tri_node_yy(ind_rr);
        tri_fluid(:,num_tri_node_find)=[];
    
    end
    rr_circ_min=min(rr_circ)*0.995;
    ind = ( rr>rr_circ_min );
    rr(ind)=[];
    rr=unique(sort(rr))';
    rr=[rr; rr_circ_min];

    %ir_delete=[];
    %for ir=2:length(rr)
    %    if (rr(ir)-rr(ir-1))<=1.e-4, ir_delete=[ir_delete; ir]; end
    %end
    %rr(ir_delete)=[];
    
    % tuning of rr(1) and rr(end)
    Nrr=length(rr);
    rr_ind=zeros(Nrr,1);
    for irr=Nrr:-1:2
        var=abs(rr(irr)-rr(irr-1))/rr(irr);
        if var<0.03 % 3%
           rr_ind(irr)=1; 
        end
    end
    rr(logical(rr_ind))=[];
    rr(1)=rr(1)*1.02;
    rr(end)=rr(end)*r_end_corr;
    Nrr=length(rr);
    ResFluid.rr=rr;
    
    %===============================================================================
    rho_fluid=FEMatricesA.PhysProp{ii_l}.rho*1000.;
    DThetaV = pi/(CompStructA.Model.DomainNth(ii_l)*ProcTE.DThetaInc); 
    ThetaV = (( - pi):DThetaV: pi )'; NThetaV=length(ThetaV);
    % end of building of polar grid
    %===============================================================================
    
    %===============================================================================
    % make memory for answer  
    ResFluid.TE_fer=zeros(length(ProcTE.neigs),Nrr); % T Energy
    % last index: 1 - for plus, 2 - for minus 
    ResFluid.cF_ur_rm_pm=zeros(length(ProcTE.neigs),Nrr,ProcTE.NHarm,2); 
    ResFluid.cF_uf_rm_pm=zeros(length(ProcTE.neigs),Nrr,ProcTE.NHarm,2); 
    ResFluid.cF_uz_rm_pm=zeros(length(ProcTE.neigs),Nrr,ProcTE.NHarm,2); 
    
    % Common parameters for all freqs
    % Rectangular grid
    xx=FEMatricesA.MeshNodes(1,FEMatricesA.DNodesComp{ii_l})';
    yy=FEMatricesA.MeshNodes(2,FEMatricesA.DNodesComp{ii_l})';

    ind_tria_fluid=find(FEMatricesA.MeshTri(11,:)==ii_l); % number of triangles in fluid
    tria_fluid=FEMatricesA.MeshTri(:,ind_tria_fluid);  % triangles with nodes
    meshprops_b=FEMatricesA.MeshProps.b(:,ind_tria_fluid); % coeff b of triangles
    meshprops_c=FEMatricesA.MeshProps.c(:,ind_tria_fluid); % coeff c of triangles
    meshprops_delta=FEMatricesA.MeshProps.delta(ind_tria_fluid); % delta of triangles

    num_tria_fluid=length(ind_tria_fluid);

    var_ind=tria_fluid(CompStructA.Advanced.N_nodes,:);
    xx_10=FEMatricesA.MeshNodes(1,var_ind)';
    yy_10=FEMatricesA.MeshNodes(2,var_ind)';
    
    % Coefficients to calculate dff/dx and dff/dy at 10 node (L_1=L_2=L_3=1/3) by interpolation polynome
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
    
    %===============================================================================
    % work massive for Discrete Fourier transform
    exp_var_p=zeros(NThetaV,ProcTE.NHarm);
    for iHarm=1:ProcTE.NHarm, exp_var_p(:,iHarm)=exp( cone.*(iHarm-1).*ThetaV(:)); end
    exp_var_m=conj(exp_var_p);

    %===============================================================================
    %===============================================================================
    % Main Caqlculations for given eigs
    %===============================================================================
    %===============================================================================
    ifreq=ProcTE.ii_f;
    ProcTE.freq=ProcTE.freq_count(ifreq);

    %===============================================================================
    % Loading the results for the specific frequency number
    %===============================================================================
    curr_dir=pwd; cd('..');
    dat_file = char(strcat(ProcTE.Dir_Name,'\Results-',num2str(ProcTE.freq_count(ifreq)),'.mat'));
    load(dat_file);
    cd(curr_dir);

    iii_eig=0; key_count_again='y'; key_count_below='y';
    while (key_count_again=='y')
    for ieig=ProcTE.neigs
        iii_eig=iii_eig+1;
        
        % Read Data of Pressure at all nodes
        eivec_var=Results.REig_vecs(1:length(Results.REig_vecs)/2,ieig);

        ivar_st=ProcTE.stfn_layers(ii_l,1);
        ivar_fn=ProcTE.stfn_layers(ii_l,2);
        eivec_fi=eivec_var(ivar_st:ivar_fn);

        if ii_l>1 && ...
           ( strcmp(CompStructA.Model.DomainType(ii_l-1),'fluid') && strcmp(CompStructA.Model.DomainType(ii_l),'fluid') )
            ivar_st_prev=ProcTE.stfn_layers(ii_l-1,1);
            ivar_fn_prev=ProcTE.stfn_layers(ii_l-1,2);
            eivec_fi_prev=eivec_var(ivar_st_prev:ivar_fn_prev);
        end

        %===============================================================================
        % Calculate of ux, uy, uz, ur, uf, te at 10 nodes
        %===============================================================================
        
        ux_10=zeros(num_tria_fluid,1);
        uy_10=zeros(num_tria_fluid,1);
        uz_10=zeros(num_tria_fluid,1);
        ur_10=zeros(num_tria_fluid,1);
        uf_10=zeros(num_tria_fluid,1);
        ff_10=zeros(num_tria_fluid,1);

        kk_10=Results.REig_vals(ieig,ieig);
        for in_10=1:num_tria_fluid

            % find potentials for all nodes (1-10) in elementary triangle
            et_nod=tria_fluid(1:CompStructA.Advanced.N_nodes,in_10);
            ff=zeros(CompStructA.Advanced.N_nodes,1);
            for inn=1:CompStructA.Advanced.N_nodes
                ii_m = (FEMatricesA.DNodesComp{ii_l}==et_nod(inn));
               
                if sum(ii_m)==1
                   ff(inn) = eivec_fi(ii_m);
                end
               
                if (length(CompStructA.Model.DomainType)==ii_l) && sum(ii_m)==0
                    ff(inn)=0.;
                end
                if ii_l>1 && ...
                   ( strcmpi(CompStructA.Model.DomainType(ii_l-1),'fluid') && strcmpi(CompStructA.Model.DomainType(ii_l),'fluid') )
                   ii_m1=(FEMatricesA.DNodesComp{ii_l-1}==et_nod(inn)); 
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
            ux_10(in_10)=-sum(CdNdx.*Cbb.*ff)./(-cone.*Results.omega_val.*CompStructA.Misc.F_conv);

            % uy displacement
            cc_2d=meshprops_c(:,in_10)./(2.*meshprops_delta(in_10));     
            Ccc=[cc_2d(1); cc_2d(2); cc_2d(3); ...
                 cc_2d(1); cc_2d(2); cc_2d(2); cc_2d(3); cc_2d(3); cc_2d(1); ...
                 (cc_2d(1)+cc_2d(2)+cc_2d(3) )];
            uy_10(in_10)=-sum(CdNdy.*Ccc.*ff)./(-cone.*Results.omega_val.*CompStructA.Misc.F_conv);

            % uz displacement
            uz_10(in_10)=(kk_10./(Results.omega_val.*CompStructA.Misc.F_conv)).*ff(CompStructA.Advanced.N_nodes);

            % ur & uf displacement
            fi_10=atan(yy_10(in_10)/abs(xx_10(in_10)));
            if (xx_10(in_10)<0. && yy_10(in_10)>0.), fi_10= pi-fi_10; end
            if (xx_10(in_10)<0. && yy_10(in_10)<0.), fi_10=-pi-fi_10; end

            %ur_10(in_10)=sqrt(abs(ux(in_10)).^2+abs(uy(in_10)).^2);
            ur_10(in_10)= cos(fi_10)*ux_10(in_10)+sin(fi_10)*uy_10(in_10);

            % uf displacement
            uf_10(in_10)=-sin(fi_10)*ux_10(in_10)+cos(fi_10)*uy_10(in_10);

            % potential f 
            ff_10(in_10)=ff(10);
        end

        %===============================================================================
        % ReCalculate of ur_fr, uf_rf, uz_rf in Polar Grid
        %===============================================================================
        xx_rf=zeros(Nrr,NThetaV); yy_rf=zeros(Nrr,NThetaV);
        ur_rf_10=zeros(Nrr,NThetaV); uf_rf_10=zeros(Nrr,NThetaV);  uz_rf_10=zeros(Nrr,NThetaV);  ff_rf_10=zeros(Nrr,NThetaV);

        %Fr_re_var=TriScatteredInterp(xx_10,yy_10,real(ur_10)); % define interpolation function in rectanglar grid
        Fr_re_var=scatteredInterpolant(xx_10,yy_10,real(ur_10)); % define interpolation function in rectanglar grid
        Fr_im_var=scatteredInterpolant(xx_10,yy_10,imag(ur_10)); % define interpolation function in rectanglar grid

        Ff_re_var=scatteredInterpolant(xx_10,yy_10,real(uf_10)); % define interpolation function in rectanglar grid
        Ff_im_var=scatteredInterpolant(xx_10,yy_10,imag(uf_10)); % define interpolation function in rectanglar grid

        Fz_re_var=scatteredInterpolant(xx_10,yy_10,real(uz_10)); % define interpolation function in rectanglar grid
        Fz_im_var=scatteredInterpolant(xx_10,yy_10,imag(uz_10)); % define interpolation function in rectanglar grid

        FF_re_var=scatteredInterpolant(xx_10,yy_10,real(ff_10)); % define interpolation function in rectanglar grid
        FF_im_var=scatteredInterpolant(xx_10,yy_10,imag(ff_10)); % define interpolation function in rectanglar grid
        
        for ir=Nrr:-1:1
            rr_ir=ones(length(ThetaV),1).*rr(ir)';
            [qxx,qyy] = pol2cart(ThetaV,rr_ir);
            xx_rf(ir,:)=qxx; yy_rf(ir,:)=qyy;

            qzz_re=Fr_re_var(qxx,qyy); qzz_im=Fr_im_var(qxx,qyy); 
            
            if ir==Nrr
               num_nan=isnan(qzz_re); % check for nan answers
               if sum(num_nan)>=1 % we have nan answer hence we'll correct r_max
                   % ivar=[1:length(rr_ir)];
                   % ind_nan=ivar(num_nan);
                   % [xx_nan,yy_nan] = pol2cart(ThetaV(ind_nan),rr_ir(1));
                   % r2_10_mas=[];
                   % for inan=1:length(ind_nan)
                   %     [~,ind]=min(((xx_nan(inan)-xx_10).^2+(yy_nan(inan)-yy_10).^2));
                   %     r2_nan=xx_nan(inan)^2+yy_nan(inan)^2;
                   %     r2_10=xx_10(ind)^2+yy_10(ind)^2;
                   %     if r2_nan>r2_10, r2_10_mas(1:end+1)=[r2_10_mas(1:end) r2_10]; end
                   % end
                   % rr(ir)=min(sqrt(r2_10_mas))*0.999;
                   % ResFluid.rr(end)=rr(ir);
                   % fprintf(1,'\tr_end_corr is corrected!\n');
                   % key_count_below='n';
                   rr(ir)=rr(ir)*0.995;
                   fprintf(1,'\tr_end_corr is corrected!\n');
                   key_count_below='n';
                   break;
               end
            end
            
            ur_rf_10(ir,:)=complex(qzz_re,qzz_im);
            qzz_re=Ff_re_var(qxx,qyy); qzz_im=Ff_im_var(qxx,qyy); 
            uf_rf_10(ir,:)=complex(qzz_re,qzz_im);
            qzz_re=Fz_re_var(qxx,qyy); qzz_im=Fz_im_var(qxx,qyy); 
            uz_rf_10(ir,:)=complex(qzz_re,qzz_im);
            qzz_re=FF_re_var(qxx,qyy); qzz_im=FF_im_var(qxx,qyy); 
            ff_rf_10(ir,:)=complex(qzz_re,qzz_im);

        end

        if key_count_below=='y'
        
            %%% Test: Compare ur in Cartesian and Polar coordinates
            %figure(ProcTE.num_fig+501); 
            %set(gca,'FontSize',ProcTE.font_size); hold on; box on; view(45,45);
            %xlim([min(xx_10) max(xx_10)]); ylim([min(yy_10) max(yy_10)]);
            %tri=delaunay(xx_10,yy_10);
            %trimesh(tri,xx_10,yy_10,abs(ur_10));
            %xlabel('x (m)'); ylabel('y (m)');
            %title_show='Fluid, Surface of T Energy';
            %title(title_show);
            % 
    %         if iii_eig==1 
    %         figure(ProcTE.num_fig+502); 
    %         set(gca,'FontSize',ProcTE.font_size); hold on; box on;
    %         for ir=1:Nrr
    %            view(45,45); plot3(xx_rf(ir,:),yy_rf(ir,:),abs(ur_rf_10(ir,:)),'r'); 
    %         end
    %         return;
    %         end
            %===============================================================================
            % Calculate T Energy 
            %===============================================================================
            varC=0.5*rho_fluid*(2.*pi*ProcTE.freq.*CompStructA.Misc.F_conv).^2;
            var_df=DThetaV/2.;
            for ir=1:Nrr
                var_te =  abs(ur_rf_10(ir,(1:NThetaV-1))).^2 ...
                        + abs(uf_rf_10(ir,(1:NThetaV-1))).^2 ...
                        + abs(uz_rf_10(ir,(1:NThetaV-1))).^2;

                ResFluid.TE_fer(iii_eig,ir)=varC.*sum(var_te).*var_df;
            end

            %===============================================================================
            % Fourier decomposition of ur_rf, uf_rf, uz_rf by angle fi for fixed r
            %===============================================================================

            %===============================================================================
            % make memory for temporary answer  
            var_ur_pm=zeros(Nrr,ProcTE.NHarm,2); 
            var_uf_pm=zeros(Nrr,ProcTE.NHarm,2); 
            var_uz_pm=zeros(Nrr,ProcTE.NHarm,2); 

            % I_i ~ 0.5*(f(x_{i-1})+f(x_i))*(x_i-x_{i-1})
            % Int f(x)=((f_0+f_n)/2+sum_{i=1}^{n-1)f_i)*dx, dx=(b-a)/n
            % we have that f_0=f_n !
            for iHarm=1:ProcTE.NHarm
            for ir=1:Nrr
                var_ur_pm(ir,iHarm,1)=sum((ur_rf_10(ir,(1:NThetaV-1)).').*exp_var_p((1:NThetaV-1),iHarm)).*var_df;
                var_uf_pm(ir,iHarm,1)=sum((uf_rf_10(ir,(1:NThetaV-1)).').*exp_var_p((1:NThetaV-1),iHarm)).*var_df;
                var_uz_pm(ir,iHarm,1)=sum((uz_rf_10(ir,(1:NThetaV-1)).').*exp_var_p((1:NThetaV-1),iHarm)).*var_df;

                var_ur_pm(ir,iHarm,2)=sum((ur_rf_10(ir,(1:NThetaV-1)).').*exp_var_m((1:NThetaV-1),iHarm)).*var_df;
                var_uf_pm(ir,iHarm,2)=sum((uf_rf_10(ir,(1:NThetaV-1)).').*exp_var_m((1:NThetaV-1),iHarm)).*var_df;
                var_uz_pm(ir,iHarm,2)=sum((uz_rf_10(ir,(1:NThetaV-1)).').*exp_var_m((1:NThetaV-1),iHarm)).*var_df;
            end
            end

            var_pm=var_ur_pm;
            ResFluid.cF_ur_rm_pm(iii_eig,:,:,:)=var_pm./(2.*pi);

            var_pm=var_uf_pm;
            ResFluid.cF_uf_rm_pm(iii_eig,:,:,:)=var_pm./(2.*pi);

            var_pm=var_uz_pm;
            ResFluid.cF_uz_rm_pm(iii_eig,:,:,:)=var_pm./(2.*pi);
            
            if ieig==length(ProcTE.neigs), key_count_again='n'; end
        else
            break;
        end % if key_count_below

        zzz=1;
    end % cycle of eigs
    end % while key_count_again
    
    %if sum(sum(sum(isnan(ResFluid.TE_fer))))==0
    %   key_count='n';
    %else
    %   r_end_corr=r_end_corr*0.999;
    %   fprintf(1,'\tr_end_corr is corrected!\n');
    %end

    % Plot polar grid for testing
    if lower(ProcTE.Polar_Grid)=='y'
        figure(ProcTE.num_fig+(ProcTE.ii_f-1)*10+ProcTE.ilayer); 
        set(gca,'FontSize',ProcTE.font_size); 
        hold on; box on;
        plot(xx,yy,'LineStyle','none','Marker','o','Markersize',5,'Color','b');
        plot(xx_10,yy_10,'LineStyle','none','Marker','x','Markersize',10,'Color','r');
        for ir=1:Nrr
            [qxx,qyy] = pol2cart(ThetaV,rr(ir));
            plot(qxx,qyy,'Color','k','LineWidth',2);
        end
        xlabel('x (m)'); ylabel('y (m)'); title('Polar Grid in Fluid');
        
        tri_fluid_ind=(FEMatricesA.MeshTri(CompStructA.Advanced.N_nodes+1,:)==ii_l);
        tri_fluid=FEMatricesA.MeshTri(1:CompStructA.Advanced.N_nodes,tri_fluid_ind);
        num_tri=size(tri_fluid,2);
        for iel=1:num_tri
            pang=tri_fluid(1:3,iel);
            pang_new=[pang; pang(1)];
            xx=FEMatricesA.MeshNodes(1,pang_new);
            yy=FEMatricesA.MeshNodes(2,pang_new);
            plot(xx,yy,'Color','b');
        end
        
    end
    
    zzz=1;
%  

  
