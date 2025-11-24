%===============================================================================
% Solid: Plot ux, uy, uz, ur, uf, T energy, U energy, T-U, Spur, angle of Spur
% poynting vector (px,py,pz,pr,pf)
%===============================================================================

function St52_Plot_Solid_sp_SAFE(CompStructA,FEMatricesA,PlotAniso,ResSolid)

    text1_add=num2str(PlotAniso.freq);
    text2_add=num2str(PlotAniso.ieig);
    ResSolid.title_add=strcat(text1_add,' kHz, k','_{',text2_add,'}');
    Result.title_add=ResSolid.title_add;
    Result.title_add_file=Result.title_add;
    
    %===============================================================================
    %===============================================================================
    % 3d Surface
    if lower(PlotAniso.plot_3d)=='y'
       key_plot_c_save=PlotAniso.plot_c;
       if lower(PlotAniso.plot_c)=='y', PlotAniso.plot_c='n'; end
   
       if lower(PlotAniso.subplot)=='y'
          text21_add=num2str(PlotAniso.neigs(1));
          text22_add=num2str(PlotAniso.neigs(length(PlotAniso.neigs)));
          ResSolid.title_add_file=strcat('3D,F=',text1_add,',Eigs=',text21_add,'-',text22_add);
          Result.title_add_file=ResSolid.title_add_file;
       end
        
       %============================================================================
       % plot at all nodes 
       Result.xx=ResSolid.xx; Result.yy=ResSolid.yy; 
       Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
       Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);
       
       Result.field='ux';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3000;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |u_x|:';
              Result.function=abs(ResSolid.ux);
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |u_x|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='uy';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3100;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |u_y|:';
              Result.function=abs(ResSolid.uy);
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |u_y|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='uz';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3200;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |u_z|:';
              Result.function=abs(ResSolid.uz);
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |u_z|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='ur';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3300;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |u_r|:';
              Result.function=abs(ResSolid.ur);
              %Result.function=imag(ResSolid.ur);              
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |u_r|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='uf';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3400;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |u_f|:';
              Result.function=abs(ResSolid.uf);
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |u_f|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='te';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3500;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of T Energy:';
              Result.function=ResSolid.TE;
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized T Energy:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       %============================================================================
       % plot at 10 nodes
  
       Result.xx=ResSolid.xx_10; Result.yy=ResSolid.yy_10; 
       Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
       Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);

       Result.field='ue';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3600;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of U Energy:';
              Result.function=ResSolid.UE_10;
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized U Energy:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='tu';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3700;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of T-U Energy:';
              Result.function=ResSolid.TmU_10;
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized T-U Energy:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='px';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4000;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |Poyn_x|:';
              Result.function=abs(ResSolid.PV_10(:,1));
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |Poyn_x|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='py';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4100;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |Poyn_y|:';
              Result.function=abs(ResSolid.PV_10(:,2));
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |Poyn_y|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pz';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4200;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |Poyn_z|:';
              Result.function=abs(ResSolid.PV_10(:,3));
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |Poyn_z|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pr';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4300;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |Poyn_r|:';
              Result.function=abs(ResSolid.PV_10(:,4));
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |Poyn_r|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pf';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4400;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |Poyn_f|:';
              Result.function=abs(ResSolid.PV_10(:,5));
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |Poyn_f|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pp';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3800;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of |Spur|:';
              Result.function=abs(ResSolid.Spur_10);
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized |Spur|:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='ap';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+3900;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end

              Result.title_show='Solid, Surface of Angle of Spur:';
              Result.function=angle(ResSolid.Spur_10);
              if lower(PlotAniso.normalization)=='y'
                  Result.title_show='Solid, Surface of Normalized Angle of Spur:';
                  Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       if lower(key_plot_c_save)=='y', PlotAniso.plot_c='y'; end
    
    end % plot 3D

    %===============================================================================
    %===============================================================================
    if lower(PlotAniso.plot_c)=='y'
       key_plot_3d_save=PlotAniso.plot_3d;
       if lower(PlotAniso.plot_3d)=='y', PlotAniso.plot_3d='n'; end
      
       if lower(PlotAniso.subplot)=='y'
          text21_add=num2str(PlotAniso.neigs(1));
          text22_add=num2str(PlotAniso.neigs(length(PlotAniso.neigs)));
          ResSolid.title_add_file=strcat('C,F=',text1_add,',Eigs=',text21_add,'-',text22_add);
          Result.title_add_file=ResSolid.title_add_file;
       end

       %============================================================================ 
       % plot at all nodes
       Result.xx=ResSolid.xx; Result.yy=ResSolid.yy; 
       Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
       Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);
       dx=(Result.xx_max-Result.xx_min)/PlotAniso.num_xy(1); dy=(Result.yy_max-Result.yy_min)/PlotAniso.num_xy(2);
       [Result.qxx,Result.qyy]=meshgrid(Result.xx_min:dx:Result.xx_max,Result.yy_min:dy:Result.yy_max);
       

       Result.field='ux';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4500;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |u_x|:';
              Result.function=abs(ResSolid.ux);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |u_x|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='uy';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4600;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |u_y|:';
              Result.function=abs(ResSolid.uy);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |u_y|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
      
       Result.field='uz';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4700;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |u_z|:';
              Result.function=abs(ResSolid.uz);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |u_z|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='ur';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4800;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |u_r|:';
              Result.function=abs(ResSolid.ur);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |u_r|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='uf';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+4900;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |u_f|:';
              Result.function=abs(ResSolid.uf);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |u_f|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end

       Result.field='te';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5000;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of T Energy:';
              Result.function=ResSolid.TE;
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized T Energy:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       %============================================================================
       % plot at 10 nodes
       Result.xx=ResSolid.xx_10; Result.yy=ResSolid.yy_10; 
       Result.xx_min=min(Result.xx); Result.xx_max=max(Result.xx);
       Result.yy_min=min(Result.yy); Result.yy_max=max(Result.yy);
       dx=(Result.xx_max-Result.xx_min)/PlotAniso.num_xy(1); dy=(Result.yy_max-Result.yy_min)/PlotAniso.num_xy(2);
       [Result.qxx,Result.qyy]=meshgrid(Result.xx_min:dx:Result.xx_max,Result.yy_min:dy:Result.yy_max);

       Result.field='ue';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5100;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of U Energy:';
              Result.function=ResSolid.UE_10;
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized U Energy:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='tu';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5200;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of T-U Energy:';
              Result.function=ResSolid.TmU_10;
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized T-U Energy:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='px';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5500;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |Poyn_x|:';
              Result.function=abs(ResSolid.PV_10(:,1));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |Poyn_x|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='py';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5600;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |Poyn_y|:';
              Result.function=abs(ResSolid.PV_10(:,2));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |Poyn_y|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pz';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5700;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |Poyn_z|:';
              Result.function=abs(ResSolid.PV_10(:,3));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |Poyn_z|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pr';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5800;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |Poyn_r|:';
              Result.function=abs(ResSolid.PV_10(:,4));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |Poyn_r|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='pf';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5900;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |Poyn_f|:';
              Result.function=abs(ResSolid.PV_10(:,5));
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |Poyn_f|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
      
       
       Result.field='pp';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5300;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of |Spur|:';
              Result.function=abs(ResSolid.Spur_10);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized |Spur|:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       Result.field='ap';
       for ifield=1:size(PlotAniso.fields,1)
           if strcmp(Result.field,PlotAniso.fields{ifield,1}) && strcmp('y',PlotAniso.fields{ifield,2})
              Result.num_fig=PlotAniso.num_fig+5400;
              if lower(PlotAniso.subplot)=='n'
                 Result.num_fig=Result.num_fig+(PlotAniso.iii_eig-1);
              end
             
              Result.title_show='Solid, Contours of Angle of Spur:';
              Result.function=angle(ResSolid.Spur_10);
              if lower(PlotAniso.normalization)=='y'
                 Result.title_show='Solid, Contours of Normalized Angle of Spur:';
                 Result.function=Result.function./max(abs(Result.function));
              end
              
              St53_Plot_1Field_sp_SAFE(PlotAniso,Result);
           end
       end
       
       if lower(key_plot_3d_save)=='y', PlotAniso.plot_3d='y'; end

    end % plot contour
    
end % program
