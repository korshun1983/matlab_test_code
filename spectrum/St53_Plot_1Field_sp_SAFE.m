%===============================================================================
% Fuild: Plot only one specified field from St52_Plot_Fluid_sp_SAFE
%===============================================================================

function St53_Plot_1Field_sp_SAFE(PlotAniso,Result)

      PlotAniso.hf=figure(Result.num_fig);
   
      if lower(PlotAniso.subplot)=='y'
         var_iii_eig=PlotAniso.iii_eig;
      else
         var_iii_eig=1;
      end
      
      subplot(PlotAniso.subplot_nxm(1),PlotAniso.subplot_nxm(2),var_iii_eig);
      set(gca,'FontSize',PlotAniso.font_size); hold on; box on; view(PlotAniso.view_3D);
      xlim([Result.xx_min Result.xx_max]); ylim([Result.yy_min Result.yy_max]);
            
      if lower(PlotAniso.plot_3d)=='y'
          %trisurf(tri,xx,yy,uxn);
          tri=delaunay(Result.xx,Result.yy);
          trimesh(tri,Result.xx,Result.yy,Result.function);
          pause(0.2); colormap('jet');
      end
      if lower(PlotAniso.plot_c)=='y'
          F_var=TriScatteredInterp(Result.xx,Result.yy,Result.function);
          qzz=F_var(Result.qxx,Result.qyy);
          contourf(Result.qxx,Result.qyy,qzz,PlotAniso.num_cont);
          pause(0.2); colormap('jet');
      end
      
      if lower(PlotAniso.subplot)=='n'
         title(strcat(Result.title_show,Result.title_add),'VerticalAlignment','bottom');
      else
         %title(strcat(Result.field,',',Result.title_add));
         title(Result.title_add,'VerticalAlignment','bottom');
      end

      if lower(PlotAniso.xlabel_show)=='y'
         xlabel('x (m)'); ylabel('y (m)');
      end
      
end % program
