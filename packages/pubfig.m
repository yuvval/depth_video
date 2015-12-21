
classdef pubfig
% A class of static functions that automates figures styling for scientific publications
% pubfig.funcname(...)

   properties(Constant)
       pathname = '~/www/figs';
   end
   
   properties(Constant)
       
   end
   
   methods(Static)
       function text_size(fontsize, fhnd)
           % function text_size(fontsize, fhnd)
           % set the text font size of all the axes in a
           % given figure

           
           if (nargin < 2)
               fhnd = gcf;
           end
           
           axAllH = findobj(fhnd, 'type', 'axes').'; 
           
           
           for h = axAllH
               set(h, 'FontSize', fontsize);
               set(get(h, 'Title'), 'FontSize', fontsize);
               set(get(h, 'Xlabel'), 'FontSize', fontsize);
               set(get(h, 'Ylabel'), 'FontSize', fontsize);
               set(get(h, 'Zlabel'), 'FontSize', fontsize);               
           end     
       end

       function mul_axis_tick_by_factor(factor, XYZind, ahnd)
           % function mul_axis_tick_by_factor(factor, XYZind, ahnd)
           % multiple axis ticks by a given factor
           % XYZind is 'X' or 'Y' or 'Z'
           % ahnd is the axes handle, default value is gca
           
           if (nargin < 3)
               ahnd = gca;
           end
           axticks = get(ahnd, [XYZind 'Tick']);           
           set(ahnd, [XYZind 'TickLabel'], factor*axticks);         
       end
       
       function make_ax_exp(XYZind,  ahnd, basis)
           % function make_ax_exp(XYZind,  ahnd, basis)
           % make an axis exponent. 
           % E.g. axis changes from -2, -3, .. to '10^-2', '10^-3', ...
           % 
           % XYZind is 'X' or 'Y' or 'Z'
           % ahnd is the axes handle, default value is gca
           % basis is 2, 10, 'exp'. default value is 10
           if (nargin < 2)
               ahnd = gca;
           end           
           if (nargin < 3)
               basis = 10;
           end           
           
           switch basis
               case 10
                   strbasis = '10';
               case 2
                   strbasis = '2';
               case 'exp'
                   strbasis = 'e';
           end
                   
          
           
           pubfig.set_axishnd_tick_spacing(5, ahnd, XYZind, true);
           axticks = get(ahnd, [XYZind 'TickLabel']);
           newticks = [];
           for cnt = 1:size(axticks,1)
               newticks{cnt} = [strbasis '^' axticks(cnt,:)];
           end
           set(ahnd, [XYZind 'TickLabel'], newticks);
           
       end
       
       
       function mirror_axis( XYZind, ahnd)
           % function mirror_axis( XYZind, ahnd)
           % mirror axis ticks. 
           % E.g. 1, 2, 3, 4, 5 changes to 5, 4, 3, 2, 1
           % XYZind is 'X' or 'Y' or 'Z'
           % ahnd is the axes handle, default value is gca
           if (nargin < 2)
               ahnd = gca;
           end
           axticks = get(ahnd, [XYZind 'TickLabel']);           
           set(ahnd, [XYZind 'TickLabel'], flipud(axticks));
       end
       
       function set_axishnd_tick_spacing(nintervals, ahnd, XYZind, toround)
           % function set_axishnd_tick_spacing(nintervals, ahnd, XYZind, toround)
           % set the spacing between a given axis values
           % XYZind is 'X' or 'Y' or 'Z'
           % ahnd is the axes handle
           % to round: true/false. true means use ticks with round values.
           %           default value is false
           
           if (nargin < 4)
               toround  = false;
           end
           inlims  = get(ahnd, [XYZind 'Lim']);
           
           interval = (diff(inlims)/nintervals);
           if(toround == true)
               interval = round(interval);
           end
           if toround == false
               set(ahnd, [XYZind 'Tick'], inlims(1):interval:inlims(2));
           else 
               set(ahnd, [XYZind 'Tick'], floor(inlims(1)):interval:ceil(inlims(2)));
           end
           
       end

       function set_axis_tick_spacing(nintervals, fhnd, toround)
           % function set_axis_tick_spacing(nintervals, fhnd, toround)
           % set the spacing between values of all figure axes
           % fhnd is the figure handle, default is gcf
           % to round: true/false. true means use ticks with round values.
           %           default value is false
           if (nargin < 2)
               fhnd = gcf;
           end
           
           if (nargin < 3)
               toround  = false;
           end
           
           axAllH = findobj(fhnd, 'type', 'axes').'; 
           axLegH = findobj(fhnd, 'tag', 'legend').'; 
           axAllH = setdiff(axAllH, axLegH); % remove the legends axes
           
           
           for h = axAllH
               
               cnt = 1;
               for ax = {'X', 'Y', 'Z'}
                   
                   inlims  = get(h, [ax{1} 'Lim']);
                   
                   interval = (diff(inlims)/nintervals);
                   if(toround == true)
                       interval = round(interval);
                   end
                   
                   set(h, [ax{1} 'Tick'], inlims(1):interval:inlims(2));
%                    set(h, [ax{1} 'TickLabel'], xyzFactors(cnt)*(inlims(1):interval:inlims(2)));
                   
                   cnt = cnt+1;
               end
           end    
           
       end

       
       function setbox(onoff, fhnd)
           % function setbox(onoff, fhnd)
           % set figure box to on/off
           % onoff values are 'on' or 'off'
           % fhnd is the figure handle, default is gcf
           if (nargin < 2)
               fhnd = gcf;
           end
           
           axAllH = findobj(fhnd, 'type', 'axes').';
           axLegH = findobj(fhnd, 'tag', 'legend').'; 
           axAllH = setdiff(axAllH, axLegH); % remove the legends axes
           
           
           for h = axAllH               
               set(h, 'Box', onoff)
           end
           
       end
       
       function print_eps(filename,  fhnd, other_format)
           % function print_eps(filename,  fhnd, other_format)
           % prints a figure to eps file
           % fhnd is the figure handle, default is gcf
           %
           % other_format, allows to export the eps as other format as well
           % e.g. 'jpeg', 'tiff', 'png'. Default value is []. 
           % To use other_format, you should install eps2xxx() from 
           % http://www.mathworks.com/matlabcentral/fileexchange/6858-eps2xxx
           % and also you need to install Ghostscript
           
           if nargin<2
               fhnd = gcf;
           end
           if nargin < 3
               other_format = [];
           end
           A3sizecm = [ 14.85, 21];
           pos = get(fhnd, 'Position');
           pos = pos(3:4);
           [mxpos, xy] = max(pos);
           paperpos = [A3sizecm(xy) A3sizecm(xy)*pos(3-xy)/mxpos];
           if (xy == 2)
               paperpos = fliplr(paperpos);
           end
           
           set(fhnd, 'PaperUnits', 'centimeters');
           set(fhnd, 'PaperPosition', [0 0 paperpos]);
           
           print(fhnd, '-depsc2', '-r300', [filename '.eps']);
           if ~isempty(other_format)
               eps2xxx([filename '.eps'], other_format);
           end
       end
       
       function print_fig(filename,  ftype, fhnd)
           % function print_fig(filename,  ftype, fhnd)
           % prints a figure to a specific format file
           % ftype is the format type. e.g. 'jpeg', 'tiff', 'png'. 
           %       Default value is 'tiff'
           % fhnd is the figure handle, default is gcf
           if nargin<3
               fhnd = gcf;
           end
           
           if nargin<2
               ftype = 'tiff';
           end    
           
           switch ftype
               case 'tiff'
                   ext = '.tif';
               case 'jpeg'
                   ext = '.jpg';
               otherwise
                   ext = ['.' ftype];
           end
           
           if nargin < 1
               
               filename = ['figure' num2str(fhnd)];
           end
           
           
           if ~ispc
               pubfig.print_eps(filename,  fhnd, {ftype})
           else
               
               A3sizecm = [ 14.85, 21];
               pos = get(fhnd, 'Position');
               pos = pos(3:4);
               [mxpos, xy] = max(pos);
               paperpos = [A3sizecm(xy) A3sizecm(xy)*pos(3-xy)/mxpos];
               if (xy == 2)
                   paperpos = fliplr(paperpos);
               end
               
               set(fhnd, 'PaperUnits', 'centimeters');
               set(fhnd, 'PaperPosition', [0 0 paperpos]);
               
               print(fhnd, ['-d' ftype], '-r300', [filename ext]);
           end
           
       end
       
       function save(fname_fig,  ftype, fhnd)
           % function print_fig(filename,  ftype, fhnd)
           % prints a figure to a specific format file
           % ftype is the format type. e.g. 'jpeg', 'tiff', 'png'. 
           %       Default value is 'tiff'
           % fhnd is the figure handle, default is gcf
           if nargin<3
               fhnd = gcf;
           end
           
           if nargin<2
               ftype = 'tiff';
           end    
           
           switch ftype
               case 'tiff'
                   ext = 'tif';
               case 'jpeg'
                   ext = 'jpg';
               otherwise
                   ext = [ftype];
           end
           
           if nargin < 1
               
               fname_fig = ['figure' num2str(fhnd)];
           end
           
           saveas(fhnd, fname_fig, ext);
       end
       function template1(fhnd)
           % function template1(fhnd)
           % template1 is a automatic template to run a set of styling
           % command on a figure
           % fhnd is the figure handle, default is gcf
           %
           
           
           if (nargin < 1)
               fhnd = gcf;
           end
           
           pubfig.text_size(20, fhnd);
           pubfig.set_axis_tick_spacing(5, fhnd);
           pubfig.setbox('off', fhnd)
       end

       function your_own_template(fhnd)
           % function your_own_template(fhnd)
           % fhnd is the figure handle, default is gcf
           % you are encouraged to add your own styling template
           
           if (nargin < 1)
               fhnd = gcf;
           end
       end
           
       function fname_fig = titlestr_to_fname(title_str)
           fname_fig = regexprep(title_str, ' ', '_' );
           fname_fig = regexprep(fname_fig, ',', '' );
           fname_fig = regexprep(fname_fig, '\\_', '_' );
           fname_fig = regexprep(fname_fig, '\n', '_' );
       end
           
       function char = newline()
           % function char = newline()
           % return a newline char (sprintf('\n'))
           char = sprintf('\n');
       end
       
       
           
       
       
       
   end
end