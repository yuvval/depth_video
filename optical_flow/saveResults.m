%save results
    cd results
    %tmp_name=strcat(name);
    mkdir(name);
    cd (name)     

    kname=sprintf('%8.8d', k);
    workspace_name=strcat(name,'_',kname);
    save(workspace_name,'uv')
    
    cd ./../..
    cd flow-code-matlab
    computeColor(uv(:,:,1), uv(:,:,2));
    figure
    h=imagesc(ans);		  % per visualizzare la matrice (e quindi il flow)
    
    cd ./..
    cd results
    cd (name)
    
    mkdir jpg
    cd jpg
    saveas(h,workspace_name, 'jpg')  
    cd ./../../..			%siamo nella cartella comparison_of
    %results_visualization
    %close all
