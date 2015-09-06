if strcmp(seq,'LiftFrontNi100')
    kk=245:304;
    if strcmp(alg,'Solari_Chessa')
        kk=245:301;
    end
    else if strcmp(seq,'PointingDiagon')
        kk=77:84;
        if strcmp(alg,'Solari_Chessa')
            kk=77:81;
        end
    else if strcmp(seq,'PointingFront')
        kk=50:59;
        if strcmp(alg,'Solari_Chessa')
            kk=50:56;
        end
    else if strcmp(seq,'TranspDiagon')
    	kk=35:49; 
        if strcmp(alg,'Solari_Chessa')
            kk=35:46; 
        end
    else if strcmp(seq,'TranspDiagonRobot')
        kk=30:41; 
        if strcmp(alg,'Solari_Chessa')
            kk=30:38;
        end
    else if strcmp(seq,'TranspFront')
        kk=30:44;
        if strcmp(alg,'Solari_Chessa')
            kk=30:41;
        end
    else if strcmp(seq,'TranspFrontRobot')
        kk=166:179; 
        if strcmp(alg,'Solari_Chessa')
            kk=166:176;
        end
	else if strcmp(seq,'GestGabFront')
        kk=150:189; 
        if strcmp(alg,'Solari_Chessa')
            kk=150:186;

        
        end
        end
        end
        end
        end
        end
        end
        end
end
    
for k=kk
    cd results
    cd(name)
    kname=sprintf('%8.8d', k);
    workspace_name=strcat(name,'_',kname);
    load (workspace_name)
    cd .\..\..
    cd flow-code-matlab
    computeColor(uv(:,:,1), uv(:,:,2));          %codifica i campi flow u,v in colori
                      % e quindi dà come risultato ans (matrice MxNx3)
    figure
    h=imagesc(ans);		  % per visualizzare la matrice (e quindi il flow)
    cd .\..
    cd results
    cd (name)
    mkdir jpg
    cd jpg
    
    saveas(h, workspace_name, 'jpg')  
    cd .\..\..\..
    close all
end