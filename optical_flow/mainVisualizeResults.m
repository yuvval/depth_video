for i=47:48
    switch i
        case 0
            seq='LiftFrontNi100';
            alg='CLG-TV';
           
        case 1
            seq='LiftFrontNi100';
            alg='Lucas_Kanade';

        case 2
            seq='LiftFrontNi100';
            alg='ba';
            
        case 3
            seq='LiftFrontNi100';
            alg='hs';
          
        case 4
            seq='LiftFrontNi100';
            alg='Solari_Chessa';
            
        case 5
            seq='PointingDiagon';
            alg='CLG-TV';
            
        case 6
            seq='PointingDiagon';
            alg='Lucas_Kanade';
            
        case 7
            seq='PointingDiagon';
            alg='ba';
                       
        case 8
            seq='PointingDiagon';
            alg='hs';
       
        case 9
            seq='PointingDiagon';
            alg='Solari_Chessa';
        case 10
            seq='PointingFront';
            alg='CLG-TV';
            name=strcat(seq,'_',alg);
            for k=50:59
            	acquisitionSequence
                main2
            end
            
        case 11
            seq='PointingFront';
            alg='Lucas_Kanade';

        case 12
            seq='PointingFront';
            alg='ba';
            
        case 13
            seq='PointingFront';
            alg='hs';

        case 14
            seq='PointingFront';
            alg='Solari_Chessa';

        case 15
            seq='TranspDiagon';
            alg='CLG-TV';
            
        case 16
            seq='TranspDiagon';
            alg='Lucas_Kanade';
            
        case 17
            seq='TranspDiagon';
            alg='ba';
            
        case 18
            seq='TranspDiagon';
            alg='hs';

        case 19
            seq='TranspDiagon';
            alg='Solari_Chessa';
            
        case 20
            seq='TranspDiagonRobot';
            alg='CLG-TV';

        case 21
            seq='TranspDiagonRobot';
            alg='Lucas_Kanade';

        case 22
            seq='TranspDiagonRobot';
            alg='ba';

        case 23
            seq='TranspDiagonRobot';
            alg='hs';

        case 24
            seq='TranspDiagonRobot';
            alg='Solari_Chessa';

        case 25
            seq='TranspFront';
            alg='CLG-TV';
            name=strcat(seq,'_',alg);

        case 26
            seq='TranspFront';
            alg='Lucas_Kanade';

        case 27
            seq='TranspFront';
            alg='ba';

        case 28
            seq='TranspFront';
            alg='hs';

        case 29
            seq='TranspFront';
            alg='Solari_Chessa';

        case 30
            seq='TranspFrontRobot';
            alg='CLG-TV';
            
        case 31
            seq='TranspFrontRobot';
            alg='Lucas_Kanade';
            
        case 32
            seq='TranspFrontRobot';
            alg='ba';
           
        case 33
            seq='TranspFrontRobot';
            alg='hs';
            
        case 34
            seq='TranspFrontRobot';
            alg='Solari_Chessa';            
            
        case 35
            seq='LiftFrontNi100';
            alg='classic++';

        case 36
            seq='PointingDiagon';
            alg='classic++';
            
        case 37
            seq='PointingFront';
            alg='classic++';

        case 38
            seq='TranspDiagon';
            alg='classic++';

        case 39
            seq='TranspDiagonRobot';
            alg='classic++';
            
        case 40
            seq='TranspFront';
            alg='classic++';
            
        case 41
            seq='TranspFrontRobot';
            alg='classic++';
        
            
            
        case 42
            seq='GestGabFront';
            alg='CLG-TV';


        case 43
            seq='GestGabFront';
            alg='Lucas_Kanade';

        case 44
            seq='GestGabFront';
            alg='ba';

        case 45
            seq='GestGabFront';
            alg='hs';

        case 46
            seq='GestGabFront';
            alg='Solari_Chessa';
            
        case 47
            seq='GestGabFront';
            alg='classic++';
            
            
    end
	name=strcat(seq,'_',alg);
    visualizeResult
    clear
    close all
    
end