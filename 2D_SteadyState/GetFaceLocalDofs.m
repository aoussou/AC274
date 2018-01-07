


function [v V_e] = GetFaceLocalDofs(face_id)

global p

switch face_id
    case 1
        if p == 0
            V_e = 1;
        elseif p ==1
            V_e = [1 2];
        elseif p==2
            V_e = [1 4 2];
        end
        v = [1 2];
    case  2
        v = [2,3];
        if p == 0
            V_e = 1;
        elseif p ==1
            V_e = [2 3];
        elseif p==2
            V_e = [2,5,3];
        end
    case 3
        v = [3,1];
        if p == 0
            V_e = 1;
        elseif p ==1
            V_e = [3 1];
        elseif p==2
            V_e = [3 6 1];
        end
end