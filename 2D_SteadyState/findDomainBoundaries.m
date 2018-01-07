% Output is the element number and the side number that 

function  [BExt] = findDomainBoundaries

Globals2D;
global initial_condition_string
K = size(EToN,1); %number of triangular elements

%% find something that will find the correct side of the element

% the two physical vertices
% the interpolated values at each of the three
% glboal degrees of freedom

% find elements with an edge on the edge of the domain = partial omega

[kb sb] = find(ElemNeighbors == -1);%the element on the boundary kb and the side of that element sb
% figure out if the boundary edge

numElemEdges = length(sb);
l = 1;
r = 1;
t = 1;
b = 1;

side2Vert = [1 2; 2 3; 3 1];
for i = 1:length(kb)
    w1 = VX(EToN(kb(i),side2Vert(sb(i),1)),:);%the x coordinate of the ith boundary element
    w2 = VX(EToN(kb(i),side2Vert(sb(i),2)),:);%the x coordinate of the ith boundary element
    
    x1 = w1(1);y1=w1(2);
    x2 = w2(1);y2=w2(2);
    
%     switch initial_condition_string
% 
%         case'Problem1a' 
            if abs(x1+0.5)<1e-6  && abs(x2+0.5)<1e-6
                BLeft(l,:)=[kb(i) sb(i)];
                l = l+1;
            elseif abs(x1-0.5)<1e-6 && abs(x2-0.5)<1e-6
                BRight(r,:)=[kb(i) sb(i)];
                r = r+1;
            elseif abs(y1-0.5)<1e-6  && abs(y2-0.5)<1e-6
                BTop(t,:)=[kb(i) sb(i)];
                t = t+1;
            elseif abs(y1+0.5)<1e-6  && abs(y2+0.5)<1e-6
                BBottom(b,:)=[kb(i) sb(i)];
                b = b+1;
            end
                BExt(i,:)=[kb(i) sb(i)];

% 
%         case 'Problem1d'
%             BLeft = [];BRight=[];BTop=[];BBottom=[];
%             BExt(i,:)=[kb(i) sb(i)];
%         case 'Problem2'
%             
% %             if abs(x1+1)<1e-6  && abs(x2+1)<1e-6
% %                 BLeft(l,:)=[kb(i) sb(i)];
% %                 l = l+1;
% %             elseif abs(x1-1)<1e-6 && abs(x2-1)<1e-6
% %                 BRight(r,:)=[kb(i) sb(i)];
% %                 r = r+1;
% %             elseif abs(y1-1)<1e-6  && abs(y2-1)<1e-6
% %                 BTop(t,:)=[kb(i) sb(i)];
% %                 t = t+1;
% %             elseif abs(y1+1)<1e-6  && abs(y2+1)<1e-6
% %                 BBottom(b,:)=[kb(i) sb(i)];
% %                 b = b+1;
% %             end
% 
%             BExt(i,:)=[kb(i) sb(i)];
% case 'Problem2b'
%             BExt(i,:)=[kb(i) sb(i)];
% 
% case 'Problem3'
%             BExt(i,:)=[kb(i) sb(i)];
% 
%     end
end







