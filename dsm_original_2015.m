
        %%%%%%%%%% STIFFNESS METHOD %%%%%%%%%
        %%% Carlo Ruiz & Elizabeth Negrete %%
        %%%%%%%%%%%%%%  May 2015 %%%%%%%%%%%%

%%%%%%%%%%%%% Definitions %%%%%%%%%%%
clc;clear;close all
B=xlsread('Stiffness Method.xlsx','BARS');          % reads excel for bar info
Y=xlsread('Stiffness Method.xlsx','JOINTS');        % reads excel for joint coordinates
F=xlsread('Stiffness Method.xlsx','FORCES');        % reads excel for force vectors
NB=size(B,1)/2;                     % number of bars
NJ=size(Y,1);                       % number of joints
M{1,NB}=[];                         % cell for each bar
m{1,NB}=[];                         % bar's vector
L{1,NB}=[];                         % Length of bar
a{1,NB}=[];                         % angle of vector with respect to x axis
T{1,NB}=[];                         % Transformation matrix
kLAA{1,NB}=[];                      % Local AA stiffness matrix
kLAB{1,NB}=[];                      % Local AB stiffness matrix
kLBA{1,NB}=[];                      % Local BA stiffness matrix
kLBB{1,NB}=[];                      % Local BB stiffness matrix
kGAA{1,NB}=[];                      % Global AA stiffness matrix
kGAB{1,NB}=[];                      % Global AB stiffness matrix
kGBA{1,NB}=[];                      % Global BA stiffness matrix
kGBB{1,NB}=[];                      % Global BB stiffness matrix
K{NJ,NJ}=[];                        % Global system stiffness matrix
d{1,NJ}=[];                         % Global displacements vector
Fa{1,NB}=[];                        % Forces at bar's start
Fb{1,NB}=[];                        % Forces at bar's end
da{1,NB}=[];                        % Displacements assigned ton bar's start
db{1,NB}=[];                        % Displacements assigned to bar's end

%%%%%%%%%%%%%% Assigns Data %%%%%%%%%%%%
for i=1:NB;
M{i}=B(2*i-1:2*i,:);                       
m{i}=M{i}(1:2,2)-M{i}(1:2,1);              
L{i}=norm(m{i});                           
a{i}=[atan2(m{i}(2,1)/m{i}(1,1),1)];
kLAA{i}=[M{i}(1,3)*M{i}(1,4)./L{i} 0 0;0 12*M{i}(2,3)*M{i}(1,4)./L{i}^3 6*M{i}(2,3)*M{i}(1,4)./L{i}^2; 0 6*M{i}(2,3)*M{i}(1,4)./L{i}^2 4*M{i}(2,3)*M{i}(1,4)./L{i}];
kLAB{i}=[-M{i}(1,3)*M{i}(1,4)./L{i} 0 0;0 -12*M{i}(2,3)*M{i}(1,4)./L{i}^3 6*M{i}(2,3)*M{i}(1,4)./L{i}^2; 0 -6*M{i}(2,3)*M{i}(1,4)./L{i}^2 2*M{i}(2,3)*M{i}(1,4)./L{i}];
kLBA{i}=kLAB{i}';
kLBB{i}=[M{i}(1,3)*M{i}(1,4)./L{i} 0 0;0 12*M{i}(2,3)*M{i}(1,4)./L{i}^3 -6*M{i}(2,3)*M{i}(1,4)./L{i}^2; 0 -6*M{i}(2,3)*M{i}(1,4)./L{i}^2 4*M{i}(2,3)*M{i}(1,4)./L{i}];
T{i}=[cos(a{i}) sin(a{i}) 0;-sin(a{i}) cos(a{i}) 0;0 0 1];
kGAA{i}=T{i}'*kLAA{i}*T{i};
kGAB{i}=T{i}'*kLAB{i}*T{i};
kGBA{i}=T{i}'*kLBA{i}*T{i};
kGBB{i}=T{i}'*kLBB{i}*T{i};
end
%%%%%%%%%% Assembles Diagonal %%%%%%%%%
for i=1:NJ;
    J{i}=Y(i,:)';
    K{i,i}=[0 0 0;0 0 0;0 0 0];
    for j=1:NB;
        if J{i}==M{j}(:,1);
            K{i,i}=K{i,i}+kGAA{j};
            else if J{i}==M{j}(:,2);
            K{i,i}=K{i,i}+kGBB{j};
        else
           K{i,i}=K{i,i};            
                end
        end
    end
end
%%%%%%% Assembles rest of Matrix %%%%%%
for k=1:NJ;
    for j=k+1:NJ;
        K{j,k}=[0 0 0;0 0 0;0 0 0];
        for l=1:NB; 
            if J{k}==M{l}(:,1) & J{j}==M{l}(:,2)
               K{j,k}=K{j,k}+kGBA{l};
            else if J{k}==M{l}(:,2) &  J{j}==M{l}(:,1);
                    K{j,k}=K{j,k}+kGAB{l};
                else
                    K{j,k}=K{j,k};
                end
            end
        end
    K{k,j}=K{j,k}';
    end
end
%%%% Calculates global joint displacements %%%%
dG=inv(cell2mat(K))*F;
for i=1:NJ;
    d{i}=dG(3*i-2:3*i,1);
end
%%%%%%%%% Calculates forces on beam's extrema %%%%%%%%%%
for i=1:NB;
    da{i}=[0 0 0]';
    db{i}=[0 0 0]';
    for l=1:NJ;
        if J{l}==M{i}(:,1);
            da{i}=da{i}+d{l};
        end   
        if J{l}==M{i}(:,2);
                db{i}=db{i}+d{l};
        end
    Fa{i}=kLAA{i}*T{i}*da{i}+kLAB{i}*T{i}*db{i};
    Fb{i}=kLBA{i}*T{i}*da{i}+kLBB{i}*T{i}*db{i};   
    end
end
%%%%%%%%%% Writes results in excel %%%%%%%%%%
xlswrite('Stiffness Method.xlsx',dG,'RESULTS','I6');
xlswrite('Stiffness Method.xlsx',cell2mat(K),'RESULTS','M6');
xlswrite('Stiffness Method.xlsx',cell2mat(Fa(:)),'RESULTS','D6');
xlswrite('Stiffness Method.xlsx',cell2mat(Fb(:)),'RESULTS','E6');
%%%%%%%%%%%%%% END OF PROGRAM %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
