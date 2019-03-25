function [Kcc,Kce,Kce_pr,Kee, Kee_pr,Fe,Fe_pr,Fc,K_cond,F_cond, compliance_truss, u_cond, coupling_dofs]=truss_stiffness_condensation_by_set(nelx,nely,E,A1,A2,A3,A4,A5,I1,I2,I3,I4,I5,L)
Pov=[0 1 0 0 0 0
     1 0 0 0 0 0
     0 0 1 0 0 0
     0 0 0 0 1 0
     0 0 0 1 0 0
     0 0 0 0 0 1];% Rotation Matrix
 
% Stiffness Matrix by Set:
Kt=E*A1/L;% axial stiffness of all the beam
Kf=E*I1/L^3;% flexural behavior

Keh_1=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam

Kt=E*A2/L;% axial stiffness of all the beam
Kf=E*I2/L^3;% flexural behavior

Keh_2=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam

Kt=E*A3/L;% axial stiffness of all the beam
Kf=E*I3/L^3;% flexural behavior

Keh_3=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam
Kev_3=Pov'*Keh_3*Pov; % stiffness matrix vertical element

Kt=E*A4/L;% axial stiffness of all the beam
Kf=E*I4/L^3;% flexural behavior

Keh_4=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam
Kev_4=Pov'*Keh_4*Pov; % stiffness matrix vertical element

Kt=E*A5/L;% axial stiffness of all the beam
Kf=E*I5/L^3;% flexural behavior

Keh_5=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam
Kev_5=Pov'*Keh_5*Pov; % stiffness matrix vertical element

node_list = 1:(2*(nelx+1)+3*(nely-1));% 
node_dofs_map=reshape((1:3*(2*(nelx+1)+3*(nely-1)))',3,[])';
ELEMENT(1:nelx,:)=[(1:nelx)',(1:nelx)'+1];%Element connectivity matrix
ELEMENT(nelx + (1:nelx),:) = [nelx+1+(1:nelx)',(1:nelx)'+nelx+2];
ELEMENT(2*nelx+(1:(nely)),:)=[1,2*nelx+2+1
                              ((2*nelx+2+1):(2*nelx+2+nely-2))',((2*nelx+2+1):(2*nelx+2+nely-2))'+1
                                                    (2*nelx+2+nely-1), nelx+2];
ELEMENT(2*nelx+nely+(1:(nely)),:)=[fix(nelx/2),2*nelx+2+nely
                                                    ((2*nelx+2+nely):(2*nelx+2+2*nely-3))',((2*nelx+2+nely):(2*nelx+2+2*nely-3))'+1
                                                    (2*nelx+2+2*nely-2), nelx+1+fix(nelx/2)];
ELEMENT(2*nelx+2*nely+(1:(nely)),:)=[nelx+1,2*nelx+2+2*nely-2+1
                                                    ((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))',((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))'+1
                                                    (2*nelx+2+3*(nely-1)), 2*nelx+2];
ELEMENT_DOF=[node_dofs_map(ELEMENT(:,1),:),node_dofs_map(ELEMENT(:,2),:)]; % Local to global DOFs

%Local indexes and Non-zero terms in the stiffness matrix of horizontal and vertical elements
[io_1,jo_1,ko_1]=find(Keh_1);

[io_2,jo_2,ko_2]=find(Keh_2);

[iv_3,jv_3,kv_3]=find(Kev_3);

[iv_4,jv_4,kv_4]=find(Kev_4);

[iv_5,jv_5,kv_5]=find(Kev_5);

%Global indexes and Non-zero terms in the stiffness matrix of hor. elements
Io_1=reshape(ELEMENT_DOF(1:nelx,io_1)',[],1);
Jo_1=reshape(ELEMENT_DOF(1:nelx,jo_1)',[],1);

Io_2=reshape(ELEMENT_DOF(nelx+1:2*nelx,io_2)',[],1);
Jo_2=reshape(ELEMENT_DOF(nelx+1:2*nelx,jo_2)',[],1);

Ko=[repmat(ko_1,nelx,1);repmat(ko_2,nelx,1)];

%Global indexes and Non-zero terms in the stiffness matrix of ver. elements
Iv_1=reshape(ELEMENT_DOF(2*nelx+(1:nely),iv_3)',[],1);
Jv_1=reshape(ELEMENT_DOF(2*nelx+(1:nely),jv_3)',[],1);

Iv_2=reshape(ELEMENT_DOF(2*nelx+(nely+1:2*nely),iv_4)',[],1);
Jv_2=reshape(ELEMENT_DOF(2*nelx+(nely+1:2*nely),jv_4)',[],1);

Iv_3=reshape(ELEMENT_DOF(2*nelx+(2*nely+1:3*nely),iv_5)',[],1);
Jv_3=reshape(ELEMENT_DOF(2*nelx+(2*nely+1:3*nely),jv_5)',[],1);

Kv=[repmat(kv_3,nely,1);repmat(kv_4,nely,1);repmat(kv_5,nely,1)];

% Global Stiffness Matrix
K = sparse([Io_1; Io_2; Iv_1; Iv_2; Iv_3], [Jo_1; Jo_2; Jv_1; Jv_2; Jv_3], [Ko; Kv], max(node_dofs_map(:)), max(node_dofs_map(:))); %stiffness matrix assembly
K=(K+K')/2; % make it exactly symmetric

F = sparse(3*((nelx+1))+1,1,-1,max(node_dofs_map(:)),1);% Load Vector Assembly
coupling_dofs_pr=[(1:3:(3*(nelx+1)-2));2:3:(3*(nelx+1)-1)]; %u and v DOFs of the 1st hor. beam
coupling_dofs_pr = coupling_dofs_pr(:);

% coupling_dofs = [1 2 3; 3*nelx+1:3*nelx+3; 3*(nelx+1)+1:3*(nelx+1)+3; 2*3*(nelx+1)-2:2*3*(nelx+1);...
%     2*3*(nelx+1)+1:2*3*(nelx+1)+3; 2*3*(nelx+1)+3*(nely-1)-2:2*3*(nelx+1)+3*(nely-1);3*(nely-1)+(2*3*(nelx+1)+1):3*(nely-1)+(2*3*(nelx+1)+3);...
%     3*(nely-1)+(2*3*(nelx+1)+3*(nely-1))-2:3*(nely-1)+(2*3*(nelx+1)+3*(nely-1));...   
%     6*(nely-1)+(2*3*(nelx+1)+1):6*(nely-1)+(2*3*(nelx+1)+3);6*(nely-1)+(2*3*(nelx+1)+3*(nely-1))-2:6*(nely-1)+(2*3*(nelx+1)+3*(nely-1))];

coupling_dofs = [1 2; 3*fix((nelx+1)/2)+1:3*fix((nelx+1)/2)+2; 3*nelx+1:3*nelx+2];%; 3*(nelx+1)+1:3*(nelx+1)+3; 3*fix(3*(nelx+1)/4)+1:3*fix(3*(nelx+1)/4)+3; 2*3*(nelx+1)-2:2*3*(nelx+1)];

coupling_dofs = coupling_dofs(:);

reduced_dofs = setdiff(node_dofs_map(:),coupling_dofs); %Slave Nodes

reduced_dofs_pr = setdiff(node_dofs_map(:),coupling_dofs_pr);

ID_Load = find(reduced_dofs == (3*((nelx+1))+1)); %Force Node

%coupling_dofs = [coupling_dofs; reduced_dofs(ID_Load)];

%coupling DOFs = master nodes; Eliminated DOFs = slave nodes
u0=zeros(size(K,1),1);
u0(reduced_dofs) = K(reduced_dofs,reduced_dofs)\F(reduced_dofs);
Kcc = K(coupling_dofs,coupling_dofs);
Kce = K(coupling_dofs,reduced_dofs);
Kce_pr = K(coupling_dofs_pr,reduced_dofs_pr);
Kee = K(reduced_dofs,reduced_dofs);
Kee_pr = K(reduced_dofs_pr,reduced_dofs_pr);
Fc = F(coupling_dofs);
Fe = F(reduced_dofs);
Fe_pr = F(reduced_dofs_pr);
B = (Kce/Kee);
%For the Condensation
K_cond = Kcc - B*Kce';
F_cond = Fc - B*Fe;
u_cond = K_cond\F_cond;
B_t = B';
compliance_truss = -1*(u0(ID_Load)+ B_t(ID_Load,:)*u_cond)/10e11; 

coupling_dofs = [1 2; 3*fix((nelx+1)/2)+1 -  fix((3*fix((nelx+1)/2)+1)/3) :3*fix((nelx+1)/2)+2 - fix((3*fix((nelx+1)/2)+1)/3); 3*nelx+1 - fix((3*nelx+1)/3):3*nelx+2 - fix((3*nelx+1)/3)];

% Plotting
Xo=[(0:nelx)';(0:nelx)'];
Yo=[zeros(nelx+1,1);repmat(-nely,nelx+1,1)];
Xv=[zeros(nely-1,1);repmat(fix(nelx/2)-1,nely-1,1);repmat(nelx,nely-1,1)];
Yv=repmat((-1:-1:-nely+1)',3,1);
X=[Xo;Xv];
Y=[Yo;Yv];
plot(X(ELEMENT)',Y(ELEMENT)','k-o','MarkerFaceColor','k')
axis equal
U=reshape(u0,3,[])';
U=U(:,1:2);
DX=U(:,1);
DY=U(:,2);
hold on
