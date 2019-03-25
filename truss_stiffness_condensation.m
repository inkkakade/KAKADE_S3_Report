function [Kcc,Kce,Kee,Fe,Fc,K_cond,F_cond, compliance_truss, u_cond,B]=truss_stiffness_condensation(nelx,nely,E,A,I,L)
Kt=E*A/L;% axial stiffness of all the beam
Kf=E*I/L^3;% flexural behavior 
Keh=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam
Pov=[0 1 0 0 0 0
     1 0 0 0 0 0
     0 0 1 0 0 0
     0 0 0 0 1 0
     0 0 0 1 0 0
     0 0 0 0 0 1];% Rotation Matrix
Kev=Pov'*Keh*Pov; % stiffness matrix vertical element
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
[io,jo,ko]=find(Keh);
[iv,jv,kv]=find(Kev);

%Global indexes and Non-zero terms in the stiffness matrix of hor. elements
Io=reshape(ELEMENT_DOF(1:2*nelx,io)',[],1);
Jo=reshape(ELEMENT_DOF(1:2*nelx,jo)',[],1);
Ko=[repmat(ko,nelx,1);repmat(ko,nelx,1)];

%Global indexes and Non-zero terms in the stiffness matrix of ver. elements
Iv=reshape(ELEMENT_DOF(2*nelx+(1:3*nely),iv)',[],1);
Jv=reshape(ELEMENT_DOF(2*nelx+(1:3*nely),jv)',[],1);
Kv=repmat(kv,3*nely,1);

% Global Stiffness Matrix
K = sparse([Io;Iv],[Jo;Jv],[Ko;Kv], max(node_dofs_map(:)), max(node_dofs_map(:))); %stiffness matrix assembly
K=(K+K')/2; % make it exactly symmetric

F = sparse(3*((nelx+1))+1,1,-1,max(node_dofs_map(:)),1);% Load Vector Assembly
coupling_dofs=[(1:3:(3*(nelx+1)-2));2:3:(3*(nelx+1)-1)]; %u and v DOFs of the 1st hor. beam
coupling_dofs = coupling_dofs(:);

reduced_dofs = setdiff(node_dofs_map(:),coupling_dofs); %Slave Nodes

ID_Load = find(reduced_dofs == (3*((nelx+1))+1)); %Force Node

%coupling_dofs = [coupling_dofs; reduced_dofs(ID_Load)];

%coupling DOFs = master nodes; Eliminated DOFs = slave nodes
u0=zeros(size(K,1),1);
u0(reduced_dofs) = K(reduced_dofs,reduced_dofs)\F(reduced_dofs);
Kcc = K(coupling_dofs,coupling_dofs);
Kce = K(coupling_dofs,reduced_dofs);
Kee = K(reduced_dofs,reduced_dofs);
Fc = F(coupling_dofs);
Fe = F(reduced_dofs);
B = (Kce/Kee);
%For the Condensation
K_cond = Kcc - B*Kce';
F_cond = Fc - B*Fe;
u_cond = K_cond\F_cond;
B_t = B';
compliance_truss = 1*(u0(ID_Load) + B_t(ID_Load,:)*u_cond)/10e12; 

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
