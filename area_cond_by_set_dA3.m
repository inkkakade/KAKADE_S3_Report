function [Kee_a3, B_a3] = area_cond_by_set_dA3(nelx,nely,E,L,Kce,Kee)

Pov=[0 1 0 0 0 0
     1 0 0 0 0 0
     0 0 1 0 0 0
     0 0 0 0 1 0
     0 0 0 1 0 0
     0 0 0 0 0 1];% Rotation Matrix
 
% Stiffness Matrix by Set:
%% For A3:
Kt=0;% axial stiffness of all the beam
Kf=0;% flexural behavior

Keh_1=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam

Kt=0;% axial stiffness of all the beam
Kf=0;% flexural behavior

Keh_2=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam

Kt=E/L;% axial stiffness of all the beam
Kf=0;% flexural behavior

Keh_3=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam
Kev_3=Pov'*Keh_3*Pov; % stiffness matrix vertical element

Kt=0;% axial stiffness of all the beam
Kf=0;% flexural behavior

Keh_4=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam
Kev_4=Pov'*Keh_4*Pov; % stiffness matrix vertical element

Kt=0;% axial stiffness of all the beam
Kf=0;% flexural behavior

Keh_5=[Kt 0 0 -Kt 0 0
    0 12*Kf 6*L*Kf 0 -12*Kf 6*L*Kf
    0 6*L*Kf 4*L^2*Kf 0 -6*L*Kf 2*L^2*Kf
    -Kt 0 0 Kt 0 0
    0 -12*Kf -6*L*Kf 0 12*Kf -6*L*Kf
    0 6*L*Kf 2*L^2*Kf 0 -6*L*Kf 4*L^2*Kf]; % Element stiffness matrix for Horizontal beam
Kev_5=Pov'*Keh_5*Pov; % stiffness matrix vertical element

node_list=1:(2*(nelx+1)+3*(nely-1));
node_dofs_map=reshape((1:3*(2*(nelx+1)+3*(nely-1)))',3,[])';
ELEMENT(1:nelx,:)=[(1:nelx)',(1:nelx)'+1];
ELEMENT(nelx+(1:nelx),:)=[nelx+1+(1:nelx)',(1:nelx)'+nelx+2];
ELEMENT(2*nelx+(1:(nely)),:)=[1,2*nelx+2+1
                              ((2*nelx+2+1):(2*nelx+2+nely-2))',((2*nelx+2+1):(2*nelx+2+nely-2))'+1
                                                    (2*nelx+2+nely-1), nelx+2];
ELEMENT(2*nelx+nely+(1:(nely)),:)=[fix(nelx/2),2*nelx+2+nely
                                                    ((2*nelx+2+nely):(2*nelx+2+2*nely-3))',((2*nelx+2+nely):(2*nelx+2+2*nely-3))'+1
                                                    (2*nelx+2+2*nely-2), nelx+1+fix(nelx/2)];
ELEMENT(2*nelx+2*nely+(1:(nely)),:)=[nelx+1,2*nelx+2+2*nely-2+1
                                                    ((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))',((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))'+1
                                                    (2*nelx+2+3*(nely-1)), 2*nelx+2];
ELEMENT_DOF=[node_dofs_map(ELEMENT(:,1),:),node_dofs_map(ELEMENT(:,2),:)];

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
K=(K+K')/2;

coupling_dofs = [1 2; 3*fix((nelx+1)/2)+1:3*fix((nelx+1)/2)+2; 3*nelx+1:3*nelx+2];%; 3*(nelx+1)+1:3*(nelx+1)+3; 3*fix(3*(nelx+1)/4)+1:3*fix(3*(nelx+1)/4)+3; 2*3*(nelx+1)-2:2*3*(nelx+1)];
coupling_dofs = coupling_dofs(:);

reduced_dofs=setdiff(node_dofs_map(:),coupling_dofs);
Kce_a = K(coupling_dofs,reduced_dofs);
Kee_a3 = K(reduced_dofs,reduced_dofs);
%B_a = Kee_a\Kce' + Kee\Kce_a';
B_a3 = Kee\Kce_a' - (Kee^-1)*transpose(Kce*(Kee^-1)); %inv?(A)B=?A^(?1)BA^(?1)
