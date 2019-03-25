function [dK_cond_dA, Kee_a, B_a] = area_cond(nelx,nely,E,L,Kce,Kee,B)

Kt_a = E/L;
Kf_a = 0;
Keo=[Kt_a 0       0        -Kt_a 0       0
     0   12*Kf_a   6*L*Kf_a   0   -12*Kf_a  6*L*Kf_a
     0   6*L*Kf_a  4*L^2*Kf_a 0   -6*L*Kf_a 2*L^2*Kf_a
     -Kt_a 0       0        Kt_a  0       0
     0   -12*Kf_a  -6*L*Kf_a  0   12*Kf_a   -6*L*Kf_a
     0   6*L*Kf_a  2*L^2*Kf_a 0   -6*L*Kf_a 4*L^2*Kf_a];
Pov=[0 1 0 0 0 0
    1 0 0 0 0 0
    0 0 1 0 0 0
    0 0 0 0 1 0
    0 0 0 1 0 0
    0 0 0 0 0 1];
Kev=Pov'*Keo*Pov;
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
[io,jo,ko]=find(Keo);
[iv,jv,kv]=find(Kev);
Io=reshape(ELEMENT_DOF(1:2*nelx,io)',[],1);
Jo=reshape(ELEMENT_DOF(1:2*nelx,jo)',[],1);
Ko =[repmat(ko,nelx,1)/4;4*repmat(ko,nelx,1)];
Iv = reshape(ELEMENT_DOF(2*nelx+(1:3*nely),iv)',[],1);
Jv = reshape(ELEMENT_DOF(2*nelx+(1:3*nely),jv)',[],1);
Kv = 4*repmat(kv,3*nely,1);
K=sparse([Io;Iv],[Jo;Jv],[Ko;Kv],max(node_dofs_map(:)),max(node_dofs_map(:)));
K=(K+K')/2;
coupling_dofs=[(1:3:(3*(nelx+1)-2));2:3:(3*(nelx+1)-1)];
coupling_dofs=coupling_dofs(:);
reduced_dofs=setdiff(node_dofs_map(:),coupling_dofs);
Kcc_a = K(coupling_dofs,coupling_dofs);
Kce_a = K(coupling_dofs,reduced_dofs);
Kee_a = K(reduced_dofs,reduced_dofs);
%B_a = Kee_a\Kce' + Kee\Kce_a';
B_a = Kee\Kce_a' - (Kee^-1)*transpose(Kce*(Kee^-1)); %inv?(A)B=?A^(?1)BA^(?1)
dK_cond_dA = Kcc_a - B_a'*Kce' + Kce_a*B';
