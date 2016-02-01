clear all; close all;

g=0.3347; delta=8.1831-10.5665; kappa=0.0024; Sigma_z=-1; Omega_c=10.5665;
rr=-0.1:0.01:1.1;
chi0=g.^2./abs(delta);

Omega_d=Omega_c+rr.*chi0;
xiC1= (abs(delta)*kappa)^(3/2)/(3^(3/4)*g^2);
AA=0.1:0.1:100;
[Omega_d,A]=meshgrid(Omega_d,AA);

chi=Sigma_z.*g.^2./(sqrt(2.*g.^2.*(A.^2 + Sigma_z) + delta.^2));
%xi=(A.^2).*((Omega_d.^2 - (Omega_c - chi).^2).^2 + kappa.^2.*Omega_d.^2);
xi=(1/Omega_c)*A.*sqrt((Omega_d.^2 - (Omega_c - chi).^2).^2 + kappa.^2.*Omega_d.^2);
%xi_Duff=(A.*sqrt(((Omega_c-(Sigma_z.*g.^2)./(abs(delta).*(1.0./delta.^2.*g.^2.*(Sigma_z+A.^2)+1.0))).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2))./Omega_c;

dxidA=sqrt(((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2)./Omega_c+(A.^2.*Sigma_z.*g.^4.*(Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).*1.0./sqrt(((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2).*1.0./(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2).^(3.0./2.0).*((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).*4.0)./Omega_c;

%dxidA_Duff=sqrt(((Omega_c-(Sigma_z.*g.^2)./(abs(delta).*(1.0./delta.^2.*g.^2.*(Sigma_z+A.^2)+1.0))).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2)./Omega_c+(A.^2.*Sigma_z.*1.0./abs(delta).^3.*g.^4.*1.0./sqrt(((Omega_c-(Sigma_z.*g.^2)./(abs(delta).*(1.0./delta.^2.*g.^2.*(Sigma_z+A.^2)+1.0))).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2).*1.0./(1.0./delta.^2.*g.^2.*(Sigma_z+A.^2)+1.0).^2.*((Omega_c-(Sigma_z.*g.^2)./(abs(delta).*(1.0./delta.^2.*g.^2.*(Sigma_z+A.^2)+1.0))).^2-Omega_d.^2).*(Omega_c-(Sigma_z.*g.^2)./(abs(delta).*(1.0./delta.^2.*g.^2.*(Sigma_z+A.^2)+1.0))).*4.0)./Omega_c;


%dxidA_C2=sqrt(Omega_d.^2.*kappa.^2+((Omega_c-(sqrt(2.0).*Sigma_z.*g.*(1.0./2.0))./A).^2-Omega_d.^2).^2)./Omega_c+(sqrt(2.0).*Sigma_z.*g.*1.0./sqrt(Omega_d.^2.*kappa.^2+((Omega_c-(sqrt(2.0).*Sigma_z.*g.*(1.0./2.0))./A).^2-Omega_d.^2).^2).*((Omega_c-(sqrt(2.0).*Sigma_z.*g.*(1.0./2.0))./A).^2-Omega_d.^2).*(Omega_c-(sqrt(2.0).*Sigma_z.*g.*(1.0./2.0))./A))./(A.*Omega_c);

%dxi2dA2=(A.^3.*Sigma_z.^2.*g.^8.*(Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2.*1.0./sqrt(((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2).*1.0./(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2).^3.*1.6e1)./Omega_c+(A.^3.*Sigma_z.^2.*g.^8.*1.0./sqrt(((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2).*1.0./(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2).^3.*((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).*8.0)./Omega_c+(A.*Sigma_z.*g.^4.*(Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).*1.0./sqrt(((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2).*1.0./(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2).^(3.0./2.0).*((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).*1.2e1)./Omega_c-(A.^3.*Sigma_z.*g.^6.*(Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).*1.0./sqrt(((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2).*1.0./(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2).^(5.0./2.0).*((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).*2.4e1)./Omega_c-(A.^3.*Sigma_z.^2.*g.^8.*(Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2.*1.0./(((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2+Omega_d.^2.*kappa.^2).^(3.0./2.0).*1.0./(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2).^3.*((Omega_c-Sigma_z.*g.^2.*1.0./sqrt(g.^2.*(Sigma_z+A.^2).*2.0+delta.^2)).^2-Omega_d.^2).^2.*1.6e1)./Omega_c;

figure; contour(Omega_d,A,xi,100);
figure; mesh(Omega_d,xi,A);
v=[0,0];
figure;contour(Omega_d,xi/xiC1,dxidA, v, 'b'); set(gca,'yscale','log');
xlabel('Drive frequency (GHz)','FontSize',15); ylabel('\epsilon/\epsilon_{C1}', 'FontSize',15); title('Bimodality Leaf','FontSize',20);
%contour(Omega_d,A,dxi2dA2, v, 'b');set(gca,'yscale','log'); hold off;

%ylim([1.1 0.7*10^2]);
