%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the electric field and ohmic loss distribution in metal-
% dielectric impedance matched metamaterials with a moth-eye surface, as shown in 
% Fig. 9c of R. Contractor, G. D’Aguanno, C. Menyuk, “Ultra-broadband, polarization-
% independent, wide-angle absorption in impedance-matched metamaterials with anti-
% reflective moth-eye surfaces”, Optics Express.

% The RCWA implementation is based on:
% M. G. Moharam and T. K. Gaylord, “Rigorous coupled-wave analysis of planar-grating 
% diffraction,” J. Opt. Soc. Am. 71, 811-818 (1995).
% L. Li, “New formulation of the Fourier modal method for crossed surface-relief 
% gratings,” J. Opt. Soc. Am. A 14, 2758-2767 (1997).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Inputs and Structure Parameters

y = 5000; % Range of wavelength in nm

% Angles
th = 0; % Incidence Angle (deg)
phi = 0; % Angle between 'x' axis and plane of incidence (deg)
psi = 90; % Polarization angle (deg) 0 -> TM, 90 -> TE

% Unit cell parameters
mlt = 2.78; % Scaling factor for the hexagonal unit cell
D_x = 300; % Length of one period in 'x'
D_x = 2*round(D_x*mlt/2); % Scaled and rounded to nearest even number
D_y = 2*round(sqrt(3)*D_x/2); % Length of one period in 'y'

% Metal-dielectric metamaterials
D_m = 5; % Thickness of metal layer in 'z'
D_d = 75; % Thickness of dielectric layer in 'z'
NN_z = 6; % =10; for uncarved % Number of uncarved unit cells (after the moth-eyes)

% Moth-eye Structures
N_m = 4; % =0; for uncarved % Number of elementary cells carved into pillars
r_1 = 85*mlt; % Radius of pillars at top
r_2 = 138*mlt; % Radius of pillars at bottom
N_zd = 3; % Number of slices of the dielectric layer (calculates radius, increase for accuracy)
D_z = N_m*(2*D_d+D_m); % Thickness of moth-eye structure in 'z'
N_z = N_m*(2*N_zd+1); % Number of layers in one carved elementary cell
h_zd = D_d/N_zd; % Resolution in 'z' for carved dielectric layers
h_z = repmat([repmat(h_zd,1,N_zd),D_m,repmat(h_zd,1,N_zd)],1,N_m); % Height of each carved layer
z_l = [0,cumsum(h_z(1:N_z-1))]; % 'z' coordinate of each carved layer
r_z = interp1([0,D_z-h_zd],[r_1,r_2],z_l); % Radius of pillars for each carved layer
N_l = N_z+3*NN_z+2; % Total layers

% Fourier modes to consider in calculations (ODD)
N_x = 7; %=1; for uncarved % Number of Fourier modes in x direction 
N_y = 7; %=1; for uncarved % Number of Fourier modes in y direction
Nt = N_x*N_y;

% Preallocations
ky = 2*pi./y; % Wavenumbers
I = eye(Nt); II = eye(2*Nt); % Identity matrices
ZZ = zeros(2*Nt); zz = zeros(2*Nt,1);
EPS = zeros(Nt,Nt,N_l);
EPSyx = zeros(Nt,Nt,N_z);
EPSxy = zeros(Nt,Nt,N_z);
W = zeros(2*Nt,2*Nt,N_l);
V = zeros(2*Nt,2*Nt,N_l);
X = zeros(2*Nt,2*Nt,N_l);
qm = zeros(2*Nt,N_l);

% Foruier mode numbers and arrays
mx = (-(N_x-1)/2:(N_x-1)/2)';
my = (-(N_y-1)/2:(N_y-1)/2);
m0x = (N_x+1)/2;
m0y = (N_y+1)/2;
m0 = (Nt+1)/2;

%% Wavelength based permittivity

% Metal (with thin layer correction)
load Tu.txt; % Taken from refractiveindex.info
Me = interp1(Tu(:,1),Tu,y/1000); % Interpolation of experimental data
Metal = (Me(:,2) + 1i*Me(:,3)).^2; Metal = Metal.';
kkk = D_m/19; % Thickness of W layer divided mean free path of W
kky = 1000*ky;
fi = (1./kkk+3/4*(1-1/12*kkk.^2).*expint(kkk)-3./(8*kkk.^2).*(1-exp(-kkk))-(5./(8*kkk)+1/16-kkk/16).*exp(-kkk)).^-1;

wp = 0.79; % Plasma frequency in micron^-1
gamma = 39.79*10^-4; gammakk = gamma*fi/kkk;
epsDrude = 1-(wp^2)./((kky/(2*pi)).^2+1i*(kky/(2*pi))*gamma); % Drude model for the Bulk

epsDrudekk = 1-(wp^2)./((kky/(2*pi)).^2+1i*(kky/(2*pi))*gammakk); % Drude Model corrected for the thin layer

epsMetal = Metal-epsDrude+epsDrudekk;

% Dielectric
load Ta2O5_n.txt; load Ta2O5_k.txt; % Taken from refractiveindex.info
De = interp1(Ta2O5_n(:,1),[Ta2O5_n(:,2),Ta2O5_k(:,2)],y/1000);
eps_d = (De(:,1) + 1i*De(:,2)).^2;

%% Index profile of moth eye structures
% Build the moth-eye structure as a binary 3D matrix

indexProfLog = false(D_x,D_y,N_z);

% Hexagon vertice positions
centre = [D_x,D_y]/2;

for nn=1:N_z
    r = round(r_z(nn));
    % Pillars
    [rr,cc] = meshgrid(1:2*r+1);
    circ = (rr-1-r).^2+(cc-1-r).^2<=r^2;
    % Positioning the pillars
    indexProfLog((centre(1)-r):(centre(1)+r),(centre(2)-r):(centre(2)+r),nn) = circ;
    indexProfLog(1:r,1:r,nn) = indexProfLog(1:r,1:r,nn) | circ(r+2:2*r+1,r+2:2*r+1);
    indexProfLog(1:r,D_y-r:D_y,nn) = indexProfLog(1:r,D_y-r:D_y,nn) | circ(r+2:2*r+1,1:r+1);
    indexProfLog(D_x-r:D_x,D_y-r:D_y,nn) = indexProfLog(D_x-r:D_x,D_y-r:D_y,nn) | circ(1:r+1,1:r+1);
    indexProfLog(D_x-r:D_x,1:r,nn) = indexProfLog(D_x-r:D_x,1:r,nn) | circ(1:r+1,r+2:2*r+1);
end

%% Structure Analysis

uy = y/1000; % Wavelength in micron

% Permittivity and refractive indices
n_1 = 1; % Refractive index of incidence medium
n_3 = 1; % Refractive index of transmission medium
e_m = epsMetal; % Permittivty of metal for current wavelength
e_d = eps_d; % Permittivity of dielectric layer (Ta2O5) 
e_f = 1; % Permittivity of material between pillars
e_p = repmat([repmat(conj(e_d),N_zd,1);conj(e_m);repmat(conj(e_d),N_zd,1)],N_m,1);

% Permittivity matrices for each layer
indexProf = ones(D_x,D_y);
for nn=1:N_z
    indexProf(indexProfLog(:,:,nn)) = e_p(nn);
    indexProf(~indexProfLog(:,:,nn)) = e_f;

    eps_fft = fftshift(fft2(indexProf))/(D_x*D_y);
    eps_mn = eps_fft(D_x/2+2-N_x:D_x/2+N_x,D_y/2+2-N_y:D_y/2+N_y);
    
    for pp = 1:N_x
        for qq = 1:N_y
            EE = rot90(eps_mn(pp:pp+N_x-1,qq:qq+N_y-1),2);
            EPS(pp+N_x*(qq-1),:,nn+1) = reshape(EE,1,[]);
        end
    end
    
    i_iepsx_mj = zeros(N_x,D_y,N_x);
    iepsx_fft = fftshift(fft(indexProf.^(-1),[],1),1)/D_x;
    
    for qq=1:D_y
        iepsx_m = iepsx_fft(D_x/2+2-N_x:D_x/2+N_x,qq);
        iepsx_mj = toeplitz(iepsx_m(N_x:2*N_x-1),flip(iepsx_m(1:N_x)));
        
        i_iepsx_mj(:,qq,:) = inv(iepsx_mj);
    end
    
    epsxy_fft = fftshift(fft(i_iepsx_mj,[],2),2)/D_y;
    epsxy_mnj = epsxy_fft(:,D_y/2+2-N_y:D_y/2+N_y,:);
    
    E4 = zeros(N_x,N_y,N_x,N_y);
    
    for pp = 1:N_x
        for qq = 1:N_x
            E4(pp,:,qq,:) = toeplitz(epsxy_mnj(pp,N_y:2*N_y-1,qq),flip(epsxy_mnj(pp,1:N_y,qq)));
        end
    end
    
    EPSxy(:,:,nn) = reshape(E4,[N_x*N_y,N_x*N_y]);
    
    i_iepsy_nl = zeros(D_x,N_y,N_y);
    iepsy_fft = fftshift(fft(indexProf.^(-1),[],2),2)/D_y;
    
    for pp=1:D_x
        iepsy_n = iepsy_fft(pp,D_y/2+2-N_y:D_y/2+N_y);
        iepsy_nl = toeplitz(iepsy_n(N_y:2*N_y-1),flip(iepsy_n(1:N_y)));
        
        i_iepsy_nl(pp,:,:) = inv(iepsy_nl);
    end
    
    epsyx_fft = fftshift(fft(i_iepsy_nl,[],1),1)/D_x;
    epsyx_mnl = epsyx_fft(D_x/2+2-N_x:D_x/2+N_x,:,:);

    E4 = zeros(N_x,N_y,N_x,N_y);

    for rr = 1:N_y
        for ss = 1:N_y
            E4(:,rr,:,ss) = toeplitz(epsyx_mnl(N_x:2*N_x-1,rr,ss),flip(epsyx_mnl(1:N_x,rr,ss)));
        end
    end
    
    EPSyx(:,:,nn) = reshape(E4,[N_x*N_y,N_x*N_y]);
end

%% Initializing variables
% Incident Field
u_x = cosd(psi)*cosd(phi)*cosd(th) - sind(psi)*sind(phi);
u_y = cosd(psi)*sind(phi)*cosd(th) + sind(psi)*cosd(phi);
inc = zz; inc(m0) = u_x; inc(Nt+m0) = u_y;

% Wavenumber Matrices
k_0 = ky;

k_xi = k_0*n_1*sind(th)*cosd(phi) + 2*pi*mx/D_x;
k_yi = k_0*n_1*sind(th)*sind(phi) + 2*pi*my/D_y;
k_x_mn = reshape(repmat(k_xi,1,N_y),[],1);
k_y_mn = reshape(repmat(k_yi,N_x,1),[],1);

k_1_zmn = conj(((n_1*k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));
k_d_zmn = conj((e_d*(k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));
k_m_zmn = conj((e_m*(k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));
k_3_zmn = conj(((n_3*k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));

K_x = 1i*diag(k_x_mn/k_0);
K_y = 1i*diag(k_y_mn/k_0);

%% Eigenmatrices
for nn=1:N_z
    % Eigenvalue formulation for the grating region
    F11 = -K_x*(EPS(:,:,nn+1)\K_y);
    F12 = I + K_x*(EPS(:,:,nn+1)\K_x);
    F21 = -I - K_y*(EPS(:,:,nn+1)\K_y);
    F22 = K_y*(EPS(:,:,nn+1)\K_x);
    F = [F11 F12;F21 F22];
    
    G11 = -K_x*K_y;
    G12 = EPSyx(:,:,nn) + K_x^2;
    G21 = -EPSxy(:,:,nn) - K_y^2;
    G22 = K_x*K_y;
    G = [G11 G12;G21 G22];
    
    [W(:,:,nn+1), Q] = eig(F*G);
    Q = (sqrt(Q));
    qm(:,nn+1) = diag(Q);
    V(:,:,nn+1) = -G*W(:,:,nn+1)/Q;
    X(:,:,nn+1) = diag(exp(-k_0*qm(:,nn+1)*h_z(nn)));
end

% Eigenmatrices for reflection and transmission layers
W(:,:,1) = II;
W(:,:,N_l) = II;

EPS(:,:,1) = n_1^2*I;
VIxx = (1i/k_0)*diag((k_x_mn.*k_y_mn)./k_1_zmn);
VIxy = (1i/k_0)*diag((k_y_mn.^2 + k_1_zmn.^2)./k_1_zmn);
VIyx = -(1i/k_0)*diag((k_x_mn.^2 + k_1_zmn.^2)./k_1_zmn);
VIyy = -VIxx;
V(:,:,1) = [VIxx VIxy;VIyx VIyy];
qm(:,1) = 1i*[k_1_zmn;k_1_zmn]/k_0;
X(:,:,1) = II;

EPS(:,:,N_l) = n_3^2*I;
VIIxx = (1i/k_0)*diag((k_x_mn.*k_y_mn)./k_3_zmn);
VIIxy = (1i/k_0)*diag((k_y_mn.^2 + k_3_zmn.^2)./k_3_zmn);
VIIyx = -(1i/k_0)*diag((k_x_mn.^2 + k_3_zmn.^2)./k_3_zmn);
VIIyy = -VIIxx;
V(:,:,N_l) = [VIIxx VIIxy;VIIyx VIIyy];
qm(:,N_l) = 1i*[k_3_zmn;k_3_zmn]/k_0;
X(:,:,N_l) = II;

% Eigenmatrices for the uncarved layers
for nn=1:3*NN_z
    if (mod(nn+1,3))==0
        k_n_zmn = k_m_zmn;
        D_n = D_m;
        EPS(:,:,1+N_z+nn) = conj(e_m)*I;
    else
        k_n_zmn = k_d_zmn;
        D_n = D_d;
        EPS(:,:,1+N_z+nn) = conj(e_d)*I;
    end
    W(:,:,1+N_z+nn) = II;
    Vnxx = (1i/k_0)*diag((k_x_mn.*k_y_mn)./k_n_zmn);
    Vnxy = (1i/k_0)*diag((k_y_mn.^2 + k_n_zmn.^2)./k_n_zmn);
    Vnyx = -(1i/k_0)*diag((k_x_mn.^2 + k_n_zmn.^2)./k_n_zmn);
    Vnyy = -Vnxx;
    V(:,:,1+N_z+nn) = [Vnxx Vnxy;Vnyx Vnyy];
    qm(:,1+N_z+nn) = 1i*[k_n_zmn;k_n_zmn]/k_0;
    X(:,:,1+N_z+nn) = diag(exp(-1i*[k_n_zmn;k_n_zmn]*D_n));
end

%% Big matrix calculation

% Preallocation
B = repmat(ZZ,2*(N_l+1),2*(N_l+1));

% Building B using b
for ll=1:N_l
    b = [W(:,:,ll) W(:,:,ll)*X(:,:,ll);
        V(:,:,ll) -V(:,:,ll)*X(:,:,ll);
        -W(:,:,ll)*X(:,:,ll) -W(:,:,ll);
        -V(:,:,ll)*X(:,:,ll) V(:,:,ll)];
    
    B(4*(ll-1)*Nt+1:4*(ll+1)*Nt,4*(ll-1)*Nt+1:4*ll*Nt) = b;
end

% The final two columns
b_l=[-W(:,:,1), ZZ;
    V(:,:,1), ZZ;
    repmat(ZZ,2*(N_l-1),2);
    ZZ, W(:,:,N_l);
    ZZ, V(:,:,N_l)];

B(:,4*N_l*Nt+1:4*(N_l+1)*Nt) = b_l; % The BIG matrix

% Incident field matrix
Ein = inc;
Hin = V(:,:,1)*Ein;
Fin = [Ein;Hin;repmat(zz,2*N_l,1)]; % Fields in each layer at incidence

C = B\Fin; % Column matrix with wave amplitudes for each layer

% Reflected and transmitted fourier modes
Cc = reshape(C,4*Nt,N_l+1); % Rearrange
Cc(2*Nt+1:4*Nt,N_l) = zz; % Forcing incidence from transmission layer to 0

R_xy = Cc(1:2*Nt,N_l+1);
T_xy = Cc(2*Nt+1:4*Nt,N_l+1);

%% Efficiency

R_z = (k_x_mn./k_1_zmn).*R_xy(1:Nt) + (k_y_mn./k_1_zmn).*R_xy(Nt+1:2*Nt);
T_z = -(k_x_mn./k_3_zmn).*T_xy(1:Nt) - (k_y_mn./k_3_zmn).*T_xy(Nt+1:2*Nt);

R_xyz = [R_xy;R_z]; T_xyz = [T_xy;T_z];

DEri = real([(k_1_zmn/k_1_zmn(m0));(k_1_zmn/k_1_zmn(m0));(k_1_zmn/k_1_zmn(m0))]).*(R_xyz.*conj(R_xyz));
DEti = real([(k_3_zmn/k_1_zmn(m0));(k_3_zmn/k_1_zmn(m0));(k_3_zmn/k_1_zmn(m0))]).*(T_xyz.*conj(T_xyz));

DEr = sum(DEri(:));
DEt = sum(DEti(:));
Abs = 1-DEr-DEt
%% Fields

% Layer hieght array
h_l = [h_z,repmat([D_d,D_m,D_d],1,NN_z)];
% Layer permittivity array
eps_l = [n_1^2;e_p;repmat(conj([e_d;e_m;e_d]),NN_z,1);n_3^2];

% Space parameters
rz = 2; % Resolution in 'z'
zStart = -20; % Start plotting from this z coordinate
zEnd = sum(h_l) - zStart; % Stop plotting at this z coordinate
z = zStart:rz:zEnd;

rx = 5; % Resolution in 'x'
x = 1:rx:D_x-rx+1; % Need to compute just a single period

ry = 2; % Resolution in 'y'
y = D_y/2;

zi = 1; % Starting z index
ll = 1; % Starting layer

% z coordinate for the end of each layer
z_l = cumsum([0,h_l,Inf]);

% Preallocating field magnitude matrices
eps_z = zeros(length(x),length(y),length(z));
E_x = zeros(length(x),length(y),length(z));
E_y = zeros(length(x),length(y),length(z));
E_z = zeros(length(x),length(y),length(z));
H_x = zeros(length(x),length(y),length(z));
H_y = zeros(length(x),length(y),length(z));
H_z = zeros(length(x),length(y),length(z));
E_v = zeros(length(x),length(y),length(z),3);
H_v = zeros(length(x),length(y),length(z),3);
S_v = zeros(length(x),length(y),length(z),3);

while zi<=length(z)
    % Propagation in z
    if z(zi)<z_l(ll)
        p_z = exp(-k_0*qm(:,ll)*(z(zi)-z_l(ll-(ll~=1))));
        Phi_z = diag([p_z;X(:,:,ll)*p_z.^(-1)]);
        
        M = [W(:,:,ll), W(:,:,ll);
             EPS(:,:,ll)\[-K_y,K_x]*V(:,:,ll), -EPS(:,:,ll)\[-K_y,K_x]*V(:,:,ll);
             V(:,:,ll), -V(:,:,ll);
             [-K_y,K_x]*W(:,:,ll), [-K_y,K_x]*W(:,:,ll)];
             
        Psi_z = M*Phi_z*Cc(:,ll); % Field amplitude coefficients at (*,*,z) in the incidence layer
        
        for yi=1:length(y)
            p_y = exp(-k_0*diag(K_y)*y(yi)); % Propagation in x
            Phi_y = diag(repmat(p_y,6,1));
    
            Psi_yz = Phi_y*Psi_z; % Field amplitude coefficients at (*,y,z)
            for xi=1:length(x)
                p_x = exp(-k_0*diag(K_x)*x(xi)); % Propagation in x
                Phi_x = diag(repmat(p_x,6,1));
    
                Psi_xyz = Phi_x*Psi_yz; % Field amplitude coefficients at (x,y,z)
    
                % Separating the field components from \Psi(x,z)
                E_x(xi,yi,zi) = sum(Psi_xyz(1:Nt));
                E_y(xi,yi,zi) = sum(Psi_xyz(Nt+1:2*Nt));
                E_z(xi,yi,zi) = sum(Psi_xyz(2*Nt+1:3*Nt));
                H_x(xi,yi,zi) = sum(Psi_xyz(3*Nt+1:4*Nt));
                H_y(xi,yi,zi) = sum(Psi_xyz(4*Nt+1:5*Nt));
                H_z(xi,yi,zi) = sum(Psi_xyz(5*Nt+1:6*Nt));
                
                E_v(xi,yi,zi,:) = [E_x(xi,yi,zi),E_y(xi,yi,zi),E_z(xi,yi,zi)];
                H_v(xi,yi,zi,:) = [H_x(xi,yi,zi),H_y(xi,yi,zi),H_z(xi,yi,zi)];
                S_v(xi,yi,zi,:) = 0.5*real(cross(E_v(xi,yi,zi,:),conj(H_v(xi,yi,zi,:))));
                if(z(zi)<D_z && z(zi)>0)
                    if(indexProfLog(round(x(xi)),round(y(yi)),ll-1))
                        eps_z(xi,yi,zi) = eps_l(ll);
                    else
                        eps_z(xi,yi,zi) = e_f;
                    end
                else
                    eps_z(xi,yi,zi) = eps_l(ll);
                end
            end
        end
        zi = zi + 1;
    else
        ll = ll + 1;
    end
end

%% Losses

absE2 = [];
imagEps = [];
absE2(:,:) = abs(E_x(:,1,:)).^2 + abs(E_y(:,1,:)).^2 + abs(E_z(:,1,:)).^2;
imagEps(:,:) = abs(imag(eps_z(:,1,:)));
ohmL = imagEps.*absE2;

%% Plot
figure('Units','inches','Position',[3 2 12 5]);
h = pcolor(z,x,ohmL);
set(h,'linestyle','none');
set(gca,'fontsize',20);
colormap('jet');
bb = colorbar;
bb.Label.String = 'Im(\epsilon)|E|^2';
bb.Label.Rotation = 270;
clims = get(gca,'clim');
bb.Label.Position = [3.5,(clims(1)+clims(2))/2,0];
xlabel('z (nm)');
ylabel('x (nm)');

figure('Units','inches','Position',[3 2 12 5]);
h = pcolor(z,x,absE2);
set(h,'linestyle','none');
set(gca,'fontsize',20);
colormap('jet');
bb = colorbar;
bb.Label.String = '|E|^2';
bb.Label.Rotation = 270;
clims = get(gca,'clim');
bb.Label.Position = [3.75,(clims(1)+clims(2))/2,0];
xlabel('z (nm)');
ylabel('x (nm)');