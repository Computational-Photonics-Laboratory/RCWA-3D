%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates the absorption spectrum of Cu-based metal-dielectric
% impedance matched metamaterials with a moth-eye surface, as shown in 
% Fig. 4 of R. Contractor, G. D’Aguanno, C. Menyuk, “Ultra-broadband, polarization-
% independent, wide-angle absorption in impedance-matched metamaterials with anti-
% reflective moth-eye surfaces”, Optics Express.

% The RCWA implementation is based on:
% M. G. Moharam and T. K. Gaylord, “Rigorous coupled-wave analysis of planar-grating 
% diffraction,” J. Opt. Soc. Am. 71, 811-818 (1995).
% L. Li, “New formulation of the Fourier modal method for crossed surface-relief 
% gratings,” J. Opt. Soc. Am. A 14, 2758-2767 (1997).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Cluster Computing
% Breakup the wavelength array to send to Maya

%{
indx = str2double(getenv('SLURM_ARRAY_TASK_ID')); 
y = 700.1:10:6250.1; 
if (indx==278)
y = y(length(y));
else
y = y((indx*2 + 1):(indx+1)*2);
end
%}

%% Inputs and Structure Parameters

y = 700.1:10:6250.1; % Range of wavelength in nm

% Angles
th = 0:0.25:89.75; % Incidence Angle (deg)
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
NN_z = 6; % Number of uncarved unit cells (after the moth-eyes)
% NN_z = 10; for uncarved % NN_z = 0; for metal moth-eyes on metal substrate

% Moth-eye Structures
N_m = 4; % = 0; for uncarved % Number of elementary cells carved into pillars
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
N_x = 9; %=1; for uncarved % Number of Fourier modes in x direction 
N_y = 9; %=1; for uncarved % Number of Fourier modes in y direction
Nt = N_x*N_y;

% Preallocations
ky = 2*pi./y; % Wavenumbers
I = eye(Nt); II = eye(2*Nt); % Identity matrices
ZZ = zeros(2*Nt);
EPS = zeros(Nt,Nt,N_z);
EPSyx = zeros(Nt,Nt,N_z);
EPSxy = zeros(Nt,Nt,N_z);
W = zeros(2*Nt,2*Nt,N_l);
V = zeros(2*Nt,2*Nt,N_l);
X = zeros(2*Nt,2*Nt,N_l);
DEr = zeros(length(y),length(th),length(psi));
DEt = zeros(length(y),length(th),length(psi));

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

%% Begin loop in wavelength

for yy=1:length(y)
	uy = y(yy)/1000; % Wavelength in micron

	%% Structure Analysis

	% Permittivity and refractive indices
	n_1 = 1; % Refractive index of incidence medium
	n_3 = 1; % Refractive index of transmission medium
        % n_3 = conj(epsMetal(yy))^2; for metal as substrate/transmission medium
	e_m = epsMetal(yy); % Permittivty of metal for current wavelength
	e_d = eps_d(yy); % Permittivity of dielectric layer (Ta2O5)
        % e_d = e_m; for all metal moth-eyes
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
				EPS(pp+N_x*(qq-1),:,nn) = reshape(EE,1,[]);
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

	%% Begin loop in polarization and incidence angle

	for ss=1:length(psi)
		for tt=1:length(th)
			% Incident Field
			u_x = cosd(psi(ss))*cosd(phi)*cosd(th(tt)) - sind(psi(ss))*sind(phi);
			u_y = cosd(psi(ss))*sind(phi)*cosd(th(tt)) + sind(psi(ss))*cosd(phi);

			% Wavenumber Matrices
			k_0 = ky(yy);

			k_xi = k_0*n_1*sind(th(tt))*cosd(phi) + 2*pi*mx/D_x;
			k_yi = k_0*n_1*sind(th(tt))*sind(phi) + 2*pi*my/D_y;
			k_x_mn = reshape(repmat(k_xi,1,N_y),[],1);
			k_y_mn = reshape(repmat(k_yi,N_x,1),[],1);

			k_1_zmn = conj(((n_1*k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));
			k_d_zmn = conj((e_d*(k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));
			k_m_zmn = conj((e_m*(k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));
			k_3_zmn = conj(((n_3*k_0)^2*ones(Nt,1) - k_x_mn.^2 - k_y_mn.^2).^(1/2));

			K_x = 1i*diag(k_x_mn/k_0);
			K_y = 1i*diag(k_y_mn/k_0);

			%% Eigenmatrices

			% Eigenvalue formulation for the grating region
			for nn=1:N_z
				F11 = -K_x*(EPS(:,:,nn)\K_y);
				F12 = I + K_x*(EPS(:,:,nn)\K_x);
				F21 = -I - K_y*(EPS(:,:,nn)\K_y);
				F22 = K_y*(EPS(:,:,nn)\K_x);
				F = [F11 F12;F21 F22];
				
				G11 = -K_x*K_y;
				G12 = EPSyx(:,:,nn) + K_x^2;
				G21 = -EPSxy(:,:,nn) - K_y^2;
				G22 = K_x*K_y;
				G = [G11 G12;G21 G22];
				
				[W(:,:,nn+1), Q] = eig(F*G);
				Q = (sqrt(Q));
				qm = diag(Q);
				V(:,:,nn+1) = -G*W(:,:,nn+1)/Q;
				X(:,:,nn+1) = diag(exp(-k_0*qm*h_z(nn)));
			end

			% Eigenmatrices for reflection and transmission layers
			W(:,:,1) = II;
			W(:,:,N_l) = II;

			VIxx = (1i/k_0)*diag((k_x_mn.*k_y_mn)./k_1_zmn);
			VIxy = (1i/k_0)*diag((k_y_mn.^2 + k_1_zmn.^2)./k_1_zmn);
			VIyx = -(1i/k_0)*diag((k_x_mn.^2 + k_1_zmn.^2)./k_1_zmn);
			VIyy = -VIxx;
			V(:,:,1) = [VIxx VIxy;VIyx VIyy];
			X(:,:,1) = II;

			VIIxx = (1i/k_0)*diag((k_x_mn.*k_y_mn)./k_3_zmn);
			VIIxy = (1i/k_0)*diag((k_y_mn.^2 + k_3_zmn.^2)./k_3_zmn);
			VIIyx = -(1i/k_0)*diag((k_x_mn.^2 + k_3_zmn.^2)./k_3_zmn);
			VIIyy = -VIIxx;
			V(:,:,N_l) = [VIIxx VIIxy;VIIyx VIIyy];
			X(:,:,N_l) = II;

			% Eigenmatrices for the uncarved layers
			for nn=1:3*NN_z
				if (mod(nn+1,3))==0
					k_n_zmn = k_m_zmn;
					D_n = D_m;
				else
					k_n_zmn = k_d_zmn;
					D_n = D_d;
				end
				W(:,:,1+N_z+nn) = II;
				Vnxx = (1i/k_0)*diag((k_x_mn.*k_y_mn)./k_n_zmn);
				Vnxy = (1i/k_0)*diag((k_y_mn.^2 + k_n_zmn.^2)./k_n_zmn);
				Vnyx = -(1i/k_0)*diag((k_x_mn.^2 + k_n_zmn.^2)./k_n_zmn);
				Vnyy = -Vnxx;
				V(:,:,1+N_z+nn) = [Vnxx Vnxy;Vnyx Vnyy];
				X(:,:,1+N_z+nn) = diag(exp(-1i*[k_n_zmn;k_n_zmn]*D_n));
			end

			%% S matrix calculation

			% Preallocating
			Suu = II;
			Sud = ZZ;
			Sdu = ZZ;
			Sdd = II;

			for ll=1:N_l-1
				S = [W(:,:,ll+1) -W(:,:,ll);V(:,:,ll+1) V(:,:,ll)]\[W(:,:,ll)*X(:,:,ll) -W(:,:,ll+1)*X(:,:,ll+1);V(:,:,ll)*X(:,:,ll) V(:,:,ll+1)*X(:,:,ll+1)];
				
				% Cascade the next interface to the current S matrix
				Auu = Suu;
				Aud = Sud;
				Adu = Sdu;
				Add = Sdd;
				Buu = S(1:2*Nt,1:2*Nt);
				Bud = S(1:2*Nt,2*Nt+1:4*Nt);
				Bdu = S(2*Nt+1:4*Nt,1:2*Nt);
				Bdd = S(2*Nt+1:4*Nt,2*Nt+1:4*Nt);
				
				Suu = Buu*((II-Aud*Bdu)\Auu);
				Sud = Bud+(Buu*Aud*((II-Bdu*Aud)\Bdd));
				Sdu = Adu+(Add*Bdu*((II-Aud*Bdu)\Auu));
				Sdd = Add*((II-Bdu*Aud)\Bdd);
			end

			%% Efficiency

			R_xy = u_x*Sdu(:,m0)+u_y*Sdu(:,Nt+m0);
			T_xy = u_x*Suu(:,m0)+u_y*Suu(:,Nt+m0);

			R_z = (k_x_mn./k_1_zmn).*R_xy(1:Nt) + (k_y_mn./k_1_zmn).*R_xy(Nt+1:2*Nt);
			T_z = -(k_x_mn./k_3_zmn).*T_xy(1:Nt) - (k_y_mn./k_3_zmn).*T_xy(Nt+1:2*Nt);

			R_xyz = [R_xy;R_z]; T_xyz = [T_xy;T_z];

			DEri = real([(k_1_zmn/k_1_zmn(m0));(k_1_zmn/k_1_zmn(m0));(k_1_zmn/k_1_zmn(m0))]).*(R_xyz.*conj(R_xyz));
			DEti = real([(k_3_zmn/k_1_zmn(m0));(k_3_zmn/k_1_zmn(m0));(k_3_zmn/k_1_zmn(m0))]).*(T_xyz.*conj(T_xyz));

			DEr(yy,tt,ss) = sum(DEri(:));
			DEt(yy,tt,ss) = sum(DEti(:));
			[num2str(yy),'/',num2str(length(y)),', ',num2str(tt),'/',num2str(length(th))] % Display progress
		end
	end
end

%% Plot

% Angular Spectrum
pp = 1; % Select the polarization for the plot in terms of psi(pp)
figure('Units','inches','Position',[5 2 6 5]);
h = pcolor(y/1000,th,(1-DEr(:,:,pp)-DEt(:,:,pp))');
set(h,'linestyle','none');
set(gca,'fontsize',20);
colormap('jet');
bb = colorbar;
bb.Label.String = 'Absorption';
bb.Label.Rotation = 270;
clims = get(gca,'clim');
bb.Label.Position = [3.75,(clims(1)+clims(2))/2,0];
caxis([0,1]);
xlabel('\lambda (\mum)');
ylabel('\theta (^o)');

% Line Plot
tt = 1; % [1,121,241]; % Incidence angles to be plotted in terms of th(tt)
figure('Units','inches','Position',[5 2 6 5]);
plot(y/1000,1-DEr(:,tt(1),pp)-DEt(:,tt(1),pp),'linewidth',2);
hold on
for aa=2:length(tt)
    plot(y/1000,1-DEr(:,tt(aa),pp)-DEt(:,tt(aa),pp),'-.r','linewidth',2);
end
xlim([min(y),max(y)]/1000);
xlabel('\lambda (\mum)');
ylabel('Absorption');
title(['\psi = ',num2str(psi(pp)),'\circ']);
set(gca,'fontsize',20);
ylim([0,1]);
labels = cell(size(tt));
for aa=1:length(tt)
    labels(aa) = {['\theta = ',num2str(th(tt(aa))),'\circ']};
end
legend(labels);

%% Save data on Maya
%{
myfile = ['W090','-',int2str(indx)];
save(myfile,'y','th','psi','DEr','DEt');
%}