%load('raw_data.mat');
%img = ifft2c(raw);
%s_ref = adaptive_est_sens(img);
%save('sref.mat','s_ref');
data = load('full_mat_for_adapt_est_sens.mat');
raw = data.raw;
img = ifft2c(raw);
s_ref=adaptive_est_sens(img);
function out = ifft2c(input)
    out = fftshift(ifft(ifftshift(input,1),[],1),1);
    out = fftshift(ifft(ifftshift(out,2),[],2),2);
end

function S = adaptive_est_sens(data)
    [Nx,Ny,Nz,Nc] = size(data);
    S = zeros(Nx,Ny,Nz,Nc);
    M = zeros(Nx,Ny,Nz);
    w = 5;
    for i = 1:Nx
        ii = max(i-w,1):min(i+w,Nx);
        for j = 1:Ny
            jj = max(j-w,1):min(j+w,Ny);
            for k = 1:Nz
                kk = max(k-w,1):min(k+w,Nz);
                kernel = reshape(data(ii,jj,kk,:),[],Nc);
                [V,D] = eigs(conj(kernel'*kernel),1);
                S(i,j,k,:) = V*exp(-1j*angle(V(1)));
                M(i,j,k) = sqrt(D);
            end
        end
    end
    S = S.*(M>0.1*max(abs(M(:))));
end