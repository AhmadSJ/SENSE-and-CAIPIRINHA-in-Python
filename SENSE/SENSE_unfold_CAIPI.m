
function [recon_im,gmaps,unmix,recon_im_temp] = SENSE_unfold_CAIPI(target_data, coilsensmaps,acc_factors,CAIPI_shifts,lambda,recvr_noise_m,choose_outpt)
%Function that combines a dataset recorded by multiple coils to one
%dataset. 
%
%Input: 
%   target_data     Data recorded by multiple coils (x,y,z,coils) in image 
%                   domain. If there are dynamics present: (x,y,
%                   z,coils,dyn). 
%   coilsensmaps    Coil sensitivity maps (x,y,z,coils)
%   acc_factors     SENSE acceleration factors (Acc_x,Acc_y,Acc_z)
%   caipi_shift     Shiftdimension,number_of voxels shifted, ie 'x1'
%Optional input:
%   lambda          Regularization factor. Default = 0 (off)
%   recvr_noise_m   Receiver noise matrix. Defined as the correlation
%                   of noise in the receiver channels (ncoils x ncoils).
%                   If not given, the identity matrix is used. 
%   choose_outpt    This option gives you a way to calculate only specific
%                   output (faster). Other output variables will be empty. 
%                       0) Calculate all outputs (default)
%                       1) Only calc the reconstructed image
%                       2) Only calc g-map
%                       3) Only calc the unmixing coefficients
%Output: 
%   recon_im        A reconstructed coil combined dataset. 
%Optional output: 
%   gmaps           g-factor maps. 
%   unmix           Unmixing coefficients
%
%                                                   02-10-2014  Arj@H

% alias_A=zeros(size(target_data,1)*size(target_data,2),3);
%  xa=1;

%Number of coils
ncoils=size(target_data,4); 

%Initial Input
if nargin<3
    acc_factors(1)=size(coilsensmaps,1)/size(target_data,1);
    acc_factors(2)=size(coilsensmaps,2)/size(target_data,2);
    acc_factors(3)=size(coilsensmaps,3)/size(target_data,3);end
if nargin<4
    CAIPI_shifts='x1';end
if nargin<5
    lambda=0; end  %Regularization
if nargin<6
    recvr_noise_m=[]; end 
if nargin<7
    choose_outpt=0;
end

%Assign matrixes
Rx=acc_factors(1);Ry=acc_factors(2);Rz=acc_factors(3);
sensmaps=coilsensmaps;
dim_xyz = [size(sensmaps,1),size(sensmaps,2),size(sensmaps,3)];

%Read CAIPIRINHA parameters
    %Read CAIPI code
    %[shiftdir,vshift] = readCAIPIcode(CAIPI_shifts);
    shiftdir = CAIPI_shifts(1);
    vshift = CAIPI_shifts(2);
    %Find opposing CAIPI direction
    oppdim_list=acc_factors;
    oppdim_list(shiftdir)=0; %remove CAIPI direction
    oppdim_list((dim_xyz==1))=0; %remove direction with length 1 (only one slice). 
    %Remove SENSE 1 directions, if opp direction is not yet found
    if sum(oppdim_list)>1
        oppdim_list((oppdim_list==1))=0; %remove direction with SENSE 1. 
    end
    shiftdir_opp = find(oppdim_list>=1);
    if length(shiftdir_opp)>1 || isempty(shiftdir_opp);
        disp(' ')
        disp('Cannot determine the opposing CAIPI direction...')
        question_prompt = 'Give the dimnumber of the dimension opposing CAIPI-shift direction: ';
        dimstr = input(question_prompt,'s');
        shiftdir_opp = str2num(dimstr);
    end
    disp(['     Opposite CAIPI direction is: ',num2str(shiftdir_opp)])

    %Check if the amount of CAIPI shift is lower then the SENSE factor in CAIPI direction. 
    %Otherwise perform correction. 
    if vshift >= acc_factors(shiftdir) %Shift is larger than SENSE factor
        disp('     The CAIPIshift is larger or equal to the SENSE-factor.')
        mfactordiff = floor(vshift/acc_factors(shiftdir));  %How many times is it larger?
        vshift = vshift - acc_factors(shiftdir)*mfactordiff;
        disp_SENSE_code=sort([shiftdir,shiftdir_opp]);
        disp(['     CAIPI shift is corrected from ',num2str(acc_factors(disp_SENSE_code(1))),'x',num2str(acc_factors(disp_SENSE_code(2))),'_',CAIPI_shifts,' to: ',num2str(acc_factors(disp_SENSE_code(1))),'x',num2str(acc_factors(disp_SENSE_code(2))),'_',CAIPI_shifts(1),num2str(vshift),' (same pattern)'])
    end

%Noise matrix
if isempty(recvr_noise_m); 
    W=eye(ncoils); %Identity matrix
else
    W=recvr_noise_m; %Noise matrix is given as input
end

%Size data
N=size(target_data);
FOVx_us=size(target_data,1);
FOVy_us=size(target_data,2);
FOVz_us=size(target_data,3);
FOV_us_opp = size(target_data,shiftdir_opp);
sensmapFOVx=size(sensmaps,1);
sensmapFOVy=size(sensmaps,2);
sensmapFOVz=size(sensmaps,3);

%Acceleration factor check, round
if Rx*FOVx_us ~= sensmapFOVx
    Rnew=sensmapFOVx/FOVx_us;
    disp(['Acceleration factor corrected from Rx ',num2str(Rx(1)), ' to: ',num2str(Rnew)])
    Rx=Rnew;
end
if Ry*FOVy_us ~= sensmapFOVy
    Rnew=sensmapFOVy/FOVy_us;
    disp(['Acceleration factor corrected from Ry ',num2str(Ry(1)), ' to: ',num2str(Rnew)])
    Ry=Rnew;
end
if Rz*FOVz_us ~= sensmapFOVz
    Rnew=sensmapFOVz/FOVz_us;
    disp(['Acceleration factor corrected from Rz ',num2str(Rz(1)), ' to: ',num2str(Rnew)])
    Rz=Rnew;
end

%Total acceleration
Rtot = Rx*Ry*Rz;

%Dimensions for reconstruction
recon_dims=[N(1)*Rx,N(2)*Ry,size(target_data,3)*Rz,N(5:end)];

% SENSE unfolding
%p= (S^H * W^-1 * S)^-1 * S^H * W^-1 * m
recon_im_temp=[];

%Initial
recon_im=zeros(recon_dims);
gmaps=zeros(recon_dims);
unmix=[];

%Keep track of reconstruction time
h = waitbar(0,'Initializing waitbar...');

%Loop over every pixel
for dynamic=1:size(target_data,5)    %Dynamic
    for z=1:size(target_data,3)       %Slice
        for x=1:size(target_data,1)    %X   
            for y=1:size(target_data,2) %Y      

                %Display the selected pixel for reconstruction
                disp(['x=',num2str(x),' y=',num2str(y),' z=',num2str(z)])

                %Select folded pixels
                xt=[]; yt=[]; zt=[];
                for n=1:ceil(Rx) %X
                    xt(n)=x + FOVx_us*(n-1);
                end
                for m=1:ceil(Ry) %Y
                    yt(m)=y + FOVy_us*(m-1);
                end
                for o=1:ceil(Rz) %Z
                    zt(o)=z + FOVz_us*(o-1);
                end

                %Center target image, due to smaller field of view. 
                %Offset of target image in respect to the sensmap
                pos_xshift=sensmapFOVx/2 - sensmapFOVx/(2*Rx);
                pos_yshift=sensmapFOVy/2 - sensmapFOVy/(2*Ry);
                pos_zshift=sensmapFOVz/2 - sensmapFOVz/(2*Rz);

                %Save pixel indices in dictionary (combine x and y coordinates)
                n=0; r=cell(1,length(xt)*length(yt)*length(zt)); %rm=[];
                for xp=1:length(xt)
                    for yp=1:length(yt)
                        for zp=1:length(zt)
                            n=n+1;
                            r{n}=[round(xt(xp)+pos_xshift),round(yt(yp)+pos_yshift),round(zt(zp)+pos_zshift)];


                            %Correct for CAIPI shift

                            %N voxels to shift = FOV_sensmap * (kspace_shift/R_total)
                            %Old: sens_map_shift =  round(size(sensmaps,shiftdir_opp) * vshift * (1/Rtot));
                            sens_map_shift =  FOV_us_opp - round(size(sensmaps,shiftdir_opp) * vshift * (1/Rtot));

                            %Give CAIPI shift to position matrix (r)
                           if vshift>0
                            if  shiftdir==1     %Shift dir is X
                                r{n}(shiftdir_opp)=r{n}(shiftdir_opp)+sens_map_shift*(xp-1);
                            elseif shiftdir==2  %Shift in Y
                                r{n}(shiftdir_opp)=r{n}(shiftdir_opp)+sens_map_shift*(yp-1);
                            elseif shiftdir==3  %Shift in Z
                                r{n}(shiftdir_opp)=r{n}(shiftdir_opp)+sens_map_shift*(zp-1);
                            end
                           end

                            %Check if the aliased pixel is more than 1x
                            %outside FOV. If yes, determine the multiplication factor. 
                            mfx=1; mfy=1; mfz=1; %Initial value 
                            if r{n}(1)>ceil(Rx)*FOVx_us
                                mfx = floor((r{n}(1)-1)/(ceil(Rx)*FOVx_us))+1;end
                            if r{n}(2)>ceil(Ry)*FOVy_us
                                mfy = floor((r{n}(2)-1)/(ceil(Ry)*FOVy_us))+1;end
                            if r{n}(3)>ceil(Rz)*FOVz_us
                                mfz = floor((r{n}(3)-1)/(ceil(Rz)*FOVz_us))+1;
                            end

                            %Correction for decimal acceleration factors (Rx=1.5)
                            %Remove index which is: FOV sensmap < index < rounded FOV
                            if r{n}(1)>(size(sensmaps,1)+ceil(Rx)*FOVx_us*(mfx-1)) && r{n}(1)<=(ceil(Rx)*FOVx_us*mfx)
                                r(n)=[]; n=n-1; end %X
                            if r{n}(2)>(size(sensmaps,2)+ceil(Ry)*FOVy_us*(mfy-1)) && r{n}(2)<=(ceil(Ry)*FOVy_us*mfy)
                                r(n)=[]; n=n-1; end %Y
                            if r{n}(3)>(size(sensmaps,3)+ceil(Rz)*FOVz_us*(mfz-1)) && r{n}(3)<=(ceil(Rz)*FOVz_us*mfz)
                                r(n)=[]; n=n-1;     %Z
                            end

                            %Fold back in image, if position is outside the
                            %field of view sensitivitymap (also possible by using the mod (modulus) function). 
                            if r{n}(1)>size(sensmaps,1) %X              %Correction when more than 1 "large sensmapsize" (=R*FOV_us) outside of FOV
                                r{n}(1)=r{n}(1)-ceil(Rx)*FOVx_us        * (mfx-1);end
                               %r{n}(1)=r{n}(1)-ceil(Rx)*FOVx_us ;end 
                            if r{n}(2)>size(sensmaps,2) %Y
                                r{n}(2)=r{n}(2)-ceil(Ry)*FOVy_us        * (mfy-1);end
                               %r{n}(2)=r{n}(2)-ceil(Ry)*FOVy_us ;end 
                            if r{n}(3)>size(sensmaps,3) %Z
                                r{n}(3)=r{n}(3)-ceil(Rz)*FOVz_us        * (mfz-1);
                               %r{n}(3)=r{n}(3)-ceil(Rz)*FOVz_us;
                            end

                             %save coordinate of aliased points.
%        alias_A(xa,:)=r{n};
%        xa=xa+1;  

                            %Display the calculated aliased pixel positions
                            %disp([num2str(xp),' ',num2str(yp),' ',num2str(zp),'   ',num2str(r{n})])

                            %rm(x,y,n,:)=[xt(xp),yt(yp)];
                        end
                    end
                end

                %Build sensitivity matrix S
                S=zeros(ncoils,length(r));
                for coil=1:ncoils
                    for n=1:length(r)
                        S(coil,n)=sensmaps(r{n}(1),r{n}(2),r{n}(3),coil);
                    end
                end

                %Regularization R = lambda^2 * I
                R = lambda.^2 * eye(size(S,2),size(S,2));

                if choose_outpt==0 || choose_outpt==1 || choose_outpt==3
                    %Build unfolding matrix
                    %warning('off')
                    U= pinv(S' * inv(W) * S + R) * S' * inv(W);  % Build unfolding matrix U(np x nc)
                    %warning('on')
                end

                %-------------Extra, unmixingcoefficients & gmaps-----------
                if choose_outpt==0 || choose_outpt==3     
                    %Save unfolding matrix to obtain unmixing coeffiecients
                    for n=1:length(r) %number of folded pixels
                        unmix(r{n}(1),r{n}(2),r{n}(3),:,dynamic)=U(n,:);
                    end
                end
                if choose_outpt==0 || choose_outpt==2     
                    %Calculate gmaps
                    %Equation g = sqrt( diag( (S^H * W^-1 * S)^-1) .* diag(S^H * W^-1 * S)) Pruessmann et al. (1999)
                    gfactors = sqrt(diag(pinv(S' * inv(W) * S  + R)) .* diag((S' * inv(W) * S)));
                    %Place g-factors in the right position. 
                    for n=1:length(r) %number of folded pixels
                        gmaps(r{n}(1),r{n}(2),r{n}(3),:,dynamic) = gfactors(n);
                    end
                end
                %-----------------------------------------------------------

                if choose_outpt==0 || choose_outpt==1  

                    %Get measured (folded) signal from target scan
                    ms=squeeze(target_data(x,y,z,:,dynamic));  %<-- wrongly selected

                    %Reconstruction
                    p=U*ms;     %Calculate new pixel values

                    for n=1:length(p)
                        recon_im(r{n}(1),r{n}(2),r{n}(3),dynamic)=p(n);
                    end

                    %Videoframes, step-by-step, moving point
%                     for n=1:length(p)
%                         if r{n}(1)<recon_dims(1) && r{n}(2)<recon_dims(2)
%                             if y~=size(target_data,2)
%                                 recon_im(r{n}(1),r{n}(2)+1,z,dynamic)=50;
%                             end
%                         end
%                     end
%                     recon_im_temp = cat(3,recon_im_temp,recon_im);

                end

            end %y
        end %x

        %Update waitbar
        process_perc = round(z/size(target_data,3) *100);
        waitbar(process_perc/100,h,sprintf('%d%% of reconstruction',process_perc))

    end %z
end %dyn
% mean_g=mean(mean(mean(gmaps)));
% var_g=var(var(var(gmaps)));
% save ('prova', 'gmaps', 'mean_g','var_g');

%Close waitbar
delete(h)