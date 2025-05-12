function [vxrms,vyrms,vxp,vyp,vxm,vym] = V_fluctuations(log,log_gl,Ls)
% Post-processing function for the droplet velocities
%
%   Takes as inputs the log file containing the individual droplet
%   properties as well as the one containing the global properties for the
%   entire computational domain as well as the domain size
%
%   Ouptuts the rms of the velocity fluctuations, the mean velocity, as well as the
%   fluctuation of the veloctiy from the mean as a function of time in both
%   the horizontal and vertical directions

%start off by turning the log file into something readable, seperate the
%different sample times and assotiate the samples to individual droplets
Drops=[];
%separate the drops based on the timestep
j=1;
for i=1:length(log)
    if(log(i,7)>0.5)
        Drops(j,:,int64((log(i,2)-log(1,2))*10+1))=log(i,:);
        if(i<length(log))
            if(log(i,1)==log(i+1,1))
                j=j+1;
            else
                j=1;
            end
        end
    end
end
clear j;
%sort the drops into a consitent order
Drops_sorted=zeros(size(Drops));
Drops_sorted(:,:,1)=Drops(:,:,1);
for i=2:size(Drops,3)
    dt=Drops(1,2,i)-Drops(1,2,i-1);
    Drops_sorted(:,:,i)=find_drops(Drops(:,:,i),Drops_sorted(:,:,i-1),dt,Ls);
end
clear i;
Drops_sorted=Drops_sorted(:,:,1:size(Drops_sorted,3)-1);

%invert dimensions 1 and 3 for added convenience

Drops_final=permute(Drops_sorted,[3,2,1]);

%t=Drops_final(:,2,1); 
vx=Drops_final(:,5,:)-log_gl(1:length(log_gl)-1,6);
vy=Drops_final(:,6,:)-log_gl(1:length(log_gl)-1,7);

%mean velocities
vxm=mean(vx,'all'); 
vym=mean(vy,'all');

% the fluctuations are the difference between the instantneous values and the mean values
vxp=(vx-vxm); 
vyp=(vy-vym);

%remove the transient phase
vxp=vxp(2000:size(vxp,1),:,:); 
vyp=vyp(2000:size(vyp,1),:,:);

% calculate the rms of the fluctuations
vxrms=sqrt(mean((vx-vxm).^2,'all')); 
vyrms=sqrt(mean((vy-vym).^2,'all'));



end

