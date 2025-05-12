function [t,t_bins,n_int_start,t_range,Dist_ref,vx,vy,vt,Dist_signal,t_cont,t_min,Dist_free,Dist_not_free,n_min,Dist] = drop_statistics(log,Ls,D,lDref,lt_range)
%Primary Post-proccessing function for the calculation of droplet
%statistics
%   Not all of the statistics that are computed are outputed so feel free
%   to either remove the parts you don't want or to add them as outputs 
%   
%   This function takes as an input a log file generated in basilisk ( see
%   sandbox/pnaanouh/rising-suspension.c), the domain size of the
%   simulation, the droplet diameter, the number of different contact
%   threshold distances to consider, and the number of bins for the
%   calculation of the histogram of the contact time pdf
%
%   This function has may outputs:
%       t: the time at every sample in the simulation
%       t_bins: the values of the bins for the histogram of the contact 
%           time pdf at different contact threshold distances
%       n_int_start: the total number of collisions counted at the start of
%           the collisions at different contact threshold distances
%       t_range: the bins of the histogram of the contact
%           time pdf
%       Dist_ref: the range of the contact threshold distance
%       vx: the instananeous horizontal velocity of the dispersed phase as 
%           a function of time
%       vy: the instananeous vertical velocity of the dispersed phase as a 
%           function of time
%       vt: the total volume of the dispersed phase as a function of time
%       Dist_signal: Distance between each droplet and its nearest
%           neigbhour as a function of time
%       t_cont: the durations of all the measured collisions for different
%           contact threshold distances
%       t_min: the durations of all the measured collisions for different
%           contact threshold distances with the assumption of binary
%           collisions
%       Dist_free: mean free path between collisions for different contact threshold
%           distances
%       Dist_not_free: mean distance travelled during collisions path for  
%           different contact threshold distances
%       n_min:the total number of collisions at different contact threshold
%           distances  with the assumption of binary collisions
%       Dist: the distance between every pair of droplets at every sample
%           time


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
%sort the drops into a consistent order
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

t=Drops_final(:,2,1);
%get the drop velocities

%Speed of each drop
%pvel=sqrt(Drops_sorted(:,5,:).^2+Drops_sorted(:,6,:).^2);
%pvel=permute(pvel,[1,3,2]);
%Average speed at every instant
%pvelt=sum(pvel.*Drops_sorted(:,7,:))./sum(Drops_sorted(:,7,:));
%pvelt=permute(pvelt,[3,2,1]);
%Average velocity at every instant
vx=sum(Drops_sorted(:,5,:).*Drops_sorted(:,7,:))./sum(Drops_sorted(:,7,:));
vx=permute(vx,[3,2,1]);

vy=sum(Drops_sorted(:,6,:).*Drops_sorted(:,7,:))./sum(Drops_sorted(:,7,:));
vy=permute(vy,[3,2,1]);

%Total volume in time
vt=sum(Drops_sorted(:,7,:));
vt=permute(vt,[3,2,1]);

%calculate the distance between all drops and thus the collision freq and
%contact time

%initialize the matrix to store the distances between all the drops at
%every timestep
Dist=zeros(size(Drops_sorted,1),size(Drops_sorted,1),size(Drops_sorted,3));
%initialize the matrix to store the distance travelled by each drop at
%every timestep
dx=zeros(size(Drops_sorted,1),size(Drops_sorted,3));
%initialize the matrix to store the distances between each drop and its
%nearest neighbour at every timestep
Dist_signal=zeros(size(Drops_sorted,1),size(Drops_sorted,3));

%set the range for the contact threshold distance between 0 and 4 diameters
%with lDref number of points
Dist_ref=linspace(0,4*D,lDref);

%initalize the vectors to store the total number of collisions n_int_start
%counted at the start of the collsions, n_int_end counted at the end of the
%collisions, and n_min calculated based on the minimum distance signal
n_int_start=zeros(size(Dist_ref));
n_int_end=zeros(size(Dist_ref));
n_min=zeros(size(Dist_ref));

%initalize the vectors to store the total contact time, the total distance
%traveled outside of collisions, and the total distance travalled during 
%collisions  
t_total=zeros(size(Dist_ref));
Dist_free=zeros(size(Dist_ref));
Dist_not_free=zeros(size(Dist_ref));

%initialize variables to store the durations of individual collisions 
t_cont=[]; % stores the collision durations
t_int=[]; % stores the time between the end of the collisions
t_min=[]; % stores the collision durations with the assumption of binary collisions

%iterate on the range of the contact reference distances for each pair of
%drops at every sample time
for n=1:length(Dist_ref)
    %initialize the matixies used to calculate the contact time 
    
    %stores the time spent in the current collision 
    t_c3=zeros(size(Drops_sorted,1),size(Drops_sorted,1));
    t_i3=zeros(size(Drops_sorted,1),size(Drops_sorted,1));
    t_m3=zeros(size(Drops_sorted,1));
    
    %stores the matrix that contains the duration of all the collisions it
    %is neccessary due to the change in the number of collisions with the
    %contact threshold distance 
    t_c2=struct;
    t_i2=struct;
    t_m2=struct;
    
    %stores the duration of all the collisions 
    t_c=[];
    t_i=[];
    t_m=[];
    
    %matix to store the state of each droplet pair( colliding 1, not
    %colliding 0) at the current sample and the previous sample
    ck=zeros(size(Drops_sorted,1),size(Drops_sorted,1),2);
    %matix to store the state of each droplet( colliding 1, not
    %colliding 0) at the current sample and the previous sample
    ck2=zeros(size(Drops_sorted,1),2);
    
    %loop on sample time
    for i=1:size(Drops_sorted,3) 
        % current sample 
        a=mod(i,2)+1; 
        %previous sample 
        b=mod(i-1,2)+1; 
        %initialize the distance to the domain size, this is unnecessary 
        Dist(:,:,i)=Ls;  
        % loop on each droplet
        for j=1:size(Drops_sorted,1)
            %the distance on the diagonal is set to the domain size, again unnecessary
            Dist(j,j,i)=Ls;   
            % matrix to store the distance travelled by the drop between the current sample 
            %and the previous sample with symmetries accross the periodic boundaries
            mv=zeros(3,3);
            %calculation of the distance travelled 
            if(i>1) 
                for l=-1:1
                    for m=-1:1
                        mv(l+2,m+2)=sqrt((Drops_sorted(j,3,i)-Drops_sorted(j,3,i-1)-Ls*l).^2+(Drops_sorted(j,4,i)-Drops_sorted(j,4,i-1)-Ls*m).^2);
                    end
                end
            end
            % the actual distance travelled is the minimum of all the symmetries
            dx(j,i)=min(mv,[],'all'); 
            
            % loop of every pair of drops starts from j+1 to get upper triangular matrix 
            for k=j+1:size(Drops_sorted,1)  
                % checks that both droplets actually exist can be relevant if there is breakup or coalescence due to how Drops_sorted is constructed
                if(Drops_sorted(k,7,i)~=0 && Drops_sorted(j,7,i)~=0)
                    % matrix to store the distance between the two drops plus the symmetries
                    DD=zeros(3,3);
                    %calculate the distance
                    for l=-1:1 
                        for m=-1:1
                            DD(l+2,m+2)=sqrt((Drops_sorted(k,3,i)-Drops_sorted(j,3,i)-Ls*l).^2+(Drops_sorted(k,4,i)-Drops_sorted(j,4,i)-Ls*m).^2);
                        end
                    end
                    % the actual distance is the minimum of the symmetries
                    Dist(j,k,i)=min(DD ,[],'All'); 
                    % adding the symmetric of the distance in the distance matrix
                    Dist(k,j,i)=Dist(j,k,i); 
                    % store if the droplets are colliding or not
                    ck(j,k,a)=Dist(j,k,i)<=Dist_ref(n); 
                    if(i>1)
                        %if the droplets were colliding and are no longer colliding 
                        if(ck(j,k,b)==1 &&ck(j,k,a)==0) 
                            % count the end of a collision 
                            n_int_end(n)=n_int_end(n)+1; 
                            % store the duration of the collision 
                            t_c=[t_c,t_c3(j,k)]; 
                            % initialize the collision duration to zero 
                            t_c3(j,k)=0; 
                            % store the time between collisions 
                            t_i=[t_i,t_i3(j,k)]; 
                            % initialize the time between collisions to zero
                            t_i3(j,k)=0; 
                            %if the droplets were not colliding and are now colliding
                        elseif(ck(j,k,b)==0 &&ck(j,k,a)==1) 
                            % count the start of a collision
                            n_int_start(n)=n_int_start(n)+1; 
                        end
                    end
                end
            end
            %computes the  minimum distance
            Dist_signal(j,i)=min(Dist(j,:,i)); 
            % checks if the droplet is colliding based on the assumption of binary collisions
            ck2(j,a)=Dist_signal(j,i)<=Dist_ref(n);
            % adds to the total distance travelled while colliding 
            Dist_free(n)=Dist_free(n)+(1-ck2(j,a))*dx(j,i); 
            % adds to the total distance travelled without colliding
            Dist_not_free(n)=Dist_not_free(n)+ck2(j,a)*dx(j,i);  
            if(i>1)
                % at the start of the collision 
                if(ck2(j,a)==1 && ck2(j,b)==0) 
                    % count the collision
                    n_min(n)=n_min(n)+1; 
                    % store the duration of the previous collision 
                    t_m=[t_m,t_m3(j)]; 
                    % initialize the collision duration to zero
                    t_m3(j)=0; 
                end
            end
        end
        %matrix computation of the collision durations, and the time
        %between collisions
        t_c3=t_c3+ck(:,:,a).*dt;  
        t_i3=t_i3+dt;
        t_m3=t_m3+ck2(:,a).*dt;
    end
    %the collision durations are stored in structures
    t_c2.n=t_c;
    t_i2.n=t_i;
    t_m2.n=t_m;
    % th structures are added to matrixies
    t_cont=[t_cont,t_c2];
    t_int=[t_int,t_i2];
    t_min=[t_min,t_m2];
end
%the mean free and mean not free paths are computed to be the total
%distance divided by the number of droplets divided by the number of
%sections 
Dist_free=Dist_free./size(Drops_sorted,1)./n_min;
Dist_not_free=Dist_not_free./size(Drops_sorted,1)./n_min;

clear i j k l m a b mv DD t_c t_c2 t_c3 t_i t_i2 t_i3 n ck ck2;

% calculation of the pdf of the contact time
% This can be replaced by a call to histogram or histcounts 

t_range=linspace(0,60,lt_range);
t_bins=zeros(length(t_range),length(t_cont));
t_bins2=zeros(length(t_range),length(t_int));
for i=1:length(t_range)-1
    for j=1:length(t_cont)
        for k=1:length(t_cont(j).n)
            if(t_cont(j).n(k)>t_range(i)&&t_cont(j).n(k)<=t_range(i+1))
                t_bins(i,j)=t_bins(i,j)+1;
            end
        end
    end
end
t_bins=t_bins./sum(t_bins);

for i=1:length(t_range)-1
    for j=1:length(t_int)
        for k=1:length(t_int(j).n)
            if(t_int(j).n(k)>t_range(i)&&t_int(j).n(k)<=t_range(i+1))
                t_bins2(i,j)=t_bins2(i,j)+1;
            end
        end
    end
end
t_bins2=t_bins2./sum(t_bins2);
clear i j k


end

