function [gav,radii] = r_distribution(log,Ls,D,rmax,dr)
%   Computes the radial distribution function of the droplets
%
%   Takes as inputs the log file generated from basilisk, the domain size,
%   Droplet diameter, maximum radius of the radial distribution function,
%   and the thickness of one radial distribution shell

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

%sets the edges of the histigram bins starting at D/100 to avoid a peak at zero
Edges=D/100:dr:rmax; 
radii=zeros(length(Edges)-1,1); 
for j=1:length(radii)
        radii(j)=(Edges(j+1)+Edges(j))/2;
end
% initialize the radial distrubution function to zero
gav=zeros(1,length(Edges)-1); 

% for every sample aside form the first 2000 that are neglected due to the transients
for i=2000:size(Drops_sorted,3)  
    %take into account the symmetries accros the periodic boundaries
    mirrorx=zeros(size(Drops_sorted,1),3,3);  
    mirrory=zeros(size(Drops_sorted,1),3,3);
    for m=-1:1
        for n=-1:1
            mirrorx(:,n+2,m+2)= Drops_sorted(:,3,i)+n.*Ls;
            mirrory(:,n+2,m+2)= Drops_sorted(:,4,i)+m.*Ls;
        end
    end
    clear n m
    % initialize the instantaneous radial distibution function to zero 
    g=zeros(1,length(Edges)-1); 
    % for every drop 
    for j=1:size(Drops_sorted,1) 
        % distance between it and every other drop with symmetries
        Distances= sqrt((Drops_sorted(j,3,i)-mirrorx).^2+(Drops_sorted(j,4,i)-mirrory).^2); 
        %radial distribution aroud the droplet
        result=histcounts(Distances,Edges);
        % instantaneous radial distribution is the average of the radial distibutions of all droplets divided by the number density of droplets
        g=g+result./size(Drops_sorted,1)./(size(Drops_sorted,1)./(Ls.^2)); 
    end
    clear result j 
    % scaled to take into account the increasing area of the shells with the increase in the radius
    for j=1:size(g,2)
        g(j)=g(j)./(pi*(Edges(j+1).^2-Edges(j).^2));
    end
    % the radial distibution function is the mean of the instantaneous radial distribution functions
    gav=gav+g./(size(Drops_sorted,3)-2000); 
end
%plots the radial distibution function
figure
plot(radii,gav);
end

