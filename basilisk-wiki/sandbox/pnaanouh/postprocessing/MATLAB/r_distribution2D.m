function [gav, X,Y] = r_distribution2D(log,Ls,D,rmax,dr)
%   Computes the 2-D pair distribution function of the droplets
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

%sets the edges of the histigram bins in the radial direction starting at D/100 to avoid a peak at zero
Edges=D/100:dr:rmax; 
%sets the edges of the histigram bins in the angular direction
Edges_a=linspace(-pi./2,pi./2,251); 
thetas=zeros(length(Edges_a)-1,1);
radii=zeros(length(Edges)-1,1);
for j=1:length(radii)
        radii(j)=(Edges(j+1)+Edges(j))/2;
end
for j=1:length(thetas)
        thetas(j)=(Edges_a(j+1)+Edges_a(j))/2;
end
% initialize the pair distrubution function to zero
gav=zeros(length(Edges)-1,length(Edges_a)-1); 

%start averaging at t=200 to avoid transients
for i=2000:size(Drops_sorted,3) 
    mirrorx=zeros(size(Drops_sorted,1),3,3);
    mirrory=zeros(size(Drops_sorted,1),3,3);
    %symetries due to periodic boundaries
    for m=-1:1
        for n=-1:1
            mirrorx(:,n+2,m+2)= Drops_sorted(:,3,i)+n.*Ls;
            mirrory(:,n+2,m+2)= Drops_sorted(:,4,i)+m.*Ls;
        end
    end
    clear n m
    % initialize the instantaneous pair distibution function to zero 
    g=zeros(length(Edges)-1,length(Edges_a)-1); 
    % for every drop
    for j=1:size(Drops_sorted,1) 
        % distance between it and every other drop with symmetries
        Distances= sqrt( (Drops_sorted(j,3,i)-mirrorx).^2+(Drops_sorted(j,4,i)-mirrory).^2);
        % Angles between it and every other drop
        Angles= asin((Drops_sorted(j,4,i)-mirrory)./Distances); 
        %pair distribution aroud the droplet
        result=histcounts2(Distances,Angles,Edges,Edges_a); 
        % instantaneous pair distribution is the average of the pair distibutions of all droplets divided by the number density of droplets
        g=g+result./size(Drops_sorted,1)./(size(Drops_sorted,1)./(Ls.^2));  
    end
    clear result j
    % scaled to take into account the area of the bins
    for j=1:size(g,1) 
        g(j,:)=g(j,:)./((pi*(Edges(j+1).^2-Edges(j).^2))./250);
    end
    % the pair distibution function is the mean of the instantaneous pair distribution functions
    gav=gav+g./(size(Drops_sorted,3)-2000);
end
%plot of the pair distibution function
    [Theta,R]=meshgrid(thetas,radii);
    X=R.*cos(Theta);
    Y=R.*sin(Theta);
    figure('Name','Surface','Position',[0 0 600 920])
    surf(X,Y,gav,'EdgeColor','none','FaceColor','flat'); 
    view(2);  
    colorbar
    colormap jet
end 