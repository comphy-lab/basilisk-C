function [sorted_state] = find_drops(curr_state,prev_state,dt,Ls)
%function that correlates the drops between two samples 
%
%   Takes as input the data for the current sample as well as the previous
%   sample, the time between samples and the domain size
%
%   Outputs a sorted version of the current state where each row in the
%   current sate and same row in the previous state correspond to the same
%   drop

%initialize the sorted state to zero and the distance between the two
%states to zero
sorted_state=zeros(size(curr_state));
Dist=ones(size(curr_state,1),size(curr_state,1))*Ls;

%loop over every drop in the current sate and ever drop 
for i=1:size(curr_state,1) in the previous state
    for j=1:size(curr_state,1)
        %predicted position of the droplet
        nx=prev_state(j,3)+prev_state(j,5)*dt; 
        ny=prev_state(j,4)+prev_state(j,6)*dt;
        % difference between the prediction and a measured position
        dx=abs(curr_state(i,3)-nx);
        dy=abs(curr_state(i,4)-ny); 
        % taking into account the periodic boundary conditions
        if(dx<Ls+prev_state(j,5)*dt && dx>Ls-prev_state(j,5)*dt)
            dx=abs(dx-Ls);
        end
        if(dy<Ls+prev_state(j,6)*dt && dy>Ls-prev_state(j,6)*dt)
            dy=abs(dy-Ls);
        end
        % Distance between the predicted position of drop j and the measured position of drop i
        Dist(i,j)=sqrt(dx^2+dy^2); 
    end
end
% vectors containing the minimums of the row and collumns of Dist
Minimums=min(Dist);
M=min(Dist,[],2);
% vector used to ensure that each drop is only correlated to one other drop
no_overlap=zeros(1,size(curr_state,1)); 
%loop on all the drops in the current and previous state
for i=1:size(curr_state,1)
    for j=1:size(curr_state,1)
    %The two drops that are the closest to each other are the same drop
        if(Dist(i,j)==Minimums(j)&&Minimums(j)==M(i)) 
            sorted_state(j,:)=curr_state(i,:);
            no_overlap(i)=1;
        end
    end
end
% If a drop isn't correlated to any other one then just assign it to any drop that is free. This only matters in the case of coalescence and breakup 
for i=1:size(curr_state,1) 
    for j=1:size(curr_state,1)
        if(sorted_state(j,7)==0 && no_overlap(i)==0)
            sorted_state(j,:)=curr_state(i,:);
            no_overlap(i)=1;
        end
    end
end
end

