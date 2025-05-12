/**
#A gnu-Octave / matlab code
To generate your own refinement attributes.

*/
% ax+by+cxy 
disp 2d
clear all
A=4*[1 0;1 -1;0 -1; -1 -1;-1 0;-1 1;0 1; 1 1]; 

x=A(:,1); 
y=A(:,2);
X(:,1)=x;
X(:,2)=y;
X(:,3)=(x.*y);
transpose(A)
C=inv(transpose(X)*X)*transpose(X);
C*192


%ax+by+cz+dxy+eyz+fxz+gxyz
disp 3dc
clear all
Locs=[ 1 0 0;1 -1 0;0 -1 0; -1 -1 0;-1 0 0;-1 1 0;0 1 0;1 1 0;
      0 0 -1; 1 0 -1;1 -1 -1;0 -1 -1; -1 -1 -1;-1 0 -1;-1 1 -1;0 1 -1;1 1 -1; 
     0 0 1; 1 0 1;1 -1 1;0 -1 1; -1 -1 1;-1 0 1;-1 1 1;0 1 1;1 1 1];
     
A=4*Locs; 
x=A(:,1); 
y=A(:,2);
z=A(:,3); 

X(:,1)=x;
X(:,2)=y;
X(:,3)=z;
X(:,4)=(x.*y);
X(:,5)=(z.*y);
X(:,6)=(x.*z);
X(:,7)=(x.*y.*z);
transpose(A)
C=inv(transpose(X)*X)*transpose(X);
C=C*8*9*8*8; 
[a,~]=size(C); 
for j=1:a
    ind=find(C(j,:)~=0);
    fprintf([num2str(round(abs(C(j,ind(1))))) '*']);
    if (length(unique(abs(C(j,ind))))==1)
        for g=1:length(ind); 
            if (C(j,ind(g))<0)
               fprintf('-');  
            end
            if (C(j,ind(g))>0 && g~=1)
                fprintf('+');
            end
        fprintf(['s[' num2str(Locs(ind(g),1)) ',' num2str(Locs(ind(g),2)) ',' num2str(Locs(ind(g),3)) ']']);
        end
        clear ind
        fprintf(';\n');
    else
        disp('Error! non unique weights');
    end
end

%ax+by+cz+dxy+eyz+fxz+gxyz+hxx+iyy+jzz+kyxx+lzxx+mxyy+nzyy+oxzz+pzzy+qxxyy+rxxzz+syyzz+txxyyzz
%20 weights!
disp ('3d third order accurate')
clear all
Locs=[0 0 0; 1 0 0;1 -1 0;0 -1 0; -1 -1 0;-1 0 0;-1 1 0;0 1 0;1 1 0;
      0 0 -1; 1 0 -1;1 -1 -1;0 -1 -1; -1 -1 -1;-1 0 -1;-1 1 -1;0 1 -1;1 1 -1; 
     0 0 1; 1 0 1;1 -1 1;0 -1 1; -1 -1 1;-1 0 1;-1 1 1;0 1 1;1 1 1;
     2 0 0; -2 0 0; 0 2 0; 0 -2 0; 0 0 2; 0 0 -2];
     
A=4*Locs; 
x=A(:,1); 
y=A(:,2);
z=A(:,3); 
% The terms
X(:,1)=x;                               %a
X(:,2)=y;                               %b
X(:,3)=z;                               %|  
X(:,4)=(x.*y);                         
X(:,5)=(x.*z);                          
X(:,6)=(y.*z);
X(:,7)=(x.*y.*z);
X(:,8)=x.*x;
X(:,9)=y.*y;
X(:,10)=z.*z;
X(:,11)=X(:,8).*y;
X(:,12)=X(:,8).*z;
X(:,13)=X(:,9).*x;
X(:,14)=X(:,9).*z;
X(:,15)=X(:,10).*x;
X(:,16)=X(:,10).*y;
X(:,17)=X(:,8).*X(:,9);
X(:,18)=X(:,8).*X(:,10);
X(:,19)=X(:,9).*X(:,10);                %|
X(:,20)=X(:,8).*X(:,9).*X(:,10);        %t
X(:,21)=ones(size(x));
S{:}=['a';'b';'c';'d';'e';'f';'g';'h';'i';'j';'k';'l';'m';'n';'o';'p';'q';'r';'s';'t';'u'];
C=inv(transpose(X)*X)*transpose(X);
 
[a,~]=size(C); 
C=round(C*100000,9);
for j=1:a
    fprintf('double %c = ',S{1}(j));
    ind=find(C(j,:)~=0);
    wei=(unique(abs(C(j,ind))));
    for h=1:length(wei); 
        ind = find(abs(C(j,:))==wei(h));
        fprintf('%.15g*(',(wei(h)));
        for g=1:length(ind); 
            if (C(j,ind(g))<0)
               fprintf('-');  
            end
            if (C(j,ind(g))>0 && g~=1)
                fprintf('+');
            end
        fprintf(['s[' num2str(Locs(ind(g),1)) ',' num2str(Locs(ind(g),2)) ',' num2str(Locs(ind(g),3)) ']']);
        end
        clear ind
        if h==length(wei)
            fprintf(');\n');
        else 
            fprintf(')+');
        end
    end
end



