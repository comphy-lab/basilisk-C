%{
# 2D differentiation matrices
Here we are going to study the 2D differential matrices for higher orders. We are going to go further with the code already developped for lower derivatives levels: [diffmat_2D.m](http://basilisk.fr/sandbox/easystab/diffmat_2D.m), [diffmat_2d_crossed_derivatives.m](http://basilisk.fr/sandbox/easystab/diffmat_2d_crossed_derivatives.m)

The first step is to define the nuber of grid points and to generate the derivative matrices. We also define a test fuction.
%}

clear all; clf

%%%% parameters and flags
Nx=20; % gridpoints in x 
Ny=20; % gridpoints in x  
Lx=2*pi % domain size in x
Ly=pi % domain size in y




maxorder = 19;

scalex=-2/Lx;
[x,DMx] = chebdif(Nx,Nx-1); 

scaley=-2/Ly;
[y,DMy] = chebdif(Ny,Ny-1); 



x=(x-1)/scalex;
y=(y-1)/scaley; 

[X,Y]=meshgrid(x,y);
f = cos(X).*sin(Y);


for i = 1:maxorder
    
%{
  
  We can then calculatethe numerical solution.
  
%}    



    
    
    for j = 1:maxorder


     
        
        dx=DMx(:,:,i)*scalex^(i+1);    
        Dx=kron(dx,eye(Ny));

  
        dy=DMy(:,:,j)*scaley^(j+1);
        Dy=kron(eye(Nx),dy);

        Df=reshape(Dy*(Dx*f(:)),Ny,Nx);
         
        
%{
  
  The following part calculate the analytical solution.
  
%}    
        
        
        i2 = i;
        while i2 > 4

        i2 = i2 - 4;

        end 
        
        
        
        j2 = j;
        while j2 > 4
 
            j2 = j2 - 4;
            
        end
        
        
        
        switch i2
            
            
            case 0 
                
                
                switch j2

                    case 0
                        df=cos(X).*sin(Y);

                    case 1
                        df=cos(X).*cos(Y);

                    case 2
                        df=cos(X).*(-sin(Y));

                    case 3
                        df=cos(X).*(-cos(Y));

                end            
            
            
            case 1
                
                
                switch j2

                    case 0
                        df=-sin(X).*sin(Y);

                    case 1
                        df=-sin(X).*cos(Y);

                    case 2
                        df=-sin(X).*(-sin(Y));

                    case 3
                        df=-sin(X).*(-cos(Y));
                        
                end
               
                
            case 2
                
                
                switch j2

                    case 0
                        df=-cos(X).*sin(Y);
                        
                    case 1
                        df=-cos(X).*cos(Y);
                        
                    case 2
                        df=-cos(X).*(-sin(Y));

                    case 3
                        df=-cos(X).*(-cos(Y));
                        
                end                
                
                
            case 3
                
                
                switch j2

                    case 0
                        df=sin(X).*sin(Y);
                        
                    case 1
                        df=sin(X).*cos(Y);
                        
                    case 2
                        df=sin(X).*(-sin(Y));

                    case 3
                        df=sin(X).*(-cos(Y));
                        
                end                    
                
                
        end
        
        
        
%{
  
  And then we compare the two solutions.
  
%}    
        
       
        
 
        erreur(i,j) = max(max(abs(Df-df)));

    end
    
end

i3 = linspace(1, maxorder, maxorder); 
j3 = linspace(1, maxorder, maxorder);  


% results

figure(1)
mesh(i3,j3,log10(erreur)); xlabel('Xorder'); ylabel('Yorder'); title('Dxiyj*f-fxiyj');


%{
![Visualisation of the error according to the derivative order](diff_2D_higer_order1.jpg)
![Visualisation of the error according to higer derivative order with an other angle of view](diff_2D_higer_order2.jpg)

The derivatives are two dimentional and we see that the error appears when we passe à separation line

![Visualisation of the error according to higer derivative order from top](diff_2D_higer_order3.jpg)

The error appears where it starts being blue and we se that it when the product of the levels of two level derivatives is near 15. That seems logical to us because we have to devide by the scale at each time we ad a level scalex^(i+1) and scaley^(j+1). That means the error appears when we get near the machine zero of variables or 10^16. Thats ad a counter-intuitive behaviour when we use a more refined grid. In fact, the error increase for the hight level dervatives because the step size between two nodes get lower and make the zero machine closer during the calculations. We can see in the following graphs the error for a twice and a half more refined grid:

![Visualisation of the error according to derivative order withe a more refined grid](diff_2D_higer_order4.jpg)
![Visualisation of the error according to derivative order withe a more refined grid from top](diff_2D_higer_order5.jpg)

#The sugestion answered to in diffmat_2D was:
Please test also higher order derivatives (third, fourth…) and show how the accuracy changes


%}
