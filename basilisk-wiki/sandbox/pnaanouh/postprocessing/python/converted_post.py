import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from math import sqrt 

def find_drops(curr_state,prev_state,dt,Ls):
    # function that correlates the drops between two samples 
    #
    #    Takes as input the data for the current sample as well as the previous
    #    sample, the time between samples and the domain size
    #
    #    Outputs a sorted version of the current state where each row in the
    #    current sate and same row in the previous state correspond to the same
    #    drop
    
    # initialize the sorted state to zero and the distance between the two
    # states to zero
    
    n_Drops=curr_state.shape[0]
 
    sorted_state=np.zeros_like(curr_state);
    Dist=np.ones((n_Drops,n_Drops))*Ls;
    
    for i in range(n_Drops): #loop over every drop in the current sate and ever drop in the previous state
        for j in range(n_Drops):
            nx=prev_state[j,2]+prev_state[j,4]*dt; #predicted position of the droplet
            ny=prev_state[j,3]+prev_state[j,5]*dt;
            dx=abs(curr_state[i,2]-nx); # difference between the prediction and a measured position
            dy=abs(curr_state[i,3]-ny);
            if(dx<Ls+prev_state[j,4]*dt and dx>Ls-prev_state[j,4]*dt): # taking into account the periodic boundary conditions
                dx=abs(dx-Ls);
            if(dy<Ls+prev_state[j,5]*dt and dy>Ls-prev_state[j,5]*dt):
                dy=abs(dy-Ls);
            Dist[i,j]=sqrt(dx**2+dy**2)

    Minimums=np.amin(Dist,0); # vectors containing the minimums of the row and collumns of Dist
    M=np.amin(Dist,1);
    no_overlap=np.zeros((n_Drops)); # vector used to ensure that each drop is only correlated to one other drop
    for i in range(n_Drops): #loop on all the drops in the current and previous state
        for j in range(n_Drops):
            if(Dist[i,j]==Minimums[j] and Minimums[j]==M[i]): #The two drops that are the closest to each other are the same drop
                sorted_state[j,:]=curr_state[i,:];
                no_overlap[i]=1;
    for i in range(n_Drops): # If a drop isn't correlated to any other one then just assign it to any drop that is free. This only matters in the case of coalescence and breakup 
        for j in range(n_Drops):
            if(sorted_state[j,6]==0 and no_overlap[i]==0):
                sorted_state[j,:]=curr_state[i,:];
                no_overlap[i]=1;
    return sorted_state

def V_fluctuations(log,log_gl,Ls) :
# Post-processing function for the droplet velocities
# 
#    Takes as inputs the log file containing the individual droplet
#    properties as well as the one containing the global properties for the
#    entire computational domain as well as the domain size
# 
#    Ouptuts the rms of the velocity fluctuations, the mean velocity, as well as the
#    fluctuation of the veloctiy from the mean as a function of time in both
#    the horizontal and vertical directions

# start off by turning the log file into something readable, seperate the
# different sample times and assotiate the samples to individual droplets


# separate the drops based on the timestep
    j=True;
    n=True;
    for i in range(0,log.shape[0]):
        if(log[i,6]>0.5):
            if(j):
                pg=log[i,:];
                j=False
            else:
                pg=np.vstack((pg,log[i,:]));
                
            #Drops[j,:,int((log[i,1]-log[0,1])*10+1)]=log[i,:];
            if(i<log.shape[0]-1):
                if(log[i,0]!=log[i+1,0]):
                    j=True;
                    if n:
                        n=False;
                        Drops=pg.copy();
                    else:
                        Drops=np.dstack((Drops,pg))
            else:
                Drops=np.dstack((Drops,pg))
                        
      
    # sort the drops into a consitent order
    Drops_sorted=np.zeros_like(Drops);
    Drops_sorted[:,:,0]=Drops[:,:,0];
    for i in range(1,Drops.shape[2]):
        dt=Drops[0,1,i]-Drops[0,1,i-1];
        Drops_sorted[:,:,i]=find_drops(Drops[:,:,i],Drops_sorted[:,:,i-1],dt,Ls);
    
    np.delete(Drops_sorted,Drops_sorted.shape[2]-1,2);
    
     
    vx=np.subtract(Drops_sorted[:,4,:],log_gl[0:log_gl.shape[0],5]).transpose();
    vy=np.subtract(Drops_sorted[:,5,:],log_gl[0:log_gl.shape[0],6]).transpose();
    
    vxm=np.mean(vx); #mean velocities
    vym=np.mean(vy);
    
    vxp=(vx-vxm); # the fluctuations are the difference between the instantneous values and the mean values
    vyp=(vy-vym);
    
    np.delete(vxp,np.arange(0,2000),0);  
    #remove the transient phase
    np.delete(vyp,np.arange(0,2000),0);
    
    
    vxrms=np.sqrt(np.mean((vx-vxm)**2)); # calculate the rms of the fluctuations
    vyrms=np.sqrt(np.mean((vy-vym)**2));

    return vxrms,vyrms,vxp,vyp,vxm,vym 

def r_distribution(log,Ls,D,rmax,dr):
    #   Computes the radial distribution function of the droplets
    #
    #   Takes as inputs the log file generated from basilisk, the domain size,
    #   Droplet diameter, maximum radius of the radial distribution function,
    #   and the thickness of one radial distribution shell
    
    #start off by turning the log file into something readable, seperate the
    #different sample times and assotiate the samples to individual droplets
    Drops=[];
    #separate the drops based on the timestep
    j=True;
    n=True;
    for i in range(0,log.shape[0]):
        if(log[i,6]>0.5):
            if(j):
                pg=log[i,:];
                j=False
            else:
                pg=np.vstack((pg,log[i,:]));
                
            #Drops[j,:,int((log[i,1]-log[0,1])*10+1)]=log[i,:];
            if(i<log.shape[0]-1):
                if(log[i,0]!=log[i+1,0]):
                    j=True;
                    if n:
                        n=False;
                        Drops=pg.copy();
                    else:
                        Drops=np.dstack((Drops,pg))
            else:
                Drops=np.dstack((Drops,pg))
                
    n_Drops=Drops.shape[0]
    n_timestep=Drops.shape[2]
      
    Edges=np.arange(D/100,rmax,dr); #sets the edges of the histigram bins starting at D/100 to avoid a peak at zero
    radii=np.zeros([len(Edges)-1,1]); 
    for j in range(len(radii)):
            radii[j]=(Edges[j+1]+Edges[j])/2;
            
    gav=np.zeros([1,len(Edges)-1]); # initialize the radial distrubution function to zero
    
    for i in range(1999,n_timestep): # for every sample aside form the first 2000 that are neglected due to the transients 
        mirrorx=np.zeros([n_Drops,3,3]); #take into account the symmetries accros the periodic boundaries 
        mirrory=np.zeros([n_Drops,3,3]);
        for m in range(-1,2):
            for n in range(-1,2):
                mirrorx[:,n+1,m+1]= Drops[:,2,i]+n*Ls;
                mirrory[:,n+1,m+1]= Drops[:,3,i]+m*Ls;
    
        g=np.zeros_like(gav) # initialize the instantaneous radial distibution function to zero 
        for j in range(n_Drops): # for every drop 
            Distances= np.sqrt((Drops[j,2,i]-mirrorx)**2+(Drops[j,3,i]-mirrory)**2); # distance between it and every other drop with symmetries
            result=np.histogram(Distances,bins=Edges); #radial distribution aroud the droplet
            g+=result[0]/n_Drops/(n_Drops/Ls**2); # instantaneous radial distribution is the average of the radial distibutions of all droplets divided by the number density of droplets
    
        for j in range(g.shape[1]): # scaled to take into account the increasing area of the shells with the increase in the radius
            g[0,j]/=np.pi*(Edges[j+1]**2-Edges[j]**2);
        gav+=g/(n_timestep-2000); # the radial distibution function is the mean of the instantaneous radial distribution functions
    gav=np.transpose(gav,[1,0])    
    
    #plots the radial distibution function
    plt.plot(radii,gav)
    plt.show()
    return [gav,radii]

def r_distribution2D(log,Ls,D,rmax,dr):
    #   Computes the 2-D pair distribution function of the droplets
    #
    #   Takes as inputs the log file generated from basilisk, the domain size,
    #   Droplet diameter, maximum radius of the radial distribution function,
    #   and the thickness of one radial distribution shell
    
    #start off by turning the log file into something readable, seperate the
    #different sample times and assotiate the samples to individual droplets
    Drops=[];
    #separate the drops based on the timestep
    j=True;
    n=True;
    for i in range(0,log.shape[0]):
        if(log[i,6]>0.5):
            if(j):
                pg=log[i,:];
                j=False
            else:
                pg=np.vstack((pg,log[i,:]));
                
            #Drops[j,:,int((log[i,1]-log[0,1])*10+1)]=log[i,:];
            if(i<log.shape[0]-1):
                if(log[i,0]!=log[i+1,0]):
                    j=True;
                    if n:
                        n=False;
                        Drops=pg.copy();
                    else:
                        Drops=np.dstack((Drops,pg))
            else:
                Drops=np.dstack((Drops,pg))
    n_Drops=Drops.shape[0]
    n_timestep=Drops.shape[2]
    
    Edges=np.arange(D/100,rmax,dr); #sets the edges of the histigram bins in the radial direction starting at D/100 to avoid a peak at zero
    Edges_a=np.linspace(-np.pi/2,np.pi/2,num=251); #sets the edges of the histigram bins in the angular direction
    thetas=np.zeros([len(Edges_a)-1,1]);
    radii=np.zeros([len(Edges)-1,1]); 
    for j in range(len(radii)):
            radii[j]=(Edges[j+1]+Edges[j])/2;
    for j in range(len(thetas)):
            thetas[j]=(Edges_a[j+1]+Edges_a[j])/2;
    gav=np.zeros([len(Edges)-1,len(Edges_a)-1]); # initialize the pair distrubution function to zero
    
    for i in range(1999,n_timestep): #start averaging at t=200 to avoid transients
        mirrorx=np.zeros([n_Drops,3,3]);
        mirrory=np.zeros([n_Drops,3,3]);
        for m in range(-1,2): #symetries due to periodic boundaries
            for n in range(-1,2):
                mirrorx[:,n+1,m+1]= Drops[:,2,i]+n*Ls;
                mirrory[:,n+1,m+1]= Drops[:,3,i]+m*Ls;
    
        g=np.zeros_like(gav); # initialize the instantaneous pair distibution function to zero 
        for j in range(n_Drops): # for every drop
            Distances= np.sqrt( (Drops[j,2,i]-mirrorx)**2+(Drops[j,3,i]-mirrory)**2);# distance between it and every other drop with symmetries
            Angles= np.arcsin((Drops[j,3,i]-mirrory)/Distances); # Angles between it and every other drop
            result=np.histogram2d(Distances.flatten(),Angles.flatten(),bins=(Edges,Edges_a)); #pair distribution aroud the droplet
            g+=result[0]/n_Drops/(n_Drops/Ls**2);  # instantaneous pair distribution is the average of the pair distibutions of all droplets divided by the number density of droplets
    
        for j in range(g.shape[0]): # scaled to take into account the area of the bins
            g[j,:]/=(np.pi*(Edges[j+1]**2-Edges[j]**2))/250;
        gav+=g/(n_timestep-2000);# the pair distibution function is the mean of the instantaneous pair distribution functions
    #plot of the pair distibution function
    Theta,R=np.meshgrid(thetas,radii);
    X=R*np.cos(Theta);
    Y=R*np.sin(Theta);
    
    
    plt.figure(figsize=(8, 16))
    plt.pcolormesh(X,Y,gav)
    plt.jet()
    plt.show()
    #figure('Name','Surface','Position',[0 0 600 920])
    #surf(X,Y,gav,'EdgeColor','none','FaceColor','flat'); 
    #view(2);  
    #colorbar
    #colormap jet
    return gav, X,Y

def drop_statistics(log,Ls,D,lDref,lt_range):
         #Primary Post-proccessing function for the calculation of droplet
        #statistics
        #   Not all of the statistics that are computed are outputed so feel free
        #   to either remove the parts you don't want or to add them as outputs 
        #   
        #   This function takes as an input a log file generated in basilisk ( see
        #   sandbox/pnaanouh/rising-suspension.c), the domain size of the
        #   simulation, the droplet diameter, the number of different contact
        #   threshold distances to consider, and the number of bins for the
        #   calculation of the histogram of the contact time pdf
        #
        #   This function has may outputs:
        #       t: the time at every sample in the simulation
        #       t_bins: the values of the bins for the histogram of the contact 
        #           time pdf at different contact threshold distances
        #       n_int_start: the total number of collisions counted at the start of
        #           the collisions at different contact threshold distances
        #       t_range: the bins of the histogram of the contact
        #           time pdf
        #       Dist_ref: the range of the contact threshold distance
        #       vx: the instananeous horizontal velocity of the dispersed phase as 
        #           a function of time
        #       vy: the instananeous vertical velocity of the dispersed phase as a 
        #           function of time
        #       vt: the total volume of the dispersed phase as a function of time
        #       Dist_signal: Distance between each droplet and its nearest
        #           neigbhour as a function of time
        #       t_cont: the durations of all the measured collisions for different
        #           contact threshold distances
        #       t_min: the durations of all the measured collisions for different
        #           contact threshold distances with the assumption of binary
        #           collisions
        #       Dist_free: mean free path between collisions for different contact threshold
        #           distances
        #       Dist_not_free: mean distance travelled during collisions path for  
        #           different contact threshold distances
        #       n_min:the total number of collisions at different contact threshold
        #           distances  with the assumption of binary collisions
        #       Dist: the distance between every pair of droplets at every sample
        #           time
        
        
    #start off by turning the log file into something readable, seperate the
    #different sample times and assotiate the samples to individual droplets
    Drops=[];
    #separate the drops based on the timestep
    j=True;
    n=True;
    for i in range(0,log.shape[0]):
        if(log[i,6]>0.5):
            if(j):
                pg=log[i,:];
                j=False
            else:
                pg=np.vstack((pg,log[i,:]));
                
            #Drops[j,:,int((log[i,1]-log[0,1])*10+1)]=log[i,:];
            if(i<log.shape[0]-1):
                if(log[i,0]!=log[i+1,0]):
                    j=True;
                    if n:
                        n=False;
                        Drops=pg.copy();
                    else:
                        Drops=np.dstack((Drops,pg))
            else:
                Drops=np.dstack((Drops,pg))
                        
      
    # sort the drops into a consitent order
    Drops_sorted=np.zeros_like(Drops);
    Drops_sorted[:,:,0]=Drops[:,:,0];
    for i in range(1,Drops.shape[2]):
        dt=Drops[0,1,i]-Drops[0,1,i-1];
        Drops_sorted[:,:,i]=find_drops(Drops[:,:,i],Drops_sorted[:,:,i-1],dt,Ls);
    
    n_Drops=Drops_sorted.shape[0]
    n_timestep=Drops_sorted.shape[2]
    
    #np.delete(Drops_sorted,n_timestep-1,2);
    
    # invert dimensions 1 and 3 for added convenience
    
    Drops_final=np.transpose(Drops_sorted.copy(),[2,1,0]);
    
    
    
    t=Drops_final[:,2,1];
    #get the drop velocities
    
    #Average velocity at every instant
    vx=np.sum(Drops_sorted[:,4,:]*Drops_sorted[:,6,:],0)/np.sum(Drops_sorted[:,6,:],0);
    
    vy=np.sum(Drops_sorted[:,5,:]*Drops_sorted[:,6,:],0)/np.sum(Drops_sorted[:,6,:],0);
    
    #Total volume in time
    vt=np.sum(Drops_sorted[:,6,:],0);
    
    
    #calculate the distance between all drops and thus the collision freq and
    #contact time
    
    #initialize the matrix to store the distances between all the drops at
    #every timestep
    Dist=np.zeros([n_Drops,n_Drops,n_timestep])
    #initialize the matrix to store the distance travelled by each drop at
    #every timestep
    dx=np.zeros([n_Drops,n_timestep])
    #initialize the matrix to store the distances between each drop and its
    #nearest neighbour at every timestep
    Dist_signal=np.zeros([n_Drops,n_timestep])
    
    #set the range for the contact threshold distance between 0 and 4 diameters
    #with lDref number of points
    Dist_ref=np.linspace(0,4*D,num=lDref);
    
    #initalize the vectors to store the total number of collisions n_int_start
    #counted at the start of the collsions, n_int_end counted at the end of the
    #collisions, and n_min calculated based on the minimum distance signal
    n_int_start=np.zeros_like(Dist_ref);
    n_int_end=np.zeros_like(Dist_ref);
    n_min=np.zeros_like(Dist_ref);
    
    #initalize the vectors to store the total contact time, the total distance
    #traveled outside of collisions, and the total distance travalled during 
    #collisions  
    t_total=np.zeros_like(Dist_ref);
    Dist_free=np.zeros_like(Dist_ref)
    Dist_not_free=np.zeros_like(Dist_ref)
    
    #initialize variables to store the durations of individual collisions 
    t_cont=[]; # stores the collision durations
    t_int=[]; # stores the time between the end of the collisions
    t_min=[]; # stores the collision durations with the assumption of binary collisions
    
    #iterate on the range of the contact reference distances for each pair of
    #drops at every sample time
    for n in range(len(Dist_ref)):
        #initialize the matixies used to calculate the contact time 
        
        #stores the time spent in the current collision 
        t_c3=np.zeros([n_Drops,n_Drops])
        t_i3=np.zeros([n_Drops,n_Drops])
        t_m3=np.zeros([n_Drops])
        
    #     %stores the matrix that contains the duration of all the collisions it
    #     %is neccessary due to the change in the number of collisions with the
    #     %contact threshold distance 
    #     t_c2=struct;
    #     t_i2=struct;
    #     t_m2=struct;
        
        #stores the duration of all the collisions 
        t_c=[]
        t_i=[]
        t_m=[]
        
        #matix to store the state of each droplet pair( colliding 1, not
        #colliding 0) at the current sample and the previous sample
        ck=np.zeros([n_Drops,n_Drops,2]);
        #matix to store the state of each droplet( colliding 1, not
        #colliding 0) at the current sample and the previous sample
        ck2=np.zeros([n_Drops,2]);
        
        for i in range(n_timestep): #loop on sample time
            a=i%2; # current sample 
            b=(i-1)%2; #previous sample 
            Dist[:,:,i]=Ls; #initialize the distance to the domain size, this is unnecessary 
            for j in range(n_Drops): # loop on each droplet
                Dist[j,j,i]=Ls;  #the distance on the diagonal is set to the domain size, again unnecessary
                mv=np.zeros([3,3]); # matrix to store the distance travelled by the drop between the current sample 
                #and the previous sample with symmetries accross the periodic boundaries
                if(i>0): #calculation of the distance travelled 
                    for l in range(-1,2):
                        for m in range(-1,2):
                            mv[l+1,m+1]=sqrt((Drops_sorted[j,2,i]-Drops_sorted[j,2,i-1]-Ls*l)**2+(Drops_sorted[j,3,i]-Drops_sorted[j,3,i-1]-Ls*m)**2);
    
                dx[j,i]=np.amin(mv); # the actual distance travelled is the minimum of all the symmetries
                
                
                for k in range(j+1,n_Drops): # loop of every pair of drops starts from j+1 to get upper triangular matrix 
                    if(Drops_sorted[k,6,i]!=0 and Drops_sorted[j,6,i]!=0): # checks that both droplets actually exist can be relevant if there is breakup or coalescence due to how Drops_sorted is constructed
                        DD=np.zeros([3,3]); # matrix to store the distance between the two drops plus the symmetries
                        for l in range(-1,2): #calculate the distance
                            for m in range(-1,2):
                                DD[l+1,m+1]=sqrt((Drops_sorted[k,2,i]-Drops_sorted[j,2,i]-Ls*l)**2+(Drops_sorted[k,3,i]-Drops_sorted[j,3,i]-Ls*m)**2);
                        Dist[j,k,i]=np.amin(DD); # the actual distance is the minimum of the symmetries
                        Dist[k,j,i]=Dist[j,k,i]; # adding the symmetric of the distance in the distance matrix
                        
                        ck[j,k,a]=Dist[j,k,i]<=Dist_ref[n]; # store if the droplets are colliding or not
                                
                        if(i>0):
                            if(ck[j,k,b] and (not ck[j,k,a])): #if the droplets were colliding and are no longer colliding 
                                n_int_end[n]+=1; # count the end of a collision 
                                t_c.append(t_c3[j,k]); # store the duration of the collision 
                                t_c3[j,k]=0; # initialize the collision duration to zero 
                                t_i.append(t_i3[j,k]); # store the time between collisions 
                                t_i3[j,k]=0; # initialize the time between collisions to zero
                            elif( (not ck[j,k,b]) and ck[j,k,a]): #if the droplets were not colliding and are now colliding
                                n_int_start[n]+=1; # count the start of a collision
                Dist_signal[j,i]=np.amin(Dist[j,:,i]); #computes the  minimum distance
                ck2[j,a]=Dist_signal[j,i]<=Dist_ref[n]; # checks if the droplet is colliding based on the assumption of binary collisions
                Dist_free[n]+=(1-ck2[j,a])*dx[j,i]; # adds to the total distance travelled while colliding 
                Dist_not_free[n]+=ck2[j,a]*dx[j,i]; # adds to the total distance travelled without colliding 
                if(i>0):
                    if(ck2[j,a] and (not ck2[j,b])): # at the start of the collision 
                        n_min[n]+=1; # count the collision
                        t_m.append(t_m3[j]); # store the duration of the previous collision 
                        t_m3[j]=0; # initialize the collision duration to zero
            #matrix computation of the collision durations, and the time
            #between collisions
            t_c3+=ck[:,:,a]*dt;  
            t_i3+=dt;
            t_m3+=ck2[:,a]*dt;
            
        #the collision durations are stored in structures
        #t_c2.n=t_c;
        #t_i2.n=t_i;
        #t_m2.n=t_m;
        #the structures are added to matrixies
        t_cont.append(t_c);
        t_int.append(t_i);
        t_min.append(t_m);
    
    #the mean free and mean not free paths are computed to be the total
    #distance divided by the number of droplets divided by the number of
    #sections 
    for i in range(len(Dist_ref)):
        if n_min[i]>0:
            Dist_free[i]/=n_Drops*n_min[i];
            Dist_not_free[i]/=n_Drops*n_min[i];
        else:
            Dist_free[i]=0
            Dist_not_free[i]=0
    
    
    # calculation of the pdf of the contact time
    
    t_range=np.linspace(0, 60, num=lt_range)
    
    t_bins=[]
    t_bins2=[]
    
    for i in range(len(Dist_ref)):
        t_bins.append(np.histogram(t_cont[i],t_range,density=True))
        t_bins2.append(np.histogram(t_int[i],t_range,density=True))
    
    return t,t_bins,n_int_start,t_range,Dist_ref,vx,vy,vt,Dist_signal,t_cont,t_min,Dist_free,Dist_not_free,n_min,Dist

Ls=2*11.4411;
D=1;

Ga10=pd.read_csv("./Results/MuR_042_PHI_15_Ga_10.out",sep='\s+',header=None,skiprows=1,skipinitialspace=True).to_numpy();
Ga10gl=pd.read_csv("./Results/MuR_042_PHI_15_Ga_10_gl.out",sep='\s+',header=None,skiprows=1,skipfooter=3).to_numpy();


[gav1d10,radii] = r_distribution(Ga10,Ls,D,3,0.01);
[gav2d10, X,Y] = r_distribution2D(Ga10,Ls,D,3,0.01);
[vxrms10,vyrms10,vxp10,vyp10,vxm10,vym10]=V_fluctuations(Ga10,Ga10gl,Ls);
[t10,t_bins10,n_int_start10,t_range10,Dist_ref10,vx10,vy10,vt10,DS10,t_cont10,t_min10,Dist_free10,Dist_not_free10,n_min10] = drop_statistics(Ga10,Ls,D,5,121);
