import numpy as np
import matplotlib.pyplot as plt
from IPython.display import clear_output
from IPython.display import display
from sklearn.linear_model import LinearRegression
import os 


######## Defining the Fragmentation Mechanism ############

##Each fragment is described by two numbers, its length that is between 0 and 1, and a neutral marker uniformly sampled between 0 and 1. 

#Whenever a fragment splits into smaller pieces, the smaller fragments inherit the same neutral marker as the original fragment. 

#Therefore, the whole Fragment List (Frag_List) will have two attributes, Frag_lens and Frag_labs for lengths and labels

#Here, we define multiple functions that are essential to calculating the rates at which immigration, fragmentation and exit events occur.

#These functions are then used to define frag_mechan which determines whether fragmentation, immigration or exit occurs for a list of fragments, and what happens to the list of fragments subsequently.

class Frag_List:
    
    def __init__(self, Frag_lens, Frag_labs):
        self.Frag_lens=Frag_lens
        self.Frag_labs=Frag_labs
        
    def __str__(self):
        information = f"Fragmentation length = {self.Frag_lens}, Fragmentation Labels = {self.Frag_labs}"
        return information

    def age(self, T):
        "age function gives the age of the fragments, T= current time"
        return np.asarray([T-lab for lab in self.Frag_labs])    
   
     
    def frag_time(self,frag_speed, alpha_f):
        "frag_time gives the next fragmentation time given list of fragments. fragments: np-array of all fragments lengths, speed: rate at which fragments split, alpha_f: self-similarity index of fragmentation"
        frag_time = np.random.exponential(1/(np.sum(frag_speed*self.Frag_lens**(alpha_f))))
        #if frag_time < 0.00001: #prevent fragmentation time that is too small
            #frag_time = 0.00001
        return frag_time
    
    def imm_time(self,x_max, imm_rate):
        "imm_time gives next immigration time, which is exponential with imm_rate"
        return np.random.exponential(1/(x_max*imm_rate))
    
    def imm(self, new_Frags):
        "Add new_Frags into existing fragment lists, fragment"
        self.Frag_lens= np.append(self.Frag_lens,new_Frags.Frag_lens)
        self.Frag_labs= np.append(self.Frag_labs,new_Frags.Frag_labs)    
    
    def exit_params(self, ex_max, ex_slope):
        "ex_params gives the time for first exit of a fragment from a Frag_list, the exit rate of each fragment according to its length for a list of fragments, and the renormalised exit probability Here we assume the exit rate linearly decreases with fragment lengths. ex_max: maximum rate, ex_slope: slope of decay"
        exit_rates=ex_max*np.ones(len(self.Frag_lens))-ex_slope*self.Frag_lens
        exit_rates[exit_rates<=0]=0.0001 #Prohibit negative exit rates
        exit_rates[exit_rates>=1000000]=1000000 #Prohibit huge exit rates
        exit_time = np.random.exponential(1/np.sum(exit_rates))
        exit_prob = exit_rates/np.linalg.norm(exit_rates,1)
        return exit_time, exit_rates, exit_prob
    
    def exit_power_par(self, B_ex, alpha_e, x_max):
        "exit_power_par gives the time for first exit of a fragment from a Frag_list, the exit rate of each fragment according to its length for a list of fragments, and the renormalised exit probability Here we assume the exit rate decreases with fragment lengths with a power alpha_e. B_ex: exit boundary, alpha_e: power of decay"
        #frag_lens = self.Frag_lens
        #frag_lens[frag_lens >= B_ex] = 0
        #exit_rates = frag_lens
        #exit_rates[frag_lens != 0] = frag_lens[frag_lens != 0]**(alpha_e) - B_ex ** alpha_e
        exit_rates = x_max** (-alpha_e) * (self.Frag_lens**(alpha_e) - B_ex ** alpha_e)
        if B_ex == np.inf:
            exit_rates = x_max** (-alpha_e) * self.Frag_lens**(alpha_e)
        #exit_rates[exit_rates>=10000]=10000 #Prohibit large exit rates
        exit_rates[exit_rates<=0]=0.000000000001
        exit_time = np.random.exponential(1/np.sum(exit_rates))
        exit_prob = exit_rates/np.linalg.norm(exit_rates,1)
        return exit_time, exit_rates, exit_prob

   
    def exit(self, ex_max, ex_slope):
        "This function deletes a fragment from the list with probability proportional to its exit rates"
        exit_prob=self.exit_params(ex_max, ex_slope)[2]
        frag_gone=np.random.choice(self.Frag_lens.size,1, p=exit_prob)[0] #the gone fragment is chosen with probability proportional to exit_prob
        #print(frag_gone)
        self.Frag_lens=np.delete(self.Frag_lens,frag_gone) #gone fragment deleted from existing system
        self.Frag_labs=np.delete(self.Frag_labs,frag_gone) #gone fragment deleted from existing system
    
    def exit_power(self, B_ex, alpha_e, x_max):
        "This function deletes a fragment from the list with probability proportional to its exit rates (given power)"
        exit_prob=self.exit_power_par(B_ex, alpha_e, x_max)[2]
        frag_gone=np.random.choice(self.Frag_lens.size,1, p=exit_prob)[0] #the gone fragment is chosen with probability proportional to exit_prob
        #print(frag_gone)
        self.Frag_lens=np.delete(self.Frag_lens,frag_gone) #gone fragment deleted from existing system
        self.Frag_labs=np.delete(self.Frag_labs,frag_gone) #gone fragment deleted from existing system
    
    def split(self, N_split, a, b, alpha_f, low_threshold=0.001):
        "This function splits a fragment into N_split many smaller fragments with beta distribution with parameters a, b"
        frag_rates=self.Frag_lens**(alpha_f)
        frag_rates[frag_rates>=100000]=100000 #prohibits large fragmentation rates 
        frag_prob = frag_rates/np.linalg.norm(frag_rates,1)
        i = np.random.choice(self.Frag_lens.size, 1, p=frag_prob)[0] #i-th fragment fragments
        #print(i)
        gone_frag_lens=self.Frag_lens[i]
        inherited_lab = self.Frag_labs[i]
        frag_interval=np.array([0,self.Frag_lens[i]]) 
        #frag_site=np.sort(np.append(np.random.uniform(0,self.Frag_lens[i],N_split-1),frag_interval))
        frag_site=np.sort(np.append(self.Frag_lens[i]*np.random.beta(a,b,N_split-1),frag_interval)) #splitting positions are distributed as Beta(a,b) 
        new_frag_lens=np.diff(frag_site)
        new_frag_lens = new_frag_lens[new_frag_lens >= low_threshold] #only retain fragments with length > low_threshold
        new_frag_labs=self.Frag_labs[i]*np.ones(len(new_frag_lens))
        self.Frag_lens = np.append(np.delete(self.Frag_lens,i),new_frag_lens)
        self.Frag_labs = np.append(np.delete(self.Frag_labs,i),new_frag_labs)
    
    def FragExIm(self, ex_max, ex_slope, frag_speed, N_split, a,b, alpha_f, Nimm, imm_rate, Initial_Time, lin_ex=True, B_ex = 1, alpha_e = 0, x_max=1, norm_imm=False, exp_imm=False, mu_i=0):
        "FragExIm picks a fragmentation / exit / immigration by sampling the next event time, then proceeds to apply frag/imm/exit process to the list of fragments. lin_ex = True if we assume linearity for exit rate."
        event_times = np.asarray([self.frag_time(frag_speed, alpha_f), self.exit_params(ex_max,ex_slope)[0], self.imm_time(x_max, imm_rate)])
        if lin_ex == False:
            event_times = np.asarray([self.frag_time(frag_speed, alpha_f), self.exit_power_par(B_ex, alpha_e, x_max)[0], self.imm_time(x_max, imm_rate)])
        
        #print(event_times)
        #print(np.where(event_times == np.amin(event_times))[0][0])
        next_event=np.where(event_times == np.amin(event_times))[0]
        if next_event[0]==0:
            Event= "fragmentation"
            self.split(N_split,a,b, alpha_f)
        elif next_event[0]==1:
            Event= "exit"
            if lin_ex == True:
                self.exit(ex_max, ex_slope)
            elif lin_ex == False:
                self.exit_power(B_ex, alpha_e, x_max)
        else:
            Event= "Immigration"
            new_Frags = Frag_List(np.random.uniform(0,x_max,Nimm),np.ones(Nimm)*(Initial_Time + np.amin(event_times))) ### New Fragments given label at birth time = Initial_Time + event_times 
            if norm_imm:
                new_Frags = Frag_List(np.random.normal(mu_i,0.1*mu_i,Nimm),np.ones(Nimm)*(Initial_Time + np.amin(event_times))) ### New Fragments given label at birth time = Initial_Time + event_times 
            if exp_imm:
                new_Frags = Frag_List(np.random.exponential(mu_i,Nimm),np.ones(Nimm)*(Initial_Time + np.amin(event_times))) ### New Fragments given label at birth time = Initial_Time + event_times 
            #new_Frags = Frag_List(np.ones(Nimm)*x_max,np.ones(Nimm)*(Initial_Time + np.amin(event_times))) ### New Fragments are of length 1 
            self.imm(new_Frags)
        #print(f"Event = {Event}, Event Time = {np.amin(event_times): .3f}")    
        return Event, np.amin(event_times)
    
    def FragEx(self, ex_max, ex_slope, frag_speed, N_split, a,b, alpha_f, lin_ex=True, B_ex = 1, alpha_e = 0, x_max=1):
        "FragExIm picks a fragmentation / exit by sampling the next event time, then proceeds to apply frag/ex process to the list of fragments. lin_ex = True if we assume linearity for exit rate."
        event_times = np.asarray([self.frag_time(frag_speed, alpha_f), self.exit_params(ex_max,ex_slope)[0]])
        if lin_ex == False:
            event_times = np.asarray([self.frag_time(frag_speed, alpha_f), self.exit_power_par(B_ex, alpha_e, x_max)[0]])
        
        #print(event_times)
        #print(np.where(event_times == np.amin(event_times))[0][0])
        next_event=np.where(event_times == np.amin(event_times))[0]
        if next_event[0]==0:
            Event= "fragmentation"
            self.split(N_split,a,b, alpha_f)
        elif next_event[0]==1:
            Event= "exit"
            if lin_ex == True:
                self.exit(ex_max, ex_slope)
            elif lin_ex == False:
                self.exit_power(B_ex, alpha_e, x_max)
        return Event, np.amin(event_times)
    
############################################################################################################################################
################Fragmentation Statistics###########
#Here we define functions that give important statistics on a list of fragmentation lengths.
#They include a function that gives the cumulative distribution (or count) of fragment lengths of the list, a function that produces the empirical density / empirical mass of the fragment distribution, and a function that produces the Kolomogorov-Smirnov Distance between two empirical distributions / empirical mass.

    def mean(self):
        return np.mean(self.Frag_lens)
    
    def cum_dist(self, N_Group, x_max, count=False): 
        "cum_dist gives the empirical cumulative distribution function of fragment lengths given a list of fragments with N_Group data points, if count=True, cumulative count is given instead"
        N_tiles=np.linspace(0,x_max,N_Group+1)
        #print(N_tiles)
        cum_count=np.asarray([np.count_nonzero(self.Frag_lens<=i) for i in N_tiles])
        cum_dist=cum_count/cum_count[-1]
        if count:
            output = cum_count
        else:
            output = cum_dist
        return  output
    
    def density(self, N_Group,x_max, count=False):
        "density gives the empirical density of a sample, if count=True, empirical mass function is given instead"
        #return np.diff(self.cum_dist(N_Group,x_max,count))/(x_max/N_Group) #To calculate density, we divide count over length of bin x_max / N_Group
        return np.diff(self.cum_dist(N_Group,x_max,count)) #To calculate density, we divide count over length of bin x_max / N_Group        

def KS_Dist(cum_dist1,cum_dist2):
    "Kolmogorov-Smirnov Distance of two empirical distributions"
    return np.max(np.abs(cum_dist1-cum_dist2))/min(cum_dist1[-1],cum_dist2[-1])

def tail_loglog_fit(density, x_max, B_ex, peak=True):
    "Returns the best fit power law decay line for the decreasing tail of an observed density function"
    x_axis = np.linspace(0,x_max,len(density)+1)
    #print(x_axis)
    ex_index = int(round(len(density)*B_ex/x_max)) #index of where exit boundary should be
    tail_density = density[ex_index: ]
    print(len(density), ex_index, tail_density)
    tail_x_axis = x_axis[ex_index: -1 ]
    #print(f'density: {density}')
    #print(f"x_axis: {x_axis}")
    #print(f"ex_index: {ex_index}")
    #print(f"tail density: {tail_density}")
    #print(f"tail_x_axis: {tail_x_axis}")
    tail_density[tail_density <= 0.001] = 0.001
    tail_x_axis[tail_x_axis <= 0.01] = 0.001
    log_den = np.log(tail_density)
    log_x_axis = np.log(tail_x_axis)
    if peak == True:
        log_den_trans = log_den - log_den[0]  #shifting the density point at B_ex to origin
        log_x_axis_trans = log_x_axis - log_x_axis[0] #shifting the point to origin
        lm = LinearRegression(fit_intercept = False) #apply a linear model with no intercept fit
        lm.fit(log_x_axis_trans.reshape(-1, 1), log_den_trans)
        [beta] = lm.coef_
        c = log_den[0] - beta * log_x_axis[0]
    else:
        beta, c = np.polyfit(log_x_axis, log_den, 1)
      
    return log_den, log_x_axis, beta, c

def tail_l_l_fit(x_axis, density, i):
    log_x_axis = np.log(x_axis + 0.01)
    log_den = np.log(density + 0.01)
    tail_x = log_x_axis[i:]
    tail_den = log_den[i:]
    beta, c = np.polyfit(tail_x,tail_den,1)
    return tail_x, tail_den, beta, c 

def ex_reg_loglog_fit(density, x_max, B_ex, frag_speed, ex_slope):
    "Returns the best fit power law decay line for an observed density function in the exit region"
    x_axis = np.linspace(0,x_max,len(density)+1)
    ex_index = int(len(density)*B_ex)+1 #index of where exit boundary should be
    exit_density = density[: ex_index ]
    exit_x_axis = x_axis[:ex_index ] + 0.001
    ex_log_den = np.log(exit_density+0.001)
    ex_log_x_axis = np.log(B_ex+ (frag_speed/ex_slope-1)* exit_x_axis)
    log_den_trans = ex_log_den - ex_log_den[-1]  #shifting the density point at B_ex to origin
    log_x_axis_trans = ex_log_x_axis - ex_log_x_axis[-1] #shifting the point to origin
    lm = LinearRegression(fit_intercept = False) #apply a linear model with no intercept fit
    lm.fit(log_x_axis_trans.reshape(-1, 1), log_den_trans)
    [gamma] = lm.coef_
    d = ex_log_den[-1] - gamma * ex_log_x_axis[-1]
    return exit_x_axis, ex_log_den, ex_log_x_axis, gamma, d

## Stationary Distribution Simulator
#The stationary distribution simulator simulates the Fragmentation with Immigration and Exit Mechanism until it reaches mixing time.
#This is achieved by comparing the current state of the Markov process with the historical average state of the Markov process.
# When they are close in the Kolmogorov-Smirnov Distance, we will terminate the process and print out the stationary distribution. 
#In the following codes, you can choose to print out the state of the Markov chain at each time step if you decide not to hide the print functions.
#You can also set weighted


def stationary(x,ex_max, ex_slope, frag_speed, N_split, a,b, alpha_f, Nimm, imm_rate, Time_scale, N_Group, KS_threshold, Step_Max, lin_ex=True, B_ex = 1, alpha_e = 0, weighted= False, plot=False, cfDNA= False, x_max=1, xlim_max=1, save_fig=False, norm_imm=False, exp_imm=False, mu_i=0):
    "This function gives the stationary distribution of FragImmEx process if it reaches so before Step_Max time steps.\n Time_scale is the time we run our Markov process before updating KS distance. \n N_Group gives the number of groups we put our fragments into according to their fragment lengths. \n KS_threshold gives the threshold for the KS distance to be close.\n weighted gives weighted mean for faster convergence. \n plot=True plots frequency distribution at each time step. \n cfDNA = True rescales x axis to represent simulated cfDNA data."
    ##################Initialisation###########################################
    T=0
    N_Event=0
    run_cum_count = x.cum_dist(N_Group,x_max, count=True)
    run_dens = x.density(N_Group, x_max, count=True)
    KS_dist = 10
    Time_steps = 0
    KS_dist_mean = 1
    Dist_mean_list = []
    if save_fig == True:
       parent_dir = r'../Simulation Output/Stationary Plots/'
       directory =  str(ex_max) +'_'+ str(ex_slope) +'_'+  str(frag_speed)+'_'+  str(N_split) +'_'+  str(a) +'_'+ str(b) +'_'+  str(alpha_f) +'_'+ str( Nimm) +'_'+  str(imm_rate) +'_'+ str(Time_scale) +'_'+ str(N_Group) +'_'+  str(KS_threshold) +'_'+  str(Step_Max)
       path = os.path.join(parent_dir, directory)
       os.mkdir(path)

    if weighted == True:
         w_run_cum_count = x.cum_dist(N_Group,x_max, count=True)
         w_run_dens = x.density(N_Group, x_max, count=True)
         w_KS_dist_mean = 1
         w_Dist_mean_list = []
    ############################# Simulation until Stationary Distribution is achieved or Time_steps > Step_Max ##########################################
    ################### Simulation Ends when (KS_distance between current state and mean state)/(KS_distance in first time step) < ratio or KS dist < threshold ####################

    while  KS_dist_mean > KS_threshold and Time_steps < Step_Max:
        t=0
        Time_steps +=1
      ################### Running FragExIm simulation for Time_scale unit of time, Tracking Event Time and Evene numbers #######################
        while t < Time_scale:
            Event, Event_time = x.FragExIm(ex_max, ex_slope, frag_speed, N_split, a,b, alpha_f, Nimm, imm_rate,T,lin_ex, B_ex, alpha_e, x_max = x_max, norm_imm= norm_imm, exp_imm=exp_imm, mu_i = mu_i)
            #print( f"{Event} occurs after time {Event_time}")
            #Time_db.append(Event_time)
            N_Event +=1
            t += Event_time
            #print(T)
            T += Event_time
                 
        #Recording new density and new cumulative count, adding it to run_cum_count and average over time to give mean_dens and mean_cum_count
        new_dens = x.density(N_Group, x_max, count=True) 
        new_cum_count = x.cum_dist(N_Group,x_max,count=True)
        run_dens += new_dens
        run_cum_count += new_cum_count 
        mean_dens = run_dens / Time_steps
        mean_cum_count = run_cum_count / Time_steps
        # We calculate the KS distance between the mean cum count and the new cum count, and update the KS_dist to a list.
        KS_dist_mean = KS_Dist(mean_cum_count,new_cum_count)
        
        # We fix the KS_dist to be 0.5 for the first 5 time steps to ensure the simulator will not terminate in the first step. 
        if Time_steps <= 5:
            KS_dist_mean = 0.5  
        Dist_mean_list.append(KS_dist_mean)

        # For weighted models we weight each new output with the time_steps so more recent fragment profile contributes more to running mean fragment profile.
        if weighted == True:
           w_run_dens += Time_steps*new_dens
           w_run_cum_count += Time_steps*new_cum_count
           w_mean_dens = 2*w_run_dens/(Time_steps*(Time_steps+1))
           w_mean_cum_count = 2*w_run_cum_count/(Time_steps*(Time_steps+1))
           w_KS_dist_mean = KS_Dist(w_mean_cum_count,new_cum_count)
           if Time_steps <= 5:
            w_KS_dist_mean = 0.5
           w_Dist_mean_list.append(w_KS_dist_mean)
  


                

        if plot == True:
            fig, ax = plt.subplots(3,figsize=(10,10),frameon=True)
            #############Plotting New Density Against Mean Density over Time ############################
            ax[0].bar(np.linspace(0,x_max,N_Group),new_dens, label=f"Density at time {T: .1f}", color = 'salmon', alpha = 0.2, width=x_max/N_Group)
            ax[0].scatter(np.linspace(0,x_max,N_Group),mean_dens, label = f"Running Mean Density at time {T: .1f}",marker='x', color = 'r', alpha = 0.5)
            ax[0].set_xlim(0,xlim_max)
            ax[0].set_xlabel('Fragment length')
            ax[0].set_ylabel('Fragment count')
            #ax[0].set_ylim(0,2)
            ax[0].set_title('Fragment Lengths Frequency')
            #############Plotting New Cumulative Count Against Mean Cumulative Count over Time #################
            #ax[1].plot(np.linspace(0,x_max,N_Group+1),new_cum_count, label=f"CDF at time {T: .1f}")
            #ax[1].plot(np.linspace(0,x_max,N_Group+1),mean_cum_count, label=f"Running Mean CDF at time {T: .1f}")
            ##############Plotting KS Distance between Existing State and Mean State over Time #################
            ax[1].plot(np.linspace(Time_scale,Time_scale*len(Dist_mean_list),len(Dist_mean_list)),Dist_mean_list, label="Distance to Running Mean")
            ax[1].set_ylim(0,0.75)
            ax[1].set_title('KS Distance to Running Mean')
            #print(N_Event)
            #############Plotting Distribution of Ancestor (Labels) ########################################
            #ax[2].hist(x.Frag_labs,19)
            ax[2].hist(x.Frag_labs,Time_steps)
            #ax[2].set_ylim(0,imm_rate*Time_scale*x_max)
            #ax[2].legend(loc=1, prop={'size': 8})
            ax[2].set_title('Immigration Time Distribution')
            if weighted == True:
                ax[0].scatter(np.linspace(0,x_max,N_Group),w_mean_dens, label = f"weighted Running Mean Density at time {T: .1f}",marker='.', color = 'k', alpha = 0.5)
                #ax[1].plot(np.linspace(0,x_max,N_Group+1),w_mean_cum_count, label=f"weighted Running Mean CDF at time {T: .1f}")
                ax[1].plot(np.linspace(Time_scale,Time_scale*len(w_Dist_mean_list),len(w_Dist_mean_list)),w_Dist_mean_list, label="Distance to Weighted Running Mean")
            ax[0].legend(loc=0, prop={'size': 8})
            ax[1].legend(loc=0, prop={'size': 8})
            ax[2].legend(loc=1, prop={'size': 8})
            #print(f"Maximum count = {np.max(new_cum_count)}")
            #print(f"new_cum_count: {new_cum_count}" )
            #print(f"run_cum_count: {run_cum_count}" )
            #print(Dist_mean_list)
            if save_fig == True:
                fig_name = str(ex_max) +'_'+ str(ex_slope) +'_'+  str(frag_speed)+'_'+  str(N_split) +'_'+  str(a) +'_'+ str(b) +'_'+  str(alpha_f) +'_'+ str( Nimm) +'_'+  str(imm_rate) +'_'+ str(Time_scale) +'_'+ str(N_Group) +'_'+  str(KS_threshold) +'_'+  str(Step_Max) +'_' + str(int(T))
                fig.savefig(path + '/'+ fig_name +'.png', transparent=False) 
            fig.tight_layout()
            display(fig)
            clear_output(wait = True)
            plt.pause(0.1)
        
        if weighted == True:
           KS_dist_mean = w_KS_dist_mean




        #KS_ratio = abs(KS_dist_mean/Dist_mean_list[0])

          
      
    ################## Recording Mixing Time, Stationary Count, Stationary Density and Age of Fragmentation #################################
      
    Mixing_Time = T
    station_count = mean_cum_count
    station_dens = np.diff(mean_cum_count)
    if weighted == True:
        station_count = w_mean_cum_count
        station_dens = np.diff(w_mean_cum_count)
    ################### Printing Results #################################################
    if Time_steps < Step_Max:
      Equilibrium = True             #Equilibrium is a Boolean variable that indicates if the process reaches equilibrium before the maximum Step size
    else:
      Equilibrium = False
    np.set_printoptions(precision=2)
    print(f"Equilibrium Reached = {Equilibrium}, Mixing Time is {Mixing_Time: .2e}, Oldest Fragment immigrated at {np.min(x.Frag_labs): .2e}")
    #print(f"Mixing Time= {Mixing_Time:.2f}, #Events={N_Event}")
    #print(f"Mean_cum_count: {mean_cum_count},\n New_cum_count: {new_cum_count}")
    #print(f"Mean difference is {d:.3f}")
    #print(f"Diff to Mean: {Diff_to_mean}")
    #print(f"Cumulative count difference ratio is {KS_dist}")
    #print(f"Dist to Run Mean: {KS_dist_mean: .4f}")
    #print(f"KS Dist to Mean History: {np.asarray(Dist_mean_list)}")
    #print(f"Fragment Age = {x.age(Mixing_Time)}")
    #clear_output(wait = True)

    return Equilibrium, Mixing_Time, station_count, x




## Simulating Multiple Samples with same Parameters

#To efficiently simulate multiple samples of FragExImm processes with the same parameters, here we define a class called Simulation. The parameters of the simulations are those needed to run one simulation.
#The function run_simulation with parameter N_sample efficiently runs the N_sample many simulations of the same input parameters. It will also aggregrate the observed count and produce the sample mean stationary distribution over all the samples with the same FragExIm parameters. 

class Simulation:
    
    def __init__(self, ex_max, ex_slope, frag_speed, N_split,  a,b, alpha_f, Nimm, imm_rate, Time_scale, N_Group, KS_threshold, Step_Max, lin_ex=True, B_ex = 1, alpha_e = 0, weighted=False, plot=False, cfDNA= False,  x_max=1, xlim_max=1, peak=True, loglog=True, norm_imm= False, exp_imm=False, mu_i=0):
        self.ex_max=ex_max
        self.ex_slope = ex_slope
        self.frag_speed = frag_speed
        self.N_split = N_split
        self.Nimm = Nimm
        self.imm_rate = imm_rate
        self.Time_scale = Time_scale
        self.N_Group = N_Group
        self.KS_threshold = KS_threshold
        self.Step_Max = Step_Max
        self.weighted = weighted
        self.plot = plot
        self.cfDNA = cfDNA
        self.a = a
        self.b = b
        self.alpha_f = alpha_f
        self.lin_ex = lin_ex
        self.B_ex = B_ex
        self.alpha_e = alpha_e
        self.x_max = x_max
        self.xlim_max = xlim_max
        self.peak = peak
        self.loglog = loglog
        self.norm_imm=norm_imm
        self.exp_imm = exp_imm
        self.mu_i=mu_i


    def __str__(self):
        information = f"Number of Immigrants = {self.Nimm}, Time_Scale = {self.Time_scale}"
        return information

    def run_simulation(self, N_sample,plot=False):
        "Produce N_sample many samples with the parameters given by self."
        print(self)
        i = 1
        Equilibrium_list  = []
        Mix_time_list = []
        Mean_age_list = []
        N_Group = self.N_Group
        imm_rate = self.imm_rate
        frag_speed = self.frag_speed
        ex_slope = self.ex_slope
        ex_max = self.ex_max
        x_max= self.x_max
        xlim_max = self.xlim_max
        a = self.a
        b = self.b
        alpha_f = self.alpha_f
        alpha_e = self.alpha_e
        lin_ex = self.lin_ex
        B_ex = self.B_ex
        loglog = self.loglog
        if lin_ex == True:
            B_ex = ex_max/ex_slope
        Station_agg = np.zeros(N_Group+1)
        Prior_Count = 1
        sim_mat = np.zeros((N_sample, N_Group))
      ################ Run N_sample many experiments with a while loop ###########################
        while i < N_sample+1:
            np.random.seed(i)
            print(f"Experiment: {i}, Seed: {i}")
            initial_condition = Frag_List(np.ones(Prior_Count)*x_max,np.zeros(Prior_Count))
            Equilibrium, Mixing_Time, station_count, x = stationary(initial_condition,self.ex_max, self.ex_slope, self.frag_speed, self.N_split, self.a, self.b, self.alpha_f, self.Nimm, self.imm_rate, self.Time_scale, self.N_Group, self.KS_threshold, self.Step_Max, lin_ex, B_ex , alpha_e , self.weighted, self.plot, self.cfDNA, self.x_max, self.xlim_max, norm_imm=self.norm_imm, exp_imm = self.exp_imm, mu_i=self.mu_i)
            mean_age=np.mean(x.age(Mixing_Time))
            Equilibrium_list.append(Equilibrium)
            Mix_time_list.append(Mixing_Time)
            sample_count = x.cum_dist(N_Group, x_max, count=True) 
            Station_agg += sample_count #Aggragating the stationary distribution across all the experiments
            sim_mat[i-1,:] = np.diff(sample_count) #Storing the sample count in the N-th row
            #print(Station_agg)
            Mean_age_list.append(mean_age)
            i+=1
        Station_mean = Station_agg/N_sample #Averaging over number of samples to produce the sample mean stationary distribution
        sample_mean_age = np.mean(Mean_age_list)
        Station_dens = np.diff(Station_mean)
        beta = 0
        if loglog and B_ex < x_max:
            log_den, log_x_axis, beta, c = tail_loglog_fit(Station_dens,x_max, B_ex, peak = self.peak)
        gamma = 0

        #if a == 1 and b == 1 and lin_ex==True:
            #exit_x_axis, ex_log_den, ex_log_x_axis, gamma, d = ex_reg_loglog_fit(Station_dens,x_max, B_ex, frag_speed, ex_slope)

        up_bound = np.percentile(sim_mat, 95, axis=0)
        low_bound = np.percentile(sim_mat, 5, axis=0)
      ##################### Plotting Sample Mean Stationary Distribution #################################################
        if plot == True:
            N_plots = 2
            #if a == 1 and b == 1 and alpha_f == 1 and lin_ex == True: #plotting exit region when fragmentation is uniform and exit mechanism is linear
                #N_plots = 3
            fig2,ax2 = plt.subplots(N_plots,figsize=(10, 10))
            #if a == 1 and b == 1 and alpha_f == 1 and lin_ex == True: #plotting exit region when fragmentation is uniform 
                #ax2[2].plot(np.log(exit_x_axis), ex_log_den, label="log  density")
                #ax2[2].plot(np.log(exit_x_axis), gamma*np.log((frag_speed/ex_slope-1)*exit_x_axis +B_ex) + d, label = f"best fit line: {gamma: .3}log({B_ex : .3f} + {frag_speed/ex_slope-1 :.3f}x) + {d: .3f}" )
                #ax2[2].set_title("exit region density log-log plot")
                #ax2[2].legend()
                #ax2[0].plot(exit_x_axis,np.exp(d)*(B_ex + (frag_speed/ex_slope-1)*exit_x_axis)**gamma, label = f"Exit Region Best Fit Curve,  {np.exp(d) :.3g}*({B_ex : .3f} + {frag_speed/ex_slope-1 :.3f}x)^{gamma :.3} ")
            ax2[0].plot(np.linspace(0,x_max,N_Group+1)[:-1],Station_dens, label = f"Running Mean Density ")
            ax2[0].plot(np.linspace(0,x_max,N_Group+1)[:-1],up_bound, linestyle='-.', color='darkorchid', label = f"95-th %")
            ax2[0].plot(np.linspace(0,x_max,N_Group+1)[:-1],low_bound, linestyle='-.', color='darkorchid', label = f"5-th %")
            ax2[0].fill_between(np.linspace(0,x_max,N_Group+1)[:-1], low_bound, up_bound, alpha=0.1, color = 'darkorchid')
            if loglog and B_ex < x_max:
                ax2[0].plot(np.exp(log_x_axis),np.exp(beta*log_x_axis + c), label = f"Tail Best Fit Curve, {np.exp(c): .3g}*x^{beta :.3f} ")
                ax2[1].plot(np.exp(log_x_axis),np.exp(beta*log_x_axis + c), label = f"Tail Best Fit Curve, {np.exp(c): .3g}*x^{beta :.3f} ")
            ax2[0].set_xlim([0,xlim_max])
            ax2[0].legend(loc=0, prop={'size': 8})
            ax2[0].set_title("Density Plot")
            ax2[0].set_xlabel('Fragment Length')
            ax2[0].set_ylabel('Fragment Count')
            ax2[1].plot(np.linspace(0,x_max,N_Group+1)[1:-1],Station_dens[1:], label = f"Running Mean Density ")
            ax2[1].plot(np.linspace(0,x_max,N_Group+1)[1:-1],up_bound[1:], linestyle='-.', color='darkorchid', label = f"95-th %")
            ax2[1].plot(np.linspace(0,x_max,N_Group+1)[1:-1],low_bound[1:], linestyle='-.', color='darkorchid', label = f"5-th %")
            ax2[1].fill_between(np.linspace(0,x_max,N_Group+1)[1:-1], low_bound[1:], up_bound[1:], alpha=0.1, color = 'darkorchid')
            #ax2[1].set_xlim([0,xlim_max])
            ax2[1].legend(loc=0, prop={'size': 8})
            ax2[1].set_title("Log -Log Density Plot")
            ax2[1].set_xlabel('Log Fragment Length')
            ax2[1].set_ylabel('Log Fragment Count')
            ax2[1].set_xscale('log')
            ax2[1].set_yscale('log')
            if lin_ex == True:
                fig2.suptitle(f"Exit Function: {self.ex_max: .3f} + {self.ex_slope: .3f} x, Imm={self.imm_rate: .3f}, B_ex: {B_ex: .3f}, Frag_Speed: {self.frag_speed}, DNA length: {x_max : .3f}")
            elif lin_ex == False:
                fig2.suptitle(f"Exit Function: x^{self.alpha_e: .3f} - {B_ex : .3f}^{self.alpha_e: .3f}, Imm={self.imm_rate: .3f}, B_ex: {B_ex: .3f}, Frag_Speed: {self.frag_speed}, Frag_index: {alpha_f: .3f}, a: {a: .3f}, b: {b: .3f}, DNA length: {x_max : .3f}")
                
            fig2.tight_layout()
            #display(fig2)
            #clear_output(wait = True)
        
        mode = np.amax(Station_dens) #The maximum count
        return sim_mat, Station_mean, Equilibrium_list, Mix_time_list, sample_mean_age, mode, beta, gamma


