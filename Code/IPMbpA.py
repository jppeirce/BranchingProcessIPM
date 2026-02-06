###########################################################################################################
# Model description provided in the paper:
# Combining branching processes and Allee effects into an integral projection model to assess invasion risk
# #########################################################################################################
import numpy as np
import scipy.stats as stats
from scipy.special import expit
import pandas as pd
from random import choices, random

class length_weight:
    '''
    Define log-10 length-weight class
    '''
    def __init__(self, intercept = 1.02, slope = 3.02):
        self.intercept = intercept
        self.slope = slope
    
    def predict_weight(self, length):
        weight = 10.0 ** (self.intercept + self.slope * np.log10(length))
        return(weight)
    
class prob_survival:
    '''
    Define surival as a function of length.
    '''
    def __init__(self, lw_func, slope = 2.7, exponent = -0.315):
        self.slope = slope
        self.exponent = exponent,
        self.lw_func = lw_func
    def prob_s(self, length):
        weight = self.lw_func.predict_weight(length) * 1000
        ## 1000 converts from kg to g
        survival = 1.0 - (1.0 - np.exp(- self.slope * weight ** (self.exponent)))
        return survival

class node:
    '''
    This is from the MetaIPM package.
    '''
    def __init__(self, vonB_Linf, vonB_K, vonB_sigma_k):
        self.vonB_Linf = vonB_Linf
        self.vonB_K = vonB_K
        self.vonB_sigma_k = vonB_sigma_k

    def growth(self, length_now, length_next):
        z = np.atleast_1d(length_now)
        z_prime = np.atleast_1d(length_next)

        project = np.zeros((len(z), len(z_prime)))
        for index in range(0, len(z)):
            location_parameter = np.exp(- self.vonB_K) * z[index] + \
                (1-np.exp(-self.vonB_K)) * self.vonB_Linf
            prob_raw = stats.norm.pdf(x=z_prime,
                                      loc=location_parameter,
                                      scale=self.vonB_sigma_k)
            project[index, :] = prob_raw / prob_raw.sum()
        return project
    
class mature():
    '''
    Define maturity as a function of length.
    '''
    def __init__(self, intercept = -7.41, slope = 24.98):
        self.mat_alpha = intercept
        self.mat_beta = slope
    def maturity_prob(self, length_in):
        return expit(self.mat_alpha + self.mat_beta * length_in)

def grow_select(growth_matrix, pop_vec_in):
    """
    Selects the growth matrix columns corresponding to the population vector
    """
    # Extract growth probs
    grow_pop_pos = growth_matrix[np.where(pop_vec_in > 0)[0], :] # selects all growth columns where pop>0
    # Create array to hold outputs
    new_pop_array = np.zeros(grow_pop_pos.shape[1])
    # Loop through growth columns
    for g_idx in range(grow_pop_pos.shape[0]):
        # select new growth bins
        new_pop_array += np.random.multinomial(pop_vec_in[pop_vec_in>0][g_idx], grow_pop_pos[g_idx])
    return new_pop_array

def ipm_bp_allee(
    p_survive,
    p_mature,
    lw_input,
    gamma,
    C,
    K,
    k,
    min_len,
    max_len,
    n_len,
    n_time,
    length_mesh = None,
    growth_matrix = None,
    yoy_mean = 0.319,
    yoy_sd = 0.040,
    init_pop = {"inital population": np.array([0.0,  50.0, 50.0, 50.0]),
                 "length (m)": np.array([0.1, 0.4, 0.6, 0.8])}):
    
    if length_mesh is None:
        length_mesh = np.linspace(min_len, max_len, n_len)

    if growth_matrix is None:
        growth_matrix = von_b_growth.growth(length_now = length_mesh, length_next = length_mesh)
        
    delta_z = np.diff(length_mesh)[0] # width of length bin
    # Convert input lengths to mesh bins
    bins = np.digitize(init_pop["length (m)"], length_mesh, right=True)
    pop_dens = np.zeros([n_len, n_time + 1])
    # for needed if multiple lengths in same bin
    for idx in range(len(bins)):
        pop_dens[bins[idx], 0] += init_pop["initial population"][idx]
    # Normalize to ensure the sum equals the initial population
    pop_dens[:, 0] *= init_pop["initial population"].sum() / (pop_dens[:, 0] * delta_z).sum()

    yoy_place_raw = stats.norm.pdf(x=length_mesh,
                               loc=yoy_mean,
                               scale=yoy_sd)
    yoy_place_den = yoy_place_raw/yoy_place_raw.sum() 
    
    for t_idx in range(n_time):
        # 1. Compute the components of A (the probability of producing k recruits)
        # current biomass
        biomass_t = (lw_input.predict_weight(length_mesh) * pop_dens[:, t_idx] * delta_z).sum()
        # predicted biomass for next year:
        biomass_tp1 =  biomass_t*np.exp(gamma*(1-biomass_t/K)*( (biomass_t-C)))
        pop_mat_dens = np.array([stats.binom.rvs(n = 1,
                            p = p_mature.maturity_prob(xx),
                            loc=0, size = int(pop_dens[xx, t_idx]),
                            random_state=None).sum() for xx in range(n_len)])
        mature_total = (pop_mat_dens*delta_z).sum() #m_t
        # survival = s(z)* n(z,t)
        surv_dens =  np.array([stats.binom.rvs(n = 1,
                     p = p_survive.prob_s(length_mesh[xx]),
                     loc=0, size = int(pop_dens[xx, t_idx]),
                     random_state=None).sum() for xx in range(n_len)])
        # and grow = G(z',z)*s(z)*n(z,t)
        surv_grow_dens = grow_select(growth_matrix = growth_matrix, pop_vec_in = surv_dens)
        survivor_biomass = (lw_input.predict_weight(length_mesh) * surv_grow_dens * delta_z).sum()
        # k C(z'):
        k_density = np.random.multinomial(int(k), yoy_place_den)/delta_z
        # biomass of (k-1) recruits
        k_biomass = (lw_input.predict_weight(length_mesh)*k_density*delta_z).sum()
        # 2. Define A
        if (biomass_tp1 - survivor_biomass) < 0: # current biomass too great, do not reproduce
            A = 0
        elif (biomass_tp1 - survivor_biomass) > (k_biomass * mature_total): # large gap, all reproduce
            A = 1
        elif k_biomass * mature_total == 0:
            A = 0
        else:
            A = (biomass_tp1 - survivor_biomass) / (k_biomass * mature_total)
        # 2. Compute BP probabilities and run update
        for i in range(len(length_mesh)):
            # Calculate the expected number of fish in the length bin
            X_i = pop_dens[i, t_idx] * delta_z
            X_int = int(X_i)  # Integer part (whole fish count)
            X_frac = X_i - X_int  # Fractional part (probabilistic fish)
            
            # Define branching process probabilities
            m = p_mature.maturity_prob(length_mesh[i])  # Probability of maturing at size z
            s = p_survive.prob_s(length_mesh[i])[0]     # Probability of surviving at size z
    
            # Define outcome probabilities for the branching process
            Xi_0 = (1 - m) * (1 - s) + m * (1 - A) * (1 - s)  
            Xi_1 = (1 - m) * s + m * (1 - A) * s
            Xi_k = m * A * (1 - s)
            Xi_kp1 = m * A * s

            # Process each whole fish
            for _ in range(X_int):
                Xi_select = choices(['0', '1', 'k', 'k+1'], [Xi_0, Xi_1, Xi_k, Xi_kp1])[0]

                # Update population density based on branching outcome
                if Xi_select == '1':  # Survives and grows
                    pop_dens[:, t_idx + 1] += growth_matrix[i, :] * (1 / delta_z)
        
                elif Xi_select == 'k':  # Spawns and dies
                    yoy = np.random.multinomial(k, yoy_place_den)
                    pop_dens[:, t_idx + 1] += yoy / delta_z
        
                elif Xi_select == 'k+1':  # Survives, grows, and spawns
                    yoy = np.random.multinomial(k, yoy_place_den)
                    pop_dens[:, t_idx + 1] += growth_matrix[i, :] * (1 / delta_z) + yoy / delta_z
                    
            # Handle probabilistic fish if fractional part exists
            if random() < X_frac:
                Xi_select = choices(['0', '1', 'k', 'k+1'], [Xi_0, Xi_1, Xi_k, Xi_kp1])[0]
        
                if Xi_select == '1':  # Survives and grows
                    pop_dens[:, t_idx + 1] += growth_matrix[i, :] * (1 / delta_z)
        
                elif Xi_select == 'k':  # Spawns and dies
                    yoy = np.random.multinomial(k, yoy_place_den)
                    pop_dens[:, t_idx + 1] += yoy / delta_z
        
                elif Xi_select == 'k+1':  # Survives, grows, and spawns
                    yoy = np.random.multinomial(k, yoy_place_den)
                    pop_dens[:, t_idx + 1] += growth_matrix[i, :] * (1 / delta_z) + yoy / delta_z    

    return pop_dens

def stoch_wrapper(
    p_survive,
    p_mature,
    lw_input,
    gamma,
    C,
    K,
    k,
    min_len,
    max_len,
    n_len,
    n_time,
    n_iter,
    length_mesh,
    growth_matrix,
    yoy_mean,
    yoy_sd,
    init_pop):
    
    total_out = pd.DataFrame()
    for iter_idx in range(n_iter):
        iter_out = ipm_bp_allee(
            p_survive = p_survive,
            p_mature = p_mature,
            lw_input = lw_input,
            gamma = gamma,
            C = C,
            K = K,
            k = k,
            min_len = min_len,
            max_len = max_len,
            n_len = n_len,
            n_time = n_time,
            length_mesh = length_mesh,
            growth_matrix = growth_matrix,
            yoy_mean = yoy_mean,
            yoy_sd = yoy_sd,
            init_pop = init_pop)    
        
        if length_mesh is None:
            length_mesh = np.linspace(min_len, max_len, n_len)
    
        delta_z = np.diff(length_mesh)[0] # width of length bin
        total_pop = pd.DataFrame(iter_out)*delta_z
        total_pop["Length"] = length_mesh
        total_pop_long = total_pop.melt(id_vars = "Length", var_name = "Year", value_name = "Population")
        total_pop_sum = total_pop_long.groupby("Year").aggregate(["sum"]).drop(axis=1, labels="Length")
        total_pop_sum.columns = total_pop_sum.columns.droplevel(1)
        total_pop_sum["iter"] = iter_idx
        total_out = pd.concat([total_out, total_pop_sum], axis = 0)

    return total_out