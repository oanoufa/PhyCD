import numpy as np
from tqdm import trange

class Cont_Detec_HMM:
    """HMM model for the detection of contamination in a genomic sequence.
    This model takes as input a sequence of depth values categorized.
    It is an assembly of several HMMs:
    - 2 states for the first HMM (No contaminant - No dropout / Dropout)
    - 4 states for each additional HMM, describing contamination with a different proportion,
      and the following states:
      With the first letter being the original sample and the second being the contaminant:
      NN / ND / DN / DD (N = No dropout, D = Dropout)
    
    The starting probabilities pi are constrained such that the model can only start in state NN of one of the HMMs.
    The transition probabilities A are constrained such that the model cannot transition from one HMM to another.
    The emission probabilities B are constrained depending on the state:
    For instance, state NN can almost only emit high depth values, while state DD can only emit very low depth values.
    States ND and DN's emission probabilities depend on the value of the contaminant proportion.
    
    
    Categories of depth values are the following:
    With median depth = 1000,
    0: 0 - 5
    1: 5 - 20
    2: 20 - 50
    3: 50 - 100
    4: 100 - 300
    5: 300+
    
    Input:
    - n_models: number of HMMs to assemble
    - cont_prop: list of ranges of proportions of contamination for each HMM
    - n_observations: number of categories for the depth values
    - random_state: random seed for reproducibility
    
    Usage:
    model = Cont_Detec_HMM(n_models=2, cont_prop=[0, 0.2, 0.5], n_observations=6, random_state=42)
    model.initialize_params(sequences)
    model.train(sequences, n_iter=10, tol=1e-4, verbose=True)
    model.viterbi(obs_seq)
    """
    def __init__(self, n_models=1, cont_prop=[0, 0.2, 0.5], n_observations=6, random_state=None):
        self.n_models = n_models
        self.n_states = 2 + n_models*4
        self.n_observations = n_observations
        self.cont_prop = cont_prop
        self.random_state = np.random.RandomState(random_state)
    
    def depth_to_category(self, depth, median_depth):
        """Convert depth values to categories.
        """
        if depth < 0.005 * median_depth:
            return 0
        elif depth < 0.02* median_depth:
            return 1
        elif depth < 0.05 * median_depth:
            return 2
        elif depth < 0.1 * median_depth:
            return 3
        elif depth < 0.3 * median_depth:
            return 4
        else:
            return 5

    def initialize_params(self, sequences):
        """Initialize the parameters of the HMM.
        """
        # Init pi for all states
        self.pi = self.random_state.dirichlet(np.ones(self.n_states))
        # Model shouldn't start in ND or DN states
        self.pi_mask = np.ones(self.n_states)
        for i in range(self.n_models):
            self.pi_mask[2 + i*4 + 1] = 0
            self.pi_mask[2 + i*4 + 2] = 0
        self.pi = self.pi_mask * self.pi
        # Normalize pi
        self.pi = self.pi / np.sum(self.pi)
        
        # Init A: Transition probabilities between states
        self.A = self.random_state.dirichlet(np.ones(self.n_states), size=self.n_states)
        self.A_min = np.zeros((self.n_states, self.n_states))
        self.A_max = np.zeros((self.n_states, self.n_states))
        # No transition between HMMs
        
        # Transition between states of the 2 states HMM
        self.A_max[:2, :2] = 1
        
        # Transition between states of the 4 states HMMs
        for i in range(self.n_models):
            self.A_max[2 + i*4:2 + (i+1)*4, 2 + i*4:2 + (i+1)*4] = 1     
            
        # Clip A and normalize
        self.A = np.clip(self.A, self.A_min, self.A_max)
        self.A = self.A / self.A.sum(axis=1, keepdims=True)
        
        # Init B: Emission probabilities
        self.B = self.random_state.dirichlet(np.ones(self.n_observations), size=self.n_states)
        self.B_min = np.zeros((self.n_states, self.n_observations))
        self.B_max = np.zeros((self.n_states, self.n_observations))
        
        # States NN (0, 2, 6 etc) have emission probabilities for categories 0, 1, 2, 3 set to 0
        # States ND and DN (3, 4, 7, 8 etc) should not be hardly constrained but clipped using cont_prop
        # States DD (1, 5, 9 etc) have emission probabilities for categories 4, 5 set to 0
        self.B_max[0, 4:] = 1
        self.B_max[1, :4] = 1
        
        median_depth = np.median(sequences)
        
        for i in range(self.n_models):
            # States NN
            self.B_max[2 + i*4, 4:] = 1
            # States DD
            self.B_max[2 + i*4 + 3, :4] = 1
            
            # cont_prop_min = self.cont_prop[i]
            # min_cont_depth_cat = self.depth_to_category(median_depth * cont_prop_min, median_depth)
            # min_main_depth_cat = self.depth_to_category(median_depth * (1 - cont_prop_min), median_depth)
            # cont_prop_max = self.cont_prop[i+1]
            # max_cont_depth_cat = self.depth_to_category(median_depth * cont_prop_max, median_depth)
            # max_main_depth_cat = self.depth_to_category(median_depth * (1 - cont_prop_max), median_depth)
            
            # print("Min cat main dropout", min_cont_depth_cat,
            #         "Max cat main dropout", max_cont_depth_cat,
            #         "Min cat cont dropout", min_main_depth_cat,
            #         "Max cat cont dropout", max_main_depth_cat,
            #         "cont_prop_min", cont_prop_min,
            #         "cont_prop_max", cont_prop_max)
            
            # # Contaminant dropout (3, 7 etc)
            # # Depth should be between cont_prop_min*typical_depth and cont_prop_max*typical_depth
            # self.B_max[2 + i*4 + 1, min_main_depth_cat:max_main_depth_cat+1] = 1
            # # Main sample dropout (4, 8 etc)
            # self.B_max[2 + i*4 + 2, min_cont_depth_cat:max_cont_depth_cat+1] = 1
            
            self.B_max[2 + i*4 + 1, :] = 1
            self.B_max[2 + i*4 + 2, :] = 1
            
        # Clip B and normalize
        self.B = np.clip(self.B, self.B_min, self.B_max)
        self.B = self.B / self.B.sum(axis=1, keepdims=True)

    def forward(self, obs_seq):
        T = len(obs_seq)
        # Separate alpha for each HMM
        alpha = np.zeros((T, self.n_states))
    
    def backward(self, obs_seq):
        T = len(obs_seq)
        beta = np.zeros((T, self.n_states))
        beta[-1] = 1
        for t in range(T - 2, -1, -1):
            for j in range(self.n_states):
                beta[t, j] = np.sum(self.A[j] * self.B[:, obs_seq[t + 1]] * beta[t + 1])
        return beta
    
    def baum_welch_step(self, sequences):
        for obs_seq in sequences:
            T = len(obs_seq)
            alpha = self.forward(obs_seq)
            beta = self.backward(obs_seq)
            
            # Compute the gamma and xi values
            gamma = np.zeros((T, self.n_states))
            xi = np.zeros((T - 1, self.n_states, self.n_states))
            
            for t in range(T):
                gamma[t] = alpha[t] * beta[t] / np.sum(alpha[t] * beta[t])
            
            for t in range(T - 1):
                xi[t] = (alpha[t][:, np.newaxis] * self.A * self.B[:, obs_seq[t + 1]] * beta[t + 1]) / np.sum(alpha[t] * beta[t])
            
        return gamma, xi
    
    def update_params(self, obs_seq, gamma, xi):
        """Update the parameters of the HMM using the gamma and xi values.
        """
        # Update pi
        self.pi = gamma[0]
        self.pi = self.pi_mask * self.pi
        self.pi = self.pi / np.sum(self.pi)
        
        # Update A
        self.A = np.sum(xi, axis=0) / np.sum(gamma[:-1], axis=0)[:, np.newaxis]
        
        # Update B
        for j in range(self.n_states):
            for k in range(self.n_observations):
                self.B[j, k] = np.sum(gamma[obs_seq == k, j]) / np.sum(gamma[:, j])
        
        # Clip and normalize
        self.A = np.clip(self.A, self.A_min, self.A_max)
        self.A = self.A / self.A.sum(axis=1, keepdims=True)
        
        self.B = np.clip(self.B, self.B_min, self.B_max)
        self.B = self.B / self.B.sum(axis=1, keepdims=True)
        
    def train(self, obs_seqs, n_iter=10, tol=1e-4, verbose=False):
        """Train the HMM using the Baum-Welch algorithm.

        Args:
            obs_seqs (list or np.ndarray): List of observation sequences or a single sequence.
            n_iter (int, optional): Number of iterations for the Baum-Welch algorithm. Defaults to 10.
            tol (float, optional): Tolerance for convergence. Defaults to 1e-4.
            verbose (bool, optional): Whether to print progress. Defaults to False.

        Returns:
            list: List of log-likelihoods for each iteration.
        """
        if isinstance(obs_seqs, np.ndarray):
            obs_seqs = [obs_seqs]
        
        log_likelihoods = []
        
        for i in trange(n_iter, desc="Training", disable=not verbose):
            total_log_likelihood = 0
            for obs_seq in obs_seqs:
                gamma, xi = self.baum_welch_step(obs_seq)
                self.update_params(obs_seq, gamma, xi)
                total_log_likelihood += np.sum(np.log(np.sum(self.forward(obs_seq), axis=1)))
            
            log_likelihoods.append(total_log_likelihood)
            
            if i > 0 and abs(log_likelihoods[-1] - log_likelihoods[-2]) < tol:
                break
        
        return log_likelihoods

    def viterbi(self, obs_seq):
        T = len(obs_seq)
        delta = np.zeros((T, self.n_states))
        psi = np.zeros((T, self.n_states), dtype=int)
        delta[0] = self.pi * self.B[:, obs_seq[0]]
        for t in range(1, T):
            for j in range(self.n_states):
                prob = delta[t - 1] * self.A[:, j]
                psi[t, j] = np.argmax(prob)
                delta[t, j] = np.max(prob) * self.B[j, obs_seq[t]]
        path = np.zeros(T, dtype=int)
        path[-1] = np.argmax(delta[-1])
        for t in range(T - 2, -1, -1):
            path[t] = psi[t + 1, path[t + 1]]
        return path
        
