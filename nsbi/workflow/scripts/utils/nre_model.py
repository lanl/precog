"""
Neural Ratio Estimation (NRE) model architecture.

Implements a binary classifier neural network that learns to distinguish
joint distribution samples (true parameter-statistic pairs) from marginal
distribution samples (randomly paired parameters and statistics).

This version supports optional masking for variable-length time series inputs.

Masking Strategy:
-----------------
When use_masking=True, masks are CONCATENATED as additional features rather than
multiplied with the input. This allows the network to learn which positions are
valid/invalid while preserving gradient flow and batch normalization stability.

Input without masking: [features(n), parameters(m)]  -> dim = n+m
Input with masking:    [features(n), mask(n), parameters(m)] -> dim = 2n+m
"""

import torch
import torch.nn as nn

class NREClassifier(nn.Module):
    """
    Neural network for ratio estimation with optional mask support.
    
    Takes concatenated summary statistics and parameters as input,
    outputs probability that the pair comes from the joint distribution.
    Supports masking for variable-length time series inputs.
    """
    
    def __init__(self, input_dim=5, hidden_dims=[64, 32, 16], 
                 dropout_rate=0.2, use_batch_norm=True, use_masking=False,
                 n_features=None):
        """
        Initialize NRE classifier.
        
        Args:
            input_dim: Input dimension (n_features + n_parameters)
                       When use_masking=True, this should be (2*n_features + n_parameters)
                       because masks are concatenated as additional features
            hidden_dims: List of hidden layer sizes
            dropout_rate: Dropout probability for regularization
            use_batch_norm: Whether to use batch normalization
            use_masking: Whether to concatenate masks as additional features for 
                         variable-length/missing data scenarios
            n_features: Number of feature dimensions (required if use_masking=True)
                        Used to know how many mask dimensions to expect
        """
        super(NREClassifier, self).__init__()
        
        self.input_dim = input_dim
        self.hidden_dims = hidden_dims
        self.dropout_rate = dropout_rate
        self.use_batch_norm = use_batch_norm
        self.use_masking = use_masking
        self.n_features = n_features
        
        # Validate configuration
        if use_masking and n_features is None:
            raise ValueError("n_features must be provided when use_masking=True")
        
        # Build network layers
        layers = []
        prev_dim = input_dim
        
        for hidden_dim in hidden_dims:
            # Linear layer
            layers.append(nn.Linear(prev_dim, hidden_dim))
            
            # Batch normalization
            if use_batch_norm:
                layers.append(nn.BatchNorm1d(hidden_dim))
            
            # Activation
            layers.append(nn.ReLU())
            
            # Dropout
            if dropout_rate > 0:
                layers.append(nn.Dropout(dropout_rate))
            
            prev_dim = hidden_dim
        
        # Output layer (single logit)
        layers.append(nn.Linear(prev_dim, 1))
        
        self.network = nn.Sequential(*layers)
    
    def forward(self, x, mask=None):
        """
        Forward pass - simplified, masks pre-concatenated in input.
        
        Args:
            x: Input tensor of shape (batch_size, input_dim)
               When masking enabled: [features, masks, parameters]
               When masking disabled: [features, parameters]
            mask: Deprecated parameter, ignored for backward compatibility
        
        Returns:
            Logits of shape (batch_size, 1)
        
        Note:
        -----
        Masks should be pre-concatenated into the input x before calling forward().
        The mask argument is kept for backward compatibility but is not used.
        """
        # Simply pass through - input already has masks concatenated if needed
        return self.network(x)
    
    def predict_proba(self, x, mask=None):
        """
        Predict probability (apply sigmoid to logits).
        
        Args:
            x: Input tensor of shape (batch_size, input_dim)
            mask: Optional mask tensor
        
        Returns:
            Probabilities of shape (batch_size, 1)
        """
        with torch.no_grad():
            logits = self.forward(x, mask)
            probs = torch.sigmoid(logits)
        return probs

