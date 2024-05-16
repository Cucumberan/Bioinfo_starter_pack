import numpy as np
from sklearn.base import BaseEstimator
from sklearn.tree import (DecisionTreeClassifier)
import random
from concurrent.futures import ProcessPoolExecutor

SEED = 111
random.seed(SEED)
np.random.seed(SEED)


def train_tree(args):
    X, y, n_features, n_samples, max_features, max_depth, seed = args
    np.random.seed(seed)
    features_choice = np.random.choice(n_features, max_features, replace=False)
    sample_choice = np.random.choice(n_samples, n_samples, replace=True)
    X_sampled = X[sample_choice][:, features_choice]
    y_sampled = y[sample_choice]
    tree = DecisionTreeClassifier(max_depth=max_depth, random_state=seed)
    tree.fit(X_sampled, y_sampled)
    return tree, features_choice

def get_tree_probas(tree_feat_ids, X):
    tree, feat_ids = tree_feat_ids
    return tree.predict_proba(X[:, feat_ids])

class RandomForestClassifierCustom(BaseEstimator):
    """
    A custom implementation of a Random Forest classifier that supports parallel processing
    for both training and predicting. This implementation leverages decision trees as the
    base estimators.

    Parameters:
        n_estimators (int): Number of trees in the forest. Defaults to 10.
        max_depth (int, optional): The maximum depth of the trees. If None, the nodes are expanded
                                   until all leaves are pure or until all leaves contain less than
                                   min_samples_split samples. Defaults to None.
        max_features (int, optional): The number of features to consider when looking for the best split.
                                      If None, all features are considered. Defaults to None.
        random_state (int): A seed used by the random number generator for reproducibility.
        n_jobs (int): The number of jobs to run in parallel for both `fit` and `predict`.
                      `n_jobs=1` means using one processor, `n_jobs=-1` means using all processors.

    Attributes:
        trees (list): A list of the fitted decision tree classifiers.
        feat_ids_by_tree (list): A list where each element is an array of indices of the features
                                 used by the corresponding tree in the forest.

    Methods:
        fit(X, y, n_jobs=1): Fits the random forest model to the training data.
        predict_proba(X, n_jobs=1): Predicts class probabilities for X using the trained model.
        predict(X, n_jobs=1): Predicts the class labels for X using the trained model.

    Example:
        >>> from sklearn.datasets import make_classification
        >>> X, y = make_classification(n_samples=100, n_features=4, random_state=42)
        >>> clf = RandomForestClassifierCustom(n_estimators=5, max_depth=10, random_state=42)
        >>> clf.fit(X, y)
        >>> print(clf.predict(X))

    Note:
        This implementation requires the function `train_tree` and `get_tree_probas` to be defined
        in the same module as helper functions, which are responsible for training each tree and 
        obtaining probabilities, respectively. These functions should be defined at the module level
        to be pickleable by the multiprocessing module.
    """
    def __init__(self, n_estimators=10, max_depth=None, max_features=None, random_state=SEED, n_jobs=1):
        self.n_estimators = n_estimators
        self.max_depth = max_depth
        self.max_features = max_features
        self.random_state = random_state
        self.n_jobs = n_jobs
        self.trees = []
        self.feat_ids_by_tree = []

    def fit(self, X, y, n_jobs=1):
        self.classes_ = sorted(np.unique(y))
        n_features = X.shape[1]
        n_samples = X.shape[0]

        args = [(X, y, n_features, n_samples, self.max_features if self.max_features is not None else n_features, self.max_depth, self.random_state + i) for i in range(self.n_estimators)]

        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            results = list(executor.map(train_tree, args))

        self.trees, self.feat_ids_by_tree = zip(*results)
        return self
    
    def predict_proba(self, X, n_jobs=1):
        # Prepare the input for parallel execution
        with ProcessPoolExecutor(max_workers=n_jobs) as executor:
            probas = list(executor.map(get_tree_probas, zip(self.trees, self.feat_ids_by_tree), [X]*len(self.trees)))
        
        # Aggregate the probabilities from each tree
        avg_probas = np.mean(probas, axis=0)
        return avg_probas
    
    def predict(self, X, n_jobs=1):
        probas = self.predict_proba(X, n_jobs=n_jobs)
        predictions = np.argmax(probas, axis=1)
        return predictions