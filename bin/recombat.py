# Author: Michael F. Adamer (https://github.com/BorgwardtLab/reComBat)
# Slightly modified by Olivier Coen (only cosmetic changes + typing)
# Only the parametric part was kept, without X nor C

import logging
import numpy as np
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet

logger = logging.getLogger(__name__)


class ReComBat:
    """
    ReComBat class

    Parameters
    ----------
    model : str
        Choose a linear model, ridge, Lasso, elastic_net.
        The default is 'linear'.
    config : dict
        A dictionary containing kwargs for the model (see sklean.linear_model for details).
        The default is None.
    conv_criterion : float, optional
        The convergence criterion for the optimization.
        The default is 1e-4.
    max_iter : int, optional
        The maximum number of steps of the parametric empirical Bayes optimization.
        The detault is 1000.
    """

    model: str
    config: dict | None
    conv_criterion: float
    max_iter: int
    is_fitted: bool
    gamma_star_hat_: np.ndarray
    delta_star_squared_hat_: np.ndarray
    alpha_: np.ndarray
    beta_x_: np.ndarray
    beta_c_: np.ndarray
    sigma_: np.ndarray

    def __init__(
        self,
        parametric: bool = True,
        model: str = 'linear',
        config: dict | None = None,
        conv_criterion: float = 1e-4,
        max_iter: int = 1000,
    ):
        self.model = model

        if config is not None:
            self.config = config
        else:
            self.config = {}

        self.conv_criterion = conv_criterion
        self.max_iter = max_iter

        # Set to True to indicate that reComBat has been fitted.
        self.is_fitted = False

        self.gamma_star_hat_ = None
        self.delta_star_squared_hat_ = None
        self.alpha_ = None
        self.beta_x_ = None
        self.beta_c_ = None
        self.sigma_ = None


    def fit_transform(self, data: np.ndarray, batches: np.ndarray) -> np.ndarray:
        '''
        Fit and transform in one go.
        '''
        self.fit(data, batches)
        return self.transform(data, batches)


    @staticmethod
    def get_dummies(batches: np.ndarray) -> np.ndarray:
        categories, codes = np.unique(batches.astype(str), return_inverse=True)
        return np.eye(len(categories), dtype=int)[codes]


    # -------------------------------------------------------
    # FIT
    # -------------------------------------------------------

    def fit(self, data: np.ndarray, batches: np.ndarray) -> None:
        """
        Fit method.

        Parameters
        ----------
        data : np.ndarray
            A numpy array containing the data matrix.
            The format is (rows x columns) = (samples x features)
        batches : pandas series
            A pandas series containing the batch of each sample in the dataframe.

        Returns
        -------
        None
        """
        logger.info("Starting to fit reComBat.")
        if np.isnan(data).any().any():
            raise ValueError("The data contains NaN values.")

        if len(np.unique(batches)) == 1:
            raise ValueError("There should be at least two batches in the dataset.")

        batches_one_hot = self.get_dummies(batches)

        if np.any(batches_one_hot.sum(axis=0) == 1):
            raise ValueError("There should be at least two values for each batch.")

        if data.shape[0] != batches.shape[0]:
            raise ValueError("The batch and data matrix have a different number of samples.")

        logger.info("Fit the linear model.")
        Z = self.fit_model_(data, batches_one_hot)

        logger.info("Starting the empirical parametric optimisation.")
        self.parametric_optimization_(Z, batches_one_hot)
        logger.info("Optimisation finished.")

        self.is_fitted = True
        logger.info("reComBat is fitted.")


    def fit_model_(self, data: np.ndarray, batches_one_hot: np.ndarray) -> np.ndarray:
        '''
        Fit the linear model.
        '''

        # Create the design matrix
        num_batches = batches_one_hot.shape[1]
        Covariates = batches_one_hot

        # Initialise the model class
        if self.model == 'linear':
            model = LinearRegression(fit_intercept=False, **self.config)
        elif self.model == 'ridge':
            model = Ridge(fit_intercept=False, **self.config)
        elif self.model == 'Lasso':
            model = Lasso(fit_intercept=False, **self.config)
        elif self.model == 'elastic_net':
            model = ElasticNet(fit_intercept=False, **self.config)
        else:
            raise ValueError('Model not implemented')

        model.fit(Covariates, data)

        # Save the fitted parameters
        # Note that alpha is computed implicitly via the contraints on the batch parameters.
        self.alpha_ = np.matmul(
            batches_one_hot.sum(axis=0, keepdims=True) / batches_one_hot.sum(),
            model.coef_.T[:num_batches]
        )

        self.beta_x_ = model.coef_.T[num_batches: num_batches]
        self.beta_c_ = model.coef_.T[num_batches:]

        # Compute the standard deviation of the reconstructed data.
        data_hat = np.matmul(Covariates, model.coef_.T)
        self.sigma_ = np.mean(
            (data - data_hat)**2,
            axis=0,
            keepdims=True
        )

        # Standardise the data.
        a = np.matmul(Covariates[:, num_batches: num_batches], self.beta_x_)
        b = np.matmul(Covariates[:, num_batches:], self.beta_c_)
        return (data.copy() - a - b - self.alpha_) / np.sqrt(self.sigma_)


    def parametric_optimization_(self, Z: np.ndarray, batches_one_hot: np.ndarray):
        '''
        Perform parametric optimization.
        '''
        gamma_hat, gamma_bar, delta_hat_squared, tau_bar_squared, lambda_bar, theta_bar = self.compute_init_values_parametric_(Z, batches_one_hot)
        n = batches_one_hot.sum(axis=0)

        gamma_star = gamma_hat.copy()
        delta_star_squared = delta_hat_squared.copy()

        for i in range(batches_one_hot.shape[1]):
            gamma_star[i], delta_star_squared[i] = self.parametric_update_(gamma_star[i],
                                                                    delta_star_squared[i],
                                                                    n[i],
                                                                    tau_bar_squared[i],
                                                                    gamma_hat[i],
                                                                    gamma_bar[i],
                                                                    theta_bar[i],
                                                                    lambda_bar[i],
                                                                    Z[batches_one_hot[:,i] == 1])

        self.gamma_star_hat_ = gamma_star
        self.delta_star_squared_hat_ = delta_star_squared


    def compute_init_values_parametric_(
        self,
        Z: np.ndarray,
        batches_one_hot: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        '''
        Compute the starting values of the Bayesian optimization.
        '''
        gamma_hat = np.array([np.mean(Z[batches_one_hot[:,i] == 1], axis=0) for i in range(batches_one_hot.shape[1])])
        gamma_bar = np.mean(gamma_hat, axis=1)
        tau_bar_squared = np.var(gamma_hat, axis=1, ddof=1)
        delta_hat_squared = np.array([np.var(Z[batches_one_hot[:,i] == 1], axis=0, ddof=1) for i in range(batches_one_hot.shape[1])])
        V_bar = np.mean(delta_hat_squared, axis=1)
        S_bar_squared = np.var(delta_hat_squared, axis=1, ddof=1)

        lambda_bar = (V_bar**2 + 2 * S_bar_squared) / S_bar_squared
        theta_bar = (V_bar**3 + V_bar * S_bar_squared) / S_bar_squared

        return gamma_hat, gamma_bar, delta_hat_squared, tau_bar_squared, lambda_bar, theta_bar


    def parametric_update_(
        self,
        gamma_star_i: np.ndarray,
        delta_star_squared_i: np.ndarray,
        n_i: int,
        tau_bar_squared_i: np.ndarray,
        gamma_hat_i: np.ndarray,
        gamma_bar_i: np.ndarray,
        theta_bar_i: np.ndarray,
        lambda_bar_i: np.ndarray,
        Z_i: np.ndarray
    ):
        '''
        Perform the optimization for one batch at a time
        '''
        gamma_star_new_i = gamma_star_i.copy()
        delta_star_squared_new_i = delta_star_squared_i.copy()

        iterations = 0
        convergence = self.conv_criterion + 1
        while (convergence > self.conv_criterion) and (iterations < self.max_iter):
            gamma_star_new_i = self.new_gamma_parametric_(n_i,
                                                          tau_bar_squared_i,
                                                          gamma_hat_i,
                                                          delta_star_squared_i,
                                                          gamma_bar_i)
            delta_star_squared_new_i = self.new_delta_star_squared_parametric_(theta_bar_i,
                                                                               Z_i,
                                                                               gamma_star_new_i,
                                                                               n_i,
                                                                               lambda_bar_i)

            convergence = np.max([
                np.max(np.abs(gamma_star_new_i - gamma_star_i) / gamma_star_i),
                np.max(np.abs(delta_star_squared_new_i - delta_star_squared_i) / delta_star_squared_i)
            ])
            gamma_star_i = gamma_star_new_i
            delta_star_squared_i = delta_star_squared_new_i
            iterations += 1

        if (iterations >= self.max_iter):
            logger.warning("Maximum number of iterations reached")

        return gamma_star_i, delta_star_squared_i


    @staticmethod
    def new_gamma_parametric_(
        n_i: int,
        tau_bar_squared_i: np.ndarray,
        gamma_hat_i: np.ndarray,
        delta_star_squared_i: np.ndarray,
        gamma_bar_i: np.ndarray,
    ) -> np.ndarray:
        numerator = n_i * tau_bar_squared_i * gamma_hat_i + delta_star_squared_i * gamma_bar_i
        denominator = n_i * tau_bar_squared_i + delta_star_squared_i
        return numerator / denominator

    @staticmethod
    def new_delta_star_squared_parametric_(
        theta_bar_i: np.ndarray,
        Z_i: np.ndarray,
        gamma_star_new_i: np.ndarray,
        n_i: int,
        lambda_bar_i: np.ndarray,
    ) -> np.ndarray:
        numerator = (theta_bar_i + 0.5 * np.sum((Z_i - gamma_star_new_i)**2, axis=0))
        denominator = (0.5 * n_i + lambda_bar_i - 1)
        return numerator / denominator


    # -------------------------------------------------------
    # TRANSFORM
    # -------------------------------------------------------

    def transform(self, data: np.ndarray, batches: np.ndarray) -> np.ndarray:
        """
        Transform method.
        ----------------

        Adjusts a dataframe. Please make sure that the number of batches,
        features and design matrix features match.

        Parameters
        ----------
        data : numpy array
            A numpy array containing the data matrix.
            The format is (rows x columns) = (samples x features)
        batches : numpy array
            A numpy array containing the batch of each sample in the dataframe.

        Returns
        -------
        A numpy array of the same shape as the input array.
        """

        logger.info("Starting to transform.")
        batches_one_hot = self.get_dummies(batches)

        if not self.is_fitted:
            raise AttributeError("reComBat has not been fitted yet.")

        if data.shape[1] != self.alpha_.shape[1]:
                raise ValueError("Wrong number of features.")

        if batches_one_hot.shape[1] != self.gamma_star_hat_.shape[0]:
            raise ValueError("Wrong number of batches.")

        Z = (data.copy() - self.alpha_) / np.sqrt(self.sigma_)

        data_adjusted = self.adjust_data_(Z, batches_one_hot)

        logger.info("Transform finished.")
        return data_adjusted


    def adjust_data_(self, Z: np.ndarray, batches_one_hot: np.ndarray):
        '''
        Perform the final adjustment step.
        '''
        tmp = np.zeros_like(Z)
        for i in range(batches_one_hot.shape[1]):
            tmp[batches_one_hot[:, i] == 1] = (Z[batches_one_hot[:, i] == 1] - self.gamma_star_hat_[i]) / np.sqrt(self.delta_star_squared_hat_[i])
        data_adjusted = np.sqrt(self.sigma_) * tmp + self.alpha_
        return data_adjusted
