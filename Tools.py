import sys
import numpy as np
import pandas as pd
# import matplotlib.pyplot as plt
from sklearn.model_selection import LeaveOneOut, KFold, cross_val_score


class CV:
    """
    Cross validation class

    :Example:
        rf_CV = CV(df_train, df_target, n_splits = 10, model = rf, params = rf_params)
        mae_mean, mae_std = rf_CV.KFold()
        rf_CV.plot_importances('filename')

    """
    def __init__(self, df_train, df_target,  model, n_splits = 10, **params):
        """

        :param params: model parameters
        :param df_train: training set
        :param df_target: vector of targets
        :param n_splits: number of KFolds
        :param model: model to get cross validation for

        """
        self.params = params
        self.df_train = df_train
        self.df_target = df_target
        self.n_splits = n_splits
        self.model = model

    def KFold(self, scoring):
        """
        Kfold cross validation
        """
        kf = KFold(n_splits = self.n_splits, random_state = 0)
        self.model.fit(self.df_train, self.df_target)
        results = cross_val_score(self.model, self.df_train, self.df_target, cv = kf, scoring = scoring)
        return results.mean(), results.std()

    def LOO(self):
        """
        Leave one out CV
        """
        loo = LeaveOneOut()
        loo_splits = [(train_index, test_index) for train_index, test_index in loo.split(self.df_train)]
        return loo_splits

    def plot_importances(self, filename = None):
        """
        Plot feature importances for decision tree based classifiers/regressors
        """
        features = self.model.feature_importances_
        cols = self.df_train.columns
        feature_dataframe = pd.DataFrame({'features': cols,
        'Random Forest feature importances': features})
        # # Feature importances bar plot
        # plt.figure(figsize = (10,10))
        # ax = feature_dataframe.plot(kind = 'barh', legend = None, figsize = (20,50))
        ax = feature_dataframe.plot(kind = 'barh', legend = None)
        ax.set_yticklabels(feature_dataframe['features'].values, fontsize = 10)
        # ax.tick_params(axis='y', which='major', pad=15)
        # ax.yaxis.set_ticks(np.arange(0,len(feature_dataframe),2))
        plt.tight_layout()
        plt.savefig('Plots/{}.png'.format(filename))

    def report(self, model, date, mean, std):
        """
        print report to stdout as well as store in a file
        """
        original = sys.stdout
        sys.stdout = open('Logs/{}.txt'.format(date), 'a')
        print('{} Scores : {} ({})\n{}'.format(model, mean, std, self.params))
        sys.stdout = original
        print('{} Scores : {} ({})'.format(model, mean, std))
