
# written by cagri for portfolio analysis
# this is an early draft and in the process of development

import pandas as pd

# 0 defining the class
class portfolio_analysis:

    #    labels_list = []

    def __init__(self, df):
        self.data = df

    # 0.0 exchange code-share code framer

    def universe_framer(self, exchange_frame, share_frame, valid_exchanges=[1, 2, 3], valid_shares=[10, 11]):

        mask_exchange = exchange_frame.isin(valid_exchanges)
        mask_share = share_frame.isin(valid_shares)

        self.framed = self.data[mask_exchange & mask_share].copy()

    # 0.1 winsorizer and truncater
    def outlier_handling(self, winsorize=True, win_limits=[.0025, .0025], truncate=True, trunc_limits=[.0025, .0025]):
        '''

        :param win_limits:
        :param winsorize:
        :param trunc_limits:
        :param truncate:
        :return:
        '''
        self.cleaned = self.framed.copy()

        # truncating
        if truncate:

            trunc_values = [self.framed.quantile(trunc_limits[0], axis=1),
                            self.framed.quantile(1 - trunc_limits[1], axis=1)]

            mask_trunc = self.framed.ge(trunc_values[0], axis=0) & self.framed.le(trunc_values[1], axis=0)
            self.cleaned = self.cleaned[mask_trunc]

        else:

            trunc_limits = [0, 0]

        if winsorize:
            win_limits = [win_limits[0] + trunc_limits[0], win_limits[1] + trunc_limits[1]]
            win_values = [self.framed.quantile(win_limits[0], axis=1), self.framed.quantile(1 - win_limits[1], axis=1)]

            self.cleaned = self.cleaned.clip(lower=win_values[0], upper=win_values[1], axis=0)

    # 0.2 zscore calculator
    def zscore(self):

        return self.cleaned.sub(self.cleaned.mean(axis=1), axis=0).divide(self.cleaned.std(axis=1), axis=0)

    # 0.4 portfolio calculator
    def portfolio_cut(df, numof_port, labels_list=None):

        if labels_list is None:
            labels_list = [i for i in range(1, numof_port + 1)]

        labels_bins = [1 / numof_port * i for i in range(numof_port + 1)]

        return df.apply(lambda x: pd.qcut(x, labels_bins, labels=labels_list), axis=1), labels_list

    def portfolio_generator(self, sorting_var=None, numof_port=10, double_sort=False):
        """this function assigns each stock to a portfolio for each date point in the time series
        """
        #        nonlocal labels_list
        self.numof_port = numof_port
        self.double_sort = double_sort

        if not self.double_sort:

            self.portfolios, self.labels_list = portfolio_analysis.portfolio_cut(self.cleaned, numof_port)

        else:

            #            self.labels_list = [10*i+j for i in range(1,numof_port[0]+1) for j in range(1,numof_port[1]+1)] # here labels are like this 11-12-13-21-22-23... etc
            self.labels_list = []

            sorting_ports, labels_sorting = portfolio_analysis.portfolio_cut(sorting_var[self.cleaned.notnull()],
                                                                             numof_port[0])
            self.portfolios = pd.DataFrame(index=self.cleaned.index, columns=self.cleaned.columns)

            for label in labels_sorting:
                mask = sorting_ports == label
                labels_list = [10 * label + i for i in
                               range(1, numof_port[1] + 1)]  # here labels are like this 11-12-13-21-22-23... etc
                self.labels_list += labels_list

                sorting_second, _ = portfolio_analysis.portfolio_cut(self.cleaned[mask], numof_port[1],
                                                                     labels_list=labels_list)

                self.portfolios.update(sorting_second)

    # 0.5 portfolio return calculator
    def port_retCalc(self, df_ret, df_mcap=None, VW=False):
        '''

        :param df_ret:
        :param df_mcap:
        :param VW:
        :return:
        '''

        self.mcaps = df_mcap

        self.port_rets_VW = pd.DataFrame(index=self.portfolios.index, columns=self.labels_list)
        self.port_rets_EW = self.port_rets_VW.copy()

        for p in self.labels_list:

            mask = self.portfolios == p

            if VW:

                weights = df_mcap[mask]
                weights = weights.divide(weights.sum(axis=1), axis=0)

                self.port_rets_VW[p] = df_ret.multiply(weights).sum(axis=1)
                self.port_rets_EW[p] = df_ret[mask].mean(axis=1)

            else:

                self.port_rets_EW[p] = df_ret[mask].mean(axis=1)

    # 0.6 factor calculator
    def factor_calc(self, ascending=True):
        '''

        :param ascending:
        :return:
        '''

        if self.double_sort:

            self.factor_EW = 0
            self.factor_VW = 0

            if ascending:

                for p in range(1, self.numof_port[0] + 1):
                    self.factor_EW += (self.port_rets_EW[10 * p + self.numof_port[1]] - self.port_rets_EW[10 * p + 1]) / \
                                      self.numof_port[0]
                    self.factor_VW += (self.port_rets_VW[10 * p + self.numof_port[1]] - self.port_rets_VW[10 * p + 1]) / \
                                      self.numof_port[0]

            else:

                for p in range(1, self.numof_port[0] + 1):
                    self.factor_EW += (self.port_rets_EW[10 * p + 1] - self.port_rets_EW[10 * p + self.numof_port[1]]) / \
                                      self.numof_port[0]
                    self.factor_VW += (self.port_rets_VW[10 * p + 1] - self.port_rets_VW[10 * p + self.numof_port[1]]) / \
                                      self.numof_port[0]
        else:

            if ascending:

                self.factor_EW = self.port_rets_EW[max(self.labels_list)] - self.port_rets_EW[min(self.labels_list)]
                self.factor_VW = self.port_rets_VW[max(self.labels_list)] - self.port_rets_VW[min(self.labels_list)]
            else:

                self.factor_EW = self.port_rets_EW[min(self.labels_list)] - self.port_rets_EW[max(self.labels_list)]
                self.factor_VW = self.port_rets_VW[min(self.labels_list)] - self.port_rets_VW[max(self.labels_list)]

