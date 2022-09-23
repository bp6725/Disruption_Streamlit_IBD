import os
import collections
import pickle
import pandas as pd

class GlobalUtils():
    BACKBONE_WIDG_FOR_DISRUPTION_MODULE_NET = ['reference', 'med', 'scd', 'med with nut', 'scd with nut',
                                               'before med with nuts', 'before scd with nuts', 'refrence with nuts']
    DIET_WIDG_FOR_DISRUPTION_MODULE_NET = ['reference', 'med', 'scd', 'med with nut', 'scd with nut',
                                           'before med with nuts', 'before scd with nuts', 'refrence with nuts']

    @staticmethod
    def cache_me(func):
        def wrapper_cache_me(*args, **kwargs):
            if not 'for_cache' in kwargs.keys():
                return func(*args)
            okw = collections.OrderedDict(sorted(kwargs['for_cache'].items()))
            cache_file_name = '../cache/' + func.__name__ + ' ' + ' '.join(
                [str(k) + '(' + str(v) + ')' for k, v in okw.items()])
            if (os.path.isfile(cache_file_name)):
                with open(cache_file_name, 'rb') as f:
                    try :
                        value = pickle.load(f)
                    except :
                        value = pd.read_pickle(cache_file_name)
            else:
                value = func(*args)
                with open(cache_file_name, 'wb') as f:
                    pickle.dump(value, f)
            return value

        return wrapper_cache_me

    @staticmethod
    def copy_df(df):
        if (type(df) is pd.DataFrame):
            return pd.DataFrame(columns=df.columns.copy(deep=True), index=df.index.copy(deep=True),
                                data=df.values.copy())
        else:
            return pd.Series(data=df.values.copy(), index=df.index.copy(deep=True))
