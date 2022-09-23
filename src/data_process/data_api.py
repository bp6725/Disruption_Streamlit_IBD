import pickle
import pandas as pd
from sklearn.neighbors import KNeighborsRegressor


class DataAPI() :
    def __init__(self):
        self.union_ref_df = self.cache_store("union_ref_df", False)
        self.union_scd_df = self.cache_store("union_scd_df", False)
        self.union_med_df = self.cache_store("union_med_df", False)
        self.rejected_med_df = self.cache_store("rejected_med_df", False)
        self.rejected_scd_df = self.cache_store("rejected_scd_df", False)
        self.rejected_before_med_df = self.cache_store("rejected_before_med_df", False)
        self.rejected_before_scd_df = self.cache_store("rejected_before_scd_df", False)
        self.rejected_refrence_med_df = self.cache_store("rejected_refrence_med_df", False)
        self.rejected_refrence_scd_df = self.cache_store("rejected_refrence_scd_df", False)
        self.imp_id_to_path = self.cache_store("imp_id_to_path", False)
        self.unified_features_df = self.cache_store("unified_features_df", False)
        self.all_features_df = self.cache_store("all_features_df", False)

        self.global_diet_to_idx = None

        # remove aggregated bioms
        with open(r"../data/dfs/biom_feautres_df.pkl", 'rb') as f:
            biom_feautres_df = pickle.load(f)
        self.unagg_bioms = biom_feautres_df.drop(columns=['PN', 'visit', 'PN_visit', 'Diet', ]).columns

        _all_features_df = self.all_features_df.reset_index().set_index('PN')
        self.med_sig_df = self.union_med_df
        self.scd_sig_df = self.union_scd_df
        self.med_nutrients_df = self.rejected_med_df
        self.scd_nutrients_df = self.rejected_scd_df

        self.before_med_nutrients_df = self.rejected_before_med_df
        self.before_scd_nutrients_df = self.rejected_before_scd_df

        med_with_nutriens = pd.concat([self.med_sig_df, self.med_nutrients_df], axis=1,sort=True)
        self.med_with_nutriens = med_with_nutriens.fillna(med_with_nutriens.mean())

        scd_with_nutriens = pd.concat([self.scd_sig_df, self.scd_nutrients_df], axis=1,sort=True)
        self.scd_with_nutriens = scd_with_nutriens.fillna(scd_with_nutriens.mean())

        _before_med_sig_df = _all_features_df[_all_features_df['Diet'] == 'Before_MED'][self.med_sig_df.columns]
        _before_scd_sig_df = _all_features_df[_all_features_df['Diet'] == 'Before_SCD'][self.scd_sig_df.columns]

        before_med_sig_df_with_nutriens = pd.concat([_before_med_sig_df, self.before_med_nutrients_df], axis=1,sort=True)
        self.before_med_sig_df_with_nutriens = before_med_sig_df_with_nutriens.fillna(
            before_med_sig_df_with_nutriens.mean())

        before_scd_sig_df_with_nutriens = pd.concat([_before_scd_sig_df, self.before_scd_nutrients_df], axis=1,sort=True)
        self.before_scd_sig_df_with_nutriens = before_scd_sig_df_with_nutriens.fillna(before_scd_sig_df_with_nutriens.mean())

        self.update_all_with_nut()

        self.feature_to_type_dict = self.build_feature_to_type_dict()
        self.feature_name_dict = self.build_feature_name_dict()

    #region INIT

    def update_all_with_nut(self):
        all_before_med_nutrients_df = self.before_med_nutrients_df.merge(self.all_features_df[self.all_features_df['Diet'] == 'Before_MED'], left_index=True, right_on=['PN'])
        all_before_scd_nutrients_df = self.before_scd_nutrients_df.merge(self.all_features_df[self.all_features_df['Diet'] == 'Before_SCD'], left_index=True, right_on=['PN'])
        rejected_reference_with_nutriens = all_before_scd_nutrients_df.columns.unique().union(all_before_med_nutrients_df.columns.unique())

        all_before_scd_with_nut = self.all_features_df[self.all_features_df['Diet'] == 'Before_SCD'].merge(
            self.rejected_refrence_scd_df, right_index=True, left_on=['PN'])
        all_before_med_with_nut = self.all_features_df[self.all_features_df['Diet'] == 'Before_MED'].merge(
            self.rejected_refrence_med_df, right_index=True, left_on=['PN'])

        all_refrence_with_nuts = pd.concat([all_before_scd_with_nut, all_before_med_with_nut])
        all_refrence_with_nuts = all_refrence_with_nuts[
            rejected_reference_with_nutriens.union(self.union_ref_df.columns)]

        cols_to_drop = all_refrence_with_nuts.columns[all_refrence_with_nuts.isna().sum() > 10].tolist()
        all_refrence_with_nuts = all_refrence_with_nuts.drop(columns=cols_to_drop)

        cols_to_impute = all_refrence_with_nuts.columns[all_refrence_with_nuts.isna().any()].tolist()
        cols_non_num = [col for col in cols_to_impute if all_refrence_with_nuts[col].dtypes != float]
        cols_to_impute = [col for col in cols_to_impute if col not in cols_non_num]
        all_refrence_with_nuts = all_refrence_with_nuts.drop(columns=cols_non_num)

        # imputer = IterativeImputer(estimator=KNeighborsRegressor())
        # all_refrence_with_nuts[cols_to_impute] = imputer.fit_transform(all_refrence_with_nuts[cols_to_impute])

        self.union_ref_with_nutrients_df = \
            pd.concat([all_before_med_nutrients_df, all_before_scd_nutrients_df], axis=1,sort=True)[
                rejected_reference_with_nutriens]
        self.all_not_sig_refrence_with_nuts = all_refrence_with_nuts

    def build_feature_name_dict(self):
        feature_name_dict_str = {str(feature): self.extract_name(feature) for feature in
                                 self.union_ref_df.columns.union(self.med_nutrients_df.columns).union(
                                     self.scd_nutrients_df.columns).tolist()}
        feature_name_dict_int = {feature: self.extract_name(feature) for feature in
                                 self.union_ref_df.columns.union(self.med_nutrients_df.columns).union(
                                     self.scd_nutrients_df.columns).tolist()}
        feature_name_dict = {**feature_name_dict_str, **feature_name_dict_int}

        for k, v in self.imp_id_to_path.items():
            if (not k in feature_name_dict.keys()):
                feature_name_dict[k] = v
                feature_name_dict[str(k)] = v

        for fea in self.unified_features_df['biom'].columns:
            if (not fea in feature_name_dict.keys()):
                feature_name_dict[fea] = self.extract_name(fea)

        return feature_name_dict

    def build_feature_to_type_dict(self) :
        feature_to_type_dict = {}
        for col in self.unified_features_df.columns:
            feature_to_type_dict[col[1]] = col[0]

        for nut in self.med_nutrients_df.columns.union(self.scd_nutrients_df.columns):
            if nut not in feature_to_type_dict.keys():
                feature_to_type_dict[nut] = 'nutrients'

        for nut in self.all_not_sig_refrence_with_nuts.columns:
            if nut not in feature_to_type_dict.keys():
                feature_to_type_dict[nut] = 'nutrients'

        return feature_to_type_dict

    def extract_name(self,feature):
        if feature in self.imp_id_to_path.keys():
            return self.imp_id_to_path[feature]

        if (type(feature) is not str):
            return feature

        if '(' in feature:
            return '(' + ','.join(feature.split(')')[0].split(',')[-2:]) + ')'

        return feature

    def feature_to_name(self,feature, feature_name_dict):
        def extract_name(feature):
            if feature in self.imp_id_to_path.keys():
                return self.imp_id_to_path[feature]

            if (type(feature) is not str):
                return feature

            if '(' in feature:
                return '(' + ','.join(feature.split(')')[0].split(',')[-2:]) + ')'

            return feature

        if (feature in feature_name_dict.keys()):
            return feature_name_dict[feature]
        if (type(feature) is str):
            return extract_name(feature)

    def return_unagg_biom_df(self,df, feature_to_type_dict):
        agg_cols = df[[col for col in df.head().columns if feature_to_type_dict[col] == 'biom']].columns.difference(
            self.unagg_bioms)
        return df.drop(columns=agg_cols)

    def  cache_store(self,file_name,save,file = None) :
        path = f"../cache/store/{file_name}.pkl"
        mode = "wb" if save else "rb"
        with open(path,mode) as f :
            if save :
                pickle.dump(file,f)
            else :
                try :
                    return pickle.load(f)
                except :
                    return pd.read_pickle(path)

    #endregion

    def return_unagg_biom_data(self):
        med_sig_df = self.return_unagg_biom_df(self.med_sig_df, self.feature_to_type_dict)
        scd_sig_df = self.return_unagg_biom_df(self.scd_sig_df, self.feature_to_type_dict)
        union_ref_df = self.return_unagg_biom_df(self.union_ref_df, self.feature_to_type_dict)
        med_with_nutriens = self.return_unagg_biom_df(self.med_with_nutriens, self.feature_to_type_dict)
        scd_with_nutriens = self.return_unagg_biom_df(self.scd_with_nutriens, self.feature_to_type_dict)
        before_med_sig_df_with_nutriens = self.return_unagg_biom_df(self.before_med_sig_df_with_nutriens,
                                                                    self.feature_to_type_dict)
        before_scd_sig_df_with_nutriens = self.return_unagg_biom_df(self.before_scd_sig_df_with_nutriens,
                                                                    self.feature_to_type_dict)
        all_refrence_with_nuts = pd.concat([before_med_sig_df_with_nutriens,
                                            before_scd_sig_df_with_nutriens],sort=True)

        return med_sig_df,scd_sig_df,union_ref_df,med_with_nutriens,scd_with_nutriens,\
               before_med_sig_df_with_nutriens,before_scd_sig_df_with_nutriens,all_refrence_with_nuts

    def return_unagg_global_diet_to_idx(self):
        if self.global_diet_to_idx is not None :
            return self.global_diet_to_idx

        med_sig_df, scd_sig_df, union_ref_df,\
        med_with_nutriens, scd_with_nutriens,\
        before_med_sig_df_with_nutriens,\
        before_scd_sig_df_with_nutriens, all_refrence_with_nuts = self.return_unagg_biom_data()

        return {0: union_ref_df, 1: self.union_med_df, 2: self.union_scd_df
            , 3: med_with_nutriens, 4: scd_with_nutriens
            , 5: before_med_sig_df_with_nutriens, 6: before_scd_sig_df_with_nutriens
            , 7: all_refrence_with_nuts}

    def return_unagg_by_index(self,idx):
        return self.return_unagg_global_diet_to_idx()[idx]

    def extract_deltas_measures(self,explanatory_features,stage_name):
        '''
        return the fold change per a diet for explanatory_features
        :param explanatory_features: the features we calculate the fold change for
        :param stage_name: the stage we calculate the feature for - always the *real* diet stage. so for dilation - diet_idx, and for creation - ref_idx
        :return:foldchange per feature per person
        '''

        full_df = self.all_features_df.copy(deep=True)
        full_df = full_df.set_index('PN', drop=False)

        after_df = full_df[full_df['Diet'] == 'After_' + stage_name]
        before_df = full_df[full_df['Diet'] == 'Before_' + stage_name]

        res = (after_df[explanatory_features] - before_df[explanatory_features]) / before_df[explanatory_features]
        res.columns = [col + '_delta' for col in res.columns]
        return res

    def get_diet_df_with_delta(self,diet_df,explanatory_features,stage_name):
        _delta_def = self.extract_deltas_measures(explanatory_features,stage_name)
        return diet_df.merge(_delta_def,left_index=True,right_index=True)


    def return_relvent_sc_imp_data(self,stage):
        unified_features_df = self.unified_features_df.copy(deep=True)
        relvent_24 = unified_features_df.loc[~unified_features_df["PN"].isin(["PN10", "PN27", "PN28", "PN29"])]
        relvent_04 = unified_features_df.loc[unified_features_df["PN"].isin(["PN10", "PN27", "PN28", "PN29"])]

        if stage == "Before_MED":
            _24 = relvent_24[relvent_24["Diet"] == "Before_MED"].set_index("PN")["ScaledImpdata"]
            _04 = relvent_04[relvent_04["Diet"] == "Before_SCD"].set_index("PN")["ScaledImpdata"]
            return _24.append(_04)

        if stage == "After_MED":
            _24 = relvent_24[relvent_24["Diet"] == "After_MED"].set_index("PN")["ScaledImpdata"]
            _04 = relvent_04[relvent_04["Diet"] == "After_MED"].set_index("PN")["ScaledImpdata"]
            return _24.append(_04)

        if stage == "After_SCD":
            _24 = relvent_24[relvent_24["Diet"] == "After_SCD"].set_index("PN")["ScaledImpdata"]
            _04 = relvent_04[relvent_04["Diet"] == "After_SCD"].set_index("PN")["ScaledImpdata"]
            return _24.append(_04)

    def rename_index(self,ind):
        if type(ind) is int:
            return f"{ind} {self.imp_id_to_path[ind]}"

        if type(ind) is str:
            if "-" in ind:
                _f_dis_ind, _s_dis_ind = ind.split("-")
                if _f_dis_ind.isdigit():
                    f_dis_ind = f"{_f_dis_ind} {self.imp_id_to_path[int(_f_dis_ind)]}"
                else:
                    f_dis_ind = _f_dis_ind

                if _s_dis_ind.isdigit():
                    s_dis_ind = f"{_s_dis_ind} {self.imp_id_to_path[int(_s_dis_ind)]}"
                else:
                    s_dis_ind = _s_dis_ind

                return f"({f_dis_ind} - {s_dis_ind})"

        dis_ind = ind[0]
        ser_ind = ind[1]

        _f_dis_ind, _s_dis_ind = dis_ind.split("-")

        if _f_dis_ind.isdigit():
            f_dis_ind = f"{_f_dis_ind} {self.imp_id_to_path[int(_f_dis_ind)]}"
        else:
            f_dis_ind = _f_dis_ind

        if _s_dis_ind.isdigit():
            s_dis_ind = f"{_s_dis_ind} {self.imp_id_to_path[int(_s_dis_ind)]}"
        else:
            s_dis_ind = _s_dis_ind

        new_dis_ind = f"({f_dis_ind} - {s_dis_ind})"
        new_ser_ind = f"{ser_ind} {self.imp_id_to_path[ser_ind]}"

        new_ind = (new_dis_ind, new_ser_ind)

        return new_ind
