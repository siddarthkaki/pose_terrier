function [data_table] = f_statistics(mobilepose_data, ptl_data)
%F_STATISTICS Summary of this function goes here

data_temp = denan(mobilepose_data.ErrorRotationInterp);
mp_ER_mean = mean(data_temp);
mp_ER_median = median(data_temp);
mp_ER_std = std(data_temp);

data_temp = denan(mobilepose_data.NormErrorTranslationInterp);
mp_ET_mean = mean(data_temp);
mp_ET_median = median(data_temp);
mp_ET_std = std(data_temp);

data_temp = denan(mobilepose_data.NormErrorTranslationNormedInterp);
mp_ETN_mean = mean(data_temp);
mp_ETN_median = median(data_temp);
mp_ETN_std = std(data_temp);

data_temp = denan(mobilepose_data.NormErrorTranslationNormedInterp + deg2rad(mobilepose_data.ErrorRotationInterp));
mp_EC_mean = mean(data_temp);
mp_EC_median = median(data_temp);
mp_EC_std = std(data_temp);


data_temp = denan(ptl_data.ErrorRotation);
fl_ER_mean = mean(data_temp);
fl_ER_median = median(data_temp);
fl_ER_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslation);
fl_ET_mean = mean(data_temp);
fl_ET_median = median(data_temp);
fl_ET_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslationNormed);
fl_ETN_mean = mean(data_temp);
fl_ETN_median = median(data_temp);
fl_ETN_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslationNormed + deg2rad(ptl_data.ErrorRotation));
fl_EC_mean = mean(data_temp);
fl_EC_median = median(data_temp);
fl_EC_std = std(data_temp);


Metric = ["E_R [deg]";"E_R [deg]";"E_R [deg]"; "E_T [m]";"E_T [m]";"E_T [m]"; "E_TN";"E_TN";"E_TN"; "E_C";"E_C";"E_C";];
Statistic = ["Mean";"Median";"Std. Dev."; "Mean";"Median";"Std. Dev."; "Mean";"Median";"Std. Dev."; "Mean";"Median";"Std. Dev."];
Unfiltered = [  mp_ER_mean ;mp_ER_median ;mp_ER_std ;   mp_ET_mean;mp_ET_median;mp_ET_std; ...
                mp_ETN_mean;mp_ETN_median;mp_ETN_std;   mp_EC_mean;mp_EC_median;mp_EC_std];
Filtered = [  fl_ER_mean ;fl_ER_median ;fl_ER_std ;   fl_ET_mean;fl_ET_median;fl_ET_std; ...
                fl_ETN_mean;fl_ETN_median;fl_ETN_std;   fl_EC_mean;fl_EC_median;fl_EC_std];

Unfiltered = round(Unfiltered,3);
Filtered = round(Filtered,3);

data_table = table(Metric, Statistic, Unfiltered, Filtered);
end

