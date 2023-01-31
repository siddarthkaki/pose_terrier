function [data_table] = f_statistics_quatera(ptl_data,iidx)
%F_STATISTICS Summary of this function goes here


data_temp = denan(ptl_data.ErrorRotation(iidx));
fl_ER_mean = mean(data_temp);
fl_ER_median = median(data_temp);
fl_ER_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslation(1:end));
fl_ET_mean = mean(data_temp);
fl_ET_median = median(data_temp);
fl_ET_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslationNormed(1:end));
fl_ETN_mean = mean(data_temp);
fl_ETN_median = median(data_temp);
fl_ETN_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslationNormed(iidx) + deg2rad(ptl_data.ErrorRotation(iidx)));
fl_EC_mean = mean(data_temp);
fl_EC_median = median(data_temp);
fl_EC_std = std(data_temp);



data_temp = denan(ptl_data.ErrorRotationNLS(iidx));
un_ER_mean = mean(data_temp);
un_ER_median = median(data_temp);
un_ER_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslationNLS(1:end));
un_ET_mean = mean(data_temp);
un_ET_median = median(data_temp);
un_ET_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslationNormedNLS(1:end));
un_ETN_mean = mean(data_temp);
un_ETN_median = median(data_temp);
un_ETN_std = std(data_temp);

data_temp = denan(ptl_data.NormErrorTranslationNormedNLS(iidx) + deg2rad(ptl_data.ErrorRotationNLS(iidx)));
un_EC_mean = mean(data_temp);
un_EC_median = median(data_temp);
un_EC_std = std(data_temp);




Metric = ["E_R [deg]";"E_R [deg]";"E_R [deg]"; "E_T [m]";"E_T [m]";"E_T [m]"; "E_TN";"E_TN";"E_TN"; "E_C";"E_C";"E_C";];
Statistic = ["Mean";"Median";"Std. Dev."; "Mean";"Median";"Std. Dev."; "Mean";"Median";"Std. Dev."; "Mean";"Median";"Std. Dev."];

Unfiltered = [  un_ER_mean ;un_ER_median ;un_ER_std ;   un_ET_mean;un_ET_median;un_ET_std; ...
                un_ETN_mean;un_ETN_median;un_ETN_std;   un_EC_mean;un_EC_median;un_EC_std];

Filtered = [  fl_ER_mean ;fl_ER_median ;fl_ER_std ;   fl_ET_mean;fl_ET_median;fl_ET_std; ...
                fl_ETN_mean;fl_ETN_median;fl_ETN_std;   fl_EC_mean;fl_EC_median;fl_EC_std];

Unfiltered = round(Unfiltered,3);
Filtered = round(Filtered,3);

data_table = table(Metric, Statistic, Unfiltered, Filtered);
end

