python nist_wct_log_spectral_scan_FIXED.py --csv Fe_lines.csv --null-n 500 --min-lines 100
C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST\nist_wct_log_spectral_scan_FIXED.py:90: FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.
  return df.applymap(clean_cell)
[raw] rows=20615 cols=['element', 'sp_num', 'obs_wl_air(nm)', 'unc_obs_wl', 'ritz_wl_air(nm)', 'unc_ritz_wl', 'wn(cm-1)', 'intens', 'Aki(s^-1)', 'Acc', 'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k', 'term_k', 'J_k', 'Type', 'tp_ref', 'line_ref', 'Unnamed: 21']
[clean] usable unique lines=19577
[scan] bins=160 k=0.5..80.0 n_k=2500
[null] 100/500
[null] 200/500
[null] 300/500
[null] 400/500
[null] 500/500

[done]
{
  "k_best": 31.358343337334933,
  "deltaD": 107.12175266766758,
  "scan_null_p": 0.001996007984031936,
  "tail_count_ge": 0,
  "null_n": 500,
  "n_obs": 10.752733902312999,
  "n_unique_lines": 19577
}

[nearest branch]
         branch         n  k_target   k_error  abs_k_error   n_error  abs_n_error
integer_1_to_40 11.000000 32.079449 -0.721105     0.721105 -0.247266     0.247266
 koide_10_15_20 10.000000 29.163135  2.195208     2.195208  0.752734     0.752734
integer_1_to_40 10.000000 29.163135  2.195208     2.195208  0.752734     0.752734
integer_1_to_40 12.000000 34.995762 -3.637419     3.637419 -1.247266     1.247266
integer_1_to_40  9.000000 26.246822  5.111522     5.111522  1.752734     1.752734
integer_1_to_40 13.000000 37.912076 -6.553733     6.553733 -2.247266     2.247266
  folded_4over9 13.333333 38.884180 -7.525837     7.525837 -2.580599     2.580599
integer_1_to_40  8.000000 23.330508  8.027835     8.027835  2.752734     2.752734
integer_1_to_40 14.000000 40.828389 -9.470046     9.470046 -3.247266     3.247266
integer_1_to_40  7.000000 20.414195 10.944149    10.944149  3.752734     3.752734
(cupy-env) PS C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST> python nist_wct_log_spectral_scan_FIXED.py --csv Fe_lines.csv --ion 2 --null-n 500 --min-lines 100

C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST\nist_wct_log_spectral_scan_FIXED.py:90: FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.
  return df.applymap(clean_cell)
[raw] rows=20615 cols=['element', 'sp_num', 'obs_wl_air(nm)', 'unc_obs_wl', 'ritz_wl_air(nm)', 'unc_ritz_wl', 'wn(cm-1)', 'intens', 'Aki(s^-1)', 'Acc', 'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k', 'term_k', 'J_k', 'Type', 'tp_ref', 'line_ref', 'Unnamed: 21']
[clean] usable unique lines=9447
[scan] bins=160 k=0.5..80.0 n_k=2500
[null] 100/500
[null] 200/500
[null] 300/500
[null] 400/500
[null] 500/500

[done]
{
  "k_best": 31.3265306122449,
  "deltaD": 315.7276606037286,
  "scan_null_p": 0.001996007984031936,
  "tail_count_ge": 0,
  "null_n": 500,
  "n_obs": 10.739649568878214,
  "n_unique_lines": 9447
}

[nearest branch]
         branch         n  k_target   k_error  abs_k_error   n_error  abs_n_error
integer_1_to_40 11.000000 32.085948 -0.759417     0.759417 -0.260350     0.260350
 koide_10_15_20 10.000000 29.169044  2.157487     2.157487  0.739650     0.739650
integer_1_to_40 10.000000 29.169044  2.157487     2.157487  0.739650     0.739650
integer_1_to_40 12.000000 35.002852 -3.676322     3.676322 -1.260350     1.260350
integer_1_to_40  9.000000 26.252139  5.074391     5.074391  1.739650     1.739650
integer_1_to_40 13.000000 37.919757 -6.593226     6.593226 -2.260350     2.260350
  folded_4over9 13.333333 38.892058 -7.565527     7.565527 -2.593684     2.593684
integer_1_to_40  8.000000 23.335235  7.991296     7.991296  2.739650     2.739650
integer_1_to_40 14.000000 40.836661 -9.510130     9.510130 -3.260350     3.260350
integer_1_to_40  7.000000 20.418330 10.908200    10.908200  3.739650     3.739650
(cupy-env) PS C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST> python nist_wct_log_spectral_scan_FIXED.py --csv Fe_lines.csv --ion 3 --null-n 500 --min-lines 100

C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST\nist_wct_log_spectral_scan_FIXED.py:90: FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.
  return df.applymap(clean_cell)
[raw] rows=20615 cols=['element', 'sp_num', 'obs_wl_air(nm)', 'unc_obs_wl', 'ritz_wl_air(nm)', 'unc_ritz_wl', 'wn(cm-1)', 'intens', 'Aki(s^-1)', 'Acc', 'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k', 'term_k', 'J_k', 'Type', 'tp_ref', 'line_ref', 'Unnamed: 21']
[clean] usable unique lines=1524
[scan] bins=160 k=0.5..80.0 n_k=2500
[null] 100/500
[null] 200/500
[null] 300/500
[null] 400/500
[null] 500/500

[done]
{
  "k_best": 35.27130852340936,
  "deltaD": 54.688241372315474,
  "scan_null_p": 0.001996007984031936,
  "tail_count_ge": 0,
  "null_n": 500,
  "n_obs": 9.451543028800112,
  "n_unique_lines": 1524
}

[nearest branch]
         branch         n  k_target    k_error  abs_k_error   n_error  abs_n_error
integer_1_to_40  9.000000 33.586238   1.685070     1.685070  0.451543     0.451543
 koide_10_15_20 10.000000 37.318043  -2.046734     2.046734 -0.548457     0.548457
integer_1_to_40 10.000000 37.318043  -2.046734     2.046734 -0.548457     0.548457
integer_1_to_40  8.000000 29.854434   5.416874     5.416874  1.451543     1.451543
integer_1_to_40 11.000000 41.049847  -5.778538     5.778538 -1.548457     1.548457
integer_1_to_40  7.000000 26.122630   9.148679     9.148679  2.451543     2.451543
integer_1_to_40 12.000000 44.781651  -9.510343     9.510343 -2.548457     2.548457
  folded_4over9  6.666667 24.878695  10.392613    10.392613  2.784876     2.784876
integer_1_to_40  6.000000 22.390826  12.880483    12.880483  3.451543     3.451543
integer_1_to_40 13.000000 48.513455 -13.242147    13.242147 -3.548457     3.548457

 python nist_wct_log_spectral_scan_FIXED.py --csv Fe_lines.csv --ion 2 --null-n 5000 --min-lines 100
C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST\nist_wct_log_spectral_scan_FIXED.py:90: FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.
  return df.applymap(clean_cell)
[raw] rows=20615 cols=['element', 'sp_num', 'obs_wl_air(nm)', 'unc_obs_wl', 'ritz_wl_air(nm)', 'unc_ritz_wl', 'wn(cm-1)', 'intens', 'Aki(s^-1)', 'Acc', 'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k', 'term_k', 'J_k', 'Type', 'tp_ref', 'line_ref', 'Unnamed: 21']
[clean] usable unique lines=9447
[scan] bins=160 k=0.5..80.0 n_k=2500
[null] 100/5000
[null] 200/5000
[null] 300/5000
[null] 400/5000
[null] 500/5000
[null] 600/5000
[null] 700/5000
[null] 800/5000
[null] 900/5000
[null] 1000/5000
[null] 1100/5000
[null] 1200/5000
[null] 1300/5000
[null] 1400/5000
[null] 1500/5000
[null] 1600/5000
[null] 1700/5000
[null] 1800/5000
[null] 1900/5000
[null] 2000/5000
[null] 2100/5000
[null] 2200/5000
[null] 2300/5000
[null] 2400/5000
[null] 2500/5000
[null] 2600/5000
[null] 2700/5000
[null] 2800/5000
[null] 2900/5000
[null] 3000/5000
[null] 3100/5000
[null] 3200/5000
[null] 3300/5000
[null] 3400/5000
[null] 3500/5000
[null] 3600/5000
[null] 3700/5000
[null] 3800/5000
[null] 3900/5000
[null] 4000/5000
[null] 4100/5000
[null] 4200/5000
[null] 4300/5000
[null] 4400/5000
[null] 4500/5000
[null] 4600/5000
[null] 4700/5000
[null] 4800/5000
[null] 4900/5000
[null] 5000/5000

[done]
{
  "k_best": 31.3265306122449,
  "deltaD": 315.7276606037286,
  "scan_null_p": 0.0001999600079984003,
  "tail_count_ge": 0,
  "null_n": 5000,
  "n_obs": 10.739649568878214,
  "n_unique_lines": 9447
}

[nearest branch]
         branch         n  k_target   k_error  abs_k_error   n_error  abs_n_error
integer_1_to_40 11.000000 32.085948 -0.759417     0.759417 -0.260350     0.260350
 koide_10_15_20 10.000000 29.169044  2.157487     2.157487  0.739650     0.739650
integer_1_to_40 10.000000 29.169044  2.157487     2.157487  0.739650     0.739650
integer_1_to_40 12.000000 35.002852 -3.676322     3.676322 -1.260350     1.260350
integer_1_to_40  9.000000 26.252139  5.074391     5.074391  1.739650     1.739650
integer_1_to_40 13.000000 37.919757 -6.593226     6.593226 -2.260350     2.260350
  folded_4over9 13.333333 38.892058 -7.565527     7.565527 -2.593684     2.593684
integer_1_to_40  8.000000 23.335235  7.991296     7.991296  2.739650     2.739650
integer_1_to_40 14.000000 40.836661 -9.510130     9.510130 -3.260350     3.260350
integer_1_to_40  7.000000 20.418330 10.908200    10.908200  3.739650     3.739650
(cupy-env) PS C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST> python nist_wct_log_spectral_scan_FIXED.py --csv Fe_lines.csv --ion 2 --bins 120 --null-n 5000 --min-lines 100
C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST\nist_wct_log_spectral_scan_FIXED.py:90: FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.
  return df.applymap(clean_cell)
[raw] rows=20615 cols=['element', 'sp_num', 'obs_wl_air(nm)', 'unc_obs_wl', 'ritz_wl_air(nm)', 'unc_ritz_wl', 'wn(cm-1)', 'intens', 'Aki(s^-1)', 'Acc', 'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k', 'term_k', 'J_k', 'Type', 'tp_ref', 'line_ref', 'Unnamed: 21']
[clean] usable unique lines=9447
[scan] bins=120 k=0.5..80.0 n_k=2500
[null] 100/5000
[null] 200/5000
[null] 300/5000
[null] 400/5000
[null] 500/5000
[null] 600/5000
[null] 700/5000
[null] 800/5000
[null] 900/5000
[null] 1000/5000
[null] 1100/5000
[null] 1200/5000
[null] 1300/5000
[null] 1400/5000
[null] 1500/5000
[null] 1600/5000
[null] 1700/5000
[null] 1800/5000
[null] 1900/5000
[null] 2000/5000
[null] 2100/5000
[null] 2200/5000
[null] 2300/5000
[null] 2400/5000
[null] 2500/5000
[null] 2600/5000
[null] 2700/5000
[null] 2800/5000
[null] 2900/5000
[null] 3000/5000
[null] 3100/5000
[null] 3200/5000
[null] 3300/5000
[null] 3400/5000
[null] 3500/5000
[null] 3600/5000
[null] 3700/5000
[null] 3800/5000
[null] 3900/5000
[null] 4000/5000
[null] 4100/5000
[null] 4200/5000
[null] 4300/5000
[null] 4400/5000
[null] 4500/5000
[null] 4600/5000
[null] 4700/5000
[null] 4800/5000
[null] 4900/5000
[null] 5000/5000

[done]
{
  "k_best": 31.3265306122449,
  "deltaD": 355.7443010501305,
  "scan_null_p": 0.0001999600079984003,
  "tail_count_ge": 0,
  "null_n": 5000,
  "n_obs": 10.717134580264222,
  "n_unique_lines": 9447
}

[nearest branch]
         branch         n  k_target   k_error  abs_k_error   n_error  abs_n_error
integer_1_to_40 11.000000 32.153355 -0.826825     0.826825 -0.282865     0.282865
 koide_10_15_20 10.000000 29.230323  2.096208     2.096208  0.717135     0.717135
integer_1_to_40 10.000000 29.230323  2.096208     2.096208  0.717135     0.717135
integer_1_to_40 12.000000 35.076388 -3.749857     3.749857 -1.282865     1.282865
integer_1_to_40  9.000000 26.307291  5.019240     5.019240  1.717135     1.717135
integer_1_to_40 13.000000 37.999420 -6.672889     6.672889 -2.282865     2.282865
  folded_4over9 13.333333 38.973764 -7.647233     7.647233 -2.616199     2.616199
integer_1_to_40  8.000000 23.384258  7.942272     7.942272  2.717135     2.717135
integer_1_to_40 14.000000 40.922452 -9.595922     9.595922 -3.282865     3.282865
integer_1_to_40  7.000000 20.461226 10.865304    10.865304  3.717135     3.717135
(cupy-env) PS C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST> python nist_wct_log_spectral_scan_FIXED.py --csv Fe_lines.csv --ion 2 --bins 200 --null-n 5000 --min-lines 100
C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST\nist_wct_log_spectral_scan_FIXED.py:90: FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.
  return df.applymap(clean_cell)
[raw] rows=20615 cols=['element', 'sp_num', 'obs_wl_air(nm)', 'unc_obs_wl', 'ritz_wl_air(nm)', 'unc_ritz_wl', 'wn(cm-1)', 'intens', 'Aki(s^-1)', 'Acc', 'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k', 'term_k', 'J_k', 'Type', 'tp_ref', 'line_ref', 'Unnamed: 21']
[clean] usable unique lines=9447
[scan] bins=200 k=0.5..80.0 n_k=2500
[null] 100/5000
[null] 200/5000
[null] 300/5000
[null] 400/5000
[null] 500/5000
[null] 600/5000
[null] 700/5000
[null] 800/5000
[null] 900/5000
[null] 1000/5000
[null] 1100/5000
[null] 1200/5000
[null] 1300/5000
[null] 1400/5000
[null] 1500/5000
[null] 1600/5000
[null] 1700/5000
[null] 1800/5000
[null] 1900/5000
[null] 2000/5000
[null] 2100/5000
[null] 2200/5000
[null] 2300/5000
[null] 2400/5000
[null] 2500/5000
[null] 2600/5000
[null] 2700/5000
[null] 2800/5000
[null] 2900/5000
[null] 3000/5000
[null] 3100/5000
[null] 3200/5000
[null] 3300/5000
[null] 3400/5000
[null] 3500/5000
[null] 3600/5000
[null] 3700/5000
[null] 3800/5000
[null] 3900/5000
[null] 4000/5000
[null] 4100/5000
[null] 4200/5000
[null] 4300/5000
[null] 4400/5000
[null] 4500/5000
[null] 4600/5000
[null] 4700/5000
[null] 4800/5000
[null] 4900/5000
[null] 5000/5000

[done]
{
  "k_best": 31.3265306122449,
  "deltaD": 259.17875676052677,
  "scan_null_p": 0.0001999600079984003,
  "tail_count_ge": 0,
  "null_n": 5000,
  "n_obs": 10.753158562046606,
  "n_unique_lines": 9447
}

[nearest branch]
         branch         n  k_target   k_error  abs_k_error   n_error  abs_n_error
integer_1_to_40 11.000000 32.045639 -0.719108     0.719108 -0.246841     0.246841
 koide_10_15_20 10.000000 29.132399  2.194132     2.194132  0.753159     0.753159
integer_1_to_40 10.000000 29.132399  2.194132     2.194132  0.753159     0.753159
integer_1_to_40 12.000000 34.958879 -3.632348     3.632348 -1.246841     1.246841
integer_1_to_40  9.000000 26.219159  5.107371     5.107371  1.753159     1.753159
integer_1_to_40 13.000000 37.872119 -6.545588     6.545588 -2.246841     2.246841
  folded_4over9 13.333333 38.843199 -7.516668     7.516668 -2.580175     2.580175
integer_1_to_40  8.000000 23.305919  8.020611     8.020611  2.753159     2.753159
integer_1_to_40 14.000000 40.785359 -9.458828     9.458828 -3.246841     3.246841
integer_1_to_40  7.000000 20.392679 10.933851    10.933851  3.753159     3.753159
(cupy-env) PS C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST> python nist_wct_log_spectral_scan_FIXED.py --csv Fe_lines.csv --null-n 5000 --min-lines 100
C:\Users\Ricky.Reyes.CST\Desktop\RESEARCH\NIST\nist_wct_log_spectral_scan_FIXED.py:90: FutureWarning: DataFrame.applymap has been deprecated. Use DataFrame.map instead.
  return df.applymap(clean_cell)
[raw] rows=20615 cols=['element', 'sp_num', 'obs_wl_air(nm)', 'unc_obs_wl', 'ritz_wl_air(nm)', 'unc_ritz_wl', 'wn(cm-1)', 'intens', 'Aki(s^-1)', 'Acc', 'Ei(cm-1)', 'Ek(cm-1)', 'conf_i', 'term_i', 'J_i', 'conf_k', 'term_k', 'J_k', 'Type', 'tp_ref', 'line_ref', 'Unnamed: 21']
[clean] usable unique lines=19577
[scan] bins=160 k=0.5..80.0 n_k=2500
[null] 100/5000
[null] 200/5000
[null] 300/5000
[null] 400/5000
[null] 500/5000
[null] 600/5000
[null] 700/5000
[null] 800/5000
[null] 900/5000
[null] 1000/5000
[null] 1100/5000
[null] 1200/5000
[null] 1300/5000
[null] 1400/5000
[null] 1500/5000
[null] 1600/5000
[null] 1700/5000
[null] 1800/5000
[null] 1900/5000
[null] 2000/5000
[null] 2100/5000
[null] 2200/5000
[null] 2300/5000
[null] 2400/5000
[null] 2500/5000
[null] 2600/5000
[null] 2700/5000
[null] 2800/5000
[null] 2900/5000
[null] 3000/5000
[null] 3100/5000
[null] 3200/5000
[null] 3300/5000
[null] 3400/5000
[null] 3500/5000
[null] 3600/5000
[null] 3700/5000
[null] 3800/5000
[null] 3900/5000
[null] 4000/5000
[null] 4100/5000
[null] 4200/5000
[null] 4300/5000
[null] 4400/5000
[null] 4500/5000
[null] 4600/5000
[null] 4700/5000
[null] 4800/5000
[null] 4900/5000
[null] 5000/5000

[done]
{
  "k_best": 31.358343337334933,
  "deltaD": 107.12175266766758,
  "scan_null_p": 0.0001999600079984003,
  "tail_count_ge": 0,
  "null_n": 5000,
  "n_obs": 10.752733902312999,
  "n_unique_lines": 19577
}

[nearest branch]
         branch         n  k_target   k_error  abs_k_error   n_error  abs_n_error
integer_1_to_40 11.000000 32.079449 -0.721105     0.721105 -0.247266     0.247266
 koide_10_15_20 10.000000 29.163135  2.195208     2.195208  0.752734     0.752734
integer_1_to_40 10.000000 29.163135  2.195208     2.195208  0.752734     0.752734
integer_1_to_40 12.000000 34.995762 -3.637419     3.637419 -1.247266     1.247266
integer_1_to_40  9.000000 26.246822  5.111522     5.111522  1.752734     1.752734
integer_1_to_40 13.000000 37.912076 -6.553733     6.553733 -2.247266     2.247266
  folded_4over9 13.333333 38.884180 -7.525837     7.525837 -2.580599     2.580599
integer_1_to_40  8.000000 23.330508  8.027835     8.027835  2.752734     2.752734
integer_1_to_40 14.000000 40.828389 -9.470046     9.470046 -3.247266     3.247266
integer_1_to_40  7.000000 20.414195 10.944149    10.944149  3.752734     3.752734