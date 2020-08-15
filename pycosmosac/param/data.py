
Mullins_2006 = {
  "a_eff" : 7.5,
  "r_av" : 0.81764,
  "f_decay" : 1.0,
  "c_hb" : 85580.0,
  "sigma_hb" : 0.0084,
  "q_0" : 79.53,
  "r_0" : 66.69,
  "z" : 10,
  "alpha_prime" : 16466.72
}


Saidi_2002 = {
  "a_eff" : 6.22,
  "f_decay" : 3.57,
  "c_hb" : 42790.0,
  "sigma_hb" : 0.00764,
  "q_0" : 79.53,
  "r_0" : 66.69,
  "z" : 10,
  "alpha_prime" : 16466.72,
}

BIOSAC_SVP_GEPOL = {
  "a_eff" : 6.28852021e+00,
  "f_decay" : 3.57,
  "c_hb" : 7.17783028e+04,
  "sigma_hb" : 7.19127849e-03,
  "q_0" : 79.53,
  "r_0" : 66.69,
  "z" : 10,
  "alpha_prime" : 1.27029075e+04,
  "omega_ring" : 2.66273110e-01,
  "eta_0" : 1.20563529e+00
}

Hsieh_2010 = {
  "a_eff" : 7.25,
  "f_decay" : 3.57,
  "sigma_hb" : 0.0084, #should this be adjustable? this number is from https://dx.doi.org/10.1021/acs.jctc.9b01016
  "c_ohoh" : 4013.78,
  "c_otot" : 932.31,
  "c_ohot" : 3016.43,
  "sigma_0" : 0.007,
  "A_es" : 6525.69,
  "B_es" : 1.4859*1e8,
  "q_0" : 79.53,
  "r_0" : 66.69,
  "z" : 10
}


disp = {
  "C(sp3)": 115.7023,
  "C(sp2)": 117.4650,
  "C(sp)": 66.0691,
  "-O-": 95.6184,
  "=O": -11.0549,
  "N(sp3)": 15.4901,
  "N(sp2)": 84.6268,
  "N(sp)": 109.6621,
  "F": 52.9318,
  "Cl": 104.2534,
  "H(OH)": 19.3477,
  "H(NH)": 141.1709,
  "H(H2O)": 58.3301,
  "H(COOH)": 58.3301, #see http://dx.doi.org/10.1016/j.fluid.2014.01.032; this number is 19.3477 in https://dx.doi.org/10.1021/acs.jctc.9b01016
  "P": 105.0,
  "S": 105.0,
  "Br": 105.0,
  "I": 105.0,
  "w" : 0.27027
}

COSMORS = {
  "a_eff" : 6.25,
  "r_av" : 0.5,
  "f_decay" : 1.0,
  "c_hb" : 8808.0,
  "sigma_hb" : 0.0085,
  "q_0" : 79.53,
  "r_0" : 66.69,
  "alpha_prime" : 1428.0,
  "lambda_0" : 1.0,
  "lambda_1" : 1.0,
  "lambda_2" : 1.0,
  "omega_ring" : 0.2136,
  "eta_0" : 5.208
}

BIOSAC_SVP_GEPOL_disp = {
  "cdisp" : 1.54042645e+00,
  "H"  : -2.86015045e-02, 
  "C"  :  6.01727630e-03,
  "N"  : -6.69178136e-02,
  "O"  : -5.21667805e-02,
  "F"  : -1.87041166e-03,
  "P"  : -1.48010041e-01,
  "S"  : -4.53578156e-02,
  "Cl" : -1.43757687e-02,
  "Br" :  7.84788951e-03,
  "I"  : -4.49383370e-02
}

disp_RS = BIOSAC_SVP_GEPOL_disp
